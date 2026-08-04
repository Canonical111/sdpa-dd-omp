#!/usr/bin/env bash
# Build, verify and optionally install an sdpa-*-omp solver from the current checkout.
#
# Usage:  install.sh [dd|gmp|qd] [--prefix DIR] [--serial] [-h|--help]
#
# Every step was verified by execution on Ubuntu x86-64 and macOS Apple Silicon. Expected
# failures (the first SPOOLES pass) are handled in explicit branches; anything else dies
# with a message naming the log. Success is printed only after every requested operation
# succeeded -- a stale binary or a failed --prefix install cannot produce a false DONE.
#
# MODIFIED (GPLv2 2a notice), 2026-08-04: build QD with hardware FMA on aarch64 only, and
# namespace the QD cache by that choice. See git log.
set -euo pipefail

QD_VER=2.3.24
QD_URL="https://www.davidhbailey.com/dhbsoftware/qd-$QD_VER.tar.gz"
QD_SHA=a47b6c73f86e6421e86a883568dd08e299b20e36c11a99bdfbe50e01bde60e38
CACHE=${XDG_CACHE_HOME:-$HOME/.cache}/sdpa-omp

usage() {
    cat <<'EOF'
Build, verify and optionally install an sdpa-*-omp solver from the current checkout.

Usage:  install.sh [dd|gmp|qd] [--prefix DIR] [--serial] [-h|--help]

  dd|gmp|qd     which solver this checkout is (auto-detected; a mismatch is refused)
  --prefix DIR  after verification, install the binary into DIR/bin
  --serial      build without OpenMP (and verify none is linked)

Success ends with DONE; any failure ends with FATAL naming the log to read.
EOF
    exit "${1:-0}"
}
note() { echo "==> $*"; }

SOLVER="" PREFIX="" OPENMP=yes
while [ $# -gt 0 ]; do case "$1" in
  dd|gmp|qd) SOLVER=$1 ;;
  --prefix) { [ $# -ge 2 ] && [ -n "$2" ] && [ "${2#-}" = "$2" ]; } \
            || { echo "--prefix needs a directory argument" >&2; usage 2; }
            PREFIX=$2; shift ;;
  --prefix=*) PREFIX=${1#--prefix=}
            [ -n "$PREFIX" ] || { echo "--prefix needs a directory argument" >&2; usage 2; } ;;
  --serial) OPENMP=no ;;
  -h|--help) usage 0 ;;
  *) echo "unknown argument: $1" >&2; usage 2 ;;
esac; shift; done

# No temporary state before we know work is actually requested.
LOG=$(mktemp -d "${TMPDIR:-/tmp}/sdpa-omp-install.XXXXXX")
die() { echo "FATAL: $*" >&2; echo "       logs: $LOG" >&2; exit 1; }

# ---- detect the checkout, and refuse a solver/repository mismatch ----------
[ -f configure.ac ] || [ -f configure.in ] || die "run from the root of an sdpa-*-omp checkout"
DETECTED=""
for s in dd gmp qd; do
    grep -qi "sdpa[-_]$s" configure.ac configure.in 2>/dev/null && DETECTED=$s && break
done
[ -n "$DETECTED" ] || die "cannot identify this checkout as sdpa-dd/gmp/qd"
if [ -n "$SOLVER" ] && [ "$SOLVER" != "$DETECTED" ]; then
    die "you asked for '$SOLVER' but this checkout is sdpa-$DETECTED -- one solver's rules must not run in another's repository"
fi
SOLVER=$DETECTED
BIN=sdpa_$SOLVER
note "solver: sdpa-$SOLVER-omp   openmp: $OPENMP   logs: $LOG"

OS=$(uname -s)
NPROC=$( (nproc || sysctl -n hw.ncpu) 2>/dev/null | head -1 )
case "$NPROC" in ''|*[!0-9]*|0) NPROC=1 ;; esac

# ---- preflight: every tool this platform+solver will need -------------------
need() { command -v "$1" >/dev/null || die "missing required tool: $1  ($2)"; }
need make "build-essential / Xcode CLT"
if [ "$OS" = Darwin ]; then
    command -v brew >/dev/null || die "Homebrew required on macOS (https://brew.sh)"
    # TRAP: /usr/bin/gcc is Apple clang -- no OpenMP; configure would SILENTLY build serial.
    GCC=$(ls "$(brew --prefix gcc 2>/dev/null)"/bin/gcc-[0-9]* 2>/dev/null | head -1) || true
    GXX=$(ls "$(brew --prefix gcc 2>/dev/null)"/bin/g++-[0-9]* 2>/dev/null | head -1) || true
    [ -n "${GCC:-}" ] && [ -n "${GXX:-}" ] || die "Homebrew GCC required: brew install gcc  (Apple clang has no OpenMP)"
    need otool "Xcode CLT"
else
    GCC=gcc GXX=g++
    need gcc "apt install build-essential"; need g++ "apt install build-essential"
    if [ "$SOLVER" = qd ]; then need nm "binutils"; else need ldd "libc-bin"; fi
fi
if [ "$SOLVER" != qd ]; then need autoreconf "autoconf automake libtool"; fi
if [ "$SOLVER" = qd ]; then need curl "curl"; need tar "tar"; fi
if [ -n "$PREFIX" ]; then need install "coreutils"; fi
note "compiler: $GXX   jobs: $NPROC"

# ---- flags: explicit, auditable ---------------------------------------------
OPT_FLAGS="-O2 -funroll-all-loops"
OMP_FLAG=""
[ "$OPENMP" = yes ] && OMP_FLAG="-fopenmp"

# ---- per-solver configure arguments ------------------------------------------
CFG=()
case "$SOLVER" in
gmp)
    if [ "$OS" = Darwin ]; then
        # TRAP: the bundled GMP 6.2.1 test suite bus-errors on Apple Silicon.
        [ -f "$(brew --prefix)/include/gmp.h" ] || die "brew install gmp  (bundled GMP fails its tests on arm64)"
        CFG+=(--with-system-gmp --with-gmp-includedir="$(brew --prefix)/include" --with-gmp-libdir="$(brew --prefix)/lib")
    fi
    CFG+=(--enable-openmp="$OPENMP")
    ;;
dd) CFG+=(--enable-openmp="$OPENMP") ;;
qd)
    # TRAP: QD must be built with the SAME compiler family. Homebrew's qd bottle is
    # clang/libc++ and fails to link against g++ (missing ostream operator<<).
    QD_DIR=${QD_DIR:-}
    if [ -z "$QD_DIR" ]; then
        if [ "$OS" != Darwin ] && [ -f /usr/include/qd/qd_real.h ]; then
            QD_DIR=/usr
        else
            # TRAP: QD's configure runs its FMA probe only under `case $host in
            # powerpc*-*-*)`, so --enable-fma=auto resolves to "none" here and two_prod --
            # the most-executed primitive in the solver -- compiles the ~17-op Dekker split
            # instead of p=a*b; err=fma(a,b,-p). Ask for it explicitly, but on aarch64 ONLY:
            # there fmadd is baseline ISA and __builtin_fma is one instruction, whereas on
            # x86-64 without -mfma GCC lowers __builtin_fma to a libm CALL that is slower
            # than the split -- and generic x86-64 is a deliberate choice here so binaries
            # stay bit-identical across AMD and Intel nodes.
            QD_FMA_ARGS=()
            case "$(uname -m)" in
                aarch64|arm64) QD_FMA_ARGS=(--enable-fma=gnu); QD_FMA_TAG=fma ;;
                *)             QD_FMA_TAG=nofma ;;
            esac
            # Namespace the cache by ABI *and* by the FMA decision: an old library built by
            # a different compiler, architecture, or two_prod implementation must never be
            # silently reused. Bumping the tag is what forces a rebuild of caches populated
            # before 2026-08-04, which contain a non-FMA QD.
            QD_DIR=$CACHE/qd-$QD_VER-$(uname -m)-gcc$("$GXX" -dumpversion | cut -d. -f1)-$QD_FMA_TAG
            # .complete is written only after make install succeeds -- a bare
            # libqd.a left by an interrupted install does not count.
            if [ ! -f "$QD_DIR/.complete" ]; then
                command -v sha256sum >/dev/null || command -v shasum >/dev/null \
                    || die "need sha256sum or shasum to verify the QD download"
                note "building QD $QD_VER from source with $GXX (cached at $QD_DIR)"
                rm -rf "$QD_DIR"
                BUILD=$(mktemp -d "${TMPDIR:-/tmp}/qd-build.XXXXXX")
                curl -fsSLo "$BUILD/qd.tar.gz" "$QD_URL" \
                    || { rm -rf "$BUILD"; die "cannot download QD; set QD_DIR to an existing GCC-built QD"; }
                GOT=$( { sha256sum "$BUILD/qd.tar.gz" 2>/dev/null || shasum -a 256 "$BUILD/qd.tar.gz"; } | awk '{print $1}')
                [ "$GOT" = "$QD_SHA" ] || { rm -rf "$BUILD"; die "QD tarball sha256 mismatch: got $GOT, expected $QD_SHA"; }
                tar xzf "$BUILD/qd.tar.gz" -C "$BUILD" \
                    || { rm -rf "$BUILD"; die "cannot extract the verified QD tarball"; }
                if ( cd "$BUILD/qd-$QD_VER" && ./configure CC="$GCC" CXX="$GXX" \
                        --prefix="$QD_DIR" --enable-fortran=no "${QD_FMA_ARGS[@]+"${QD_FMA_ARGS[@]}"}" \
                        && make -j"$NPROC" && make install
                   ) >"$LOG/qd-build.log" 2>&1; then
                    touch "$QD_DIR/.complete"
                    rm -rf "$BUILD"
                else
                    rm -rf "$QD_DIR"   # never leave a half-installed cache to be reused
                    die "QD build failed; see $LOG/qd-build.log"
                fi
            fi
        fi
    fi
    LIBSUB=lib; [ -d "$QD_DIR/lib/x86_64-linux-gnu" ] && LIBSUB=lib/x86_64-linux-gnu
    CFG+=(--with-qd-includedir="$QD_DIR/include" --with-qd-libdir="$QD_DIR/$LIBSUB")
    # TRAP: the 2009 config.guess predates arm64.
    if [ "$OS" = Darwin ]; then
        AM=$(ls -d "$(brew --prefix)"/share/automake-* 2>/dev/null | head -1)
        [ -n "$AM" ] || die "automake required: brew install automake  (arm64-aware config.guess)"
        cp "$AM"/config.guess "$AM"/config.sub . && chmod u+w config.guess config.sub
    fi
    ;;
esac

# ---- configure ----------------------------------------------------------------
note "configure"
if [ "$SOLVER" = qd ]; then
    # TRAP: sdpa-qd's --enable-openmp is vestigial; -fopenmp must go in the flags.
    ./configure CC="$GCC" CXX="$GXX" \
        CXXFLAGS="$OPT_FLAGS $OMP_FLAG" CFLAGS="$OPT_FLAGS $OMP_FLAG" LDFLAGS="$OMP_FLAG" \
        "${CFG[@]}" >"$LOG/configure.log" 2>&1 \
        || die "configure failed; see $LOG/configure.log"
else
    if [ ! -x configure ]; then
        note "autoreconf (upstream ships only configure.ac)"
        autoreconf -fi >"$LOG/autoreconf.log" 2>&1 || die "autoreconf failed; see $LOG/autoreconf.log"
    fi
    ./configure CC="$GCC" CXX="$GXX" "${CFG[@]}" >"$LOG/configure.log" 2>&1 \
        || die "configure failed; see $LOG/configure.log"
    if [ "$OPENMP" = yes ] && grep -q "support OpenMP... unsupported" "$LOG/configure.log"; then
        die "configure found NO OpenMP support -- this would be a SILENT serial build"
    fi
fi

# ---- build --------------------------------------------------------------------
# Objects from an earlier run must not survive into this one: make does not rebuild
# on a compiler/flag/OpenMP-mode change, and mixed objects fail to link (or worse).
note "make clean (a mode or compiler switch must not reuse old objects)"
make clean >"$LOG/clean.log" 2>&1 || die "make clean failed; see $LOG/clean.log"
# A binary from an earlier run must not be able to vouch for this one.
rm -f "$BIN"

spooles_rescue() {
    local spd
    case "$SOLVER" in qd) spd=spooles/build ;; *) spd=external/spooles/work/internal ;; esac
    [ -f "$spd/Make.inc" ] || die "build failed before SPOOLES was extracted; see $LOG/make1.log"
    note "SPOOLES rescue (its .POSIX Make.inc selects 'c99', which rejects the flags)"
    # -D_GNU_SOURCE: some toolchains hide struct timezone even under gcc.
    sed -i.bak -e "s|^# CC = gcc|  CC = $GCC|" \
        -e "s|^  CFLAGS = |  CFLAGS = -D_GNU_SOURCE |" \
        -e "s|^  CFLAGS += -O2 -funroll-all-loops|  CFLAGS += -O2 -funroll-all-loops -D_GNU_SOURCE -Wno-error=int-conversion -Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types|" \
        "$spd/Make.inc"
    rm -f "$spd/Make.inc.bak"
    ( cd "$spd" && find . -name '*.o' -delete && rm -f spooles.a libspooles.a \
        && make global -f makefile ) >"$LOG/spooles.log" 2>&1 \
        || die "SPOOLES rebuild failed; see $LOG/spooles.log"
    if [ "$SOLVER" = qd ]; then
        # TRAP: a retried make re-untars SPOOLES over any fix; satisfying the
        # libspooles.a target prevents the re-extraction.
        cp "$spd/spooles.a" "$spd/libspooles.a"
    else
        mkdir -p external/i/SPOOLES/lib && cp "$spd/spooles.a" external/i/SPOOLES/lib/libspooles.a
    fi
}

note "make (a first-pass stop inside SPOOLES is expected on many systems)"
if make -j"$NPROC" >"$LOG/make1.log" 2>&1; then
    :
else
    spooles_rescue
    make -j"$NPROC" >"$LOG/make2.log" 2>&1 || die "build failed after SPOOLES rescue; see $LOG/make2.log"
fi
[ -x "$BIN" ] || die "make reported success but produced no $BIN; see $LOG"

# ---- verify --------------------------------------------------------------------
note "verifying"
if [ "$OS" = Darwin ]; then HAVE_OMP=$(otool -L "$BIN" | grep -c gomp || true)
elif [ "$SOLVER" = qd ]; then HAVE_OMP=$(nm "$BIN" | grep -c GOMP || true)   # static: ldd is blind
else HAVE_OMP=$(ldd "$BIN" | grep -c gomp || true); fi
[ "$OPENMP" = yes ] && [ "$HAVE_OMP" -eq 0 ] && die "binary has NO OpenMP runtime -- silent serial build"
[ "$OPENMP" = no ] && [ "$HAVE_OMP" -gt 0 ] && die "serial requested but OpenMP is linked"

./"$BIN" -ds example1.dat-s -o "$LOG/example.result" -p param.sdpa >"$LOG/example.log" 2>&1 \
    || die "solver failed on the bundled example; see $LOG/example.log"
OBJ=$(grep -m1 objValPrimal "$LOG/example.result" 2>/dev/null | sed 's/.*= *//' || true)
[ -n "$OBJ" ] || die "no objValPrimal in the example output; see $LOG/example.log and $LOG/example.result"
case "$OBJ" in -4.19*|-4.1899*) : ;; *) die "example objective $OBJ != expected -4.19e+01" ;; esac
note "example solved: objValPrimal = $OBJ   (OpenMP runtime markers: $HAVE_OMP)"

if [ -n "$PREFIX" ]; then
    install -d "$PREFIX/bin" || die "cannot create $PREFIX/bin"
    install -m 0755 "$BIN" "$PREFIX/bin/$BIN" || die "cannot install to $PREFIX/bin/$BIN"
    note "installed $PREFIX/bin/$BIN"
fi

rm -rf "$LOG"
note "DONE. Run with OMP_NUM_THREADS=<physical cores>; pin with taskset/OMP_PLACES=cores on Linux."
