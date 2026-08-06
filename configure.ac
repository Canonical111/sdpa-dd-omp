dnl
dnl Configure for SDPA-DD 7.1.3
dnl
dnl MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: QD_FMA_ARG -- enable QD's
dnl hardware fused-multiply-add path on aarch64 only. See git log.

AC_INIT(sdpa-dd, 7.1.3)
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_CC
AC_PROG_CXX
AC_CANONICAL_HOST
AC_OPENMP
AC_LANG([C++])
AC_OPENMP

highest_abs_top_srcdir=`cd "$srcdir" && pwd`

AC_ARG_WITH(system-spooles,
[  --with-system-spooles      Use spooles already on system],
WITH_SYSTEM_SPOOLES=$withval, WITH_SYSTEM_SPOOLES=no)
AM_CONDITIONAL(BUILD_SPOOLES, test x$WITH_SYSTEM_SPOOLES = xno)

AC_ARG_WITH(spooles-includedir,
[  --with-spooles-includedir    Gives spooles include dir you use(e.g., /usr/local/include/spooles).],
SPOOLES_INCLUDEDIR=$withval
)

AC_ARG_WITH(spooles-libdir,
[  --with-spooles-libdir    Gives spooles lib dir you use(e.g., /usr/local/lib).],
SPOOLES_LIBDIR=$withval
)

AC_ARG_WITH(system-qd,
[  --with-system-qd      Use qd already on system],
WITH_SYSTEM_QD=$withval, WITH_SYSTEM_QD=no)
AM_CONDITIONAL(BUILD_QD, test x$WITH_SYSTEM_QD = xno)

AC_ARG_WITH(qd-includedir,
[  --with-qd-includedir    Gives qd include dir you use(e.g., /usr/local/include).],
QD_INCLUDEDIR=$withval
)

AC_ARG_WITH(qd-libdir,
[  --with-qd-libdir    Gives qd lib dir you use(e.g., /usr/local/lib).],
QD_LIBDIR=$withval
)

AC_ARG_WITH(system-mplapack,
[  --with-system-mplapack      Use mplapack already on system],
WITH_SYSTEM_MPLAPACK=$withval, WITH_SYSTEM_MPLAPACK=no)
AM_CONDITIONAL(BUILD_MPLAPACK, test x$WITH_SYSTEM_MPLAPACK = xno)

AC_ARG_WITH(mplapack-includedir,
[  --with-mplapack-includedir    Gives mplapack include dir you use(e.g., /usr/local/include).],
MPLAPACK_INCLUDEDIR=$withval
)

AC_ARG_WITH(mplapack-libdir,
[  --with-mplapack-libdir    Gives mplapack lib dir you use(e.g., /usr/local/lib).],
MPLAPACK_LIBDIR=$withval
)

AC_ARG_ENABLE(metis,
[  --enable-metis    Build with metis],
[case "${enableval}" in
  yes) with_metis=yes ;;
  no)  with_metis=no  ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-metis) ;;
esac],[with_metis=no])
AM_CONDITIONAL(ENABLE_METIS, test x$with_metis = xyes)

dnl spooles requires -lm when linking.
LIBS="-lm"

IS_LINUX=0
case "$host" in
    *linux*)
        IS_LINUX=1
        ;;
    *)
esac
AM_CONDITIONAL(IS_LINUX, test x$IS_LINUX = x1)

IS_MINGW=0
case "$host" in
        *-mingw32*) # mingw
          IS_MINGW=1
        ;;
        *)
esac
AM_CONDITIONAL(IS_MINGW, test x$IS_MINGW = x1)

IS_MACOS=0
case "$host" in
         *apple-darwin*) # MacOSX
          IS_MACOS=1
        ;;
        *)
esac
AM_CONDITIONAL(IS_MACOS, test x$IS_MACOS = x1)

dnl Hardware FMA for the bundled QD build -- aarch64 ONLY.
dnl
dnl QD's own configure runs its FMA probe only under `case $host in powerpc*-*-*)`, so the
dnl default --enable-fma=auto silently resolves to "none" on aarch64 and x86-64 alike:
dnl QD_FMA/QD_FMS stay undefined and two_prod -- the most-executed primitive in the solver --
dnl compiles the Dekker two-way-split fallback (~17 double ops) instead of p=a*b; err=fma(a,b,-p).
dnl
dnl The gate is aarch64-only on purpose, and it is NOT a portability oversight:
dnl   - on aarch64, fmadd is baseline ISA, so __builtin_fma is a single instruction;
dnl   - on x86-64 WITHOUT -mfma, GCC lowers __builtin_fma to a libm CALL, which is far slower
dnl     than the Dekker split. This project builds generic x86-64 deliberately (no -march,
dnl     no -mfma) so binaries stay bit-identical across AMD and Intel nodes.
dnl On every non-aarch64 host QD_FMA_ARG expands empty and the arguments handed to QD's
dnl configure are unchanged, so published x86-64 state hashes cannot move.
dnl --enable-qd-fma=auto|gnu|no overrides the host gate. The DEFAULT is auto, which is
dnl exactly the host gate described above, so the shipped build is unchanged. The explicit
dnl values exist so that "same build, QD's explicit __builtin_fma turned off" can be
dnl configured without hand-editing this file -- hand-editing, or overriding QD_FMA_ARG on a
dnl make command line, does NOT invalidate external/qd's configure stamp and would silently
dnl reuse the previous library.
AC_ARG_ENABLE(qd-fma,
[  --enable-qd-fma=MODE       QD hardware FMA: auto (aarch64 only), gnu, no],
[qd_fma_mode=$enableval], [qd_fma_mode=auto])

case "$qd_fma_mode" in
  auto)
    case "$host_cpu" in
      aarch64|arm64) QD_FMA_ARG="--enable-fma=gnu" ;;
      *)             QD_FMA_ARG="" ;;
    esac ;;
  gnu)     QD_FMA_ARG="--enable-fma=gnu" ;;
  no|none) QD_FMA_ARG="--enable-fma=no" ;;
  *) AC_MSG_ERROR([bad value ${qd_fma_mode} for --enable-qd-fma (use auto, gnu or no)]) ;;
esac
AC_SUBST([QD_FMA_ARG])

dnl What the generated include/qd/qd_config.h must then look like. external/qd/Makefile.am
dnl checks this after the sub-configure runs: the FMA decision is invisible in the SDPA
dnl build otherwise, and a QD built without it has already once been mistaken for one built
dnl with it.
if test x"$QD_FMA_ARG" = x"--enable-fma=gnu"; then
  QD_EXPECT_FMA=yes
else
  QD_EXPECT_FMA=no
fi
AC_SUBST([QD_EXPECT_FMA])

dnl How thread we link thread library esp. needed for spooles.
case "$host" in
        *-*-linux*|*-*-uclinux*)
            pthread_cflags="-D_REENTRANT"
            pthread_lib="-lpthread"
            ;;
        *-*-freebsd*|*-*-dragonfly*)
            pthread_cflags=""
            pthread_lib="-pthread"
        ;;
esac
LIBS="$LIBS $pthread_lib"
CFLAGS="$CFLAGS $pthread_cflags"
CXXFLAGS="$CXXFLAGS $pthread_cflags"
AC_SUBST(pthread_lib)
dnl 2026-08-05: pthread_cflags was referenced by Makefile.am (sdpa_dd_CXXFLAGS) and by the
dnl SPOOLES patch-Make.inc, but was never AC_SUBST'ed, so both expanded to nothing. Harmless
dnl here -- CXXFLAGS is preserved by this fork and already carries the flag -- but the
dnl SPOOLES sub-build reads neither, and external/spooles/Makefile.am now forwards this
dnl variable on the SPOOLES make line, so it has to actually exist.
AC_SUBST(pthread_cflags)

CFLAGS="$CFLAGS -O2 -funroll-all-loops"
CXXFLAGS="$CXXFLAGS -O2 -funroll-all-loops"
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)

dnl ---------------------------------------------------------------------------
dnl Floating-point contraction pin (--enable-fp-contract=off|on|fast).
dnl
dnl Nothing in this build pinned -ffp-contract before this option existed, and GCC/clang
dnl default to =fast for C++. dd_real arithmetic is header-inline, so the COMPILER decides
dnl whether the compensation terms of the double-double algorithms get folded into an FMA.
dnl On baseline x86-64 there is no FMA instruction to fold into, so those binaries behave as
dnl if contraction were off; on aarch64 fmadd is baseline ISA, so they do not, and the two
dnl platforms then execute different arithmetic. This option makes the semantics explicit.
dnl
dnl It appends the flag LAST, after -O2 -funroll-all-loops, so it wins over anything earlier,
dnl and it is applied to CFLAGS and CXXFLAGS, i.e. to every SDPA and MPLAPACK translation
dnl unit. It is also forwarded to the bundled QD sub-build -- see external/qd/Makefile.am;
dnl a half-applied pin (SDPA TUs pinned, libqd.a not) would look like a numerical result
dnl instead of the build bug it is.
dnl
dnl DEFAULT IS "off", i.e. this fork ships with contraction pinned. 2026-08-05.
dnl
dnl That default is a deliberate change, and it is what makes the arm64 and x86-64 builds of
dnl this fork the same solver numerically. Measured on the complete published Mac problem set,
dnl an M1 build with the pin matches the published x86-64 reference on 20 of 20 problems --
dnl iteration counts AND full-precision objectives -- where the unpinned build matches on none
dnl of them at the state-hash level. Nine iteration counts and three objectives move, and every
dnl one of them moves ONTO the x86 value; none moves away and none lands anywhere else. See
dnl results/b4_mac_rebaseline.tsv.
dnl
dnl Two of those disagreements were invisible to the comparison this project had been using to
dnl check for cross-platform agreement: mcp100 changed 66 -> 54 iterations at an UNCHANGED
dnl objective, and qap9 kept 19 iterations while its objective moved at the 8th significant
dnl digit. An iteration check misses qap9; an objective check misses mcp100.
dnl
dnl Cost: 1.0% per iteration at 1 thread, 1.7% at 8, against a per-cell repeat spread of
dnl 0.1-3.1% on the measuring machine -- at or below resolution. Total iterations over 12
dnl problems went 630 -> 626, a small net saving. (The audit's estimate of 9-10% does not
dnl reproduce: gpp124-1 measured -0.1%, arch0 +2.3%.) A1's hardware FMA is retained in full:
dnl -ffp-contract=off does not disable an explicit __builtin_fma, verified at the instruction
dnl level -- 143 fused instructions survive in libqd under the pin, versus 0 with FMA disabled,
dnl and the pinned-no-FMA build is 1.69x SLOWER.
dnl
dnl x86-64 is provably unaffected: all 76 translation units compile to byte-identical .text
dnl with and without the flag, and the unpinned build contains ZERO fused FMA instructions --
dnl there is nothing on baseline x86-64 for contraction to fold into. So the published x86
dnl benchmark and state-hash tables are untouched by this default.
dnl
dnl Pass --enable-fp-contract=fast to recover the old unpinned behaviour, or =default to add no
dnl flag at all.
AC_ARG_ENABLE(fp-contract,
[  --enable-fp-contract=MODE  pin FP contraction: off|on|fast|default (default: off)],
[fp_contract_mode=$enableval], [fp_contract_mode=off])

case "$fp_contract_mode" in
  default)      FPCONTRACT_FLAG="" ;;
  off|on|fast)  FPCONTRACT_FLAG="-ffp-contract=$fp_contract_mode" ;;
  *) AC_MSG_ERROR([bad value ${fp_contract_mode} for --enable-fp-contract (use off, on or fast)]) ;;
esac

if test x"$FPCONTRACT_FLAG" != x; then
  dnl A flag the compiler rejects must stop configure here, not silently vanish.
  AC_MSG_CHECKING([whether $CXX accepts $FPCONTRACT_FLAG])
  save_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS $FPCONTRACT_FLAG"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[double f(double a,double b,double c){return a*b+c;}]],[[]])],
    [AC_MSG_RESULT([yes])],
    [CXXFLAGS="$save_CXXFLAGS"
     AC_MSG_RESULT([no])
     AC_MSG_ERROR([$CXX does not accept $FPCONTRACT_FLAG])])
  CFLAGS="$CFLAGS $FPCONTRACT_FLAG"
fi
AC_SUBST([FPCONTRACT_FLAG])

dnl QD's configure takes its C/C++ flags from autoconf's AC_PROG_CC/AC_PROG_CXX default
dnl ("-g -O2"); QD's own configure.ac adds no optimisation flags of its own. external/qd's
dnl sub-configure is therefore handed exactly that default PLUS the pin, so the only
dnl difference between a pinned and an unpinned libqd is the contraction flag itself.
dnl QD_FLAGS_ENV is EMPTY when nothing is pinned, which leaves the sub-configure environment
dnl exactly as it was (autoconf then supplies "-g -O2" itself).
QD_DEFAULT_FLAGS="-g -O2"
if test x"$FPCONTRACT_FLAG" = x; then
  QD_FLAGS_ENV=""
else
  QD_FLAGS_ENV="\"CFLAGS=$QD_DEFAULT_FLAGS $FPCONTRACT_FLAG\" \"CXXFLAGS=$QD_DEFAULT_FLAGS $FPCONTRACT_FLAG\""
fi
AC_SUBST([QD_FLAGS_ENV])

dnl some dummy settings
if test x$SPOOLES_INCLUDEDIR = x; then
SPOOLES_INCLUDEDIR="."
fi
if test x$SPOOLES_LIBDIR = x; then
SPOOLES_LIBDIR="."
fi

dnl

if test x$WITH_SYSTEM_SPOOLES != xno; then
  dnl Check for SPOOLES include and library files
  HAVE_SPOOLES_HEADER=""
  HAVE_SPOOLES=""
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$SPOOLES_INCLUDEDIR"
  AC_CHECK_HEADER(SPOOLES.h, HAVE_SPOOLES_HEADER=yes)
  CPPFLAGS=$save_CPPFLAGS

  if test x$HAVE_SPOOLES_HEADER = xyes; then
    save_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS -L$SPOOLES_LIBDIR"
    AC_CHECK_LIB(spooles, InpMtx_init, [SPOOLES_LIB="-lspooles"; HAVE_SPOOLES=yes])
    LDFLAGS=$save_LDFLAGS
  fi

  if test x$HAVE_SPOOLES != xyes; then
    echo "Please install SPOOLES"
    echo "SPOOLES http://www.netlib.org/linalg/spooles/spooles.2.2.html"
    AC_MSG_ERROR([No SPOOLES found])
  fi
else
SPOOLES_LIBDIR=$highest_abs_top_srcdir/external/i/SPOOLES/lib
SPOOLES_INCLUDEDIR=$highest_abs_top_srcdir/external/spooles/work/internal #just too messy
fi
AC_SUBST(SPOOLES_LIBDIR)
AC_SUBST(SPOOLES_INCLUDEDIR)

dnl Check for QD include and library files
if test x$WITH_SYSTEM_QD != xno; then
  HAVE_QD_HEADER=""
  HAVE_QD=""
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS -I$QD_INCLUDEDIR"
  AC_CHECK_HEADER(qd/qd_real.h, HAVE_QD_HEADER=yes)
  CPPFLAGS=$save_CPPFLAGS

  if test x$HAVE_QD_HEADER = xyes; then
    save_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS -L$QD_LIBDIR"
    AC_CHECK_LIB(qd, c_qd_sqrt, [QD_LIBS="-lqd"; HAVE_QD=yes])
    LDFLAGS=$save_LDFLAGS
  fi

  if test x$HAVE_QD != xyes; then
    echo "Please install QD"
    echo "You can download Double-Double and Quad-Double Arithmetic at https://www.davidhbailey.com/dhbsoftware/"
    AC_MSG_ERROR([No QD found])
  fi
else
  QD_LIBDIR=$highest_abs_top_srcdir/external/i/QD/lib
  QD_INCLUDEDIR=$highest_abs_top_srcdir/external/i/QD/include
fi
AC_SUBST(QD_LIBDIR)
AC_SUBST(QD_INCLUDEDIR)

dnl Check for MPLAPACK include and library files
dnl some dummy settings
if test x$MPLAPACK_INCLUDEDIR = x; then
MPLAPACK_INCLUDEDIR="."
fi
if test x$MPLAPACK_LIBDIR = x; then
MPLAPACK_LIBDIR="."
fi

if test x$WITH_SYSTEM_MPLAPACK != xno; then
HAVE_MPLAPACK_HEADER=""
HAVE_MPLAPACK=""
save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS -I$MPLAPACK_INCLUDEDIR"
AC_LANG_PUSH([C++])
AC_CHECK_HEADER(mpblas_dd.h, HAVE_MPLAPACK_HEADER=yes)
CPPFLAGS=$save_CPPFLAGS

if test x$HAVE_MPLAPACK_HEADER = xyes; then
  save_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS -L$MPLAPACK_LIBDIR"
  LIBS="$LDFLAGS -L$MPLAPACK_LIBDIR -lmplapack_dd -lmpblas_dd_opt"
 # AC_CHECK_FUNC(Rsyev, [HAVE_MPLAPACK=yes])
  HAVE_MPLAPACK=yes
  if test x$IS_LINUX = x1; then
    LDFLAGS="$save_LDFLAGS -Wl,--rpath=$MPLAPACK_LIBDIR"
  fi
fi

if test x$HAVE_MPLAPACK != xyes; then
  echo "Please install MPLAPACK (MPBLAS/MPLAPACK)"
  echo "https://github.com/nakatamaho/mplapack/"
  AC_MSG_ERROR([No MPLAPACK found])
fi
else
MPLAPACK_LIBDIR=$highest_abs_top_srcdir/mplapack
MPLAPACK_INCLUDEDIR=$highest_abs_top_srcdir/mplapack
fi
AM_CONDITIONAL(SYSTEM_MPLAPACK, test x$with_system_mplapack = xyes)
AC_SUBST(MPLAPACK_LIBDIR)
AC_SUBST(MPLAPACK_INCLUDEDIR)

AC_CONFIG_FILES([
external/Makefile
external/spooles/Makefile
external/qd/Makefile
])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
