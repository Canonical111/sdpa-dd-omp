# MODIFIED from upstream (GPLv2 2a notice), 2026-08-04: pass $(QD_FMA_ARG) to QD's configure
# so the hardware fused-multiply-add two_prod is used on aarch64. See git log.
#
# QD's configure probes for FMA only on powerpc, so its --enable-fma=auto default leaves
# QD_FMA/QD_FMS undefined here and two_prod compiles the ~17-op Dekker split. QD_FMA_ARG is
# set by our configure.ac to --enable-fma=gnu when $host_cpu is aarch64/arm64 and to the
# EMPTY string on every other host -- notably x86-64, where __builtin_fma without -mfma
# becomes a libm call and would be slower than the split. When it is empty the argument
# vector handed to ./configure below is exactly what it was before this change.
QDVERSION=2.3.23

if IS_MACOS
CONFIGURE_ARGS=--prefix=$(abs_builddir)/../i/QD/ --host=$(host) $(QD_FMA_ARG)
else
CONFIGURE_ARGS=--prefix=$(abs_builddir)/../i/QD/ --host=$(host) --disable-shared $(QD_FMA_ARG)
endif

EXTRACT_STAMP=.qd_extract_done
PATCH_STAMP=.qd_patch_done
CONFIGURE_STAMP=.qd_configure_done
BUILD_STAMP=.qd_build_done
INSTALL_STAMP=.qd_install_done
INSTALL_INTERNAL_STAMP=.qd_install_internal_done
CHECK_STAMP=.qd_check_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal
	cd $(abs_builddir)/work/internal; tar xvfz ../../download/qd-$(QDVERSION).tar.gz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	touch $(abs_builddir)/work/$(PATCH_STAMP)

# MODIFIED (GPLv2 2a notice), 2026-08-04: forward the FP-contraction pin into this
# sub-build. Until now NO C/C++ flags reached QD's configure at all -- the env below listed
# CC/CXX/CPPFLAGS/F77/LDFLAGS/NM/RANLIB and nothing else -- so libqd was compiled with
# autoconf's bare "-g -O2" and kept the compiler's DEFAULT -ffp-contract=fast even when the
# SDPA translation units were built with -ffp-contract=off. dd_real arithmetic is
# header-inline and libqd.a also carries out-of-line dd/qd routines (sqrt(dd_real) on the
# Lanczos/Rpotf2 path), so a pin that stops at the SDPA TUs is only half applied: the
# resulting numbers would look like a numerical finding rather than the build bug they are.
#
# QD_FLAGS_ENV comes from our configure.ac. It is EMPTY unless --enable-fp-contract was
# given, in which case it is exactly
#     "CFLAGS=-g -O2 <pin>" "CXXFLAGS=-g -O2 <pin>"
# i.e. autoconf's own default flags -- which is what QD would have chosen for itself, since
# QD's configure.ac adds no optimisation flags of its own -- plus the pin. Only the pin
# differs between a pinned and an unpinned libqd; the optimisation level does not move.
# When it is empty the environment handed to ./configure is byte-for-byte the old one and
# autoconf supplies "-g -O2" itself, so an unpinned build is bit-for-bit the old build.
#
# Depend on Makefile as well as the patch stamp: CONFIGURE_ARGS and QD_FLAGS_ENV live here,
# so without this an existing build tree keeps a stamp made with the OLD arguments and
# silently produces a library configured the old way. That is exactly how a QD built WITHOUT
# hardware FMA can be mistaken for one built with it.
#
# "$(MAKE) clean" before reconfiguring is not belt-and-braces: automake does NOT make object
# files depend on the Makefile, so a reconfigured tree with new flags would happily keep the
# .lo/.o files compiled with the old ones. Removing the install directory too means a build
# that fails after this point can never leave a stale libqd for the SDPA link to pick up.
work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP) Makefile
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; test ! -f Makefile || ${MAKE} clean
	rm -rf $(abs_builddir)/../i/QD
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" $(QD_FLAGS_ENV) "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS) --enable-fortran=no
	@qdmk=$(abs_builddir)/work/internal/qd-$(QDVERSION)/src/Makefile; \
	qdflags=`sed -n 's/^CXXFLAGS = //p' $$qdmk`; \
	echo "QD sub-build CXXFLAGS = $$qdflags"; \
	case " $$qdflags " in *" -O2 "*) ;; \
	  *) echo "*** QD sub-build lost -O2 (CXXFLAGS = $$qdflags in $$qdmk)" >&2; exit 1 ;; esac; \
	if test -n "$(FPCONTRACT_FLAG)"; then \
	  case " $$qdflags " in *" $(FPCONTRACT_FLAG) "*) ;; \
	    *) echo "*** $(FPCONTRACT_FLAG) did NOT reach the QD sub-build (CXXFLAGS = $$qdflags in $$qdmk)" >&2; exit 1 ;; esac; \
	fi
	@qdcfg=$(abs_builddir)/work/internal/qd-$(QDVERSION)/include/qd/qd_config.h; \
	if grep -q '^#define QD_FMA(' $$qdcfg; then got=yes; else got=no; fi; \
	echo "QD hardware FMA in $$qdcfg: $$got (expected $(QD_EXPECT_FMA))"; \
	if test "$(QD_EXPECT_FMA)" = unknown; then \
	  echo "    (--enable-qd-fma=auto on a non-aarch64 host: QD ran its own probe, so this"; \
	  echo "     is REPORTED, not asserted. QD's probe can enable FMA on powerpc*, and"; \
	  echo "     aborting there would fail a legitimate build on an untested host.)"; \
	else \
	  test "$$got" = "$(QD_EXPECT_FMA)" || \
	    { echo "*** QD FMA is '$$got' but $(QD_FMA_ARG) asked for '$(QD_EXPECT_FMA)'" >&2; exit 1; }; \
	fi
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE}
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
#	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE} check
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; ${MAKE} install
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

install: work/$(INSTALL_STAMP)

distdir:

clean:
	rm -rf work
distclean:
	rm -rf work
