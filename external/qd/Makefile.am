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

# Depend on Makefile as well as the patch stamp: CONFIGURE_ARGS lives here, so without
# this an existing build tree keeps a stamp made with the OLD arguments and silently
# produces a library configured the old way. That is exactly how a QD built WITHOUT
# hardware FMA can be mistaken for one built with it.
work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP) Makefile
	cd $(abs_builddir)/work/internal/qd-$(QDVERSION) ; env "CC=$(CC)" "CXX=$(CXX)" "CPPFLAGS=$(CPPFLAGS)" "F77=$(F77)" "LDFLAGS=$(LDFLAGS)" "NM=$(NM)" "RANLIB=$(RANLIB)" ./configure $(CONFIGURE_ARGS) --enable-fortran=no
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
