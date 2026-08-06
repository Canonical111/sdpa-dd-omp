SPOOLESVERSION=2.2

EXTRACT_STAMP=.spooles_extract_done
PATCH_STAMP=.spooles_patch_done
CONFIGURE_STAMP=.spooles_configure_done
BUILD_STAMP=.spooles_build_done
INSTALL_STAMP=.spooles_install_done
INSTALL_INTERNAL_STAMP=.spooles_install_internal_done
CHECK_STAMP=.spooles_check_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal
	cd $(abs_builddir)/work/internal; tar xvfz ../../download/spooles.$(SPOOLESVERSION).tgz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-A2+makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-ETree+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-IVL+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MPI+makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MPI+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MT+drivers+AllInOneMT.c
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MT+makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-Make.inc
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-Tree+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-Utilities+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-timings.h
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

# 2026-08-05: CC is supplied on the SPOOLES make line, and the build log is then checked.
#
# WHY THE MAKE LINE.  SPOOLES's Make.inc opens with ".POSIX:", which makes GNU make default
# CC to c99 rather than cc, and every CC assignment in Make.inc is commented out.  Apple's
# c99 rejects every -f option, so the shipped recipe died on its first compile
# ("c99: illegal option -- f") and NO fork could be built by `make` alone -- only by
# install.sh's spooles_rescue(), which is why the working recipe had drifted out of the
# tree and into CI YAML and a skill script.  A command-line assignment outranks every
# assignment inside the makefiles, so this is also the only way to guarantee that the
# compiler ./configure selected is the compiler that builds SPOOLES.
#
# WHY THE ASSERTIONS.  A silently-empty flag variable has already cost this project a CI
# job that passed while testing nothing (configure-time CXXFLAGS discarded; 0 of 3671
# compile lines carried the flag under test).  So the log is kept, the compile lines are
# required to carry the patched CFLAGS, and spooles.a is required to be non-empty -- a
# failed SPOOLES build leaves a short but non-zero archive behind, so `test -s` alone is
# not enough and the exit status is checked first.
work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/ ; ${MAKE} global -f makefile CC='$(CC)' \
	  pthread_cflags='$(pthread_cflags)' \
	  > $(abs_builddir)/work/spooles-build.log 2>&1 \
	  || { cat $(abs_builddir)/work/spooles-build.log >&2 ; \
	       echo "SPOOLES build FAILED; full log: $(abs_builddir)/work/spooles-build.log" >&2 ; \
	       exit 1 ; }
	@grep -q -- '-funroll-all-loops' $(abs_builddir)/work/spooles-build.log \
	  || { echo "ERROR: no SPOOLES compile line carries the patched CFLAGS -- patch-Make.inc did not take effect" >&2 ; exit 1 ; }
	@grep -q -- '-Wno-error=int-conversion' $(abs_builddir)/work/spooles-build.log \
	  || { echo "ERROR: no SPOOLES compile line carries the GCC-14 -Wno-error flags" >&2 ; exit 1 ; }
	@test -s $(abs_builddir)/work/internal/spooles.a \
	  || { echo "ERROR: SPOOLES reported success but produced no spooles.a" >&2 ; exit 1 ; }
	@echo "SPOOLES built with CC=$(CC); log $(abs_builddir)/work/spooles-build.log"
#	cd $(abs_builddir)/work/internal/MT/src ; ${MAKE} global -f makeGlobalLib
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
#	cd $(abs_builddir)/work/internal/ ; ${MAKE} check
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	rm -rf ../i/SPOOLES/lib
	mkdir -p ../i/SPOOLES/lib
	cd $(abs_builddir)/work/internal/    ; cp spooles.a ../../../i/SPOOLES/lib/libspooles.a
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

distdir:

clean:
	rm -rf work
distclean:
	rm -rf work
