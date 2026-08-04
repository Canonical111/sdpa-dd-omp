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
case "$host_cpu" in
  aarch64|arm64) QD_FMA_ARG="--enable-fma=gnu" ;;
  *)             QD_FMA_ARG="" ;;
esac
AC_SUBST([QD_FMA_ARG])

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
CXXFLAGS="$CXXFLAGS $pthread_cflags"
AC_SUBST(pthread_lib)

CFLAGS="$CFLAGS -O2 -funroll-all-loops"
CXXFLAGS="$CXXFLAGS -O2 -funroll-all-loops"
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)

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
