--- Make.inc.orig	1999-01-19 06:48:55.000000000 +0900
+++ Make.inc	2026-08-05 00:00:00.000000000 +0900
@@ -12,7 +12,15 @@
 # for solaris
 #
 # CC = gcc
-  CC = /usr/lang-4.0/bin/cc
+# Every CC assignment in this file is left commented out ON PURPOSE.  The sdpa build
+# passes CC on the SPOOLES make line (see external/spooles/Makefile.am, or spooles/Makefile
+# in the qd fork), and a command-line assignment beats anything written here, so the
+# compiler ./configure selected is the one that builds SPOOLES.  If you drive this makefile
+# YOURSELF you must supply it too --  make global -f makefile CC=gcc  -- because the
+# ".POSIX:" on line 1 makes GNU make default CC to c99, not cc, and c99 rejects every -f
+# option, so the build dies on the -funroll-all-loops below.  That failure is the one
+# install.sh's spooles_rescue() exists to paper over.
+# CC = cc
 #
 # for sgi
 #
@@ -35,7 +43,23 @@
 # CFLAGS = -Wall -g
 # CFLAGS = -Wall -pg
 # CFLAGS = $(OPTLEVEL) -D_POSIX_C_SOURCE=199506L
-  CFLAGS = $(OPTLEVEL)
+# SPOOLES is 1999 C.  GCC 14 promoted -Wint-conversion, -Wimplicit-function-declaration and
+# -Wincompatible-pointer-types from warnings to errors, and SPOOLES trips the first of them
+# (ETree/src/transform.c:294, IVinit(nfront, NULL)).  Measured on gcc-16: without the three
+# -Wno-error= flags the build stops with a truncated 674 KB spooles.a; with them it completes
+# at 2545752 bytes.  -D_GNU_SOURCE: some toolchains hide struct timezone even under gcc.
+# These belong here, in the shipped patch, not only in CI YAML and a skill script.
+#
+# `override` (GNU make) is deliberate.  Jobs that instrument the tree pass CFLAGS on the make
+# line -- CI builds SPOOLES with -ftrivial-auto-var-init and with ASan -- and a plain `+=`
+# loses to a command-line assignment, which would silently take these flags away again and
+# put the GCC-14 errors back.  With `override` the instrumentation flags and these are BOTH
+# applied.  The trade is that CFLAGS on the make line can no longer remove -O2; if you need
+# that, edit this line.
+  override CFLAGS += -O2 -funroll-all-loops $(pthread_cflags) -D_GNU_SOURCE \
+                     -Wno-error=int-conversion \
+                     -Wno-error=implicit-function-declaration \
+                     -Wno-error=incompatible-pointer-types
 # CFLAGS = -Wall $(OPTLEVEL)
 #
 #----------------------------------------
@@ -79,8 +103,8 @@
 #  set the ranlib environment
 #  (if ranlib is not needed, we echo the library name)
 #
-# RANLIB = ranlib
-  RANLIB = echo
+  RANLIB = ranlib
+# RANLIB = echo
 #
 #---------------------------------------------------------------------
 #
