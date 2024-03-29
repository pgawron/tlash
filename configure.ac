dnl
dnl  libflame
dnl  An object-based infrastructure for developing high-performance
dnl  dense linear algebra libraries.
dnl
dnl  Copyright (C) 2011, The University of Texas
dnl
dnl  libflame is free software; you can redistribute it and/or modify
dnl  it under the terms of the GNU Lesser General Public License as
dnl  published by the Free Software Foundation; either version 2.1 of
dnl  the License, or (at your option) any later version.
dnl
dnl  libflame is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
dnl  Lesser General Public License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public
dnl  License along with libflame; if you did not receive a copy, see
dnl  http://www.gnu.org/licenses/.
dnl
dnl  For more information, please contact us at flame@cs.utexas.edu or
dnl  send mail to:
dnl
dnl  Field G. Van Zee and/or
dnl  Robert A. van de Geijn
dnl  The University of Texas at Austin
dnl  Department of Computer Sciences
dnl  1 University Station C0500
dnl  Austin TX 78712
dnl

dnl 
dnl configure.ac
dnl 
dnl Process this file with autoconf to produce a configure script.
dnl 
dnl Notice that this file has distinct sections. While the ordering of these
dnl sections (and by proxy, their constituent macros) is not strictly required,
dnl the order given here is encouraged by the autoconf manual published at 
dnl http://www.gnu.org/.
dnl
dnl Friendly reminder: All macro arguments must be quoted with `[' and `]'!
dnl Also, all literal strings should be double quoted!
dnl


dnl ----------------------------------------------------------------------------
dnl Initialization and configuration
dnl ----------------------------------------------------------------------------

dnl Initialize autoconf. The arguments here are autoconf macros defined to
dnl their appropriate values elsewhere in the build directory.
dnl AC_INIT([libflame],
dnl         [r]m4_esyscmd_s([cat revision]),
dnl         [flame@cs.utexas.edu],
dnl         [libflame],
dnl         [http://www.cs.utexas.edu/users/flame/])
AC_INIT([libflame],
        [r`cat revision`],
        [flame@cs.utexas.edu],
        [libflame],
        [http://www.cs.utexas.edu/users/flame/])

dnl Check for a good autoconf version. We need to do more research into the
dnl functionality provided by various versions of autoconf.
AC_PREREQ(2.15)

dnl Copyright notice for autoconf-generated configure script.
AC_COPYRIGHT([[Copyright 2011 The University of Texas
This configure script may be copied, distributed and modified under the terms
of the LGPL. See LICENSE for more details.]])

dnl Specify some file in the source directory to ensure so that the configure
dnl script can verify that the directory in which the file resides actually 
dnl does exist. This is a safety check.
AC_CONFIG_SRCDIR([src/base/flamec/include/FLAME.h])

dnl Specify the directory containing our auxiliary autoconf scripts. These
dnl include missing, mkinstalldirs, install-sh, config.guess, and config.sub.
AC_CONFIG_AUX_DIR([build/ac-utils])

dnl Specify the directory containing additional autoconf macros.
dnl (Not needed for now since we invoke aclocal manually with -I argument.)
dnl AC_CONFIG_MACRO_DIR([build/ac-macros])

dnl Specify the name of the global config header.
AC_CONFIG_HEADERS([build/FLA_config.h])

dnl Register commands to create directories for configure byproducts and
dnl object files.
AC_CONFIG_COMMANDS_PRE(
[
	dnl Update the revision file if we have a working copy.
	echo "configure: updating/checking revision file."
	build/update-check-rev-file.sh -v

	dnl Create the config/$host directory.
	echo "configure: creating config/$host"
	mkdir -p config/$host

	dnl Create the obj/$host directory.
	echo "configure: creating obj/$host"
	mkdir -p obj/$host

	dnl Create the lib/$host directory.
	echo "configure: creating lib/$host"
	mkdir -p lib/$host

	dnl Save the host string to a file so it can be retrieved from
	dnl the config.status commands defined below.
	echo "configure: creating config.sys_type"
	echo $host > config.sys_type
])

dnl Register commands to be executed by config.status after the main configure
dnl script is finished executing.
AC_CONFIG_COMMANDS([config.status],
[
	dnl Retreive the host_cpu value that was left by the AC_CONFIG_COMMANDS_PRE
	dnl commands above.
	echo "config.status: retreiving host string from config.sys_type"
	host=$(cat config.sys_type)

	dnl Move the config.status byproducts into the appropriate config
	dnl sub-directory for the current architecture.
	echo "config.status: moving config.mk to config/$host"
	mv build/config.mk         config/$host/
	echo "config.status: moving post-configure.sh to config/$host"
	mv build/post-configure.sh config/$host/
	echo "config.status: moving FLA_config.h to config/$host"
	mv build/FLA_config.h      config/$host/

	dnl Mirror the source tree directory in obj/$host.
	echo "config.status: mirroring source directory into obj/$host"
	build/mirror-tree.sh src obj/$host

	dnl Generate makefile fragments throughout the source tree.
	echo "config.status: running gen-make-frags.sh"
	build/gen-make-frags.sh

	dnl Before running the post-configure script, make sure its permissions
	dnl are set to be executable.
	echo "config.status: running post-configure.sh"
	chmod 0744 config/$host/post-configure.sh
	config/$host/post-configure.sh
])

dnl Determine the canonical build and host system types, which has the form
dnl `cpu-vendor-os' (or `cpu-vendor-kernel-os'). Ex: `i686-pc-linux-gnu'.
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

dnl Set the @fla_host_cpu@ output variable to something simple like p4 or it2
dnl or opteron based on the canonical host. Also, set the search paths for 
dnl AC_PROG_CC and AC_PROG_F77 marcos in the variables fla_c_compiler_list and
dnl fla_f77_compiler_list, respectively.
FLA_OBSERVE_HOST_CPU_TYPE

dnl Set the default value of the $prefix value, normally set by the --prefix
dnl option to configure.
AC_PREFIX_DEFAULT([$HOME/flame])

dnl Check for default set of C header files.
dnl AC_DEFAULT_INCLUDES




dnl ----------------------------------------------------------------------------
dnl Checks for programs.
dnl ----------------------------------------------------------------------------

dnl Check for GNU make.
FLA_REQUIRE_GNU_MAKE 

dnl Check for bash 2.0 or higher.
FLA_REQUIRE_GNU_BASH
dnl FLA_REQUIRE_GNU_BASH_VERSION(2.0)

dnl Check whether the user requested a specific C compiler. If so, add it to the
dnl front of the compiler search list. Otherwise do nothing.
FLA_CHECK_WITH_CC

dnl Determine the C compiler.
FLA_REQUIRE_CC

dnl Check whether the user requested a specific library archiver. If so, assign
dnl it to the requested archiver variable. Otherwise do nothing.
FLA_CHECK_WITH_AR

dnl Check for ar.
FLA_REQUIRE_AR

dnl Check for ranlib.
dnl AC_PROG_RANLIB
FLA_CHECK_WITH_RANLIB

dnl Require ranlib.
FLA_REQUIRE_RANLIB

dnl Check for install.
AC_PROG_INSTALL

dnl Check for find.
FLA_REQUIRE_FIND

dnl Check for xargs.
FLA_REQUIRE_XARGS

dnl Check for ln -s.
AC_PROG_LN_S



dnl ----------------------------------------------------------------------------
dnl Checks for libraries.
dnl ----------------------------------------------------------------------------

dnl Observe whether the user wants us to auto-detect linker flags that may be
dnl necessary for linking against Fortran object code.
dnl FLA_CHECK_ENABLE_AUTODETECT_F77_LDFLAGS

dnl Observe whether the user wants us to auto-detect the name-mangling strategy
dnl for caling Fortran routines from C (and C-compiled routines from Fortran).
dnl FLA_CHECK_ENABLE_AUTODETECT_F77_NAME_MANGLING

dnl NOTE: These conditionals may not be spun off into a separate macro because
dnl autoconf places dependent macros BEFORE the body of the macro that requires
dnl it, which means the dependent Fortran macros would get called regardless.
dnl NOTE: These conditionals are ordered very carefully! Don't muck with them
dnl unless you know what you're doing! That goes for you too, Field!
dnl if test "$fla_enable_autodetect_f77_ldflags"       = "yes" -o \
dnl         "$fla_enable_autodetect_f77_name_mangling" = "yes" ; then
dnl 
dnl	AC_F77_LIBRARY_LDFLAGS
dnl fi
dnl
dnl if test "$fla_enable_autodetect_f77_name_mangling" = "yes" ; then
dnl
dnl	FLA_OBSERVE_F77_NAME_MANGLING
dnl
dnl fi
dnl
dnl if test "$fla_enable_autodetect_f77_name_mangling" = "no"  -a \
dnl        "$fla_enable_autodetect_f77_ldflags"       = "no"  ; then

dnl	fla_f77_name_mangling="lower case, underscore, no extra underscore (static default)"
dnl	AC_SUBST(fla_f77_name_mangling)

dnl	FLIBS="# <Please-add-Fortran-77-linker-flags-here!>"
dnl	AC_SUBST(FLIBS)
dnl fi

dnl if test "$fla_enable_autodetect_f77_name_mangling" = "yes"  -a \
dnl        "$fla_enable_autodetect_f77_ldflags"       = "no"  ; then
dnl
dnl	AC_MSG_ERROR([User requested that we disable autodetection of Fortran linker flags but enable autodetection of Fortran name-mangling. However, we can't proceed since the former must first be detected in order to detect the latter.])
dnl
dnl fi



dnl ----------------------------------------------------------------------------
dnl Checks for header files.
dnl ----------------------------------------------------------------------------

AC_HEADER_STDC

dnl Itanium architectures seem to need this header, at least when compiling with
dnl Intel compilers.
AC_CHECK_HEADERS([ia64intrin.h])

dnl Other important headers.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h signal.h assert.h math.h])

dnl Check for the existence of the header file sys/time.h.
AC_CHECK_HEADERS([sys/time.h])

dnl If a program may include both time.h and sys/time.h, define the preprocessor
dnl macro TIME_WITH_SYS_TIME. On some older systems, sys/time.h includes time.h,
dnl but time.h is not protected against multiple inclusions, so programs should
dnl not explicitly include both files.
AC_HEADER_TIME




dnl ----------------------------------------------------------------------------
dnl Checks for typedefs and structures.
dnl ----------------------------------------------------------------------------





dnl ----------------------------------------------------------------------------
dnl Checks for standard compiler characteristics.
dnl ----------------------------------------------------------------------------

dnl Checks whether the C compiler can be run as a C preprocessor, which is
dnl usually invoked as "$CC -E".
dnl AC_PROG_CPP

dnl Check whether the C compiler found above accepts -c and -o options
dnl simultaneously.
AC_PROG_CC_C_O

dnl Check for type 'size_t'
AC_TYPE_SIZE_T

dnl Check for 'const' keyword support
AC_C_CONST

dnl Check for support for the inline keyword
AC_C_INLINE

dnl Check for C prototypes
AC_C_PROTOTYPES



dnl ----------------------------------------------------------------------------
dnl Checks for build-related options.
dnl ----------------------------------------------------------------------------

dnl Observe whether we're going to enable verbose output from make.
FLA_CHECK_ENABLE_VERBOSE_MAKE_OUTPUT

dnl Observe whether we're going to build a static library.
FLA_CHECK_ENABLE_STATIC_BUILD

dnl Observe whether we're going to build a dynamic library.
FLA_CHECK_ENABLE_DYNAMIC_BUILD

dnl Observe whether we should archive object files from a separate flat object
dnl directory so that systems with kernels with small ARG_MAX can still archive
dnl the library.
FLA_CHECK_ENABLE_MAX_ARG_LIST_HACK



dnl Make sure the user requested at least one type of build (static and/or
dnl dynamic).
FLA_REQUIRE_STATIC_OR_DYNAMIC_BUILD

dnl Check whether the user requested the ability to specify extra C compiler
dnl flags.
FLA_CHECK_WITH_EXTRA_CFLAGS

dnl Set any generally recommended C preprocessor flags for the chosen compiler.
FLA_SET_C_PREPROC_FLAGS

dnl Set any generally recommended C language flags for the chosen compiler.
FLA_SET_C_LANG_FLAGS



dnl ----------------------------------------------------------------------------
dnl Checks for various configuration options.
dnl ----------------------------------------------------------------------------

dnl Observe whether we're going to enable various non-critical code segments.
FLA_CHECK_ENABLE_NON_CRITICAL_CODE

dnl Observe whether we're going to enable a built-in implementation of the BLAS.
FLA_CHECK_ENABLE_BUILTIN_BLAS

dnl Observe whether we're going to build and include the lapack2flame
dnl compatibility layer.
FLA_CHECK_ENABLE_LAPACK2FLAME

dnl Observe whether we're going to enable external LAPACK for computing
dnl subproblems.
FLA_CHECK_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS

dnl Observe whether we should build code that provides FLAME-like interfaces to
dnl certain LAPACK routines.
FLA_CHECK_ENABLE_EXTERNAL_LAPACK_INTERFACES

dnl Make sure the user did not enable both lapack2flame and external LAPACK
dnl interfaces.
FLA_REQUIRE_NO_LAPACK_NAMESPACE_CONFLICT

dnl Observe how we should implement the level-3 BLAS front-ends.
FLA_CHECK_ENABLE_BLAS3_FRONT_END_CNTL_TREES

dnl Observe the threading switch.
FLA_CHECK_ENABLE_MULTITHREADING

dnl Observe the SuperMatrix switch.
FLA_CHECK_ENABLE_SUPERMATRIX

dnl Observe the GPU extension switch.
FLA_CHECK_ENABLE_GPU

dnl Make sure the user did not request a weird combination of parallelization
dnl options.
FLA_REQUIRE_SUPERMATRIX_ENABLED

dnl Observe the vector intrinsics switch and type.
FLA_CHECK_ENABLE_VECTOR_INTRINSICS

dnl Observe whether we're being given a specific memory alignment boundary.
FLA_CHECK_ENABLE_MEMORY_ALIGNMENT

dnl Observe whether we're being given a specific memory alignment boundary.
FLA_CHECK_ENABLE_LDIM_ALIGNMENT

dnl Set the appropriate compiler optimization flags.
FLA_CHECK_ENABLE_OPTIMIZATIONS

dnl Set the appropriate compiler warning flags, if requested.
FLA_CHECK_ENABLE_WARNINGS

dnl Set the appropriate compiler debugging flags, if requested.
FLA_CHECK_ENABLE_DEBUG

dnl Set the appropriate compiler profiling flags, if requested.
FLA_CHECK_ENABLE_PROFILING

dnl Observe whether we're going to perform internal error checking.
FLA_CHECK_ENABLE_INTERNAL_ERROR_CHECKING

dnl Observe whether we're going to enable the libflame memory leak detector.
FLA_CHECK_ENABLE_MEMORY_LEAK_COUNTER

dnl Observe whether we're going to enable usage of FLA_malloc() within
dnl bli_malloc().
FLA_CHECK_ENABLE_BLIS_USE_OF_FLA_MALLOC

dnl Observe whether we're going to enable interfaces to internal/low-level
dnl libgoto symbols.
FLA_CHECK_ENABLE_GOTO_INTERFACES

dnl Observe whether we're going to enable the CBLAS interface to the BLAS.
FLA_CHECK_ENABLE_CBLAS_INTERFACES

dnl Observe whether we're being given alternate default blocksizes.
FLA_CHECK_ENABLE_DEFAULT_M_BLOCKSIZE
FLA_CHECK_ENABLE_DEFAULT_K_BLOCKSIZE
FLA_CHECK_ENABLE_DEFAULT_N_BLOCKSIZE

dnl Check for a request to define and compile FLA_Clock() using the more
dnl portable gettimeofday() instead of using architecture-specific inline
dnl assembly.
FLA_CHECK_ENABLE_PORTABLE_TIMER

dnl Observe whether we're building for Windows. (Note this actually never gets
dnl invoked, because the Windows build doesn't use a configure script.)
FLA_CHECK_ENABLE_WINDOWS_BUILD

dnl Observe the SCC extension switch.
FLA_CHECK_ENABLE_SCC

dnl Observe the TI DSP switch.
FLA_CHECK_ENABLE_TIDSP

dnl Unconditionally define BLIS_FROM_LIBFLAME, which signals to BLIS that
dnl it is being used from within libflame.
dnl AC_DEFINE(BLIS_FROM_LIBFLAME,1,
dnl           [Unconditionally use libflame types (such as complex types) in BLIS.])

dnl Set a preprocessor constant so that the BLIS includes FLAME.h.
fla_blis_flags="-DBLIS_FROM_LIBFLAME"
AC_SUBST(fla_blis_flags)


dnl ----------------------------------------------------------------------------
dnl Checks for library functions.
dnl ----------------------------------------------------------------------------

dnl Check for error_on_line() function (used by #error directive?).
AC_FUNC_ERROR_AT_LINE

dnl Check for GNU-compatible malloc().
dnl AC_FUNC_MALLOC

dnl Check for working memcmp().
AC_FUNC_MEMCMP

dnl Check for functions that appear in FLAME source (?)
dnl AC_CHECK_FUNCS([fdatasync])

dnl Check for sqrt in the math library
AC_CHECK_LIB([m],[sqrt])
AC_CHECK_LIB([m],[fabs])
AC_CHECK_LIB([m],[pow])



dnl ----------------------------------------------------------------------------
dnl Output
dnl ----------------------------------------------------------------------------

dnl Specify files to be instantiated (by default from filename.in), substituting
dnl the output variable values for each file. This macro is used by AC_OUTPUT.
AC_CONFIG_FILES([build/config.mk build/post-configure.sh])

dnl Write the config.status file and launch it, which will create the individual
dnl Makefiles specified by AC_CONFIG_FILES.
AC_OUTPUT

