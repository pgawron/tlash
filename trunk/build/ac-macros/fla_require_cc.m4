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
AC_DEFUN([FLA_REQUIRE_CC],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])

	dnl Save the value of CFLAGS. This will come in useful later in determining
	dnl whether the user provided his own definition of CFLAGS.
	fla_userdef_cflags=$CFLAGS
	
	dnl Find a C compiler. 
	dnl If the CC environment variable is not already set, search for the
	dnl compiler defined by fla_requested_cc (which may be empty) and then
	dnl continue searching for the compilers in $fla_c_compiler_list, if
	dnl necessary. Also, if the C compiler found is not in ANSI mode, then
	dnl try to add an option to make it so. If the GNU gcc was found, then
	dnl GCC shell variable is set to `yes'. 
	AC_PROG_CC([$fla_requested_cc $fla_c_compiler_list])

	if test "$CC" = "" ; then
		AC_MSG_ERROR([Could not locate any of the following C compilers: $CC $fla_requested_cc $fla_c_compiler_list. Cannot continue without a C compiler.],[1])
	fi
	
	dnl Substitute the user-defined CFLAGS into the autoconf output files.
	AC_SUBST(fla_userdef_cflags)
])
