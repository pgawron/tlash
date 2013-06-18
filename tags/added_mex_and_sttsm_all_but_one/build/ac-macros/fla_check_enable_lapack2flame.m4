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
AC_DEFUN([FLA_CHECK_ENABLE_LAPACK2FLAME],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested the lapack2flame compatibility layer])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([lapack2flame],
	              AC_HELP_STRING([--enable-lapack2flame],[Compile and build into libflame a compatibility layer that maps LAPACK invocations to their corresponding FLAME/C implementations. Note that erroneous input parameters are reported according to libflame conventions, NOT LAPACK conventions. That is, if libflame error checking is disabled, no error checking is performed, and if erroneous input parameters are detected, the library aborts. Also, if this option is enabled, then external-lapack-for-subproblems MUST be disabled. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_lapack2flame=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_lapack2flame=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_LAPACK2FLAME!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_lapack2flame=no
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_lapack2flame" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_LAPACK2FLAME,1,
		          [Determines whether the LAPACK compatibility layer is included in libflame.])
		
	elif test "$fla_enable_lapack2flame" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_LAPACK2FLAME!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_lapack2flame)

])
