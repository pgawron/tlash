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
AC_DEFUN([FLA_CHECK_ENABLE_BLIS_USE_OF_FLA_MALLOC],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested enabling code that defines bli_malloc() in terms of FLA_malloc()])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([blis-use-of-fla-malloc],
	              AC_HELP_STRING([--enable-blis-use-of-fla-malloc],[Enable code that defines bli_malloc() in terms of FLA_malloc(). One benefit of this is that BLIS memory allocations can be tracked, along with other libflame memory allocations, if the memory leak counter is enabled. A second benefit is that BLIS memory allocations can be aligned to boundaries if libflame memory alignment is enabled. Note this option may only be set at configure-time. (Enabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_blis_use_of_fla_malloc=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_blis_use_of_fla_malloc=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BLIS_USE_OF_FLA_MALLOC!]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to enable the option.
		fla_enable_blis_use_of_fla_malloc=yes
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_blis_use_of_fla_malloc" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_BLIS_USE_OF_FLA_MALLOC,1,
		          [Determines whether to define bli_malloc() in terms of FLA_malloc().])
		
	elif test "$fla_enable_blis_use_of_fla_malloc" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_BLIS_USE_OF_FLA_MALLOC!]])
	fi
	
	dnl Substitute the output variable.
	AC_SUBST(fla_enable_blis_use_of_fla_malloc)

])
