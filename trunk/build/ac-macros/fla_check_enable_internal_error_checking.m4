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
AC_DEFUN([FLA_CHECK_ENABLE_INTERNAL_ERROR_CHECKING],
[
	dnl Initialize some variables.
	fla_enable_internal_error_checking=no
	fla_internal_error_checking_level=none
	
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested support for internal error checking])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([internal-error-checking],
	              AC_HELP_STRING([--enable-internal-error-checking=level],[Enable various internal runtime checks of function parameters and object properties to prevent functions from executing with unexpected values. Note that this option determines the default level, which may be changed at runtime. Valid values for level are "full", "minimal", and "none". (Enabled by default to "full".)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "full" ; then
			
			dnl Enable with OpenMP support.
			fla_enable_internal_error_checking=yes
			fla_internal_error_checking_level=full

		elif test "$enableval" = "minimal" ; then
			
			dnl Enable with POSIX threads support.
			fla_enable_internal_error_checking=yes
			fla_internal_error_checking_level=minimal

		elif test "$enableval" = "none" ; then
			
			dnl Disable internal error checking.
			fla_enable_internal_error_checking=no
			fla_internal_error_checking_level=none

		else
			
			dnl Invalid option.
			AC_MSG_ERROR([[Invalid option to --enable-internal-error-checking. Valid options are "full", "minimal", and "none".]])
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_internal_error_checking=yes
		fla_internal_error_checking_level=full
	]
	)
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_enable_internal_error_checking])
		
	dnl Now act according to whether the option was requested.
	if test "$fla_enable_internal_error_checking" = "yes" ; then
		
		dnl Tell the user we're checking the value given.
		AC_MSG_CHECKING([user-requested internal error checking level])
		AC_MSG_RESULT([$fla_internal_error_checking_level])

		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_INTERNAL_ERROR_CHECKING,1,
		          [Determines whether to enable internal runtime consistency checks of function parameters and return values.])

		dnl Now we set cpp macros related to each internal error checking level.
		if test "$fla_internal_error_checking_level" = "full" ; then
		
			dnl Encode the C prepropcessor value for full internal error checking.
			fla_level_val=2
			
		elif test "$fla_internal_error_checking_level" = "minimal" ; then

			dnl Encode the C prepropcessor value for minimal internal error checking.
			fla_level_val=1

		fi

	else

		dnl Encode the C prepropcessor value for no internal error checking.
		fla_level_val=0

	fi

	dnl Define the preprocessor macro INTERNAL_ERROR_CHECKING_LEVEL to the value corresponding
	dnl to full, minimal, or no internal error checking, depending on how fla_level_val
	dnl was set above.
	AC_DEFINE_UNQUOTED(FLA_INTERNAL_ERROR_CHECKING_LEVEL,$fla_level_val,
	                   [Encodes the default level of internal error checking chosen at configure-time.])

	dnl Substitute the output variables.
	AC_SUBST(fla_enable_internal_error_checking)
	AC_SUBST(fla_internal_error_checking_level)
])
