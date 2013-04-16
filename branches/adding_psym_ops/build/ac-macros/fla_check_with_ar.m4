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
AC_DEFUN([FLA_CHECK_WITH_AR],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a specific library archiver])
	
	dnl Determine whether the user gave the --enable-<option> or
    dnl --disable-<option>. If so, then run the first snippet of code;
    dnl otherwise, run the second code block.
	AC_ARG_WITH([ar],
	            AC_HELP_STRING([--with-ar=ar],[ Search for and use a library archiver named <ar>. If <ar> is not found, then use the first library archiver found from the default search list for the detected build architecture. Note: the library archiver search list usually consists only of "ar".]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_specific_ar=without
			fla_requested_ar=''

		elif test "$withval" = "yes" ; then
			
			dnl User provided --with-<option>=yes or --with-<option>.
			fla_with_specific_ar=no
			fla_requested_ar=''
		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_specific_ar=yes
			fla_requested_ar=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
        dnl Default behavior is to disable the option.
		fla_with_specific_ar=no
		fla_requested_ar=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_specific_ar" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_requested_ar)])
		
	elif test "$fla_with_specific_ar" = "without" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl The user requested --with-ar=no or --without-ar. Scold him.
		AC_MSG_ERROR([[Detected --without-ar. Cannot continue with library archiver disabled!]])
		
	else
		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
	
	dnl Check for AR environment variable, which would override everything else.
	if test "$AR" != "" ; then
		AC_MSG_NOTICE([[AR environment variable is set to $AR, which will override --with-ar option and default search list for library archiver.]])
	fi

])