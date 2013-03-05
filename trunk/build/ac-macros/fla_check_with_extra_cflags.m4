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
AC_DEFUN([FLA_CHECK_WITH_EXTRA_CFLAGS],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested using extra C compiler flags])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_WITH([extra-cflags],
	            AC_HELP_STRING([--with-extra-cflags=flagstring],[When compiling C code, use the flags in flagstring in addition to the flags that configure would normally use. This is useful when the user wants some extra flags passed to the compiler but does not want to manually set the CFLAGS environment variable and thus override all of the default compiler flags.]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$withval" = "no" ; then
			
			dnl User provided --with-<option>=no or --without-<option>.
			fla_with_extra_cflags=no
			fla_extra_cflags=''

		else
			
			dnl User provided argument value: --with-<option>=value.
			fla_with_extra_cflags=yes
			fla_extra_cflags=$withval
		fi
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_with_extra_cflags=no
		fla_extra_cflags=''
	]
	)
	
	dnl Now act according to whether a specific value was given.
	if test "$fla_with_extra_cflags" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
	else

		dnl Output the result.
		AC_MSG_RESULT([yes ($fla_extra_cflags)])
		
	fi
	
	dnl Check for CC environment variable, which would override everything else.
	if test "$EXTRA_CFLAGS" != "" ; then
		AC_MSG_NOTICE([[EXTRA_CFLAGS environment variable is set to $EXTRA_CFLAGS, which will override --with-extra-cflags option.]])
        fla_with_extra_cflags=yes
        fla_extra_cflags="$EXTRA_CFLAGS"
	fi

	dnl Declare EXTRA_CFLAGS as precious
	AC_ARG_VAR([EXTRA_CFLAGS],[extra C compiler flags to be used in addition to those determined automatically by configure])
	
	dnl Substitute the output variable.
	AC_SUBST(fla_with_extra_cflags)
	AC_SUBST(fla_extra_cflags)
])
