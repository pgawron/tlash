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
AC_DEFUN([FLA_CHECK_ENABLE_DEFAULT_K_BLOCKSIZE],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested specific default blocksize in k dimension])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([default_k_blocksize],
	              [  --enable-default-k-blocksize=kb],
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" || test "$enableval" = "0" ; then
			
			dnl User provided --enable-<option>=no or --enable-<option>=0.
			AC_MSG_ERROR([[Invalid option to --enable-default-k-blocksize]])

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes.
			AC_MSG_ERROR([[Invalid option to --enable-default-k-blocksize]])
		else
			
			dnl User provided a valid argument (hopefully).
			fla_enable_default_k_blocksize=yes
			fla_default_k_blocksize=$enableval
		fi
		
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_default_k_blocksize=no
	]
	)

	dnl Now act according to whether the option was requested.
	if test "$fla_enable_default_k_blocksize" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Enable option by setting a corresponding preprocessor directive
		dnl to the requested value.
		AC_DEFINE_UNQUOTED(FLA_DEFAULT_K_BLOCKSIZE,$fla_default_k_blocksize,
		                   [Sets the default blocksize in the k dimension.])

		dnl Tell the user we're checking the value given.
		AC_MSG_CHECKING([user-requested k dimension blocksize])
		AC_MSG_RESULT([$fla_default_k_blocksize])
		
		dnl Substitute output variable values.
		AC_SUBST(fla_enable_default_k_blocksize)
		AC_SUBST(fla_default_k_blocksize)

	else

		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
])
