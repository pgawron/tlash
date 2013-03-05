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
AC_DEFUN([FLA_CHECK_ENABLE_MEMORY_ALIGNMENT],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested support for memory alignment])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([memory-alignment],
	              AC_HELP_STRING([--enable-memory-alignment=N],[Enable code that aligns dynamically allocated memory regions at N-byte boundaries. Note: N must be a power of two and multiple of sizeof(void*), which is usually 4 on 32-bit architectures and 8 on 64-bit architectures. (Disabled by default.)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "yes" || test "$enableval" = "0" ; then
			
			dnl Disallow disabling.
			AC_MSG_ERROR([[Invalid option to --enable-memory-alignment]])

		elif test "$enableval" = "no" ; then
			
			dnl Disable.
			fla_enable_memory_alignment=no

		else
			
			dnl User provided a valid argument (hopefully). Enable.
			fla_enable_memory_alignment=yes
			fla_memory_alignment_boundary=$enableval
		fi
		
	],
	[
		dnl User did not specify whether to enable or disable the option.
		dnl Default behavior is to disable the option.
		fla_enable_memory_alignment=no
	]
	)

	dnl Now act according to whether the option was requested.
	if test "$fla_enable_memory_alignment" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Verify that we have posix_memalign() on hand.
		AC_CHECK_FUNC([posix_memalign])
		
		dnl Define the preprocessor macro to enable the option.
		AC_DEFINE(FLA_ENABLE_MEMORY_ALIGNMENT,1,
		          [Determines whether memory is aligned to user-requested boundaries.])

		dnl Define an additional preprocessor directive that specifies the
		dnl requested value.
		AC_DEFINE_UNQUOTED(FLA_MEMORY_ALIGNMENT_BOUNDARY,$fla_memory_alignment_boundary,
		                   [Sets the byte boundary used to align the starting address of all memory allocated dynamically through libflame.])

		dnl Tell the user we're checking the value given.
		AC_MSG_CHECKING([user-requested memory alignment boundary])
		AC_MSG_RESULT([$fla_memory_alignment_boundary])
		
		dnl Substitute output variable values.
		AC_SUBST(fla_enable_memory_alignment)
		AC_SUBST(fla_memory_alignment_boundary)

	else

		dnl Output the result.
		AC_MSG_RESULT([no])
	fi
])

