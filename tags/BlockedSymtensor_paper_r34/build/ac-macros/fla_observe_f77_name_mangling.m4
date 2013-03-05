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
AC_DEFUN([FLA_OBSERVE_F77_NAME_MANGLING],
[
	dnl This is important. We have to invoke the AC_F77_WRAPPERS macro first so
	dnl it will invoke the name-mangling macro.
	AC_REQUIRE([AC_F77_WRAPPERS])
	
	fla_f77_foobar_unmangled='foobar'
	fla_f77_foobar_mangled=''

	dnl Ask the system to mangle a test name for us using the mangling that
	dnl will be used by the current Fortran environment.
	AC_F77_FUNC($fla_f77_foobar_unmangled,fla_f77_foobar_mangled)

	dnl Substitute the output variables.
	AC_SUBST(fla_f77_foobar_unmangled)
	AC_SUBST(fla_f77_foobar_mangled)
])
