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
AC_DEFUN([FLA_REQUIRE_GNU_MAKE],
[
	AC_CACHE_CHECK( for GNU make,_cv_gnu_make_command,
		_cv_gnu_make_command='';
		dnl Search all the common names for GNU make
		for a in "$MAKE" make gmake gnumake ; do
			if test -z "$a" ; then continue ; fi ;
			if ( sh -c "$a --version" 2> /dev/null | grep GNU  2>&1 > /dev/null ) ;  then
				_cv_gnu_make_command=$a;
				break;
			fi
		done;
	);
	dnl If there was a GNU version, then set @fla_gnu_make_found@ to "true"
	if test  "x$_cv_gnu_make_command" != "x"  ; then
		fla_gnu_make_found=yes
	else
		fla_gnu_make_found=no
		AC_MSG_RESULT("not found!");
		AC_MSG_ERROR([[Could not locate GNU make! Bailing out...]],[1])
	fi
	AC_SUBST(fla_gnu_make_found)
])
