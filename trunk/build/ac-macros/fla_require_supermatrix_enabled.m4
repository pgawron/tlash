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
AC_DEFUN([FLA_REQUIRE_SUPERMATRIX_ENABLED],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_SUPERMATRIX])
	AC_REQUIRE([FLA_CHECK_ENABLE_MULTITHREADING])
	AC_REQUIRE([FLA_CHECK_ENABLE_GPU])

	dnl Make sure the user did not request an invalid combination of
	dnl parallelization options.

	dnl Scenario 1: User wants GPU support but forgot to enable SM.
	if test "$fla_enable_gpu" = "yes" ; then
		if test "$fla_enable_supermatrix" = "no" ; then
			AC_MSG_ERROR([Configuring libflame to enable GPU support without also enabling SuperMatrix is not allowed. GPU utilization requires that SuperMatrix be enabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi
])
