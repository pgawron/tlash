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
AC_DEFUN([FLA_REQUIRE_NO_LAPACK_NAMESPACE_CONFLICT],
[
	AC_REQUIRE([FLA_CHECK_ENABLE_LAPACK2FLAME])
	AC_REQUIRE([FLA_CHECK_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS])
	AC_REQUIRE([FLA_CHECK_ENABLE_EXTERNAL_LAPACK_INTERFACES])

	dnl Make sure the user did not request an invalid/dangerous combination of
	dnl LAPACK compatibility/interfacing options.

	dnl Scenarios 1 and 2
	if test "$fla_enable_lapack2flame" = "yes" ; then
		if test "$fla_enable_external_lapack_for_subproblems" = "yes" ; then
			AC_MSG_ERROR([Configuring libflame to enable both lapack2flame and external-lapack-for-subproblems is not allowed. lapack2flame requires that external-lapack-for-subproblems be disabled. Please adjust your configure options accordingly and then re-run configure.],[1])
		fi
	fi

	dnl Scenario 3
	if test "$fla_enable_lapack2flame" = "no" ; then
		if test "$fla_enable_external_lapack_for_subproblems" = "yes" ; then
			if test "$fla_enable_external_lapack_interfaces" = "no" ; then
				AC_MSG_ERROR([Configuring libflame to enable external-lapack-for-subproblems without external-lapack-interfaces is not allowed. external-lapack-for-subproblems requires that external-lapack-interfaces be enabled. Please adjust your configure options accordingly and then re-run configure.],[1])
			fi
		fi
	fi
])
