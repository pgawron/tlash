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
AC_DEFUN([FLA_SET_C_OPENMP_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
		
	dnl Echo C OpenMP flags to user.
	AC_MSG_CHECKING([for (guessing) OpenMP flags for $CC])
		
	dnl Set the OpenMP compiler flags based on which C compiler we're going to use.
	case $CC in
		dnl Intel icc.
		icc)
			fla_c_openmp_flags='-openmp'
		;;
		dnl GNU gcc.
		gcc)
			fla_c_openmp_flags='-fopenmp'
		;;
		dnl PathScale pathcc
		pathcc)
			fla_c_openmp_flags='-mp'
		;;
		dnl PGI pgcc.
		pgcc)
			fla_c_openmp_flags='-mp'
		;;
		dnl NEC sxcc.
		sxcc)
			fla_c_openmp_flags='-P openmp'
		;;
		dnl IBM xlc
		xlc)
			fla_c_openmp_flags='-qsmp=omp'
		;;
		dnl for all other C compilers.
		*)
			fla_c_openmp_flags='unknown'
		;;
	esac
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_c_openmp_flags])
	
	dnl Check string in case C OpenMP compiler is unknown.
	if test "$fla_c_openmp_flags" = "unknown" ; then
		
		dnl Tell the user we can't continue unless we know what flags
		dnl to pass to the C compiler to enable OpenMP.
		AC_MSG_ERROR([configure doesn't know what flag to give $CC in order to enable OpenMP support. Please submit a bug report to the FLAME developers at FLA_BUG_REPORT_ADDRESS.])
	fi

	dnl Output the C OpenMP flags variable.
	AC_SUBST(fla_c_openmp_flags)
])
