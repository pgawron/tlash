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
AC_DEFUN([FLA_SET_C_PROF_FLAGS],
[
	AC_REQUIRE([FLA_OBSERVE_HOST_CPU_TYPE])
	
	AC_MSG_CHECKING([for (guessing) appropriate $CC profiling flags])

	if test "$1" == "yes" ; then
		
		dnl Set C compiler flags assuming we found...
		case $CC in
			dnl GNU gcc.
			gcc)
				fla_c_prof_flags='-pg'
			;;
			dnl Intel cc.
			icc)
				fla_c_prof_flags='-p'
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_prof_flags='-pg'
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_prof_flags='-pg'
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_prof_flags='-p'
			;;
			dnl IBM xlc.
			xlc)
				fla_c_prof_flags='-pg'
			;;
			dnl ambiguous cc.
			cc)
				fla_c_prof_flags='-pg'
			;;
			dnl for all other C compilers.
			*)
				fla_c_prof_flags=''
			;;
		esac
	else

		dnl Set C compiler flags assuming we found...
		case $CC in
			dnl GNU gcc.
			gcc)
				fla_c_prof_flags=''
			;;
			dnl Intel cc.
			icc)
				fla_c_prof_flags=''
			;;
			dnl PathScale pathcc.
			pathcc)
				fla_c_prof_flags=''
			;;
			dnl PGI pgcc.
			pgcc)
				fla_c_prof_flags=''
			;;
			dnl NEC sxcc.
			sxcc)
				fla_c_prof_flags='-Np'
			;;
			dnl IBM xlc.
			xlc)
				fla_c_prof_flags=''
			;;
			dnl ambiguous cc.
			cc)
				fla_c_prof_flags=''
			;;
			dnl for all other C compilers.
			*)
				fla_c_prof_flags=''
			;;
		esac
	fi
	
	dnl Output the result.
	AC_MSG_RESULT([$fla_c_prof_flags])
	
	dnl Substitute the prof flags into the autoconf output files
	AC_SUBST(fla_c_prof_flags)

])
