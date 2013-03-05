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
AC_DEFUN([FLA_REQUIRE_GNU_BASH],
[
	dnl Output a message saying we're going to check for GNU bash.
	AC_MSG_CHECKING([[for GNU bash]])

	dnl Check the cache for the result. If it's not result is not yet cached,
	dnl execute the shell code to set the cache-id (ie: look for bash).
	AC_CACHE_VAL([_cv_gnu_bash_command],
	[
		dnl Initialize the cache-id to null.
		_cv_gnu_bash_command='';
	
		dnl Check that some version of bash is present.
		dnl If bash is not present, then output an error message.
		if ( sh -c "bash --version" 2> /dev/null | grep GNU 2>&1 > /dev/null ); then
			_cv_gnu_bash_command='bash';
		fi
	])
	
	dnl Now that we've checked for GNU bash, let's determine the result of the
	dnl macro and report and error if bash was not found.
	if test "$_cv_gnu_bash_command" = "bash" ; then
		AC_MSG_RESULT([[bash]])
	else
		AC_MSG_RESULT([[not found!]])
		AC_MSG_ERROR([[Could not locate GNU bash! Bailing out...]],[1])
	fi
])
