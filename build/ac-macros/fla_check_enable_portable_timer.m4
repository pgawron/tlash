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
AC_DEFUN([FLA_CHECK_ENABLE_PORTABLE_TIMER],
[
	dnl Tell the user we're checking whether to enable the option.
	AC_MSG_CHECKING([whether user requested a portable FLA_Clock() timer])
	
	dnl Determine whether the user gave the --enable-<option> or
	dnl --disable-<option>. If so, then run the first snippet of code;
	dnl otherwise, run the second code block.
	AC_ARG_ENABLE([portable-timer],
	              AC_HELP_STRING([--enable-portable-timer],[Define the FLA_Clock() timer function using clock_gettime(). If that function is not available, then getttimeofday() is used. If neither function is available, FLA_Clock() is will return a static value. (By default, a portable timer is used (if it exists).)]),
	[
		dnl If any form of the option is given, handle each case.
		if test "$enableval" = "no" ; then
			
			dnl User provided --enable-<option>=no or --disable-<option>.
			fla_enable_portable_timer=no

		elif test "$enableval" = "yes" ; then
			
			dnl User provided --enable-<option>=yes or --enable-<option>.
			fla_enable_portable_timer=yes
		else

			dnl We don't need an else branch because the configure script
			dnl should detect whether the user provided an unexpected argument
			dnl with the option.
			AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_PORTABLE_TIMER!]])
		fi
	],
	[
			dnl User did not specify whether to enable or disable the portable
			dnl implementation of FLA_Clock(). Default behavior is to enable
			dnl the portable timer.
			fla_enable_portable_timer=yes
	]
	)
	
	dnl Now act according to whether the option was requested.
	if   test "$fla_enable_portable_timer" = "yes" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([yes])
		
		dnl Define the macro.
		AC_DEFINE(FLA_ENABLE_PORTABLE_TIMER,1,
		          [Determines whether to define a portable FLA_Clock() in terms of clock_gettime() or gettimeofday() from time.h.])
		
		dnl Look for clock_gettime() and gettimeofday().
		AC_CHECK_FUNC([clock_gettime])
		AC_CHECK_FUNC([gettimeofday])

		fla_portable_timer_function=''

		dnl Define the right cpp macro depending on which function was
		dnl present in the environment.
		if   test "$ac_cv_func_clock_gettime" = "yes" ; then

			fla_portable_timer_function='clock_gettime()'
			AC_DEFINE(FLA_PORTABLE_TIMER_IS_CLOCK_GETTIME,1,
			          [Determines whether clock_gettime() was present on the system (via time.h).])

		elif test "$ac_cv_func_clock_gettime" = "no" ; then

			if   test "$ac_cv_func_gettimeofday" = "yes" ; then

				fla_portable_timer_function='gettimeofday()'
				AC_DEFINE(FLA_PORTABLE_TIMER_IS_GETTIMEOFDAY,1,
				          [Determines whether gettimeofday() was present on the system (via time.h).])

			elif test "$ac_cv_func_gettimeofday" = "no" ; then

				AC_MSG_ERROR([[Neither clock_gettime() nor gettimeofday() were found! FLA_Clock() will be broken!]])

				fla_portable_timer_function='none found!'
				AC_DEFINE(FLA_PORTABLE_TIMER_IS_UNKNOWN,1,
				          [Determines whether a timer was found at all.])
			fi
		fi
		
	elif test "$fla_enable_portable_timer" = "no" ; then
		
		dnl Output the result.
		AC_MSG_RESULT([no])
		
		dnl We'll need strchr(), which is used by detect_clocks() in
		dnl FLA_Clock.c. Verify that we have it.
		AC_CHECK_FUNC([strchr],[],
		[
			AC_MSG_ERROR([Failed to find a working version of strchr()! Try enabling the portable timer, which does not need strchr(), and rerun configure.])
		])
		
	else
		dnl Only "yes" and "no" are accepted, so this block is empty.
		AC_MSG_ERROR([[Reached unreachable branch in FLA_CHECK_ENABLE_PORTABLE_TIMER!]])
	fi
	
	dnl Substitute the output variables.
	AC_SUBST(fla_enable_portable_timer)
	AC_SUBST(fla_portable_timer_function)

])
