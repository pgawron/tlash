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
AC_DEFUN([FLA_REQUIRE_GNU_BASH_VERSION],
[
	dnl Require that GNU bash already be present.
	AC_REQUIRE([FLA_REQUIRE_GNU_BASH])

	dnlfla_tmp_file="fla_require_gnu_bash_version_$1.tmp";
	fla_tmp_file=fla_bash.tmp;

	dnl Output a message saying we're going to check that the version of
	dnl GNU bash is equal to or later than the argument.
	AC_MSG_CHECKING([[for GNU bash >= $1]])
		
	dnlbash --version | grep 'version' | sed -e 's;.*version.;;' | sed -e 's;\..*$;;' 2> /dev/null 1> $fla_tmp_file
	bash --version | grep bash | dd skip=18 bs=1 count=3 2> /dev/null 1> $fla_tmp_file
	dnlecho -n "3.0" > $fla_tmp_file
	sh -c "echo ' >= '$1" >> $fla_tmp_file
	sh -c "echo 'halt'" >> $fla_tmp_file
	dnl dnlbroken fla_bash_version=`bash --version | grep 'ver' | sed -e 's;.*version.\([0-9\.][0-9\.]*\).*;\1;'`
	fla_result=`bc -q $fla_tmp_file`

	if test $fla_result = "1" ; then
		_cv_gnu_bash_version=yes;
		AC_MSG_RESULT([[yes]]);
	else
		AC_MSG_RESULT([[no]]);
		AC_MSG_ERROR([[Could not find GNU bash version >= $1! Bailing out...]],[1])
	fi
	rm -f $fla_tmp_file 

])
