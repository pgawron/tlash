#!/bin/bash
#
#  libflame
#  An object-based infrastructure for developing high-performance
#  dense linear algebra libraries.
#
#  Copyright (C) 2011, The University of Texas
#
#  libflame is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation; either version 2.1 of
#  the License, or (at your option) any later version.
#
#  libflame is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with libflame; if you did not receive a copy, see
#  http://www.gnu.org/licenses/.
#
#  For more information, please contact us at flame@cs.utexas.edu or
#  send mail to:
#
#  Field G. Van Zee and/or
#  Robert A. van de Geijn
#  The University of Texas at Austin
#  Department of Computer Sciences
#  1 University Station C0500
#  Austin TX 78712
#

#
# query-revision-file.sh
#
# Field G. Van Zee
#

print_usage()
{
	local script_name
	
	# Get the script name
	script_name=${0##*/}
	
	# Echo usage info
	echo " "
	echo " "$script_name
	echo " "
	echo " Field G. Van Zee"
	echo " "
	echo " Query the revision string for the HEAD of the repository."
	echo " "
	echo " Usage:"
	echo "   ${script_name}"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	# Assume we are being run at the top-level directory.
	host=flame.cs.utexas.edu
	repos_dirpath=/u/field/repos/flame
	test_filepath=trunk/configure

	# Check the number of arguments after command line option processing.
	if [ $# != "0" ]; then
		print_usage
	fi
	
	# Construct the subversion URL to the FLAME repository.
	svn_url_revtest="svn+ssh://${host}/${repos_dirpath}/${test_filepath}"

	# Acquire the revision number for the test file of the HEAD.
	rev_num_str=$(svn info ${svn_url_revtest} | grep 'Revision' \
	                                          | sed 's/Revision: //g')

	# Output the revision number.
	echo "${rev_num_str}"

	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
