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
# svn-info-revision.sh
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
	echo " "${script_name}
	echo " "
	echo " Field G. Van Zee"
	echo " "
	echo " Query the subversion revision number for a file using 'svn info'."
	echo " This script is designed for use with doxygen."
	echo " "
	echo " Usage:"
	echo "   ${script_name} filepath"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	local filepath
	
	# The name of the script
	script_name=${0##*/}

	# Check the number of arguments.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	# Extract arguments.
	filepath="$1"
	
	# Check that old_dir actually exists and is a directory
	if [ ! -f ${filepath} ]; then
		echo "Hmm, ${filepath} does not seem to be a valid regular file."
		exit 1
	fi
	
	# Acquire the revision number for the file given.
	revision=$(svn info ${filepath} | grep 'Revision' | sed 's/Revision: //g')
	
	# Output the revision number.
	echo "r${revision}"

	# Exit peacefully.
	return 0
}

main "$@"
