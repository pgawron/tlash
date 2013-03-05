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
# lntexinputs.sh
#
# Field G. Van Zee
#
# Creates symbolic links in the current directory of the files passed in as 
# arguments. This is because it seems nearly impossible to set an environment
# variable from within the top-level Makefile. (In this case, TEXINPUTS.)
#

print_usage()
{
	script_name="$0"
	
	echo "Usage: $script_name src_list"
	
	exit 1
}

main()
{
	local src_list
	
	
	# Check the number of arguments after command line option processing.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	
	# Extract arguments.
	src_list="$1"
	
	
	# Soft link each file into the current directory.
	for src_file in $src_list; do
		
		# Create a soft link of the $src_file in the current directory.
		ln -s $src_file .
	done
	
	
	# Exit peacefully.
	return 0
}


main "$@"
