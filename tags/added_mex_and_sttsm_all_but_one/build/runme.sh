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
# runme.sh
# Field G. Van Zee
#

print_usage()
{
	echo "Usage: $0 quoted_dir_list" 
	exit 1
}

main()
{
	local var_dirs
	local input_files
	local output_filepath
	local input_suffix
	local driver
	local output_prefix
	
	# Check the number of arguments
	if [ $# != "1" ]; then
		print_usage
	fi

	# Extract arguments.
	var_dirs="$1"
	
	# Establish some constants.
	driver=$(ls test*.x)
	output_prefix="output"

	
	for var_dir in ${var_dirs}; do
		
		input_files=$(ls ${var_dir}/input*)
		
		for input_filepath in ${input_files}; do
			
			input_filename=${input_filepath##*/}
			input_suffix=${input_filename##input_}
			output_filepath="${var_dir}/${output_prefix}_${input_suffix}.m"
			
			echo ">>> ${driver} < ${input_filepath} > ${output_filepath}"
			./${driver} < ${input_filepath} > ${output_filepath}
		
		done
		
	done
	
	# Exit peacefully.
	return 0
}

main "$@"
