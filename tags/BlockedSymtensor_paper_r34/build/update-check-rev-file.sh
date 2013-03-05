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
# update-check-rev-file.sh
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
	echo " Create a file that contains the revision number for the current"
	echo " working copy. The revision number is obtained from the subversion"
	echo " file located in .svn/entries of the top-level directory of the"
	echo " source tree. It is assumed that this script is being run at the"
	echo " top-level directory of a checked-out working copy, since an"
	echo " exported tarball will already have the revision file created. If"
	echo " the entries file does not exist, it will verify that the revision"
	echo " file exists. If it does not, it is created with a dummy revision"
	echo " string. Otherwise, if it does exist, the script will leave the"
	echo " the file untouched."
	echo " "
	echo " This script is run by configure to ensure that the revision file"
	echo " exists before it is needed."
	echo " "
	echo " Usage:"
	echo "   ${script_name}"
	echo " "
	echo " The following arguments are accepted/required:"
	echo " "
	echo "   -v          verbose"
	echo "                 Be verbose. Output what's happening."
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	# Assume we are being run at the top-level directory.
	toplevel_dirpath=.
	svn_admin_dirpath=${toplevel_dirpath}/.svn
	revision_filename=revision
	revision_filepath=${toplevel_dirpath}/${revision_filename}
	dummy_rev_string=unknown

	# Variables set by getopts.
	verbose_flag=""

	# Process our command line options.
	while getopts ":v" opt; do
		case $opt in
			v  ) verbose_flag=1 ;;
			\? ) print_usage
		esac
	done
	shift $(($OPTIND - 1))
	
	# Check the number of arguments after command line option processing.
	if [ $# != "0" ]; then
		print_usage
	fi
	
	# Extract our arguments.
	script_name=${0##*/}

	# If we are in a working copy, we can overwrite the revision file with a
	# potentially new value.
	if [ -d "${svn_admin_dirpath}" ]; then

		# Acquire the revision number for the file given.
		rev_num=$(svn info . | grep 'Revision' | sed 's/Revision: //g')

		# Be verbose if requested.
		if [ -n "${verbose_flag}" ]; then
			echo "${script_name}: Found working copy; writing revision string \"${rev_num}\" to ${revision_filepath}"
		fi

		# Echo the revision number to the revision file.
		echo "${rev_num}" > ${revision_filepath}

 	# If we can't find the svn admin directory, we probably are in an exported
	# copy: either an official snapshot, or a copy that someone exported
	# manually--hopefully (and likely) the former.
	else

		# Be verbose if requested.
		if [ -n "${verbose_flag}" ]; then
			echo "${script_name}: Found export. Checking for revision file..."
		fi

		# Make sure a revision file exists. If it does not, then create
		# a dummy file so configure has something to work with.
		if [ ! -f "${revision_filepath}" ]; then

			# Be verbose if requested.
			if [ -n "${verbose_flag}" ]; then
				echo "${script_name}: Revision file not found. Writing file ${revision_filepath} with dummy revision string \"${dummy_rev_string}\"."
			fi

			# Echo the revision number to the revision file.
			echo "${dummy_rev_string}" > ${revision_filepath}

		else

			# Get the revision number from the file just for the purposes of
			# being verbose, if it was requested.
			rev_num=$(cat ${revision_filepath})

			# Be verbose if requested.
			if [ -n "${verbose_flag}" ]; then
				echo "${script_name}: Revision file found containing revision string \"${rev_num}\". Snapshot is valid export."
			fi


		fi
	fi
	
	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
