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
# rebuild-libflame-book.sh
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
	echo " Rebuild the libflame doxygen documentation."
	echo " Note:"
	echo "  o The toplevel_dir argument should be the absolute path to the"
	echo "    top-level directory where the libflame distribution is located."
	echo "  o The toplevel_dir directory will have \'svn update\' run on it"
	echo "    before the build is performed."
	echo " "
	echo " Usage:"
	echo "   ${script_name} toplevel_dir"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	# Make sure we have the svn binary.
	export PATH=$PATH:/lusr/bin
	
	# Check the number of arguments.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	# Extract our arguments.
	top_dirpath=$1

	# Change to the source directory.
	cd ${top_dirpath}
	
	# Run 'svn update' to make sure we have the latest source.
	svn update

	# Update the revision file in case it is out-of-date.
	./build/update-check-rev-file.sh
	
	# Update the revision number in the doxygen config file.
	revnum=$(cat revision)
	cat Doxyfile | sed "s/revision_anchor/$revnum/g" > Doxyfile.temp
	mv Doxyfile Doxyfile.orig
	mv Doxyfile.temp Doxyfile

	# Make the doxygen directory so we don't cause doxygen to state that
	# it's creating the directory for us.
	mkdir doxygen
	
	# Run doxygen.
	doxygen

	# Publish the new documentation.
	rm -rf /u/www/users/flame/libflame/docs/html
	mv doxygen/html /u/www/users/flame/libflame/docs
	
	# Restore the original file for next time.
	mv Doxyfile.orig Doxyfile

	# Remove the doxygen docs directory.
	rm -rf doxygen

	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
