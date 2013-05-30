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
	echo " Rebuild the libflame book from LaTeX source."
	echo " "
	echo " Note:"
	echo "  o The src_dir argument should be the absolute path to the"
	echo "    top-level libflame directory."
	echo "  o The docs/libflame subdirectory of src_dir will have \'svn update\'"
	echo "    run on it before the build is performed."
	echo " "
	echo " Usage:"
	echo "   ${script_name} src_dir"
	echo " "
	
	# Exit with non-zero exit status
	exit 1
}


main()
{
	# Make sure we have the svn binary.
	export PATH=$PATH:/lusr/bin
	
	# Make sure we have the pdflatex binary.
	export PATH=$PATH:/lusr/tex/bin
	
	# Check the number of arguments.
	if [ $# != "1" ]; then
		print_usage
	fi
	
	# Extract our arguments.
	toplevel_dirpath=$1
	libflame_dirpath=docs/libflame

	# Change to the top-level directory.
	cd ${toplevel_dirpath}
	
	# Run 'svn update' on docs directory to make sure we have the latest
	# source.
	svn update docs
	svn update -N

	# Update the revision file in case it is out-of-date.
	#./build/svn-info-revision.sh Makefile > revision
	./build/update-check-rev-file.sh

	# Change into libflame directory.
	cd ${libflame_dirpath}

	# Clean out the directory and then rebuild the document via make.
	make clean
	make 2> make.err
	
	# Exit peacefully.
	return 0
}


# The script's main entry point, passing all parameters given.
main "$@"
