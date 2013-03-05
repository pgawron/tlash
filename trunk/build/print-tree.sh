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
# print-tree.sh
#
# Field G. Van Zee
#
# Usage:
#   print-tree.sh [-d] dir1 [dir2 ...]
#
# Recursively descends through the directory arguments given and prints to
# standard output all items found (files and directories). Each item printed
# is indented with an additional '\t' tab character according to its depth
# in the directory hierarchy. 
# 


main()
{
	# Process our command line options. We only respond to the -d flag,
	# in which case we'll only echo directories as we go along.
	while getopts ":d" opt; do
		case $opt in
			d  ) dirs_only_flag="1" ;;
			\? ) echo 'Usage: print-tree.sh [-d] dir1 [dir2 ...]'
				exit 1
		esac
	done
	shift $(($OPTIND - 1))
	
	
	# Define our tab character
	one_tab="\t"
	
	
	# Process each command line argument given
	for top_thing in "$@"; do
		
		# Echo the item
		echo $top_thing
		
		
		# If top_thing is a directory, then descend recursively into it.
		if [ -d "$top_thing" ]; then
			this_thing=$top_thing
			print_tree $(command ls $top_thing)
		fi
	done
	
	
	# Exit peacefully.
	return 0
}

print_tree()
{
	# Add a tab character to our tab string
	tab_str=$tab_str$one_tab
	
	
	# Process each item in our argument list (ie: each item in this_thing)
	for thing in "$@"; do
		
		# Append the current item to the current directory
		this_thing=$this_thing/$thing
		
		
		# Echo the tabbed item. If the dirs_only_flag was given,
		# and $thing is a directory, then output. Otherwise, do not
		# output. If the dirs_only_flag was not given, then output
		# regardless of file type.
		if [ -n "$dirs_only_flag" ]; then
			if [ -d $this_thing ]; then
				echo -e "$tab_str$thing"
			fi
		else
			echo -e "$tab_str$thing"
		fi
		
		
		# If this_thing is a directory, then descent recursively into it.
		if [ -d "$this_thing" ]; then
			print_tree $(command ls $this_thing)
		fi
		
		
		# Delete the end of the path, up to the first / character to
		# prepare for the next "file" in $@
		this_thing=${this_thing%/*}
	done
	
	
	# Remove one tab from the end of our tab_str string, given that
	# we're about to return up the function stack one level.
	tab_str=${tab_str%"\t"}
}

# The script's main entry point, passing all parameters given.
main "$@"

