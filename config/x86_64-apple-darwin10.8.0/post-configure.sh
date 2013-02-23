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

#echo "################################################################################"
echo ""
echo "libflame configuration summary"

echo ""
echo "revision........................................ : r`cat revision`"

echo ""
echo "build system type............................... : x86_64-apple-darwin10.8.0"
echo "host system type................................ : x86_64-apple-darwin10.8.0"

echo ""
echo "Enable verbose make output...................... : no"

echo ""
echo "Enable maximum argument list hack............... : no"

echo ""
echo "C compiler...................................... : gcc"
echo "Library archiver................................ : ar"
echo "Library archive indexer......................... : ranlib"

echo ""
echo "Enable Windows build (experimental)............. : no"

echo ""
echo "Create static library........................... : yes"
echo "Create shared (dynamically-linked) library...... : no"

echo ""
echo "Enable non-critical code........................ : yes"

echo ""
echo "Enable built-in BLAS implementation............. : no"

echo ""
echo "Enable lapack2flame............................. : no"
echo "Enable external LAPACK for subproblems.......... : no"
echo "Enable external LAPACK interfaces............... : no"

echo ""
echo "Enable multithreading support................... : no"
if   [ "none" = "openmp" ] ; then
echo "   Threading implementation..................... : OpenMP"
echo "      OpenMP C compiler/linker flags............ : "
elif [ "none" = "pthreads" ] ; then
echo "   Threading implementation..................... : POSIX threads"
fi

echo ""
echo "Enable SuperMatrix support...................... : no"

echo ""
echo "Enable GPU support.............................. : no"

echo ""
echo "Enable SCC support.............................. : no"

echo ""
echo "Enable support for Texas Instruments' DSP....... : no"

echo ""
echo "Enable vector intrinsics........................ : no"
if   [ "none" = "sse" ] ; then
echo "   Vector intrinsic type........................ : SSE"
echo "      SSE C compiler flags...................... : "
fi

echo ""
echo "Enable memory alignment......................... : no"
if [ "no" = "yes" ] ; then
echo "   User-requested memory alignment boundary..... : "
echo "   Enable leading dimesion alignment............ : no"
fi

echo ""
echo "C compiler language flags....................... : -std=c99"

echo ""
echo "Enable compiler optimizations................... : yes"
echo "   C compiler optimization flags................ : -O3"

echo "Enable compiler warnings........................ : yes"
echo "   C compiler warning flags..................... : -Wall -Wno-comment"

echo "Enable compiler debugging symbols............... : no"
echo "   C compiler debug flags....................... : -g0"

echo "Enable compiler profiling symbols............... : no"
echo "   C compiler profiling flags................... : "

echo ""
echo "Enable internal error checking.................. : yes"
if   [ "yes" = "yes" ] ; then
echo "   Internal error checking level................ : full"
fi

echo ""
echo "Enable memory leak counter...................... : no"

echo ""
echo "Enable level-3 BLAS front-end control trees..... : yes"

echo ""
echo "Enable BLIS use of FLA_malloc()................. : yes"

echo ""
echo "Enable interfaces to internal libgoto symbols... : no"

echo ""
echo "Enable interfaces to CBLAS...................... : no"

echo ""
echo "Enable user-defined default m blocksize......... : no"
if [ "no" = "yes" ] ; then
echo "   User-requested default m blocksize........... : "
fi
echo "Enable user-defined default k blocksize......... : no"
if [ "no" = "yes" ] ; then
echo "   User-requested default k blocksize........... : "
fi
echo "Enable user-defined default n blocksize......... : no"
if [ "no" = "yes" ] ; then
echo "   User-requested default n blocksize........... : "
fi

echo ""
echo "Enable portable timer........................... : yes"
if [ "yes" = "yes" ] ; then
echo "   Portable timer function...................... : gettimeofday()"
fi

echo ""
echo "Autodetect Fortran linker flags................. : @fla_enable_autodetect_f77_ldflags@"
if [ "@fla_enable_autodetect_f77_ldflags@" = "yes" ] ; then
echo "   Fortran linker flags......................... : @FLIBS@"
fi
echo "Autodetect Fortran name-mangling................ : @fla_enable_autodetect_f77_name_mangling@"
if [ "@fla_enable_autodetect_f77_name_mangling@" = "yes" ] ; then
echo "   Unmangled name............................... : @fla_f77_foobar_unmangled@"
echo "   Mangled name................................. : @fla_f77_foobar_mangled@"
fi

echo ""
echo "Compile with extra C compiler flags............. : no"
if [ "no" = "yes" ] ; then
echo "   User-requested extra C compiler flags........ : "
fi

echo ""
echo "libflame install directory prefix............... : /Users/martin/flame"
echo ""
echo "Configuration complete!"
echo ""
if [ "@fla_enable_autodetect_f77_ldflags@" = "yes" ] ; then
echo "NOTE: Autodetection of Fortran linker flags was enabled. The configure"
echo "script thinks that the flags listed above are necessary to successfully"
echo "link a program to Fortran object code. If your program uses any Fortran"
echo "libraries, you will probably need to link with these flags."
echo ""
fi
echo "You may now run 'make' to build all libflame libraries and then 'make install'"
echo "to install the libraries." 
echo ""
#echo "################################################################################"

