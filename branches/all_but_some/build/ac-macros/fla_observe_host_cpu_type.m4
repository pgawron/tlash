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
AC_DEFUN([FLA_OBSERVE_HOST_CPU_TYPE],
[
	AC_REQUIRE([AC_CANONICAL_HOST])
	
	case $host in
		dnl Intel Pentium-based class of processors, as well as those processors
		dnl such as AMD Athlons, Durons, etc that fall under the i?86 family.
		i386*-*-* | i586*-*-* | i686*-*-*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="gcc icc cc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl Intel EM64T or AMD Opteron/Athlon64 processors.
		x86_64*-*-*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="gcc icc pathcc cc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl Intel Itanium processors.
		ia64*-*-*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="icc gcc cc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl NEC SX systems.
		sx*-nec-superux*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="sxcc"
			fla_mpicc_compiler_list="mpicc"
			fla_ar=sxar
			fla_ranlib=ranlib
		;;
		dnl IBM POWER/AIX systems.
		powerpc*-ibm-aix*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="xlc"
			fla_mpicc_compiler_list="mpcc mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
		dnl PowerPC/Cell systems.
		powerpc64-*-linux-gnu)
			if test "$fla_enable_cell_spu_parallelism" = "yes" ; then
				fla_host_cpu=$host_cpu
				fla_c_compiler_list="ppu-gcc"
				fla_mpicc_compiler_list="mpicc"
				fla_ar=ppu-ar
				fla_ranlib=ppu-ranlib
			else
				fla_host_cpu=$host_cpu
				fla_c_compiler_list="gcc xlc"
				fla_mpicc_compiler_list="mpicc"
				fla_ar=ar
				fla_ranlib=ranlib
			fi
		;;
		dnl For all other proessors, use a basic search path.
		*)
			fla_host_cpu=$host_cpu
			fla_c_compiler_list="gcc cc" 
			fla_mpicc_compiler_list="mpicc"
			fla_ar=ar
			fla_ranlib=ranlib
		;;
	esac
	
	dnl Substitute the cpu type into the autoconf output files
	AC_SUBST(fla_host_cpu)

])
