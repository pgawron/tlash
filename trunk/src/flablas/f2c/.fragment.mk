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
# fragment.mk 
#
# This is an automatically-generated makefile fragment and will likely get
# overwritten or deleted if the user is not careful. Modify at your own risk.
#

# These two mmakefile variables need to be set in order for the recursive
# include process to work!
CURRENT_DIR_NAME := f2c
CURRENT_SUB_DIRS := static util

# Source files local to this fragment
LOCAL_SRC_FILES  := caxpy.c ccopy.c cdotc.c cdotu.c cgbmv.c cgemm.c cgemv.c cgerc.c cgeru.c chbmv.c chemm.c chemv.c cher.c cher2.c cher2k.c cherk.c chpmv.c chpr.c chpr2.c crotg.c cscal.c csrot.c csscal.c cswap.c csymm.c csyr2k.c csyrk.c ctbmv.c ctbsv.c ctpmv.c ctpsv.c ctrmm.c ctrmv.c ctrsm.c ctrsv.c dasum.c daxpy.c dcabs1.c dcopy.c ddot.c dgbmv.c dgemm.c dgemv.c dger.c dnrm2.c drot.c drotg.c drotm.c drotmg.c dsbmv.c dscal.c dsdot.c dspmv.c dspr.c dspr2.c dswap.c dsymm.c dsymv.c dsyr.c dsyr2.c dsyr2k.c dsyrk.c dtbmv.c dtbsv.c dtpmv.c dtpsv.c dtrmm.c dtrmv.c dtrsm.c dtrsv.c dzasum.c dznrm2.c icamax.c idamax.c isamax.c izamax.c lsame.c sasum.c saxpy.c scasum.c scnrm2.c scopy.c sdot.c sdsdot.c sgbmv.c sgemm.c sgemv.c sger.c snrm2.c srot.c srotg.c srotm.c srotmg.c ssbmv.c sscal.c sspmv.c sspr.c sspr2.c sswap.c ssymm.c ssymv.c ssyr.c ssyr2.c ssyr2k.c ssyrk.c stbmv.c stbsv.c stpmv.c stpsv.c strmm.c strmv.c strsm.c strsv.c zaxpy.c zcopy.c zdotc.c zdotu.c zdrot.c zdscal.c zgbmv.c zgemm.c zgemv.c zgerc.c zgeru.c zhbmv.c zhemm.c zhemv.c zher.c zher2.c zher2k.c zherk.c zhpmv.c zhpr.c zhpr2.c zrotg.c zscal.c zswap.c zsymm.c zsyr2k.c zsyrk.c ztbmv.c ztbsv.c ztpmv.c ztpsv.c ztrmm.c ztrmv.c ztrsm.c ztrsv.c

# Add the fragment's local source files to the _global_variable_ variable.
MK_FLABLAS_F2C_SRC += $(addprefix $(PARENT_PATH)/$(CURRENT_DIR_NAME)/, $(LOCAL_SRC_FILES))




# -----------------------------------------------------------------------------
# NOTE: The code below is generic and should remain in all fragment.mk files!
# -----------------------------------------------------------------------------

# Add the current fragment to the global list of fragments so the top-level
# Makefile knows which directories are participating in the build.
FRAGMENT_DIR_PATHS  += $(PARENT_PATH)/$(CURRENT_DIR_NAME)

# Recursively descend into other subfragments' local makefiles and include them.
ifneq ($(strip $(CURRENT_SUB_DIRS)),)
key                 := $(key).x
stack_$(key)        := $(PARENT_PATH)
PARENT_PATH         := $(PARENT_PATH)/$(CURRENT_DIR_NAME)
FRAGMENT_SUB_DIRS   := $(addprefix $(PARENT_PATH)/, $(CURRENT_SUB_DIRS))
-include  $(addsuffix /$(FRAGMENT_MK), $(FRAGMENT_SUB_DIRS))
PARENT_PATH         := $(stack_$(key))
key                 := $(basename $(key))
endif
