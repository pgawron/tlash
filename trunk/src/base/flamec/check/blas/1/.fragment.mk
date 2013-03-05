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
CURRENT_DIR_NAME := 1
CURRENT_SUB_DIRS := 

# Source files local to this fragment
LOCAL_SRC_FILES  := FLA_Amax_check.c FLA_Asum_check.c FLA_Axpy_check.c FLA_Axpy_internal_check.c FLA_Axpyrt_check.c FLA_Axpys_check.c FLA_Axpyt_check.c FLA_Axpyt_internal_check.c FLA_Copy_check.c FLA_Copy_internal_check.c FLA_Copyr_check.c FLA_Copyr_internal_check.c FLA_Copyrt_check.c FLA_Copyt_check.c FLA_Copyt_internal_check.c FLA_Dot2cs_check.c FLA_Dot2s_check.c FLA_Dot_check.c FLA_Dotc_check.c FLA_Dotcs_check.c FLA_Dots_check.c FLA_Inv_scal_check.c FLA_Inv_scalc_check.c FLA_Nrm2_check.c FLA_Scal_check.c FLA_Scal_internal_check.c FLA_Scalc_check.c FLA_Scalr_check.c FLA_Scalr_internal_check.c FLA_Swap_check.c FLA_Swapt_check.c

# Add the fragment's local source files to the _global_variable_ variable.
MK_BASE_FLAMEC_SRC += $(addprefix $(PARENT_PATH)/$(CURRENT_DIR_NAME)/, $(LOCAL_SRC_FILES))




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
