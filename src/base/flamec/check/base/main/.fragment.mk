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
CURRENT_DIR_NAME := main
CURRENT_SUB_DIRS := 

# Source files local to this fragment
LOCAL_SRC_FILES  := FLA_Axpy_buffer_to_object_check.c FLA_Axpy_object_to_buffer_check.c FLA_Cont_with_1x3_to_1x2_check.c FLA_Cont_with_3x1_to_2x1_check.c FLA_Cont_with_3x3_to_2x2_check.c FLA_Copy_buffer_to_object_check.c FLA_Copy_object_to_buffer_check.c FLA_Merge_1x2_check.c FLA_Merge_2x1_check.c FLA_Merge_2x2_check.c FLA_Obj_attach_buffer_check.c FLA_Obj_buffer_at_view_check.c FLA_Obj_copy_view_check.c FLA_Obj_create_buffer_check.c FLA_Obj_create_complex_constant_check.c FLA_Obj_create_conf_to_check.c FLA_Obj_create_constant_check.c FLA_Obj_create_constant_ext_check.c FLA_Obj_create_ext_check.c FLA_Obj_create_without_buffer_check.c FLA_Obj_datatype_check.c FLA_Obj_datatype_proj_to_real_check.c FLA_Obj_datatype_size_check.c FLA_Obj_elem_size_check.c FLA_Obj_elemtype_check.c FLA_Obj_equals_check.c FLA_Obj_extract_complex_scalar_check.c FLA_Obj_extract_imag_part_check.c FLA_Obj_extract_real_part_check.c FLA_Obj_extract_real_scalar_check.c FLA_Obj_free_buffer_check.c FLA_Obj_free_check.c FLA_Obj_free_without_buffer_check.c FLA_Obj_fshow_check.c FLA_Obj_ge_check.c FLA_Obj_gt_check.c FLA_Obj_le_check.c FLA_Obj_lt_check.c FLA_Obj_set_imag_part_check.c FLA_Obj_set_real_part_check.c FLA_Obj_show_check.c FLA_Part_1x2_check.c FLA_Part_2x1_check.c FLA_Part_2x2_check.c FLA_Repart_1x2_to_1x3_check.c FLA_Repart_2x1_to_3x1_check.c FLA_Repart_2x2_to_3x3_check.c FLA_Submatrix_at_check.c

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
