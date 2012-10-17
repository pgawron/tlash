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
CURRENT_DIR_NAME := base
CURRENT_SUB_DIRS := transpose

# Source files local to this fragment
LOCAL_SRC_FILES  := FLA_Absolute_square.c FLA_Absolute_value.c FLA_Add_to_diag.c FLA_Adjust_2D_Info.c FLA_Binomial.c FLA_Clock.c FLA_Compare_Pairwise_Sort.c FLA_Conjugate.c FLA_Conjugate_r.c FLA_Fill_with_cluster_dist.c FLA_Fill_with_geometric_dist.c FLA_Fill_with_inverse_dist.c FLA_Fill_with_linear_dist.c FLA_Fill_with_logarithmic_dist.c FLA_Fill_with_random_dist.c FLA_Hermitianize.c FLA_Inv_scal_elemwise.c FLA_Invert.c FLA_Max_abs_value.c FLA_Max_abs_value_herm.c FLA_Max_elemwise_diff.c FLA_Mult_add.c FLA_Negate.c FLA_Norm1.c FLA_Norm_frob.c FLA_Norm_inf.c FLA_Pow.c FLA_Random_herm_matrix.c FLA_Random_matrix.c FLA_Random_spd_matrix.c FLA_Random_symm_matrix.c FLA_Random_symm_tensor.c FLA_Random_tri_matrix.c FLA_Random_unitary_matrix.c FLA_Scal_elemwise.c FLA_Scale_diag.c FLA_Set.c FLA_Set_diag.c FLA_Set_offdiag.c FLA_Set_to_identity.c FLA_Setr.c FLA_Shift_diag.c FLA_Sort.c FLA_Sqrt.c FLA_Symmetrize.c FLA_Triangularize.c FLA_random_number.c

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
