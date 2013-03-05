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
CURRENT_DIR_NAME := util
CURRENT_SUB_DIRS := 

# Source files local to this fragment
LOCAL_SRC_FILES  := FLA_Absolute_square_check.c FLA_Absolute_value_check.c FLA_Add_to_diag_check.c FLA_Apply_GTG_check.c FLA_Apply_G_1x2_check.c FLA_Apply_G_check.c FLA_Apply_G_mx2_check.c FLA_Conjugate_check.c FLA_Conjugate_r_check.c FLA_Fill_with_cluster_dist_check.c FLA_Fill_with_geometric_dist_check.c FLA_Fill_with_inverse_dist_check.c FLA_Fill_with_linear_dist_check.c FLA_Fill_with_logarithmic_dist_check.c FLA_Fill_with_random_dist_check.c FLA_Form_perm_matrix_check.c FLA_Givens1_check.c FLA_Givens2_check.c FLA_Hermitianize_check.c FLA_Househ2_UT_check.c FLA_Househ2s_UT_check.c FLA_Househ3UD_UT_check.c FLA_Introduce_bulge_check.c FLA_Inv_scal_elemwise_check.c FLA_Invert_check.c FLA_LU_find_zero_on_diagonal_check.c FLA_Max_abs_value_check.c FLA_Max_abs_value_herm_check.c FLA_Max_elemwise_diff_check.c FLA_Mult_add_check.c FLA_Negate_check.c FLA_Norm1_check.c FLA_Norm_frob_check.c FLA_Norm_inf_check.c FLA_Pow.c FLA_Random_herm_matrix_check.c FLA_Random_matrix_check.c FLA_Random_spd_matrix_check.c FLA_Random_symm_matrix_check.c FLA_Random_tri_matrix_check.c FLA_Random_unitary_matrix_check.c FLA_Scal_elemwise_check.c FLA_Scale_diag_check.c FLA_Set_check.c FLA_Set_diag_check.c FLA_Set_to_identity_check.c FLA_Setr_check.c FLA_Shift_diag_check.c FLA_Shift_pivots_to_check.c FLA_Sort_check.c FLA_Sort_evd_check.c FLA_Sort_svd_check.c FLA_Sqrt_check.c FLA_Symmetrize_check.c FLA_Transpose_check.c FLA_Triangularize_check.c FLA_Wilkshift_tridiag_check.c

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
