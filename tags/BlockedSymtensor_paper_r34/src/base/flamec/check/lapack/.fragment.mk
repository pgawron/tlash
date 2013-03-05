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
CURRENT_DIR_NAME := lapack
CURRENT_SUB_DIRS := util

# Source files local to this fragment
LOCAL_SRC_FILES  := FLA_Apply_CAQ2_UT_internal_check.c FLA_Apply_CAQ_UT_inc_check.c FLA_Apply_CAQ_UT_inc_internal_check.c FLA_Apply_Q2_UT_check.c FLA_Apply_Q2_UT_internal_check.c FLA_Apply_QUD_UT_check.c FLA_Apply_QUD_UT_inc_check.c FLA_Apply_QUD_UT_inc_internal_check.c FLA_Apply_QUD_UT_internal_check.c FLA_Apply_Q_UT_check.c FLA_Apply_Q_UT_inc_check.c FLA_Apply_Q_UT_inc_internal_check.c FLA_Apply_Q_UT_internal_check.c FLA_Apply_Q_check.c FLA_Apply_diag_matrix_check.c FLA_Apply_pivots_check.c FLA_Bidiag_UT_check.c FLA_Bidiag_UT_extract_diagonals_check.c FLA_Bidiag_UT_form_U_check.c FLA_Bidiag_UT_form_V_check.c FLA_Bidiag_UT_internal_check.c FLA_Bidiag_UT_realify_check.c FLA_Bidiag_UT_recover_tau_check.c FLA_Bidiag_check.c FLA_Bidiag_form_U_check.c FLA_Bidiag_form_V_check.c FLA_CAQR2_UT_internal_check.c FLA_CAQR_UT_inc_check.c FLA_CAQR_UT_inc_solve_check.c FLA_Chol_check.c FLA_Chol_internal_check.c FLA_Chol_solve_check.c FLA_Eig_gest_check.c FLA_Eig_gest_internal_check.c FLA_FS_incpiv_check.c FLA_Hess_UT_check.c FLA_Hess_UT_internal_check.c FLA_Hess_UT_recover_tau_check.c FLA_Hess_check.c FLA_Hevd_check.c FLA_Hevd_compute_scaling_check.c FLA_Hevdd_check.c FLA_Hevdr_check.c FLA_LQ_UT_check.c FLA_LQ_UT_form_Q_check.c FLA_LQ_UT_internal_check.c FLA_LQ_UT_recover_tau_check.c FLA_LQ_UT_solve_check.c FLA_LQ_check.c FLA_LU_incpiv_check.c FLA_LU_incpiv_solve_check.c FLA_LU_nopiv_check.c FLA_LU_nopiv_internal_check.c FLA_LU_nopiv_solve_check.c FLA_LU_piv_check.c FLA_LU_piv_solve_check.c FLA_Lyap_check.c FLA_Lyap_internal_check.c FLA_QR2_UT_check.c FLA_QR2_UT_internal_check.c FLA_QR_UT_check.c FLA_QR_UT_copy_internal_check.c FLA_QR_UT_form_Q_check.c FLA_QR_UT_inc_check.c FLA_QR_UT_inc_solve_check.c FLA_QR_UT_internal_check.c FLA_QR_UT_recover_tau_check.c FLA_QR_UT_solve_check.c FLA_QR_check.c FLA_QR_form_Q_check.c FLA_SPDinv_check.c FLA_SPDinv_internal_check.c FLA_Svd_check.c FLA_Svd_compute_scaling_check.c FLA_Svdd_check.c FLA_Sylv_check.c FLA_Sylv_internal_check.c FLA_Tridiag_UT_check.c FLA_Tridiag_UT_extract_diagonals_check.c FLA_Tridiag_UT_form_Q_check.c FLA_Tridiag_UT_internal_check.c FLA_Tridiag_UT_realify_check.c FLA_Tridiag_UT_recover_tau_check.c FLA_Tridiag_UT_shift_U_check.c FLA_Tridiag_apply_Q_check.c FLA_Tridiag_check.c FLA_Tridiag_form_Q_check.c FLA_Trinv_check.c FLA_Trinv_internal_check.c FLA_Ttmm_check.c FLA_Ttmm_internal_check.c FLA_UDdate_UT_check.c FLA_UDdate_UT_inc_check.c FLA_UDdate_UT_inc_solve_check.c FLA_UDdate_UT_inc_update_rhs_check.c FLA_UDdate_UT_internal_check.c FLA_UDdate_UT_solve_check.c FLA_UDdate_UT_update_rhs_check.c

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
