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
# Makefile
#
# Field G. Van Zee
# 
# Makefile for libflame test suite driver.
#



#
# --- Include libflame config makefile fragment --------------------------------
#

# Determine the path to the libflame config makefile fragment. We'll use
# several variables defined there.
BUILD_DIRPATH       := ../build
CONFIG_DIRPATH      := ../config
HOST                := $(shell sh $(BUILD_DIRPATH)/ac-utils/config.guess)
CONFIG_MK_FRAGMENT  := $(CONFIG_DIRPATH)/$(HOST)/config.mk

# Include the definitions in the config makefile fragment.
-include $(CONFIG_MK_FRAGMENT)



#
# --- Optional overrides -------------------------------------------------------
#

# Uncomment and modify these definitions if you wish to override the values
# present in the master config makefile fragment.
 CC             := gcc
# LINKER         := $(CC)
 CFLAGS         := -ggdb -O0 -Wall -Wno-comment
 LDFLAGS        := 
 INSTALL_PREFIX := $(HOME)/flame



#
# --- BLAS and LAPACK implementations ------------------------------------------
#

# BLAS implementation path. A BLAS library must be given in order to run
# the libflame test suite. Modify these definitions if needed.
#LIBBLAS        := $(LIBBLAS_PATH)/libblas.a
#LIBBLAS        := $(LIBBLAS_PATH)/libgoto.a

#My Mac
LIBBLAS_PATH   := $(INSTALL_PREFIX)/lib
LIBBLAS        := -framework accelerate

#Clarksville
#LIBBLAS_PATH    := $(HOME)/OpenBLAS/lib
#LIBBLAS         := $(LIBBLAS_PATH)/libopenblas_single.a

# LAPACK implementation path. These values only matter if libflame was
# configured with the external-lapack-interfaces option enabled. Modify
# these definitions if needed.
#LIBLAPACK      := $(LIBLAPACK_PATH)/liblapack.a

#My Mac
LIBLAPACK_PATH := $(INSTALL_PREFIX)/lib
LIBLAPACK      := -framework accelerate

#
# --- General build definitions ------------------------------------------------
#

TEST_SRC_PATH  := src
TEST_OBJ_PATH  := obj

FLA_LIB_PATH   := $(INSTALL_PREFIX)/lib
FLA_INC_PATH   := $(INSTALL_PREFIX)/include
LIBFLAME       := $(FLA_LIB_PATH)/libflame.a
LIBFLAME_PROF  := $(HOME)/flame_prof/lib/libflame.a

CFLAGS         += -I$(FLA_INC_PATH) -I$(TEST_SRC_PATH)

FNAME          := libflame

TEST_OBJS      := $(patsubst $(TEST_SRC_PATH)/%.c, \
                             $(TEST_OBJ_PATH)/%.o, \
                             $(wildcard $(TEST_SRC_PATH)/*.c))
TEST_BIN       := test_$(FNAME).x

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.c
		$(CC) $(CFLAGS) -c $< -o $@

test_$(FNAME): $(TEST_OBJS)
		$(LINKER) $(TEST_OBJS) $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o $(TEST_BIN)

sttsm: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/test_sttsm.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o test_sttsm

sttsm_but_one: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/test_sttsm_but_one.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o test_sttsm_but_one
		
tensor_print: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/test_tensor_print.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_print

tensor_part: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/test_tensor_part.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_part

tensor_permute: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/tensor_permute.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_permute

tensor_ttm: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/tensor_ttm.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_ttm

tensor_sym: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/tensor_sym.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_sym

tensor_sym_view: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/tensor_sym_view.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_sym_view

tensor_psym: $(TEST_OBJS) $(LIBFLAME)
		$(LINKER) $(TEST_OBJ_PATH)/tensor_psym.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_psym

tensor_psttm: $(TEST_OBJS) $(LIBFLAME)
		$(LINKER) $(TEST_OBJ_PATH)/tensor_psttm.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o tensor_psttm

sttsm_dense: $(TEST_OBJS)
		$(LINKER) $(TEST_OBJ_PATH)/test_sttsm_dense.o $(LDFLAGS) $(LIBFLAME) $(LIBLAPACK) $(LIBBLAS) -o test_sttsm_dense

sttsm_prof: $(TEST_OBJS)
		$(LINKER) -pg $(TEST_OBJ_PATH)/test_sttsm.o $(LDFLAGS) $(LIBFLAME_PROF) $(LIBLAPACK) $(LIBBLAS) -o test_sttsm_prof

sttsm_dense_prof: $(TEST_OBJS)
		$(LINKER) -pg $(TEST_OBJ_PATH)/test_sttsm_dense.o $(LDFLAGS) $(LIBFLAME_PROF) $(LIBLAPACK) $(LIBBLAS) -o test_sttsm_dense_prof

clean:
		$(RM_F) $(TEST_OBJS) $(TEST_BIN)

all: sttsm sttsm_but_one tensor_print tensor_part tensor_permute tensor_ttm tensor_sym tensor_sym_view tensor_psym tensor_psttm sttsm_dense
