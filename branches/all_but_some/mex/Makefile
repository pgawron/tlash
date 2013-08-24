CC = gcc
MEX = mex

FLAME_TOP_DIR = /Users/martin/flame
MATLAB_TOP_DIR = /Applications/MATLAB_R2013a.app

BLAS_CFLAGS = 
BLAS_LDFLAGS = -framework accelerate

#-----------------------------------------------------------------
FLAME_INC_DIR = $(FLAME_TOP_DIR)/include
FLAME_LIB_DIR = $(FLAME_TOP_DIR)/lib
FLAME_LIB = flame

MEX_INC_DIR = $(MATLAB_TOP_DIR)/extern/include

SRC_TOP_DIR = $(PWD)/src
SRC_BASE_MAIN_DIR = $(SRC_TOP_DIR)/base/flamec/tensor/main
SRC_BASE_INC_DIR = $(SRC_TOP_DIR)/base/flamec/tensor/include
SRC_BASE_UTIL_BASE_DIR = $(SRC_TOP_DIR)/base/flamec/tensor/util/base
SRC_BTAS_TOP_DIR = $(SRC_TOP_DIR)/btas
SRC_BTAS_TTM_DIR = $(SRC_BTAS_TOP_DIR)/ttm/flamec
SRC_BTAS_STTSM_DIR = $(SRC_BTAS_TOP_DIR)/sttsm/flamec
SRC_MATLAB_DIR = $(SRC_TOP_DIR)/matlab

MEX_CONVERT_PRIMITIVES_FILE = $(SRC_BASE_UTIL_BASE_DIR)/TLA_MEX_Convert_primitives.c
MEX_OBJ_FILE = $(SRC_BASE_MAIN_DIR)/TLA_MEX_Obj.c
TTM_FILE = $(SRC_BTAS_TTM_DIR)/TLA_ttm.c
STTSM_FILE = $(SRC_BTAS_STTSM_DIR)/TLA_sttsm.c

OBJ_TOP_DIR = $(PWD)/obj
OBJ_MEX_OBJ_FILE = $(OBJ_TOP_DIR)/TLA_MEX_Obj.o
OBJ_MEX_CONVERT_PRIMITIVES_FILE = $(OBJ_TOP_DIR)/TLA_MEX_Convert_primitives.o  
OBJS = $(OBJ_MEX_OBJ_FILE) $(OBJ_MEX_CONVERT_PRIMITIVES_FILE)

LDFLAGS = -L$(FLAME_LIB_DIR) -l$(FLAME_LIB) $(BLAS_LDFLAGS)
CFLAGS = -g -I$(FLAME_INC_DIR) -I$(MEX_INC_DIR) -I$(SRC_BASE_INC_DIR) $(BLAS_CFLAGS)

MEX_CFLAGS = CFLAGS='$$CFLAGS $(CFLAGS)'
MEX_LDFLAGS = LDFLAGS='$$LDFLAGS $(LDFLAGS)'
MEX_FLAGS = $(MEX_CFLAGS) $(MEX_LDFLAGS)

$(OBJ_MEX_OBJ_FILE) : $(MEX_OBJ_FILE)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_MEX_CONVERT_PRIMITIVES_FILE) : $(MEX_CONVERT_PRIMITIVES_FILE)
	$(CC) $(CFLAGS) -c $< -o $@
	
ttm: $(OBJS)
	$(MEX) $(MEX_FLAGS) $(TTM_FILE) $(OBJS) 

sttsm: $(OBJS)
	$(MEX) $(MEX_FLAGS) $(STTSM_FILE) $(OBJS)
	
tensor: $(OBJS)
	$(MEX) $(MEX_FLAGS) $(SRC_MATLAB_DIR)/TLA_tensor.c $(OBJS)

blockedtensor: $(OBJS)
	$(MEX) $(MEX_FLAGS) $(SRC_MATLAB_DIR)/TLA_blockedtensor.c $(OBJS)

blockedpsymtensor: $(OBJS)
	$(MEX) $(MEX_FLAGS) $(SRC_MATLAB_DIR)/TLA_blockedpsymtensor.c $(OBJS)

all: ttm sttsm
