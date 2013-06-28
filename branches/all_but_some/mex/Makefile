CC = gcc

FLAME_DIR = /Users/martin/flame/lib
FLAME_LIB = flame
FLAME_INCLUDE = /Users/martin/flame/include

MEX_INCLUDE = /Applications/MATLAB_R2013a.app/extern/include

UTIL_INCLUDE = $(PWD)/src/base/flamec/tensor/include
UTIL_FILE = $(PWD)/src/base/flamec/tensor/util/base/TLA_MEX_Convert_primitives.c

MEX_OBJ_FILE = $(PWD)/src/base/flamec/tensor/main/TLA_MEX_Obj.c

TTM_FILE = $(PWD)/src/btas/ttm/flamec/TLA_ttm.c
STTSM_FILE = $(PWD)/src/btas/sttsm/flamec/TLA_sttsm.c

LDFLAGS = -L$(FLAME_DIR) -l$(FLAME_LIB)
CFLAGS = -g -I$(FLAME_INCLUDE) -I$(MEX_INCLUDE) -I$(UTIL_INCLUDE)

OBJ_PATH = $(PWD)/obj

OBJS = $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o

$(OBJ_PATH)/TLA_MEX_Obj.o : $(MEX_OBJ_FILE)
	$(CC) $(CFLAGS) -c $(MEX_OBJ_FILE) -o $(OBJ_PATH)/TLA_MEX_Obj.o

$(OBJ_PATH)/TLA_MEX_Convert_primitives.o : $(UTIL_FILE)
	$(CC) $(CFLAGS) -c $(UTIL_FILE) -o $(OBJ_PATH)/TLA_MEX_Convert_primitives.o
	
ttm: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex $(CFLAGS) $(TTM_FILE) $(OBJS) $(LDFLAGS) 

sttsm: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex $(CFLAGS) $(STTSM_FILE) $(OBJS) $(LDFLAGS)
	
tensor: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex $(CFLAGS) $(PWD)/src/TLA_tensor.c $(OBJS) $(LDFLAGS)

blockedtensor: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex $(CFLAGS) $(PWD)/src/TLA_blockedtensor.c $(OBJS) $(LDFLAGS)

blockedpsymtensor: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex $(CFLAGS) $(PWD)/src/TLA_blockedpsymtensor.c $(OBJS) $(LDFLAGS)

all: ttm sttsm
