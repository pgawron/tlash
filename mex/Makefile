CC = /usr/bin/gcc

FLAME_DIR = /Users/martin/flame/lib
FLAME_LIB = flame
FLAME_INCLUDE = /Users/martin/flame/include

MEX_INCLUDE = /Applications/MATLAB_R2012a.app/extern/include

UTIL_INCLUDE = $(PWD)/src/base/flamec/tensor/include
UTIL_FILE = $(PWD)/src/base/flamec/tensor/util/base/TLA_MEX_Convert_primitives.c

OBJ_FILE = $(PWD)/src/base/flamec/tensor/main/TLA_MEX_Obj.c

FILES = $(PWD)/src/btas/ttm/flamec/TLA_ttm.c

LDFLAGS = -L$(FLAME_DIR) -l$(FLAME_LIB)
CFLAGS = -I$(FLAME_INCLUDE) -I$(MEX_INCLUDE) -I$(UTIL_INCLUDE)

OBJ_PATH = $(PWD)/obj

OBJS = $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o

$(OBJ_PATH)/TLA_MEX_Obj.o : $(OBJ_FILE)
	$(CC) -ggdb $(CFLAGS) -c $(OBJ_FILE) -o $(OBJ_PATH)/TLA_MEX_Obj.o

$(OBJ_PATH)/TLA_MEX_Convert_primitives.o : $(UTIL_FILE)
	$(CC) -ggdb $(CFLAGS) -c $(UTIL_FILE) -o $(OBJ_PATH)/TLA_MEX_Convert_primitives.o
	
all: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex -g $(CFLAGS) $(FILES) $(OBJS) $(LDFLAGS) 

tensor: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex -g $(CFLAGS) $(PWD)/src/TLA_tensor.c $(OBJS) $(LDFLAGS)

blockedtensor: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex -g $(CFLAGS) $(PWD)/src/TLA_blockedtensor.c $(OBJS) $(LDFLAGS)

blockedpsymtensor: $(OBJ_PATH)/TLA_MEX_Convert_primitives.o $(OBJ_PATH)/TLA_MEX_Obj.o
	mex -g $(CFLAGS) $(PWD)/src/TLA_blockedpsymtensor.c $(OBJS) $(LDFLAGS)