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
# --- Include variables determined at configure-time --------------------------
#
CONFIGURE_DEFS = config\config.mk

!if exist ( $(CONFIGURE_DEFS) )
!include $(CONFIGURE_DEFS)
!else
!error nmake: $(CONFIGURE_DEFS) does not exist! Run configure.cmd first.
!endif



#
# --- Include environment- and build-specific definitions ----------------------
#

MAKE_DEFS = build\defs.mk

# Include build definitions
!if exist ( $(MAKE_DEFS) )
!include $(MAKE_DEFS)
!else
!error nmake: $(MAKE_DEFS) does not exist! Your libflame distribution may be incomplete.
!endif



#
# --- TLASH specific definitions

CFLAGS = $(CFLAGS) /I"$(FLAME_INC_DIR)"
LDFLAGS = $(LDFLAGS) /LIBPATH:"$(FLAME_LIB_DIR)" $(FLAME_LIB)


#
# --- BLAS definitions
#

LDFLAGS = $(LDFLAGS) /LIBPATH:"$(BLAS_LIB_DIR)" $(BLAS_LIB)


#
# --- MEX specific definitions
#

MEX_FLAGS = COMPFLAGS="$$COMPFLAGS $(CFLAGS)" LINKFLAGS="$$LINKFLAGS $(LDFLAGS)"

TTM_FILE = $(PWD)/../src/btas/ttm/flamec/TLA_ttm.c
STTSM_FILE = $(PWD)/../src/btas/sttsm/flamec/TLA_sttsm.c
OBJS = $(OBJ_FLA_DIRPATH)\TLA_MEX_Convert_primitives.obj $(OBJ_FLA_DIRPATH)\TLA_MEX_Obj.obj


#
# --- Variable modifications ---------------------------------------------------
#

LIBFLAME_OBJS = $(FLAMEC_OBJS)


#
# --- High-level rules ---------------------------------------------------------
#


all: sttsm ttm

sttsm: $(OBJ_FLA_DIRPATH)\TLA_MEX_Convert_primitives.obj $(OBJ_FLA_DIRPATH)\TLA_MEX_Obj.obj
	"$(MEX) $(MEX_FLAGS) $(STTSM_FILE) $(OBJS)" 
	
ttm: $(OBJ_FLA_DIRPATH)/TLA_MEX_Convert_primitives.obj $(OBJ_FLA_DIRPATH)/TLA_MEX_Obj.obj
	"$(MEX) $(MEX_FLAGS) $(TTM_FILE) $(OBJS)"
	 
tensor: $(OBJ_DIRPATH)/TLA_MEX_Convert_primitives.o $(OBJ_DIRPATH)/TLA_MEX_Obj.o
	"$(MEX) $(MEX_FLAGS) $(PWD)/src/TLA_tensor.c $(OBJS)"

blockedtensor: $(OBJ_DIRPATH)/TLA_MEX_Convert_primitives.o $(OBJ_DIRPATH)/TLA_MEX_Obj.o
	"$(MEX) $(MEX_FLAGS) $(PWD)/src/TLA_blockedtensor.c $(OBJS)"

blockedpsymtensor: $(OBJ_DIRPATH)/TLA_MEX_Convert_primitives.o $(OBJ_DIRPATH)/TLA_MEX_Obj.o
	"$(MEX) $(MEX_FLAGS) $(PWD)/src/TLA_blockedpsymtensor.c $(OBJS)"
	
#
# --- Source code (inference) rules --------------------------------------------
#

# --- C source files in flamec directory ---
{$(SRC_FLA_DIRPATH)}.c{$(OBJ_FLA_DIRPATH)}.obj:
!ifdef VERBOSE
	if not exist $(OBJ_FLA_DIRPATH) \
	   ( $(MKDIR) $(OBJ_FLA_DIRPATH) )
	   echo $(OBJ_FLA_DIRPATH)
	$(CC) $(MEX_INCLUDE_FLAGS) $(CFLAGS) /c $< /Fo$@ >> $(CC_LOG_FILE)
!else
	@if not exist $(OBJ_FLA_DIRPATH) \
	 (  ( $(ECHO) nmake: Creating $(OBJ_FLA_DIRPATH) directory ) & \
	    ( $(MKDIR) $(OBJ_FLA_DIRPATH) ) )
	@$(ECHO) nmake: Compiling $<
	@$(CC) $(MEX_INCLUDE_FLAGS) $(CFLAGS) /c $< /Fo$@ >> $(CC_LOG_FILE)
!endif

	

#
# --- Clean rules --------------------------------------------------------------
#

clean-log:
!ifdef VERBOSE
	if exist $(CC_LOG_FILE) \
	   ( $(DEL) $(CC_LOG_FILE) )
	if exist $(FC_LOG_FILE) \
	   ( $(DEL) $(FC_LOG_FILE) )
	if exist $(COPY_LOG_FILE) \
	   ( $(DEL) $(COPY_LOG_FILE) )
!else
	@if exist $(CC_LOG_FILE) \
	 (  ( $(ECHO) nmake: Deleting $(CC_LOG_FILE) ) & \
	    ( $(DEL) $(CC_LOG_FILE) ) )
	@if exist $(FC_LOG_FILE) \
	 (  ( $(ECHO) nmake: Deleting $(FC_LOG_FILE) ) & \
	    ( $(DEL) $(FC_LOG_FILE) ) )
	@if exist $(COPY_LOG_FILE) \
	 (  ( $(ECHO) nmake: Deleting $(COPY_LOG_FILE) ) & \
	    ( $(DEL) $(COPY_LOG_FILE) ) )
!endif

clean-config:
!ifdef VERBOSE
	if exist $(CNF_DIRNAME) \
	   ( $(RMDIR) $(CNF_DIRNAME) )
	if exist $(INC_DIRNAME) \
	   ( $(RMDIR) $(INC_DIRNAME) )
	if exist $(SRC_DIRNAME) \
	   ( $(RMDIR) $(SRC_DIRNAME) )
!else
	@if exist $(CNF_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(CNF_DIRNAME) directory ) & \
	    ( $(RMDIR) $(CNF_DIRNAME) ) )
	@if exist $(INC_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(INC_DIRNAME) directory ) & \
	    ( $(RMDIR) $(INC_DIRNAME) ) )
	@if exist $(SRC_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(SRC_DIRNAME) directory ) & \
	    ( $(RMDIR) $(SRC_DIRNAME) ) )
!endif

clean-build:
!ifdef VERBOSE
	if exist $(OBJ_DIRNAME) \
	   ( $(RMDIR) $(OBJ_DIRNAME) )
	if exist $(LIB_DIRNAME) \
	   ( $(RMDIR) $(LIB_DIRNAME) )
	if exist $(DLL_DIRNAME) \
	   ( $(RMDIR) $(DLL_DIRNAME) )
!else
	@if exist $(OBJ_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(OBJ_DIRNAME) directory ) & \
	    ( $(RMDIR) $(OBJ_DIRNAME) ) )
	@if exist $(LIB_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(LIB_DIRNAME) directory ) & \
	    ( $(RMDIR) $(LIB_DIRNAME) ) )
	@if exist $(DLL_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(DLL_DIRNAME) directory ) & \
	    ( $(RMDIR) $(DLL_DIRNAME) ) )
!endif

# Useful for developing when all we want to do is remove the library products.
clean-lib:
!ifdef VERBOSE
	if exist $(LIB_DIRNAME) \
	   ( $(RMDIR) $(LIB_DIRNAME) )
	if exist $(DLL_DIRNAME) \
	   ( $(RMDIR) $(DLL_DIRNAME) )
!else
	@if exist $(LIB_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(LIB_DIRNAME) directory ) & \
	    ( $(RMDIR) $(LIB_DIRNAME) ) )
	@if exist $(DLL_DIRNAME) \
	 (  ( $(ECHO) nmake: Deleting $(DLL_DIRNAME) directory ) & \
	    ( $(RMDIR) $(DLL_DIRNAME) ) )
!endif



#
# --- Help target --------------------------------------------------------------
#

help:
	@$(NMAKE_HELP)

