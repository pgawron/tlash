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
# --- Include environment- and build-specific definitions ----------------------
#

CONFIG_MK = ..\..\windows\config\config.mk
DEFS_MK   = ..\..\windows\build\defs.mk

# Include build definitions
!if exist ( $(CONFIG_MK) )
!include $(CONFIG_MK)
!else
!error nmake: $(CONFIG_MK) does not exist! Run configure.cmd first.
!endif

# Include build definitions
!if exist ( $(DEFS_MK) )
!include $(DEFS_MK)
!else
!error nmake: $(DEFS_MK) does not exist! Your libflame distribution may be incomplete.
!endif



#
# --- General build definitions ------------------------------------------------
#

# Path to test suite source code.
SRC_PATH         = ..\src

# Path to directory where test suite object files will be stored.
OBJ_PATH         = obj

# Path to libflame and local headers.
LIBFLAME_HEADERS = $(INSTALL_PREFIX_INC)
LOCAL_HEADERS    = $(SRC_PATH)

# The installation paths for the static and dynamic libflame library builds.
LIB_LIBPATH      = C:\tlash_trunk\windows\lib\x86\debug
DLL_LIBPATH      = $(INSTALL_PREFIX_DLL)

# We have to reconstruct the C compiler flags with the new include paths.
CPPROCFLAGS      = /I$(LIBFLAME_HEADERS) /I$(LOCAL_HEADERS)
CFLAGS           = $(CMISCFLAGS) $(CLANGFLAGS) $(CPPROCFLAGS) $(CWARNFLAGS) \
                   $(CDBGFLAGS) $(COPTFLAGS) $(CRTIMEFLAGS) $(CMTHREADFLAGS)

# We use the C compiler as our linker.
LINKER           = $(CC)

# The linker flags. This might need to be modified depending on your build
# environment.
LDFLAGS          = /nologo /link \
                   /LIBPATH:"C:\Program Files\Microsoft SDKs\Windows\v7.0A\Lib" \
                   /LIBPATH:"C:\Program Files\Microsoft Visual Studio 10.0\VC\lib" \
                   /nodefaultlib:libcmt /nodefaultlib:libc /nodefaultlib:libmmt \
                   msvcrt.lib \
				   /LIBPATH:"C:\tlash_trunk\windows\lib\x86\debug" \
				   libblas.lib

TEST_OBJS        = $(OBJ_PATH)\test_sttsm.obj

# Binary path. We put it one directory up because that is where the input
# paremeter files reside.
TEST_BIN         = ..\test_sttsm.exe

# The inference rule to compile source to object files.
{$(SRC_PATH)}.c{$(OBJ_PATH)}.obj:
	$(CC) $(CFLAGS) /c $< /Fo$@

# The rule to build the test suite driver by linking against a static
# build of libflame. This is the default.
with-lib: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) /Fe$(TEST_BIN) $(LDFLAGS) $(LIB_LIBPATH)\$(LIBFLAME_LIB)

# The rule to build the test suite driver by linking against a dynamic
# build of libflame.
with-dll: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) /Fe$(TEST_BIN) $(LDFLAGS) $(DLL_LIBPATH)\$(LIBFLAME_LIB)

clean:
	del /F /Q $(OBJ_PATH)\*.obj
	del /F /Q $(TEST_BIN)

