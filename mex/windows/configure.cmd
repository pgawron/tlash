::
:: libflame
:: An object-based infrastructure for developing high-performance
:: dense linear algebra libraries.
::
:: Copyright (C) 2011, The University of Texas
::
:: libflame is free software; you can redistribute it and/or modify
:: it under the terms of the GNU Lesser General Public License as
:: published by the Free Software Foundation; either version 2.1 of
:: the License, or (at your option) any later version.
::
:: libflame is distributed in the hope that it will be useful, but
:: WITHOUT ANY WARRANTY; without even the implied warranty of
:: MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
:: Lesser General Public License for more details.
::
:: You should have received a copy of the GNU Lesser General Public
:: License along with libflame; if you did not receive a copy, see
:: http://www.gnu.org/licenses/.
::
:: For more information, please contact us at flame@cs.utexas.edu or
:: send mail to:
::
:: Field G. Van Zee and/or
:: Robert A. van de Geijn
:: The University of Texas at Austin
:: Department of Computer Sciences
:: 1 University Station C0500
:: Austin TX 78712
::

@echo off

:ENVIRONMENT
	set GEN_CHECK_REV_FILE=.\build\gen-check-rev-file.py
	set GATHER_SRC=.\build\gather-src-for-windows.py
	set GEN_CONFIG_FILE=.\build\gen-config-file.py
	set CONFIG_DEFS_TEMPL=.\build\config.mk.in
	set SRC_TREE_DIR=..\src
	set TOP_BUILD_DIR=.

:PARAMS
	if "%1"=="" (goto USAGE)
	if "%2"=="" (goto USAGE)
	if "%3"=="" (goto USAGE)
	if %4==""   (goto USAGE)
	if %5==""   (goto USAGE)
	if %6==""   (goto USAGE)

	set ARCH=%1
	set BUILD=%2
	set CCOMPILER=%3
	set MATLAB_TOP_DIR=%4
	set FLAME_LIB_PATH=%5
	set BLAS_LIB_PATH=%6
	
:TASK_UNIT
	if "%ARCH%"=="x64" (goto ARCHERROR)
	
	echo %0: Checking/updating revision file.
	%GEN_CHECK_REV_FILE% -v
	echo %0: Gathering source files into local flat directories.
	%GATHER_SRC% %SRC_TREE_DIR% %TOP_BUILD_DIR%
	echo %0: Creating configure definitions file.
	%GEN_CONFIG_FILE% %TOP_BUILD_DIR% %ARCH% %BUILD% %CCOMPILER% %MATLAB_TOP_DIR% %FLAME_LIB_PATH% %BLAS_LIB_PATH% %CONFIG_DEFS_TEMPL%
	echo %0: Configuration and setup complete. You may now run nmake. 

	goto END

:USAGE
	echo. 
	echo  configure.cmd
	echo. 
	echo  A wrapper script for various configuration and setup scripts that need
	echo  to be run before nmake when building tlash MATLAB interface for 
	echo  Microsoft Windows.
	echo. 
	echo  USAGE:
	echo     %0 [arch] [build] [cc] [MATLAB_dir] [flame_path] [blas_path]
	echo.
	echo        arch       -- The architecture string to build.
	echo                      Supported values: {x86,x64}
	echo        build      -- The kind of build.
	echo                      Supported values: {debug,release}
	echo        cc         -- The C compiler to use.
	echo                      Supported values: {icl,cl}
	echo        MATLAB_dir -- The directory (in quotes) where MATLAB
	echo                      is installed (including MATLAB version)
	echo        flame_path -- The directory (in quotes) where
	echo                      libFLAME is installed
	echo        blas_path  -- The full path (in quotes) where
	echo                      BLAS library is installed
	echo. 
	echo  examples:
	echo     %0 x86 debug icl "C:\Program Files (x86)\MATLAB\R2013a" "C:\field\C:\field\lib\libflame\lib\libflame-x86-runknown.lib" "C:\lib\blas\x86\libblas.lib"
	echo     %0 x64 release cl "C:\Program Files\MATLAB\R2013a" "C:\field\C:\field\lib\libflame\lib\libflame-x64-runknown.lib" "C:\lib\blas\x64\libblas.lib"
	echo.
	
	goto END

:ARCHERROR
	echo x64 build for Windows unsupported at this time
	
	goto END

:END
