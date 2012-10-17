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

	set ARCH=%1
	set BUILD=%2
	set CCOMPILER=%3
	
:TASK_UNIT
	echo %0: Checking/updating revision file.
	%GEN_CHECK_REV_FILE% -v
	echo %0: Gathering source files into local flat directories.
	%GATHER_SRC% %SRC_TREE_DIR% %TOP_BUILD_DIR%
	echo %0: Creating configure definitions file.
	%GEN_CONFIG_FILE% %TOP_BUILD_DIR% %ARCH% %BUILD% %CCOMPILER% %CONFIG_DEFS_TEMPL%
	echo %0: Configuration and setup complete. You may now run nmake. 

	goto END

:USAGE
	echo. 
	echo  configure.cmd
	echo. 
	echo  A wrapper script for various configuration and setup scripts that need
	echo. to be run before nmake when building libflame for Microsoft Windows.
	echo. 
	echo  USAGE:
	echo     %0 [arch] [build] [cc]
	echo.
	echo        arch     -- The architecture string to build.
	echo                    Supported values: {x86,x64}
	echo        build    -- The kind of build.
	echo                    Supported values: {debug,release}
	echo        cc       -- The C compiler to use.
	echo                    Supported values: {icl,cl}
	echo. 
	echo  examples:
	echo     %0 x86 debug icl
	echo     %0 x64 release cl
	echo.

:END
