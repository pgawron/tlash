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

echo. 
echo  Makefile
echo. 
echo  Field G. Van Zee
echo.  
echo  nmake Makefile for building libflame for Microsoft Windows. nmake targets
echo  may be invoked after running the configure.cmd script. Valid targets are:
echo. 
echo    all          - Invoke the lib and dll targets.
echo    lib          - Build libflame as a static library.
echo    dll          - Build libflame as a dynamically-linked library.
echo    help         - Output help and usage information.
echo    clean        - Invoke clean-log and clean-build targets.
echo    clean-log    - Remove any log files present.
echo    clean-config - Remove all products of configure.cmd. Namely, remove the
echo                   config, include, and src directories.
echo    clean-build  - Remove all products of the compilation portion of the build
echo                   process. Namely, remove the obj, lib, and dll directories.
echo    distclean    - Invoke clean-log, clean-config, and clean-build targets.
echo.
echo  The Makefile also recognizes configuration options corresponding to the
echo  following Makefile variables:
echo.
echo    VERBOSE               - When defined, nmake outputs the actual commands
echo                            executed instead of more concise one-line progress
echo                            indicators. (Undefined by default.)
echo.
echo  Typically, these options are specified by commenting or uncommenting the
echo  corresponding lines in the Makefile. However, if the Makefile currently does
echo  not define one of the options, and you wish to enable the corresponding
echo  feature without editing the Makefile, you may define the variable at the
echo  command line when nmake is invoked. For example, you may enable verboseness
echo  while invoking the lib target as follows:
echo.
echo    nmake lib VERBOSE=1
echo.
