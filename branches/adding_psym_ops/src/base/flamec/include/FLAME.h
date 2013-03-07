/*
   libflame
   An object-based infrastructure for developing high-performance
   dense linear algebra libraries.

   Copyright (C) 2011, The University of Texas

   libflame is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   libflame is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with libflame; if you did not receive a copy, see
   http://www.gnu.org/licenses/.

   For more information, please contact us at flame@cs.utexas.edu or
   send mail to:

   Field G. Van Zee and/or
   Robert A. van de Geijn
   The University of Texas at Austin
   Department of Computer Sciences
   1 University Station C0500
   Austin TX 78712
*/

#ifndef FLAME_H
#define FLAME_H

// Allow C++ users to include this header file in their source code. However,
// we make the extern "C" conditional on whether we're using a C++ compiler,
// since regular C compilers don't understand the extern "C" construct.
#ifdef __cplusplus
extern "C" {
#endif 

  // Include autoconf-related preprocessor defines.
  #include "FLA_config.h"
  #include "FLA_config_check.h"

  // Include standard C header files.
  #include <stdio.h>
  #include <stdlib.h>
  #include <stdarg.h>
  #include <string.h>
  #ifdef FLA_ENABLE_WINDOWS_BUILD
    #include <windows.h>
  #else
  #ifndef FLA_ENABLE_TIDSP
    // TI CG does not support POSIX
    #include <unistd.h>
    #include <fcntl.h>
    #include <sys/types.h>
  #endif
  #endif
  #include <math.h>
  #include <float.h>
  #include <signal.h>

  // Include prototypes for BLAS-like interfaces.
  #ifndef BLIS_FROM_LIBFLAME
    #define BLIS_FROM_LIBFLAME
  #endif
  #include "blis.h"

  // Include f2c definitions.
  #include "FLA_f2c.h"

  // Include general FLAME macro and _PTR macro definitions.
  #include "FLA_macro_defs.h"
  #include "FLA_macro_ptr_defs.h"

  // Include general FLAME type definitions, including those for FLA_Obj.
  #include "FLA_type_defs.h"

  // Include "extern" definitions for global FLAME scalar constants.
  #include "FLA_extern_defs.h"

  // Include control tree structure definitions, utility prototypes, and
  // initialization prototypes.
  #include "FLA_Cntl.h"
  #include "FLA_Cntl_init.h"

  // Include prototypes for base FLAME routines.
  #include "FLA_main_prototypes.h"
  #include "FLA_util_base_prototypes.h"

  // Include prototypes for FLAME interfaces to BLAS and LAPACK operations.
  #include "FLA_blas1_prototypes.h"
  #include "FLA_blas2_prototypes.h"
  #include "FLA_blas3_prototypes.h"

  // Include prototypes for FLAME implementations of BLAS and LAPACK operations.
  #include "FLA_blas_var_prototypes.h"

  // Include FLASH headers.
  #include "FLASH.h"

  // Include SuperMatrix headers.
  #include "FLASH_Queue.h"

  // Include Fortran name-mangling macro (if not already defined).
  #include "FLA_f77_name_mangling.h"

  // Include tensor headers.
  #include "TLASH.h"

// End extern "C" construct block.
#ifdef __cplusplus
}
#endif 

#endif

