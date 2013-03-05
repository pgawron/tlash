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

#ifndef BLIS_H
#define BLIS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Determine whether or not we are using BLIS from libflame.
//#define BLIS_FROM_LIBFLAME

#ifdef BLIS_FROM_LIBFLAME

  // If using libflame, pull in its header files so that
  // vector intrinsics-related macro constants are set properly.
  //#include "FLAME.h"
  #include "FLA_config.h"
  #include "FLA_macro_defs.h"
  #include "FLA_type_defs.h"

  // --- Pass-through macros for BLIS ---
  #ifdef FLA_ENABLE_CBLAS_INTERFACES
    #define BLIS_ENABLE_CBLAS_INTERFACES
  #endif
  #ifdef FLA_ENABLE_WINDOWS_BUILD
    #define BLIS_ENABLE_WINDOWS_BUILD
  #endif
  #ifdef FLA_ENABLE_UPPERCASE_F77
    #define BLIS_ENABLE_UPPERCASE_F77
  #endif
  #ifdef FLA_ENABLE_VECTOR_INTRINSICS
    #define BLIS_ENABLE_VECTOR_INTRINSICS
  #endif

  #define BLIS_VECTOR_INTRINSIC_TYPE FLA_VECTOR_INTRINSIC_TYPE

#else

  // --- BLIS configuration options ---

  // #define BLIS_ENABLE_USE_OF_FLA_MALLOC
  // #define BLIS_ENABLE_CBLAS_INTERFACES
  // #define BLIS_ENABLE_WINDOWS_BUILD
  // #define BLIS_ENABLE_UPPERCASE_F77
  // #define BLIS_ENABLE_VECTOR_INTRINSICS
  //   #define BLIS_VECTOR_INTRINSIC_TYPE BLIS_NO_INTRINSICS
  //   #define BLIS_VECTOR_INTRINSIC_TYPE BLIS_SSE_INTRINSICS

#endif

#include "blis_macro_defs.h"
#include "blis_type_defs.h"

#include "blis_prototypes_util.h"
#include "blis_prototypes_query.h"
#include "blis_prototypes_misc.h"

#include "blis_prototypes_level1.h"
#include "blis_prototypes_level2.h"
#include "blis_prototypes_level3.h"

#include "blis_prototypes_fused1.h"

#include "blis_f77_name_mangling.h"

#ifdef BLIS_ENABLE_CBLAS_INTERFACES
  #include "blis_prototypes_cblas.h"
#else
  #include "blis_prototypes_blas.h"
#endif

#endif
