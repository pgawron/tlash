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

#ifndef BLIS_TYPE_DEFS_H
#define BLIS_TYPE_DEFS_H

// --- Basic type definitions -------------------------------------------------

/*
// --- trans ---

#define BLIS_NO_TRANSPOSE      'n'
#define BLIS_TRANSPOSE         't'
#define BLIS_CONJ_NO_TRANSPOSE 'c'
#define BLIS_CONJ_TRANSPOSE    'h'

// --- conj ---

#define BLIS_NO_CONJUGATE      'n'
#define BLIS_CONJUGATE         'c'

// --- uplo ---

#define BLIS_LOWER_TRIANGULAR  'l'
#define BLIS_UPPER_TRIANGULAR  'u'

// --- side ---

#define BLIS_LEFT              'l'
#define BLIS_RIGHT             'r'

// --- diag ---

#define BLIS_NONUNIT_DIAG      'n'
#define BLIS_UNIT_DIAG         'u'
#define BLIS_ZERO_DIAG         'z'
*/

#define BLIS_TRANS_BEGIN 100
#define BLIS_UPLO_BEGIN  200
#define BLIS_SIDE_BEGIN  300
#define BLIS_DIAG_BEGIN  400
#define BLIS_CONJ_BEGIN  500

typedef enum
{
	BLIS_NO_TRANSPOSE = BLIS_TRANS_BEGIN,
	BLIS_TRANSPOSE,
	BLIS_CONJ_NO_TRANSPOSE,
	BLIS_CONJ_TRANSPOSE
} trans_t;

typedef enum
{
	BLIS_LOWER_TRIANGULAR = BLIS_UPLO_BEGIN,
	BLIS_UPPER_TRIANGULAR
} uplo_t;

typedef enum
{
	BLIS_LEFT = BLIS_SIDE_BEGIN,
	BLIS_RIGHT
} side_t;

typedef enum
{
	BLIS_NONUNIT_DIAG = BLIS_DIAG_BEGIN,
	BLIS_UNIT_DIAG,
	BLIS_ZERO_DIAG
} diag_t;

typedef enum
{
	BLIS_NO_CONJUGATE = BLIS_CONJ_BEGIN,
	BLIS_CONJUGATE
} conj_t;





// --- Intrinsic/assembly definitions ----------------------------------------

/*
#ifndef BLIS_FROM_LIBFLAME
    #error "NOT using blis from libflame"
#else
    #error "using blis from libflame"
#endif
*/

/*
#if BLIS_VECTOR_INTRINSIC_TYPE == BLIS_SSE_INTRINSICS
    #error "using sse in blis"
#elif  BLIS_VECTOR_INTRINSIC_TYPE == BLIS_NO_INTRINSICS
    #error "NOT using sse in blis"
#else
    #error "undefined case!"
#endif
*/

// Only define vector intrinsics types if they are not already provided by
// libflame.
#ifndef BLIS_FROM_LIBFLAME

#if BLIS_VECTOR_INTRINSIC_TYPE == BLIS_SSE_INTRINSICS

#include "pmmintrin.h"
typedef union
{
    __m128d v; 
    double  d[2];
} v2df_t;
#endif

#endif


// --- Complex type definitions -----------------------------------------------

// Only define complex types if they are not already provided by libflame.
//#ifndef BLIS_ENABLE_USE_OF_LIBFLAME_TYPES
#ifndef BLIS_FROM_LIBFLAME

typedef struct scomplex
{
  float real, imag;
} scomplex;

typedef struct dcomplex
{
  double real, imag;
} dcomplex;

#endif


#endif // BLIS_TYPE_DEFS_H
