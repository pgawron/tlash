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

#include "FLAME.h"


static FLA_Bool FLA_initialized = FALSE;

FLA_Obj FLA_THREE;
FLA_Obj FLA_TWO;
FLA_Obj FLA_ONE;
FLA_Obj FLA_ONE_HALF;
FLA_Obj FLA_ZERO;
FLA_Obj FLA_MINUS_ONE_HALF;
FLA_Obj FLA_MINUS_ONE;
FLA_Obj FLA_MINUS_TWO;
FLA_Obj FLA_MINUS_THREE;




/* *************************************************************************

   FLA_Init()

 *************************************************************************** */

void FLA_Init()
{
  if ( FLA_initialized == TRUE ) return;
  
  FLA_initialized = TRUE;

  FLA_Error_messages_init();

  FLA_Memory_leak_counter_init();

  FLA_Init_constants();

  FLA_Cntl_init();

#if FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS
  _MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
#endif

#ifdef FLA_ENABLE_SUPERMATRIX
  FLASH_Queue_init();
#endif
}

/* *************************************************************************

  FLA_Finalize()

 *************************************************************************** */

void FLA_Finalize()
{
  if ( FLA_initialized == FALSE ) return;

  FLA_initialized = FALSE;

  FLA_Finalize_constants();

  FLA_Cntl_finalize();

#ifdef FLA_ENABLE_SUPERMATRIX
  FLASH_Queue_finalize();
#endif

  FLA_Memory_leak_counter_finalize();
}

/* *************************************************************************

  FLA_Init_safe()

 *************************************************************************** */

void FLA_Init_safe( FLA_Error* init_result )
{
  if ( FLA_Initialized() )
  {
    *init_result = FLA_FAILURE;
  }
  else
  {
    FLA_Init();
    *init_result = FLA_SUCCESS;
  }
}

/* *************************************************************************

  FLA_Finalize_safe()

 *************************************************************************** */

void FLA_Finalize_safe( FLA_Error init_result )
{
  if ( init_result == FLA_SUCCESS )
    FLA_Finalize();
}

/* *************************************************************************

   FLA_Initialized()

 *************************************************************************** */

FLA_Bool FLA_Initialized( void )
{
  return FLA_initialized;
}

/* *************************************************************************

   FLA_Init_constants()

 *************************************************************************** */

void FLA_Init_constants()
{
  FLA_Obj_create_constant(  3.0, &FLA_THREE );
  FLA_Obj_create_constant(  2.0, &FLA_TWO );
  FLA_Obj_create_constant(  1.0, &FLA_ONE );
  FLA_Obj_create_constant(  0.5, &FLA_ONE_HALF );
  FLA_Obj_create_constant(  0.0, &FLA_ZERO );
  FLA_Obj_create_constant( -0.5, &FLA_MINUS_ONE_HALF );
  FLA_Obj_create_constant( -1.0, &FLA_MINUS_ONE );
  FLA_Obj_create_constant( -2.0, &FLA_MINUS_TWO );
  FLA_Obj_create_constant( -3.0, &FLA_MINUS_THREE );
}

/* *************************************************************************

   FLA_Finalize_constants()

 *************************************************************************** */

void FLA_Finalize_constants()
{
  FLA_Obj_free( &FLA_THREE );
  FLA_Obj_free( &FLA_TWO );
  FLA_Obj_free( &FLA_ONE );
  FLA_Obj_free( &FLA_ONE_HALF );
  FLA_Obj_free( &FLA_ZERO );
  FLA_Obj_free( &FLA_MINUS_ONE_HALF );
  FLA_Obj_free( &FLA_MINUS_ONE );
  FLA_Obj_free( &FLA_MINUS_TWO );
  FLA_Obj_free( &FLA_MINUS_THREE );
}

