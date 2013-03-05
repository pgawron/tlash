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


#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1


FLA_Error REF_Lyap( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale );

void time_Lyap(
                int param_combo, int type, int nrepeats, int m,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale,
                double *dtime, double *diff, double *gflops );


void time_Lyap(
                int param_combo, int type, int nrepeats, int m,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale,
                double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_save, norm;

  if ( param_combo == 0 && type == FLA_ALG_FRONT )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_save );
  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( C ), 1, 1, 0, 0, &norm );

  FLASH_Copy( C, C_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( C_save, C );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:{
      switch( type ){
      //case FLA_ALG_REFERENCE:
      //  REF_Lyap( FLA_NO_TRANSPOSE, isgn, A_flat, C_flat, scale );
      //  break;
      case FLA_ALG_FRONT:
        FLASH_Lyap( FLA_NO_TRANSPOSE, isgn, A, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      //case FLA_ALG_REFERENCE:
      //  REF_Lyap( FLA_CONJ_TRANSPOSE, isgn, A_flat, C_flat, scale );
      //  break;
      case FLA_ALG_FRONT:
        FLASH_Lyap( FLA_CONJ_TRANSPOSE, isgn, A, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );
  }

/*
  if ( type == FLA_ALG_REFERENCE )
  {
    FLASH_Obj_hierarchify( C_flat, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLASH_Max_elemwise_diff( C, C_ref );
  }
*/

  {
    FLA_Obj X, W;

    FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &X );
    FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &W );

    FLASH_Copy( C, X );
    FLASH_Hermitianize( FLA_UPPER_TRIANGULAR, X );

    if ( param_combo == 0 )
    {
      FLASH_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,   FLA_ONE, A, X, FLA_ZERO, W );
      FLASH_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, X, A, FLA_ONE,  W );
    }
    else if ( param_combo == 1 )
    {
      FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, X, FLA_ZERO, W );
      FLASH_Gemm( FLA_NO_TRANSPOSE,   FLA_NO_TRANSPOSE, FLA_ONE, X, A, FLA_ONE,  W );
    }

    FLASH_Scal( isgn, W );

    FLASH_Axpy( FLA_MINUS_ONE, C_save, W );
    FLASH_Norm1( W, norm );
    FLA_Obj_extract_real_scalar( norm, diff );

    FLASH_Obj_free( &X );
    FLASH_Obj_free( &W );
  }

  *gflops = ( 2.0 / 3.0 ) * ( m * m * m ) / 
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( C ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( C_save, C );

  FLASH_Obj_free( &C_save );
  FLA_Obj_free( &norm );
}

