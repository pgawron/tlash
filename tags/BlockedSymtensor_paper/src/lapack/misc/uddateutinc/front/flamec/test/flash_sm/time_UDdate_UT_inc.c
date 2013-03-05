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


void time_UDdate_UT_inc(
                 int variant, int type, int n_repeats, int mB, int mC, int mD, int n,
                 FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops );


void time_UDdate_UT_inc(
                 int variant, int type, int n_repeats, int mB, int mC, int mD, int n,
                 FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    RR, EE, R_save, C_save, D_save;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, &RR );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, E, &EE );

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, &R_save );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_save );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, D, &D_save );

  FLASH_Copy( R, R_save );
  FLASH_Copy( C, C_save );
  FLASH_Copy( D, D_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLASH_Copy( R_save, R );
    FLASH_Copy( C_save, C );
    FLASH_Copy( D_save, D );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:{
      switch( type ){
      case FLA_ALG_FRONT:
        FLASH_UDdate_UT_inc( R, C, D, T, W );

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

  {
    FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, R, R, FLA_ZERO, RR );
    FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, E, E, FLA_ZERO, EE );

    *diff = FLASH_Max_elemwise_diff( RR, EE );
  }

  *gflops = 2 * ( ( mC + mD ) * n * n +
                  ( mC + mD ) * n * 6 ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( R ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( R_save, R );
  FLASH_Copy( C_save, C );
  FLASH_Copy( D_save, D );

  FLASH_Obj_free( &R_save );
  FLASH_Obj_free( &C_save );
  FLASH_Obj_free( &D_save );

  FLASH_Obj_free( &RR );
  FLASH_Obj_free( &EE );
}

