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


FLA_Error REF_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );

void time_Sylv(
                int param_combo, int type, int nrepeats, int m, int n,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
                double *dtime, double *diff, double *gflops );


void time_Sylv(
                int param_combo, int type, int nrepeats, int m, int n,
                FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref, FLA_Obj scale,
                double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old, A_flat, B_flat, C_flat;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, A, &A_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, B, &B_flat );
  FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, C, &C_flat );

  FLASH_Copy( C, C_old );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( C_old, C );
    FLASH_Obj_flatten( A, A_flat );
    FLASH_Obj_flatten( B, B_flat );
    FLASH_Obj_flatten( C, C_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A_flat, B_flat, C_flat, scale );
        break;
      case FLA_ALG_FRONT:
        FLASH_Sylv( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, isgn, A_flat, B_flat, C_flat, scale );
        break;
      case FLA_ALG_FRONT:
        FLASH_Sylv( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, isgn, A, B, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A_flat, B_flat, C_flat, scale );
        break;
      case FLA_ALG_FRONT:
        FLASH_Sylv( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 3:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Sylv( FLA_TRANSPOSE, FLA_TRANSPOSE, isgn, A_flat, B_flat, C_flat, scale );
        break;
      case FLA_ALG_FRONT:
        FLASH_Sylv( FLA_TRANSPOSE, FLA_TRANSPOSE, isgn, A, B, C, scale );
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

  if ( type == FLA_ALG_REFERENCE )
  {
    FLASH_Obj_hierarchify( C_flat, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLASH_Max_elemwise_diff( C, C_ref );
  }

  *gflops = ( m * m * n + n * n * m ) / 
            dtime_old / 1e9;

  *dtime = dtime_old;

  FLASH_Copy( C_old, C );

  FLASH_Obj_free( &C_old );
  FLASH_Obj_free( &A_flat );
  FLASH_Obj_free( &B_flat );
  FLASH_Obj_free( &C_flat );
}
