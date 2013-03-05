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


void time_Apply_Q2_UT(
               int param_combo, int type, int nrepeats, int m,
               FLA_Obj A, FLA_Obj A_ref, FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj E,
               FLA_Obj T, FLA_Obj TB, FLA_Obj TD, FLA_Obj TE, FLA_Obj W,
               double *dtime, double *diff, double *gflops );


void time_Apply_Q2_UT(
               int param_combo, int type, int nrepeats, int m,
               FLA_Obj A, FLA_Obj A_ref, FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj E,
               FLA_Obj T, FLA_Obj TB, FLA_Obj TD, FLA_Obj TE, FLA_Obj W,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj A_save;

  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLASH_Copy( A, A_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( A_save, A );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        FLASH_QR2_UT( B,
                      D, TD );
        FLASH_Apply_Q2_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                             D, TD, W, C, E );
        break;
      case FLA_ALG_FRONT:
        FLASH_QR2_UT( B,
                      D, TD );
        FLASH_Apply_Q2_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                             D, TD, W, C, E );
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
    FLASH_Copy( A, A_ref );

    *diff = 0.0;
  }
  else
  {
    *diff = FLASH_Max_elemwise_diff( A, A_ref );
  }

  *gflops = 2.0 * m * m * (m - m/3.0) /
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Copy( A_save, A );

  FLASH_Obj_free( &A_save );
}

