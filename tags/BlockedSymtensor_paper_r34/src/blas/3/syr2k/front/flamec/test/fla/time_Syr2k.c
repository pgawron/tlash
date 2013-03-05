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


FLA_Error REF_Syr2k( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C );
void time_Syr2k(
               int param_combo, int type, int nrepeats, int m, int k,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Syr2k( 
               int param_combo, int type, int nrepeats, int m, int k,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    C_old;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );


  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Syr2k( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Syr2k( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 1
    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Syr2k( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Syr2k( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 2
    case 2:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Syr2k( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Syr2k( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 3
    case 3:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Syr2k( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Syr2k( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, C );
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
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }

  *gflops = 2.0 * m * m * k /
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( C ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

