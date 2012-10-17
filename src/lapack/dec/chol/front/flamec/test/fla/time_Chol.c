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


FLA_Error REF_Chol( FLA_Trans trans, FLA_Obj A );

void time_Chol(
                int param_combo, int type, int nrepeats, int m,
                FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                double *dtime, double *diff, double *gflops );


void time_Chol(
                int param_combo, int type, int nrepeats, int m,
                FLA_Obj A, FLA_Obj b, FLA_Obj b_orig, FLA_Obj norm,
                double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_save = 1.0e9;

  FLA_Obj
    A_save, b_save, b_orig_save;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b_orig, &b_orig_save );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( b, b_save );
  FLA_Copy_external( b_orig, b_orig_save );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Chol( FLA_LOWER_TRIANGULAR, A );
        break;
      case FLA_ALG_FRONT:
        FLA_Chol( FLA_LOWER_TRIANGULAR, A );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 1:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Chol( FLA_UPPER_TRIANGULAR, A );
        break;
      case FLA_ALG_FRONT:
        FLA_Chol( FLA_UPPER_TRIANGULAR, A );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_save = min( *dtime, dtime_save );
  }

  if ( param_combo == 0 )
  {
    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_UNIT_DIAG, A, b );
    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );

    FLA_Hemv_external( FLA_LOWER_TRIANGULAR,
                       FLA_ONE, A_save, b, FLA_MINUS_ONE, b_orig );

    FLA_Nrm2_external( b_orig, norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, norm,
                               1, 1, diff, 1, 1 );
  }
  else if ( param_combo == 1 )
  {
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                       FLA_UNIT_DIAG, A, b );
    FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, A, b );

    FLA_Hemv_external( FLA_UPPER_TRIANGULAR,
                       FLA_ONE, A_save, b, FLA_MINUS_ONE, b_orig );

    FLA_Nrm2_external( b_orig, norm );
    FLA_Copy_object_to_buffer( FLA_NO_TRANSPOSE, 0, 0, norm,
                               1, 1, diff, 1, 1 );
  }

  *gflops = 1.0 / 3.0  * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) / 
            dtime_save / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_save;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( b_save, b );
  FLA_Copy_external( b_orig_save, b_orig );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &b_save );
  FLA_Obj_free( &b_orig_save );
}

