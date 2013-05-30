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


FLA_Error REF_Trmm( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
void time_Trmm(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Trmm( 
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A, FLA_Obj C, FLA_Obj C_ref,
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
        REF_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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
        REF_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 4
    case 4:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 5
    case 5:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 6
    case 6:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 7
    case 7:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 8
    case 8:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 9
    case 9:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 10
    case 10:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    // Time parameter combination 11
    case 11:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        REF_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
        break;
      case FLA_ALG_FRONT:
        FLA_Trmm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_TWO, A, C );
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

  *gflops = 1.0 * 
            FLA_Obj_length( C ) * 
            FLA_Obj_width( C ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1.0e9;

  if ( param_combo == 0 ||
       param_combo == 3 ||
       param_combo == 6 ||
       param_combo == 9 )
  *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

