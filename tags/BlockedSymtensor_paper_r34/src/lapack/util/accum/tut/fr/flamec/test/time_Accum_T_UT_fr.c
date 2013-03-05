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
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_UNB_OPT   2


FLA_Error REF_Accum_T_UT_fr( FLA_Obj A, FLA_Obj t );
void time_Accum_T_UT_fr(
                 int variant, int type, int nrepeats, int m, int n,
                 FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W, FLA_Obj b, FLA_Obj b_ref,
                 double *dtime, double *diff, double *gflops );


void time_Accum_T_UT_fr(
                 int variant, int type, int nrepeats, int m, int n,
                 FLA_Obj A, FLA_Obj t, FLA_Obj T, FLA_Obj W, FLA_Obj b, FLA_Obj b_ref,
                 double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, b_save, norm;


  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, b, &b_save );

  if ( FLA_Obj_is_single_precision( A ) )
    FLA_Obj_create( FLA_FLOAT, 1, 1, 0, 0, &norm );
  else
    FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( b, b_save );


  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:{
      switch( type ){
      case FLA_ALG_REFERENCE:
        FLA_LQ_UT( A, T );
        break;
      case FLA_ALG_UNBLOCKED:
        FLA_LQ_UT( A, T );
        FLA_LQ_UT_recover_tau( T, t );
        FLA_Set( FLA_ZERO, T );
        FLA_Accum_T_UT_fr_unb_var1( A, t, T );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_LQ_UT( A, T );
        FLA_LQ_UT_recover_tau( T, t );
        FLA_Set( FLA_ZERO, T );
        FLA_Accum_T_UT_fr_opt_var1( A, t, T );
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
    FLA_Copy_external( b, b_ref );
    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A, b );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE, A, T, W, b );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, b, FLA_ONE, b_ref );
    FLA_Nrm2_external( b_ref, norm );
    if ( FLA_Obj_is_single_precision( A ) )
      *diff = *(FLA_FLOAT_PTR(norm));
    else
      *diff = *(FLA_DOUBLE_PTR(norm));
  }
  else
  {
    FLA_Copy_external( b, b_ref );
    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                       FLA_NONUNIT_DIAG, FLA_ONE, A, b );
    FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE, A, T, W, b );
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A_save, b, FLA_ONE, b_ref );
    FLA_Nrm2_external( b_ref, norm );
    if ( FLA_Obj_is_single_precision( A ) )
      *diff = *(FLA_FLOAT_PTR(norm));
    else
      *diff = *(FLA_DOUBLE_PTR(norm));
  }

  *gflops = 2.0 * n * n * (m - n/3.0) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( b_save, b );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &b_save );
  FLA_Obj_free( &norm );
}

