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

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_UNB_OPT_FUSED 3
#define FLA_ALG_BLOCKED       4
#define FLA_ALG_BLOCKED_FUSED 5

FLA_Error REF_Hess_UT( FLA_Obj A, FLA_Obj t );
void time_Hess(
               int variant, int type, int n_repeats, int m, int nb_alg,
               FLA_Obj A, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj TT, FLA_Obj t,
               double *dtime, double *diff, double *gflops );


void time_Hess(
               int variant, int type, int n_repeats, int m, int nb_alg,
               FLA_Obj A, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj TT, FLA_Obj t,
               double *dtime, double *diff, double *gflops )
{
  int irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_save, norm;

  if ( 
       ( variant == 1 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       ( variant == 5 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       ( variant == 1 && type == FLA_ALG_BLOCKED_FUSED ) ||
       ( variant == 5 && type == FLA_ALG_BLOCKED_FUSED ) ||
       FALSE
     )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }

  if (
       ( variant == 0 ) ||

       ( variant == 1 && type == FLA_ALG_UNBLOCKED ) ||
       ( variant == 2 && type == FLA_ALG_UNBLOCKED ) ||
       ( variant == 3 && type == FLA_ALG_UNBLOCKED ) ||
       //( variant == 4 && type == FLA_ALG_UNBLOCKED ) ||
       ( variant == 5 && type == FLA_ALG_UNBLOCKED ) ||

       ( variant == 1 && type == FLA_ALG_UNB_OPT ) ||
       ( variant == 2 && type == FLA_ALG_UNB_OPT ) ||
       ( variant == 3 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 4 && type == FLA_ALG_UNB_OPT ) ||
       ( variant == 5 && type == FLA_ALG_UNB_OPT ) ||

       ( variant == 2 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       ( variant == 3 && type == FLA_ALG_UNB_OPT_FUSED ) ||
       //( variant == 4 && type == FLA_ALG_UNB_OPT_FUSED ) ||

       ( variant == 1 && type == FLA_ALG_BLOCKED ) ||
       ( variant == 2 && type == FLA_ALG_BLOCKED ) ||
       ( variant == 3 && type == FLA_ALG_BLOCKED ) ||
       //( variant == 4 && type == FLA_ALG_BLOCKED ) ||
       ( variant == 5 && type == FLA_ALG_BLOCKED ) ||

       ( variant == 2 && type == FLA_ALG_BLOCKED_FUSED ) ||
       ( variant == 3 && type == FLA_ALG_BLOCKED_FUSED ) ||
       //( variant == 4 && type == FLA_ALG_BLOCKED_FUSED ) ||

       FALSE
  )
  {
    *gflops = 0.0;
    *diff   = 0.0;
    return;
  }


  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

  FLA_Copy_external( A, A_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Hess_UT( A, t );
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Hess_UT_unb_var1( A, T );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Hess_UT_opt_var1( A, T );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        break;
      case FLA_ALG_BLOCKED:
        FLA_Hess_UT_blk_var1( A, TT );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        break;
      }
      break;
    }

    // Time variant 2
    case 2:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Hess_UT_unb_var2( A, T );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Hess_UT_opt_var2( A, T );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        FLA_Hess_UT_ofu_var2( A, T );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Hess_UT_blk_var2( A, TT );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        FLA_Hess_UT_blf_var2( A, TT );
        break;
      }
      break;
    }

    // Time variant 3
    case 3:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Hess_UT_unb_var3( A, T );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Hess_UT_opt_var3( A, T );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        FLA_Hess_UT_ofu_var3( A, T );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Hess_UT_blk_var3( A, TT );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        FLA_Hess_UT_blf_var3( A, TT );
        break;
      }
      break;
    }

    // Time variant 4
    case 4:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Hess_UT_unb_var4( A, T );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Hess_UT_opt_var4( A, T );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        FLA_Hess_UT_ofu_var4( A, T );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Hess_UT_blk_var4( A, TT );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        FLA_Hess_UT_blf_var4( A, TT );
        break;
      }
      break;
    }

    // Time variant 5
    case 5:
    {
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Hess_UT_unb_var5( A, T );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_Hess_UT_opt_var5( A, T );
        break;
      case FLA_ALG_UNB_OPT_FUSED:
        break;
      case FLA_ALG_BLOCKED:
        FLA_Hess_UT_blk_var5( A, TT );
        break;
      case FLA_ALG_BLOCKED_FUSED:
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }

  if ( variant > 0 )
  {

/*
if ( type == FLA_ALG_UNBLOCKED && variant <= 2 )
{
FLA_Obj_show( "A", A, "%10.3e + %10.3e ", "" );
FLA_Obj_show( "T", T, "%10.3e + %10.3e ", "" );
}
*/


    FLA_Obj    AT, AB;
    FLA_Obj Q, QT, QB;
    FLA_Obj E, ET, EB;
    FLA_Obj F;
    FLA_Obj W, WW, eye;
    dim_t   m_A, n_Q, m_T, m_TT;

//FLA_Obj_show( "A", A, "%10.3e", "" );
//FLA_Obj_show( "T", T, "%10.3e", "" );

    m_A = FLA_Obj_length( A );
    m_T = FLA_Obj_length( T );
    m_TT = FLA_Obj_length( TT );

    FLA_Obj_create( FLA_Obj_datatype( A ), m_A,  m_A, 0, 0, &Q );
    FLA_Obj_create( FLA_Obj_datatype( A ), m_T,  m_A, 0, 0, &W );
    FLA_Obj_create( FLA_Obj_datatype( A ), m_TT, m_A, 0, 0, &WW );
    FLA_Set_to_identity( Q );

    FLA_Part_2x1( Q,   &QT,
                       &QB,   1, FLA_TOP );
    FLA_Part_2x1( A,   &AT,
                       &AB,   1, FLA_TOP );

    if ( type == FLA_ALG_BLOCKED || type == FLA_ALG_BLOCKED_FUSED )
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, AB, TT, WW, QB );
    else
      FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, AB, T, W, QB );

/*
    n_Q = FLA_Obj_width( QB );

    FLA_Obj_create( FLA_Obj_datatype( A ), n_Q, n_Q, 0, 0, &eye );
    FLA_Set_to_identity( eye );

    FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, Q, Q, FLA_MINUS_ONE, eye );

    FLA_Norm_frob( eye, norm );
    FLA_Obj_extract_real_scalar( norm, diff );
*/

    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &E );     
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &F );

    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, A_save, Q, FLA_ZERO, E );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, Q, E, FLA_ZERO, F );
/*
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, A_save, Q, FLA_ZERO, E );
    FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, Q, E, FLA_ZERO, F );
*/

    FLA_Copy( A, E );
    FLA_Part_2x1( E,    & ET,
                        & EB,      1, FLA_TOP );
    FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, EB );

    *diff = FLA_Max_elemwise_diff( E, F );

    //FLA_Obj_free( &T );
    FLA_Obj_free( &W );
    FLA_Obj_free( &WW );
    FLA_Obj_free( &E );
    FLA_Obj_free( &F );
    FLA_Obj_free( &Q );
    //FLA_Obj_free( &eye );
  }
  else
  {
    *diff = 0.0;
  }

  *gflops = 10.0 / 3.0 * m * m * m /
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &norm );
}

