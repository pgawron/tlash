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

FLA_Error FLA_Apply_Q_UT_rnbr_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
/*
  Apply a unitary matrix Q to a matrix B from the right,

    B :=  B Q

  where Q is the backward product of Householder transformations:

    Q  =  H(k-1) ... H(1) H(0)

  where H(i) corresponds to the Householder vector stored above the diagonal
  in the ith row of A. Thus, the operation becomes:

    B :=  B Q
       =  B ( H(k-1) ... H(1) H(0) )
       =  B ( H(k-1)' ... H(1)' H(0)' )
       =  B ( H(0) H(1) ... H(k-1) )'
       =  B H(k-1)' ... H(1)' H(0)'

  From this, we can see that we must move through A from bottom-right to top-
  left, since the Householder vector for H(k-1) was stored in the last row
  of A. We intend to apply blocks of reflectors at a time, where a block
  reflector H of b consecutive Householder transforms may be expressed as:

    H  =  ( H(i) H(i+1) ... H(i+b-1) )'
       =  ( I - U inv(T) U' )'

  where:
    - U^T is the strictly upper trapezoidal (with implicit unit diagonal) matrix
      of Householder vectors, stored above the diagonal of A in rows i through
      i+b-1, corresponding to H(i) through H(i+b-1).
    - T is the upper triangular block Householder matrix corresponding to
      Householder vectors i through i+b-1.

  Consider applying H to B as an intermediate step towards applying all of Q:

    B  :=  B H
        =  B ( I - U inv(T) U' )'
        =  B ( I - U inv(T)' U' )
        =  B - B U inv(T)' U'

  We must move from bottom-right to top-left. So, we partition:

    U^T -> ( U11 U12 )  B -> ( B1 B2 )  T -> ( T2 T1 )

  where:
    - U11 is stored in strictly upper triangle of A11 with implicit unit
      diagonal.
    - U12 is stored in A12.
    - T1 is an upper triangular block of row-panel matrix T.

  Substituting repartitioned U, B, and T, we have:

    ( B1 B2 ) :=  ( B1 B2 ) - ( B1 B2 ) ( U11 U12 )^T inv(T1)' conj( U11 U12 )
               =  ( B1 B2 ) - ( B1 B2 ) / U11^T \ inv(T1)' conj( U11 U12 )
                                        \ U12^T /
               =  ( B1 B2 ) - ( B1 U11^T + B2 U12^T ) inv(T1)' conj( U11 U12 )

  Thus, B1 is updated as:

      B1    :=      B1   -   ( B1 U11^T + B2 U12^T ) inv(T1)' conj(U11)

  And B2 is updated as:

      B2    :=      B2   -   ( B1 U11^T + B2 U12^T ) inv(T1)' conj(U12)

  Note that:

    ( B1 U11^T + B2 U12^T ) inv(T1)'

  is common to both updates, and thus may be computed and stored in
  workspace, and then re-used.

  -FGVZ
*/
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj T1T,
          T2B;

  FLA_Obj WTL,  WTR,
          WBL,  WBR;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  dim_t   b_alg, b;
  dim_t   m_BR, n_BR;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  // If m < n, then we have to initialize our partitionings carefully so
  // that we begin in the proper location in A and B (since we traverse
  // matrix A from BR to TL).
  if ( FLA_Obj_length( A ) < FLA_Obj_width( A ) )
  {
    m_BR = 0;
    n_BR = FLA_Obj_width( A ) - FLA_Obj_length( A );
  }
  else
  {
    m_BR = 0;
    n_BR = 0;
  }

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     m_BR, n_BR, FLA_BR );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_RIGHT );

  FLA_Part_1x2( B,    &BL,  &BR,      n_BR, FLA_RIGHT );

  while ( FLA_Obj_min_dim( ATL ) > 0 ){

    b = min( b_alg, FLA_Obj_min_dim( ATL ) );

    // Since T was filled from left to right, and since we need to access them
    // in reverse order, we need to handle the case where the last block is
    // smaller than the other b x b blocks.
    if ( FLA_Obj_width( TR ) == 0 && FLA_Obj_width( T ) % b_alg > 0 )
      b = FLA_Obj_width( T ) % b_alg;

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, &T1, /**/ &T2,
                           b, FLA_LEFT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &B1, /**/ &B2,
                           b, FLA_LEFT );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T, 
                         &T2B,     b, FLA_TOP );

    FLA_Part_2x2( W,     &WTL, &WTR,
                         &WBL, &WBR,     b, FLA_Obj_length( B1 ), FLA_TL );

    // WTL = B1^T;

    FLA_Copyt_internal( FLA_TRANSPOSE, B1, WTL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // U11 = triuu( A11 );
    // U12 = A12;
    // Let WTL^T be conformal to B1.
    //
    // WTL^T = ( B1 * U11^T + B2 * U12^T ) * inv( triu(T1T)' );
    // WTL   = inv( conj(triu(T1T)) ) * ( U11 * B1^T + U12 * B2^T );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm1( cntl ) );

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, 
                       FLA_ONE, A12, B2, FLA_ONE, WTL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1T, WTL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // B2 = B2 - WTL^T * conj(U12);
    // B1 = B1 - WTL^T * conj(U11);
    //    = B1 - ( U11' * WTL )^T;

    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE,
                       FLA_MINUS_ONE, WTL, A12, FLA_ONE, B2,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_MINUS_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm2( cntl ) );

    FLA_Axpyt_internal( FLA_TRANSPOSE, FLA_ONE, WTL, B1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, /**/ T1, T2,
                              FLA_RIGHT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ B1, B2,
                              FLA_RIGHT );
  }

  return FLA_SUCCESS;
}

