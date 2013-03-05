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

FLA_Error FLA_QR_UT_form_Q( FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
	FLA_Error r_val = FLA_SUCCESS;
	FLA_Obj   QTL, QTR,
	          QBL, QBR;
	FLA_Obj   QL;
	FLA_Obj   W;
	dim_t     n_A;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_QR_UT_form_Q_check( A, T, Q );

	// Zero out the upper triangle of Q.
	FLA_Setr( FLA_UPPER_TRIANGULAR, FLA_ZERO, Q );

	// If A and Q are different objects, copy the Householder vectors
	// from A to QT, and zero out the lower triangle of QBR. If they
	// are the same object, we don't need to do the copy, and don't
	// need to zero anything out since the user should only have A and
	// Q be the same object if A is square, since Q needs to be square
	// (specifically, dim(Q) needs to equal m(A)).
	if ( FLA_Obj_is( A, Q ) == FALSE )
	{
		n_A = FLA_Obj_width( A );
		FLA_Part_2x2( Q,   &QTL, &QTR,
		                   &QBL, &QBR,   n_A, n_A, FLA_TL );
		FLA_Merge_2x1( QTL,
		               QBL,  &QL );

		// Copy the Householder vectors in A to QL.
		FLA_Copyr( FLA_LOWER_TRIANGULAR, A, QL );

		// Zero out the lower triangle of QBR.
		FLA_Setr( FLA_LOWER_TRIANGULAR, FLA_ZERO, QBR );
	}

	// Set the digaonal to one.
	FLA_Set_diag( FLA_ONE, Q );

	// Create workspace for applying the block Householder transforms.
	FLA_Apply_Q_UT_create_workspace( T, Q, &W );

	// Overwrite Q, which currently contains Householder vectors in the
	// strictly lower triangle and identity in the upper triangle, with
	// the unitary matrix associated with those Householder transforms.
	r_val = FLA_QR_UT_form_Q_blk_var1( Q, T, W );
	
	// Free the temporary workspace.
	FLA_Obj_free( &W );

/*
FLA_Apply_Q_UT_create_workspace( T, Q, &W );
FLA_Set_to_identity( Q );
FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                A, T, W, Q );
FLA_Obj_free( &W );
FLA_Obj_show( "Q", Q, "%8.1e %8.1e ", "" );
*/

	return r_val;
}

FLA_Error FLA_QR_UT_form_Q_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj T1T,
          T2B;

  FLA_Obj WTL,  WTR,
          WBL,  WBR;

  FLA_Obj AB1,   AB2;

  dim_t   b, b_alg;
  dim_t   m_BR, n_BR;

  b_alg = FLA_Obj_length( T );

  // If A is wider than T, then we need to position ourseves carefully
  // within the matrix for the initial partitioning.
  if ( FLA_Obj_width( A ) > FLA_Obj_width( T ) )
  {
    m_BR = FLA_Obj_length( A ) - FLA_Obj_width( T );
    n_BR = FLA_Obj_width( A ) - FLA_Obj_width( T );
  }
  else
  {
    m_BR = 0;
    n_BR = 0;
  }
//printf( "\nn(A) n(T) = %d %d\n", FLA_Obj_width( A ), FLA_Obj_width( T ) );
//printf(   "m_BR n_BR = %d %d\n", m_BR, n_BR );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     m_BR, n_BR, FLA_BR );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_RIGHT );

  while ( FLA_Obj_length( ATL ) > 0 )
  {
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

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T,
                         &T2B,     b, FLA_TOP );

    FLA_Part_2x2( W,     &WTL, &WTR,
                         &WBL, &WBR,     b, FLA_Obj_width( A12 ), FLA_TL );

    // Use an unblocked algorithm for the first (or only) block.
    if ( FLA_Obj_length( ABR ) == 0 )
    {
      FLA_QR_UT_form_Q_opt_var1( A11, T1T );
    }
    else
    {
      FLA_Merge_2x1( A11,
                     A21,   &AB1 );
      FLA_Merge_2x1( A12,
                     A22,   &AB2 );

      // Apply the block Householder transforms to A12 and A22.
      FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                      AB1, T1T, WTL, AB2 );

      // Apply H to the current block panel consisting of A11 and A21.
      FLA_QR_UT_form_Q_opt_var1( AB1, T1T );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, /**/ T1, T2,
                              FLA_RIGHT );
  }

  return FLA_SUCCESS;
}


FLA_Error FLA_QR_UT_form_Q_opt_var1( FLA_Obj A, FLA_Obj T )
{
	FLA_Datatype datatype;
	int          m_A, n_A;
	int          rs_A, cs_A;
	int          rs_T, cs_T;

	datatype = FLA_Obj_datatype( A );

	m_A      = FLA_Obj_length( A );
	n_A      = FLA_Obj_width( A );
	rs_A     = FLA_Obj_row_stride( A );
	cs_A     = FLA_Obj_col_stride( A );

	rs_T     = FLA_Obj_row_stride( T );
	cs_T     = FLA_Obj_col_stride( T );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_A = ( float* ) FLA_FLOAT_PTR( A );
			float*    buff_T = ( float* ) FLA_FLOAT_PTR( T );

			FLA_QR_UT_form_Q_ops_var1( m_A,
			                           n_A,
			                           buff_A, rs_A, cs_A,
			                           buff_T, rs_T, cs_T );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_A = ( double* ) FLA_DOUBLE_PTR( A );
			double*   buff_T = ( double* ) FLA_DOUBLE_PTR( T );

			FLA_QR_UT_form_Q_opd_var1( m_A,
			                           n_A,
			                           buff_A, rs_A, cs_A,
			                           buff_T, rs_T, cs_T );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );
			scomplex* buff_T = ( scomplex* ) FLA_COMPLEX_PTR( T );

			FLA_QR_UT_form_Q_opc_var1( m_A,
			                           n_A,
			                           buff_A, rs_A, cs_A,
			                           buff_T, rs_T, cs_T );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
			dcomplex* buff_T = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( T );

			FLA_QR_UT_form_Q_opz_var1( m_A,
			                           n_A,
			                           buff_A, rs_A, cs_A,
			                           buff_T, rs_T, cs_T );

			break;
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_QR_UT_form_Q_ops_var1( int       m_A,
                                     int       n_A,
                                     float*    buff_A, int rs_A, int cs_A,
                                     float*    buff_T, int rs_T, int cs_T )
{
	float    one     = bli_d1();
	int      min_m_n = min( m_A, n_A );
	int      i;

	for ( i = min_m_n - 1; i >= 0; --i )
	{
		//float*    a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
		float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
		float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
		float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
		float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

		float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

		float     minus_inv_tau11;

		//int       m_behind = i;
		int       n_ahead  = n_A - i - 1;
		int       m_ahead  = m_A - i - 1;

		FLA_Apply_H2_UT_l_ops_var1( m_ahead,
		                            n_ahead,
		                            tau11,
		                            a21,  rs_A,
		                            a12t, cs_A,
		                            A22,  rs_A, cs_A );

		minus_inv_tau11 = -one / *tau11;

		*alpha11 = one + minus_inv_tau11;

		bli_sscalv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            &minus_inv_tau11,
		            a21, rs_A );

		// Not necessary if upper triangle of A is initialized to identity.
		//bli_ssetv( m_behind,
		//           &zero,
		//           a01, rs_A );
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_QR_UT_form_Q_opd_var1( int       m_A,
                                     int       n_A,
                                     double*   buff_A, int rs_A, int cs_A,
                                     double*   buff_T, int rs_T, int cs_T )
{
	double   one     = bli_d1();
	int      min_m_n = min( m_A, n_A );
	int      i;

	for ( i = min_m_n - 1; i >= 0; --i )
	{
		//double*   a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
		double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
		double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
		double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
		double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

		double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

		double    minus_inv_tau11;

		//int       m_behind = i;
		int       n_ahead  = n_A - i - 1;
		int       m_ahead  = m_A - i - 1;

		FLA_Apply_H2_UT_l_opd_var1( m_ahead,
		                            n_ahead,
		                            tau11,
		                            a21,  rs_A,
		                            a12t, cs_A,
		                            A22,  rs_A, cs_A );

		minus_inv_tau11 = -one / *tau11;

		*alpha11 = one + minus_inv_tau11;

		bli_dscalv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            &minus_inv_tau11,
		            a21, rs_A );

		// Not necessary if upper triangle of A is initialized to identity.
		//bli_dsetv( m_behind,
		//           &zero,
		//           a01, rs_A );
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_QR_UT_form_Q_opc_var1( int       m_A,
                                     int       n_A,
                                     scomplex* buff_A, int rs_A, int cs_A,
                                     scomplex* buff_T, int rs_T, int cs_T )
{
	scomplex zero    = bli_c0();
	scomplex one     = bli_c1();
	int      min_m_n = min( m_A, n_A );
	int      i;

	for ( i = min_m_n - 1; i >= 0; --i )
	{
		//scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
		scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
		scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
		scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
		scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

		scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

		scomplex  minus_inv_tau11;

		//int       m_behind = i;
		int       n_ahead  = n_A - i - 1;
		int       m_ahead  = m_A - i - 1;

		FLA_Apply_H2_UT_l_opc_var1( m_ahead,
		                            n_ahead,
		                            tau11,
		                            a21,  rs_A,
		                            a12t, cs_A,
		                            A22,  rs_A, cs_A );

		minus_inv_tau11.real = -one.real / tau11->real;
		minus_inv_tau11.imag = zero.imag;

		alpha11->real = one.real + minus_inv_tau11.real;
		alpha11->imag = zero.imag;

		bli_cscalv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            &minus_inv_tau11,
		            a21, rs_A );

		// Not necessary if upper triangle of A is initialized to identity.
		//bli_csetv( m_behind,
		//           &zero,
		//           a01, rs_A );
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_QR_UT_form_Q_opz_var1( int       m_A,
                                     int       n_A,
                                     dcomplex* buff_A, int rs_A, int cs_A,
                                     dcomplex* buff_T, int rs_T, int cs_T )
{
	dcomplex zero    = bli_z0();
	dcomplex one     = bli_z1();
	int      min_m_n = min( m_A, n_A );
	int      i;

	for ( i = min_m_n - 1; i >= 0; --i )
	{
		//dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
		dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
		dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
		dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
		dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

		dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

		dcomplex  minus_inv_tau11;

		//int       m_behind = i;
		int       n_ahead  = n_A - i - 1;
		int       m_ahead  = m_A - i - 1;

		FLA_Apply_H2_UT_l_opz_var1( m_ahead,
		                            n_ahead,
		                            tau11,
		                            a21,  rs_A,
		                            a12t, cs_A,
		                            A22,  rs_A, cs_A );

		minus_inv_tau11.real = -one.real / tau11->real;
		minus_inv_tau11.imag = zero.imag;

		alpha11->real = one.real + minus_inv_tau11.real;
		alpha11->imag = zero.imag;

		bli_zscalv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            &minus_inv_tau11,
		            a21, rs_A );

		// Not necessary if upper triangle of A is initialized to identity.
		//bli_zsetv( m_behind,
		//           &zero,
		//           a01, rs_A );
	}

	return FLA_SUCCESS;
}

