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

FLA_Error FLA_Tridiag_UT_shift_U( FLA_Uplo uplo, FLA_Obj A )
{
	FLA_Datatype datatype;
	int          m_A;
	int          rs_A, cs_A;

	datatype = FLA_Obj_datatype( A );
	m_A      = FLA_Obj_length( A );
	rs_A     = FLA_Obj_row_stride( A );
	cs_A     = FLA_Obj_col_stride( A );

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Tridiag_UT_shift_U_check( uplo, A );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_A = ( float* ) FLA_FLOAT_PTR( A );

			if ( uplo == FLA_LOWER_TRIANGULAR )
				FLA_Tridiag_UT_shift_U_l_ops( m_A,
				                              buff_A, rs_A, cs_A );
			//else // if ( uplo == FLA_UPPER_TRIANGULAR )
			//	FLA_Tridiag_UT_shift_U_u_ops( m_A,
			//	                              buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_A = ( double* ) FLA_DOUBLE_PTR( A );

			if ( uplo == FLA_LOWER_TRIANGULAR )
				FLA_Tridiag_UT_shift_U_l_opd( m_A,
				                              buff_A, rs_A, cs_A );
			//else // if ( uplo == FLA_UPPER_TRIANGULAR )
			//	FLA_Tridiag_UT_shift_U_u_opd( m_A,
			//	                              buff_A, rs_A, cs_A );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

			if ( uplo == FLA_LOWER_TRIANGULAR )
				FLA_Tridiag_UT_shift_U_l_opc( m_A,
				                              buff_A, rs_A, cs_A );
			//else // if ( uplo == FLA_UPPER_TRIANGULAR )
			//	FLA_Tridiag_UT_shift_U_u_opc( m_A,
			//	                              buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

			if ( uplo == FLA_LOWER_TRIANGULAR )
				FLA_Tridiag_UT_shift_U_l_opz( m_A,
				                              buff_A, rs_A, cs_A );
			//else // if ( uplo == FLA_UPPER_TRIANGULAR )
			//	FLA_Tridiag_UT_shift_U_u_opz( m_A,
			//	                              buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_shift_U_l_ops( int       m_A,
                                        float*    buff_A, int rs_A, int cs_A )
{
	float*  a00  = buff_A;
	float*  a10  = buff_A + rs_A;
	float   zero = bli_s0();
	float   one  = bli_s1();
	int     j;

	for ( j = m_A - 1; j > 0; --j )
	{
		float*    alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		float*    alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		float*    a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		float*    a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bli_scopyv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bli_ssetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_shift_U_l_opd( int       m_A,
                                        double*   buff_A, int rs_A, int cs_A )
{
	double* a00  = buff_A;
	double* a10  = buff_A + rs_A;
	double  zero = bli_d0();
	double  one  = bli_d1();
	int     j;

	for ( j = m_A - 1; j > 0; --j )
	{
		double*   alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		double*   alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		double*   a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		double*   a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bli_dcopyv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bli_dsetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_shift_U_l_opc( int       m_A,
                                        scomplex* buff_A, int rs_A, int cs_A )
{
	scomplex* a00  = buff_A;
	scomplex* a10  = buff_A + rs_A;
	scomplex  zero = bli_c0();
	scomplex  one  = bli_c1();
	int       j;

	for ( j = m_A - 1; j > 0; --j )
	{
		scomplex* alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		scomplex* alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		scomplex* a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		scomplex* a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bli_ccopyv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bli_csetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_shift_U_l_opz( int       m_A,
                                        dcomplex* buff_A, int rs_A, int cs_A )
{
	dcomplex* a00  = buff_A;
	dcomplex* a10  = buff_A + rs_A;
	dcomplex  zero = bli_z0();
	dcomplex  one  = bli_z1();
	int       j;

	for ( j = m_A - 1; j > 0; --j )
	{
		dcomplex* alpha01  = buff_A + (j  )*cs_A + (0  )*rs_A;
		dcomplex* alpha11  = buff_A + (j  )*cs_A + (j  )*rs_A;
		dcomplex* a20      = buff_A + (j-1)*cs_A + (j+1)*rs_A;
		dcomplex* a21      = buff_A + (j  )*cs_A + (j+1)*rs_A;

		int       m_ahead  = m_A - j - 1;

		*alpha01 = zero;
		*alpha11 = one;
		bli_zcopyv( BLIS_NO_CONJUGATE,
		            m_ahead,
		            a20, rs_A,
		            a21, rs_A );
	}

	*a00 = one;
	bli_zsetv( m_A - 1,
	           &zero,
	           a10, rs_A );

	return FLA_SUCCESS;
}

