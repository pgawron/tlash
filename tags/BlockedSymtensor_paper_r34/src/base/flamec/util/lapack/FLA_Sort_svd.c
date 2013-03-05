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

FLA_Error FLA_Sort_svd( FLA_Direct direct, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
	FLA_Datatype datatype;
	dim_t        m_U, n_V;
	dim_t        rs_U, cs_U;
	dim_t        rs_V, cs_V;
	dim_t        inc_s;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Sort_svd_check( direct, s, U, V );

	datatype = FLA_Obj_datatype( U );

	m_U      = FLA_Obj_length( U );
	n_V      = FLA_Obj_length( V );

	rs_U     = FLA_Obj_row_stride( U );
	cs_U     = FLA_Obj_col_stride( U );

	rs_V     = FLA_Obj_row_stride( V );
	cs_V     = FLA_Obj_col_stride( V );

	inc_s    = FLA_Obj_vector_inc( s );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float* s_p = ( float* ) FLA_FLOAT_PTR( s );
			float* U_p = ( float* ) FLA_FLOAT_PTR( U );
			float* V_p = ( float* ) FLA_FLOAT_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_svd_f_ops( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_svd_b_ops( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );

			break;
		}

		case FLA_DOUBLE:
		{
			double* s_p = ( double* ) FLA_DOUBLE_PTR( s );
			double* U_p = ( double* ) FLA_DOUBLE_PTR( U );
			double* V_p = ( double* ) FLA_DOUBLE_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_svd_f_opd( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_svd_b_opd( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    s_p = ( float*    ) FLA_FLOAT_PTR( s );
			scomplex* U_p = ( scomplex* ) FLA_COMPLEX_PTR( U );
			scomplex* V_p = ( scomplex* ) FLA_COMPLEX_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_svd_f_opc( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_svd_b_opc( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   s_p = ( double*   ) FLA_DOUBLE_PTR( s );
			dcomplex* U_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( U );
			dcomplex* V_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_svd_f_opz( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_svd_b_opz( m_U,
				                    n_V,
				                    s_p, inc_s,
				                    U_p, rs_U, cs_U,
				                    V_p, rs_V, cs_V );

			break;
		}

	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Sort_svd_f_ops( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              float*    U, int rs_U, int cs_U,
                              float*    V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_svd_b_ops( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              float*    U, int rs_U, int cs_U,
                              float*    V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_svd_f_opd( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              double*   U, int rs_U, int cs_U,
                              double*   V, int rs_V, int cs_V )
{
	int    min_m_n = min( m_U, n_V );
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < min_m_n; ++ii )
	{
		i = ii - 1;
		k = i;

		p = s[ i*inc_s ];

		for ( j = ii; j < min_m_n; ++j )
		{
			if ( s[ j*inc_s ] < p )
			{
				k = j;
				p = s[ j*inc_s ];
			}
		}

		if ( k != i )
		{
			s[ k*inc_s ] = s[ i ];
			s[ i       ] = p;
			bli_dswapv( m_U,
			            U + i*cs_U, rs_U,
			            U + k*cs_U, rs_U );
			bli_dswapv( n_V,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_svd_b_opd( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              double*   U, int rs_U, int cs_U,
                              double*   V, int rs_V, int cs_V )
{
	int    min_m_n = min( m_U, n_V );
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < min_m_n; ++ii )
	{
		i = ii - 1;
		k = i;

		p = s[ i*inc_s ];

		for ( j = ii; j < min_m_n; ++j )
		{
			if ( s[ j*inc_s ] > p )
			{
				k = j;
				p = s[ j*inc_s ];
			}
		}

		if ( k != i )
		{
			s[ k*inc_s ] = s[ i ];
			s[ i       ] = p;
			bli_dswapv( m_U,
			            U + i*cs_U, rs_U,
			            U + k*cs_U, rs_U );
			bli_dswapv( n_V,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_svd_f_opc( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              scomplex* U, int rs_U, int cs_U,
                              scomplex* V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_svd_b_opc( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              scomplex* U, int rs_U, int cs_U,
                              scomplex* V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_svd_f_opz( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              dcomplex* U, int rs_U, int cs_U,
                              dcomplex* V, int rs_V, int cs_V )
{
	int    min_m_n = min( m_U, n_V );
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < min_m_n; ++ii )
	{
		i = ii - 1;
		k = i;

		p = s[ i*inc_s ];

		for ( j = ii; j < min_m_n; ++j )
		{
			if ( s[ j*inc_s ] < p )
			{
				k = j;
				p = s[ j*inc_s ];
			}
		}

		if ( k != i )
		{
			s[ k*inc_s ] = s[ i ];
			s[ i       ] = p;
			bli_zswapv( m_U,
			            U + i*cs_U, rs_U,
			            U + k*cs_U, rs_U );
			bli_zswapv( n_V,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_svd_b_opz( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              dcomplex* U, int rs_U, int cs_U,
                              dcomplex* V, int rs_V, int cs_V )
{
	int    min_m_n = min( m_U, n_V );
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < min_m_n; ++ii )
	{
		i = ii - 1;
		k = i;

		p = s[ i*inc_s ];

		for ( j = ii; j < min_m_n; ++j )
		{
			if ( s[ j*inc_s ] > p )
			{
				k = j;
				p = s[ j*inc_s ];
			}
		}

		if ( k != i )
		{
			s[ k*inc_s ] = s[ i ];
			s[ i       ] = p;
			bli_zswapv( m_U,
			            U + i*cs_U, rs_U,
			            U + k*cs_U, rs_U );
			bli_zswapv( n_V,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}
