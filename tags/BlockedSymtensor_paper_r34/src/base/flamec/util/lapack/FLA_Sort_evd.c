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

FLA_Error FLA_Sort_evd( FLA_Direct direct, FLA_Obj l, FLA_Obj V )
{
	FLA_Datatype datatype;
	dim_t        m_A;
	dim_t        rs_V, cs_V;
	dim_t        inc_l;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Sort_evd_check( direct, l, V );

	datatype = FLA_Obj_datatype( V );

	m_A      = FLA_Obj_length( V );

	rs_V     = FLA_Obj_row_stride( V );
	cs_V     = FLA_Obj_col_stride( V );

	inc_l    = FLA_Obj_vector_inc( l );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float* l_p = ( float* ) FLA_FLOAT_PTR( l );
			float* V_p = ( float* ) FLA_FLOAT_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_evd_f_ops( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_evd_b_ops( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );

			break;
		}

		case FLA_DOUBLE:
		{
			double* l_p = ( double* ) FLA_DOUBLE_PTR( l );
			double* V_p = ( double* ) FLA_DOUBLE_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_evd_f_opd( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_evd_b_opd( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    l_p = ( float*    ) FLA_FLOAT_PTR( l );
			scomplex* V_p = ( scomplex* ) FLA_COMPLEX_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_evd_f_opc( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_evd_b_opc( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   l_p = ( double*   ) FLA_DOUBLE_PTR( l );
			dcomplex* V_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( V );

			if ( direct == FLA_FORWARD )
				FLA_Sort_evd_f_opz( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_evd_b_opz( m_A,
				                    l_p, inc_l,
				                    V_p, rs_V, cs_V );

			break;
		}

	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Sort_evd_f_ops( int       m_A,
                              float*    l, int inc_l,
                              float*    V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_evd_b_ops( int       m_A,
                              float*    l, int inc_l,
                              float*    V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_evd_f_opd( int       m_A,
                              double*   l, int inc_l,
                              double*   V, int rs_V, int cs_V )
{
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < m_A; ++ii )
	{
		i = ii - 1;
		k = i;

		p = l[ i*inc_l ];

		for ( j = ii; j < m_A; ++j )
		{
			if ( l[ j*inc_l ] < p )
			{
				k = j;
				p = l[ j*inc_l ];
			}
		}

		if ( k != i )
		{
			l[ k*inc_l ] = l[ i ];
			l[ i       ] = p;
			bli_dswapv( m_A,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_evd_b_opd( int       m_A,
                              double*   l, int inc_l,
                              double*   V, int rs_V, int cs_V )
{
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < m_A; ++ii )
	{
		i = ii - 1;
		k = i;

		p = l[ i*inc_l ];

		for ( j = ii; j < m_A; ++j )
		{
			if ( l[ j*inc_l ] > p )
			{
				k = j;
				p = l[ j*inc_l ];
			}
		}

		if ( k != i )
		{
			l[ k*inc_l ] = l[ i ];
			l[ i       ] = p;
			bli_dswapv( m_A,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_evd_f_opc( int       m_A,
                              float*    l, int inc_l,
                              scomplex* V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_evd_b_opc( int       m_A,
                              float*    l, int inc_l,
                              scomplex* V, int rs_V, int cs_V )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_evd_f_opz( int       m_A,
                              double*   l, int inc_l,
                              dcomplex* V, int rs_V, int cs_V )
{
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < m_A; ++ii )
	{
		i = ii - 1;
		k = i;

		p = l[ i*inc_l ];

		for ( j = ii; j < m_A; ++j )
		{
			if ( l[ j*inc_l ] < p )
			{
				k = j;
				p = l[ j*inc_l ];
			}
		}

		if ( k != i )
		{
			l[ k*inc_l ] = l[ i ];
			l[ i       ] = p;
			bli_zswapv( m_A,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_evd_b_opz( int       m_A,
                              double*   l, int inc_l,
                              dcomplex* V, int rs_V, int cs_V )
{
	int    i, ii, j, k;
	double p;

	for ( ii = 1; ii < m_A; ++ii )
	{
		i = ii - 1;
		k = i;

		p = l[ i*inc_l ];

		for ( j = ii; j < m_A; ++j )
		{
			if ( l[ j*inc_l ] > p )
			{
				k = j;
				p = l[ j*inc_l ];
			}
		}

		if ( k != i )
		{
			l[ k*inc_l ] = l[ i ];
			l[ i       ] = p;
			bli_zswapv( m_A,
			            V + i*cs_V, rs_V,
			            V + k*cs_V, rs_V );
		}
	}

	return FLA_SUCCESS;
}
