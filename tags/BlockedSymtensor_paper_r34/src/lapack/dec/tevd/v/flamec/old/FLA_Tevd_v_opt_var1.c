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

FLA_Error FLA_Tevd_v_opt_var1( FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U )
{
	FLA_Datatype datatype;
	int          m_A, m_U, n_G;
	int          inc_d;
	int          inc_e;
	int          rs_G, cs_G;
	int          rs_U, cs_U;

	datatype = FLA_Obj_datatype( U );

	m_A      = FLA_Obj_vector_dim( d );
	m_U      = FLA_Obj_length( U );
	n_G      = FLA_Obj_width( G );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	
	rs_G     = FLA_Obj_row_stride( G );
	cs_G     = FLA_Obj_col_stride( G );

	rs_U     = FLA_Obj_row_stride( U );
	cs_U     = FLA_Obj_col_stride( U );


	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			float*    buff_U = FLA_FLOAT_PTR( U );

			FLA_Tevd_v_ops_var1( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_U = FLA_DOUBLE_PTR( U );

			FLA_Tevd_v_opd_var1( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			scomplex* buff_U = FLA_COMPLEX_PTR( U );

			FLA_Tevd_v_opc_var1( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );

			FLA_Tevd_v_opz_var1( m_A,
			                     m_U,
			                     n_G,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_ops_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opd_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opc_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opz_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_U, int rs_U, int cs_U )
{
	double  gamma, sigma;
	int     i, k;
	int     k_total = 0;
	int     k_weight = 0;

	for ( i = m_A - 1; i > 1; --i )
	{
		int m_ATL = i + 1;

		/*------------------------------------------------------------*/

		// Find an eigenvalue of the top-left m_ATL-by-m_ATL matrix.
		FLA_Tevd_eigval_v_opd_var1( m_ATL,
		                            n_G,
		                            buff_G, rs_G, cs_G,
		                            buff_d, inc_d,
		                            buff_e, inc_e,
		                            &k );

		k_total += k;
		k_weight += i * k;
//printf( "FLA_Tevd_v_opz_var1: found eig  %14.11f in col %3d after %3d iterations\n", buff_d[ i*inc_d ], i, k );

		// Apply the Givens rotations to update the eigenvectors.
		FLA_Apply_G_rf_opz_var1( k,
		                         m_U,
		                         m_ATL,
		                         buff_G, rs_G, cs_G,
		                         buff_U, rs_U, cs_U );

		/*------------------------------------------------------------*/
	}


//printf( "FLA_Tevd_v_opz_var1: total iter:        %d\n", k_total );
//printf( "FLA_Tevd_v_opz_var1: weighted avg iter: %.3f\n", ( double ) k_weight / ( m_A * m_A / 2 ) );

	// Find the eigenvalue decomposition of the remaining (or only) 2x2
	// submatrix.
	FLA_Hevv_2x2_opd( buff_d + (0  )*inc_d,
	                  buff_e + (0  )*inc_e,
	                  buff_d + (1  )*inc_d,
	                  buff_d + (0  )*inc_d,
	                  buff_d + (1  )*inc_d,
	                  &gamma,
	                  &sigma );

	// Update the eigenvectors.
	FLA_Apply_G_mx2_opz( m_U,
	                     &gamma,
	                     &sigma,
	                     buff_U + (0  )*cs_U, rs_U,
	                     buff_U + (1  )*cs_U, rs_U );

	return FLA_SUCCESS;
}
