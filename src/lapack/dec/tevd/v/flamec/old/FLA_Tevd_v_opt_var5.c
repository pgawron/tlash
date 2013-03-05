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

FLA_Error FLA_Tevd_v_opt_var5( FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U )
{
	FLA_Datatype datatype;
	int          m_A, m_U;
	int          inc_d;
	int          inc_e;
	int          rs_G, cs_G;
	int          rs_U, cs_U;

	datatype = FLA_Obj_datatype( U );

	m_A      = FLA_Obj_vector_dim( d );
	m_U      = FLA_Obj_length( U );

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
			float*    buff_G = FLA_FLOAT_PTR( G );
			float*    buff_U = FLA_FLOAT_PTR( U );

			FLA_Tevd_v_ops_var5( m_A,
			                     m_U,
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
			double*   buff_G = FLA_DOUBLE_PTR( G );
			double*   buff_U = FLA_DOUBLE_PTR( U );

			FLA_Tevd_v_opd_var5( m_A,
			                     m_U,
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
			float*    buff_G = FLA_FLOAT_PTR( G );
			scomplex* buff_U = FLA_COMPLEX_PTR( U );

			FLA_Tevd_v_opc_var5( m_A,
			                     m_U,
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
			double*   buff_G = FLA_DOUBLE_PTR( G );
			dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );

			FLA_Tevd_v_opz_var5( m_A,
			                     m_U,
			                     buff_d, inc_d,
			                     buff_e, inc_e,
			                     buff_G, rs_G, cs_G,
			                     buff_U, rs_U, cs_U );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_ops_var5( int       m_A,
                               int       m_U,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_G, int rs_G, int cs_G,
                               float*    buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_opd_var5( int       m_A,
                               int       m_U,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_G, int rs_G, int cs_G,
                               double*   buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opc_var5( int       m_A,
                               int       m_U,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_G, int rs_G, int cs_G,
                               scomplex* buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_v_opz_var5( int       m_A,
                               int       m_U,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_G, int rs_G, int cs_G,
                               dcomplex* buff_U, int rs_U, int cs_U )
{
	double  one  = bli_d1();
	double  zero = bli_d0();
	double* buff_Ur;
	double* buff_Uc;
	int     rs_Ur, cs_Ur;
	int     rs_Uc, cs_Uc;
	int     i;
	dgiv_t  rots;
	
/*
	double* buff_d_copy;
	double* buff_e_copy;
	int     n_iter;

	buff_d_copy = bli_dallocv( m_A );
	buff_e_copy = bli_dallocv( m_A - 1 );
	bli_dcopyv( BLIS_NO_CONJUGATE,
	            m_A,
	            buff_d,      inc_d,
	            buff_d_copy, 1 );
	bli_dcopyv( BLIS_NO_CONJUGATE,
	            m_A - 1,
	            buff_e,      inc_e,
	            buff_e_copy, 1 );
	n_iter = FLA_Tevd_n_opd_var1( m_A,
	                              buff_d_copy, 1,
	                              buff_e_copy, 1 );
*/

//printf( "num iterations: %d\n", n_iter );

/*
typedef struct fla_dgivens
{
	int        n_g;
	int        n_g_alloc;
	dcomplex** g;
	int*       j_off;
	int*       m_g;
	
} dgiv_t;
*/
	// Allocate and initialize the Givens structure.
	rots.n_g       = 0;
	rots.n_g_alloc = 20 * m_A;
	rots.g         = bli_vallocv( rots.n_g_alloc, sizeof( dcomplex* ) );
	rots.j_off     = bli_iallocv( rots.n_g_alloc );
	rots.m_g       = bli_iallocv( rots.n_g_alloc );

// Goal: to save up all k_total iterations worth of C and S and apply them to
// G all at once, and then apply G to U.

	// Find the eigenvalues of the tridiagonal matrix and accumulate the
	// Givens rotations in G.
	FLA_Tevd_v_opd_var5r( m_A,
	                      0,
	                      buff_d, inc_d,
	                      buff_e, inc_e,
	                      &rots );

	FLA_Apply_G_rf_opd_var5( m_U,
	                         &rots,
	                         buff_G, rs_G, cs_G );



	// Free the Givens structure and all of its elements.
	for ( i = 0; i < rots.n_g; ++i )
		bli_zfree( rots.g[i] );
	bli_ifree( rots.m_g );
	bli_ifree( rots.j_off );
	bli_vfree( rots.g );

	// Allocate temporary matrices for real and imaginary parts of U.
	buff_Ur = bli_dallocm( m_U, m_U );
	rs_Ur   = 1;
	cs_Ur   = m_U;
	buff_Uc = bli_dallocm( m_U, m_U );
	rs_Uc   = 1;
	cs_Uc   = m_U;

	// Copy the real and imaginary parts of U into separate contiguous
	// matrices.
	bli_dcopymt( BLIS_NO_TRANSPOSE,
	             m_U,
	             m_U,
	             ((double*)(buff_U))+0, rs_U*2, cs_U*2,
	             buff_Ur,               rs_Ur,  cs_Ur );
	bli_dcopymt( BLIS_NO_TRANSPOSE,
	             m_U,
	             m_U,
	             ((double*)(buff_U))+1, rs_U*2, cs_U*2,
	             buff_Uc,               rs_Uc,  cs_Uc );

	// Apply G to U by computing with the real and imaginary parts separately.
	bli_dgemm( BLIS_NO_TRANSPOSE,
	           BLIS_NO_TRANSPOSE, 
	           m_U,
	           m_U, 
	           m_U, 
	           &one,
	           buff_Ur, rs_Ur, cs_Ur,
	           buff_G,  rs_G,  cs_G,
	           &zero,
	           ((double*)(buff_U))+0, rs_U*2, cs_U*2 );
	bli_dgemm( BLIS_NO_TRANSPOSE,
	           BLIS_NO_TRANSPOSE, 
	           m_U,
	           m_U, 
	           m_U, 
	           &one,
	           buff_Uc, rs_Uc, cs_Uc,
	           buff_G,  rs_G,  cs_G,
	           &zero,
	           ((double*)(buff_U))+1, rs_U*2, cs_U*2 );

	// Free temporary matrices.
	bli_dfree( buff_Ur );
	bli_dfree( buff_Uc );
/*
	bli_dfree( buff_d_copy );
	bli_dfree( buff_e_copy );
*/

	return FLA_SUCCESS;
}

