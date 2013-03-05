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

FLA_Error FLA_Bsvd_sinval_v_opt_var1( FLA_Obj tol, FLA_Obj thresh, FLA_Obj G, FLA_Obj H, FLA_Obj d, FLA_Obj e, FLA_Obj k )
{
	FLA_Datatype datatype;
	int          m_A, n_GH;
	int          rs_G, cs_G;
	int          rs_H, cs_H;
	int          inc_d;
	int          inc_e;

	datatype = FLA_Obj_datatype( d );

	m_A      = FLA_Obj_vector_dim( d );
	n_GH     = FLA_Obj_width( G );

	rs_G     = FLA_Obj_row_stride( G );
	cs_G     = FLA_Obj_col_stride( G );

	rs_H     = FLA_Obj_row_stride( H );
	cs_H     = FLA_Obj_col_stride( H );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_tol    = FLA_FLOAT_PTR( tol );
			float*    buff_thresh = FLA_FLOAT_PTR( thresh );
			scomplex* buff_G      = FLA_COMPLEX_PTR( G );
			scomplex* buff_H      = FLA_COMPLEX_PTR( H );
			float*    buff_d      = FLA_FLOAT_PTR( d );
			float*    buff_e      = FLA_FLOAT_PTR( e );
			int*      buff_k      = FLA_INT_PTR( k );

			FLA_Bsvd_sinval_v_ops_var1( m_A,
			                            n_GH,
			                            9,
			                            *buff_tol,
			                            *buff_thresh,
			                            buff_G, rs_G, cs_G,
			                            buff_H, rs_H, cs_H,
			                            buff_d, inc_d,
			                            buff_e, inc_e,
			                            buff_k );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_tol    = FLA_DOUBLE_PTR( tol );
			double*   buff_thresh = FLA_DOUBLE_PTR( thresh );
			dcomplex* buff_G      = FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_H      = FLA_DOUBLE_COMPLEX_PTR( H );
			double*   buff_d      = FLA_DOUBLE_PTR( d );
			double*   buff_e      = FLA_DOUBLE_PTR( e );
			int*      buff_k      = FLA_INT_PTR( k );

			FLA_Bsvd_sinval_v_opd_var1( m_A,
			                            n_GH,
			                            9,
			                            *buff_tol,
			                            *buff_thresh,
			                            buff_G, rs_G, cs_G,
			                            buff_H, rs_H, cs_H,
			                            buff_d, inc_d,
			                            buff_e, inc_e,
			                            buff_k );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_sinval_v_ops_var1( int       m_A,
                                      int       n_GH,
                                      int       n_iter_allowed,
                                      float     tol,
                                      float     thresh,
                                      scomplex* buff_G, int rs_G, int cs_G,
                                      scomplex* buff_H, int rs_H, int cs_H,
                                      float*    buff_d, int inc_d, 
                                      float*    buff_e, int inc_e,
                                      int*      n_iter )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Bsvd_sinval_v_opd_var1( int       m_A,
                                      int       n_GH,
                                      int       n_iter_allowed,
                                      double    tol,
                                      double    thresh,
                                      dcomplex* buff_G, int rs_G, int cs_G,
                                      dcomplex* buff_H, int rs_H, int cs_H,
                                      double*   buff_d, int inc_d, 
                                      double*   buff_e, int inc_e,
                                      int*      n_iter )
{
	FLA_Error r_val;
	double    one = bli_d1();
	//double*   d_first;
	//double*   d_last_m1;
	double*   e_last;
	//double*   d_last;
	double    smax;
	double    smin;
	double    sminl;
	double    shift;
	int       k;

	// Initialize pointers to some diagonal and superdiagonal elements
	// that we will refer to later.
	e_last    = buff_e + (m_A-2)*inc_e;
	//d_last_m1 = buff_d + (m_A-2)*inc_d;
	//d_last    = buff_d + (m_A-1)*inc_d;
	//d_first   = buff_d + (0    )*inc_d;

	// Find the largest element of the diagonal or superdiagonal.
	// This is used later when checking the shift.
	FLA_Bsvd_find_max_min_opd( m_A,
	                           buff_d, inc_d,
	                           buff_e, inc_e,
	                           &smax,
	                           &smin );

	// Perform some iterations.
	for ( k = 0; k < n_iter_allowed; ++k )
	{
		dcomplex* g1 = buff_G + (k  )*cs_G;
		dcomplex* h1 = buff_H + (k  )*cs_H;

		/*------------------------------------------------------------*/

		// Before we perform any rotations, check for pre-existing deflation.
		r_val = FLA_Bsvd_find_converged_opd( m_A,
		                                     tol,
		                                     buff_d, inc_d,
		                                     buff_e, inc_e,
		                                     &sminl );

		// If r_val is positive, then deflation was found.
		if ( 0 <= r_val )
		{
#ifdef PRINTF
			printf( "FLA_Bsvd_sinval_v_opt_var1: Deflation detected in col %d, sval %d\n", r_val, m_A - 1 );
			printf( "FLA_Bsvd_sinval_v_opt_var1: alpha11 alpha12 = %23.19e %23.19e\n", buff_d[r_val*inc_d], buff_e[r_val*inc_e] );
			printf( "FLA_Bsvd_sinval_v_opt_var1:         alpha22 =         %43.19e\n", buff_d[(r_val+1)*inc_d] );
#endif

			// Set the off-diagonal element to zero.
			buff_e[ (r_val)*inc_e ] = 0.0;

			*n_iter = k;
			return r_val;
		}


		// Compute a shift with the last 2x2 matrix.
		FLA_Bsvd_compute_shift_opd( m_A,
		                            tol,
		                            sminl,
		                            smax,
		                            buff_d, inc_d,
		                            buff_e, inc_e,
		                            &shift );

		// Perform a Francis step.
		r_val = FLA_Bsvd_francis_v_opd_var1( m_A,
		                                     shift,
		                                     g1,     rs_G,
		                                     h1,     rs_H,
		                                     buff_d, inc_d,
		                                     buff_e, inc_e );

		// Check for convergence using thresh.
		if ( MAC_Bsvd_sinval_is_converged_opd( thresh, one, *e_last ) )
		{
			*e_last = 0.0;
			*n_iter = k + 1;
			return m_A - 1;
		}

		/*------------------------------------------------------------*/
	}

	*n_iter = n_iter_allowed;
	return FLA_FAILURE;
}

