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

FLA_Error FLA_Bsvd_compute_shift( FLA_Obj tol, FLA_Obj sminl, FLA_Obj smax, FLA_Obj d, FLA_Obj e, FLA_Obj shift )
{
	FLA_Datatype datatype;
	int          m_A;
	int          inc_d;
	int          inc_e;

	datatype = FLA_Obj_datatype( d );

	m_A      = FLA_Obj_vector_dim( d );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_tol       = FLA_FLOAT_PTR( tol );
			float*    buff_sminl     = FLA_FLOAT_PTR( sminl );
			float*    buff_smax      = FLA_FLOAT_PTR( smax );
			float*    buff_d         = FLA_FLOAT_PTR( d );
			float*    buff_e         = FLA_FLOAT_PTR( e );
			float*    buff_shift     = FLA_FLOAT_PTR( shift );

			FLA_Bsvd_compute_shift_ops( m_A,
			                            *buff_tol,
			                            *buff_sminl,
			                            *buff_smax,
			                            buff_d, inc_d,
			                            buff_e, inc_e,
			                            buff_shift );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_tol       = FLA_DOUBLE_PTR( tol );
			double*   buff_sminl     = FLA_DOUBLE_PTR( sminl );
			double*   buff_smax      = FLA_DOUBLE_PTR( smax );
			double*   buff_d         = FLA_DOUBLE_PTR( d );
			double*   buff_e         = FLA_DOUBLE_PTR( e );
			double*   buff_shift     = FLA_DOUBLE_PTR( shift );

			FLA_Bsvd_compute_shift_opd( m_A,
			                            *buff_tol,
			                            *buff_sminl,
			                            *buff_smax,
			                            buff_d, inc_d,
			                            buff_e, inc_e,
			                            buff_shift );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_compute_shift_ops( int       m_A,
                                      float     tol,
                                      float     sminl,
                                      float     smax,
                                      float*    buff_d, int inc_d,
                                      float*    buff_e, int inc_e,
                                      float*    shift )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_compute_shift_opd( int       m_A,
                                      double    tol,
                                      double    sminl,
                                      double    smax,
                                      double*   buff_d, int inc_d,
                                      double*   buff_e, int inc_e,
                                      double*   shift )
{
	double  hndrth = 0.01;
	double  eps;
	double* d_first;
	double* e_last;
	double* d_last_m1;
	double* d_last;
	double  sll, temp;

	eps = FLA_Mach_params_opd( FLA_MACH_EPS );

	d_first   = buff_d + (0    )*inc_d;
	e_last    = buff_e + (m_A-2)*inc_e;
	d_last_m1 = buff_d + (m_A-2)*inc_d;
	d_last    = buff_d + (m_A-1)*inc_d;

	// If the shift would ruin relative accuracy, set it to zero.
	if ( m_A * tol * ( sminl / smax ) <= max( eps, hndrth * tol ) )
	{
#ifdef PRINTF
printf( "FLA_Bsvd_compute_shift_opd: shift would ruin accuracy; setting shift to 0.\n" );
printf( "                   m_A = %d     \n", m_A );
printf( "                   tol = %20.15e\n", tol );
printf( "                 sminl = %20.15e\n", sminl );
printf( "                  smax = %20.15e\n", smax );
printf( "                   LHS = %20.15e\n", m_A * tol * ( sminl / smax ) );
printf( "      max(eps,0.01*tol)= %20.15e\n", max( eps, hndrth * tol ) );
#endif
		*shift = 0.0;
	}
	else
	{
		// Compute the shift from the last 2x2 matrix.
		FLA_Sv_2x2_opd( d_last_m1,
		                e_last,
		                d_last,
		                shift,
		                &temp );

		sll = fabs( *d_first );

		// Check if the shift is negligible; if so, set it to zero.
		if ( sll > 0.0 )
		{
			temp = ( *shift / sll );
			if ( temp * temp < eps )
			{
#ifdef PRINTF
printf( "FLA_Bsvd_compute_shift_opd: shift is negligible; setting shift to 0.\n" );
#endif
				*shift = 0.0;
			}
		}
	}

	return FLA_SUCCESS;
}

