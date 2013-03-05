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

FLA_Error FLA_Bsvd_find_max( FLA_Obj d, FLA_Obj e, FLA_Obj smax, FLA_Obj smin )
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
			float*    buff_d    = FLA_FLOAT_PTR( d );
			float*    buff_e    = FLA_FLOAT_PTR( e );
			float*    buff_smax = FLA_FLOAT_PTR( smax );
			float*    buff_smin = FLA_FLOAT_PTR( smin );

			FLA_Bsvd_find_max_min_ops( m_A,
			                           buff_d, inc_d,
			                           buff_e, inc_e,
			                           buff_smax,
			                           buff_smin );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d    = FLA_DOUBLE_PTR( d );
			double*   buff_e    = FLA_DOUBLE_PTR( e );
			double*   buff_smax = FLA_DOUBLE_PTR( smax );
			double*   buff_smin = FLA_DOUBLE_PTR( smin );

			FLA_Bsvd_find_max_min_opd( m_A,
			                           buff_d, inc_d,
			                           buff_e, inc_e,
			                           buff_smax,
			                           buff_smin );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_max_min_ops( int       m_A,
                                     float*    buff_d, int inc_d, 
                                     float*    buff_e, int inc_e, 
                                     float*    smax,
                                     float*    smin )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_max_min_opd( int       m_A,
                                     double*   buff_d, int inc_d, 
                                     double*   buff_e, int inc_e, 
                                     double*   smax,
                                     double*   smin )
{
	double smax_cand;
	double smin_cand;
	int    i;

	smax_cand = fabs( buff_d[ (m_A-1)*inc_d ] );
	smin_cand = smax_cand;

	for ( i = 0; i < m_A - 1; ++i )
	{
		double abs_di = fabs( buff_d[ i*inc_d ] );
		double abs_ei = fabs( buff_e[ i*inc_e ] );

		// Track the minimum element.
		smin_cand = min( smin_cand, abs_di );

		// Track the maximum element.
		smax_cand = max( smax_cand, abs_di );
		smax_cand = max( smax_cand, abs_ei );
	}

	// Save the results of the search.
	*smax = smax_cand;
	*smin = smin_cand;

	return FLA_SUCCESS;
}

