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


FLA_Error FLA_Tevd_find_perfshift_ops( int       m_d,
                                       int       m_l,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e, 
                                       float*    buff_l, int inc_l, 
                                       int*      buff_ls, int inc_ls, 
                                       float*    buff_pu, int inc_pu, 
                                       int*      ij_shift )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_find_perfshift_opd( int       m_d,
                                       int       m_l,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e, 
                                       double*   buff_l, int inc_l, 
                                       int*      buff_ls, int inc_ls, 
                                       double*   buff_pu, int inc_pu, 
                                       int*      ij_shift )
{
	double* d1p;
	double* e1p;
	double* d2p;
	double  wilkshift;
	int     i;
	int     ij_cand;
	double  dist_cand;
	double  pshift_cand;
	
	d1p = buff_d + (m_d-2)*inc_d;
	e1p = buff_e + (m_d-2)*inc_e;
	d2p = buff_d + (m_d-1)*inc_d;

	if ( *buff_ls == -1 )
	{
		*ij_shift = -1;
		return FLA_FAILURE;
	}

	FLA_Wilkshift_tridiag_opd( *d1p,
	                           *e1p,
	                           *d2p,
	                           &wilkshift );

/*
	// If we have shifted here previously, use a Wilkinson shfit.
	prev_shift = buff_pu[ (m_d-1)*inc_pu ];

	if ( prev_shift != 0.0 )
	{
		// *shift = prev_shift;
		*shift = wilkshift;
		return FLA_SUCCESS;
	}
*/
	ij_cand = -1;

	// Find an available (unused) shift.
	for ( i = 0; i < m_l; ++i )
	{
		int* status = buff_ls + (i  )*inc_ls;

		if ( *status == 0 )
		{
			double* lambda1 = buff_l + (i  )*inc_l;
			ij_cand     = i;
			pshift_cand = *lambda1;
			dist_cand   = fabs( wilkshift - pshift_cand );
		}
	}

	if ( ij_cand == -1 )
	{
		*ij_shift = -1;
		*buff_ls  = -1;
		return FLA_FAILURE;
	}

	// Now try to find a shift closer to wilkshift than the
	// first one we found.
	for ( i = 0; i < m_l; ++i )
	{
		double* lambda1 = buff_l  + (i  )*inc_l;
		int*    status  = buff_ls + (i  )*inc_ls;
		double  dist    = fabs( wilkshift - *lambda1 );

		if ( *status == 1 ) continue;

		if ( dist < dist_cand )
		{
			ij_cand = i;
			pshift_cand = *lambda1;
			dist_cand = dist;
		}
	}

	*ij_shift = ij_cand;

	return FLA_SUCCESS;
}

