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


FLA_Error FLA_Bsvd_find_submatrix_ops( int       mn_A,
                                       int       ij_begin,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Bsvd_find_submatrix_opd( int       mn_A,
                                       int       ij_begin,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR )
{
	double rzero = bli_d0();
	int    ij_tl;
	int    ij_br;

	// Search for the first non-zero superdiagonal element starting at
	// the index specified by ij_begin.
	for ( ij_tl = ij_begin; ij_tl < mn_A - 1; ++ij_tl )
	{
		double* e1 = buff_e + (ij_tl  )*inc_e;

		// If we find a non-zero element, record it and break out of this
		// loop.
		if ( *e1 != rzero )
		{
#ifdef PRINTF
printf( "FLA_Bsvd_find_submatrix_opd: found non-zero superdiagonal element\n" );
printf( "                             e[%3d] = %22.19e\n", ij_tl, *e1 );
#endif
			*ijTL = ij_tl;
			break;
		}
	}

	// If ij_tl was incremented all the way up to mn_A - 1, then we didn't
	// find any non-zeros.
	if ( ij_tl == mn_A - 1 )
	{
#ifdef PRINTF
printf( "FLA_Bsvd_find_submatrix_opd: no submatrices found.\n" );
#endif
		return FLA_FAILURE;
	}

	// If we've gotten this far, then a non-zero superdiagonal element was
	// found. Now we must walk the remaining portion of the superdiagonal
	// to find the first zero element, or if one is not found, we simply
	// use the last element of the superdiagonal.
	for ( ij_br = ij_tl; ij_br < mn_A - 1; ++ij_br )
	{
		double* e1 = buff_e + (ij_br  )*inc_e;

		// If we find a zero element, record it and break out of this
		// loop.
		if ( *e1 == rzero )
		{
#ifdef PRINTF
printf( "FLA_Bsvd_find_submatrix_opd: found zero superdiagonal element\n" );
printf( "                             e[%3d] = %22.19e\n", ij_br, *e1 );
#endif
			break;
		}
	}

	// If a zero element was found, then ij_br should hold the index of
	// that element. If a zero element was not found, then ij_br should
	// hold mn_A - 1. Either way, we save the value and return success.
	*ijBR = ij_br;

	return FLA_SUCCESS;
}
