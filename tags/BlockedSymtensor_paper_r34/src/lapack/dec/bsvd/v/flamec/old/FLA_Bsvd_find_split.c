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

FLA_Error FLA_Bsvd_find_split( FLA_Obj d, FLA_Obj e )
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
			float*    buff_d      = FLA_FLOAT_PTR( d );
			float*    buff_e      = FLA_FLOAT_PTR( e );

			FLA_Bsvd_find_split_ops( m_A,
			                         buff_d, inc_d,
			                         buff_e, inc_e );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d      = FLA_DOUBLE_PTR( d );
			double*   buff_e      = FLA_DOUBLE_PTR( e );

			FLA_Bsvd_find_split_opd( m_A,
			                         buff_d, inc_d,
			                         buff_e, inc_e );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_split_ops( int       m_A,
                                   float*    buff_d, int inc_d, 
                                   float*    buff_e, int inc_e )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_split_opd( int       m_A,
                                   double*   buff_d, int inc_d, 
                                   double*   buff_e, int inc_e )
{
	int i;

	for ( i = 0; i < m_A - 1; ++i )
	{
		double* epsilon1 = buff_e + (i  )*inc_e;

		if ( *epsilon1 == 0.0 )
		{
			// Return index of split as i+1 since e_i is in the same
			// column as d_(i+1).
			return i + 1;
		}
	}

	// Return with no split found found.
	return FLA_FAILURE;
}

