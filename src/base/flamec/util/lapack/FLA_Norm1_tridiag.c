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

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm )
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
			float*    buff_norm = FLA_FLOAT_PTR( norm );

			FLA_Norm1_tridiag_ops( m_A,
			                       buff_d, inc_d,
			                       buff_e, inc_e,
			                       buff_norm );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d    = FLA_DOUBLE_PTR( d );
			double*   buff_e    = FLA_DOUBLE_PTR( e );
			double*   buff_norm = FLA_DOUBLE_PTR( norm );

			FLA_Norm1_tridiag_opd( m_A,
			                       buff_d, inc_d,
			                       buff_e, inc_e,
			                       buff_norm );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Norm1_tridiag_ops( int       m_A,
                                 float*    buff_d, int inc_d, 
                                 float*    buff_e, int inc_e,
                                 float*    norm )
{
	float*  d  = buff_d;
	float*  e  = buff_e;
	float   nm;
	int     i;

	if ( m_A == 1 )
	{
		nm = fabs( *d );
	}
	else
	{
		float  d_first = d[ (0    )*inc_d ];
		float  e_first = e[ (0    )*inc_e ];
		float  e_last  = e[ (m_A-2)*inc_e ];
		float  d_last  = d[ (m_A-1)*inc_d ];

		// Record the maximum of the absolute row/column sums for the
		// first and last row/columns.
		nm = max( fabs( d_first ) + fabs( e_first ),
		          fabs( e_last  ) + fabs( d_last  ) );

		for ( i = 1; i < m_A - 2; ++i )
		{
			float  e0 = e[ (i-1)*inc_e ];
			float  e1 = e[ (i  )*inc_e ];
			float  d1 = d[ (i  )*inc_d ];

			// Update nm with the absolute row/column sum for the ith
			// row/column.
			nm = max( nm, fabs( e0 ) +
			              fabs( d1 ) +
			              fabs( e1 ) );
		}
	}

	*norm = nm;

	return FLA_SUCCESS;
}



FLA_Error FLA_Norm1_tridiag_opd( int       m_A,
                                 double*   buff_d, int inc_d, 
                                 double*   buff_e, int inc_e,
                                 double*   norm )
{
	double* d  = buff_d;
	double* e  = buff_e;
	double  nm;
	int     i;

	if ( m_A == 1 )
	{
		nm = fabs( *d );
	}
	else
	{
		double d_first = d[ (0    )*inc_d ];
		double e_first = e[ (0    )*inc_e ];
		double e_last  = e[ (m_A-2)*inc_e ];
		double d_last  = d[ (m_A-1)*inc_d ];

		// Record the maximum of the absolute row/column sums for the
		// first and last row/columns.
		nm = max( fabs( d_first ) + fabs( e_first ),
		          fabs( e_last  ) + fabs( d_last  ) );

		for ( i = 1; i < m_A - 2; ++i )
		{
			double e0 = e[ (i-1)*inc_e ];
			double e1 = e[ (i  )*inc_e ];
			double d1 = d[ (i  )*inc_d ];

			// Update nm with the absolute row/column sum for the ith
			// row/column.
			nm = max( nm, fabs( e0 ) +
			              fabs( d1 ) +
			              fabs( e1 ) );
		}
	}

	*norm = nm;

	return FLA_SUCCESS;
}

