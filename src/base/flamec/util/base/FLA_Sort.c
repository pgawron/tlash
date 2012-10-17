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

int fla_scomp_f( const void* a, const void* b );
int fla_scomp_b( const void* a, const void* b );
int fla_dcomp_f( const void* a, const void* b );
int fla_dcomp_b( const void* a, const void* b );

FLA_Error FLA_Sort( FLA_Direct direct, FLA_Obj x )
{
	FLA_Datatype datatype;
	FLA_Obj      x_use;
	dim_t        m_x;
	dim_t        inc_x;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Sort_check( direct, x );

	datatype = FLA_Obj_datatype( x );

	m_x      = FLA_Obj_vector_dim( x );
	inc_x    = FLA_Obj_vector_inc( x );

	// If the vector does not have unit stride, copy it to a temporary vector
	// that does have unit stride.
	if ( inc_x != 1 )
	{
		FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, x, &x_use );
		inc_x = FLA_Obj_vector_inc( x_use );
	}
	else
	{
		x_use = x;
	}

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float* x_p = ( float* ) FLA_FLOAT_PTR( x_use );

			if ( direct == FLA_FORWARD )
				FLA_Sort_f_ops( m_x,
				                x_p, inc_x );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_b_ops( m_x,
				                x_p, inc_x );

			break;
		}

		case FLA_DOUBLE:
		{
			double* x_p = ( double* ) FLA_DOUBLE_PTR( x_use );

			if ( direct == FLA_FORWARD )
				FLA_Sort_f_opd( m_x,
				                x_p, inc_x );
			else // if ( direct == FLA_BACKWARD )
				FLA_Sort_b_opd( m_x,
				                x_p, inc_x );

			break;
		}

	}

	if ( inc_x != 1 )
	{
		FLA_Copy( x_use, x );
		FLA_Obj_free( &x_use );
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Sort_f_ops( int     m_x,
                          float*  x, int inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( float ),
	       fla_scomp_f );

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_b_ops( int     m_x,
                          float*  x, int inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( float ),
	       fla_scomp_b );

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_f_opd( int     m_x,
                          double* x, int inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( double ),
	       fla_dcomp_f );

	return FLA_SUCCESS;
}

FLA_Error FLA_Sort_b_opd( int     m_x,
                          double* x, int inc_x )
{
	qsort( x,
	       m_x,
	       sizeof( double ),
	       fla_dcomp_b );

	return FLA_SUCCESS;
}





int fla_scomp_f( const void* a, const void* b )
{
	float*  da = ( float*  ) a;
	float*  db = ( float*  ) b;
	int     r_val;

	if      ( *da < *db ) r_val = -1;
	else if ( *da > *db ) r_val =  1;
	else                  r_val =  0;

	return r_val;
}

int fla_scomp_b( const void* a, const void* b )
{
	float*  da = ( float*  ) a;
	float*  db = ( float*  ) b;
	int     r_val;

	if      ( *da < *db ) r_val =  1;
	else if ( *da > *db ) r_val = -1;
	else                  r_val =  0;

	return r_val;
}

int fla_dcomp_f( const void* a, const void* b )
{
	double* da = ( double* ) a;
	double* db = ( double* ) b;
	int     r_val;

	if      ( *da < *db ) r_val = -1;
	else if ( *da > *db ) r_val =  1;
	else                  r_val =  0;

	return r_val;
}

int fla_dcomp_b( const void* a, const void* b )
{
	double* da = ( double* ) a;
	double* db = ( double* ) b;
	int     r_val;

	if      ( *da < *db ) r_val =  1;
	else if ( *da > *db ) r_val = -1;
	else                  r_val =  0;

	return r_val;
}

