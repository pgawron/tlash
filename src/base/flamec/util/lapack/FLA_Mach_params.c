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

FLA_Error FLA_Mach_params( FLA_Machval machval, FLA_Obj val )
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( val );

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Mach_params_check( machval, val );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*  val_p = ( float* ) FLA_FLOAT_PTR( val );

			*val_p = FLA_Mach_params_ops( machval );

			break;
		}

		case FLA_DOUBLE:
		{
			double* val_p = ( double* ) FLA_DOUBLE_PTR( val );

			*val_p = FLA_Mach_params_opd( machval );

			break;
		}
	}

	return FLA_SUCCESS;
}


float FLA_Mach_params_ops( FLA_Machval machval )
{
	static int    first_time = TRUE;
	static float  vals[FLA_MACH_N_VALS];

	if ( first_time )
	{
		char lapack_machval;
		int  i;

		for( i = 0; i < FLA_MACH_N_VALS - 1; ++i )
		{
			FLA_Param_map_flame_to_netlib_machval( FLA_MACH_START + i, &lapack_machval );
//printf( "querying %d %c\n", FLA_MACH_START + i, lapack_machval );
			vals[i] = fla_slamch( &lapack_machval, 1 );
//printf( "got back  %34.29e\n", vals[i] );
		}

		// Store epsilon^2 in the last element.
		vals[i] = vals[0] * vals[0];

		first_time = FALSE;
	}

	return vals[ machval - FLA_MACH_START ];
}

double FLA_Mach_params_opd( FLA_Machval machval )
{
	static int    first_time = TRUE;
	static double vals[FLA_MACH_N_VALS];

	if ( first_time )
	{
		char lapack_machval;
		int  i;

		for( i = 0; i < FLA_MACH_N_VALS - 1; ++i )
		{
			FLA_Param_map_flame_to_netlib_machval( FLA_MACH_START + i, &lapack_machval );
//printf( "querying %d %c\n", FLA_MACH_START + i, lapack_machval );
			vals[i] = fla_dlamch( &lapack_machval, 1 );
//printf( "got back  %34.29e\n", vals[i] );
		}

		// Store epsilon^2 in the last element.
		vals[i] = vals[0] * vals[0];

		first_time = FALSE;
	}

	return vals[ machval - FLA_MACH_START ];
}

