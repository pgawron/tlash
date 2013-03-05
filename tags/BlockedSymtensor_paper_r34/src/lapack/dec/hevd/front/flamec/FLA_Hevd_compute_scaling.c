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

FLA_Error FLA_Hevd_compute_scaling( FLA_Uplo uplo, FLA_Obj A, FLA_Obj sigma )
{
	FLA_Datatype dt_real;
	FLA_Obj      norm;
	FLA_Obj      safmin;
	FLA_Obj      prec;
	FLA_Obj      rmin;
	FLA_Obj      rmax;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Hevd_compute_scaling_check( uplo, A, sigma );

	dt_real = FLA_Obj_datatype_proj_to_real( A );

	FLA_Obj_create( dt_real, 1, 1, 0, 0, &norm );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &safmin );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &prec );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &rmin );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &rmax );

	// Query safmin, precision.
	FLA_Mach_params( FLA_MACH_SFMIN, safmin );
	FLA_Mach_params( FLA_MACH_PREC,  prec );

//FLA_Obj_show( "safmin", safmin, "%20.12e", "" );
//FLA_Obj_show( "prec", prec, "%20.12e", "" );

	// rmin = sqrt( safmin / prec );
	FLA_Copy( safmin, rmin );
	FLA_Inv_scal( prec, rmin );
	FLA_Copy( rmin, rmax );
	FLA_Sqrt( rmin );

	// rmax = sqrt( 1 / ( safmin / prec ) );
	FLA_Invert( FLA_NO_CONJUGATE, rmax );
	FLA_Sqrt( rmax );

//FLA_Obj_show( "rmin", rmin, "%20.12e", "" );
//FLA_Obj_show( "rmax", rmax, "%20.12e", "" );

	// Find the maximum absolute value of A.
	FLA_Max_abs_value_herm( uplo, A, norm );

	if ( FLA_Obj_gt( norm, FLA_ZERO ) && FLA_Obj_lt( norm, rmin ) )
	{
		// sigma = rmin / norm;
		FLA_Copy( rmin, sigma );
		FLA_Inv_scal( norm, sigma );
	}
	else if ( FLA_Obj_gt( norm, rmax ) )
	{
		// sigma = rmax / norm;
		FLA_Copy( rmax, sigma );
		FLA_Inv_scal( norm, sigma );
	}
	else
	{
		// sigma = 1.0;
		FLA_Copy( FLA_ONE, sigma );
	}

	FLA_Obj_free( &norm );
	FLA_Obj_free( &safmin );
	FLA_Obj_free( &prec );
	FLA_Obj_free( &rmin );
	FLA_Obj_free( &rmax );

	return FLA_SUCCESS;
}

