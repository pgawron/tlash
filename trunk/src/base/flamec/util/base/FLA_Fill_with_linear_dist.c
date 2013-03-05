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

FLA_Error FLA_Fill_with_linear_dist( FLA_Obj shift, FLA_Obj delta, FLA_Obj x )
{
	FLA_Obj      lT,              l0,
	             lB,              lambda1,
	                              l2;
	FLA_Obj      l, k, delta2;
	FLA_Datatype dt_real;
	dim_t        n_x;


	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Fill_with_linear_dist_check( shift, delta, x );

	dt_real = FLA_Obj_datatype_proj_to_real( x );
	n_x     = FLA_Obj_vector_dim( x );

	// Create a local counter to increment as we create the distribution.
	FLA_Obj_create( dt_real, 1,   1, 0, 0, &k );

	// Create a local vector l. We will work with this vector, which is
	// the same length as x, so that we can use vertical partitioning.
	FLA_Obj_create( dt_real, n_x, 1, 0, 0, &l );

	// Create a local real scalar alpha2 of the same precision as
	// alpha. Then copy alpha to alpha2, which will convert the
	// complex value to real, if necessary (ie: if alpha is complex).
	FLA_Obj_create( dt_real, 1,   1, 0, 0, &delta2 );
	FLA_Copy( delta, delta2 );

	// Initialize k to shift + delta2.
	FLA_Set( shift, k );
	FLA_Mult_add( FLA_ONE, delta2, k );

	FLA_Part_2x1( l,    &lT,
	                    &lB,            0, FLA_TOP );

	while ( FLA_Obj_length( lB ) > 0 )
	{
		FLA_Repart_2x1_to_3x1( lT,                &l0,
		                    /* ** */            /* ******* */
		                                          &lambda1,
		                       lB,                &l2,        1, FLA_BOTTOM );

		/*------------------------------------------------------------*/

		// lambda1 = k;
		FLA_Copy( k, lambda1 );

		// k = k + delta2;
		FLA_Mult_add( FLA_ONE, delta2, k );

		/*------------------------------------------------------------*/

		FLA_Cont_with_3x1_to_2x1( &lT,                l0,
		                                              lambda1,
		                        /* ** */           /* ******* */
		                          &lB,                l2,     FLA_TOP );
	}

	// Normalize by last element.
	//FLA_Part_2x1( l,    &lT,
	//                    &lB,            1, FLA_BOTTOM );
	//FLA_Inv_scal( lB, l );

	// Overwrite x with the distribution we created in l.
	// If x is complex, then this is where the conversion between
	// datatypes happens.
	FLA_Copy( l, x );

	FLA_Obj_free( &l );
	FLA_Obj_free( &k );
	FLA_Obj_free( &delta2 );

	return FLA_SUCCESS;
}

