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

FLA_Error FLASH_Axpy_buffer_to_hier( FLA_Obj alpha, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj H )
{
	FLA_Obj      flat_matrix;
	FLA_Datatype datatype;
	FLA_Error    e_val;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_if_scalar( alpha );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_consistent_object_datatype( alpha, H );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_matrix_strides( m, n, rs, cs );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_submatrix_dims_and_offset( m, n, i, j, H );
		FLA_Check_error_code( e_val );
	}

	// Acquire the datatype from the hierarchical matrix object.
	datatype = FLASH_Obj_datatype( H );

	// Create a temporary conventional matrix object of the requested datatype
	// and dimensions and attach the given buffer containing the incoming data.
	FLA_Obj_create_without_buffer( datatype, m, n, &flat_matrix );
	FLA_Obj_attach_buffer( buffer, rs, cs, &flat_matrix );

	// Recurse through H, adding in the corresponding elements of flat_matrix,
	// starting at the (i,j) element offset.
	FLASH_Axpy_flat_to_hier( alpha, flat_matrix, i, j, H );

	// Free the object (but don't free the buffer!).
	FLA_Obj_free_without_buffer( &flat_matrix );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Axpy_hier_to_buffer( FLA_Obj alpha, dim_t i, dim_t j, FLA_Obj H, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs )
{
	FLA_Obj      flat_matrix;
	FLA_Datatype datatype;
	FLA_Error    e_val;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_if_scalar( alpha );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_consistent_object_datatype( alpha, H );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_matrix_strides( m, n, rs, cs );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_submatrix_dims_and_offset( m, n, i, j, H );
		FLA_Check_error_code( e_val );
	}

	// Acquire the datatype from the hierarchical matrix object.
	datatype = FLASH_Obj_datatype( H );

	// Create a temporary conventional matrix object of the requested datatype
	// and dimensions and attach the given buffer containing the incoming data.
	FLA_Obj_create_without_buffer( datatype, m, n, &flat_matrix );
	FLA_Obj_attach_buffer( buffer, rs, cs, &flat_matrix );

	// Recurse through H, adding in the corresponding elements of flat_matrix,
	// starting at the (i,j) element offset.
	FLASH_Axpy_hier_to_flat( alpha, i, j, H, flat_matrix );

	// Free the object (but don't free the buffer!).
	FLA_Obj_free_without_buffer( &flat_matrix );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Axpy_flat_to_hier( FLA_Obj alpha, FLA_Obj F, dim_t i, dim_t j, FLA_Obj H )
{
	FLA_Obj HTL, HTR,
	        HBL, HBR;
	FLA_Obj HBR_tl, HBR_tr,
	        HBR_bl, HBR_br;
	dim_t   m, n;

	m = FLA_Obj_length( F );
	n = FLA_Obj_width( F );

	FLASH_Part_create_2x2( H,   &HTL, &HTR,
	                            &HBL, &HBR,    i, j, FLA_TL );

	FLASH_Part_create_2x2( HBR,   &HBR_tl, &HBR_tr,
	                              &HBR_bl, &HBR_br,    m, n, FLA_TL );

	FLASH_Axpy_hierarchy( FLA_FLAT_TO_HIER, alpha, F, &HBR_tl );

	FLASH_Part_free_2x2( &HBR_tl, &HBR_tr,
	                     &HBR_bl, &HBR_br );

	FLASH_Part_free_2x2( &HTL, &HTR,
	                     &HBL, &HBR );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Axpy_hier_to_flat( FLA_Obj alpha, dim_t i, dim_t j, FLA_Obj H, FLA_Obj F )
{
	FLA_Obj HTL, HTR,
	        HBL, HBR;
	FLA_Obj HBR_tl, HBR_tr,
	        HBR_bl, HBR_br;
	dim_t   m, n;

	m = FLA_Obj_length( F );
	n = FLA_Obj_width( F );

	FLASH_Part_create_2x2( H,   &HTL, &HTR,
	                            &HBL, &HBR,    i, j, FLA_TL );

	FLASH_Part_create_2x2( HBR,   &HBR_tl, &HBR_tr,
	                              &HBR_bl, &HBR_br,    m, n, FLA_TL );

	FLASH_Axpy_hierarchy( FLA_HIER_TO_FLAT, alpha, F, &HBR_tl );

	FLASH_Part_free_2x2( &HBR_tl, &HBR_tr,
	                     &HBR_bl, &HBR_br );

	FLASH_Part_free_2x2( &HTL, &HTR,
	                     &HBL, &HBR );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Axpy_hierarchy( int direction, FLA_Obj alpha, FLA_Obj F, FLA_Obj* H )
{
	// Once we get down to a submatrix whose elements are scalars, we are down
	// to our base case.
	if ( FLA_Obj_elemtype( *H ) == FLA_SCALAR )
	{
		// Depending on which top-level function invoked us, we either axpy
		// the source data in the flat matrix to the leaf-level submatrix of
		// the hierarchical matrix, or axpy the data in the hierarchical
		// submatrix to the flat matrix.
		if      ( direction == FLA_FLAT_TO_HIER )
		{
#ifdef FLA_ENABLE_SCC
			if ( FLA_is_owner() )
#endif
			FLA_Axpy_external( alpha, F, *H );
		}
		else if ( direction == FLA_HIER_TO_FLAT )
		{
#ifdef FLA_ENABLE_SCC
			if ( FLA_is_owner() )
#endif
			FLA_Axpy_external( alpha, *H, F );
		}
	}
	else
	{
		FLA_Obj HL,  HR,       H0,  H1,  H2;
		FLA_Obj FL,  FR,       F0,  F1,  F2;

		FLA_Obj H1T,           H01,
		        H1B,           H11,
		                       H21;
		FLA_Obj F1T,           F01,
		        F1B,           F11,
		                       F21;

		dim_t b_m;
		dim_t b_n;

		FLA_Part_1x2( *H,    &HL,  &HR,      0, FLA_LEFT );
		FLA_Part_1x2(  F,    &FL,  &FR,      0, FLA_LEFT );

		while ( FLA_Obj_width( HL ) < FLA_Obj_width( *H ) )
		{
			FLA_Repart_1x2_to_1x3( HL,  /**/ HR,        &H0, /**/ &H1, &H2,
			                       1, FLA_RIGHT );

			// Get the scalar width of H1 and use that to determine the
			// width of F1.
			b_n = FLASH_Obj_scalar_width( H1 );

			FLA_Repart_1x2_to_1x3( FL,  /**/ FR,        &F0, /**/ &F1, &F2,
			                       b_n, FLA_RIGHT );

			// -------------------------------------------------------------

			FLA_Part_2x1( H1,    &H1T,
			                     &H1B,       0, FLA_TOP );
			FLA_Part_2x1( F1,    &F1T,
			                     &F1B,       0, FLA_TOP );

			while ( FLA_Obj_length( H1T ) < FLA_Obj_length( H1 ) )
			{
				FLA_Repart_2x1_to_3x1( H1T,               &H01,
				                    /* ** */            /* *** */
				                                          &H11,
				                       H1B,               &H21,        1, FLA_BOTTOM );

				// Get the scalar length of H11 and use that to determine the
				// length of F11.
				b_m = FLASH_Obj_scalar_length( H11 );

				FLA_Repart_2x1_to_3x1( F1T,               &F01,
				                    /* ** */            /* *** */
				                                          &F11,
				                       F1B,               &F21,        b_m, FLA_BOTTOM );
				// -------------------------------------------------------------

				// Recursively axpy between F11 and H11.
				FLASH_Axpy_hierarchy( direction, alpha, F11,
				                      FLASH_OBJ_PTR_AT( H11 ) );

				// -------------------------------------------------------------

				FLA_Cont_with_3x1_to_2x1( &H1T,               H01,
				                                              H11,
				                        /* ** */           /* *** */
				                          &H1B,               H21,     FLA_TOP );
				FLA_Cont_with_3x1_to_2x1( &F1T,               F01,
				                                              F11,
				                        /* ** */           /* *** */
				                          &F1B,               F21,     FLA_TOP );
			}

			// -------------------------------------------------------------

			FLA_Cont_with_1x3_to_1x2( &HL,  /**/ &HR,        H0, H1, /**/ H2,
			                          FLA_LEFT );
			FLA_Cont_with_1x3_to_1x2( &FL,  /**/ &FR,        F0, F1, /**/ F2,
			                          FLA_LEFT );
		}
	}

	return FLA_SUCCESS;
}

