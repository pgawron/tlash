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

#include "blis.h"

void bli_sfree_saved_contigmr( uplo_t uplo, int m, int n, float* a_save, int a_rs_save, int a_cs_save, float** a, int* a_rs, int* a_cs )
{
	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bli_scopymr( uplo,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bli_sfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bli_dfree_saved_contigmr( uplo_t uplo, int m, int n, double* a_save, int a_rs_save, int a_cs_save, double** a, int* a_rs, int* a_cs )
{
	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bli_dcopymr( uplo,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bli_dfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bli_cfree_saved_contigmr( uplo_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs )
{
	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bli_ccopymr( uplo,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bli_cfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

void bli_zfree_saved_contigmr( uplo_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs )
{
	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Copy the contents of the temporary matrix back to the original.
		bli_zcopymr( uplo,
		             m,
		             n,
		             *a,     *a_rs,     *a_cs,
		             a_save, a_rs_save, a_cs_save );

		// Free the temporary contiguous storage for the matrix.
		bli_zfree( *a );

		// Restore the original matrix address.
		*a = a_save;

		// Restore the original row and column strides.
		*a_rs = a_rs_save;
		*a_cs = a_cs_save;
	}
}

