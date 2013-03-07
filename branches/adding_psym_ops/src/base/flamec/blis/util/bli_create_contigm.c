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

void bli_screate_contigm( int m, int n, float* a_save, int a_rs_save, int a_cs_save, float** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bli_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bli_sallocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bli_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bli_scopymt( BLIS_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

void bli_dcreate_contigm( int m, int n, double* a_save, int a_rs_save, int a_cs_save, double** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bli_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bli_dallocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bli_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bli_dcopymt( BLIS_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

void bli_ccreate_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bli_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bli_callocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bli_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bli_ccopymt( BLIS_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

void bli_zcreate_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs )
{
	int m_contig, n_contig;

	if ( bli_is_gen_storage( a_rs_save, a_cs_save ) )
	{
		// Initialize dimensions assuming no transposition needed during copy.
		m_contig = m;
		n_contig = n;

/*
		// Transpose the dimensions of the contiguous matrix, if requested.
		if ( bli_does_trans( trans_copy ) )
		{
			m_contig = n;
			n_contig = m;
		}
*/

		// Allocate temporary contiguous storage for the matrix.
		*a = bli_zallocm( m_contig, n_contig );

		// Set the row and column strides for the temporary matrix.
		bli_set_contig_strides( m_contig, n_contig, a_rs, a_cs );

		// Initialize the contiguous matrix with the contents of the original.
		bli_zcopymt( BLIS_NO_TRANSPOSE,
		             m_contig,
		             n_contig,
		             a_save, a_rs_save, a_cs_save,
		             *a,     *a_rs,     *a_cs );
	}
}

