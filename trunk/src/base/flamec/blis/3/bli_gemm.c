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

void bli_sgemm( trans_t transa, trans_t transb, int m, int k, int n, float* alpha, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs, float* beta, float* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	float*    a_save    = a;
	float*    b_save    = b;
	float*    c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	float     zero = bli_s0();
	float     one  = bli_s1();
	float*    a_unswap;
	float*    b_unswap;
	float*    c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;

	// Return early if possible.
	if ( bli_zero_dim3( m, k, n ) )
	{
		bli_sscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bli_screate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bli_screate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bli_screate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bli_is_col_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transa );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );

				bli_sswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bli_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bli_is_row_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bli_swap_ints( ldc, incc );
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transa );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_sswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bli_swap_ints( ldc, incc );
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transb );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_sswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_sswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bli_sallocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bli_sgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a,       lda,
		                b,       ldb,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bli_sscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bli_saxpymt( BLIS_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bli_sfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bli_sgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a, lda,
		                b, ldb,
		                beta,
		                c, ldc );
	}

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bli_sfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bli_sfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bli_sfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bli_dgemm( trans_t transa, trans_t transb, int m, int k, int n, double* alpha, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs, double* beta, double* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	double*   a_save    = a;
	double*   b_save    = b;
	double*   c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	double    zero = bli_d0();
	double    one  = bli_d1();
	double*   a_unswap;
	double*   b_unswap;
	double*   c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;

	// Return early if possible.
	if ( bli_zero_dim3( m, k, n ) )
	{
		bli_dscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bli_dcreate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bli_dcreate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bli_dcreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bli_is_col_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transa );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );

				bli_dswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bli_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bli_is_row_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bli_swap_ints( ldc, incc );
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transa );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_dswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bli_swap_ints( ldc, incc );
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transb );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_dswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_dswap_pointers( a, b );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bli_dallocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bli_dgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a,       lda,
		                b,       ldb,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bli_dscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bli_daxpymt( BLIS_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bli_dfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bli_dgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a, lda,
		                b, ldb,
		                beta,
		                c, ldc );
	}

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bli_dfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bli_dfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bli_dfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bli_cgemm( trans_t transa, trans_t transb, int m, int k, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	scomplex* a_save    = a;
	scomplex* b_save    = b;
	scomplex* c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	scomplex  zero = bli_c0();
	scomplex  one  = bli_c1();
	scomplex* a_unswap;
	scomplex* b_unswap;
	scomplex* a_conj;
	scomplex* b_conj;
	scomplex* c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       lda_conj, inca_conj;
	int       ldb_conj, incb_conj;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;
	int       a_was_copied;
	int       b_was_copied;

	// Return early if possible.
	if ( bli_zero_dim3( m, k, n ) )
	{
		bli_cscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bli_ccreate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bli_ccreate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bli_ccreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Figure out whether A and/or B was copied to contiguous memory. This
	// is used later to prevent redundant copying.
	a_was_copied = ( a != a_save );
	b_was_copied = ( b != b_save );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bli_is_col_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transa );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );

				bli_cswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bli_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bli_is_row_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bli_swap_ints( ldc, incc );
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transa );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_cswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bli_swap_ints( ldc, incc );
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transb );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_cswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_cswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
	}

	// We need a temporary matrix for the case where A is conjugated.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;

	// If transa indicates conjugate-no-transpose and A was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transa indicates conjugate-no-transpose and A was already copied,
	// just conjugate it.
	if ( bli_is_conjnotrans( transa ) && !a_was_copied )
	{
		a_conj    = bli_callocm( m_gemm, k );
		lda_conj  = m_gemm;
		inca_conj = 1;

		bli_ccopymt( BLIS_CONJ_NO_TRANSPOSE,
		             m_gemm,
		             k,
		             a,      inca,      lda,
		             a_conj, inca_conj, lda_conj );
	}
	else if ( bli_is_conjnotrans( transa ) && a_was_copied )
	{
		bli_cconjm( m_gemm,
		            k,
		            a_conj, inca_conj, lda_conj );
	}

	// We need a temporary matrix for the case where B is conjugated.
	b_conj    = b;
	ldb_conj  = ldb;
	incb_conj = incb;

	// If transb indicates conjugate-no-transpose and B was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transb indicates conjugate-no-transpose and B was already copied,
	// just conjugate it.
	if ( bli_is_conjnotrans( transb ) && !b_was_copied )
	{
		b_conj    = bli_callocm( k, n_gemm );
		ldb_conj  = k;
		incb_conj = 1;

		bli_ccopymt( BLIS_CONJ_NO_TRANSPOSE,
		             k,
		             n_gemm,
		             b,      incb,      ldb,
		             b_conj, incb_conj, ldb_conj );
	}
	else if ( bli_is_conjnotrans( transb ) && b_was_copied )
	{
		bli_cconjm( k,
		            n_gemm,
		            b_conj, incb_conj, ldb_conj );
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bli_callocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bli_cgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj,  lda_conj,
		                b_conj,  ldb_conj,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bli_cscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bli_caxpymt( BLIS_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bli_cfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bli_cgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj, lda_conj,
		                b_conj, ldb_conj,
		                beta,
		                c,      ldc );
	}

	if ( bli_is_conjnotrans( transa ) && !a_was_copied )
		bli_cfree( a_conj );

	if ( bli_is_conjnotrans( transb ) && !b_was_copied )
		bli_cfree( b_conj );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bli_cfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bli_cfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bli_cfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bli_zgemm( trans_t transa, trans_t transb, int m, int k, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	dcomplex* a_save    = a;
	dcomplex* b_save    = b;
	dcomplex* c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	dcomplex  zero = bli_z0();
	dcomplex  one  = bli_z1();
	dcomplex* a_unswap;
	dcomplex* b_unswap;
	dcomplex* a_conj;
	dcomplex* b_conj;
	dcomplex* c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       lda_conj, inca_conj;
	int       ldb_conj, incb_conj;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;
	int       a_was_copied;
	int       b_was_copied;

	// Return early if possible.
	if ( bli_zero_dim3( m, k, n ) )
	{
		bli_zscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bli_zcreate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bli_zcreate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bli_zcreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Figure out whether A and/or B was copied to contiguous memory. This
	// is used later to prevent redundant copying.
	a_was_copied = ( a != a_save );
	b_was_copied = ( b != b_save );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bli_is_col_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transa );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );

				bli_zswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bli_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bli_is_row_storage( c_rs, c_cs ) )
	{
		if ( bli_is_col_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bli_swap_ints( ldc, incc );
				bli_swap_ints( ldb, incb );

				bli_toggle_trans( transa );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_zswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
		else // if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			if ( bli_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bli_swap_ints( ldc, incc );
				bli_swap_ints( lda, inca );

				bli_toggle_trans( transb );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_zswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
			else // if ( bli_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
				bli_swap_ints( ldc, incc );

				bli_swap_ints( m, n );
				bli_swap_ints( m_gemm, n_gemm );
				bli_zswap_pointers( a, b );
				bli_swap_ints( a_was_copied, b_was_copied );
				bli_swap_ints( lda, ldb );
				bli_swap_ints( inca, incb );
				bli_swap_trans( transa, transb );
			}
		}
	}

	// We need a temporary matrix for the case where A is conjugated.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;

	// If transa indicates conjugate-no-transpose and A was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transa indicates conjugate-no-transpose and A was already copied,
	// just conjugate it.
	if ( bli_is_conjnotrans( transa ) && !a_was_copied )
	{
		a_conj    = bli_zallocm( m_gemm, k );
		lda_conj  = m_gemm;
		inca_conj = 1;

		bli_zcopymt( BLIS_CONJ_NO_TRANSPOSE,
		             m_gemm,
		             k,
		             a,      inca,      lda,
		             a_conj, inca_conj, lda_conj );
	}
	else if ( bli_is_conjnotrans( transa ) && a_was_copied )
	{
		bli_zconjm( m_gemm,
		            k,
		            a_conj, inca_conj, lda_conj );
	}

	// We need a temporary matrix for the case where B is conjugated.
	b_conj    = b;
	ldb_conj  = ldb;
	incb_conj = incb;

	// If transb indicates conjugate-no-transpose and B was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transb indicates conjugate-no-transpose and B was already copied,
	// just conjugate it.
	if ( bli_is_conjnotrans( transb ) && !b_was_copied )
	{
		b_conj    = bli_zallocm( k, n_gemm );
		ldb_conj  = k;
		incb_conj = 1;

		bli_zcopymt( BLIS_CONJ_NO_TRANSPOSE,
		             k,
		             n_gemm,
		             b,      incb,      ldb,
		             b_conj, incb_conj, ldb_conj );
	}
	else if ( bli_is_conjnotrans( transb ) && b_was_copied )
	{
		bli_zconjm( k,
		            n_gemm,
		            b_conj, incb_conj, ldb_conj );
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bli_zallocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bli_zgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj,  lda_conj,
		                b_conj,  ldb_conj,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bli_zscalm( BLIS_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bli_zaxpymt( BLIS_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bli_zfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bli_zgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj, lda_conj,
		                b_conj, ldb_conj,
		                beta,
		                c,      ldc );
	}

	if ( bli_is_conjnotrans( transa ) && !a_was_copied )
		bli_zfree( a_conj );

	if ( bli_is_conjnotrans( transb ) && !b_was_copied )
		bli_zfree( b_conj );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bli_zfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bli_zfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bli_zfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

// --- Classic routine wrappers ---

void bli_sgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, float* alpha, float* a, int lda, float* b, int ldb, float* beta, float* c, int ldc )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );
	bli_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_sgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             *alpha,
	             a, lda,
	             b, ldb,
	             *beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bli_param_map_to_netlib_trans( transa, &blas_transa );
	bli_param_map_to_netlib_trans( transb, &blas_transb );

	F77_sgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bli_dgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, double* alpha, double* a, int lda, double* b, int ldb, double* beta, double* c, int ldc )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );
	bli_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_dgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             *alpha,
	             a, lda,
	             b, ldb,
	             *beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bli_param_map_to_netlib_trans( transa, &blas_transa );
	bli_param_map_to_netlib_trans( transb, &blas_transb );

	F77_dgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bli_cgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );
	bli_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_cgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             alpha,
	             a, lda,
	             b, ldb,
	             beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bli_param_map_to_netlib_trans( transa, &blas_transa );
	bli_param_map_to_netlib_trans( transb, &blas_transb );

	F77_cgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bli_zgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );
	bli_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_zgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             alpha,
	             a, lda,
	             b, ldb,
	             beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bli_param_map_to_netlib_trans( transa, &blas_transa );
	bli_param_map_to_netlib_trans( transb, &blas_transb );

	F77_zgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

