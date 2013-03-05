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

FLA_Error FLA_Apply_G_rf_asm_var5b( FLA_Obj G, FLA_Obj A )
/*
  Apply k sets of Givens rotations to a matrix A from the right,
  where each set takes the form:

    A := A ( G(n-1,k) ... G(1,k) G(0,k) )'
       = A G(0,k)' G(1,k)' ... G(n-1,k)'

  where Gik is the ith Givens rotation formed from the kth set,
  stored in the (i,k) entries of of C and S:

    Gik  =  / gamma_ik  -sigma_ik \
            \ sigma_ik   gamma_ik /

  -FGVZ
*/
{
	FLA_Datatype datatype;
	int          k_G, m_A, n_A;
	int          rs_G, cs_G;
	int          rs_A, cs_A;

	datatype = FLA_Obj_datatype( A );

	k_G      = FLA_Obj_width( G );
	m_A      = FLA_Obj_length( A );
	n_A      = FLA_Obj_width( A );

	rs_G     = FLA_Obj_row_stride( G );
	cs_G     = FLA_Obj_col_stride( G );

	rs_A     = FLA_Obj_row_stride( A );
	cs_A     = FLA_Obj_col_stride( A );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			float*    buff_A = ( float*    ) FLA_FLOAT_PTR( A );

			FLA_Apply_G_rf_ass_var5b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_A = ( double*   ) FLA_DOUBLE_PTR( A );

			FLA_Apply_G_rf_asd_var5b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

			FLA_Apply_G_rf_asc_var5b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

			FLA_Apply_G_rf_asz_var5b( k_G,
			                         m_A,
			                         n_A,
			                         0,
			                         0,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_ass_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asd_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A )
{
	double             one  = bli_d1();
	double             zero = bli_d0();
	double             gamma12, sigma12;
	double             gamma23, sigma23;
	double*   restrict a1;
	double*   restrict a2;
	double*   restrict a3;
	dcomplex* restrict g12;
	dcomplex* restrict g23;
	int                j, g, k;
	int                nG, nG_app;
	int                k_minus_1;
	int                is_ident12;
	int                is_ident23;

	int                n_run  = ( n_A - 1 ) / 2;
	int                n_left = ( n_A - 1 ) % 2;
	int                m_app;
	int                m_base = i_k + 4 - iTL;


	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 3(k - 1).

	if ( nG < 3*k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_asd_var1( k_G,
		                         m_A,
		                         n_A,
		                         buff_G, rs_G, cs_G,
		                         buff_A, rs_A, cs_A );
		return FLA_SUCCESS;
	}


	// Start-up phase.

	for ( j = 0; j < k_minus_1; ++j )
	{
		nG_app = j + 1;

		for ( k = 0, g = nG_app - 1; k < nG_app; ++k, --g )
		{
			g12   = buff_G + (2*g    )*rs_G + (k  )*cs_G;
			g23   = buff_G + (2*g + 1)*rs_G + (k  )*cs_G;
			a1    = buff_A + (2*g    )*cs_A;
			a2    = buff_A + (2*g + 1)*cs_A;
			a3    = buff_A + (2*g + 2)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;
			gamma23 = g23->real;
			sigma23 = g23->imag;

			is_ident12 = ( gamma12 == one && sigma12 == zero );
			is_ident23 = ( gamma23 == one && sigma23 == zero );

			//m_app = min( i_k + 4 + 2*j - iTL, m_A );
			m_app = min( m_base + 2*j, m_A );

			if      ( !is_ident12 && is_ident23 )
			{
				// Apply only to columns 1 and 2.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 2,
				                     a2, 2 );
			}
			else if ( is_ident12 && !is_ident23 )
			{
				// Apply only to columns 2 and 3.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma23,
				                     &sigma23,
				                     a2, 2,
				                     a3, 2 );
			}
			else if ( !is_ident12 && !is_ident23 )
			{
				// Apply to all three columns.

				MAC_Apply_G_mx3_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     &gamma23,
				                     &sigma23,
				                     a1, 2,
				                     a2, 2,
				                     a3, 2 );
			}
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < n_run; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g12   = buff_G + (2*g    )*rs_G + (k  )*cs_G;
			g23   = buff_G + (2*g + 1)*rs_G + (k  )*cs_G;
			a1    = buff_A + (2*g    )*cs_A;
			a2    = buff_A + (2*g + 1)*cs_A;
			a3    = buff_A + (2*g + 2)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;
			gamma23 = g23->real;
			sigma23 = g23->imag;

			is_ident12 = ( gamma12 == one && sigma12 == zero );
			is_ident23 = ( gamma23 == one && sigma23 == zero );

//printf( "%03d %03d %03d: m_app m_left %d %d %d\n", j, k, g, m_app, m_app % 2, m_app / 2 );
			//m_app = min( i_k + 4 + 2*j - iTL, m_A );
			m_app = min( m_base + 2*j, m_A );

			if      ( !is_ident12 && is_ident23 )
			{
				// Apply only to columns 1 and 2.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 2,
				                     a2, 2 );
			}
			else if ( is_ident12 && !is_ident23 )
			{
				// Apply only to columns 2 and 3.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma23,
				                     &sigma23,
				                     a2, 2,
				                     a3, 2 );
			}
			else if ( !is_ident12 && !is_ident23 )
			{
				// Apply to all three columns.

				MAC_Apply_G_mx3_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     &gamma23,
				                     &sigma23,
				                     a1, 2,
				                     a2, 2,
				                     a3, 2 );
			}
		}
	}


	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG + 1; ++j )
	{
		nG_app = nG - j;

		if ( n_left == 1 )
		{
			k     = k_G - 1 - nG_app;
			g12   = buff_G + (n_A-2)*rs_G + (k  )*cs_G;
			a1    = buff_A + (n_A-2)*cs_A;
			a2    = buff_A + (n_A-1)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;

			is_ident12 = ( gamma12 == one && sigma12 == zero );

			m_app = m_A;

			if ( !is_ident12 )
			{
				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 2,
				                     a2, 2 );
			}
		}

		for ( k = k_G - nG_app, g = n_run - 1; k < k_G; ++k, --g )
		{
			g12   = buff_G + (2*g    )*rs_G + (k  )*cs_G;
			g23   = buff_G + (2*g + 1)*rs_G + (k  )*cs_G;
			a1    = buff_A + (2*g    )*cs_A;
			a2    = buff_A + (2*g + 1)*cs_A;
			a3    = buff_A + (2*g + 2)*cs_A;

			gamma12 = g12->real;
			sigma12 = g12->imag;
			gamma23 = g23->real;
			sigma23 = g23->imag;

			//m_app = min( i_k + 4 + 2*j - iTL, m_A );
			//m_app = min( m_base + 2*j, m_A );
			m_app = m_A;

			is_ident12 = ( gamma12 == one && sigma12 == zero );
			is_ident23 = ( gamma23 == one && sigma23 == zero );

			if      ( !is_ident12 && is_ident23 )
			{
				// Apply only to columns 1 and 2.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     a1, 2,
				                     a2, 2 );
			}
			else if ( is_ident12 && !is_ident23 )
			{
				// Apply only to columns 2 and 3.

				MAC_Apply_G_mx2_asd( m_app,
				                     &gamma23,
				                     &sigma23,
				                     a2, 2,
				                     a3, 2 );
			}
			else if ( !is_ident12 && !is_ident23 )
			{
				// Apply to all three columns.

				MAC_Apply_G_mx3_asd( m_app,
				                     &gamma12,
				                     &sigma12,
				                     &gamma23,
				                     &sigma23,
				                     a1, 2,
				                     a2, 2,
				                     a3, 2 );
			}
		}
	}
	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asc_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asz_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A )
{
	return FLA_SUCCESS;
}

