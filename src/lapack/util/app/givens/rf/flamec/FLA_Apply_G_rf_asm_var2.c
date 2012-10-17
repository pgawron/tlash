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

FLA_Error FLA_Apply_G_rf_asm_var2( FLA_Obj G, FLA_Obj A )
/*
  Apply k sets of Givens rotations to a matrix A from the right,
  where each set takes the form:

    A := A ( G(n-1,k) ... G(1,k) G(0,k) )'
       = A G(0,k)' G(1,k)' ... G(n-1,k)'

  where Gik is the ith Givens rotation formed from the kth set,
  stored in the (i,k) entries of of G:

    Gik  =  / gamma_ik  -sigma_ik \
            \ sigma_ik   gamma_ik /

  This variant iterates in pipelined, overlapping fashion and
  applies rotations to two columns at a time.

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

			FLA_Apply_G_rf_ass_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_A = ( double*   ) FLA_DOUBLE_PTR( A );

			FLA_Apply_G_rf_asd_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
			scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

			FLA_Apply_G_rf_asc_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

			FLA_Apply_G_rf_asz_var2( k_G,
			                         m_A,
			                         n_A,
			                         buff_G, rs_G, cs_G,
			                         buff_A, rs_A, cs_A );

			break;
		}
	}

	return FLA_SUCCESS;
}


FLA_Error FLA_Apply_G_rf_ass_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A )
{
	float     one  = bli_s1();
	float     zero = bli_s0();
	float     gamma;
	float     sigma;
	float*    a1;
	float*    a2;
	scomplex* g11;
	int       j, g, k;
	int       nG, nG_app;
	int       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_ass_var1( k_G,
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
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_ass( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_ass( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_ass( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asd_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A )
{
	double    one  = bli_d1();
	double    zero = bli_d0();
	double    gamma;
	double    sigma;
	double*   a1;
	double*   a2;
	dcomplex* g11;
	int       j, g, k;
	int       nG, nG_app;
	int       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
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
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asd( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asd( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asd( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asc_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A )
{
	float     one  = bli_s1();
	float     zero = bli_s0();
	float     gamma;
	float     sigma;
	scomplex* a1;
	scomplex* a2;
	scomplex* g11;
	int       j, g, k;
	int       nG, nG_app;
	int       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_asc_var1( k_G,
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
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asc( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asc( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asc( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Apply_G_rf_asz_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A )
{
	double    one  = bli_d1();
	double    zero = bli_d0();
	double    gamma;
	double    sigma;
	dcomplex* a1;
	dcomplex* a2;
	dcomplex* g11;
	int       j, g, k;
	int       nG, nG_app;
	int       k_minus_1;

	k_minus_1 = k_G - 1;
	nG        = n_A - 1;

	// Use the simple variant for nG < 2(k - 1).
	if ( nG < k_minus_1 || k_G == 1 )
	{
		FLA_Apply_G_rf_asz_var1( k_G,
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
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Pipeline stage

	for ( j = k_minus_1; j < nG; ++j )
	{
		nG_app = k_G;

		for ( k = 0, g = j; k < nG_app; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	// Shutdown stage

	for ( j = nG - k_minus_1; j < nG; ++j )
	{
		nG_app = nG - j;

		for ( k = k_G - nG_app, g = nG - 1; k < k_G; ++k, --g )
		{
			g11   = buff_G + (g    )*rs_G + (k  )*cs_G;
			a1    = buff_A + (g    )*cs_A;
			a2    = buff_A + (g + 1)*cs_A;

			gamma = g11->real;
			sigma = g11->imag;

			// Skip the current iteration if the rotation is identity.
			if ( gamma == one && sigma == zero ) continue;

			MAC_Apply_G_mx2_asz( m_A,
			                     &gamma,
			                     &sigma,
			                     a1, 1,
			                     a2, 1 );
		}
	}

	return FLA_SUCCESS;
}

