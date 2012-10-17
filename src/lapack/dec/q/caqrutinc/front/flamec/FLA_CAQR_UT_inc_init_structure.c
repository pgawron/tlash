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

void FLA_CAQR_UT_inc_init_structure( dim_t p, dim_t nb_part, FLA_Obj R )
{
	dim_t    m, n;
	dim_t    rs, cs;
	dim_t    i, j, ip;
	FLA_Obj* buff_R;

	m      = FLA_Obj_length( R );
	n      = FLA_Obj_width( R );
	rs     = FLA_Obj_row_stride( R );
	cs     = FLA_Obj_col_stride( R );
	buff_R = FLA_Obj_buffer_at_view( R );

	// Fill in R by row panels.
	for ( ip = 0; ip < p; ++ip )
	{
		FLA_Obj* buff_R1 = buff_R + (ip*nb_part)*rs;

		int  m_behind   = ip*nb_part;
		int  m_ahead    = m - m_behind;

		int  m_cur      = min( nb_part, m_ahead );
		int  n_cur      = n;

		// Iterate across columns for the current panel.
		for ( j = 0; j < n_cur; ++j )
		{
			FLA_Obj* rho = buff_R1 + j*cs;

			// Mark the above-diagonal blocks as full.
			for ( i = 0; i < j; ++i )
			{
				rho->base->uplo = FLA_FULL_MATRIX;
				rho += rs;
			}

			// Mark the diagonal block as triangular.
			rho->base->uplo = FLA_UPPER_TRIANGULAR;
			rho += rs;
			
			// Mark the below-diagonal blocks as zero.
			for ( i = j + 1; i < m_cur; ++i )
			{
				rho->base->uplo = FLA_ZERO_MATRIX;
				rho += rs;
			}
		}
	}
}

