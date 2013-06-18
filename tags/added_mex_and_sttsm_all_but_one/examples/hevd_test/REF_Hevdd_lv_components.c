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

FLA_Error REF_Hevdd_lv_components( FLA_Obj A, FLA_Obj l,
                                   double* dtime_tred, double* dtime_tevd, double* dtime_appq )
/*
{
  return FLA_Hevdd_external( FLA_EVD_WITH_VECTORS, FLA_LOWER_TRIANGULAR, A, l );
}
*/
{
  FLA_Datatype dt_A;
  FLA_Datatype dt_A_real;
  FLA_Uplo     uplo;
  dim_t        m_A, n_A;
  FLA_Obj      t, d, e, W;
  double       dtime_temp;

  uplo      = FLA_LOWER_TRIANGULAR;
  dt_A      = FLA_Obj_datatype( A );
  dt_A_real = FLA_Obj_datatype_proj_to_real( A );
  m_A       = FLA_Obj_length( A );
  n_A       = FLA_Obj_width( A );

  FLA_Obj_create( dt_A,      m_A,   1, 0, 0, &t );
  FLA_Obj_create( dt_A_real, m_A,   1, 0, 0, &d );
  FLA_Obj_create( dt_A_real, m_A-1, 1, 0, 0, &e );
  FLA_Obj_create( dt_A,      m_A, n_A, 0, 0, &W );


  dtime_temp = FLA_Clock();
  {
    // Reduce to tridiagonal form.
    FLA_Tridiag_blk_external( uplo, A, t );
    FLA_Tridiag_UT_extract_diagonals( uplo, A, d, e );
  }
  *dtime_tred = FLA_Clock() - dtime_temp;


  dtime_temp = FLA_Clock();
  {
    // Divide-and-conquor algorithm.
    FLA_Tevdd_external( FLA_EVD_OF_TRIDIAG_WITH_VECTORS, d, e, W );
  }
  *dtime_tevd = FLA_Clock() - dtime_temp;


  dtime_temp = FLA_Clock();
  {
    // Apply Q.
    FLA_Tridiag_apply_Q_external( FLA_LEFT, uplo, FLA_NO_TRANSPOSE, A, t, W );
    FLA_Copy( W, A );
  }
  *dtime_appq = FLA_Clock() - dtime_temp;


  // Copy eigenvalues to output vector.
  FLA_Copy( d, l );

  // Sort eigenvalues and eigenvectors.
  FLA_Sort_evd( FLA_FORWARD, l, A );

  FLA_Obj_free( &t );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &W );

  return FLA_SUCCESS;
}

