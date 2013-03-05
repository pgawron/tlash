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

FLA_Error REF_Gemm_nn( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Datatype datatype;
  int          m, k, n, ldim_A, ldim_B, ldim_C;

  datatype = FLA_Obj_datatype( A );
  ldim_A   = FLA_Obj_ldim( A );
  ldim_B   = FLA_Obj_ldim( B );
  ldim_C   = FLA_Obj_ldim( C );
  m        = FLA_Obj_length( A );
  k        = FLA_Obj_width( A );
  n        = FLA_Obj_width( B );
  
  switch( datatype ){
    case FLA_DOUBLE:
    {
      double *buff_A, *buff_B, *buff_C, d_one=1.0;

      buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
      buff_B = ( double * ) FLA_Obj_buffer_at_view( B );
      buff_C = ( double * ) FLA_Obj_buffer_at_view( C );
    
      FLA_C2F( dgemm )( "N", "N", &m, &n, &k,
                        &d_one, buff_A, &ldim_A, buff_B, &ldim_B,
                        &d_one, buff_C, &ldim_C );
    } break;
  }
  
  return 0;
}

