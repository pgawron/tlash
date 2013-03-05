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

//TODO: THIS FUNCTION PROKEN!!!
FLA_Error FLA_Random_unitary_matrix( FLA_Obj A )
{
  printf("FLA_Random_unitary_matrix is non-functioning as LAPACK routines have been removed\n");
  return FLA_SUCCESS;

//Original libFLAME lines
/*

  FLA_Obj B, T;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Random_unitary_matrix_check( A );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &B );

  FLA_Random_matrix( B );

  FLA_QR_UT_create_T( B, &T );

  FLA_QR_UT( B, T );

  FLA_QR_UT_form_Q( B, T, A );
  //FLA_Apply_Q_UT_create_workspace( A, T, &W );
  //FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, B, T, W, A );

  //FLA_Obj_free( &W );
  FLA_Obj_free( &T );
  FLA_Obj_free( &B );
*/

/*
  FLA_Datatype datatype;
  FLA_Obj      v, tau;
  FLA_Obj      aT,
               AB;
  int          i, mn;
  int          k;

  datatype = FLA_Obj_datatype( A );
  mn       = FLA_Obj_length( A );
  k        = 1;

  FLA_Obj_create( datatype, mn-1, 1, 0, 0, &v );
  FLA_Obj_create( datatype, 1,    1, 0, 0, &tau );

  FLA_Obj_set_to_identity( A );

  FLA_Part_2x1( A,   &aT,
                     &AB,    1, FLA_TOP );

  for ( i = 0; i < k; ++i )
  {
    FLA_Random_matrix( tau );
    FLA_Random_matrix( v );

    FLA_Apply_H2_UT( FLA_LEFT, tau, v, aT,
                                       AB );
  }

  FLA_Obj_free( &tau );
  FLA_Obj_free( &v );
*/

  return FLA_SUCCESS;
}

