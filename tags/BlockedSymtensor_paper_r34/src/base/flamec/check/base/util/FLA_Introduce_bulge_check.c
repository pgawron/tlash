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

FLA_Error FLA_Introduce_bulge_check( FLA_Obj shift, FLA_Obj gamma, FLA_Obj sigma, FLA_Obj delta1, FLA_Obj epsilon1, FLA_Obj delta2, FLA_Obj beta, FLA_Obj epsilon2 )
{
  FLA_Error e_val;

  e_val = FLA_Check_nonconstant_object( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, epsilon1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, delta2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, epsilon2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( epsilon1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( epsilon2 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

