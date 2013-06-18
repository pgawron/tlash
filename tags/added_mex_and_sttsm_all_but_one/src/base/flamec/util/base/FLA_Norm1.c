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

FLA_Error FLA_Norm1( FLA_Obj A, FLA_Obj norm )
{
  FLA_Obj AL,   AR,       A0,  a1,  A2;

  FLA_Obj b;
  FLA_Obj bL,   bR,       b0,  beta1,  b2;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Norm1_check( A, norm );

  FLA_Obj_create( FLA_Obj_datatype( A ), 1, FLA_Obj_width( A ), 0, 0, &b );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( b,    &bL,  &bR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( bL,  /**/ bR,        &b0, /**/ &beta1, &b2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Asum( a1, beta1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &bL,  /**/ &bR,        b0, beta1, /**/ b2,
                              FLA_LEFT );
  }

  FLA_Max_abs_value( b, norm );

  FLA_Obj_free( &b );

  return FLA_SUCCESS;
}

