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

FLA_Error FLA_Norm_inf( FLA_Obj A, FLA_Obj norm )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj bT,              b0,
          bB,              beta1,
                           b2;
  FLA_Obj b;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Norm_inf_check( A, norm );

  FLA_Obj_create( FLA_Obj_datatype( A ), FLA_Obj_length( A ), 1, 0, 0, &b );

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( b,    &bT, 
                      &bB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* *** */
                                              &a1t, 
                           AB,                &A2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( bT,                &b0, 
                        /* ** */            /* ***** */
                                              &beta1, 
                           bB,                &b2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Asum( a1t, beta1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  a1t, 
                            /* ** */           /* *** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &bT,                b0, 
                                                  beta1, 
                            /* ** */           /* ***** */
                              &bB,                b2,     FLA_TOP );

  }

  FLA_Max_abs_value( b, norm );

  FLA_Obj_free( &b );

  return FLA_SUCCESS;
}

