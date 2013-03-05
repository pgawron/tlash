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

FLA_Error FLA_CAQR_UT_inc_factorize_panels( dim_t nb_part, FLA_Obj A, FLA_Obj TW )
{
  FLA_Obj AT,              A0, 
          AB,              A1,
                           A2;

  FLA_Obj TWT,             TW0, 
          TWB,             TW1,
                           TW2;

  dim_t b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( TW,   &TWT, 
                      &TWB,           0, FLA_TOP );

  while ( FLA_Obj_length( AB ) > 0 ){

    b = min( nb_part, FLA_Obj_length( AB ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( TWT,               &TW0, 
                        /* ** */            /* ** */
                                              &TW1, 
                           TWB,               &TW2,       b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // Perform an incremental QR factorization on A1, writing triangular
    // block Householder factors to T in TW1.
    FLASH_QR_UT_inc( A1, TW1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,               A0, 
                                                 A1, 
                            /* ** */          /* ** */
                              &AB,               A2,      FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &TWT,              TW0, 
                                                 TW1, 
                            /* ** */          /* ** */
                              &TWB,              TW2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

