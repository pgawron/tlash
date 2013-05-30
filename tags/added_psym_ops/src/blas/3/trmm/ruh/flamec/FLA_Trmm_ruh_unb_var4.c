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

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Trmm_ruh_unb_var4( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj BT,              B0,
          BB,              b1t,
                           B2;

  FLA_Scal_external( alpha, B );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( BB ) < FLA_Obj_length( B ) ){

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                                              &b1t, 
                        /* ** */            /* *** */
                           BB,                &B2,        1, FLA_TOP );

    /*------------------------------------------------------------*/

    /* b1t = b1t * triu( A )'; */
    FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diagA, A, b1t );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                            /* ** */           /* *** */
                                                  b1t, 
                              &BB,                B2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
