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

FLA_Error FLA_Apply_H2_UT_r_unb_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 )
/*
  Apply a single Householder transform H  from the right to a column vector a1
  and a matrix A2:

    ( a1  A2 ) := ( a1  A2 ) H                                              (1)

  H is defined as:

    H  =  / I - inv(tau) / 1  \ ( 1  u2' ) \                                (2)
          \              \ u2 /            /

  Typically, a1 and A2 are horizontally adjacent views into a larger matrix,
  though this is not always the case.

  Substituting (2) into (1), we have:

    ( a1  A2 ) := ( a1  A2 ) H
                = ( a1  A2 ) / I - inv(tau) / 1  \ ( 1  u2' ) \
                             \              \ u2 /            /
                = ( a1  A2 ) / I - inv(tau) / 1     u2'  \ \
                             \              \ u2  u2 u2' / /
                = ( a1  A2 ) - inv(tau) ( a1  A2 ) / 1     u2'  \
                                                   \ u2  u2 u2' /
                = ( a1  A2 ) - inv(tau) ( a1 + A2 u2   a1 u2' + A2 u2 u2' )
                = ( a1  A2 ) - ( inv(tau) ( a1 + A2 u2 )   inv(tau) ( a1 u2' + A2 u2 u2' ) )

  Thus, a1t is updated as:

     a1        :=   a1  -  inv(tau) ( a1 + A2 u2 )

  And A2 is updated as:

     A2        :=   A2  -  inv(tau) ( a1 u2' + A2 u2 u2' )
                =   A2  -  inv(tau) ( a1 + A2 u2 ) u2'

  Note that:

    inv(tau) ( a1 + A2 u2 )

  is common to both updates, and thus may be computed and stored in
  workspace, and then re-used.
 
  -FGVZ
*/
{
  FLA_Obj w1;

  if ( FLA_Obj_has_zero_dim( a1 ) ) return FLA_SUCCESS;

  // w1 = a1;
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  FLA_Copy_external( a1, w1 );

  // w1 = w1 + A2 * u2;
  FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );

  // w1 = w1 / tau;
  FLA_Inv_scal_external( tau, w1 );

  // a1 = a1 - w1;
  FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );

  // A2 = A2 - w1 * u2';
  FLA_Gerc_external( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );

  FLA_Obj_free( &w1 );

  return FLA_SUCCESS;
}
