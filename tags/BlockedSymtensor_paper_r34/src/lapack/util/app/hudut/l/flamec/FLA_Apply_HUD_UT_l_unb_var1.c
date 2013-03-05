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

FLA_Error FLA_Apply_HUD_UT_lh_unb_var1( FLA_Obj tau, FLA_Obj w12t,
                                                     FLA_Obj r12t,
                                        FLA_Obj u1,  FLA_Obj C2,
                                        FLA_Obj v1,  FLA_Obj D2 )
/*
  Apply a single up-and-downdating Householder transform H from the left to a
  row vector r12t and matrices C2 and D2:

    / r12t \        / r12t \
    |  C2  |  :=  H |  C2  |                                                (1)
    \  D2  /        \  D2  /

  H is defined as:

          / / 1 0 0 \            / 1 0  0 \  / 1  \ ( 1  u1' v1' ) \
    H  =  | | 0 I 0 | - inv(tau) | 0 I  0 |  | u1 |                |        (2)
          \ \ 0 0 I /            \ 0 0 -I /  \ v1 /                /

  Substituting (2) into (1), we have:

    / r12t \        / r12t \
    |  C2  |  :=  H |  C2  |
    \  D2  /        \  D2  /
                  / / 1 0 0 \            / 1 0  0 \  / 1  \ ( 1  u1' v1' ) \' / r12t \
               =  | | 0 I 0 | - inv(tau) | 0 I  0 |  | u1 |                |  |  C2  |
                  \ \ 0 0 I /            \ 0 0 -I /  \ v1 /                /  \  D2  /
                  / / 1 0 0 \            / 1 0  0 \  / 1    u1'     v1'   \ \  / r12t \
               =  | | 0 I 0 | - inv(tau) | 0 I  0 |  | u1  u1 u1'  u1 v1' | |  |  C2  |
                  \ \ 0 0 I /            \ 0 0 -I /  \ v1  v1 u1'  v1 v1' / /  \  D2  /
                  / / 1 0 0 \            /  1     u1'      v1'   \ \  / r12t \
               =  | | 0 I 0 | - inv(tau) |  u1   u1 u1'   u1 v1' | |  |  C2  |
                  \ \ 0 0 I /            \ -v1  -v1 u1'  -v1 v1' / /  \  D2  /
                  / r12t \            /  1     u1'      v1'   \ / r12t \
               =  |  C2  | - inv(tau) |  u1   u1 u1'   u1 v1' | |  C2  |
                  \  D2  /            \ -v1  -v1 u1'  -v1 v1' / \  D2  /
                  / r12t \   /  inv(tau)      inv(tau) u1'      inv(tau) v1'    \ / r12t \
               =  |  C2  | - |  inv(tau) u1   inv(tau) u1 u1'   inv(tau) u1 v1' | |  C2  |
                  \  D2  /   \ -inv(tau) v1  -inv(tau) v1 u1'  -inv(tau) v1 v1' / \  D2  /
                  / r12t \   /  inv(tau) r12t    + inv(tau) u1' C2    + inv(tau) v1' D2    \
               =  |  C2  | - |  inv(tau) u1 r12t + inv(tau) u1 u1' C2 + inv(tau) u1 v1' D2 |
                  \  D2  /   \ -inv(tau) v1 r12t - inv(tau) v1 u1' C2 - inv(tau) v1 v1' D2 /

  Thus, r12t is updated as:

      r12t    :=    r12t  -  inv(tau) r12t  +  inv(tau) u1' C2  +  inv(tau) v1' D2
               =    r12t  -  inv(tau) ( r12t  +  u1' C2  +  v1' D2 )

  C2 is updated as:

      C2      :=     C2   -  inv(tau) u1 r12t  +  inv(tau) u1 u1' C2  +  inv(tau) u1 v1' D2
               =     C2   -  inv(tau) u1 ( r12t  +  u1' C2  +  v1' D2 )

  And D2 is updated as:

      D2      :=     D2   -  ( -inv(tau) v1 r12t  -  inv(tau) v1 u1' C2  -  inv(tau) v1 v1' D2 )
               =     D2   +  inv(tau) v1 r12t  +  inv(tau) v1 u1' C2  +  inv(tau) v1 v1' D2
               =     D2   +  inv(tau) v1 ( r12t  +  u1' C2  +  v1' D2 )

  Note that:

    inv(tau) ( r12t  +  u1' C2  +  v1' D2 )

  is common to both updates, and thus may be computed and stored in
  workspace, and then re-used.
 
  -FGVZ
*/
{
  if ( FLA_Obj_has_zero_dim( r12t ) ) return FLA_SUCCESS;

  // w12t = r12t;
  FLA_Copy_external( r12t, w12t );

  // w12t = w12t + u1' * C2;
  //      = w12t + C2^T * conj(u1);
  FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, C2, u1, FLA_ONE, w12t );

  // w12t = w12t + v1' * D2;
  //      = w12t + D2^T * conj(v1);
  FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, D2, v1, FLA_ONE, w12t );

  // w12t = w12t / tau;
  FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w12t );

  // r12t = - w12t + r12t;
  FLA_Axpy_external( FLA_MINUS_ONE, w12t, r12t );

  // C2 = - u1 * w12t + C2;
  FLA_Ger_external( FLA_MINUS_ONE, u1, w12t, C2 );

  // D2 = v1 * w12t + D2;
  FLA_Ger_external( FLA_ONE, v1, w12t, D2 );

  return FLA_SUCCESS;
}

