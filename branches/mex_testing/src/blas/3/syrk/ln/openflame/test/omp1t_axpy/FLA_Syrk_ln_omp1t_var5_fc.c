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
#include "FLA_Syrk_ln_omp.h"

FLA_Error FLA_Syrk_ln_omp1t_var5_fc( FLA_Obj A, FLA_Obj C, int nb_alg )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;
  FLA_Obj MyC;

  int b;
  
  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    nb_alg = FLA_Obj_width( A )/omp_get_num_threads() + 1;
    b = min( FLA_Obj_width( AR ), nb_alg );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    #pragma intel omp task captureprivate(A1) private(MyC)
    {
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &MyC );
    FLA_Obj_set_to_zero( MyC );
    
    /* MyC := A1 * A1' */
    FLA_Syrk( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_ONE, A1, FLA_ZERO, MyC );

    /* C := MyC */
    FLA_Axpy_sync_circular( FLA_ONE, MyC, C );

    FLA_Obj_free( &MyC );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );
  }
  }

  return FLA_SUCCESS;
}

