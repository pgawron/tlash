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

extern FLA_Bool FLA_Gemm_initialized;

FLA_Error FLA_Gemm_pp_nn_var1( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj C, int nb_alg )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj packed_C1;

  int b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  /* Initialize the FLA_Gemm() interface to the kernel environment
     Note: the blocksize given to the kernel environment can be non-
     square. We pass the m and n dimensions of the blocksize here. */

  FLA_Gemm_init( nb_alg, FLA_Obj_width( A ) );

  /* Pack B */
  /* Note: the idea here is that, optionally,
     -  B is packed, and/or
     -  B is scaled                                    
     If B needs not be packed, it is not packed.  If the multiplication
     by alpha happens elsewhere, no scaling occurs.
     The "NoTranspose" means that in the version of gemm being updated
     B is not transposed.  In packing, B could be transposed, if there
     is an advantage to this.  So, the "NoTranspose" means that input
     B is not transposed in the FLA_Gemm call.  */
 
  FLA_Gemm_pack_andor_scale_B( FLA_NO_TRANSPOSE, alpha, B );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    b = min( FLA_Obj_length( AB ), nb_alg );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &C1, 
                           CB,                &C2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /* C1 = alpha * A1 * B + C1; */

    /* Pack C1 */
    /* Note: the idea here is that, optionally,
       -  packed space is provided for computing packed_C1 = alpha * A1 * B 
       If C1 needs not be packed, then C can just be returned by this routine. */
    FLA_Gemm_pack_C( FLA_NO_TRANSPOSE, C1 );

    /* Pack A */
    /* Note: the idea here is that, optionally,
        -  A is packed, and/or
        -  A is scaled                                    
	If A needs not be packed, it is not packed.  If the multiplication
	by alpha happens elsewhere, no scaling occurs. */
    FLA_Gemm_pack_andor_scale_A( FLA_NO_TRANSPOSE, alpha, A1 );

    /* Call the kernel routine */
    FLA_Gemm_kernel( alpha, A1, B, C1 );

    /* Unpack C1 */
    /* Note: the idea here is that, optionally,
       -  packed_C1 is added to C1, possibly scaled at this point. */
    FLA_Gemm_unpack_andor_scale_C( FLA_NO_TRANSPOSE, alpha, C1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  A1, 
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );

  }

  /* Release the space used to pack A1 */
  /* Note: notice that the space provided for A1 can be recycled 
     everytime through the loop, which is why this call is outside the loop. 
     If the space is statically allocated, or A1 was not packed,
     this could be a no-op. */
  FLA_Gemm_release_pack_A( FLA_NO_TRANSPOSE, A1 );

  /* Release the space used to pack B */
  /* Note: If the space is statically allocated, or B was not packed,
     this could be a no-op. */
  FLA_Gemm_release_pack_B( FLA_NO_TRANSPOSE, B );

  FLA_Gemm_finish();

  return FLA_SUCCESS;
}
