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

//NOTE: FIX for breaking/retaining symmetry
//psttm: Partial symmetric tensor times matrix
FLA_Error FLA_Psttm( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
    FLA_Obj BT, BB;
    FLA_Obj B0, B1, B2;

    FLA_Obj CT, CB;
    FLA_Obj C0, C1, C2;

    FLA_Part_1xmode2(B, &BT,
                     &BB, 0, 0, FLA_TOP);
    FLA_Part_1xmode2(C, &CT,
                     &CB, mode, 0, FLA_TOP);

    while(FLA_Obj_dimsize(CT, mode) < FLA_Obj_dimsize(C, mode)){
        dim_t b = 1;
        FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
                                      /**/ /**/
                                          &B1,
                                      BB, &B2, 0, b, FLA_BOTTOM);
        FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
                                      /**/ /**/
                                          &C1,
                                      CB, &C2, mode, b, FLA_BOTTOM);
        /*********************************/
        FLA_Psttv(alpha, A, mode, beta, B1, C1);
        /*********************************/
        FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
                                               C1,
                                        /********/
                                          &CB, C2, mode, FLA_TOP);
        FLA_Cont_with_1xmode3_to_1xmode2( &BT, B0,
                                               B1,
                                        /********/
                                          &BB, B2, 0, FLA_TOP);
    }

	return FLA_SUCCESS;
}

