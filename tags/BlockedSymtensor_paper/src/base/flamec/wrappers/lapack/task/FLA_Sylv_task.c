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

extern fla_sylv_t* fla_sylv_cntl_leaf;

FLA_Error FLA_Sylv_task( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  return FLA_Sylv_internal( transa, transb,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_nn_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_nh_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_hn_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

FLA_Error FLA_Sylv_hh_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  //return FLA_Sylv_unb_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, isgn, A, B, C, scale );
  return FLA_Sylv_internal( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                            isgn, A, B, C, scale,
                            fla_sylv_cntl_leaf );
}

