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

FLA_Error FLA_Trinv_lu_blk_var1( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_lu_blk_var2( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_lu_blk_var3( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_lu_blk_var4( FLA_Obj A, fla_trinv_t* cntl );

FLA_Error FLA_Trinv_lu_unb_var1( FLA_Obj A );
FLA_Error FLA_Trinv_lu_unb_var2( FLA_Obj A );
FLA_Error FLA_Trinv_lu_unb_var3( FLA_Obj A );
FLA_Error FLA_Trinv_lu_unb_var4( FLA_Obj A );

FLA_Error FLA_Trinv_lu_opt_var1( FLA_Obj A );
FLA_Error FLA_Trinv_lu_ops_var1( int mn_A,
                                 float*    A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opd_var1( int mn_A,
                                 double*   A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opc_var1( int mn_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opz_var1( int mn_A,
                                 dcomplex* A, int rs_A, int cs_A );

FLA_Error FLA_Trinv_lu_opt_var2( FLA_Obj A );
FLA_Error FLA_Trinv_lu_ops_var2( int mn_A,
                                 float*    A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opd_var2( int mn_A,
                                 double*   A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opc_var2( int mn_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opz_var2( int mn_A,
                                 dcomplex* A, int rs_A, int cs_A );

FLA_Error FLA_Trinv_lu_opt_var3( FLA_Obj A );
FLA_Error FLA_Trinv_lu_ops_var3( int mn_A,
                                 float*    A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opd_var3( int mn_A,
                                 double*   A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opc_var3( int mn_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opz_var3( int mn_A,
                                 dcomplex* A, int rs_A, int cs_A );

FLA_Error FLA_Trinv_lu_opt_var4( FLA_Obj A );
FLA_Error FLA_Trinv_lu_ops_var4( int mn_A,
                                 float*    A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opd_var4( int mn_A,
                                 double*   A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opc_var4( int mn_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_Trinv_lu_opz_var4( int mn_A,
                                 dcomplex* A, int rs_A, int cs_A );
