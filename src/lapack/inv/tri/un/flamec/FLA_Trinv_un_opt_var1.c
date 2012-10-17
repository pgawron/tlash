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

FLA_Error FLA_Trinv_un_opt_var1( FLA_Obj A )
{
  FLA_Datatype datatype;
  int          mn_A;
  int          rs_A, cs_A;

  datatype = FLA_Obj_datatype( A );

  mn_A     = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );

      FLA_Trinv_un_ops_var1( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_un_opd_var1( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_un_opc_var1( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_un_opz_var1( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_ops_var1( int mn_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float     alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_strmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a01 );
    // FLA_Inv_scal_external( alpha11, a01 );
    bli_sneg2( alpha11, &alpha11_m1 );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   mn_behind,
                   &alpha11_m1,
                   a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_sinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opd_var1( int mn_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double    alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_dtrmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a01 );
    // FLA_Inv_scal_external( alpha11, a01 );
    bli_dneg2( alpha11, &alpha11_m1 );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   mn_behind,
                   &alpha11_m1,
                   a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_dinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opc_var1( int mn_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex  alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_ctrmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a01 );
    // FLA_Inv_scal_external( alpha11, a01 );
    bli_cneg2( alpha11, &alpha11_m1 );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   mn_behind,
                   &alpha11_m1,
                   a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_cinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opz_var1( int mn_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex  alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_ztrmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a01 );
    // FLA_Inv_scal_external( alpha11, a01 );
    bli_zneg2( alpha11, &alpha11_m1 );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   mn_behind,
                   &alpha11_m1,
                   a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_zinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
