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

FLA_Error FLA_Trinv_un_opt_var2( FLA_Obj A )
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

      FLA_Trinv_un_ops_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_un_opd_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_un_opc_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_un_opz_var2( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_ops_var2( int mn_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float     alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_strsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Inv_scal_external( alpha11, a12t );
    bli_sneg2( alpha11, &alpha11_m1 );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a12t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_sinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opd_var2( int mn_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double    alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_dtrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Inv_scal_external( alpha11, a12t );
    bli_dneg2( alpha11, &alpha11_m1 );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a12t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_dinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opc_var2( int mn_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex  alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_ctrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Inv_scal_external( alpha11, a12t );
    bli_cneg2( alpha11, &alpha11_m1 );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a12t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_cinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opz_var2( int mn_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex  alpha11_m1;
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_ztrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Inv_scal_external( alpha11, a12t );
    bli_zneg2( alpha11, &alpha11_m1 );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   &alpha11_m1,
                   a12t, cs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_zinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
