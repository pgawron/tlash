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

FLA_Error FLA_Chol_u_opt_var2( FLA_Obj A )
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

      FLA_Chol_u_ops_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Chol_u_opd_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Chol_u_opc_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Chol_u_opz_var2( mn_A,
                           buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_u_ops_var2( int mn_A,
                               float* buff_A, int rs_A, int cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bli_sdots( BLIS_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bli_sgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // value = FLA_Sqrt( alpha11 );
    // if ( value != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bli_ssqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_u_opd_var2( int mn_A,
                               double* buff_A, int rs_A, int cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bli_ddots( BLIS_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bli_dgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // value = FLA_Sqrt( alpha11 );
    // if ( value != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bli_dsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_u_opc_var2( int mn_A,
                               scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bli_cdots( BLIS_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bli_cgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // value = FLA_Sqrt( alpha11 );
    // if ( value != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bli_csqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_u_opz_var2( int mn_A,
                               dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bli_zdots( BLIS_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bli_zgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // value = FLA_Sqrt( alpha11 );
    // if ( value != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bli_zsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

