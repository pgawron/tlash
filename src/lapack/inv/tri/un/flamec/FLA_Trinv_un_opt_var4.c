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

FLA_Error FLA_Trinv_un_opt_var4( FLA_Obj A )
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

      FLA_Trinv_un_ops_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_Trinv_un_opd_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_Trinv_un_opc_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_Trinv_un_opz_var4( mn_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_ops_var4( int mn_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_sscalv( BLIS_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );
    bli_strsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a01, a12t, A02 );
    bli_sger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              mn_behind,
              mn_ahead,
              buff_m1,
              a01, rs_A,
              a12t, cs_A,
              A02, rs_A, cs_A );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_strmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_sinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opd_var4( int mn_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_dscalv( BLIS_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );
    bli_dtrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a01, a12t, A02 );
    bli_dger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              mn_behind,
              mn_ahead,
              buff_m1,
              a01, rs_A,
              a12t, cs_A,
              A02, rs_A, cs_A );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_dtrmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_dinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opc_var4( int mn_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_cscalv( BLIS_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );
    bli_ctrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a01, a12t, A02 );
    bli_cger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              mn_behind,
              mn_ahead,
              buff_m1,
              a01, rs_A,
              a12t, cs_A,
              A02, rs_A, cs_A );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_ctrmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_cinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Trinv_un_opz_var4( int mn_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    int       mn_ahead  = mn_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Scal_external( FLA_MINUS_ONE, a12t );
    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A22, a12t );
    bli_zscalv( BLIS_NO_CONJUGATE,
                mn_ahead,
                buff_m1,
                a12t, cs_A );
    bli_ztrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_ahead,
               A22, rs_A, cs_A,
               a12t, cs_A );

    // FLA_Ger_external( FLA_MINUS_ONE, a01, a12t, A02 );
    bli_zger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              mn_behind,
              mn_ahead,
              buff_m1,
              a01, rs_A,
              a12t, cs_A,
              A02, rs_A, cs_A );

    // FLA_Trmv_external( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a01 );
    bli_ztrmv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Invert( FLA_NO_CONJUGATE, alpha11 );
    bli_zinverts( BLIS_NO_CONJUGATE,
                  alpha11 );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
