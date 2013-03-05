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

FLA_Error FLA_SA_LU_unb( FLA_Obj U, FLA_Obj D, FLA_Obj p, FLA_Obj L )
{
  FLA_Datatype datatype;
  int          m_U, cs_U;
  int          m_D, cs_D;
  int               cs_L;
  // int               rs_U;
  int               rs_D;
  // int               rs_L;
  int          m_U_min_j, m_U_min_j_min_1; 
  int          j, ipiv;
  int*         buff_p;

  if ( FLA_Obj_has_zero_dim( U ) ) return FLA_SUCCESS;
  
  datatype = FLA_Obj_datatype( U );

  m_U      = FLA_Obj_length( U );
  // rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );

  m_D      = FLA_Obj_length( D );
  rs_D     = FLA_Obj_row_stride( D );
  cs_D     = FLA_Obj_col_stride( D );
  
  // rs_L     = FLA_Obj_row_stride( L );
  cs_L     = FLA_Obj_col_stride( L );

  FLA_Copy_external( U, L );
  FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, L );

  buff_p     = ( int * ) FLA_INT_PTR( p );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_U      = ( float * ) FLA_FLOAT_PTR( U );
    float* buff_D      = ( float * ) FLA_FLOAT_PTR( D );
    float* buff_L      = ( float * ) FLA_FLOAT_PTR( L );
    float* buff_minus1 = ( float * ) FLA_FLOAT_PTR( FLA_MINUS_ONE );
    float  L_tmp;
    float  D_tmp;
    float  d_inv_Ljj;

    for ( j = 0; j < m_U; ++j )
    {
      bli_samax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabsf( L_tmp ) < fabsf( D_tmp ) )
      {
        bli_sswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      d_inv_Ljj = 1.0F / buff_L[ j*cs_L + j ];

      bli_sscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bli_sger( BLIS_NO_CONJUGATE,
                  BLIS_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bli_scopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_U      = ( double * ) FLA_DOUBLE_PTR( U );
    double* buff_D      = ( double * ) FLA_DOUBLE_PTR( D );
    double* buff_L      = ( double * ) FLA_DOUBLE_PTR( L );
    double* buff_minus1 = ( double * ) FLA_DOUBLE_PTR( FLA_MINUS_ONE );
    double  L_tmp;
    double  D_tmp;
    double  d_inv_Ljj;

    for ( j = 0; j < m_U; ++j )
    {
      bli_damax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabs( L_tmp ) < fabs( D_tmp ) )
      {
        bli_dswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      d_inv_Ljj = 1.0 / buff_L[ j*cs_L + j ];

      bli_dscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bli_dger( BLIS_NO_CONJUGATE,
                  BLIS_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bli_dcopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex* buff_U      = ( scomplex * ) FLA_COMPLEX_PTR( U );
    scomplex* buff_D      = ( scomplex * ) FLA_COMPLEX_PTR( D );
    scomplex* buff_L      = ( scomplex * ) FLA_COMPLEX_PTR( L );
    scomplex* buff_minus1 = ( scomplex * ) FLA_COMPLEX_PTR( FLA_MINUS_ONE );
    scomplex  L_tmp;
    scomplex  D_tmp;
    scomplex  d_inv_Ljj;
    scomplex  Ljj;
    float     temp;

    for ( j = 0; j < m_U; ++j )
    {
      bli_camax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabsf( L_tmp.real + L_tmp.imag ) < fabsf( D_tmp.real + D_tmp.imag ) )
      {
        bli_cswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      Ljj = buff_L[ j*cs_L + j ];

      // d_inv_Ljj = 1.0 / Ljj
      temp = 1.0F / ( Ljj.real * Ljj.real +
                      Ljj.imag * Ljj.imag );
      d_inv_Ljj.real = Ljj.real *  temp;
      d_inv_Ljj.imag = Ljj.imag * -temp;

      bli_cscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bli_cger( BLIS_NO_CONJUGATE,
                  BLIS_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bli_ccopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_U      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( U );
    dcomplex* buff_D      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( D );
    dcomplex* buff_L      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( L );
    dcomplex* buff_minus1 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
    dcomplex  L_tmp;
    dcomplex  D_tmp;
    dcomplex  d_inv_Ljj;
    dcomplex  Ljj;
    double    temp;

    for ( j = 0; j < m_U; ++j )
    {
      bli_zamax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabs( L_tmp.real + L_tmp.imag ) < fabs( D_tmp.real + D_tmp.imag ) )
      {
        bli_zswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      Ljj = buff_L[ j*cs_L + j ];

      // d_inv_Ljj = 1.0 / Ljj
      temp = 1.0  / ( Ljj.real * Ljj.real +
                      Ljj.imag * Ljj.imag );
      d_inv_Ljj.real = Ljj.real *  temp;
      d_inv_Ljj.imag = Ljj.imag * -temp;

      bli_zscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bli_zger( BLIS_NO_CONJUGATE,
                  BLIS_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bli_zcopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  }

  return FLA_SUCCESS;
}
