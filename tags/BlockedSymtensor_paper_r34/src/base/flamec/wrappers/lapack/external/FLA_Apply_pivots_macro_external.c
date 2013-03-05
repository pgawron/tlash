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

FLA_Error FLA_Apply_pivots_macro_external( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A )
{
   int          i, j;
   int          ipiv;
   int*         buf_p    = ( int* ) FLA_Obj_buffer_at_view( p );
   FLA_Obj*     blocks   = FLASH_OBJ_PTR_AT( A );
   int          m_blocks = FLA_Obj_length( A );
   int          m_A      = FLA_Obj_length( *blocks );
   int          n_A      = FLA_Obj_width( *blocks );
   FLA_Datatype datatype = FLA_Obj_datatype( A );

#ifdef FLA_ENABLE_WINDOWS_BUILD
   int* m  = ( int* ) _alloca( m_blocks * sizeof( int ) );
   int* cs = ( int* ) _alloca( m_blocks * sizeof( int ) );
#else
   int* m  = ( int* ) malloc( m_blocks * sizeof( int ) );
   int* cs = ( int* ) malloc( m_blocks * sizeof( int ) );
   //int m[m_blocks];
   //int cs[m_blocks];
#endif

   if ( side != FLA_LEFT || trans != FLA_NO_TRANSPOSE )
      FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

   switch ( datatype )
   {
      case FLA_FLOAT:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         float** buffer = ( float** ) _alloca( m_blocks * sizeof( float* ) );
#else
         float** buffer = ( float** ) malloc( m_blocks * sizeof( float* ) );
         //float*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( float* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bli_sswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
      case FLA_DOUBLE:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         double** buffer = ( double** ) _alloca( m_blocks * sizeof( double* ) );
#else
         double** buffer = ( double** ) malloc( m_blocks * sizeof( double* ) );
         //double*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( double* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bli_dswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
      case FLA_COMPLEX:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         scomplex** buffer = ( scomplex** ) _alloca( m_blocks * sizeof( scomplex* ) );
#else
         scomplex** buffer = ( scomplex** ) malloc( m_blocks * sizeof( scomplex* ) );
         //scomplex*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( scomplex* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bli_cswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
      case FLA_DOUBLE_COMPLEX:
      {
#ifdef FLA_ENABLE_WINDOWS_BUILD
         dcomplex** buffer = ( dcomplex** ) _alloca( m_blocks * sizeof( dcomplex* ) );
#else
         dcomplex** buffer = ( dcomplex** ) malloc( m_blocks * sizeof( dcomplex* ) );
         //dcomplex*  buffer[m_blocks];
#endif
         for ( i = 0; i < m_blocks; i++ )
         {
            buffer[i] = ( dcomplex* ) FLA_Obj_buffer_at_view( blocks[i] );
            
            m[i] = FLA_Obj_length( blocks[i] );
            cs[i] = FLA_Obj_col_stride( blocks[i] );
         }
         
         for ( j = 0; j < m_A; j++ )
         {
            ipiv = buf_p[j] + j;
            
            if ( ipiv != j )
            {
               i = 0;
               
               while ( ipiv >= m[i] )
               {
                  ipiv = ipiv - m[i];
                  i++;
               }

               bli_zswapv( n_A,
                           buffer[0] + j,    cs[0],
                           buffer[i] + ipiv, cs[i] );
            }
         }
#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
         free( buffer ); 
#endif
         break;
      }
   }

#ifdef FLA_ENABLE_WINDOWS_BUILD
#else
   free( m );
   free( cs );
#endif

   return FLA_SUCCESS;
}

