#include "FLAME.h"

#define AA( i,j ) buff_A[ (j)*ldim_A + (i) ]
#define xx( i )   buff_x[ i ]
#define vv( i )   buff_v[ i ]
#define ww( i )   buff_w[ i ]

/*
       Effective computation:
       v = A' * x;
       w = A * x;
 */

int NoFLA_Atx_Ax_var1_opt2(
       int m, int n, double * buff_A, int ldim_A, 
       double * buff_x, double * buff_v, double * buff_w ) {
  int     j, i_one = 1;
  double  ddot_();

  //// int m, n, ldA, j, i_one=1;
  //// double *buff_A, *buff_x, *buff_v, *buff_w, ddot_();
  //// m = FLA_Obj_length( A );
  //// n = FLA_Obj_width ( A );
  //// ldA = FLA_Obj_ldim( A );
  //// buff_A = ( double * ) FLA_Obj_buffer_at_view( A );  
  //// buff_x = ( double * ) FLA_Obj_buffer_at_view( x );
  //// buff_v = ( double * ) FLA_Obj_buffer_at_view( v );
  //// buff_w = ( double * ) FLA_Obj_buffer_at_view( w );

  //// printf( "NoFLA_Atx_Ax_var1_opt2\n" );
 
  //// MyFLA_Obj_set_to_zero( v );
  //// MyFLA_Obj_set_to_zero( w );
  for ( j = 0; j < n; j++ ) {
    vv( j ) = 0.0;
    ww( j ) = 0.0;
  }

  for ( j = 0; j < n; j++ ) {
    vv( j ) =  ddot_( &m, &AA( 0,j ), &i_one, buff_x, &i_one );
    daxpy_( &m, &xx( j ), &AA( 0,j ), &i_one, buff_w, &i_one );
  }  

  return FLA_SUCCESS;
}

