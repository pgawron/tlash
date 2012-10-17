#include "FLAME.h"


#define AA( i,j ) buff_A[ (j)*ldA + (i) ]
#define uu( i ) buff_u[ i ]
#define vv( i ) buff_v[ i ]
#define ww( i ) buff_w[ i ]
#define xx( i ) buff_x[ i ]
#define yy( i ) buff_y[ i ]
#define zz( i ) buff_z[ i ]


/*
       Effective computation:
       A = A + beta * ( u * y' + z * u' );
       v = A' * x;
       w = A * x;
 */

int NoFLA_Ger2_Atx_Ax_var1_opt3( 
        int m, int n, double * buff_A, int ldA, double * beta,
        double * buff_u, double * buff_y, double * buff_z, double * buff_x, 
        double * buff_v, double * buff_w ) {
  int     i, j, i_one = 1;
  double  temp1, temp2, temp3, ddot_();
  double  * buff_A_col_j, * pA, * pu, * pz;
  /*
  double *buff_A, *buff_u, *buff_y, *buff_z, *buff_x, *buff_v, *buff_w, 
         ddot_();
  */
  
  //// MyFLA_Obj_set_to_zero( w );
  for ( i = 0; i < m; i++ ) {
    ww( i ) = 0.0;
  }

  //// m = FLA_Obj_length( A );
  //// n = FLA_Obj_width ( A );
  //// ldA = FLA_Obj_ldim( A );
  
  //// buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  //// buff_u = ( double * ) FLA_Obj_buffer_at_view( u );
  //// buff_y = ( double * ) FLA_Obj_buffer_at_view( y );
  //// buff_z = ( double * ) FLA_Obj_buffer_at_view( z );
  //// buff_x = ( double * ) FLA_Obj_buffer_at_view( x );
  //// buff_v = ( double * ) FLA_Obj_buffer_at_view( v );
  //// buff_w = ( double * ) FLA_Obj_buffer_at_view( w );

  //// printf( " pepe\n" );
 
  for ( j = 0; j < n; j++ ) {
    temp1 = *beta * yy( j );
    temp2 = *beta * buff_u[ j ];
    // daxpy_( &m, &temp1, buff_u, &i_one, &AA( 0,j ), &i_one );
    // daxpy_( &m, &temp2, buff_z, &i_one, &AA( 0, j ), &i_one );
    buff_A_col_j = &( buff_A[ 0 + ldA * j ] );
    pA = buff_A_col_j;
    pu = buff_u;
    pz = buff_z;
    for( i = 0; i < m; i++ ) {
      temp3 = temp1 * *pu;
      *pA += temp3 + temp2 * *pz;
      pA++;
      pu++;
      pz++;
    }
    //// daxpy_( &m, &uu( j ), buff_z, &i_one, &AA( 0,j ), &i_one );

    vv( j ) =  ddot_( &m, &AA( 0, j ), &i_one, buff_x, &i_one );
    //bli_ddot( BLIS_CONJUGATE,
    //          m,
    //          &AA(0,j), i_one,
    //          buff_x,   i_one,
    //          &vv(j) );
    daxpy_( &m, &xx( j ), &AA( 0, j ), &i_one, buff_w, &i_one );
    //bli_daxpyv( BLIS_NO_CONJUGATE,
    //            m,
    //            &xx(j),
    //            &AA(0,j), i_one,
    //            buff_w,   i_one );

  }  

  return FLA_SUCCESS;
}

