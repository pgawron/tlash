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

int main( int argc, char *argv[] )
{
  int
    m_input, n_input,
    m, n, rs, cs,
    i,
    datatype;

  int blocksize[3];
  int depth;
  double buffer[64];
  double buffer2[64];

  FLA_Obj Af, Ah, Bh;

  FLA_Init();

  fprintf( stdout, "%c Enter hierarchy depth:", '%' );
  scanf( "%d", &depth );
  fprintf( stdout, "%c %d\n", '%', depth );

  for ( i = 0; i < depth; ++i )
  {
    fprintf( stdout, "%c Enter blocksize %d:", '%', i );
    scanf( "%d", &blocksize[i] );
    fprintf( stdout, "%c %d\n", '%', blocksize[i] );
  }

  fprintf( stdout, "%c enter m n: ", '%' );
  scanf( "%d%d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );

  datatype      = FLA_DOUBLE;
  m             = m_input;
  n             = n_input;
  rs            = 1;
  cs            = m_input;

  for( i = 0; i < 64; i++ ) buffer[i] = ( double ) i;
  for( i = 0; i < 64; i++ ) buffer2[i] = ( double ) 0;

  //FLASH_Obj_create( datatype, m, n, depth, blocksize, &Ah );
  FLASH_Obj_create_without_buffer( datatype, m, n, depth, blocksize, &Ah );
  FLASH_Obj_attach_buffer( buffer, rs, cs, &Ah );

  //FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, Ah, &Af );
  //FLASH_Obj_create_hier_conf_to_flat( FLA_NO_TRANSPOSE, Af, depth, blocksize, &Bh );
  //FLASH_Obj_create_flat_copy_of_hier( Ah, &Af );
  //FLASH_Obj_create_hier_copy_of_flat( Af, depth, blocksize, &Bh );
  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, Ah, &Bh );

  //FLASH_Axpy( FLA_TWO, Ah, Bh );
  FLASH_Copy( Ah, Bh );

  //FLA_Obj_create_without_buffer( datatype, 4, 4, &Af );
  //FLA_Obj_attach_buffer( buffer2, 4, &Af );

  //FLASH_Axpy_flat_to_hier( FLA_TWO, Af, 1, 1, Ah );
  //FLASH_Axpy_hier_to_flat( FLA_TWO, 1, 1, Ah, Af );
  //FLASH_Axpy_buffer_to_hier( FLA_ONE, 4, 4, buffer, 4, 1, 1, Ah );

  //FLASH_Axpy_hier_to_buffer( FLA_ONE, 2, 2, Ah, 4, 4, buffer2, 4 );

  //fprintf( stderr, "T: Am An = %d %d\n", FLASH_Obj_scalar_length( Ah ), 
  //                                       FLASH_Obj_scalar_width( Ah ) );

  //FLASH_Random_matrix( Ah );

  //fprintf( stderr, "depth = %d\n", FLASH_Obj_depth( Ah ) );;

/*
  {
    int depth;
    int b_m[4];
    int b_n[4];

    depth = FLASH_Obj_blocksizes( Bh, b_m, b_n );
    fprintf( stderr, "depth = %d\n", depth );;
    fprintf( stderr, "b_m[0] = %d\n", b_m[0] );;
    fprintf( stderr, "b_n[0] = %d\n", b_n[0] );;
  }
*/

  FLASH_Obj_show( "", Ah, "%11.4e", "" );
  FLASH_Obj_show( "", Bh, "%11.4e", "" );

  //FLA_Obj_show( "", Af, "%11.4e", "" );
  //FLASH_print_struct( Ah );

  //fprintf( stderr, "max_diff = %e\n", FLASH_Max_elemwise_diff( Ah, Bh ) );;

  //FLASH_Obj_free_without_buffer( &Ah );
  //FLASH_Obj_free( &Af );
  //FLA_Obj_free( &Af );


  FLA_Finalize();

  return 0;
}

