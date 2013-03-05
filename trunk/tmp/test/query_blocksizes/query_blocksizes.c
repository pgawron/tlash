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


int main()
{
    fla_blocksize_t* bp_m;
    fla_blocksize_t* bp_k;
    fla_blocksize_t* bp_n;
    fla_blocksize_t* bp_min;

    FLA_Init();

	bp_m   = FLA_Query_blocksizes( FLA_DIMENSION_M );
	bp_k   = FLA_Query_blocksizes( FLA_DIMENSION_K );
	bp_n   = FLA_Query_blocksizes( FLA_DIMENSION_N );
    bp_min = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	fprintf( stdout, "             m      k      n    min\n" );
	fprintf( stdout, "float    %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_m ), 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_k ), 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_n ), 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_min ) ); 
	fprintf( stdout, "double   %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_m ), 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_k ), 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_n ), 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_min ) ); 
	fprintf( stdout, "complex  %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_m ), 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_k ), 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_n ), 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_min ) ); 
	fprintf( stdout, "dcomplex %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_m ), 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_k ), 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_n ), 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_min ) ); 

    FLA_Blocksize_free( bp_m );
    FLA_Blocksize_free( bp_k );
    FLA_Blocksize_free( bp_n );
    FLA_Blocksize_free( bp_min );

    FLA_Finalize();

	return 0;
}

