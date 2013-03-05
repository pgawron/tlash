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

#include "FLA_type_defs.h"

// --- Pointer-accessing FLAME macro definitions ------------------------------------

#define FLA_CONSTANT_I_OFFSET  0
#define FLA_CONSTANT_S_OFFSET  ( sizeof(double) )
#define FLA_CONSTANT_D_OFFSET  ( sizeof(double) + sizeof(double) )
#define FLA_CONSTANT_C_OFFSET  ( sizeof(double) + sizeof(double) + sizeof(double) )
#define FLA_CONSTANT_Z_OFFSET  ( sizeof(double) + sizeof(double) + sizeof(double) + sizeof( scomplex ) )
#define FLA_CONSTANT_SIZE      ( sizeof(double) + sizeof(double) + sizeof(double) + sizeof( scomplex ) + sizeof( dcomplex ) )

#define FLA_INT_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( int * )      ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_I_OFFSET             ) ) : \
                     ( ( ( int * )      ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_FLOAT_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( float * )    ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_S_OFFSET             ) ) : \
                     ( ( ( float * )    ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_DOUBLE_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( double * )   ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_D_OFFSET             ) ) : \
                     ( ( ( double * )   ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_COMPLEX_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( scomplex * ) ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_C_OFFSET             ) ) : \
                     ( ( ( scomplex * ) ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

#define FLA_DOUBLE_COMPLEX_PTR( x ) \
  ( ((x).base)->datatype == FLA_CONSTANT ? \
    ( ( dcomplex * ) ( ( ( char * )     ((x).base)->buffer ) + FLA_CONSTANT_Z_OFFSET             ) ) : \
                     ( ( ( dcomplex * ) ((x).base)->buffer ) + ( size_t ) (x).offn * ((x).base)->cs + \
                                                               ( size_t ) (x).offm * ((x).base)->rs ) )

