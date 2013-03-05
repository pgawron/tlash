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

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_BLOCKED   2
#define FLA_ALG_RECURSIVE 3
#define FLA_ALG_OPTIMIZED 4
#define FLA_ALG_OPENMP_1VAR     5
#define FLA_ALG_OPENMP_2VARS    6




void time_Gemm_nn(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
               double *dtime, double *diff, double *gflops );


void time_Gemm_nn( 
	       int variant, int type, int nrepeats, int n, int nb_alg,
	       FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj Cref,
	       double *dtime, double *diff, double *gflops )
{
  int
    irep,
    info, lwork;

  double
    dtime_old,
    d_minus_one = -1.0, d_one = 1.0;

  FLA_Obj
    Cold;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &Cold );

  FLA_Copy_external( C, Cold );

  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( Cold, C );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:
      REF_Gemm_nn( A, B, C );
      break;

    default:
	 printf("trouble\n");
      break;
    }

    if ( irep == 0 )
      dtime_old = FLA_Clock() - *dtime;
    else{
      *dtime = FLA_Clock() - *dtime;
      dtime_old = min( *dtime, dtime_old );
    }
  }


  if ( variant == 0 ){
    //FLA_Copy_external( C, Cref );
    *diff = 0.0;
  }
  else{
    //*diff = FLA_Max_elemwise_diff( C, Cref );
  }

  *gflops = 2.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_width( A ) * 
            FLA_Obj_width( B ) / 
            dtime_old / 
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( Cold, C );

  FLA_Obj_free( &Cold );
}

