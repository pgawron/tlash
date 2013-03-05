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

#define FLA_ALG_REFERENCE   0
#define FLA_ALG_OPENMP_BVAR 1
#define FLA_ALG_OPENMP_CVAR 2




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
      // Time reference implementation
      REF_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                FLA_ONE, A, B, FLA_ONE, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var1( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 2:{
      // Time variant 2
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var2( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 3:{
      // Time variant 3
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var3( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 4:{
      // Time variant 4
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var4( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 5:{
      // Time variant 5
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var5( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 6:{
      // Time variant 6
      switch( type ){
      case FLA_ALG_OPENMP_BVAR:
        FLA_Gemm_nn_omp_var6( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 13:{
      // Time variant 1->3
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var13( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 15:{
      // Time variant 1->5
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var15( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 31:{
      // Time variant 3->1 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var31( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 35:{
      // Time variant 3->5 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var35( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 51:{
      // Time variant 5->1 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var51( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 53:{
      // Time variant 5->3 
      switch( type ){
      case FLA_ALG_OPENMP_CVAR:
        FLA_Gemm_nn_omp_var53( FLA_ONE, A, B, C, nb_alg );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    }

    if ( irep == 0 )
      dtime_old = FLA_Clock() - *dtime;
    else{
      *dtime = FLA_Clock() - *dtime;
      dtime_old = min( *dtime, dtime_old );
    }
  }


  if ( variant == 0 ){
    FLA_Copy_external( C, Cref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( C, Cref );
    //FLA_Obj_show( "C:", C, "%f", "\n");
  }

  *gflops = 2.0 * 
            FLA_Obj_length( C ) * 
            FLA_Obj_width( C ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( Cold, C );

  FLA_Obj_free( &Cold );
}

