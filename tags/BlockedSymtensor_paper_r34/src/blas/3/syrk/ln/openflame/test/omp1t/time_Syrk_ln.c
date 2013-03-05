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
#define FLA_ALG_OPENMP_1TASK     5
#define FLA_ALG_OPENMP_2TASKS    6
#define FLA_ALG_OPENMP_2LOOPS    7




void time_Syrk_ln(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Syrk_ln( 
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old;

  FLA_Obj
    C_old;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );

  for ( irep = 0 ; irep < nrepeats; irep++ ){

    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:
      // Time reference implementation
      REF_Syrk_ln( FLA_ONE, A, FLA_ONE, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var1( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var1( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var1( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 2:{
      // Time variant 2
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var2( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var2( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var2( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    } 
    case 3:{
      // Time variant 3 
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var3( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var3( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var3( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 4:{
      // Time variant 4
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var4( A, C );
        break;
      case FLA_ALG_OPENMP_2TASKS:
        FLA_Syrk_ln_omp2t_var4( A, C );
        break;
      case FLA_ALG_OPENMP_2LOOPS:
        FLA_Syrk_ln_omp2l_var4( A, C );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }
    case 5:{
      // Time variant 5
      switch( type ){
      case FLA_ALG_OPENMP_1TASK:
        FLA_Syrk_ln_omp1t_var5( A, C );
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
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( C, C_ref );
    //FLA_Obj_show( "C:", C, "%f", "\n");
  }

  *gflops = 1.0 * 
            FLA_Obj_length( A ) * 
            FLA_Obj_length( A ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

