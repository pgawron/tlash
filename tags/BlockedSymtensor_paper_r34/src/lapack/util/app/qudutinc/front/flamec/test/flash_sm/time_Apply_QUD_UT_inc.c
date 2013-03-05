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
#define FLA_ALG_FRONT     1


void time_Apply_QUD_UT_inc(
                 int n_repeats, int mB, int mC, int mD, int n, int n_rhs, dim_t b_alg,
                 FLA_Obj R_BC, FLA_Obj R_BD, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W,
                 FLA_Obj bR_BC, FLA_Obj bR_BD, FLA_Obj bC, FLA_Obj bD,
                 double *dtime, double *diff, double *gflops );

void time_Apply_QUD_UT_inc(
                 int n_repeats, int mB, int mC, int mD, int n, int n_rhs, dim_t b_alg,
                 FLA_Obj R_BC, FLA_Obj R_BD, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W,
                 FLA_Obj bR_BC, FLA_Obj bR_BD, FLA_Obj bC, FLA_Obj bD,
                 double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    bR_BD_save, bC_save, bD_save;

  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bR_BD, &bR_BD_save );
  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bC, &bC_save );
  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bD, &bD_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLASH_Copy( bR_BD_save, bR_BD );
    FLASH_Copy( bC_save, bC );
    FLASH_Copy( bD_save, bD );

    *dtime = FLA_Clock();

    FLASH_Apply_QUD_UT_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                            T, W,
                               bR_BD,
                            C, bC,
                            D, bD );

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }

  {

    FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                FLA_ONE, R_BD, bR_BD );

    FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                FLA_ONE, R_BC, bR_BC );

    *diff = FLASH_Max_elemwise_diff( bR_BD, bR_BC );
  }

  *gflops = n * n_rhs * ( 2.0 * mC + 2.0 * mD + 0.5 * b_alg + 0.5 ) /
            dtime_old / 1e9;
  if ( FLA_Obj_is_complex( R_BD ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLASH_Obj_free( &bR_BD_save );
  FLASH_Obj_free( &bC_save );
  FLASH_Obj_free( &bD_save );
}

