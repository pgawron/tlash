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

#define FLA_ALG_FRONT_OPT0 0
#define FLA_ALG_FRONT_OPT1 1


FLA_Error REF_LU_incpiv( FLA_Obj A, FLA_Obj p );
void time_LU(
              int pivot_combo, int type, int nrepeats, int m, int n, dim_t nb_alg, dim_t nb_flash,
              FLA_Obj A, FLA_Obj p, FLA_Obj x, FLA_Obj b, FLA_Obj norm,
              double *dtime, double *diff, double *gflops );


void time_LU(
              int pivot_combo, int type, int nrepeats, int m, int n, dim_t nb_alg, dim_t nb_flash,
              FLA_Obj A, FLA_Obj p, FLA_Obj x, FLA_Obj b, FLA_Obj norm,
              double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj AH_save, b_save;
  FLA_Obj AH, pH, bH, LH;

  FLASH_LU_incpiv_create_hier_matrices( A, 1, &nb_flash, nb_alg,
                                        &AH, &pH, &LH );
  FLASH_Obj_create_hier_copy_of_flat( b, 1, &nb_flash, &bH );

  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, AH, &AH_save );

  FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, b, &b_save );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( AH_save, AH );

    *dtime = FLA_Clock();

    switch( pivot_combo ){

    case 0:
    {
      switch( type )
      {
      case FLA_ALG_FRONT_OPT0:
        FLASH_LU_incpiv_noopt( AH, pH, LH );
        break;
      case FLA_ALG_FRONT_OPT1:
        FLASH_LU_incpiv_opt1( AH, pH, LH );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );
  }

  {
    FLASH_FS_incpiv( AH, pH, LH, bH );
    FLASH_Trsv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                AH, bH );

    FLASH_Obj_flatten( bH, x );

    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE,
                       A, x, FLA_MINUS_ONE, b );

    FLA_Nrm2_external( b, norm );
    FLA_Obj_extract_real_scalar( norm, diff );
  }

  *gflops = 2.0 / 3.0 * m * m * n /
            dtime_old / 1e9;

  if ( FLA_Obj_is_complex( A ) ) *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy( b_save, b );

  FLASH_Obj_free( &AH );
  FLASH_Obj_free( &pH );
  FLASH_Obj_free( &bH );
  FLASH_Obj_free( &LH );

  FLA_Obj_free( &b_save );
  FLASH_Obj_free( &AH_save );
}

