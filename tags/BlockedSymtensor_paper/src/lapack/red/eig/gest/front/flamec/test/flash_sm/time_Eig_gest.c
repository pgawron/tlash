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


FLA_Error REF_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

void time_Eig_gest(
                int param_combo, int type, int n_repeats, int n,
                FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj norm,
                double *dtime, double *diff, double *gflops );


void time_Eig_gest(
                int param_combo, int type, int nrepeats, int m,
                FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj norm,
                double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_save = 1.0e9;

  FLA_Obj
    A_save, B_save, A_flat, B_flat;

  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );
  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, B, &B_save );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLASH_Copy( A_save, A );
    FLASH_Copy( B_save, B );

    *dtime = FLA_Clock();

    switch( param_combo ){

    case 0:
    {
      switch( type ){
      case FLA_ALG_FRONT:
        FLASH_Eig_gest( FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    case 1:
    {
      switch( type ){
      case FLA_ALG_FRONT:
        FLASH_Eig_gest( FLA_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    case 2:
    {
      switch( type ){
      case FLA_ALG_FRONT:
        FLASH_Eig_gest( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    case 3:
    {
      switch( type ){
      case FLA_ALG_FRONT:
        FLASH_Eig_gest( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_save = min( *dtime, dtime_save );
  }

  FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );
  FLASH_Obj_create_flat_copy_of_hier( B, &B_flat );

  FLA_Hermitianize( uplo, A_flat );

  // Recover A.
  if ( inv == FLA_NO_INVERSE )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      // A = L' * A_orig * L
      // A_orig = inv(L') * A * inv(L)
      FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A_flat );
      FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
      FLA_Trsm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      // A = U * A_orig * U'
      // A_orig = inv(U) * A * inv(U')
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, A_flat );
      FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
      FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
    }
  }
  else // if ( inv == FLA_INVERSE )
  {
    if ( uplo == FLA_LOWER_TRIANGULAR )
    {
      // A = inv(L) * A_orig * inv(L')
      // A_orig = L * A * L'
      FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A_flat );
      FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
      FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
    }
    else // if ( uplo == FLA_UPPER_TRIANGULAR )
    {
      // A = inv(U') * A_orig * inv(U)
      // A_orig = U' * A * U
      FLA_Hermitianize( FLA_UPPER_TRIANGULAR, A_flat );
      FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
      FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                         FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                         FLA_ONE, B_flat, A_flat );
    }
  }

  FLASH_Obj_hierarchify( A_flat, A );
  *diff = FLASH_Max_elemwise_diff( A, A_save );

  *gflops = 1.0 * m * m * m /
            dtime_save / 1e9;

  if ( FLA_Obj_is_complex( A ) )
    *gflops *= 4.0;

  *dtime = dtime_save;

  FLASH_Copy( A_save, A );
  FLASH_Copy( B_save, B );

  FLA_Obj_free( &A_flat );
  FLA_Obj_free( &B_flat );
  FLASH_Obj_free( &A_save );
  FLASH_Obj_free( &B_save );
}

