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

#define FLA_ALG_UNBLOCKED 0
#define FLA_ALG_UNB_OPT   1


void time_QR2_UT(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A_flat, FLA_Obj A_flat_ref, FLA_Obj B_flat, FLA_Obj D_flat, FLA_Obj T_flat,
               double *dtime, double *diff, double *gflops );


void time_QR2_UT(
               int param_combo, int type, int nrepeats, int m, int n,
               FLA_Obj A_flat, FLA_Obj A_flat_ref, FLA_Obj B_flat, FLA_Obj D_flat, FLA_Obj T_flat,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_flat_save;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A_flat, &A_flat_save );

  FLA_Copy( A_flat, A_flat_save );

  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy( A_flat_save, A_flat );

    *dtime = FLA_Clock();

    switch( param_combo ){

    // Time parameter combination 0
    case 0:{
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_QR2_UT_unb_var1( B_flat,
                                       D_flat, T_flat );
        break;
      case FLA_ALG_UNB_OPT:
        FLA_QR2_UT_opt_var1( B_flat,
                                       D_flat, T_flat );
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


  if ( type == FLA_ALG_UNBLOCKED )
  {
    FLA_Copy( A_flat, A_flat_ref );
    //FLA_Copy( T_flat, A_flat_ref );

    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( A_flat, A_flat_ref );
    //*diff = FLA_Max_elemwise_diff( T_flat, A_flat_ref );
  }

  *gflops = 2.0 * n * n * (m - n/3.0) /
            dtime_old / 
            1.0e9;

  if ( FLA_Obj_is_complex( A_flat ) )
    *gflops *= 4.0;

  *dtime = dtime_old;

  FLA_Copy( A_flat_save, A_flat );

  FLA_Obj_free( &A_flat_save );
}

