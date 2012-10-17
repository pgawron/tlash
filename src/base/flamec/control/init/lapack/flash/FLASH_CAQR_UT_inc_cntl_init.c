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

extern fla_caqr2ut_t*  flash_caqr2ut_cntl;
extern fla_apcaq2ut_t* flash_apcaq2ut_cntl;

fla_caqrutinc_t*     flash_caqrutinc_cntl;
fla_blocksize_t*     flash_caqrutinc_var1_bsize;

void FLASH_CAQR_UT_inc_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_caqrutinc_var1_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_caqrutinc_cntl = FLA_Cntl_caqrutinc_obj_create( FLA_HIER,
	                                                      FLA_BLOCKED_VARIANT1, 
	                                                      flash_caqrutinc_var1_bsize,
	                                                      flash_caqr2ut_cntl,
	                                                      flash_apcaq2ut_cntl );
}

void FLASH_CAQR_UT_inc_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_caqrutinc_cntl );

	FLA_Blocksize_free( flash_caqrutinc_var1_bsize );
}

