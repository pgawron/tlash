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


//
// Level-1 BLAS
//

struct fla_axpy_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_axpy_s* sub_axpy;
};
typedef struct fla_axpy_s fla_axpy_t;


struct fla_axpyt_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_axpyt_s* sub_axpyt;
};
typedef struct fla_axpyt_s fla_axpyt_t;


struct fla_copy_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_copy_s* sub_copy;
};
typedef struct fla_copy_s fla_copy_t;


struct fla_copyt_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_copyt_s* sub_copyt;
};
typedef struct fla_copyt_s fla_copyt_t;


struct fla_copyr_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_copyr_s* sub_copyr;
	struct fla_copy_s*  sub_copy;
};
typedef struct fla_copyr_s fla_copyr_t;


struct fla_scal_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_scal_s* sub_scal;
};
typedef struct fla_scal_s fla_scal_t;


struct fla_scalr_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_scalr_s* sub_scalr;
	struct fla_scal_s*  sub_scal;
};
typedef struct fla_scalr_s fla_scalr_t;


struct fla_swap_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_swap_s* sub_swap;
};
typedef struct fla_swap_s fla_swap_t;


struct fla_tpose_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_tpose_s* sub_trans;
	struct fla_swap_s*  sub_swap;
};
typedef struct fla_tpose_s fla_tpose_t;


#define FLA_Cntl_sub_axpy( cntl )     cntl->sub_axpy
#define FLA_Cntl_sub_axpy1( cntl )    cntl->sub_axpy1
#define FLA_Cntl_sub_axpy2( cntl )    cntl->sub_axpy2
#define FLA_Cntl_sub_axpy3( cntl )    cntl->sub_axpy3
#define FLA_Cntl_sub_axpyt( cntl )    cntl->sub_axpyt
#define FLA_Cntl_sub_copy( cntl )     cntl->sub_copy
#define FLA_Cntl_sub_copyt( cntl )    cntl->sub_copyt
#define FLA_Cntl_sub_copyr( cntl )    cntl->sub_copyr
#define FLA_Cntl_sub_scal( cntl )     cntl->sub_scal
#define FLA_Cntl_sub_scalr( cntl )    cntl->sub_scalr
#define FLA_Cntl_sub_swap( cntl )     cntl->sub_swap
#define FLA_Cntl_sub_trans( cntl )    cntl->sub_trans


fla_axpy_t* FLA_Cntl_axpy_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_axpy_t*      sub_axpy );
fla_axpyt_t* FLA_Cntl_axpyt_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_axpyt_t*     sub_axpyt );
fla_copy_t* FLA_Cntl_copy_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_copy_t*      sub_copy );
fla_copyt_t* FLA_Cntl_copyt_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_copyt_t*     sub_copyt );
fla_copyr_t* FLA_Cntl_copyr_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_copyr_t*     sub_copyr,
                                        fla_copy_t*      sub_copy );
fla_scal_t* FLA_Cntl_scal_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal );
fla_scalr_t* FLA_Cntl_scalr_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_scalr_t*     sub_scalr,
                                        fla_scal_t*      sub_scal );
fla_swap_t* FLA_Cntl_swap_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_swap_t*      sub_swap );
fla_tpose_t* FLA_Cntl_tpose_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_tpose_t*     sub_trans,
                                        fla_swap_t*      sub_swap );

