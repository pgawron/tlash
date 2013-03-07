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



FLA_Error FLA_Obj_copy_view( FLA_Obj A, FLA_Obj* B )
{
  FLA_Obj  A_view;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_copy_view_check( A, B );

  // Set the m_inner and n_inner fields of a temporary copy of A.
  A_view         = A;
  A_view.m_inner = FLASH_Obj_scalar_length( A );
  A_view.n_inner = FLASH_Obj_scalar_width( A ); 

  // Copy the modified view into B. 
  *B = A_view; 

  return FLA_SUCCESS;
}



void FLA_Obj_extract_real_scalar( FLA_Obj alpha, double* alpha_value )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_real_scalar_check( alpha, alpha_value );

  if ( FLA_Obj_is_single_precision( alpha ) )
    *alpha_value = ( double ) *FLA_FLOAT_PTR( alpha );
  else
    *alpha_value = *FLA_DOUBLE_PTR( alpha );
}



void FLA_Obj_extract_complex_scalar( FLA_Obj alpha, dcomplex* alpha_value )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_complex_scalar_check( alpha, alpha_value );

  if ( FLA_Obj_is_single_precision( alpha ) )
  {
    scomplex temp = *FLA_COMPLEX_PTR( alpha );
    alpha_value->real = ( double ) temp.real;
    alpha_value->imag = ( double ) temp.imag;
  }
  else
    *alpha_value = *FLA_DOUBLE_COMPLEX_PTR( alpha );
}



void FLA_Obj_extract_real_part( FLA_Obj alpha, FLA_Obj beta )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_real_part_check( alpha, beta );

  if ( FLA_Obj_is_real( alpha ) )
    FLA_Copy( alpha, beta );
  else // if ( FLA_Obj_is_complex( alpha ) )
  {
    if      ( FLA_Obj_datatype( alpha ) == FLA_COMPLEX )
    {
      scomplex* buff_alpha = FLA_COMPLEX_PTR( alpha );
      float*    buff_beta  = FLA_FLOAT_PTR( beta );

      *buff_beta = buff_alpha->real;
    }
    else if ( FLA_Obj_datatype( alpha ) == FLA_DOUBLE_COMPLEX )
    {
      dcomplex* buff_alpha = FLA_DOUBLE_COMPLEX_PTR( alpha );
      double*   buff_beta  = FLA_DOUBLE_PTR( beta );

      *buff_beta = buff_alpha->real;
    }
  }
}



void FLA_Obj_extract_imag_part( FLA_Obj alpha, FLA_Obj beta )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_imag_part_check( alpha, beta );

  if ( FLA_Obj_is_real( alpha ) )
    FLA_Set( FLA_ZERO, beta );
  else // if ( FLA_Obj_is_complex( alpha ) )
  {
    if      ( FLA_Obj_datatype( alpha ) == FLA_COMPLEX )
    {
      scomplex* buff_alpha = FLA_COMPLEX_PTR( alpha );
      float*    buff_beta  = FLA_FLOAT_PTR( beta );

      *buff_beta = buff_alpha->imag;
    }
    else if ( FLA_Obj_datatype( alpha ) == FLA_DOUBLE_COMPLEX )
    {
      dcomplex* buff_alpha = FLA_DOUBLE_COMPLEX_PTR( alpha );
      double*   buff_beta  = FLA_DOUBLE_PTR( beta );

      *buff_beta = buff_alpha->imag;
    }
  }
}



void FLA_Obj_set_real_part( FLA_Obj alpha, FLA_Obj B )
{
  dim_t m_B;
  dim_t n_B;
  dim_t rs_B;
  dim_t cs_B;
  dim_t i, j;

  m_B  = FLA_Obj_length( B );
  n_B  = FLA_Obj_width( B );
  rs_B = FLA_Obj_row_stride( B );
  cs_B = FLA_Obj_col_stride( B );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_set_real_part_check( alpha, B );

  if ( FLA_Obj_is_complex( B ) )
  {
    if      ( FLA_Obj_datatype( B ) == FLA_COMPLEX )
    {
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      scomplex* buff_B     = FLA_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          scomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->real = *buff_alpha;
        }
      }
    }
    else if ( FLA_Obj_datatype( B ) == FLA_DOUBLE_COMPLEX )
    {
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      dcomplex* buff_B     = FLA_DOUBLE_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          dcomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->real = *buff_alpha;
        }
      }
    }
  }
}



void FLA_Obj_set_imag_part( FLA_Obj alpha, FLA_Obj B )
{
  dim_t m_B;
  dim_t n_B;
  dim_t rs_B;
  dim_t cs_B;
  dim_t i, j;

  m_B  = FLA_Obj_length( B );
  n_B  = FLA_Obj_width( B );
  rs_B = FLA_Obj_row_stride( B );
  cs_B = FLA_Obj_col_stride( B );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_set_imag_part_check( alpha, B );

  if ( FLA_Obj_is_complex( B ) )
  {
    if      ( FLA_Obj_datatype( B ) == FLA_COMPLEX )
    {
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      scomplex* buff_B     = FLA_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          scomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->imag = *buff_alpha;
        }
      }
    }
    else if ( FLA_Obj_datatype( B ) == FLA_DOUBLE_COMPLEX )
    {
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      dcomplex* buff_B     = FLA_DOUBLE_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          dcomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->imag = *buff_alpha;
        }
      }
    }
  }
}



FLA_Error FLA_Obj_fshow( FILE* file, char *s1, FLA_Obj A, char *format, char *s2 )
{
  FLA_Datatype datatype;
  dim_t        i, j, m, n;
  dim_t        rs, cs;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_fshow_check( file, s1, A, format, s2 );

  datatype = FLA_Obj_datatype( A );
  m        = FLA_Obj_length( A );
  n        = FLA_Obj_width( A );
  rs       = FLA_Obj_row_stride( A );
  cs       = FLA_Obj_col_stride( A );

  fprintf( file, "%s\n", s1 );

  switch ( datatype ){

  case FLA_CONSTANT:
  {
    int*      consti = FLA_INT_PTR( A );
    float*    consts = FLA_FLOAT_PTR( A );
    double*   constd = FLA_DOUBLE_PTR( A );
    scomplex* constc = FLA_COMPLEX_PTR( A );
    dcomplex* constz = FLA_DOUBLE_COMPLEX_PTR( A );

    fprintf( file, "int      = %d\n", *consti );
    fprintf( file, "float    = %e\n", *consts );
    fprintf( file, "double   = %e\n", *constd );
    fprintf( file, "scomplex = %e + %e\n", constc->real, constc->imag );
    fprintf( file, "dcomplex = %e + %e\n", constz->real, constc->imag );

    break;
  }

  case FLA_FLOAT:
  {
    float *buffer = ( float * ) FLA_FLOAT_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_DOUBLE:
  {
    double *buffer = ( double * ) FLA_DOUBLE_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buffer = ( scomplex * ) FLA_COMPLEX_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        //fprintf( file, format, buffer[ j*cs + i*rs ].real, buffer[ j*cs + i*rs ].imag );
        //fprintf( file, " " );
		fprintf( file, format, buffer[ j*cs + i*rs ].real );
		fprintf( file, " + " );
		fprintf( file, format, buffer[ j*cs + i*rs ].imag );
		fprintf( file, "  " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buffer = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        //fprintf( file, format, buffer[ j*cs + i*rs ].real, buffer[ j*cs + i*rs ].imag );
        //fprintf( file, " " );
		fprintf( file, format, buffer[ j*cs + i*rs ].real );
		fprintf( file, " + " );
		fprintf( file, format, buffer[ j*cs + i*rs ].imag );
		fprintf( file, "  " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_INT:
  {
    int *buffer = ( int * ) FLA_INT_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  }

  fprintf( file, "%s\n", s2 );

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_print_tensor_under_permutation(FLA_Obj A, dim_t permutation[]){
    dim_t i;
    dim_t order = FLA_Obj_order(A);
	dim_t* size = FLA_Obj_size(A);
	dim_t* stride = FLA_Obj_stride(A);
	dim_t ipermutation[order];

	for(i = 0; i < order; i++)
		ipermutation[permutation[i]] = i;


    if(FLA_Obj_elemtype(A) == FLA_TENSOR || FLA_Obj_elemtype(A) == FLA_MATRIX){
        FLA_Obj* buf = FLA_Obj_base_buffer(A);

		dim_t curIndex[order];
		dim_t updatePtr = 0;
		memset(&(curIndex[0]), 0, order * sizeof(dim_t));
		while(TRUE){
			dim_t linIndex;

			dim_t ipermutedIndex[order];
			FLA_Permute_array(order, curIndex, ipermutation, ipermutedIndex);
			FLA_TIndex_to_LinIndex(order, stride, ipermutedIndex, &linIndex);

			printf("blkIndex: ");
			for(i = 0; i < order; i++)
				printf(" %d", curIndex[i]);
			printf(" ");
			printf("datablk: ");
			FLA_Obj_print_tensor_under_permutation(buf[linIndex], FLA_Obj_permutation(buf[linIndex]));

			//update
			curIndex[updatePtr]++;
			while(updatePtr < order && curIndex[updatePtr] == size[updatePtr]){
				updatePtr++;
				if(updatePtr < order)
					curIndex[updatePtr]++;
			}
			if(updatePtr >= order)
				break;
			for(i = updatePtr - 1; i < order; i--)
				curIndex[i] = 0;
			updatePtr = 0;
		}
		
    }
    else{
        double* buf = (double*)FLA_Obj_base_buffer(A);

		dim_t curIndex[order];
		dim_t updatePtr = 0;
		memset(&(curIndex[0]), 0, order * sizeof(dim_t));
		while(TRUE){
			dim_t linIndex;

			dim_t ipermutedIndex[order];
			FLA_Permute_array(order, curIndex, ipermutation, ipermutedIndex);
			FLA_TIndex_to_LinIndex(order, stride, ipermutedIndex, &linIndex);

			printf("\n\tcurIndex:");
			for(i = 0; i < order; i++)
				printf(" %d", curIndex[i]);
			printf(" %.3f", buf[linIndex]);
			//update
			curIndex[updatePtr]++;
			while(updatePtr < order && curIndex[updatePtr] == size[updatePtr]){
				updatePtr++;
				if(updatePtr < order)
					curIndex[updatePtr]++;
			}
			if(updatePtr >= order)
				break;
			for(i = updatePtr - 1; i < order; i--)
				curIndex[i] = 0;
			updatePtr = 0;
		}
        printf("\n");
    }

	FLA_free(size);
	FLA_free(stride);
    return FLA_SUCCESS;
}

FLA_Error FLA_Obj_print_tensor(FLA_Obj A){
	dim_t* permutation = FLA_Obj_permutation(A);
	FLA_Obj_print_tensor_under_permutation(A, permutation);

	FLA_free(permutation);
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_show( char *s1, FLA_Obj A, char *format, char *s2 )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_show_check( s1, A, format, s2 );

  return FLA_Obj_fshow( stdout, s1, A, format, s2 );
}

FLA_Error FLA_Obj_print_flat_tensor(FLA_Obj A){
	dim_t i;
	FLA_Elemtype elemtype = FLA_Obj_elemtype(A);
	dim_t order = FLA_Obj_order(A);
	dim_t* size = FLA_Obj_size(A);
	dim_t* stride = FLA_Obj_stride(A);

	if(elemtype == FLA_TENSOR || elemtype == FLA_MATRIX){
//		printf("data:");
		dim_t* blkSize = FLA_Obj_size(((FLA_Obj*)FLA_Obj_base_buffer(A))[0]);
		dim_t* blkStride = FLA_Obj_stride(((FLA_Obj*)FLA_Obj_base_buffer(A))[0]);
		dim_t flatSize[order];
		for(i = 0; i < order; i++)
			flatSize[i] = size[i] * blkSize[i];

		dim_t curIndex[order];
		memset(&(curIndex[0]), 0, order * sizeof(dim_t));

		dim_t iblkPermutation[order];
		dim_t ipermutedIndex[order];
		dim_t blkIndex[order];
		dim_t blkLinOffset = 0;
		dim_t innerBlkIndex[order];
		dim_t innerBlkLinOffset = 0;
		dim_t update_ptr = 0;
		while(TRUE){
			for(i = 0; i < order; i++)
				blkIndex[i] = curIndex[i] / blkSize[i];

			FLA_TIndex_to_LinIndex(order, stride, blkIndex, &blkLinOffset);
			FLA_Obj blk = ((FLA_Obj*)FLA_Obj_base_buffer(A))[blkLinOffset];

			memcpy(&(iblkPermutation[0]), &(blk.permutation[0]), order * sizeof(dim_t));
			for(i = 0; i < order; i++)
				iblkPermutation[blk.permutation[i]] = i;

			for(i = 0; i < order; i++)
				innerBlkIndex[i] = curIndex[i] - (blkIndex[i] * blkSize[i]);

			FLA_Permute_array(order, innerBlkIndex, iblkPermutation, ipermutedIndex);
//			FLA_Permute_array(order, innerBlkIndex, blk.permutation, ipermutedIndex);

			FLA_TIndex_to_LinIndex(order, blkStride, ipermutedIndex, &innerBlkLinOffset);

			printf(" %.6f", ((double*)FLA_Obj_base_buffer(blk))[innerBlkLinOffset]);
			
			//update
			curIndex[update_ptr]++;
			while(curIndex[update_ptr] == flatSize[update_ptr] && update_ptr < order){
				update_ptr++;
				if(update_ptr < order)
					curIndex[update_ptr]++;
			}
			if(update_ptr >= order)
				break;

			for(i = update_ptr - 1; i < order; i--)
				curIndex[i] = 0;
			update_ptr = 0;
		}
//		printf("\n");
		FLA_free(blkStride);
		FLA_free(blkSize);
	}else{
//		printf("data:");

		dim_t curIndex[order];
		memset(&(curIndex[0]), 0, order * sizeof(dim_t));

		dim_t ipermutation[order];
		memcpy(&(ipermutation[0]), &(A.permutation[0]), order * sizeof(dim_t));

		for(i = 0; i < order; i++)
			ipermutation[A.permutation[i]] = i;


		dim_t ipermutedSize[order];
		FLA_Permute_array(order, size, ipermutation, ipermutedSize);
		dim_t ipermutedStride[order];
		FLA_Set_tensor_stride(order, ipermutedSize, ipermutedStride);
		dim_t ipermutedIndex[order];
		dim_t linIndex = 0;
		dim_t update_ptr = 0;
		while(TRUE){
			FLA_Permute_array(order, curIndex, ipermutation, ipermutedIndex);
			FLA_TIndex_to_LinIndex(order, ipermutedStride, ipermutedIndex, &linIndex);

			printf(" %.6f", ((double*)FLA_Obj_base_buffer(A))[linIndex]);
			
			//update
			curIndex[update_ptr]++;
			while(curIndex[update_ptr] == size[update_ptr] && update_ptr < order){
				update_ptr++;
				if(update_ptr < order)
					curIndex[update_ptr]++;
			}
			if(update_ptr >= order)
				break;

			for(i = update_ptr - 1; i < order; i--)
				curIndex[i] = 0;
			update_ptr = 0;
		}
//		printf("\n");

	}

	FLA_free(size);
	FLA_free(stride);

	return FLA_SUCCESS;
}
