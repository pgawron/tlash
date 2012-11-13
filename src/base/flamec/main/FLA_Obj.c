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


#ifdef FLA_ENABLE_SCC
typedef volatile unsigned char* t_vcharp;
t_vcharp RCCE_shmalloc(size_t);
void     RCCE_shfree(t_vcharp);
int      RCCE_ue(void);


void* FLA_shmalloc( size_t size )
{
  return ( void * ) RCCE_shmalloc( size );
}


void FLA_shfree( void* ptr )
{
  RCCE_shfree( ( t_vcharp ) ptr );
}


FLA_Bool FLA_is_owner( void )
{
  if ( RCCE_ue() == 0 )
    return TRUE;
  return FALSE;
}
#endif

FLA_Error FLA_Obj_create_blocked_tensor( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], dim_t blkSize[order], FLA_Obj *obj){

  FLA_Obj_create_blocked_tensor_ext( datatype, FLA_SCALAR, order, size, size, stride, blkSize, obj);

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_tensor( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], FLA_Obj *obj)
{

  FLA_Obj_create_tensor_ext( datatype, FLA_SCALAR, order, size, size, stride, obj );

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create( FLA_Datatype datatype, dim_t m, dim_t n, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  FLA_Obj_create_ext( datatype, FLA_SCALAR, m, n, m, n, rs, cs, obj );

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_blocked_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, dim_t size[order], dim_t size_inner[order], dim_t stride[order], dim_t blkSize[order], FLA_Obj *obj )
{
  dim_t i,j;
  dim_t sizeObj[order];
  dim_t strideObj[order];
  dim_t nBlks = 1;
  strideObj[0] = 1;
  sizeObj[0] = size[0] / blkSize[0];

  for(i = 1; i < order; i++){
	sizeObj[i] = size[i] / blkSize[i];
	strideObj[i] = strideObj[i-1] * sizeObj[i-1];
  }
  //First set up the obj to store the FLA_Objs
  FLA_Obj_create_tensor_ext( datatype, FLA_TENSOR, order, sizeObj, sizeObj, strideObj, obj);

  FLA_Obj* buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

  for(i = 0; i < order; i++)
	nBlks *= sizeObj[i];

  //Create the blocks of obj
  for(i = 0; i < nBlks; i++){
    dim_t strideBlk[order];
    strideBlk[0] = 1;
    for(j = 1; j < order; j++)
		strideBlk[j] = strideBlk[j-1] * blkSize[j-1];
	FLA_Obj_create_tensor( datatype, order, blkSize, strideBlk, &(buf[i]));
  }

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, dim_t size[order], dim_t size_inner[order], dim_t stride[order], FLA_Obj *obj )
{
  dim_t i;
  dim_t nSecondDim = 1;
  for(i = 1; i < order; i++)
	nSecondDim *= size[i];

  //First set the 2-D object info
  FLA_Obj_create_ext( datatype, elemtype, size[0], nSecondDim, size[0], nSecondDim, stride[0], stride[0]*size[0], obj );

  //Update the tensor info
  obj->order = order;
  memcpy(&((obj->size)[0]), &(size[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->size_inner)[0]), &(size_inner[0]), order * sizeof( dim_t ) );
  memset(&((obj->offset)[0]), 0, order * sizeof( dim_t ) );

  obj->base->order = order;
  memcpy(&((obj->base->size)[0]), &(size[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->base->size_inner)[0]), &(size_inner[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof( dim_t ) );
  memset(&((obj->base->index)[0]), 0, order * sizeof( dim_t ) );

  obj->base->n_elem_alloc = size[0] * nSecondDim;

  for(i = 0; i < order; i++)
	obj->permutation[i] = i;
  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t m, dim_t n, dim_t m_inner, dim_t n_inner, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  size_t buffer_size;
  size_t n_elem;

  // Adjust the strides, if necessary.
  FLA_adjust_strides( m, n, &rs, &cs );

/**///Not sure how to extend yet
//  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
//    FLA_Obj_create_ext_check( datatype, elemtype, m, n, m_inner, n_inner, rs, cs, obj );

  // Populate the fields in the view object.
  obj->m                = m;
  obj->n                = n;
  obj->offm             = 0;
  obj->offn             = 0;
  obj->m_inner          = m_inner;
  obj->n_inner          = n_inner;

  // Allocate the base object field.
  obj->base             = ( FLA_Base_obj * ) FLA_malloc( sizeof( FLA_Base_obj ) );

  // Populate the fields in the base object.
  obj->base->datatype   = datatype;
  obj->base->elemtype   = elemtype;
  obj->base->m          = m;
  obj->base->n          = n;
  obj->base->m_inner    = m_inner;
  obj->base->n_inner    = n_inner;
  obj->base->id         = ( unsigned long ) obj->base;
  obj->base->m_index    = 0;
  obj->base->n_index    = 0;

  // Compute the number of elements needed for the buffer, adjusting
  // the strides for alignment if needed.
  n_elem = FLA_compute_num_elem( FLA_Obj_elem_size( *obj ),
                                 m, n, &rs, &cs );

  // Compute the buffer size in bytes.
  buffer_size = ( size_t ) n_elem *
                ( size_t ) FLA_Obj_elem_size( *obj );

  // Allocate the base object's element buffer.
#ifdef FLA_ENABLE_SCC
  obj->base->buffer = ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_malloc( buffer_size ) : FLA_shmalloc( buffer_size ) );
#else
  obj->base->buffer = FLA_malloc( buffer_size );
#endif
  obj->base->buffer_info = 0;

  // Just in case this is a FLASH object, save the number of elements
  // allocated so that we can more easily free the elements later on.
  obj->base->n_elem_alloc = n_elem;

  // Save the row and column strides used in the memory allocation.
  obj->base->rs     = rs;
  obj->base->cs     = cs;

#ifdef FLA_ENABLE_SUPERMATRIX
  // Initialize SuperMatrix fields.
  obj->base->n_read_tasks   = 0;
  obj->base->read_task_head = NULL;
  obj->base->read_task_tail = NULL;
  obj->base->write_task     = NULL;
#endif

  return FLA_SUCCESS;
}


dim_t FLA_compute_num_elem( dim_t elem_size, dim_t m, dim_t n, dim_t* rs, dim_t* cs )
{
  dim_t n_elem;

  // Determine the amount of space we need to allocate based on the values of
  // the row and column strides.
  if ( m == 0 || n == 0 )
  {
    // For empty objects, set the length of the buffer to 0. Row and column
    // strides should remain unchanged (because alignment is not needed).
    n_elem = 0;
  }
  else if ( *rs == 1 )
  {
    // For column-major storage, use cs for computing the length of the buffer
    // to allocate.

    // Align the leading dimension to some user-defined address multiple,
    // if requested at configure-time.
    *cs = FLA_align_ldim( *cs, elem_size );

    // Compute the length of the buffer needed for the object we're creating.
    n_elem = ( size_t ) *cs *
             ( size_t ) n;
  }
  else if ( *cs == 1 )
  {
    // For row-major storage, use rs for computing the length of the buffer
    // to allocate.

    // Align the leading dimension to some user-defined address multiple,
    // if requested at configure-time.
    *rs = FLA_align_ldim( *rs, elem_size );

    // Compute the length of the buffer needed for the object we're creating.
    n_elem = ( size_t ) m *
             ( size_t ) *rs;
  }
  else
  {
    // For general storage, use rs and cs to compute the length of the buffer
    // to allocate.

    // Compute the size of the buffer needed for the object we're creating.
    if ( *rs < *cs )
    {
      *cs = FLA_align_ldim( *cs, elem_size );

      n_elem = ( size_t ) *cs *
               ( size_t ) n;
    }
    else if ( *rs > *cs )
    {
      *rs = FLA_align_ldim( *rs, elem_size );

      n_elem = ( size_t ) m *
               ( size_t ) *rs;
    }
    else // if ( rs == cs )
    {
      //rs = FLA_align_ldim( rs, FLA_Obj_elem_size( *obj ) );
      *cs = FLA_align_ldim( *cs, elem_size );

      // Note that if rs == cs, then we must be creating either a 1-by-n matrix
      // or a m-by-1 matrix. This constraint is enforced in
      // FLA_Check_matrix_strides(). Thus, we can compute the buffer length:
      // m * n * (rs|cs).
      n_elem = ( size_t ) m *
               ( size_t ) n *
               ( size_t ) *cs;
    }
  }

  return n_elem;
}


dim_t FLA_align_ldim( dim_t ldim, dim_t elem_size )
{
#ifdef FLA_ENABLE_MEMORY_ALIGNMENT
  #ifdef FLA_ENABLE_LDIM_ALIGNMENT
    // Increase ldim so that ( ldim * elem_size ) is a multiple of the desired
    // alignment.
    ldim = ( ( ldim * elem_size + FLA_MEMORY_ALIGNMENT_BOUNDARY - 1 ) / 
             FLA_MEMORY_ALIGNMENT_BOUNDARY ) *
           FLA_MEMORY_ALIGNMENT_BOUNDARY /
           elem_size;
  #endif
#endif

  return ldim;
}


void FLA_adjust_strides( dim_t m, dim_t n, dim_t* rs, dim_t* cs )
{
  // Check the strides, and modify them if needed.
  if ( *rs == 0 && *cs == 0 )
  {
    // Default values induce column-major storage, except when m == 1,
    // because we dont want both strides to be unit.
    if ( m == 1 && n > 1 )
    {
      *rs = n;
      *cs = 1;
    }
    else
    {
      *rs = 1;
      *cs = m;
    }
  }
  else if ( *rs == 1 && *cs == 1 )
  {
    // If both strides are unit, this is probably a "lazy" request for a
    // single vector (but could also be a request for a 1xn matrix in column-
    // major order or an mx1 matrix in row-major order). In libflame, we have
    // decided to "reserve" the case where rs == cs == 1 for scalars only, as
    // having unit strides can upset the BLAS error checking when attempting
    // to induce a row-major operation. Also, there is another special case
    // where rs == cs == 1 and one or both of m and n equal zero. This last
    // case is supported to allow creating "empty" objects.

    if ( m == 0 || n == 0 )
    {
      // Nothing needs to be done for the "empty" case where m and/or n
      // equal zero.
    }
    else if ( m == 1 && n == 1 )
    {
      // Nothing needs to be done for the scalar case where m == n == 1.
    }
    else if ( m > 1 && n == 1 )
    {
      // Set the column stride to indicate that this is a column vector stored
      // in column-major order. This is necessary because, in some cases, we
      // have to satisify the error checking in the underlying BLAS library,
      // which expects the leading dimension to be set to at least m, even if
      // it will never be used for indexing since it is a vector and thus only
      // has one column of data.
      *cs = m;
    }
    else if ( m == 1 && n > 1 )
    {
      // Set the row stride to indicate that this is a row vector stored
      // in row-major order.
      *rs = n;
    }
  }
}


FLA_Error FLA_Obj_create_conf_to( FLA_Trans trans, FLA_Obj obj_cur, FLA_Obj *obj_new )
{
  FLA_Datatype datatype;
  FLA_Elemtype elemtype;
  dim_t        m, n;
  dim_t        rs, cs;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_conf_to_check( trans, obj_cur, obj_new );

  datatype = FLA_Obj_datatype( obj_cur );
  elemtype = FLA_Obj_elemtype( obj_cur );

  // Query the dimensions of the existing object.
  if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
  {
    m = FLA_Obj_length( obj_cur );
    n = FLA_Obj_width( obj_cur );
  }
  else // if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
  {
    m = FLA_Obj_width( obj_cur );
    n = FLA_Obj_length( obj_cur );
  }

  // Query the row and column strides of the existing object. We don't care
  // about the actual leading dimension of the existing object, only whether
  // it is in row- or column-major format.
  rs = FLA_Obj_row_stride( obj_cur );
  cs = FLA_Obj_col_stride( obj_cur );

  if ( ( rs == 1 && cs == 1 ) )
  {
    // Do nothing. This special case will be handled by FLA_adjust_strides().
    ;
  }
  else if ( rs == 1 )
  {
    // For column-major storage, use the m dimension as the column stride.
    // Row stride is already set to 1.
    cs = m;
  }
  else if ( cs == 1 )
  {
    // For row-major storage, use the n dimension as the row stride.
    // Column stride is already set to 1.
    rs = n;
  }

  // Handle empty views.
  if ( m == 0 ) cs = 1;
  if ( n == 0 ) rs = 1;

  FLA_Obj_create_ext( datatype, elemtype, m, n, m, n, rs, cs, obj_new );

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_copy_of( FLA_Trans trans, FLA_Obj obj_cur, FLA_Obj *obj_new )
{
  // Create a new object conformal to the current object.
  FLA_Obj_create_conf_to( trans, obj_cur, obj_new );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  // Copy the contents of the current object to the new object.
  FLA_Copyt_external( trans, obj_cur, *obj_new );

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_tensor_without_buffer( FLA_Datatype datatype, dim_t order, dim_t size[order], FLA_Obj *obj ){
	dim_t i;
	dim_t nSecondDim = 1;
	for(i = 1; i < order; i++)
		nSecondDim *= size[i];

	//First set 2-D object info
	FLA_Obj_create_without_buffer(datatype, size[0], nSecondDim, obj);

	//Update tensor info
	obj->order = order;
    memcpy(&((obj->size)[0]), &(size[0]), order * sizeof( dim_t ) );
    memcpy(&((obj->size_inner)[0]), &(size[0]), order * sizeof( dim_t ) );
    memset(&((obj->offset)[0]), 0, order * sizeof( dim_t ) );

    obj->base->order = order;
    memcpy(&((obj->base->size)[0]), &(size[0]), order * sizeof( dim_t ) );
    memcpy(&((obj->base->size_inner)[0]), &(size[0]), order * sizeof( dim_t ) );
    memset(&((obj->base->index)[0]), 0, order * sizeof( dim_t ) );
    memset(&((obj->base->stride)[0]), 0, order * sizeof( dim_t ) );

	obj->base->n_elem_alloc = size[0] * nSecondDim;
	for(i = 0; i < order; i++)
		(obj->permutation)[i] = i;

	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_without_buffer( FLA_Datatype datatype, dim_t m, dim_t n, FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_without_buffer_check( datatype, m, n, obj );

  // Populate the fields in the view object.
  obj->m                = m;
  obj->n                = n;
  obj->offm             = 0;
  obj->offn             = 0;
  obj->m_inner          = m;
  obj->n_inner          = n;

  // Allocate the base object field.
  obj->base             = ( FLA_Base_obj * ) FLA_malloc( sizeof( FLA_Base_obj ) );

  // Populate the fields in the base object.
  obj->base->datatype   = datatype;
  obj->base->elemtype   = FLA_SCALAR;
  obj->base->m          = m;
  obj->base->n          = n;
  obj->base->m_inner    = m;
  obj->base->n_inner    = n;
  obj->base->id         = ( unsigned long ) obj->base;
  obj->base->m_index    = 0;
  obj->base->n_index    = 0;

  // Set the row and column strides to invalid values.
  obj->base->rs         = 0;
  obj->base->cs         = 0;

  // Initialize the base object's element buffer to NULL.
  obj->base->buffer       = NULL;
  obj->base->buffer_info  = 0;
  obj->base->n_elem_alloc = 0;

#ifdef FLA_ENABLE_SUPERMATRIX
  // Initialize SuperMatrix fields.
  obj->base->n_read_tasks   = 0;
  obj->base->read_task_head = NULL;
  obj->base->read_task_tail = NULL;
  obj->base->write_task     = NULL;
#endif

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_constant( double const_real, FLA_Obj *obj )
{
  int*      temp_i;
  float*    temp_s;
  double*   temp_d;
  scomplex* temp_c;
  dcomplex* temp_z;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_constant_check( const_real, obj );

  FLA_Obj_create( FLA_CONSTANT, 1, 1, 0, 0, obj );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  temp_i       = FLA_INT_PTR( *obj );
  temp_s       = FLA_FLOAT_PTR( *obj );
  temp_d       = FLA_DOUBLE_PTR( *obj );
  temp_c       = FLA_COMPLEX_PTR( *obj );
  temp_z       = FLA_DOUBLE_COMPLEX_PTR( *obj );

  *temp_i      = ( int   ) const_real;
  *temp_s      = ( float ) const_real;
  *temp_d      =           const_real;
  temp_c->real = ( float ) const_real;
  temp_c->imag = ( float ) 0.0;
  temp_z->real =           const_real;
  temp_z->imag =           0.0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_constant_ext( float const_s, double const_d, FLA_Obj *obj )
{
  int*      temp_i;
  float*    temp_s;
  double*   temp_d;
  scomplex* temp_c;
  dcomplex* temp_z;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_constant_ext_check( const_s, const_d, obj );

  FLA_Obj_create( FLA_CONSTANT, 1, 1, 0, 0, obj );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  temp_i       = FLA_INT_PTR( *obj );
  temp_s       = FLA_FLOAT_PTR( *obj );
  temp_d       = FLA_DOUBLE_PTR( *obj );
  temp_c       = FLA_COMPLEX_PTR( *obj );
  temp_z       = FLA_DOUBLE_COMPLEX_PTR( *obj );

  *temp_i      = ( int   ) const_s;
  *temp_s      =           const_s;
  *temp_d      =           const_d;
  temp_c->real =           const_s;
  temp_c->imag =           0.0F;
  temp_z->real =           const_d;
  temp_z->imag =           0.0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_complex_constant( double const_real, double const_imag, FLA_Obj *obj )
{
  int*      temp_i;
  float*    temp_s;
  double*   temp_d;
  scomplex* temp_c;
  dcomplex* temp_z;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_complex_constant_check( const_real, const_imag, obj );

  FLA_Obj_create( FLA_CONSTANT, 1, 1, 0, 0, obj );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  temp_i       = FLA_INT_PTR( *obj );
  temp_s       = FLA_FLOAT_PTR( *obj );
  temp_d       = FLA_DOUBLE_PTR( *obj );
  temp_c       = FLA_COMPLEX_PTR( *obj );
  temp_z       = FLA_DOUBLE_COMPLEX_PTR( *obj );

  *temp_i      = ( int   ) const_real;
  *temp_s      = ( float ) const_real;
  *temp_d      =           const_real;
  temp_c->real = ( float ) const_real;
  temp_c->imag = ( float ) const_imag;
  temp_z->real =           const_real;
  temp_z->imag =           const_imag;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_attach_buffer_to_tensor( void *buffer, dim_t order, dim_t stride[order], FLA_Obj *obj ){
	//Omitting some things attach_buffer does because not sure how to extend yet
	obj->base->buffer = buffer;
	memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof(dim_t));

	FLA_Adjust_2D_info(obj);
	
	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_attach_buffer( void *buffer, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  dim_t m, n;

  m = FLA_Obj_length( *obj );
  n = FLA_Obj_width( *obj );

  // Adjust the strides, if necessary.
  FLA_adjust_strides( m, n, &rs, &cs );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_attach_buffer_check( buffer, rs, cs, obj );

  obj->base->buffer      = buffer;
  obj->base->rs          = rs;
  obj->base->cs          = cs;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_buffer( dim_t rs, dim_t cs, FLA_Obj *obj )
{
  size_t buffer_size;
  size_t n_elem;
  dim_t  m, n;

  m = FLA_Obj_length( *obj );
  n = FLA_Obj_width( *obj );

  // Adjust the strides, if necessary.
  FLA_adjust_strides( m, n, &rs, &cs );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_buffer_check( rs, cs, obj );

  // Compute the number of elements needed for the buffer, adjusting
  // the strides for alignment if needed.
  n_elem = FLA_compute_num_elem( FLA_Obj_elem_size( *obj ),
                                 m, n, &rs, &cs );

  // Compute the buffer size in bytes.
  buffer_size = ( size_t ) n_elem *
                ( size_t ) FLA_Obj_elem_size( *obj );

  // Allocate the base object's element buffer.
#ifdef FLA_ENABLE_SCC
  obj->base->buffer = ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_malloc( buffer_size ) : FLA_shmalloc( buffer_size ) );
#else
  obj->base->buffer = FLA_malloc( buffer_size );
#endif
  obj->base->buffer_info = 0;

  // Save the number of elements allocated (for use with FLASH).
  obj->base->n_elem_alloc = n_elem;

  // Save the row and column strides used in the memory allocation.
  obj->base->rs     = rs;
  obj->base->cs     = cs;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_free( FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_free_check( obj );

#ifdef FLA_ENABLE_SCC
  ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_free( obj->base->buffer ) : FLA_shfree( obj->base->buffer ) );
#else
//printf( "freeing buff %p\n", obj->base->buffer ); fflush( stdout );
  FLA_free( obj->base->buffer );
#endif
//printf( "freeing base %p\n", obj->base ); fflush( stdout );
  FLA_free( ( void * ) obj->base );

  obj->offm = 0;
  obj->offn = 0;
  obj->m    = 0;
  obj->n    = 0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_free_without_buffer( FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_free_without_buffer_check( obj );

  FLA_free( ( void * ) obj->base );

  obj->offm = 0;
  obj->offn = 0;
  obj->m    = 0;
  obj->n    = 0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_free_buffer( FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_free_buffer_check( obj );

#ifdef FLA_ENABLE_SCC
  ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_free( obj->base->buffer ) : FLA_shfree( obj->base->buffer ) );
#else
  FLA_free( obj->base->buffer );
#endif
  obj->base->buffer = NULL;

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_blocked_free_buffer( FLA_Obj *obj)
{
	if(FLA_Obj_elemtype(*obj) == FLA_TENSOR || FLA_Obj_elemtype(*obj) == FLA_MATRIX){
		dim_t i;
		FLA_Obj* buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);
		for(i = 0; i < FLA_Obj_num_elem_alloc(*obj); i++){
			FLA_Obj_free_buffer(&(buf[i]));
			FLA_Obj_free_without_buffer(&(buf[i]));
		}
		FLA_Obj_free_buffer(obj);
	}
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_blocked_symm_free_buffer( FLA_Obj *obj)
{
	dim_t order = FLA_Obj_order(*obj);
	dim_t* endIndex = FLA_Obj_size(*obj);
	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	dim_t* stride = FLA_Obj_stride(*obj);
	FLA_Obj* buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

	dim_t update_ptr = order - 1;
	while(TRUE){
		dim_t linIndex;
		FLA_TIndex_to_LinIndex(order, curIndex, stride, &linIndex);

		FLA_Obj_free_buffer(&(buf[linIndex]));
		FLA_Obj_free_without_buffer(&(buf[linIndex]));
		
		//Update
		curIndex[update_ptr]++;
		while(update_ptr < order && curIndex[update_ptr] == endIndex[update_ptr]){
			update_ptr--;
			if(update_ptr < order)
				curIndex[update_ptr]++;
		}
		if(update_ptr >= order)
			break;
		for(dim_t i = update_ptr+1; i < order; i++)
			curIndex[i] = curIndex[update_ptr];
		update_ptr = order - 1;
	}
	FLA_Obj_free_buffer(obj);
	FLA_free(endIndex);
	FLA_free(stride);
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_symm_tensor_without_buffer(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t blkSize, FLA_Obj *obj){
	dim_t i;
	dim_t nTBlks;
	dim_t nBlksPerMode;
	FLA_Obj* t_blks;
	dim_t curIndex[order];
	dim_t endIndex[order];
	dim_t updateIndex;
	dim_t objLinIndex;
	dim_t sizeBlock[order];

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Obj_create_symm_tensor_without_buffer_check( datatype, order, size, blkSize, obj );

	//Figure out how many FLA_Objs we need for intermediate level
	nTBlks = 1;
	nBlksPerMode = size[0] / blkSize;
	for(i = 0; i < order; i++)
		nTBlks *= nBlksPerMode;

	for(i = 0; i < order; i++)
		sizeBlock[i] = blkSize;

	//Create the FLA_Objs
	t_blks = (FLA_Obj*)FLA_malloc(nTBlks * sizeof(FLA_Obj));

	//curIndex is our counter for each block's logical index
	//We loop over this until we hit endIndex  and create the apropriate FLA_Objs
	//updateIndex tells us which index in curIndex we need to update
	//objLinIndex tells us which linear object we are dealing with
	memset(curIndex, 0, order * sizeof(dim_t));
	for(i = 0; i < order; i++)
		endIndex[i] = nBlksPerMode;
	updateIndex = order - 1;
	objLinIndex = 0;

	while(TRUE){
		//Set up the FLA_Obj at this index
		FLA_Obj *curObj = &(t_blks[objLinIndex]);

		FLA_Obj_create_tensor_without_buffer( datatype, order, sizeBlock, curObj);
		//Set the offset array (we will use as an index identifier)
		memcpy(&((curObj->offset)[0]), &(curIndex[0]), order * sizeof(dim_t));

		//Loop update
		//Update current index
		curIndex[updateIndex]++;
		objLinIndex++;
		//If we hit the end, loop until we find the index to update
		while(updateIndex < order && curIndex[updateIndex] == endIndex[updateIndex]){
			updateIndex--;
			if(updateIndex < order)
				curIndex[updateIndex]++;
		}
		//If we run off the edge, we know we are at the end, so break out
		if(updateIndex >= order)
			break;
		//Otherwise, update current index, and reset all others
		for(i = updateIndex + 1; i < order; i++)
			curIndex[i] = 0;
		updateIndex = order - 1;
	}

	//Buffer of tensor blocks created, set the main obj to represent this hierarchy
	
	FLA_Obj_create_tensor_without_buffer( datatype, order, endIndex, obj);
	obj->base->elemtype = FLA_TENSOR;

	dim_t* size_obj = FLA_Obj_size(*obj);
	dim_t stride_obj[order];
	stride_obj[0] = 1;
	for(i = 1; i < order; i++)
		stride_obj[i] = stride_obj[i-1]*size_obj[i-1];
	
	FLA_Obj_attach_buffer_to_tensor(t_blks, order, stride_obj, obj);

	FLA_Adjust_2D_info(obj);

	FLA_free(size_obj);
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_attach_buffer_to_symm_tensor( void *buffer[], dim_t order, dim_t stride[order], FLA_Obj *obj ){
	dim_t i, j, updateIndex;
	dim_t curIndex[order];
	dim_t endIndex[order];
	dim_t* size_obj;
	dim_t countBuffer;
	dim_t objLinIndex;
	dim_t sortedIndex[order];
	dim_t permutation[order];
	dim_t* stride_obj;
	dim_t stride_dup[order];
	FLA_Obj *buffer_obj;

	size_obj = FLA_Obj_size(*obj);
	buffer_obj = (FLA_Obj*)FLA_Obj_base_buffer(*obj);
	stride_obj = FLA_Obj_stride(*obj);

	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	memcpy(&(endIndex[0]), &(size_obj[0]), order * sizeof(dim_t));
	countBuffer = 0;
	objLinIndex = 0;

	FLA_Paired_Sort index_pairs[order];
	//Loop over indices
	//If we hit a unique, set its buffer to the next in the list
	//Otherwise, find the unique index and set the buffer to that buffer (adjusting strides)
	//By looping correctly we will hit the unique before any dupes (I think)
	updateIndex = order - 1;
	while(TRUE){
		objLinIndex = 0;
		for(i = 0; i < order; i++)
			objLinIndex += curIndex[i] * stride_obj[i];

		for(i = 0; i < order; i++){
			index_pairs[i].index = i;
			index_pairs[i].val = curIndex[i];
		}

		qsort(index_pairs, order, sizeof(FLA_Paired_Sort), compare_pairwise_sort);

		for(i = 0; i < order; i++){
			permutation[i] = index_pairs[i].index;
			sortedIndex[i] = index_pairs[i].val;
		}

		//Check if this is unique or not
		dim_t uniqueIndex = TRUE;
		for(i = 0; i < order; i++){
			if(sortedIndex[i] != curIndex[i]){
				uniqueIndex = FALSE;
				break;
			}
		}
		if(uniqueIndex){
			(buffer_obj[objLinIndex].base)->buffer = buffer[countBuffer];
			//Update stride - TODO: MOVE THIS ELSEWHERE
			((buffer_obj[objLinIndex].base)->stride)[0] = 1;
			for(i = 1; i < order; i++)
				((buffer_obj[objLinIndex].base)->stride)[i] = ((buffer_obj[objLinIndex].base)->stride)[i-1] * ((buffer_obj[objLinIndex].base)->size)[i-1];
			countBuffer++;
		}else{
			dim_t ipermutation[order];
			dim_t uniqueLinIndex;
			FLA_TIndex_to_LinIndex(order, stride_obj, sortedIndex, &(uniqueLinIndex));

			//point this non-unique FLA_Obj to the correct base
			//WARNING: HACK
			FLA_free(buffer_obj[objLinIndex].base);
			(buffer_obj[objLinIndex]).base = (buffer_obj[uniqueLinIndex]).base;
			//Set the right permutation
			memcpy(&(ipermutation[0]), &(permutation[0]), order * sizeof(dim_t));
			for(i = 0; i < order; i++)
				ipermutation[permutation[i]] = i;

//			memcpy(&(((buffer_obj[objLinIndex]).permutation)[0]), &(permutation[0]), order * sizeof(dim_t));
			memcpy(&(((buffer_obj[objLinIndex]).permutation)[0]), &(ipermutation[0]), order * sizeof(dim_t));
		}
		FLA_Adjust_2D_info(&(buffer_obj[objLinIndex]));
		//Loop update
		//Update current index
		curIndex[updateIndex]++;
		//If we hit the end, loop until we find the index to update
		while(updateIndex < order && curIndex[updateIndex] == endIndex[updateIndex]){
			updateIndex--;
			if(updateIndex < order)
				curIndex[updateIndex]++;
		}
		//If we run off the edge, we know we are at the end, so break out
		if(updateIndex >= order)
			break;
		//Otherwise, update current index, and reset all others
		for(i = updateIndex+1; i < order; i++)
			curIndex[i] = 0;
		updateIndex = order - 1;
	}

	//Omitting some things attach_buffer does because not sure how to extend yet
	//obj->base->buffer = buffer;
	memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof(dim_t));

	FLA_free(size_obj);
	FLA_free(stride_obj);
	return FLA_SUCCESS;
}
