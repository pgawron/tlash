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

dim_t FLA_Obj_order( FLA_Obj obj )
{
  return obj.order;
}


dim_t FLA_Obj_dimsize( FLA_Obj obj, dim_t mode )
{
  return obj.size[mode];
}


dim_t* FLA_Obj_size( FLA_Obj obj )
{
  dim_t* tmp = FLA_malloc(FLA_MAX_ORDER * sizeof(dim_t));
  memcpy(&(tmp[0]), &(obj.size[0]), FLA_Obj_order( obj ) * sizeof( dim_t ) );
  return tmp;
}


dim_t* FLA_Obj_permutation( FLA_Obj obj )
{
  dim_t* tmp = FLA_malloc(FLA_MAX_ORDER * sizeof(dim_t));
  memcpy(&(tmp[0]), &(obj.permutation[0]), FLA_Obj_order( obj ) * sizeof( dim_t ) );
  return tmp;
}


dim_t* FLA_Obj_stride( FLA_Obj obj )
{
  dim_t* tmp = FLA_malloc(FLA_MAX_ORDER * sizeof(dim_t));
  memcpy(&(tmp[0]), &(((obj.base)->stride)[0]), FLA_Obj_order( obj ) * sizeof( dim_t ) );
  return tmp;
}


dim_t FLA_Obj_dimstride( FLA_Obj obj, dim_t mode )
{
  return ((obj.base)->stride)[mode];
}


dim_t* FLA_Obj_offset( FLA_Obj obj )
{
	dim_t* tmp = FLA_malloc(FLA_MAX_ORDER * sizeof(dim_t));
	memcpy(&(tmp[0]), &(obj.offset[0]), FLA_Obj_order( obj ) * sizeof( dim_t ) );
	return tmp;
}


dim_t FLA_Obj_mode_offset( FLA_Obj obj, dim_t mode )
{
  return obj.offset[mode];
}


dim_t* FLA_Obj_base_size( FLA_Obj obj )
{
  dim_t* tmp = FLA_malloc(FLA_MAX_ORDER * sizeof(dim_t));
  memcpy(&(tmp[0]), &(((obj.base)->size)[0]), FLA_Obj_order( obj ) * sizeof( dim_t ) );
  return tmp;
}


dim_t FLA_Obj_base_dimsize( FLA_Obj obj, dim_t mode )
{
	return ((obj.base)->size)[mode];
}


void* FLA_Obj_tensor_buffer_at_view( FLA_Obj obj )
{
	dim_t i;
	dim_t order;
	char*  buffer;
	size_t elem_size;
	dim_t* offset;
	dim_t* stride;
	size_t byte_offset;
	
//	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
//		FLA_Obj_buffer_at_view_check( obj );
	
	order	    = FLA_Obj_order( obj ); 
	elem_size   = ( size_t ) FLA_Obj_elem_size( obj );
	stride      = ( dim_t* ) FLA_Obj_stride( obj );
	offset      = ( dim_t* ) FLA_Obj_offset( obj );
	
	byte_offset = 0;

	for(i = 0; i < order; i++)
		byte_offset += offset[i] * stride[i];

	byte_offset *= elem_size;
	
	buffer      = ( char * ) (obj.base)->buffer;
	
	FLA_free(stride);
	FLA_free(offset);
	
	return ( void* ) ( buffer + byte_offset );
}


dim_t* FLA_Obj_base_scalar_size(FLA_Obj A){
	FLA_Elemtype elemtype = FLA_Obj_elemtype(A);
	dim_t order = FLA_Obj_order(A);
	dim_t* size = FLA_malloc(order * sizeof(dim_t));
	dim_t i;
	
	if(elemtype == FLA_SCALAR){
		return FLA_Obj_base_size(A);
	}else{
		for(i = 0; i < order; i++)
			size[i] = FLA_Obj_base_scalar_dimsize(A, i);
		return size;
	}
}


dim_t FLA_Obj_base_scalar_dimsize(FLA_Obj A, dim_t mode){
	FLA_Elemtype elemtype = FLA_Obj_elemtype(A);
	dim_t size_base = 0;
	
	if(elemtype == FLA_SCALAR){
		return FLA_Obj_base_dimsize(A, mode);
	}else{
		FLA_Obj* buffer;
		dim_t mode_dimension;
		dim_t stride;
		dim_t i;
		
		buffer = FLA_Obj_base_buffer(A);
		mode_dimension = FLA_Obj_dimsize(A, mode);
		stride = FLA_Obj_dimstride(A, mode);
		
		for(i = 0; i < mode_dimension; i++){
			FLA_Obj obj = buffer[stride * i];
			
			size_base += FLA_Obj_base_dimsize(obj, mode);
		}
	}
	return size_base;
}

dim_t TLA_mode_at_sym_pos( TLA_sym S, dim_t pos ){
	return S.symModes[pos];
}

dim_t TLA_sym_group_of_pos( TLA_sym S, dim_t pos ){
	dim_t i;
	dim_t passedModes = 0;

	for(i = 0; i < S.nSymGroups; i++){
		passedModes += S.symGroupLens[i];
		if(pos < passedModes)
			return i;
	}
	return -1;
}

dim_t TLA_sym_pos_of_mode(TLA_sym S, dim_t mode){
	dim_t i;
	for(i = 0; i < S.order; i++)
	    if(S.symModes[i] == mode)
	        return i;
	return -1;
}

dim_t TLA_sym_group_of_mode( TLA_sym S, dim_t mode){
	return TLA_sym_group_of_pos(S, TLA_sym_pos_of_mode(S, mode));
}

dim_t TLA_sym_group_mode_offset(TLA_sym S, dim_t symGroup){
    dim_t i;
    dim_t offset = 0;
    for(i = 0; i < symGroup; i++)
        offset += S.symGroupLens[i];
    return offset;
}

dim_t TLA_sym_group_of_mode_size(TLA_sym S, dim_t mode){
    dim_t i,j;
    dim_t count = 0;
    for(i = 0; i < S.nSymGroups; i++){
        for(j = 0; j < S.symGroupLens[i]; j++){
            if(S.symModes[count++] == mode){
                return S.symGroupLens[i];
            }
        }
    }
    return -1;
}

dim_t TLA_sym_group_size(TLA_sym S, dim_t symGroup){
	return S.symGroupLens[symGroup];
}
