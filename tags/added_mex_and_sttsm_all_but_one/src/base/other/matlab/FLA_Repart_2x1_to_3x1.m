%
%  libflame
%  An object-based infrastructure for developing high-performance
%  dense linear algebra libraries.
%
%  Copyright (C) 2011, The University of Texas
%
%  libflame is free software; you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as
%  published by the Free Software Foundation; either version 2.1 of
%  the License, or (at your option) any later version.
%
%  libflame is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%  Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public
%  License along with libflame; if you did not receive a copy, see
%  http://www.gnu.org/licenses/.
%
%  For more information, please contact us at flame@cs.utexas.edu or
%  send mail to:
%
%  Field G. Van Zee and/or
%  Robert A. van de Geijn
%  The University of Texas at Austin
%  Department of Computer Sciences
%  1 University Station C0500
%  Austin TX 78712
%
function [ A0,...
           A1,...
           A2 ] = FLA_Repart_2x1_to_3x1( AT,...
                                         AB,...
                                         mb, side )
%
% Purpose: Repartition a 2 x 1 partitioning of matrix A into
% a 3 x 1 partitioning where submatrix A1 with mb rows is split from
% the side indicated by side
%
  [ mt, nt ] = size( AT );
  [ mbt, nbt ] = size( AB );
  [ mside, nside ] = size( side );

%
% Check input parameters
%
  if( nt ~= nbt )
    error('input matrices must have the same number of columns');
  elseif( ( mside ~= 1 )|( ( nside ~= 7 )&( nside ~= 10 ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( ( nside == 7 )&( ~strcmp( side(1:7), 'FLA_TOP' ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( ( nside == 10 )&( ~strcmp( side(1:10), 'FLA_BOTTOM' ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  end
%
% Repartitioning...
%
  if( strcmp( side(1:7), 'FLA_TOP' ) )
    A0 = AT( 1:mt-mb, : );
    A1 = AT( mt-mb+1:mt, : );
    A2 = AB;
  else
    A0 = AT;
    A1 = AB( 1:mb, : );
    A2 = AB( mb+1:mbt, : );
  end
%
  return;
%
% End of FLA_Repart_2x1_to_3x1
%
