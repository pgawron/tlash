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
  function [ AT,...
             AB ]   = FLA_Part_2x1( A,...
                                    mb, side )
%
% function [ AT,...
%            AB ]   = FLA_Part_2x1( A,...
%                                   mb, side )
%
% Purpose: Partition matrix A into a top and a bottom side
% where the side indicated by {\tt side} has $ {\tt mb} $ rows
%
  m = size( A, 1 );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ( mside ~= 1 )|( ( nside ~= 7 )&( nside ~= 10 ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( ( nside == 7 )&( ~strcmp( side(1:7), 'FLA_TOP' ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( ( nside == 10 )&( ~strcmp( side(1:10), 'FLA_BOTTOM' ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  end
%
% Partitioning...
%
  if( strcmp( side(1:7), 'FLA_TOP' ) )
    AT = A(1:mb,:); 
    AB = A(mb+1:m,:);
  else
    AT = A(1:m-mb,:); 
    AB = A(m-mb+1:m,:);
  end
%
  return;
%
% End of FLA_Part_2x1
%
