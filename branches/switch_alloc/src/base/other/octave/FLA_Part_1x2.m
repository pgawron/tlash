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
  function [ AL, AR ] = FLA_Part_1x2( A,...
                                      nb, side )
%
% function [ AL, AR ] = FLA_Part_1x2( A,...
%                                     nb, side )
% Purpose: Partition matrix A into a left and a right side
% where the side indicated by side has nb columns
%
  n = size( A, 2 );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ( mside ~= 1 )|( nside < 8 )|( nside > 9 ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( ( nside == 8 )&( ~strcmp( side(1:8), 'FLA_LEFT' ) ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( nside == 9 )
    if ( ~strcmp( side(1:9), 'FLA_RIGHT' ) )
      error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
    end
  end
%
% Partitioning...
%
  if( strcmp( side(1:8), 'FLA_LEFT' ) )
    AL = A(:,1:nb); 
    AR = A(:,nb+1:n);
  else
    AL = A(:,1:n-nb); 
    AR = A(:,n-nb+1:n);
  end
%
  return;
%
% End of FLA_Part_1x2
%
