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
  function [ A0, A1, A2 ] = FLA_Repart_1x2_to_1x3( AL, AR,...
                                                   nb, side )
%
% function [ A0, A1, A2 ] = FLA_Repart_1x2_to_1x3( AL, AR,...
%                                                  nb, side )
%
% Purpose: Repartition a 1 x 2 partitioning of matrix A into
% a 1 x 3 partitioning where submatrix A1 with nb columns is split from
% the side indicated by side
%
  [ ml, nl ] = size( AL );
  [ mr, nr ] = size( AR );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ml ~= mr )
    error('input matrices must have the same number of rows');
  elseif( ( mside ~= 1 )|( nside < 8 )|( nside > 9 ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( ( nside == 8 )&( ~strcmp( side(1:8), 'FLA_LEFT' ) ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( ( nside == 9 )&( ~strcmp( side(1:9), 'FLA_RIGHT' ) ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  end
%
% Repartitioning...
%
  if( strcmp( side(1:8), 'FLA_LEFT' ) )
    A0 = AL( :, 1:nl-nb );
    A1 = AL( :, nl-nb+1:nl );
    A2 = AR;
  else
    A0 = AL;
    A1 = AR( :, 1:nb );
    A2 = AR( :, nb+1:nr );
  end
%
  return;
%
% End of FLA_Repart_1x2_to_1x3
%
