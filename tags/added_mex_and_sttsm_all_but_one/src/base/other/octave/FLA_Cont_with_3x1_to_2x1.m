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
             AB ] = FLA_Cont_with_3x1_to_2x1( A0,...
                                              A1,...
                                              A2,...
                                              side )
%
% function [ AT,...
%            AB ] = FLA_Cont_with_3x1_to_2x1( A0,...
%                                             A1,...
%                                             A2,...
%                                             side )
%
% Purpose: Update the 2 x 1 partitioning of matrix A by moving the
% boundaries so that A1 is added to the side indicated by side
%
  [ m0, n0 ] = size( A0 );
  [ m1, n1 ] = size( A1 );
  [ m2, n2 ] = size( A2 );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ( n0 ~= n1 )|( n1 ~= n2 ) )
    error('input matrices must have the same number of columns');
  elseif( ( mside ~= 1 )|( ( nside ~= 7 )&( nside ~= 10 ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( ( nside == 7 )&( ~strcmp( side(1:7), 'FLA_TOP' ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( nside == 10 )
    if ( ~strcmp( side(1:10), 'FLA_BOTTOM' ) )
      error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
    end
  end
%
% Continue with... 
%
  if( strcmp( side(1:7), 'FLA_TOP' ) )
    if( ( m0+m1 == 0 )|( n0 == 0 ) )
      AT = zeros( m0+m1, n0 );
    else
      AT = [ A0;...
             A1 ];
    end
    AB = A2;
  else
    AT = A0;
    if( ( m1+m2 == 0 )|( n0 == 0 ) )
      AB = zeros( m1+m2, n0 );
    else
      AB = [ A1;...
             A2 ];
    end
  end
%
  return;
%
% End of FLA_Cont_with_3x1_to_2x1
%
