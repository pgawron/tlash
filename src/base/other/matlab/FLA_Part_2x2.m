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
  function [ ATL, ATR,...
             ABL, ABR     ] = FLA_Part_2x2( A,...
                                            mb, nb, quadrant )
%
% function [ ATL, ATR,...
%            ABL, ABR     ] = FLA_Part_2x2( A,...
%                                           mb, nb, quadrant )
% Purpose: Partition matrix A into four quadrants
% where the quadrant indicated by quadrant is mb x nb
%
  [ m, n  ] = size( A );
  [ mquadrant, nquadrant ] = size( quadrant );
%
% Check input parameters
%
  if( ( mquadrant ~= 1 )|( nquadrant ~= 6 ) )
    error('quadrant must be a string with contents equal to FLA_TL, FLA_TR, FLA_BL, or FLA_BR');
  elseif( ( ~strcmp( quadrant(1:6), 'FLA_TL' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_TR' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_BL' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_BR' ) ) )
    error('quadrant must be a string with contents equal to FLA_TL, FLA_TR, FLA_BL, or FLA_BR');
  end
%
% Partitioning...
%
  if( strcmp( quadrant(1:6), 'FLA_TL' ) )
    ATL = A(1:mb,1:nb);   
        ATR = A(1:mb,nb+1:n);
    ABL = A(mb+1:m,1:nb); 
        ABR = A(mb+1:m,nb+1:n);
  elseif( strcmp( quadrant(1:6), 'FLA_TR' ) )
    ATL = A(1:mb,1:n-nb);   
        ATR = A(1:mb,n-nb+1:n);
    ABL = A(mb+1:m,1:n-nb); 
        ABR = A(mb+1:m,n-nb+1:n);
  elseif( strcmp( quadrant(1:6), 'FLA_BL' ) )
    ATL = A(1:m-mb,1:nb);   
        ATR = A(1:m-mb,nb+1:n);
    ABL = A(m-mb+1:m,1:nb); 
        ABR = A(m-mb+1:m,nb+1:n);
  else
    ATL = A(1:m-mb,1:n-nb);   
        ATR = A(1:m-mb,n-nb+1:n);
    ABL = A(m-mb+1:m,1:n-nb); 
        ABR = A(m-mb+1:m,n-nb+1:n);
  end
%
  return;
%
% End of FLA_Part_2x2
%
