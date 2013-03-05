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
function [ X ] = Symm_ll_unb_var9( alpha, A, B, C )
%
% Invariant: [ CL, CR ] = [ hatCL + alpha + A * BL, hatCR ]
%
  [ BL, BR ] = FLA_Part_1x2( B, 0, 'FLA_LEFT' );

  [ CL, CR ] = FLA_Part_1x2( C, 0, 'FLA_LEFT' );

  while( size( CL, 2 ) ~= size( C, 2 ) )
    [ C0, c1, C2 ] = FLA_Repart_1x2_to_1x3( CL, CR,
					     1, 'FLA_RIGHT' );

    [ B0, b1, B2 ] = FLA_Repart_1x2_to_1x3( BL, BR,
 				            1, 'FLA_RIGHT' );
%* ********************************************************************** *%
    c1 = alpha * ( tril( A ) + tril( A, -1 )' ) * b1 + c1;
%* ********************************************************************** *%
    [ CL, CR ] = FLA_Cont_with_1x3_to_1x2( C0, c1, C2, 'FLA_LEFT' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, b1, B2, 'FLA_LEFT' );
  end

  X = CL;
  return;
