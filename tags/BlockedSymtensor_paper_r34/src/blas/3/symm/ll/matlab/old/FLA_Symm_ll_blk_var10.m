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
function [ X ] = Symm_ll_blk_var10( alpha, A, B, C, nb )
%
% Invariant: [ CL, CR ] = [ hatCL, hatCR  + alpha + A * BR]
%
  [ BL, BR ] = FLA_Part_1x2( B, 0, 'FLA_RIGHT' );

  [ CL, CR ] = FLA_Part_1x2( C, 0, 'FLA_RIGHT' );

  while( size( CR, 2 ) ~= size( C, 2 ) )
     b = min( size( CL, 2 ), nb );

    [ C0, C1, C2 ] = FLA_Repart_1x2_to_1x3( CL, CR,
					    b, 'FLA_LEFT' );

    [ B0, B1, B2 ] = FLA_Repart_1x2_to_1x3( BL, BR,
 				            b, 'FLA_LEFT' );
%* ********************************************************************** *%
    C1 = Symm_ll_unb_var10( alpha, A, B1, C1 );
%* ********************************************************************** *%
    [ CL, CR ] = FLA_Cont_with_1x3_to_1x2( C0, C1, C2, 'FLA_RIGHT' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, 'FLA_RIGHT' );
  end

  X = CR;
  return;
