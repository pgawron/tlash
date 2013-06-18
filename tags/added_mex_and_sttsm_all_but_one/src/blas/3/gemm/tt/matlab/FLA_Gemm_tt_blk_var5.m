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

function [ C_out ] = FLA_Gemm_tt_blk_var5( alpha, A, B, C, nb_alg )

  [ AT, ...
    AB ] = FLA_Part_2x1( A, ...
                         0, 'FLA_TOP' );

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_LEFT' );

  while ( size( AT, 1 ) < size( A, 1 ) )

    b = min( size( AB, 1 ), nb_alg );

    [ A0, ...
      A1, ...
      A2 ] = FLA_Repart_2x1_to_3x1( AT, ...
                                    AB, ...
                                    b, 'FLA_BOTTOM' );

    [ B0, B1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         b, 'FLA_RIGHT' );

    %------------------------------------------------------------%

    C = alpha * A1' * B1' + C;

    %------------------------------------------------------------%

    [ AT, ...
      AB ] = FLA_Cont_with_3x1_to_2x1( A0, ...
                                       A1, ...
                                       A2, ...
                                       'FLA_TOP' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, ...
                                           'FLA_LEFT' );

  end

  C_out = C;


return


