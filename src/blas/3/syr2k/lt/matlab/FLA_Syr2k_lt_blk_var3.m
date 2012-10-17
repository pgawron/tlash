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

function [ C_out ] = FLA_Syr2k_lt_blk_var3( A, B, C, nb_alg )

  [ AL, AR ] = FLA_Part_1x2( A, ...
                               0, 'FLA_LEFT' );

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_LEFT' );

  [ CTL, CTR, ...
    CBL, CBR ] = FLA_Part_2x2( C, ...
                               0, 0, 'FLA_TL' );

  while ( size( AL, 2 ) < size( A, 2 ) )

    b = min( size( AR, 2 ), nb_alg );

    [ A0, A1, A2 ]= FLA_Repart_1x2_to_1x3( AL, AR, ...
                                         b, 'FLA_RIGHT' );

    [ B0, B1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         b, 'FLA_RIGHT' );

    [ C00, C01, C02, ...
      C10, C11, C12, ...
      C20, C21, C22 ] = FLA_Repart_2x2_to_3x3( CTL, CTR, ...
                                               CBL, CBR, ...
                                               b, b, 'FLA_BR' );

    %------------------------------------------------------------%

    C10 = C10 + A1' * B0;
    C21 = C21 + B2' * A1;
    C11 = C11 + A1' * B1 + B1' * A1;

    %------------------------------------------------------------%

    [ AL, AR ] = FLA_Cont_with_1x3_to_1x2( A0, A1, A2, ...
                                           'FLA_LEFT' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, ...
                                           'FLA_LEFT' );

    [ CTL, CTR, ...
      CBL, CBR ] = FLA_Cont_with_3x3_to_2x2( C00, C01, C02, ...
                                             C10, C11, C12, ...
                                             C20, C21, C22, ...
                                             'FLA_TL' );

  end

  C_out = [ CTL, CTR
            CBL, CBR ];

return

