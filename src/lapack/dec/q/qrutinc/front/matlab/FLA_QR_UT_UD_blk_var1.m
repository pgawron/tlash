% Copyright 2003 The University of Texas at Austin
%
% For licensing information see
%                http://www.cs.utexas.edu/users/flame/license.html 
%                                                                                 
% Programmed by: Name of author
%                Email of author

function [ U_out, D_out, T_out ] = FLA_QR_UT_UD_blk_var1( U, ...
                                                          D, T, nb_alg )

  [ UTL, UTR, ...
    UBL, UBR ] = FLA_Part_2x2( U, ...
                               0, 0, 'FLA_TL' );

  [ DL, DR ] = FLA_Part_1x2( D, ...
                               0, 'FLA_LEFT' );

  [ TL, TR ] = FLA_Part_1x2( T, ...
                               0, 'FLA_LEFT' );

  while ( size( UTL, 1 ) < size( U, 1 ) )

    b = min( size( UBR, 1 ), nb_alg );

    [ U00, U01, U02, ...
      U10, U11, U12, ...
      U20, U21, U22 ] = FLA_Repart_2x2_to_3x3( UTL, UTR, ...
                                               UBL, UBR, ...
                                               b, b, 'FLA_BR' );

    [ D0, D1, D2 ]= FLA_Repart_1x2_to_1x3( DL, DR, ...
                                         b, 'FLA_RIGHT' );

    [ T0, T1, T2 ]= FLA_Repart_1x2_to_1x3( TL, TR, ...
                                         b, 'FLA_RIGHT' );

    %------------------------------------------------------------%

    [ T1T, ...
      T2B ] = FLA_Part_2x1( T1,  b, 'FLA_TOP' );

    [ U11, ...
      D1, T1T ] = FLA_QR_UT_UD_Accum_T_unb_var1( U11,
                                                 D1, T1T );

    if ( size( U12, 2 ) > 0 )

      [ W2T, ...
        W2B ] = FLA_Part_2x1( T2, b, 'FLA_TOP' );

      W2T = inv( triu( T1T ) )' * ( U12 + D1' * D2 );

      U12 = U12 - W2T;
      D2  = D2  - D1 * W2T;

    end
    
    T1 = [ T1T
           T2B ];

    %------------------------------------------------------------%

    [ UTL, UTR, ...
      UBL, UBR ] = FLA_Cont_with_3x3_to_2x2( U00, U01, U02, ...
                                             U10, U11, U12, ...
                                             U20, U21, U22, ...
                                             'FLA_TL' );

    [ DL, DR ] = FLA_Cont_with_1x3_to_1x2( D0, D1, D2, ...
                                           'FLA_LEFT' );

    [ TL, TR ] = FLA_Cont_with_1x3_to_1x2( T0, T1, T2, ...
                                           'FLA_LEFT' );

  end

  U_out = [ UTL, UTR
            UBL, UBR ];

  D_out = [ DL, DR ];

  T_out = [ TL, TR ];

return
