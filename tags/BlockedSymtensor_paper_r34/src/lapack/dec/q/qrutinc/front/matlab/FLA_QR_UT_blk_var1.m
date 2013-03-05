function [ A_out, T_out ] = FLA_QR_UT_blk_var1( A, T, nb_alg )
% 
%  UT blocked QR factorization
%
  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );

  [ TL, TR ] = FLA_Part_1x2( T, ...
                             0, 'FLA_LEFT' );

  while ( size( ATL, 1 ) < size( A, 1 ) & ...
          size( ATL, 2 ) < size( A, 2 ) )

    b = min( size( ABR, 1 ), nb_alg );

    [ A00, A01, A02, ...
      A10, A11, A12, ...
      A20, A21, A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                               ABL, ABR, ...
                                               b, b, 'FLA_BR' );

    [ T0, T1, T2 ]= FLA_Repart_1x2_to_1x3( TL, TR, ...
                                           b, 'FLA_RIGHT' );

    %------------------------------------------------------------%

    [ T1T, ...
      T2B ] = FLA_Part_2x1( T1,  b, 'FLA_TOP' );

    [ U1, T1T ] = FLA_QR_UT_Accum_T_unb_var1( [ A11
                                                A21 ], T1T );

    [ A11, ...
      A21 ] = FLA_Part_2x1( U1,  b, 'FLA_TOP' );

    if ( size( A12, 2 ) > 0 )

      [ W12T, ...
        W12B ] = FLA_Part_2x1( T2, b, 'FLA_TOP' );

      U11 = trilu( A11 );
      U21 = A21;
      
      W12T = inv( triu( T1T ) )' * ( U11' * A12 + U21' * A22 );

      A12 = A12 - U11 * W12T;
      A22 = A22 - U21 * W12T;

    end
    
    T1 = [ T1T
           T2B ];

    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00, A01, A02, ...
                                             A10, A11, A12, ...
                                             A20, A21, A22, ...
                                             'FLA_TL' );

    [ TL, TR ] = FLA_Cont_with_1x3_to_1x2( T0, T1, T2, ...
                                           'FLA_LEFT' );

  end

  A_out = [ ATL, ATR
            ABL, ABR ];

  T_out = [ TL, TR ];

return
