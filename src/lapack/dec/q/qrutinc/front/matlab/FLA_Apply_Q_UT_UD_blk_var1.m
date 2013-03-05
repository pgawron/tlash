function [ C_out, E_out ] = FLA_Apply_Q_UT_UD_blk_var1( D, T, C, ...
                                                              E )
  nb_alg = size( T, 1 );

  [ DL, DR ] = FLA_Part_1x2( D, ...
                               0, 'FLA_LEFT' );

  [ TL, TR ] = FLA_Part_1x2( T, ...
                               0, 'FLA_LEFT' );

  [ CT, ...
    CB ] = FLA_Part_2x1( C, ...
                         0, 'FLA_TOP' );

  while ( size( DL, 2 ) < size( D, 2 ) )

    b = min( size( DR, 2 ), nb_alg );

    [ D0, D1, D2 ]= FLA_Repart_1x2_to_1x3( DL, DR, ...
                                         b, 'FLA_RIGHT' );

    [ T0, T1, T2 ]= FLA_Repart_1x2_to_1x3( TL, TR, ...
                                         b, 'FLA_RIGHT' );

    [ C0, ...
      C1, ...
      C2 ] = FLA_Repart_2x1_to_3x1( CT, ...
                                    CB, ...
                                    b, 'FLA_BOTTOM' );

    %------------------------------------------------------------%

    [ T1T, ...
      T2B ] = FLA_Part_2x1( T1,  b, 'FLA_TOP' );

    W2T = inv( triu( T1T ) )' * ( C1 + D1' * E );

    C1 = C1 - W2T;
    E  = E  - D1 * W2T;

    %------------------------------------------------------------%

    [ DL, DR ] = FLA_Cont_with_1x3_to_1x2( D0, D1, D2, ...
                                           'FLA_LEFT' );

    [ TL, TR ] = FLA_Cont_with_1x3_to_1x2( T0, T1, T2, ...
                                           'FLA_LEFT' );

    [ CT, ...
      CB ] = FLA_Cont_with_3x1_to_2x1( C0, ...
                                       C1, ...
                                       C2, ...
                                       'FLA_TOP' );

  end

  C_out = [ CT
            CB ];

  E_out = E;


return
