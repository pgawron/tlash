clear

B = rand( 4, 4 );
C = rand( 4, 4 );
D = rand( 4, 4 );
E = rand( 4, 4 );

T = rand( 2, 4 );
S = rand( 2, 4 );

TB = zeros( 2, 4 );
TD = zeros( 2, 4 );
TE = zeros( 2, 4 );

[ Q, R ] = qr( [ B C
                 D E ] );

[ Bout, Cout, ...
  Dout, Eout, ...
  TBout, ...
  TDout, TEout ] = FLA_QR_2x2( B, C, ...
                               D, E, ...
                               TB, ...
                               TD, TE, 2 );

R

[ Bout Cout
  Dout Eout ]

