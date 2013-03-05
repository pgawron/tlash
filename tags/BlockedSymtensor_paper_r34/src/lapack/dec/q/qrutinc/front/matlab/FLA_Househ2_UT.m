function [ chi1out, ...
           x2out, tau ] = FLA_Househ2_UT( chi1, ...
                                          x2 )
%
% Compute the Householder transform
% / I - 1/tau / 1  \ ( 1  u2^H ) \
% \           \ u2 /             /
% by computing u2 and tau such that
% / I - 1/tau / 1  \ ( 1  u2^H ) \ / chi1 \  = / +- || x ||_2 \
% \           \ u2 /             / \  x2  /    \       0      /
% 
% where x = / chi1 \
%           \  x2  /
%    tau = ( 1 + u2^H u2 )/2
%
% Upon completion +- || x ||_2 has overwritten chi1 and u2 has overwritten x2
%
% Compute || x ||_2
%
  normx = sqrt( chi1 * chi1 + x2' * x2 );

  chi1out = - sign( chi1 ) * normx;
  x2out = x2 / ( chi1 + sign( chi1 ) * normx );

  tau = ( 1 + x2out' * x2out ) / 2

return

