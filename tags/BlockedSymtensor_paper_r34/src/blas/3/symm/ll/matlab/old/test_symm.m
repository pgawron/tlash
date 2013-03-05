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
nvariants = 10
m = 20
n = 10
nb = 4

% create a random matrix A

A = rand( m, m );

% create a random matrix B 

B = rand( m, n );

% create a random matrix C 

C = rand( m, n );

% Create a table to report results

table_ll = zeros( nvariants, 2 );

% Compute D = symm( A ) * B + C with matlab

D = ( tril( A ) + tril( A, -1 )' ) * B + C;

for variant=1:nvariants
% Compute E = symm( A ) * B + C with unblocked variant

  E = Symm_ll( variant, 1.0, A, B, C, 1 );

  table_ll( variant, 1 ) = norm( D - E, 1 );

% Compute E = symm( A ) * B + C with blocked variant

  E = Symm_ll( variant, 1.0, A, B, C, nb );

  table_ll( variant, 2 ) = norm( D - E, 1 );
end

printf(" Results for Symm_ll \n\n");
printf("        |      Difference       \n");
printf("        +--------------------------\n");
printf("variant |  unblocked  |   blocked   \n");
printf("========+=============+=============\n");
for variant=1:nvariants
  printf("  %2d    | %5.2e    |  %5.2e\n", ...
                           variant, ...
                           table_ll( variant, 1 ), ...
                           table_ll( variant, 2 ) );
end


  
