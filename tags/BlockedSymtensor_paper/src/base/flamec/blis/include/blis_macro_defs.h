/*
   libflame
   An object-based infrastructure for developing high-performance
   dense linear algebra libraries.

   Copyright (C) 2011, The University of Texas

   libflame is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   libflame is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with libflame; if you did not receive a copy, see
   http://www.gnu.org/licenses/.

   For more information, please contact us at flame@cs.utexas.edu or
   send mail to:

   Field G. Van Zee and/or
   Robert A. van de Geijn
   The University of Texas at Austin
   Department of Computer Sciences
   1 University Station C0500
   Austin TX 78712
*/

#ifndef BLIS_MACRO_DEFS_H
#define BLIS_MACRO_DEFS_H

// --- Constants ---------------------------------------------------------------

#define BLIS_NO_INTRINSICS  0
#define BLIS_SSE_INTRINSICS 3

// --- boolean ---

#undef FALSE
#define FALSE 0

#undef TRUE
#define TRUE 1

/*
// --- trans ---

#define BLIS_NO_TRANSPOSE      'n'
#define BLIS_TRANSPOSE         't'
#define BLIS_CONJ_NO_TRANSPOSE 'c'
#define BLIS_CONJ_TRANSPOSE    'h'

// --- conj ---

#define BLIS_NO_CONJUGATE      'n'
#define BLIS_CONJUGATE         'c'

// --- uplo ---

#define BLIS_LOWER_TRIANGULAR  'l'
#define BLIS_UPPER_TRIANGULAR  'u'

// --- side ---

#define BLIS_LEFT              'l'
#define BLIS_RIGHT             'r'

// --- diag ---

#define BLIS_NONUNIT_DIAG      'n'
#define BLIS_UNIT_DIAG         'u'
#define BLIS_ZERO_DIAG         'z'
*/

// --- Functional macros -------------------------------------------------------

// --- Type-agnostic ---

// --- min ---

#define bli_min( x, y ) \
( (x) < (y) ? (x) : (y) )

// --- max ---

#define bli_max( x, y ) \
( (x) > (y) ? (x) : (y) )

// --- Type-dependent ---

// --- neg1 ---

// void bli_sneg1( float* x );
#define bli_sneg1( x ) \
*(x)     *= -1.0F;

// void bli_dneg1( double* x );
#define bli_dneg1( x ) \
*(x)     *= -1.0;

// void bli_cneg1( scomplex* x );
#define bli_cneg1( x ) \
(x)->real *= -1.0F; \
(x)->imag *= -1.0F;

// void bli_zneg1( dcomplex* x );
#define bli_zneg1( x ) \
(x)->real *= -1.0; \
(x)->imag *= -1.0;

// --- neg2 ---

// void bli_sneg2( float* x, float* y );
#define bli_sneg2( x, y ) \
*(y)      = -1.0F * *(x);

// void bli_dneg2( double* x, double* y );
#define bli_dneg2( x, y ) \
*(y)      = -1.0  * *(x);

// void bli_cneg2( scomplex* x, scomplex* y );
#define bli_cneg2( x, y ) \
(y)->real = -1.0F * (x)->real; \
(y)->imag = -1.0F * (x)->imag;

// void bli_zneg2( dcomplex* x, dcomplex* y );
#define bli_zneg2( x, y ) \
(y)->real = -1.0  * (x)->real; \
(y)->imag = -1.0  * (x)->imag;

// --- sqrte ---

// void bli_ssqrte( float* alpha, int* error );
#define bli_ssqrte( alpha, error ) \
if ( *(alpha)      < 0.0F ) {                            *(error) = FLA_FAILURE; } \
else { *(alpha)      =  ( float ) sqrt( *(alpha)      ); *(error) = FLA_SUCCESS; }

// void bli_dsqrte( double* alpha, int* error );
#define bli_dsqrte( alpha, error ) \
if ( *(alpha)      < 0.0  ) {                            *(error) = FLA_FAILURE; } \
else { *(alpha)      = ( double ) sqrt( *(alpha)      ); *(error) = FLA_SUCCESS; }

// void bli_csqrte( scomplex* alpha, int* error );
#define bli_csqrte( alpha, error ) \
if ( (alpha)->real < 0.0F ) \
{                     *(error) = FLA_FAILURE; } \
else { \
(alpha)->real =  ( float ) sqrt( (alpha)->real ); \
(alpha)->imag = 0.0F; *(error) = FLA_SUCCESS; }

// void bli_zsqrte( dcomplex* alpha, int* error );
#define bli_zsqrte( alpha, error ) \
if ( (alpha)->real < 0.0  ) \
{                     *(error) = FLA_FAILURE; } \
else { \
(alpha)->real = ( double ) sqrt( (alpha)->real ); \
(alpha)->imag = 0.0;  *(error) = FLA_SUCCESS; }

// --- absval2 ---

// void bli_sabsval2( float* alpha, float* absval );
#define bli_sabsval2( alpha, absval ) \
*(absval) = ( float ) fabs( ( double ) *(alpha) );

// void bli_dabsval2( double* alpha, double* absval );
#define bli_dabsval2( alpha, absval ) \
*(absval) = fabs( *(alpha) );

// void bli_cabsval2( scomplex* alpha, scomplex* absval );
#define bli_cabsval2( alpha, absval ) \
(absval)->real = ( float ) sqrt( ( double ) (alpha)->real * (alpha)->real + \
                                            (alpha)->imag * (alpha)->imag ); \
(absval)->imag = 0.0F;

// void bli_csabsval2( scomplex* alpha, float* absval );
#define bli_csabsval2( alpha, absval ) \
*(absval)      = ( float ) sqrt( ( double ) (alpha)->real * (alpha)->real + \
                                            (alpha)->imag * (alpha)->imag ); \

// void bli_zabsval2( dcomplex* alpha, dcomplex* absval );
#define bli_zabsval2( alpha, absval ) \
(absval)->real = sqrt( (alpha)->real * (alpha)->real + \
                       (alpha)->imag * (alpha)->imag ); \
(absval)->imag = 0.0;

// void bli_zdabsval2( dcomplex* alpha, double* absval );
#define bli_zdabsval2( alpha, absval ) \
*(absval)      = sqrt( (alpha)->real * (alpha)->real + \
                       (alpha)->imag * (alpha)->imag ); \


// --- absqr ---

// void bli_sabsqr( float* alpha );
#define bli_sabsqr( alpha ) \
*(alpha) = *(alpha) * *(alpha);

// void bli_dabsqr( double* alpha );
#define bli_dabsqr( alpha ) \
*(alpha) = *(alpha) * *(alpha);

// void bli_cabsqr( scomplex* alpha );
#define bli_cabsqr( alpha ) \
(alpha)->real = (alpha)->real * (alpha)->real + (alpha)->imag * (alpha)->imag; \
(alpha)->imag = 0.0F;

// void bli_zabsqr( dcomplex* alpha );
#define bli_zabsqr( alpha ) \
(alpha)->real = (alpha)->real * (alpha)->real + (alpha)->imag * (alpha)->imag; \
(alpha)->imag = 0.0;

// --- invscals ---

// void bli_sinvscals( float* a, float* y );
#define bli_sinvscals( a, y ) \
*(y) = *(y) / *(a);

// void bli_dinvscals( double* a, double* y );
#define bli_dinvscals( a, y ) \
*(y) = *(y) / *(a);

// void bli_csinvscals( float* a, scomplex* y );
#define bli_csinvscals( a, y ) \
{ \
(y)->real = (y)->real / *(a); \
(y)->imag = (y)->imag / *(a); \
}

// void bli_cinvscals( scomplex* a, scomplex* y );
#define bli_cinvscals( a, y ) \
{ \
float temp  = (a)->real * (a)->real + (a)->imag * (a)->imag; \
float zreal = ( (y)->real * (a)->real + (y)->imag * (a)->imag ) / temp; \
float zimag = ( (y)->imag * (a)->real - (y)->real * (a)->imag ) / temp; \
(y)->real = zreal; \
(y)->imag = zimag; \
}

// void bli_zdinvscals( double* a, dcomplex* y );
#define bli_zdinvscals( a, y ) \
{ \
(y)->real = (y)->real / *(a); \
(y)->imag = (y)->imag / *(a); \
}

// void bli_zinvscals( dcomplex* a, dcomplex* y );
#define bli_zinvscals( a, y ) \
{ \
double temp  = (a)->real * (a)->real + (a)->imag * (a)->imag; \
double zreal = ( (y)->real * (a)->real + (y)->imag * (a)->imag ) / temp; \
double zimag = ( (y)->imag * (a)->real - (y)->real * (a)->imag ) / temp; \
(y)->real = zreal; \
(y)->imag = zimag; \
}

// --- div3 ---

// void bli_sdiv3( float* x, float* y, float* a );
#define bli_sdiv3( x, y, a ) \
*(a) = *(x) / *(y);

// void bli_ddiv3( double* x, double* y, double* a );
#define bli_ddiv3( x, y, a ) \
*(a) = *(x) / *(y);

// void bli_cdiv3( scomplex* x, scomplex* y, scomplex* a );
#define bli_cdiv3( x, y, a ) \
{ \
float temp  = (y)->real * (y)->real + (y)->imag * (y)->imag; \
float areal = ( (x)->real * (y)->real + (x)->imag * (y)->imag ) / temp; \
float aimag = ( (x)->imag * (y)->real - (x)->real * (y)->imag ) / temp; \
(a)->real = areal; \
(a)->imag = aimag; \
}

// void bli_zdiv3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bli_zdiv3( x, y, a ) \
{ \
double temp  = (y)->real * (y)->real + (y)->imag * (y)->imag; \
double areal = ( (x)->real * (y)->real + (x)->imag * (y)->imag ) / temp; \
double aimag = ( (x)->imag * (y)->real - (x)->real * (y)->imag ) / temp; \
(a)->real = areal; \
(a)->imag = aimag; \
}

// --- add3 ---

// void bli_sadd3( float* x, float* y, float* a );
#define bli_sadd3( x, y, a ) \
*(a) = *(x) + *(y);

// void bli_dadd3( double* x, double* y, double* a );
#define bli_dadd3( x, y, a ) \
*(a) = *(x) + *(y);

// void bli_cadd3( scomplex* x, scomplex* y, scomplex* a );
#define bli_cadd3( x, y, a ) \
{ \
(a)->real = (x)->real + (y)->real; \
(a)->imag = (x)->imag + (y)->imag; \
}

// void bli_zadd3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bli_zadd3( x, y, a ) \
{ \
(a)->real = (x)->real + (y)->real; \
(a)->imag = (x)->imag + (y)->imag; \
}

// --- copys ---

// void bli_scopys( conj_t conj, float* x, float* y );
#define bli_scopys( conj, x, y ) \
*(y) = *(x);

// void bli_dcopys( conj_t conj, double* x, double* y );
#define bli_dcopys( conj, x, y ) \
*(y) = *(x);

// void bli_ccopys( conj_t conj, scomplex* x, scomplex* y );
#define bli_ccopys( conj, x, y ) \
*(y) = *(x); \
if ( bli_is_conj( conj ) ) (y)->imag *= -1.0F;

// void bli_zcopys( conj_t conj, dcomplex* x, dcomplex* y );
#define bli_zcopys( conj, x, y ) \
*(y) = *(x); \
if ( bli_is_conj( conj ) ) (y)->imag *= -1.0;

// --- scals ---

// void bli_sscals( float* a, float* y );
#define bli_sscals( a, y ) \
*(y) = *(a) * *(y);

// void bli_dscals( double* a, double* y );
#define bli_dscals( a, y ) \
*(y) = *(a) * *(y);

// void bli_csscals( float* a, scomplex* y );
#define bli_csscals( a, y ) \
{ \
(y)->real = *(a) * (y)->real; \
(y)->imag = *(a) * (y)->imag; \
}

// void bli_cscals( scomplex* a, scomplex* y );
#define bli_cscals( a, y ) \
{ \
float tempr = (a)->real * (y)->real - (a)->imag * (y)->imag; \
float tempi = (a)->imag * (y)->real + (a)->real * (y)->imag; \
(y)->real = tempr; \
(y)->imag = tempi; \
}

// void bli_zdscals( double* a, dcomplex* y );
#define bli_zdscals( a, y ) \
{ \
(y)->real = *(a) * (y)->real; \
(y)->imag = *(a) * (y)->imag; \
}

// void bli_zscals( dcomplex* a, dcomplex* y );
#define bli_zscals( a, y ) \
{ \
double tempr = (a)->real * (y)->real - (a)->imag * (y)->imag; \
double tempi = (a)->imag * (y)->real + (a)->real * (y)->imag; \
(y)->real = tempr; \
(y)->imag = tempi; \
}

// --- mult3 ---

// void bli_smult3( float* x, float* y, float* a );
#define bli_smult3( x, y, a ) \
*(a) = *(x) * *(y);

// void bli_dmult3( double* x, double* y, double* a );
#define bli_dmult3( x, y, a ) \
*(a) = *(x) * *(y);

// void bli_cmult3( scomplex* x, scomplex* y, scomplex* a );
#define bli_cmult3( x, y, a ) \
{ \
float tempr = (x)->real * (y)->real - (x)->imag * (y)->imag; \
float tempi = (x)->imag * (y)->real + (x)->real * (y)->imag; \
(a)->real = tempr; \
(a)->imag = tempi; \
}

// void bli_zmult3( dcomplex* x, dcomplex* y, dcomplex* a );
#define bli_zmult3( x, y, a ) \
{ \
double tempr = (x)->real * (y)->real - (x)->imag * (y)->imag; \
double tempi = (x)->imag * (y)->real + (x)->real * (y)->imag; \
(a)->real = tempr; \
(a)->imag = tempi; \
}

// --- mult4 ---

// void bli_smult4( float* alpha, float* x, float* y1, float* y2 );
#define bli_smult4( alpha, x, y1, y2 ) \
*(y2) = *(y1) + *(alpha) * *(x);

// void bli_dmult4( double* alpha, double* x, double* y1, double* y2 );
#define bli_dmult4( alpha, x, y1, y2 ) \
*(y2) = *(y1) + *(alpha) * *(x);

// void bli_cmult4( scomplex* alpha, scomplex* x, scomplex* y1, scomplex* y2 );
#define bli_cmult4( alpha, x, y1, y2 ) \
{ \
(y2)->real = (y1)->real + (alpha)->real * (x)->real - (alpha)->imag * (x)->imag; \
(y2)->imag = (y1)->imag + (alpha)->imag * (x)->real + (alpha)->real * (x)->imag; \
}

// void bli_zmult4( dcomplex* alpha, dcomplex* x, dcomplex* y1, dcomplex* y2 );
#define bli_zmult4( alpha, x, y1, y2 ) \
{ \
(y2)->real = (y1)->real + (alpha)->real * (x)->real - (alpha)->imag * (x)->imag; \
(y2)->imag = (y1)->imag + (alpha)->imag * (x)->real + (alpha)->real * (x)->imag; \
}

// --- conjs ---

// void bli_sconjs( float* a );
#define bli_sconjs( a ) \
;

// void bli_dconjs( double* a );
#define bli_dconjs( a ) \
;

// void bli_cconjs( scomplex* a );
#define bli_cconjs( a ) \
(a)->imag *= -1.0F;

// void bli_zconjs( dcomplex* a );
#define bli_zconjs( a ) \
(a)->imag *= -1.0;

// --- copyconj ---

// void bli_scopyconj( float* x, float* y );
#define bli_scopyconj( x, y ) \
*(y) = *(x);

// void bli_dcopyconj( double* x, double* y );
#define bli_dcopyconj( x, y ) \
*(y) = *(x);

// void bli_ccopyconj( scomplex* x, scomplex* y );
#define bli_ccopyconj( x, y ) \
(y)->real =         (x)->real; \
(y)->imag = -1.0F * (x)->imag;

// void bli_zcopyconj( dcomplex* x, dcomplex* y );
#define bli_zcopyconj( x, y ) \
(y)->real =         (x)->real; \
(y)->imag = -1.0  * (x)->imag;

// --- eq1 ---

// void bli_seq1( float* alpha );
#define bli_seq1( alpha ) \
  ( *alpha == 1.0F )

// void bli_deq1( double* alpha );
#define bli_deq1( alpha ) \
  ( *alpha == 1.0 )

// void bli_ceq1( scomplex* alpha );
#define bli_ceq1( alpha ) \
  ( (alpha)->real == 1.0F && (alpha)->imag == 0.0F )

// void bli_zeq1( dcomplex* alpha );
#define bli_zeq1( alpha ) \
  ( (alpha)->real == 1.0 && (alpha)->imag == 0.0 )

// --- Swapping/toggle macros --------------------------------------------------

// --- swap_pointers ---

#define bli_sswap_pointers( a, b ) \
{ \
float* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bli_dswap_pointers( a, b ) \
{ \
double* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bli_cswap_pointers( a, b ) \
{ \
void* temp = (a); \
(a) = (b); \
(b) = temp; \
}

#define bli_zswap_pointers( a, b ) \
{ \
void* temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_ints ---

#define bli_swap_ints( a, b ) \
{ \
int temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_trans ---

#define bli_swap_trans( a, b ) \
{ \
trans_t temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- swap_conj ---

#define bli_swap_conj( a, b ) \
{ \
conj_t temp = (a); \
(a) = (b); \
(b) = temp; \
}

// --- toggle_side ---

#define bli_toggle_side( side ) \
{ \
if ( bli_is_left( side ) ) side = BLIS_RIGHT; \
else                       side = BLIS_LEFT; \
}

// --- toggle_uplo ---

#define bli_toggle_uplo( uplo ) \
{ \
if ( bli_is_lower( uplo ) ) uplo = BLIS_UPPER_TRIANGULAR; \
else                        uplo = BLIS_LOWER_TRIANGULAR; \
}

// --- toggle_trans ---
#define bli_toggle_trans( trans ) \
{ \
if      ( bli_is_notrans( trans ) )     trans = BLIS_TRANSPOSE; \
else if ( bli_is_trans( trans ) )       trans = BLIS_NO_TRANSPOSE; \
else if ( bli_is_conjnotrans( trans ) ) trans = BLIS_CONJ_TRANSPOSE; \
else                                    trans = BLIS_CONJ_NO_TRANSPOSE; \
}

// --- toggle_conjtrans ---
#define bli_toggle_conjtrans( trans ) \
{ \
if      ( bli_is_notrans( trans ) )     trans = BLIS_CONJ_TRANSPOSE; \
else                                    trans = BLIS_NO_TRANSPOSE; \
}

// --- toggle_conj ---

#define bli_toggle_conj( conj ) \
{ \
if ( bli_is_conj( conj ) ) conj = BLIS_NO_CONJUGATE; \
else                       conj = BLIS_CONJUGATE; \
}

#endif // #ifndef BLIS_MACRO_DEFS_H
