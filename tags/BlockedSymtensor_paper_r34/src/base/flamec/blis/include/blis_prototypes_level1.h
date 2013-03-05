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

// --- Level-1 BLAS-like prototypes --------------------------------------------

// --- amax ---

void bli_samax( int n, float*    x, int incx, int* index );
void bli_damax( int n, double*   x, int incx, int* index );
void bli_camax( int n, scomplex* x, int incx, int* index );
void bli_zamax( int n, dcomplex* x, int incx, int* index );

// --- asum ---

void bli_sasum( int n, float*    x, int incx, float*  norm );
void bli_dasum( int n, double*   x, int incx, double* norm );
void bli_casum( int n, scomplex* x, int incx, float*  norm );
void bli_zasum( int n, dcomplex* x, int incx, double* norm );

// --- axpy ---

void bli_saxpy( int n, float*    alpha, float*    x, int incx, float*    y, int incy );
void bli_daxpy( int n, double*   alpha, double*   x, int incx, double*   y, int incy );
void bli_caxpy( int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy );
void bli_zaxpy( int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy );

// --- axpyv ---

void bli_saxpyv( conj_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy );
void bli_daxpyv( conj_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy );
void bli_caxpyv( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy );
void bli_zaxpyv( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy );

// --- axpymt ---

void bli_saxpymt( trans_t trans, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_daxpymt( trans_t trans, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_caxpymt( trans_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zaxpymt( trans_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- axpymrt ---

void bli_saxpymrt( uplo_t uplo, trans_t trans, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_daxpymrt( uplo_t uplo, trans_t trans, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_caxpymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zaxpymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- axpysv ---

void bli_saxpysv( int n, float*    alpha0, float*    alpha1, float*    x, int incx, float*    beta, float*    y, int incy );
void bli_daxpysv( int n, double*   alpha0, double*   alpha1, double*   x, int incx, double*   beta, double*   y, int incy );
void bli_caxpysv( int n, scomplex* alpha0, scomplex* alpha1, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_zaxpysv( int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- axpysmt ---

void bli_saxpysmt( trans_t trans, int m, int n, float*    alpha0, float*    alpha1, float*    a, int a_rs, int a_cs, float*    beta, float*    b, int b_rs, int b_cs );
void bli_daxpysmt( trans_t trans, int m, int n, double*   alpha0, double*   alpha1, double*   a, int a_rs, int a_cs, double*   beta, double*   b, int b_rs, int b_cs );
void bli_caxpysmt( trans_t trans, int m, int n, scomplex* alpha0, scomplex* alpha1, scomplex* a, int a_rs, int a_cs, scomplex* beta, scomplex* b, int b_rs, int b_cs );
void bli_zaxpysmt( trans_t trans, int m, int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* a, int a_rs, int a_cs, dcomplex* beta, dcomplex* b, int b_rs, int b_cs );

// --- conjv ---

void bli_sconjv( int m, float* x, int incx );
void bli_dconjv( int m, double* x, int incx );
void bli_cconjv( int m, scomplex* x, int incx );
void bli_zconjv( int m, dcomplex* x, int incx );

// --- conjm ---

void bli_sconjm( int m, int n, float*    a, int a_rs, int a_cs );
void bli_dconjm( int m, int n, double*   a, int a_rs, int a_cs );
void bli_cconjm( int m, int n, scomplex* a, int a_rs, int a_cs );
void bli_zconjm( int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- conjmr ---

void bli_sconjmr( uplo_t uplo, int m, int n, float*    a, int a_rs, int a_cs );
void bli_dconjmr( uplo_t uplo, int m, int n, double*   a, int a_rs, int a_cs );
void bli_cconjmr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs );
void bli_zconjmr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- copy ---

void bli_scopy( int m, float*    x, int incx, float*    y, int incy );
void bli_dcopy( int m, double*   x, int incx, double*   y, int incy );
void bli_ccopy( int m, scomplex* x, int incx, scomplex* y, int incy );
void bli_zcopy( int m, dcomplex* x, int incx, dcomplex* y, int incy );

// --- copyv ---

void bli_icopyv( conj_t conj, int m, int*      x, int incx, int*      y, int incy );
void bli_scopyv( conj_t conj, int m, float*    x, int incx, float*    y, int incy );
void bli_dcopyv( conj_t conj, int m, double*   x, int incx, double*   y, int incy );
void bli_ccopyv( conj_t conj, int m, scomplex* x, int incx, scomplex* y, int incy );
void bli_zcopyv( conj_t conj, int m, dcomplex* x, int incx, dcomplex* y, int incy );

void bli_sdcopyv( conj_t conj, int m, float*    x, int incx, double*   y, int incy );
void bli_dscopyv( conj_t conj, int m, double*   x, int incx, float*    y, int incy );
void bli_sccopyv( conj_t conj, int m, float*    x, int incx, scomplex* y, int incy );
void bli_cscopyv( conj_t conj, int m, scomplex* x, int incx, float*    y, int incy );
void bli_szcopyv( conj_t conj, int m, float*    x, int incx, dcomplex* y, int incy );
void bli_zscopyv( conj_t conj, int m, dcomplex* x, int incx, float*    y, int incy );
void bli_dccopyv( conj_t conj, int m, double*   x, int incx, scomplex* y, int incy );
void bli_cdcopyv( conj_t conj, int m, scomplex* x, int incx, double*   y, int incy );
void bli_dzcopyv( conj_t conj, int m, double*   x, int incx, dcomplex* y, int incy );
void bli_zdcopyv( conj_t conj, int m, dcomplex* x, int incx, double*   y, int incy );
void bli_czcopyv( conj_t conj, int m, scomplex* x, int incx, dcomplex* y, int incy );
void bli_zccopyv( conj_t conj, int m, dcomplex* x, int incx, scomplex* y, int incy );

// --- copymr ---

void bli_scopymr( uplo_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dcopymr( uplo_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_ccopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zcopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bli_sscopymr( uplo_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_sdcopymr( uplo_t uplo, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_dscopymr( uplo_t uplo, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_sccopymr( uplo_t uplo, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_cscopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_szcopymr( uplo_t uplo, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zscopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_ddcopymr( uplo_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_dccopymr( uplo_t uplo, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_cdcopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_dzcopymr( uplo_t uplo, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zdcopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_cccopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_czcopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zccopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zzcopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- copymrt ---

void bli_scopymrt( uplo_t uplo, trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dcopymrt( uplo_t uplo, trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_ccopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zcopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bli_sscopymrt( uplo_t uplo, trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_sdcopymrt( uplo_t uplo, trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_sccopymrt( uplo_t uplo, trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_szcopymrt( uplo_t uplo, trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_dscopymrt( uplo_t uplo, trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_ddcopymrt( uplo_t uplo, trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_dccopymrt( uplo_t uplo, trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_dzcopymrt( uplo_t uplo, trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_cscopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_cdcopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_cccopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_czcopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zscopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_zdcopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_zccopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zzcopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- copymt ---

void bli_icopymt( trans_t trans, int m, int n, int*      a, int a_rs, int a_cs, int*      b, int b_rs, int b_cs );
void bli_scopymt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dcopymt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_ccopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zcopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bli_sscopymt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_sdcopymt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_dscopymt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_sccopymt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_cscopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_szcopymt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zscopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_ddcopymt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_dccopymt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_cdcopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_dzcopymt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zdcopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_cccopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_czcopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zccopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zzcopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- dot ---

void bli_cdot_in( conj_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho );
void bli_zdot_in( conj_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho );

void bli_sdot( conj_t conj, int n, float*    x, int incx, float*    y, int incy, float*    rho );
void bli_ddot( conj_t conj, int n, double*   x, int incx, double*   y, int incy, double*   rho );
void bli_cdot( conj_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho );
void bli_zdot( conj_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho );

// --- dots ---

void bli_sdots( conj_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    beta, float*    rho );
void bli_ddots( conj_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   beta, double*   rho );
void bli_cdots( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho );
void bli_zdots( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho );

// --- dot2s ---

void bli_sdot2s( conj_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    beta, float*    rho );
void bli_ddot2s( conj_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   beta, double*   rho );
void bli_cdot2s( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho );
void bli_zdot2s( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho );

// --- fnorm ---

void bli_sfnorm( int m, int n, float*    a, int a_rs, int a_cs, float*  norm );
void bli_dfnorm( int m, int n, double*   a, int a_rs, int a_cs, double* norm );
void bli_cfnorm( int m, int n, scomplex* a, int a_rs, int a_cs, float*  norm );
void bli_zfnorm( int m, int n, dcomplex* a, int a_rs, int a_cs, double* norm );

// --- invscalv ---

void bli_sinvscalv(  conj_t conj, int n, float*    alpha, float*    x, int incx );
void bli_dinvscalv(  conj_t conj, int n, double*   alpha, double*   x, int incx );
void bli_csinvscalv( conj_t conj, int n, float*    alpha, scomplex* x, int incx );
void bli_cinvscalv(  conj_t conj, int n, scomplex* alpha, scomplex* x, int incx );
void bli_zdinvscalv( conj_t conj, int n, double*   alpha, dcomplex* x, int incx );
void bli_zinvscalv(  conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx );

// --- invscalm ---

void bli_sinvscalm(  conj_t conj, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bli_dinvscalm(  conj_t conj, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bli_csinvscalm( conj_t conj, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bli_cinvscalm(  conj_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bli_zdinvscalm( conj_t conj, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bli_zinvscalm(  conj_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- nrm2 ---

void bli_snrm2( int n, float*    x, int incx, float*  norm );
void bli_dnrm2( int n, double*   x, int incx, double* norm );
void bli_cnrm2( int n, scomplex* x, int incx, float*  norm );
void bli_znrm2( int n, dcomplex* x, int incx, double* norm );

// --- scal ---

void bli_sscal(  int n, float*    alpha, float*    x, int incx );
void bli_dscal(  int n, double*   alpha, double*   x, int incx );
void bli_csscal( int n, float*    alpha, scomplex* x, int incx );
void bli_cscal(  int n, scomplex* alpha, scomplex* x, int incx );
void bli_zdscal( int n, double*   alpha, dcomplex* x, int incx );
void bli_zscal(  int n, dcomplex* alpha, dcomplex* x, int incx );

// --- scalv ---

void bli_sscalv(  conj_t conj, int n, float*    alpha, float*    x, int incx );
void bli_dscalv(  conj_t conj, int n, double*   alpha, double*   x, int incx );
void bli_csscalv( conj_t conj, int n, float*    alpha, scomplex* x, int incx );
void bli_cscalv(  conj_t conj, int n, scomplex* alpha, scomplex* x, int incx );
void bli_zdscalv( conj_t conj, int n, double*   alpha, dcomplex* x, int incx );
void bli_zscalv(  conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx );

// --- scalm ---

void bli_sscalm(  conj_t conj, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bli_dscalm(  conj_t conj, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bli_csscalm( conj_t conj, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bli_cscalm(  conj_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bli_zdscalm( conj_t conj, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bli_zscalm(  conj_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- scalmr ---

void bli_sscalmr(  uplo_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bli_dscalmr(  uplo_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bli_csscalmr( uplo_t uplo, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bli_cscalmr(  uplo_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bli_zdscalmr( uplo_t uplo, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bli_zscalmr(  uplo_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- swap ---

void bli_sswap( int n, float*    x, int incx, float*    y, int incy );
void bli_dswap( int n, double*   x, int incx, double*   y, int incy );
void bli_cswap( int n, scomplex* x, int incx, scomplex* y, int incy );
void bli_zswap( int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- swapv ---

void bli_sswapv( int n, float*    x, int incx, float*    y, int incy );
void bli_dswapv( int n, double*   x, int incx, double*   y, int incy );
void bli_cswapv( int n, scomplex* x, int incx, scomplex* y, int incy );
void bli_zswapv( int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- swapmt ---

void bli_sswapmt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dswapmt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_cswapmt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zswapmt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

