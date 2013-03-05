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

// --- Name-mangling macro definitions -----------------------------------------

// --- Define Fortran name-mangling macro --------------------------

#define F77_spotrf F77_FUNC( spotrf , SPOTRF )
#define F77_dpotrf F77_FUNC( dpotrf , DPOTRF )
#define F77_cpotrf F77_FUNC( cpotrf , CPOTRF )
#define F77_zpotrf F77_FUNC( zpotrf , ZPOTRF )
      
#define F77_spotf2 F77_FUNC( spotf2 , SPOTF2 )
#define F77_dpotf2 F77_FUNC( dpotf2 , DPOTF2 )
#define F77_cpotf2 F77_FUNC( cpotf2 , CPOTF2 )
#define F77_zpotf2 F77_FUNC( zpotf2 , ZPOTF2 )
      
      
#define F77_sgetrf F77_FUNC( sgetrf , SGETRF )
#define F77_dgetrf F77_FUNC( dgetrf , DGETRF )
#define F77_cgetrf F77_FUNC( cgetrf , CGETRF )
#define F77_zgetrf F77_FUNC( zgetrf , ZGETRF )
      
#define F77_sgetf2 F77_FUNC( sgetf2 , SGETF2 )
#define F77_dgetf2 F77_FUNC( dgetf2 , DGETF2 )
#define F77_cgetf2 F77_FUNC( cgetf2 , CGETF2 )
#define F77_zgetf2 F77_FUNC( zgetf2 , ZGETF2 )
      
      
#define F77_sgeqrf F77_FUNC( sgeqrf , SGEQRF )
#define F77_dgeqrf F77_FUNC( dgeqrf , DGEQRF )
#define F77_cgeqrf F77_FUNC( cgeqrf , CGEQRF )
#define F77_zgeqrf F77_FUNC( zgeqrf , ZGEQRF )
      
#define F77_sgeqr2 F77_FUNC( sgeqr2 , SGEQR2 )
#define F77_dgeqr2 F77_FUNC( dgeqr2 , DGEQR2 )
#define F77_cgeqr2 F77_FUNC( cgeqr2 , CGEQR2 )
#define F77_zgeqr2 F77_FUNC( zgeqr2 , ZGEQR2 )
      
      
#define F77_sgelqf F77_FUNC( sgelqf , SGELQF )
#define F77_dgelqf F77_FUNC( dgelqf , DGELQF )
#define F77_cgelqf F77_FUNC( cgelqf , CGELQF )
#define F77_zgelqf F77_FUNC( zgelqf , ZGELQF )
      
#define F77_sgelq2 F77_FUNC( sgelq2 , SGELQ2 )
#define F77_dgelq2 F77_FUNC( dgelq2 , DGELQ2 )
#define F77_cgelq2 F77_FUNC( cgelq2 , CGELQ2 )
#define F77_zgelq2 F77_FUNC( zgelq2 , ZGELQ2 )
      
      
#define F77_slauum F77_FUNC( slauum , SLAUUM )
#define F77_dlauum F77_FUNC( dlauum , DLAUUM )
#define F77_clauum F77_FUNC( clauum , CLAUUM )
#define F77_zlauum F77_FUNC( zlauum , ZLAUUM )
      
#define F77_slauu2 F77_FUNC( slauu2 , SLAUU2 )
#define F77_dlauu2 F77_FUNC( dlauu2 , DLAUU2 )
#define F77_clauu2 F77_FUNC( clauu2 , CLAUU2 )
#define F77_zlauu2 F77_FUNC( zlauu2 , ZLAUU2 )
      
      
#define F77_strtri F77_FUNC( strtri , STRTRI )
#define F77_dtrtri F77_FUNC( dtrtri , DTRTRI )
#define F77_ctrtri F77_FUNC( ctrtri , CTRTRI )
#define F77_ztrtri F77_FUNC( ztrtri , ZTRTRI )
      
#define F77_strti2 F77_FUNC( strti2 , STRTI2 )
#define F77_dtrti2 F77_FUNC( dtrti2 , DTRTI2 )
#define F77_ctrti2 F77_FUNC( ctrti2 , CTRTI2 )
#define F77_ztrti2 F77_FUNC( ztrti2 , ZTRTI2 )
      
      
#define F77_strsyl F77_FUNC( strsyl , STRSYL )
#define F77_dtrsyl F77_FUNC( dtrsyl , DTRSYL )
#define F77_ctrsyl F77_FUNC( ctrsyl , CTRSYL )
#define F77_ztrsyl F77_FUNC( ztrsyl , ZTRSYL )
      
      
#define F77_sgehrd F77_FUNC( sgehrd , SGEHRD )
#define F77_dgehrd F77_FUNC( dgehrd , DGEHRD )
#define F77_cgehrd F77_FUNC( cgehrd , CGEHRD )
#define F77_zgehrd F77_FUNC( zgehrd , ZGEHRD )
      
#define F77_sgehd2 F77_FUNC( sgehd2 , SGEHD2 )
#define F77_dgehd2 F77_FUNC( dgehd2 , DGEHD2 )
#define F77_cgehd2 F77_FUNC( cgehd2 , CGEHD2 )
#define F77_zgehd2 F77_FUNC( zgehd2 , ZGEHD2 )
      
      
#define F77_ssytrd F77_FUNC( ssytrd , SSYTRD )
#define F77_dsytrd F77_FUNC( dsytrd , DSYTRD )
#define F77_chetrd F77_FUNC( chetrd , CHETRD )
#define F77_zhetrd F77_FUNC( zhetrd , ZHETRD )
      
#define F77_ssytd2 F77_FUNC( ssytd2 , SSYTD2 )
#define F77_dsytd2 F77_FUNC( dsytd2 , DSYTD2 )
#define F77_chetd2 F77_FUNC( chetd2 , CHETD2 )
#define F77_zhetd2 F77_FUNC( zhetd2 , ZHETD2 )
      
      
#define F77_sgebrd F77_FUNC( sgebrd , SGEBRD )
#define F77_dgebrd F77_FUNC( dgebrd , DGEBRD )
#define F77_cgebrd F77_FUNC( cgebrd , CGEBRD )
#define F77_zgebrd F77_FUNC( zgebrd , ZGEBRD )
      
#define F77_sgebd2 F77_FUNC( sgebd2 , SGEBD2 )
#define F77_dgebd2 F77_FUNC( dgebd2 , DGEBD2 )
#define F77_cgebd2 F77_FUNC( cgebd2 , CGEBD2 )
#define F77_zgebd2 F77_FUNC( zgebd2 , ZGEBD2 )
      
      
#define F77_ssygst F77_FUNC( ssygst , SSYGST )
#define F77_dsygst F77_FUNC( dsygst , DSYGST )
#define F77_chegst F77_FUNC( chegst , CHEGST )
#define F77_zhegst F77_FUNC( zhegst , ZHEGST )
      
#define F77_ssygs2 F77_FUNC( ssygs2 , SSYGS2 )
#define F77_dsygs2 F77_FUNC( dsygs2 , DSYGS2 )
#define F77_chegs2 F77_FUNC( chegs2 , CHEGS2 )
#define F77_zhegs2 F77_FUNC( zhegs2 , ZHEGS2 )
      
      
#define F77_slarft F77_FUNC( slarft , SLARFT )
#define F77_dlarft F77_FUNC( dlarft , DLARFT )
#define F77_clarft F77_FUNC( clarft , CLARFT )
#define F77_zlarft F77_FUNC( zlarft , ZLARFT )
      
      
#define F77_slarfg F77_FUNC( slarfg , SLARFG )
#define F77_dlarfg F77_FUNC( dlarfg , DLARFG )
#define F77_clarfg F77_FUNC( clarfg , CLARFG )
#define F77_zlarfg F77_FUNC( zlarfg , ZLARFG )
      
      
#define F77_sorgqr F77_FUNC( sorgqr , SORGQR )
#define F77_dorgqr F77_FUNC( dorgqr , DORGQR )
#define F77_cungqr F77_FUNC( cungqr , CUNGQR )
#define F77_zungqr F77_FUNC( zungqr , ZUNGQR )

      
#define F77_sormqr F77_FUNC( sormqr , SORMQR )
#define F77_dormqr F77_FUNC( dormqr , DORMQR )
#define F77_cunmqr F77_FUNC( cunmqr , CUNMQR )
#define F77_zunmqr F77_FUNC( zunmqr , ZUNMQR )
      
#define F77_sorm2r F77_FUNC( sorm2r , SORM2R )
#define F77_dorm2r F77_FUNC( dorm2r , DORM2R )
#define F77_cunm2r F77_FUNC( cunm2r , CUNM2R )
#define F77_zunm2r F77_FUNC( zunm2r , ZUNM2R )
      
      
#define F77_sormlq F77_FUNC( sormlq , SORMLQ )
#define F77_dormlq F77_FUNC( dormlq , DORMLQ )
#define F77_cunmlq F77_FUNC( cunmlq , CUNMLQ )
#define F77_zunmlq F77_FUNC( zunmlq , ZUNMLQ )
      
#define F77_sorml2 F77_FUNC( sorml2 , SORML2 )
#define F77_dorml2 F77_FUNC( dorml2 , DORML2 )
#define F77_cunml2 F77_FUNC( cunml2 , CUNML2 )
#define F77_zunml2 F77_FUNC( zunml2 , ZUNML2 )
      
      
#define F77_sorgtr F77_FUNC( sorgtr , SORGTR )
#define F77_dorgtr F77_FUNC( dorgtr , DORGTR )
#define F77_cungtr F77_FUNC( cungtr , CUNGTR )
#define F77_zungtr F77_FUNC( zungtr , ZUNGTR )
      
      
#define F77_sormtr F77_FUNC( sormtr , SORMTR )
#define F77_dormtr F77_FUNC( dormtr , DORMTR )
#define F77_cunmtr F77_FUNC( cunmtr , CUNMTR )
#define F77_zunmtr F77_FUNC( zunmtr , ZUNMTR )
      
      
#define F77_sorgbr F77_FUNC( sorgbr , SORGBR )
#define F77_dorgbr F77_FUNC( dorgbr , DORGBR )
#define F77_cungbr F77_FUNC( cungbr , CUNGBR )
#define F77_zungbr F77_FUNC( zungbr , ZUNGBR )
      
      
#define F77_sormbr F77_FUNC( sormbr , SORMBR )
#define F77_dormbr F77_FUNC( dormbr , DORMBR )
#define F77_cunmbr F77_FUNC( cunmbr , CUNMBR )
#define F77_zunmbr F77_FUNC( zunmbr , ZUNMBR )
      
      
#define F77_ssteqr F77_FUNC( ssteqr , SSTEQR )
#define F77_dsteqr F77_FUNC( dsteqr , DSTEQR )
#define F77_csteqr F77_FUNC( csteqr , CSTEQR )
#define F77_zsteqr F77_FUNC( zsteqr , ZSTEQR )
      
      
#define F77_sstedc F77_FUNC( sstedc , SSTEDC )
#define F77_dstedc F77_FUNC( dstedc , DSTEDC )
#define F77_cstedc F77_FUNC( cstedc , CSTEDC )
#define F77_zstedc F77_FUNC( zstedc , ZSTEDC )
      
      
#define F77_sstemr F77_FUNC( sstemr , SSTEMR )
#define F77_dstemr F77_FUNC( dstemr , DSTEMR )
#define F77_cstemr F77_FUNC( cstemr , CSTEMR )
#define F77_zstemr F77_FUNC( zstemr , ZSTEMR )
      
      
#define F77_ssyev  F77_FUNC( ssyev  , SSYEV  )
#define F77_dsyev  F77_FUNC( dsyev  , DSYEV  )
#define F77_cheev  F77_FUNC( cheev  , CHEEV  )
#define F77_zheev  F77_FUNC( zheev  , ZHEEV  )
      
      
#define F77_ssyevd F77_FUNC( ssyevd , SSYEVD )
#define F77_dsyevd F77_FUNC( dsyevd , DSYEVD )
#define F77_cheevd F77_FUNC( cheevd , CHEEVD )
#define F77_zheevd F77_FUNC( zheevd , ZHEEVD )
      
      
#define F77_ssyevr F77_FUNC( ssyevr , SSYEVR )
#define F77_dsyevr F77_FUNC( dsyevr , DSYEVR )
#define F77_cheevr F77_FUNC( cheevr , CHEEVR )
#define F77_zheevr F77_FUNC( zheevr , ZHEEVR )
      
      
#define F77_sbdsqr F77_FUNC( sbdsqr , SBDSQR )
#define F77_dbdsqr F77_FUNC( dbdsqr , DBDSQR )
#define F77_cbdsqr F77_FUNC( cbdsqr , CBDSQR )
#define F77_zbdsqr F77_FUNC( zbdsqr , ZBDSQR )
      
      
#define F77_sbdsdc F77_FUNC( sbdsdc , SBDSDC )
#define F77_dbdsdc F77_FUNC( dbdsdc , DBDSDC )
      
      
#define F77_sgesvd F77_FUNC( sgesvd , SGESVD )
#define F77_dgesvd F77_FUNC( dgesvd , DGESVD )
#define F77_cgesvd F77_FUNC( cgesvd , CGESVD )
#define F77_zgesvd F77_FUNC( zgesvd , ZGESVD )
      
      
#define F77_sgesdd F77_FUNC( sgesdd , SGESDD )
#define F77_dgesdd F77_FUNC( dgesdd , DGESDD )
#define F77_cgesdd F77_FUNC( cgesdd , CGESDD )
#define F77_zgesdd F77_FUNC( zgesdd , ZGESDD )
      
      
#define F77_slaswp F77_FUNC( slaswp , SLASWP )
#define F77_dlaswp F77_FUNC( dlaswp , DLASWP )
#define F77_claswp F77_FUNC( claswp , CLASWP )
#define F77_zlaswp F77_FUNC( zlaswp , ZLASWP )
      
      
#define F77_slaset F77_FUNC( slaset , SLASET )
#define F77_dlaset F77_FUNC( dlaset , DLASET )
#define F77_claset F77_FUNC( claset , CLASET )
#define F77_zlaset F77_FUNC( zlaset , ZLASET )
      

// --- Cholesky factorization ---

void F77_spotrf( char* uplo, int* n, float*    a, int* lda, int* info );
void F77_dpotrf( char* uplo, int* n, double*   a, int* lda, int* info );
void F77_cpotrf( char* uplo, int* n, scomplex* a, int* lda, int* info );
void F77_zpotrf( char* uplo, int* n, dcomplex* a, int* lda, int* info );

void F77_spotf2( char* uplo, int* n, float*    a, int* lda, int* info );
void F77_dpotf2( char* uplo, int* n, double*   a, int* lda, int* info );
void F77_cpotf2( char* uplo, int* n, scomplex* a, int* lda, int* info );
void F77_zpotf2( char* uplo, int* n, dcomplex* a, int* lda, int* info );

// --- LU factorization with partial pivoting ---

void F77_sgetrf( int* m, int* n, float*    a, int* lda, int* ipiv, int* info );
void F77_dgetrf( int* m, int* n, double*   a, int* lda, int* ipiv, int* info );
void F77_cgetrf( int* m, int* n, scomplex* a, int* lda, int* ipiv, int* info );
void F77_zgetrf( int* m, int* n, dcomplex* a, int* lda, int* ipiv, int* info );

void F77_sgetf2( int* m, int* n, float*    a, int* lda, int* ipiv, int* info );
void F77_dgetf2( int* m, int* n, double*   a, int* lda, int* ipiv, int* info );
void F77_cgetf2( int* m, int* n, scomplex* a, int* lda, int* ipiv, int* info );
void F77_zgetf2( int* m, int* n, dcomplex* a, int* lda, int* ipiv, int* info );

// --- QR factorization (classic) ---

void F77_sgeqrf( int* m, int* n, float*    a, int* lda, float*    tau, float*    work, int* lwork, int* info );
void F77_dgeqrf( int* m, int* n, double*   a, int* lda, double*   tau, double*   work, int* lwork, int* info );
void F77_cgeqrf( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info );
void F77_zgeqrf( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info );

void F77_sgeqr2( int* m, int* n, float*    a, int* lda, float*    tau, float*    work, int* info );
void F77_dgeqr2( int* m, int* n, double*   a, int* lda, double*   tau, double*   work, int* info );
void F77_cgeqr2( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info );
void F77_zgeqr2( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info );

// --- LQ factorization (classic) ---

void F77_sgelqf( int* m, int* n, float*    a, int* lda, float*    tau, float*    work, int* lwork, int* info );
void F77_dgelqf( int* m, int* n, double*   a, int* lda, double*   tau, double*   work, int* lwork, int* info );
void F77_cgelqf( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info );
void F77_zgelqf( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info );

void F77_sgelq2( int* m, int* n, float*    a, int* lda, float*    tau, float*    work, int* info );
void F77_dgelq2( int* m, int* n, double*   a, int* lda, double*   tau, double*   work, int* info );
void F77_cgelq2( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info );
void F77_zgelq2( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info );

// --- Triangular-transpose matrix multiply ---

void F77_slauum( char* uplo, int* n, float*    a, int* lda, int* info );
void F77_dlauum( char* uplo, int* n, double*   a, int* lda, int* info );
void F77_clauum( char* uplo, int* n, scomplex* a, int* lda, int* info );
void F77_zlauum( char* uplo, int* n, dcomplex* a, int* lda, int* info );

void F77_slauu2( char* uplo, int* n, float*    a, int* lda, int* info );
void F77_dlauu2( char* uplo, int* n, double*   a, int* lda, int* info );
void F77_clauu2( char* uplo, int* n, scomplex* a, int* lda, int* info );
void F77_zlauu2( char* uplo, int* n, dcomplex* a, int* lda, int* info );

// --- Triangular matrix inversion ---

void F77_strtri( char* uplo, char* diag, int* n, float*    a, int* lda, int* info );
void F77_dtrtri( char* uplo, char* diag, int* n, double*   a, int* lda, int* info );
void F77_ctrtri( char* uplo, char* diag, int* n, scomplex* a, int* lda, int* info );
void F77_ztrtri( char* uplo, char* diag, int* n, dcomplex* a, int* lda, int* info );

void F77_strti2( char* uplo, char* diag, int* n, float*    a, int* lda, int* info );
void F77_dtrti2( char* uplo, char* diag, int* n, double*   a, int* lda, int* info );
void F77_ctrti2( char* uplo, char* diag, int* n, scomplex* a, int* lda, int* info );
void F77_ztrti2( char* uplo, char* diag, int* n, dcomplex* a, int* lda, int* info );

// --- Triangular Sylvester equation solve ---

void F77_strsyl( char* transa, char* transb, int* isgn, int* m, int* n, float*    a, int* lda, float*    b, int* ldb, float*    c, int* ldc, float*    scale, int* info );
void F77_dtrsyl( char* transa, char* transb, int* isgn, int* m, int* n, double*   a, int* lda, double*   b, int* ldb, double*   c, int* ldc, double*   scale, int* info );
void F77_ctrsyl( char* transa, char* transb, int* isgn, int* m, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* c, int* ldc, float*    scale, int* info );
void F77_ztrsyl( char* transa, char* transb, int* isgn, int* m, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* c, int* ldc, double*   scale, int* info );

// --- Reduction to upper Hessenberg form ---

void F77_sgehrd( int* n, int* ilo, int* ihi, float*    a, int* lda, float*    tau, float*    work, int* lwork, int* info );
void F77_dgehrd( int* n, int* ilo, int* ihi, double*   a, int* lda, double*   tau, double*   work, int* lwork, int* info );
void F77_cgehrd( int* n, int* ilo, int* ihi, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info );
void F77_zgehrd( int* n, int* ilo, int* ihi, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info );

void F77_sgehd2( int* n, int* ilo, int* ihi, float*    a, int* lda, float*    tau, float*    work, int* info );
void F77_dgehd2( int* n, int* ilo, int* ihi, double*   a, int* lda, double*   tau, double*   work, int* info );
void F77_cgehd2( int* n, int* ilo, int* ihi, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info );
void F77_zgehd2( int* n, int* ilo, int* ihi, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info );

// --- Reduction to tridiagonal form ---

void F77_ssytrd( char* uplo, int* n, float*    a, int* lda, float*  d, float*  e, float*    tau, float*    work, int* lwork, int* info );
void F77_dsytrd( char* uplo, int* n, double*   a, int* lda, double* d, double* e, double*   tau, double*   work, int* lwork, int* info );
void F77_chetrd( char* uplo, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tau, scomplex* work, int* lwork, int* info );
void F77_zhetrd( char* uplo, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tau, dcomplex* work, int* lwork, int* info );

void F77_ssytd2( char* uplo, int* n, float*    a, int* lda, float*  d, float*  e, float*    tau, int* info );
void F77_dsytd2( char* uplo, int* n, double*   a, int* lda, double* d, double* e, double*   tau, int* info );
void F77_chetd2( char* uplo, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tau, int* info );
void F77_zhetd2( char* uplo, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tau, int* info );

// --- Reduction to bidiagonal form ---

void F77_sgebrd( int* m, int* n, float*    a, int* lda, float*  d, float*  e, float*    tauq, float*    taup, float*    work, int* lwork, int* info );
void F77_dgebrd( int* m, int* n, double*   a, int* lda, double* d, double* e, double*   tauq, double*   taup, double*   work, int* lwork, int* info );
void F77_cgebrd( int* m, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, int* lwork, int* info );
void F77_zgebrd( int* m, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, int* lwork, int* info );

void F77_sgebd2( int* m, int* n, float*    a, int* lda, float*  d, float*  e, float*    tauq, float*    taup, float*    work, int* info );
void F77_dgebd2( int* m, int* n, double*   a, int* lda, double* d, double* e, double*   tauq, double*   taup, double*   work, int* info );
void F77_cgebd2( int* m, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, int* info );
void F77_zgebd2( int* m, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, int* info );

// --- Reduce Hermitian-definite generalized eigenproblem to standard form ---

void F77_ssygst( int* itype, char* uplo, int* n, float*    a, int* lda, float*    b, int* ldb, int* info );
void F77_dsygst( int* itype, char* uplo, int* n, double*   a, int* lda, double*   b, int* ldb, int* info );
void F77_chegst( int* itype, char* uplo, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, int* info );
void F77_zhegst( int* itype, char* uplo, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, int* info );

void F77_ssygs2( int* itype, char* uplo, int* n, float*    a, int* lda, float*    b, int* ldb, int* info );
void F77_dsygs2( int* itype, char* uplo, int* n, double*   a, int* lda, double*   b, int* ldb, int* info );
void F77_chegs2( int* itype, char* uplo, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, int* info );
void F77_zhegs2( int* itype, char* uplo, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, int* info );

// --- Accumulate block Householder matrix T (classic) ---

void F77_slarft( char* direct, char* storev, int* n, int* k, float*    v, int* ldv, float*    tau, float*    t, int* ldt );
void F77_dlarft( char* direct, char* storev, int* n, int* k, double*   v, int* ldv, double*   tau, double*   t, int* ldt );
void F77_clarft( char* direct, char* storev, int* n, int* k, scomplex* v, int* ldv, scomplex* tau, scomplex* t, int* ldt );
void F77_zlarft( char* direct, char* storev, int* n, int* k, dcomplex* v, int* ldv, dcomplex* tau, dcomplex* t, int* ldt );

// --- Generate a Householder vector (classic) ---

void F77_slarfg( int* n, float*    alpha, float*    x, int* incx, float*    tau );
void F77_dlarfg( int* n, double*   alpha, double*   x, int* incx, double*   tau );
void F77_clarfg( int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* tau );
void F77_zlarfg( int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* tau );

// --- Form Q from QR factorization ---

void F77_sorgqr( int* m, int* n, int* k, float*    a, int* lda, float*    tau, float*    work, int* lwork, int* info );
void F77_dorgqr( int* m, int* n, int* k, double*   a, int* lda, double*   tau, double*   work, int* lwork, int* info );
void F77_cungqr( int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info );
void F77_zungqr( int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info );

// --- Apply Q or Q' from QR factorization ---

void F77_sormqr( char* side, char* trans, int* m, int* n, int* k, float*    a, int* lda, float*    tau, float*    c, int* ldc, float*    work, int* lwork, int* info );
void F77_dormqr( char* side, char* trans, int* m, int* n, int* k, double*   a, int* lda, double*   tau, double*   c, int* ldc, double*   work, int* lwork, int* info );
void F77_cunmqr( char* side, char* trans, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info );
void F77_zunmqr( char* side, char* trans, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info );

void F77_sorm2r( char* side, char* trans, int* m, int* n, int* k, float*    a, int* lda, float*    tau, float*    c, int* ldc, float*    work, int* info );
void F77_dorm2r( char* side, char* trans, int* m, int* n, int* k, double*   a, int* lda, double*   tau, double*   c, int* ldc, double*   work, int* info );
void F77_cunm2r( char* side, char* trans, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* info );
void F77_zunm2r( char* side, char* trans, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* info );

// --- Apply Q or Q' from LQ factorization ---

void F77_sormlq( char* side, char* trans, int* m, int* n, int* k, float*    a, int* lda, float*    tau, float*    c, int* ldc, float*    work, int* lwork, int* info );
void F77_dormlq( char* side, char* trans, int* m, int* n, int* k, double*   a, int* lda, double*   tau, double*   c, int* ldc, double*   work, int* lwork, int* info );
void F77_cunmlq( char* side, char* trans, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info );
void F77_zunmlq( char* side, char* trans, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info );

void F77_sorml2( char* side, char* trans, int* m, int* n, int* k, float*    a, int* lda, float*    tau, float*    c, int* ldc, float*    work, int* info );
void F77_dorml2( char* side, char* trans, int* m, int* n, int* k, double*   a, int* lda, double*   tau, double*   c, int* ldc, double*   work, int* info );
void F77_cunml2( char* side, char* trans, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* info );
void F77_zunml2( char* side, char* trans, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* info );

// --- Form Q from tridiagonal reduction ---

void F77_sorgtr( char* uplo, int* m, float*    a, int* lda, float*    tau, float*    work, int* lwork, int* info );
void F77_dorgtr( char* uplo, int* m, double*   a, int* lda, double*   tau, double*   work, int* lwork, int* info );
void F77_cungtr( char* uplo, int* m, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info );
void F77_zungtr( char* uplo, int* m, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info );

// --- Apply Q or Q' from tridiagonal reduction ---

void F77_sormtr( char* side, char* uplo, char* trans, int* m, int* n, float*    a, int* lda, float*    tau, float*    c, int* ldc, float*    work, int* lwork, int* info );
void F77_dormtr( char* side, char* uplo, char* trans, int* m, int* n, double*   a, int* lda, double*   tau, double*   c, int* ldc, double*   work, int* lwork, int* info );
void F77_cunmtr( char* side, char* uplo, char* trans, int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info );
void F77_zunmtr( char* side, char* uplo, char* trans, int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info );

// --- Form Q from bidiagonal reduction ---

void F77_sorgbr( char* vect, int* m, int* n, int* k, float*    a, int* lda, float*    tau, float*    work, int* lwork, int* info );
void F77_dorgbr( char* vect, int* m, int* n, int* k, double*   a, int* lda, double*   tau, double*   work, int* lwork, int* info );
void F77_cungbr( char* vect, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info );
void F77_zungbr( char* vect, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info );

// --- Apply Q or Q' from bidiagonal reduction ---

void F77_sormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, float*    a, int* lda, float*    tau, float*    c, int* ldc, float*    work, int* lwork, int* info );
void F77_dormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, double*   a, int* lda, double*   tau, double*   c, int* ldc, double*   work, int* lwork, int* info );
void F77_cunmbr( char* vect, char* side, char* trans, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info );
void F77_zunmbr( char* vect, char* side, char* trans, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info );

// --- Tridiagonal QR algorithm ---

void F77_ssteqr( char* jobz, int* n, float*    d, float*    e, float*    z, int* ldz, float*  work, int* info ); 
void F77_dsteqr( char* jobz, int* n, double*   d, double*   e, double*   z, int* ldz, double* work, int* info ); 
void F77_csteqr( char* jobz, int* n, float*    d, float*    e, scomplex* z, int* ldz, float*  work, int* info ); 
void F77_zsteqr( char* jobz, int* n, double*   d, double*   e, dcomplex* z, int* ldz, double* work, int* info ); 

// --- Tridiagonal divide-and-conquer algorithm ---

void F77_sstedc( char* compz, int* n, float*    d, float*    e, float*    z, int* ldz, float*    work, int* lwork,                             int* iwork, int* liwork, int* info );
void F77_dstedc( char* compz, int* n, double*   d, double*   e, double*   z, int* ldz, double*   work, int* lwork,                             int* iwork, int* liwork, int* info );
void F77_cstedc( char* compz, int* n, float*    d, float*    e, scomplex* z, int* ldz, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info );
void F77_zstedc( char* compz, int* n, double*   d, double*   e, dcomplex* z, int* ldz, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info );

// --- Tridiagonal MRRR algorithm ---

void F77_sstemr( char* jobz, char* range, int* n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, float*    z, int* ldz, int* nzc, int* isuppz, int* tryrac, float*  work, int* lwork, int* iwork, int* liwork, int* info );
void F77_dstemr( char* jobz, char* range, int* n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, double*   z, int* ldz, int* nzc, int* isuppz, int* tryrac, double* work, int* lwork, int* iwork, int* liwork, int* info );
void F77_cstemr( char* jobz, char* range, int* n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, scomplex* z, int* ldz, int* nzc, int* isuppz, int* tryrac, float*  work, int* lwork, int* iwork, int* liwork, int* info );
void F77_zstemr( char* jobz, char* range, int* n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, dcomplex* z, int* ldz, int* nzc, int* isuppz, int* tryrac, double* work, int* lwork, int* iwork, int* liwork, int* info );

// --- Hermitian eigenvalue decomposition (QR algorithm) ---

void F77_ssyev( char* jobz, char* uplo, int* n, float*    a, int* lda, float*  w, float*    work, int* lwork, float*  rwork, int* info ); 
void F77_dsyev( char* jobz, char* uplo, int* n, double*   a, int* lda, double* w, double*   work, int* lwork, double* rwork, int* info ); 
void F77_cheev( char* jobz, char* uplo, int* n, scomplex* a, int* lda, float*  w, scomplex* work, int* lwork, float*  rwork, int* info ); 
void F77_zheev( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* info ); 

// --- Hermitian eigenvalue decomposition (divide-and-conquer) ---

void F77_ssyevd( char* jobz, char* uplo, int* n, float*    a, int* lda, float*  w, float*    work, int* lwork,                             int* iwork, int* liwork, int* info ); 
void F77_dsyevd( char* jobz, char* uplo, int* n, double*   a, int* lda, double* w, double*   work, int* lwork,                             int* iwork, int* liwork, int* info ); 
void F77_cheevd( char* jobz, char* uplo, int* n, scomplex* a, int* lda, float*  w, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info ); 
void F77_zheevd( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info ); 

// --- Hermitian eigenvalue decomposition (MRRR) ---

void F77_ssyevr( char* jobz, char* range, char* uplo, int* n, float*    a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, float*    z, int* ldz, int* isuppz, float*    work, int* lwork,                             int* iwork, int* liwork, int* info ); 
void F77_dsyevr( char* jobz, char* range, char* uplo, int* n, double*   a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, double*   z, int* ldz, int* isuppz, double*   work, int* lwork,                             int* iwork, int* liwork, int* info ); 
void F77_cheevr( char* jobz, char* range, char* uplo, int* n, scomplex* a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, scomplex* z, int* ldz, int* isuppz, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info ); 
void F77_zheevr( char* jobz, char* range, char* uplo, int* n, dcomplex* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, dcomplex* z, int* ldz, int* isuppz, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info ); 

// --- Bidiagonal QR algorithm ---

void F77_sbdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, float*    d, float*    e, float*    vt, int* ldvt, float*    u, int* ldu, float*    c, int* ldc, float*  rwork, int* info ); 
void F77_dbdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, double*   d, double*   e, double*   vt, int* ldvt, double*   u, int* ldu, double*   c, int* ldc, double* rwork, int* info ); 
void F77_cbdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, float*    d, float*    e, scomplex* vt, int* ldvt, scomplex* u, int* ldu, scomplex* c, int* ldc, float*  rwork, int* info ); 
void F77_zbdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, double*   d, double*   e, dcomplex* vt, int* ldvt, dcomplex* u, int* ldu, dcomplex* c, int* ldc, double* rwork, int* info ); 

// --- Bidiagonal divide-and-conquor algorithm ---

void F77_sbdsdc( char* uplo, char* compq, int* n, float*  d, float*  e, float*  u, int* ldu, float*  vt, int* ldvt, float*  q, float*  iq, float*  work, int* iwork, int* info ); 
void F77_dbdsdc( char* uplo, char* compq, int* n, double* d, double* e, double* u, int* ldu, double* vt, int* ldvt, double* q, double* iq, double* work, int* iwork, int* info ); 

// --- General matrix singular value decomposition (QR algorithm) ---

void F77_sgesvd( char* jobu, char* jobv, int* m, int* n, float*    a, int* lda, float*  s, float*    u, int* ldu, float*    vt, int* ldvt, float*    work, int* lwork,                int* info );
void F77_dgesvd( char* jobu, char* jobv, int* m, int* n, double*   a, int* lda, double* s, double*   u, int* ldu, double*   vt, int* ldvt, double*   work, int* lwork,                int* info );
void F77_cgesvd( char* jobu, char* jobv, int* m, int* n, scomplex* a, int* lda, float*  s, scomplex* u, int* ldu, scomplex* vt, int* ldvt, scomplex* work, int* lwork, float*  rwork, int* info );
void F77_zgesvd( char* jobu, char* jobv, int* m, int* n, dcomplex* a, int* lda, double* s, dcomplex* u, int* ldu, dcomplex* vt, int* ldvt, dcomplex* work, int* lwork, double* rwork, int* info );

// --- General matrix singular value decomposition (divide-and-conquer) ---

void F77_sgesdd( char* jobz, int* m, int* n, float*    a, int* lda, float*  s, float*    u, int* ldu, float*    vt, int* ldvt, float*    work, int* lwork,                int* iwork, int* info );
void F77_dgesdd( char* jobz, int* m, int* n, double*   a, int* lda, double* s, double*   u, int* ldu, double*   vt, int* ldvt, double*   work, int* lwork,                int* iwork, int* info );
void F77_cgesdd( char* jobz, int* m, int* n, scomplex* a, int* lda, float*  s, scomplex* u, int* ldu, scomplex* vt, int* ldvt, scomplex* work, int* lwork, float*  rwork, int* iwork, int* info );
void F77_zgesdd( char* jobz, int* m, int* n, dcomplex* a, int* lda, double* s, dcomplex* u, int* ldu, dcomplex* vt, int* ldvt, dcomplex* work, int* lwork, double* rwork, int* iwork, int* info );

// --- Swap rows ---

void F77_slaswp( int* n, float*    a, int* lda, int* k1, int* k2, int* ipiv, int* incx );
void F77_dlaswp( int* n, double*   a, int* lda, int* k1, int* k2, int* ipiv, int* incx );
void F77_claswp( int* n, scomplex* a, int* lda, int* k1, int* k2, int* ipiv, int* incx );
void F77_zlaswp( int* n, dcomplex* a, int* lda, int* k1, int* k2, int* ipiv, int* incx );

// --- Initialize a matrix ---

void F77_slaset( char* uplo, int* m, int* n, float*    alpha, float*    beta, float*    a, int* lda );
void F77_dlaset( char* uplo, int* m, int* n, double*   alpha, double*   beta, double*   a, int* lda );
void F77_claset( char* uplo, int* m, int* n, scomplex* alpha, scomplex* beta, scomplex* a, int* lda );
void F77_zlaset( char* uplo, int* m, int* n, dcomplex* alpha, dcomplex* beta, dcomplex* a, int* lda );

