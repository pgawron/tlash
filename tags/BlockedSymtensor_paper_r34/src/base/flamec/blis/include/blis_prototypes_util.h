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

// --- Utility-level BLAS-like prototypes --------------------------------------

// --- constant-generating functions ---

float    bli_s2( void );
double   bli_d2( void );
scomplex bli_c2( void );
dcomplex bli_z2( void );
float    bli_s1( void );
double   bli_d1( void );
scomplex bli_c1( void );
dcomplex bli_z1( void );
float    bli_s1h( void );
double   bli_d1h( void );
scomplex bli_c1h( void );
dcomplex bli_z1h( void );
float    bli_s0( void );
double   bli_d0( void );
scomplex bli_c0( void );
dcomplex bli_z0( void );
float    bli_sm1h( void );
double   bli_dm1h( void );
scomplex bli_cm1h( void );
dcomplex bli_zm1h( void );
float    bli_sm1( void );
double   bli_dm1( void );
scomplex bli_cm1( void );
dcomplex bli_zm1( void );
float    bli_sm2( void );
double   bli_dm2( void );
scomplex bli_cm2( void );
dcomplex bli_zm2( void );

// --- allocv ---

void*     bli_vallocv( unsigned int n_elem, unsigned int elem_size );
int*      bli_iallocv( unsigned int n_elem );
float*    bli_sallocv( unsigned int n_elem );
double*   bli_dallocv( unsigned int n_elem );
scomplex* bli_callocv( unsigned int n_elem );
dcomplex* bli_zallocv( unsigned int n_elem );

// --- allocm ---

void*     bli_vallocm( unsigned int m, unsigned int n, unsigned int elem_size );
int*      bli_iallocm( unsigned int m, unsigned int n );
float*    bli_sallocm( unsigned int m, unsigned int n );
double*   bli_dallocm( unsigned int m, unsigned int n );
scomplex* bli_callocm( unsigned int m, unsigned int n );
dcomplex* bli_zallocm( unsigned int m, unsigned int n );

// --- apdiagmv ---

void bli_sapdiagmv( side_t side, conj_t conj, int m, int n, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bli_dapdiagmv( side_t side, conj_t conj, int m, int n, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bli_csapdiagmv( side_t side, conj_t conj, int m, int n, float*    x, int incx, scomplex* a, int a_rs, int a_cs );
void bli_capdiagmv( side_t side, conj_t conj, int m, int n, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bli_zdapdiagmv( side_t side, conj_t conj, int m, int n, double*   x, int incx, dcomplex* a, int a_rs, int a_cs );
void bli_zapdiagmv( side_t side, conj_t conj, int m, int n, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

// --- create_contigm ---

void bli_screate_contigm( int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dcreate_contigm( int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_ccreate_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zcreate_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmt ---

void bli_screate_contigmt( trans_t trans_dims, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dcreate_contigmt( trans_t trans_dims, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_ccreate_contigmt( trans_t trans_dims, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zcreate_contigmt( trans_t trans_dims, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmr ---

void bli_screate_contigmr( uplo_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dcreate_contigmr( uplo_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_ccreate_contigmr( uplo_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zcreate_contigmr( uplo_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmsr ---

void bli_screate_contigmsr( side_t side, uplo_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dcreate_contigmsr( side_t side, uplo_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_ccreate_contigmsr( side_t side, uplo_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zcreate_contigmsr( side_t side, uplo_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_contigm ---

void bli_sfree_contigm( float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dfree_contigm( double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_cfree_contigm( scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zfree_contigm( dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigm ---

void bli_sfree_saved_contigm( int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dfree_saved_contigm( int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_cfree_saved_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zfree_saved_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigmr ---

void bli_sfree_saved_contigmr( uplo_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dfree_saved_contigmr( uplo_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_cfree_saved_contigmr( uplo_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zfree_saved_contigmr( uplo_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigmsr ---

void bli_sfree_saved_contigmsr( side_t side, uplo_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bli_dfree_saved_contigmsr( side_t side, uplo_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bli_cfree_saved_contigmsr( side_t side, uplo_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bli_zfree_saved_contigmsr( side_t side, uplo_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- ewinvscalv ---

void bli_sewinvscalv( conj_t conj, int n, float*    x, int incx, float*    y, int incy );
void bli_dewinvscalv( conj_t conj, int n, double*   x, int incx, double*   y, int incy );
void bli_csewinvscalv( conj_t conj, int n, float*    x, int incx, scomplex* y, int incy );
void bli_cewinvscalv( conj_t conj, int n, scomplex* x, int incx, scomplex* y, int incy );
void bli_zdewinvscalv( conj_t conj, int n, double*   x, int incx, dcomplex* y, int incy );
void bli_zewinvscalv( conj_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- ewscalmt ---

void bli_sewinvscalmt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dewinvscalmt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_csewinvscalmt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_cewinvscalmt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zdewinvscalmt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zewinvscalmt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- ewscalv ---

void bli_sewscalv( conj_t conj, int n, float*    x, int incx, float*    y, int incy );
void bli_dewscalv( conj_t conj, int n, double*   x, int incx, double*   y, int incy );
void bli_csewscalv( conj_t conj, int n, float*    x, int incx, scomplex* y, int incy );
void bli_cewscalv( conj_t conj, int n, scomplex* x, int incx, scomplex* y, int incy );
void bli_zdewscalv( conj_t conj, int n, double*   x, int incx, dcomplex* y, int incy );
void bli_zewscalv( conj_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- ewscalmt ---

void bli_sewscalmt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dewscalmt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_csewscalmt( trans_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_cewscalmt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_zdewscalmt( trans_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bli_zewscalmt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- free ---

void bli_vfree( void*     p );
void bli_ifree( int*      p );
void bli_sfree( float*    p );
void bli_dfree( double*   p );
void bli_cfree( scomplex* p );
void bli_zfree( dcomplex* p );

// --- inverts ---

void bli_sinverts( conj_t conj, float*    alpha );
void bli_dinverts( conj_t conj, double*   alpha );
void bli_cinverts( conj_t conj, scomplex* alpha );
void bli_zinverts( conj_t conj, dcomplex* alpha );

// --- invert2s ---

void bli_sinvert2s( conj_t conj, float*    alpha, float*    beta );
void bli_dinvert2s( conj_t conj, double*   alpha, double*   beta );
void bli_cinvert2s( conj_t conj, scomplex* alpha, scomplex* beta );
void bli_zinvert2s( conj_t conj, dcomplex* alpha, dcomplex* beta );

// --- invertv ---

void bli_sinvertv( conj_t conj, int n, float*    x, int incx );
void bli_dinvertv( conj_t conj, int n, double*   x, int incx );
void bli_cinvertv( conj_t conj, int n, scomplex* x, int incx );
void bli_zinvertv( conj_t conj, int n, dcomplex* x, int incx );

// --- ident ---

void bli_sident( int m, float*    a, int a_rs, int a_cs );
void bli_dident( int m, double*   a, int a_rs, int a_cs );
void bli_cident( int m, scomplex* a, int a_rs, int a_cs );
void bli_zident( int m, dcomplex* a, int a_rs, int a_cs );

// --- maxabsv ---

void bli_smaxabsv( int n, float*    x, int incx, float*  maxabs );
void bli_dmaxabsv( int n, double*   x, int incx, double* maxabs );
void bli_cmaxabsv( int n, scomplex* x, int incx, float*  maxabs );
void bli_zmaxabsv( int n, dcomplex* x, int incx, double* maxabs );

// --- maxabsm ---

void bli_smaxabsm( int m, int n, float*    a, int a_rs, int a_cs, float*  maxabs );
void bli_dmaxabsm( int m, int n, double*   a, int a_rs, int a_cs, double* maxabs );
void bli_cmaxabsm( int m, int n, scomplex* a, int a_rs, int a_cs, float*  maxabs );
void bli_zmaxabsm( int m, int n, dcomplex* a, int a_rs, int a_cs, double* maxabs );

// --- maxabsmr ---

void bli_smaxabsmr( uplo_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*  maxabs );
void bli_dmaxabsmr( uplo_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double* maxabs );
void bli_cmaxabsmr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float*  maxabs );
void bli_zmaxabsmr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double* maxabs );

// --- rands ---

void bli_srands( float*    alpha );
void bli_drands( double*   alpha );
void bli_crands( scomplex* alpha );
void bli_zrands( dcomplex* alpha );

// --- randv ---

void bli_srandv( int n, float*    x, int incx );
void bli_drandv( int n, double*   x, int incx );
void bli_crandv( int n, scomplex* x, int incx );
void bli_zrandv( int n, dcomplex* x, int incx );

// --- randm ---

void bli_srandm( int m, int n, float*    a, int a_rs, int a_cs );
void bli_drandm( int m, int n, double*   a, int a_rs, int a_cs );
void bli_crandm( int m, int n, scomplex* a, int a_rs, int a_cs );
void bli_zrandm( int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- randmr ---
void bli_srandmr( uplo_t uplo, diag_t diag, int m, int n, float*    a, int a_rs, int a_cs );
void bli_drandmr( uplo_t uplo, diag_t diag, int m, int n, double*   a, int a_rs, int a_cs );
void bli_crandmr( uplo_t uplo, diag_t diag, int m, int n, scomplex* a, int a_rs, int a_cs );
void bli_zrandmr( uplo_t uplo, diag_t diag, int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- set_contig_strides ---

void bli_set_contig_strides( int m, int n, int* rs, int* cs );

// --- set_dims_with_side ---

void bli_set_dim_with_side( side_t side, int m, int n, int* dim_new );

// --- set_dims_with_trans ---

void bli_set_dims_with_trans( trans_t trans, int m, int n, int* m_new, int* n_new );

// --- setv ---

void bli_isetv( int m, int*      sigma, int*      x, int incx );
void bli_ssetv( int m, float*    sigma, float*    x, int incx );
void bli_dsetv( int m, double*   sigma, double*   x, int incx );
void bli_csetv( int m, scomplex* sigma, scomplex* x, int incx );
void bli_zsetv( int m, dcomplex* sigma, dcomplex* x, int incx );

// --- setm ---

void bli_isetm( int m, int n, int*      sigma, int*      a, int a_rs, int a_cs );
void bli_ssetm( int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bli_dsetm( int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bli_csetm( int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bli_zsetm( int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- setmr ---

void bli_ssetmr( uplo_t uplo, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bli_dsetmr( uplo_t uplo, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bli_csetmr( uplo_t uplo, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bli_zsetmr( uplo_t uplo, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- setdiag ---

void bli_isetdiag( int offset, int m, int n, int*      sigma, int*      a, int a_rs, int a_cs );
void bli_ssetdiag( int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bli_dsetdiag( int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bli_csetdiag( int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bli_zsetdiag( int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- scalediag ---

void bli_sscalediag( conj_t conj, int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bli_dscalediag( conj_t conj, int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bli_cscalediag( conj_t conj, int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bli_zscalediag( conj_t conj, int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );
void bli_csscalediag( conj_t conj, int offset, int m, int n, float*    sigma, scomplex* a, int a_rs, int a_cs );
void bli_zdscalediag( conj_t conj, int offset, int m, int n, double*   sigma, dcomplex* a, int a_rs, int a_cs );

// --- shiftdiag ---

void bli_sshiftdiag( conj_t conj, int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bli_dshiftdiag( conj_t conj, int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bli_cshiftdiag( conj_t conj, int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bli_zshiftdiag( conj_t conj, int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );
void bli_csshiftdiag( conj_t conj, int offset, int m, int n, float*    sigma, scomplex* a, int a_rs, int a_cs );
void bli_zdshiftdiag( conj_t conj, int offset, int m, int n, double*   sigma, dcomplex* a, int a_rs, int a_cs );

// --- symmize ---

void bli_ssymmize( conj_t conj, uplo_t uplo, int m, float*    a, int a_rs, int a_cs );
void bli_dsymmize( conj_t conj, uplo_t uplo, int m, double*   a, int a_rs, int a_cs );
void bli_csymmize( conj_t conj, uplo_t uplo, int m, scomplex* a, int a_rs, int a_cs );
void bli_zsymmize( conj_t conj, uplo_t uplo, int m, dcomplex* a, int a_rs, int a_cs );

