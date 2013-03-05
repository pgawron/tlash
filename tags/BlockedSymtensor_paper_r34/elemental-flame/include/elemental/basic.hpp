/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental-FLAME.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_BASIC_HPP
#define ELEMENTAL_BASIC_HPP 1

#include "elemental/partitioning.hpp"

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

namespace elemental {
namespace basic {

//----------------------------------------------------------------------------//
// Tuning parameters                                                          //
//----------------------------------------------------------------------------//

template<typename T> void SetLocalHemvBlocksize( int blocksize );
template<> void SetLocalHemvBlocksize<float>( int blocksize );
template<> void SetLocalHemvBlocksize<double>( int blocksize );
#ifndef WITHOUT_COMPLEX
template<> void SetLocalHemvBlocksize< std::complex<float> >( int blocksize );
template<> void SetLocalHemvBlocksize< std::complex<double> >( int blocksize );
#endif // WITHOUT_COMPLEX

template<typename T> void SetLocalSymvBlocksize( int blocksize );
template<> void SetLocalSymvBlocksize<float>( int blocksize );
template<> void SetLocalSymvBlocksize<double>( int blocksize );
#ifndef WITHOUT_COMPLEX
template<> void SetLocalSymvBlocksize< std::complex<float> >( int blocksize );
template<> void SetLocalSymvBlocksize< std::complex<double> >( int blocksize );
#endif // WITHOUT_COMPLEX

template<typename T> void SetLocalTriangularRankKBlocksize( int blocksize );
template<> void SetLocalTriangularRankKBlocksize<float>( int blocksize );
template<> void SetLocalTriangularRankKBlocksize<double>( int blocksize );
#ifndef WITHOUT_COMPLEX
template<> void 
SetLocalTriangularRankKBlocksize< std::complex<float> >( int blocksize );
template<> void 
SetLocalTriangularRankKBlocksize< std::complex<double> >( int blocksize );
#endif // WITHOUT_COMPLEX

template<typename T> void SetLocalTriangularRank2KBlocksize( int blocksize );
template<> void SetLocalTriangularRank2KBlocksize<float>( int blocksize );
template<> void SetLocalTriangularRank2KBlocksize<double>( int blocksize );
#ifndef WITHOUT_COMPLEX
template<> void 
SetLocalTriangularRank2KBlocksize< std::complex<float> >( int blocksize );
template<> void 
SetLocalTriangularRank2KBlocksize< std::complex<double> >( int blocksize );
#endif // WITHOUT_COMPLEX

template<typename T> int LocalHemvBlocksize();
template<> int LocalHemvBlocksize<float>();
template<> int LocalHemvBlocksize<double>();
#ifndef WITHOUT_COMPLEX
template<> int LocalHemvBlocksize<scomplex>();
template<> int LocalHemvBlocksize<dcomplex>();
#endif // WITHOUT_COMPLEX

template<typename T> int LocalSymvBlocksize();
template<> int LocalSymvBlocksize<float>();
template<> int LocalSymvBlocksize<double>();
#ifndef WITHOUT_COMPLEX
template<> int LocalSymvBlocksize<scomplex>();
template<> int LocalSymvBlocksize<dcomplex>();
#endif // WITHOUT_COMPLEX

template<typename T> int LocalTriangularRankKBlocksize();
template<> int LocalTriangularRankKBlocksize<float>();
template<> int LocalTriangularRankKBlocksize<double>();
#ifndef WITHOUT_COMPLEX
template<> int LocalTriangularRankKBlocksize<scomplex>();
template<> int LocalTriangularRankKBlocksize<dcomplex>();
#endif // WITHOUT_COMPLEX

template<typename T> int LocalTriangularRank2KBlocksize();
template<> int LocalTriangularRank2KBlocksize<float>();
template<> int LocalTriangularRank2KBlocksize<double>();
#ifndef WITHOUT_COMPLEX
template<> int LocalTriangularRank2KBlocksize<scomplex>();
template<> int LocalTriangularRank2KBlocksize<dcomplex>();
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Level 1 BLAS-like functionality                                            //
//----------------------------------------------------------------------------//

//
// Axpy (Alpha X Plus Y):
//
// Y := alpha X + Y
//

// Serial version
template<typename T>
void
Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y );

// Parallel version
template<typename T, Distribution U, Distribution V>
void
Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );

//
// Copy:
//
// Y := X
//

// Serial version
template<typename T>
void
Copy( const Matrix<T>& X, Matrix<T>& Y );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z >
void
Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// Dot: 
// 
// Returns (x,y) = x^H y.
//
// Though the standard BLAS interface only defines DOT for real 
// datatypes, it is naturally generalized to an inner product over the
// complex field. Recall that the conjugate symmetry of inner products 
// requires that (x,y) = conj(y,x), so that (x,x) = conj( (x,x) ) => 
// (x,x) is real. This requires that we choose (x,x) = conj(x)^T * x.
//

// Serial version
template<typename T>
T
Dot( const Matrix<T>& x, const Matrix<T>& y );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z >
T
Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

//
// Dotc:
//
// Returns (x,y) = x^H y.
//
// This is the sister routine to 'Dot'; while 'Dot' is originally defined 
// only over the reals, 'Dotc' was defined only over the complex field. 
// They both have been extended to the same function, so from our point of 
// view they are identical.
//

// Serial version
template<typename T>
T
Dotc( const Matrix<T>& x, const Matrix<T>& y );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z >
T
Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

//
// Dotu: 
//
// Returns x^T y.
//
// Note: in the complex case, this is NOT an inner product.
//

// Serial version
template<typename T>
T
Dotu( const Matrix<T>& x, const Matrix<T>& y );

// Parallel version
template<typename T, Distribution U, Distribution V,
                     Distribution W, Distribution Z >
T
Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

//
// Nrm2 (2-norm):
//
// || x ||_2 = sqrt( x^H x ).
//

// Serial version for real datatypes
template<typename R>
R
Nrm2( const Matrix<R>& x ); 

#ifndef WITHOUT_COMPLEX
// Serial version for complex datatypes
template<typename R>
R
Nrm2( const Matrix< std::complex<R> >& x );
#endif

// Parallel version for real datatypes
template<typename R>
R
Nrm2( const DistMatrix<R,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
// Parallel version for complex datatypes
template<typename R>
R
Nrm2( const DistMatrix< std::complex<R>, MC, MR >& x );
#endif

// 
// Scal:
//
// X := alpha X
//

// Serial version
template<typename T>
void
Scal( T alpha, Matrix<T>& X );
    
// Parallel version
template<typename T, Distribution U, Distribution V>
void
Scal
( T alpha, DistMatrix<T,U,V>& A );
    
//----------------------------------------------------------------------------//
// Level 1 BLAS-like extensions                                               //
//----------------------------------------------------------------------------//

//
// Conj: 
//
// Conjugates a matrix. The in-place version performs A := Conj(A), while the 
// out-of-place sets B := Conj(A).
//

// In-place serial version for real datatypes. 
// Note: this is a no-op.
template<typename Z>
void
Conj( Matrix<Z>& A );

#ifndef WITHOUT_COMPLEX
// In-place serial version for complex datatypes.
template<typename Z>
void
Conj( Matrix< std::complex<Z> >& A );
#endif

// In-place parallel version
template<typename T, Distribution U, Distribution V>
void
Conj( DistMatrix<T,U,V>& A );

// Out-of-place serial version.
template<typename T>
void
Conj( const Matrix<T>& A, Matrix<T>& B );

// Out-of-place parallel version.
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void
Conj( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// ConjTrans:
//
// B := A^H
//

// Serial version
template<typename T>
void
ConjTrans( const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void
ConjTrans( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//
// Trans:
//
// B := A^T
//

// Serial version
template<typename T>
void
Trans( const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T, 
         Distribution U, Distribution V,
         Distribution W, Distribution Z>
void
Trans( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

//----------------------------------------------------------------------------//
// Level 2 BLAS-like routines                                                 //
//----------------------------------------------------------------------------//

//
// Gemv (GEneral Matrix-Vector multiply):
//
// y := alpha orientation( A ) x + beta y,
// where orientation( A ) is determined by 'orientation'.
//

// Serial version
template<typename T>
void
Gemv
( Orientation orientation, 
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void
Gemv
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Ger (GEneral Rank-one update):
//
// A := alpha x y^H + A
//

// Serial version
template<typename T>
void
Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Ger
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Gerc (GEneral Rank-one Conjugated update):
//
// A := alpha x y^H + A
//
// This is identical to Ger because both have been extended to work for both
// real and complex datatypes.
//

// Serial version
template<typename T>
void
Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Gerc
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Geru (GEneral Rank-one Unconjugated update):
//
// A := alpha x y^T + A
//

// Serial version
template<typename T>
void
Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Geru
( T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Hemv (HErmitian Matrix-Vector multiply):
//
// Implicitly performs
//   y := alpha A x + beta y,
// where only the triangle specified by 'shape' is referenced and the other
// triangle is implied by the Hermitian assumption.

// Serial version
template<typename T>
void
Hemv( Shape shape,
      T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void
Hemv
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Her (HErmitian Rank-one update):
//
// Implicitly performs
//   A := alpha x x^H + A,
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Her( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A );

// Parallel version
template<typename T>
void
Her
( Shape shape, T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A );

//
// Her2 (HErmitian Rank-2 update):
//
// Implicitly performs
//   A := alpha ( x y^H + y x^H ) + A,
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Her2
( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Her2
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Symv (SYmmetric Matrix-Vector multiply):
//
// Implicitly performs
//   y := alpha A x + beta y,
// where only the triangle specified by 'shape' is referenced and the other
// triangle is implied by the symmetry assumption.
//

// Serial version
template<typename T>
void
Symv
( Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

// Parallel version
template<typename T>
void
Symv
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

//
// Syr (SYmmetric Rank-one update):
//
// Implicitly performs the update
//   A := alpha x x^T + A,
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syr( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A );

// Parallel version
template<typename T>
void
Syr
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x, DistMatrix<T,MC,MR>& A );

//
// Syr2 (SYmmetric Rank-2 update):
//
// Implicitly perform the update
//   A := alpha ( x y^T + y x^T ) + A
// where only the triangle specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syr2
( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

// Parallel version
template<typename T>
void
Syr2
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x, const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A );

//
// Trmv (TRiangular Matrix-Vector multiply):
//
// Performs the update
//   x := orientation( A ) x,
// where 'shape' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diagonal' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename T>
void
Trmv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x );

// Parallel version
template<typename T>
void
Trmv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& x );

//
// Trsv (TRiangular Solve with a Vector):
//
// Performs the update
//   x := orientation( A )^-1 x,
// where 'shape' determines whether or not A is to be implicitly treated as 
// lower or upper triangular, and 'diagonal' specifies whether it has an 
// implicit unit diagonal.
//

// Serial version
template<typename F>
void
Trsv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<F>& A, Matrix<F>& x );

// Parallel version
template<typename F>
void
Trsv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& x );

//----------------------------------------------------------------------------//
// Level 3 BLAS-like routines                                                 //
//----------------------------------------------------------------------------//

//
// Gemm (GEneral Matrix-Matrix multiplication):
//
// C := alpha orientationOfA( A ) orientationOfB( B ) + beta C
//

// Serial version
template<typename T>
void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

template<typename T>
void
Strassen
( const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  DistMatrix<T,MC,MR>& C );


//
// Hemm (HErmitian Matrix-Matrix multiplication):
//
// Performs the update
//   C := alpha A B + beta C,  { side = Left }
// or
//   C := alpha B A + beta C,  { side = Right }
// where only the triangle of 'A' specified by 'shape' is referenced, and the
// other triangle is implied by the Hermitian assumption.
//

// Serial version
template<typename T>
void
Hemm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Hemm
( Side side, Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Her2k (HErmitian Rank-2K update):
//
// Performs the update
//   C := alpha ( A B^H + B A^H ) + beta C, { orientation = Normal }
// or
//   C := alpha ( A^H B + B^H A ) + beta C, { orientation = ConjugateTranspose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Her2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Her2k
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Herk (HErmitian Rank-K update):
//
// Performs the update
//   C := alpha A B^H + beta C,  { orientation = Normal }
// or
//   C := alpha A^H B + beta C,  { orientation = ConjugateTranspose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Herk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Herk
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

//
// Symm (SYmmetric Matrix-Matrix multiplication):
//
// Performs the update
//   C := alpha A B + beta C,  { side = Left }
// or
//   C := alpha B A + beta C,  { side = Right }
// where only the triangle of A specified by 'shape' is referenced, and the 
// other triangle is implied by the symmetry assumption.
//

// Serial version
template<typename T>
void
Symm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C ); 

// Parallel version
template<typename T>
void
Symm
( Side side, Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Syr2k (SYmmetric Rank-2K update):
//
// Performs the update
//   C := alpha ( A B^H + B A^H ) + beta C,  { orientation = Normal }
// or
//   C := alpha ( A^H B + B^H A ) + beta C,  { orientation = Transpose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syr2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Syr2k
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

//
// Syrk (SYmmetric Rank-K update):
//
// Performs the update
//   C := alpha A B^H + beta C,  { orientation = Normal }
// or
//   C := alpha A^H B + beta C,  { orientation = Transpose }
// where only the triangle of C specified by 'shape' is updated.
//

// Serial version
template<typename T>
void
Syrk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C );

// Parallel version
template<typename T>
void
Syrk
( Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

//
// Trmm (TRiangular Matrix-Matrix multiplication):
//
// Performs the update
//   B := alpha orientation( A ) B,  { side = Left }
// or
//   B := alpha B orientation( A ),  { side = Right }
// where 'shape' determines whether A is assumed to be upper or lower 
// triangular and 'diagonal' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename T>
void
Trmm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B );

// Parallel version
template<typename T>
void
Trmm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

//
// Trsm (TRiangular Solve with Multiple right-hand sides):
//
// Performs the update
//   B := alpha orientation( A )^-1 B,  { side = Left }
// or
//   B := alpha B orientation( A )^-1,  { side = Right }
// where 'shape' determines whether A is assumed to be upper or lower
// triangular and 'diagonal' determines whether A has an implicit unit
// diagonal.
//

// Serial version
template<typename F>
void
Trsm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  F alpha, const Matrix<F>& A, Matrix<F>& B ); 
        
// Parallel version
template<typename F>
void
Trsm
( Side side, Shape shape, Orientation orientation, Diagonal diagonal,
  F alpha, const DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B );

} // basic
} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 1                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::basic::Axpy
( T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("basic::Axpy");
#endif
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( (X.Height()==1 || X.Width()==1) && (Y.Height()==1 || Y.Width()==1) )
    {
        const unsigned XLength = ( X.Width()==1 ? X.Height() : X.Width() );
#ifndef RELEASE
        const unsigned YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
        if( XLength != YLength )
            throw std::logic_error( "Nonconformal Axpy." );
#endif
        if( X.Width()==1 && Y.Width()==1 )
        {
            elemental::imports::blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), 1 );
        }
        else if( X.Width()==1 )
        {
            elemental::imports::blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), Y.LDim() );
        }
        else if( Y.Width()==1 )
        {
            elemental::imports::blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), 1 );
        }
        else
        {
            elemental::imports::blas::Axpy
            ( XLength, alpha, 
              X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), Y.LDim() );
        }
    }
    else 
    {
#ifndef RELEASE
        if( X.Height() != Y.Height() || X.Width() != Y.Width() )
            throw std::logic_error( "Nonconformal Axpy." );
#endif
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                elemental::imports::blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                elemental::imports::blas::Axpy
                ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                    Y.Buffer(i,0),       Y.LDim() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Copy
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
elemental::basic::Dot
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Dot");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error( "Expected vector inputs." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error( "x and y must be the same length." );
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = imports::blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = imports::blas::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = imports::blas::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = imports::blas::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
elemental::basic::Dotc
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Dotc");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error( "Expected vector inputs." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error( "x and y must be the same length." );
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = imports::blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = imports::blas::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = imports::blas::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = imports::blas::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
elemental::basic::Dotu
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Dotu");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1) )
        throw std::logic_error( "Expected vector inputs." );
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
        throw std::logic_error( "x and y must be the same length." );
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = imports::blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = imports::blas::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = imports::blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = imports::blas::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename R>
inline R
elemental::basic::Nrm2
( const Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("basic::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "Expected vector input." );
#endif
    R norm;
    if( x.Width() == 1 )
    {
        norm = imports::blas::Nrm2
               ( x.Height(), x.LockedBuffer(), 1 );
    }
    else
    {
        norm = imports::blas::Nrm2
               ( x.Width(), x.LockedBuffer(), x.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
elemental::basic::Nrm2
( const Matrix< std::complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("basic::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "Expected vector input." );
#endif
    R norm;
    if( x.Width() == 1 )
    {
        norm = imports::blas::Nrm2
               ( x.Height(), x.LockedBuffer(), 1 );
    }
    else
    {
        norm = imports::blas::Nrm2
               ( x.Width(), x.LockedBuffer(), x.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif

template<typename T>
inline void
elemental::basic::Scal
( T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("basic::Scal");
#endif
    if( alpha != (T)1 )
    {
        if( alpha == (T)0 )
        {
            for( int j=0; j<X.Width(); ++j )
                for( int i=0; i<X.Height(); ++i )
                    X.Set(i,j,0);
        }
        else
        {
            for( int j=0; j<X.Width(); ++j )
            {
                imports::blas::Scal
                ( X.Height(), alpha, X.Buffer(0,j), 1 );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 1 (extensions)                             //
//----------------------------------------------------------------------------//

// Default case is for real datatypes
template<typename Z>
inline void
elemental::basic::Conj
( Matrix<Z>& A )
{ }

#ifndef WITHOUT_COMPLEX
// Specialization is to complex datatypes
template<typename Z>
inline void
elemental::basic::Conj
( Matrix< std::complex<Z> >& A )
{
#ifndef RELEASE
    PushCallStack("basic::Conj (in-place)");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set(i,j,elemental::Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template<typename T>
inline void
elemental::basic::Conj
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Conj");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(i,j,elemental::Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::ConjTrans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("basic::ConjTrans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( !B.Viewing() )
    {
        B.ResizeTo( n, m );
    }
    else if( B.Height() != n || B.Width() != m )
    {
        throw std::logic_error
              ("If ConjTrans'ing into a view, it must be the right size.");
    }
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(j,i,elemental::Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Trans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Trans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( !B.Viewing() )
    {
        B.ResizeTo( n, m );
    }
    else if( B.Height() != n || B.Width() != m )
    {
        throw std::logic_error
              ("If Transposing into a view, it must be the right size.");
    }
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(j,i,A.Get(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 2                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::basic::Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
    {
        std::ostringstream msg;
        msg << "x and y must be vectors: " << std::endl
            << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
            << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
        throw std::logic_error( msg.str() );
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == Normal )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:" << std::endl
                << "  A ~ " << A.Height() << " x " << A.Width() << std::endl
                << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
                << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
            throw std::logic_error( msg.str() );
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:" << std::endl
                << "  A ~ " << A.Height() << " x " << A.Width() << std::endl
                << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
                << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
            throw std::logic_error( msg.str() );
        }
    }
#endif
    const char transChar = OrientationToChar( orientation );
    const int m = A.Height();
    const int n = A.Width();
    const int k = ( transChar == 'N' ? n : m );
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( k != 0 )
    {
        imports::blas::Gemv
        ( transChar, m, n, 
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
          beta,  y.Buffer(), incy );
    }
    else
    {
        basic::Scal( beta, y );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Ger
( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Ger");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal Ger: " << std::endl
            << "  x ~ " << x.Height() << " x " << x.Width() << std::endl
            << "  y ~ " << y.Height() << " x " << y.Width() << std::endl
            << "  A ~ " << A.Height() << " x " << A.Width() << std::endl;
        throw std::logic_error( msg.str() );
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    imports::blas::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Gerc
( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Gerc");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error( "Nonconformal Gerc." );
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    imports::blas::Gerc
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Geru
( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Geru");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
        throw std::logic_error( "Nonconformal Geru." );
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    imports::blas::Geru
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Hemv
( Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Hemv");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error( "A must conform with x and y." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    imports::blas::Hemv
    ( uploChar, m,
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Her
( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Her");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error( "x must be a vector." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    imports::blas::Her
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Her2
( Shape shape,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Her2");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error( "x and y must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    imports::blas::Her2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Symv
( Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Symv");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error( "A must conform with x and y." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    imports::blas::Symv
    ( uploChar, m, 
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Syr
( Shape shape, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Syr");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error( "x must be a vector." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    imports::blas::Syr
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Syr2
( Shape shape, T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Syr2");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1) )
        throw std::logic_error( "x and y must be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
        throw std::logic_error( "x and y must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    imports::blas::Syr2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Trmv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x )
{
#ifndef RELEASE
    PushCallStack("basic::Trmv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "x must be a vector." );
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    imports::blas::Trmv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::basic::Trsv
( Shape shape, Orientation orientation, Diagonal diagonal,
  const Matrix<F>& A, Matrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("basic::Trsv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error( "x must be a vector." );
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error( "x must conform with A." );
#endif
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    imports::blas::Trsv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS-like routines: Level 3                                          //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::basic::Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Gemm");
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height() )
            throw std::logic_error( "Nonconformal GemmNN." );
    }
    else if( orientationOfA == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width() )
            throw std::logic_error( "Nonconformal GemmN(T/C)." );
    }
    else if( orientationOfB == Normal )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height() )
            throw std::logic_error( "Nonconformal Gemm(T/C)N." );
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width() )
            throw std::logic_error( "Nonconformal Gemm(T/C)(T/C)." );
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == Normal ? A.Width() : A.Height() );
    if( k != 0 )
    {
        imports::blas::Gemm
        ( transA, transB, m, n, k, 
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        basic::Scal( beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Hemm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Hemm");
#endif
    const char sideChar = SideToChar( side );
    const char shapeChar = ShapeToChar( shape );
    imports::blas::Hemm
    ( sideChar, shapeChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Her2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Her2k");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width() )
            throw std::logic_error( "Nonconformal Her2k." );
    }
    else if( orientation == ConjugateTranspose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width() )
            throw std::logic_error( "Nonconformal Her2k." );
    }
    else
        throw std::logic_error
        ( "Her2k only accepts Normal and ConjugateTranspose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    imports::blas::Her2k
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Herk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Herk");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error( "Nonconformal Herk." );
    }
    else if( orientation == ConjugateTranspose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error( "Nonconformal Herk." );
    }
    else
        throw std::logic_error
        ( "Herk only accepts Normal and ConjugateTranpose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    imports::blas::Herk
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Symm
( Side side, Shape shape,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Symm");
#endif
    const char sideChar = SideToChar( side );
    const char shapeChar = ShapeToChar( shape );
    imports::blas::Symm
    ( sideChar, shapeChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Syr2k
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Syr2k");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width()    )
            throw std::logic_error( "Nonconformal Syr2k." );
    }
    else if( orientation == Transpose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width()   )
            throw std::logic_error( "Nonconformal Syr2k." );
    }
    else
        throw std::logic_error
        ( "Syr2k only accepts Normal and Tranpose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    imports::blas::Syr2k
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Syrk
( Shape shape, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Syrk");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error( "Nonconformal Syrk." );
    }
    else if( orientation == Transpose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error( "Nonconformal Syrk." );
    }
    else
        throw std::logic_error
        ( "Syrk only accepts Normal and Tranpose options." );
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    imports::blas::Syrk
    ( uplo, trans, C.Height(), k, 
      alpha, A.LockedBuffer(), A.LDim(), 
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::Trmm
( Side side, Shape shape, 
  Orientation orientation, Diagonal diagonal,
  T alpha, const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Trmm");
    if( A.Height() != A.Width() )
        throw std::logic_error( "Triangular matrix must be square." );
    if( side == Left )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error( "Nonconformal Trmm." );
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error( "Nonconformal Trmm." );
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    imports::blas::Trmm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
elemental::basic::Trsm
( Side side, Shape shape,
  Orientation orientation,Diagonal diagonal,
  F alpha, const Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Trsm");
    if( A.Height() != A.Width() )
        throw std::logic_error( "Triangular matrix must be square." );
    if( side == Left )
    {
        if( A.Height() != B.Height() )
            throw std::logic_error( "Nonconformal Trsm." );
    }
    else
    {
        if( A.Height() != B.Width() )
            throw std::logic_error( "Nonconformal Trsm." );
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    imports::blas::Trsm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 1                                    //
//----------------------------------------------------------------------------//

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::basic::Axpy
( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("basic::Axpy");
    if( X.Grid() != Y.Grid() )
        throw std::logic_error
        ( "X and Y must be distributed over the same grid." );
    if( X.ColAlignment() != Y.ColAlignment() ||
        X.RowAlignment() != Y.RowAlignment() )
        throw std::logic_error( "Axpy requires X and Y be aligned." );
#endif
    basic::Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::basic::Copy
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

// Our extended Dotc is equivalent to our extended Dot, 
// but we are burdened with consistency
template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline T
elemental::basic::Dotc
( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Dotc");
#endif
    T globalDot = basic::Dot( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::basic::Scal
( T alpha, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Scal");
#endif
    basic::Scal( alpha, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 1 (extensions)                       //
//----------------------------------------------------------------------------//

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::basic::Conj
( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Conj (in-place)");
#endif
    basic::Conj( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::basic::Conj
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Conj");
#endif
    B = A;
    basic::Conj( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::basic::ConjTrans
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("basic::ConjTrans");
#endif
    DistMatrix<T,Z,W> C( B.Grid() );
    if( B.ConstrainedColAlignment() )
        C.AlignRowsWith( B );
    if( B.ConstrainedRowAlignment() )
        C.AlignColsWith( B );

    C = A;

    if( !B.Viewing() )
    {
        if( !B.ConstrainedColAlignment() )
            B.AlignColsWith( C );
        if( !B.ConstrainedRowAlignment() )
            B.AlignRowsWith( C );

        B.ResizeTo( A.Width(), A.Height() );
    }
    else if( B.Height() != A.Width() || B.Width() != A.Height() )
    {
        throw std::logic_error
              ("If ConjTrans'ing into a view, it must be the right size.");
    }
    basic::ConjTrans( C.LockedLocalMatrix(), B.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V,
                     elemental::Distribution W, elemental::Distribution Z >
inline void
elemental::basic::Trans
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("basic::Trans");
#endif
    DistMatrix<T,Z,W> C( B.Grid() );
    if( B.ConstrainedColAlignment() )
        C.AlignRowsWith( B );
    if( B.ConstrainedRowAlignment() )
        C.AlignColsWith( B );

    C = A;

    if( !B.Viewing() )
    {
        if( !B.ConstrainedColAlignment() )
            B.AlignColsWith( C );
        if( !B.ConstrainedRowAlignment() )
            B.AlignRowsWith( C );

        B.ResizeTo( A.Width(), A.Height() );
    }
    else if( B.Height() != A.Width() || B.Width() != A.Height() )
    {
        throw std::logic_error
              ("If Transposing into a view, it must be the right size.");
    }
    basic::Trans( C.LockedLocalMatrix(), B.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS-like routines: Level 2                                    //
//----------------------------------------------------------------------------//

// Our extended Ger and Gerc are equivalent
template<typename T>
inline void
elemental::basic::Gerc
( T alpha, const DistMatrix<T,MC,MR>& x,
           const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Gerc");
#endif
    basic::Ger( alpha, x, y, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif /* ELEMENTAL_BASIC_HPP */

