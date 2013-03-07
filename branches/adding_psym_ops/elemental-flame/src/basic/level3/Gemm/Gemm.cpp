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
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;

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

template<typename T>
void
elemental::basic::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::Gemm");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        basic::internal::GemmNN( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        basic::internal::GemmNT( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        basic::internal::GemmTN
        ( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        basic::internal::GemmTT
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmA");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        basic::internal::GemmNNA( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        basic::internal::GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        basic::internal::GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        basic::internal::GemmTTA
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmB");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        basic::internal::GemmNNB( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        basic::internal::GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        basic::internal::GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        basic::internal::GemmTTB
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmC");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        basic::internal::GemmNNC( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        basic::internal::GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        basic::internal::GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        basic::internal::GemmTTC
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmDot");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
        basic::internal::GemmNNDot( alpha, A, B, beta, C );
    else
        throw logic_error( "GemmDot currently only implemented for NN case." );
    // This code will be enabled when the routines are implemented
    /*
    else if( orientationOfA == Normal )
    {
        basic::internal::GemmNTDot( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        basic::internal::GemmTNDot( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        basic::internal::GemmTTDot( orientationOfA, orientationOfB,
                                   alpha, A, B, beta, C );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::basic::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::basic::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::basic::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::basic::internal::GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::basic::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::basic::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::basic::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::basic::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::basic::internal::GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::basic::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::basic::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::basic::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::basic::internal::GemmDot
( Orientation orientationOfA,
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::basic::Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::basic::internal::GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::basic::internal::GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::basic::internal::GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::basic::internal::GemmDot
( Orientation orientationOfA,
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
