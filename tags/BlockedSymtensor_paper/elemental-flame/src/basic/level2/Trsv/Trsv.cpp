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

template<typename F>
void
elemental::basic::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<F,MC,MR>& A,
        DistMatrix<F,MC,MR>& x   )
{
#ifndef RELEASE
    PushCallStack("basic::Trsv");
#endif
    if( shape == Lower )
    {
        if( orientation == Normal )
            basic::internal::TrsvLN( diagonal, A, x );
        else
            basic::internal::TrsvLT( orientation, diagonal, A, x );
    }
    else
    {
        if( orientation == Normal )
            basic::internal::TrsvUN( diagonal, A, x );
        else
            basic::internal::TrsvUT( orientation, diagonal, A, x );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::basic::Trsv
( Shape shape, 
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& x );

template void
elemental::basic::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void
elemental::basic::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& x );

template void
elemental::basic::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

