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
#include "elemental/advanced_internal.hpp"
using namespace elemental;

template<typename F> // represents a real or complex number
void
elemental::advanced::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<F,MC,MR>& A  )
{
#ifndef RELEASE
    PushCallStack("advanced::Trinv");
#endif
    advanced::internal::TrinvVar3( shape, diagonal, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
elemental::advanced::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<F,MC,MR>& A  )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::TrinvVar3");
#endif
    if( shape == Lower )
        advanced::internal::TrinvLVar3( diagonal, A );
    else
        advanced::internal::TrinvUVar3( diagonal, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::advanced::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<float,MC,MR>& A );

template void elemental::advanced::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<double,MC,MR>& A );

template void elemental::advanced::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<scomplex,MC,MR>& A );

template void elemental::advanced::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal,
  DistMatrix<scomplex,MC,MR>& A );

template void elemental::advanced::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<dcomplex,MC,MR>& A );

template void elemental::advanced::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal,
  DistMatrix<dcomplex,MC,MR>& A );
#endif

