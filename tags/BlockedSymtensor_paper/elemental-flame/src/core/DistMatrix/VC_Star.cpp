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
#include "elemental/dist_matrix.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;
using namespace elemental::utilities;

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

template<typename Z>
void
elemental::DistMatrix<Z,VC,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int p           = this->Grid().Size();
    const int colShift    = this->ColShift();

    this->SetToRandom();

    Z* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*p;
        if( i < width )
            thisLocalBuffer[iLocal+i*thisLDim] += width;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename Z>
void
elemental::DistMatrix<complex<Z>,VC,Star>::SetToRandomHPD()
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetToRandomHPD");
    this->AssertNotLockedView();
    if( this->Height() != this->Width() )
        throw logic_error( "Positive-definite matrices must be square." );
#endif
    const int width       = this->Width();
    const int localHeight = this->LocalHeight();
    const int p           = this->Grid().Size();
    const int colShift    = this->ColShift();

    this->SetToRandom();

    complex<Z>* thisLocalBuffer = this->LocalBuffer();
    const int thisLDim = this->LocalLDim();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*p;
        if( i < width )
        {
            const Z value = real(thisLocalBuffer[iLocal+i*thisLDim]);
            thisLocalBuffer[iLocal+i*thisLDim] = value + width;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,VC,Star>::GetReal
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetReal");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetRealLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
Z
elemental::DistMatrix<complex<Z>,VC,Star>::GetImag
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::GetImag");
    this->AssertValidEntry( i, j );
#endif
    // We will determine the owner rank of entry (i,j) and broadcast from that
    // process over the entire g
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    Z u;
    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        u = this->GetImagLocalEntry(iLoc,j);
    }
    mpi::Broadcast( &u, 1, ownerRank, g.VCComm() );

#ifndef RELEASE
    PopCallStack();
#endif
    return u;
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,VC,Star>::SetReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetReal");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        this->SetRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,VC,Star>::SetImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::SetImag");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        this->SetImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,VC,Star>::UpdateReal
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateReal");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateRealLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
void
elemental::DistMatrix<complex<Z>,VC,Star>::UpdateImag
( int i, int j, Z u )
{
#ifndef RELEASE
    PushCallStack("[VC,* ]::UpdateImag");
    this->AssertValidEntry( i, j );
#endif
    const elemental::Grid& g = this->Grid();
    const int ownerRank = (i + this->ColAlignment()) % g.Size();

    if( g.VCRank() == ownerRank )
    {
        const int iLoc = (i-this->ColShift()) / g.Size();
        this->UpdateImagLocalEntry(iLoc,j,u);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template class elemental::DistMatrix<int,   VC,Star>;
template class elemental::DistMatrix<float, VC,Star>;
template class elemental::DistMatrix<double,VC,Star>;
#ifndef WITHOUT_COMPLEX
template class elemental::DistMatrix<scomplex,VC,Star>;
template class elemental::DistMatrix<dcomplex,VC,Star>;
#endif

