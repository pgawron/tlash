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
using namespace std;
using namespace elemental;
using namespace elemental::imports;

// Algorithmic controls
namespace {
advanced::internal::TridiagApproach tridiagApproach = 
    advanced::internal::TRIDIAG_NORMAL;
advanced::internal::GridOrder gridOrder = 
    advanced::internal::ROW_MAJOR;
}
void 
elemental::advanced::internal::SetTridiagApproach
( advanced::internal::TridiagApproach approach )
{ ::tridiagApproach = approach; }
void 
elemental::advanced::internal::SetTridiagSquareGridOrder
( advanced::internal::GridOrder order )
{ ::gridOrder = order; }

template<typename R> // representation of a real number
void
elemental::advanced::Tridiag
( Shape shape, DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::Tridiag");
#endif
    const Grid& g = A.Grid();
    if( ::tridiagApproach == advanced::internal::TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( shape == Lower )
            advanced::internal::TridiagL( A );
        else 
            advanced::internal::TridiagU( A );
    }
    else if( ::tridiagApproach == advanced::internal::TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh
        int p = g.Size();
        int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

        std::vector<int> squareRanks(pSqrt*pSqrt);
        if( ::gridOrder == advanced::internal::COL_MAJOR )
        {
            for( int j=0; j<pSqrt; ++j )
                for( int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = i+j*pSqrt;
        }
        else
        {
            for( int j=0; j<pSqrt; ++j )
                for( int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = j+i*pSqrt;
        }

        mpi::Group owningGroup = g.OwningGroup();
        mpi::Group squareGroup;
        mpi::GroupIncl
        ( owningGroup, squareRanks.size(), &squareRanks[0], squareGroup );

        mpi::Comm viewingComm = g.ViewingComm();
        const Grid squareGrid( viewingComm, squareGroup, pSqrt, pSqrt );
        DistMatrix<R,MC,MR> ASquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( shape == Lower )
            advanced::internal::TridiagLSquare( ASquare );
        else
            advanced::internal::TridiagUSquare( ASquare ); 
        A = ASquare;

        mpi::GroupFree( squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( shape == Lower )
                advanced::internal::TridiagLSquare( A );
            else
                advanced::internal::TridiagUSquare( A );
        }
        else
        {
            if( shape == Lower )
                advanced::internal::TridiagL( A );
            else
                advanced::internal::TridiagU( A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
void
elemental::advanced::Tridiag
( Shape shape, 
  DistMatrix<complex<R>,MC,  MR  >& A,
  DistMatrix<complex<R>,Star,Star>& t )
{
#ifndef RELEASE
    PushCallStack("advanced::Tridiag");
#endif
    typedef complex<R> C;

    const Grid& g = A.Grid();
    if( ::tridiagApproach == advanced::internal::TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( shape == Lower )
            advanced::internal::TridiagL( A, t );
        else
            advanced::internal::TridiagU( A, t );
    }
    else if( ::tridiagApproach == advanced::internal::TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh 
        int p = g.Size();
        int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

        std::vector<int> squareRanks(pSqrt*pSqrt);
        if( ::gridOrder == advanced::internal::COL_MAJOR )
        {
            for( int j=0; j<pSqrt; ++j )
                for( int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = i+j*pSqrt;
        }
        else
        {
            for( int j=0; j<pSqrt; ++j )
                for( int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = j+i*pSqrt;
        }

        mpi::Group owningGroup = g.OwningGroup();
        mpi::Group squareGroup;
        mpi::GroupIncl
        ( owningGroup, squareRanks.size(), &squareRanks[0], squareGroup );

        mpi::Comm viewingComm = g.ViewingComm();
        const Grid squareGrid( viewingComm, squareGroup, pSqrt, pSqrt );
        DistMatrix<C,MC,MR> ASquare(squareGrid);
        DistMatrix<C,Star,Star> tSquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( shape == Lower )
            advanced::internal::TridiagLSquare( ASquare, tSquare );
        else
            advanced::internal::TridiagUSquare( ASquare, tSquare ); 
        A = ASquare;
        t = tSquare;

        mpi::GroupFree( squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( shape == Lower )
                advanced::internal::TridiagLSquare( A, t );
            else
                advanced::internal::TridiagUSquare( A, t ); 
        }
        else
        {
            if( shape == Lower )
                advanced::internal::TridiagL( A, t );
            else
                advanced::internal::TridiagU( A, t );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::advanced::Tridiag
( Shape shape, 
  DistMatrix<float,MC,MR>& A );

template void elemental::advanced::Tridiag
( Shape shape, 
  DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::Tridiag
( Shape shape,
  DistMatrix<scomplex,MC,  MR  >& A,
  DistMatrix<scomplex,Star,Star>& t );

template void elemental::advanced::Tridiag
( Shape shape,
  DistMatrix<dcomplex,MC,  MR  >& A,
  DistMatrix<dcomplex,Star,Star>& t );
#endif

