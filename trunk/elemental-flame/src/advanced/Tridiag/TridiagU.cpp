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
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;

template<typename R> // representation of a real number
void
elemental::advanced::internal::TridiagU
( DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::TridiagU");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
#endif
    const Grid& g = A.Grid();

    if( g.InGrid() )
    {
        // Matrix views 
        DistMatrix<R,MC,MR> 
            ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
            ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                             A20(g), A21(g), A22(g);

        // Temporary distributions
        DistMatrix<R,Star,Star> A11_Star_Star(g);
        DistMatrix<R,MC,  Star> APan_MC_Star(g);
        DistMatrix<R,MR,  Star> APan_MR_Star(g);
        DistMatrix<R,MC,  Star> A11_MC_Star(g);
        DistMatrix<R,MR,  Star> A11_MR_Star(g);
        DistMatrix<R,MC,  Star> A01_MC_Star(g);
        DistMatrix<R,MR,  Star> A01_MR_Star(g);
        DistMatrix<R,MC,  MR  > WPan(g);
        DistMatrix<R,MC,  Star> WPan_MC_Star(g);
        DistMatrix<R,MR,  Star> WPan_MR_Star(g);
        DistMatrix<R,MC,  Star> W11_MC_Star(g);
        DistMatrix<R,MR,  Star> W11_MR_Star(g);
        DistMatrix<R,MC,  Star> W01_MC_Star(g);
        DistMatrix<R,MR,  Star> W01_MR_Star(g);

        PartitionUpDiagonal
        ( A, ATL, ATR,
             ABL, ABR, 0 );
        while( ABR.Height() < A.Height() )
        {
            RepartitionUpDiagonal
            ( ATL, /**/ ATR,  A00, A01, /**/ A02,
                   /**/       A10, A11, /**/ A12,
             /*************/ /******************/
              ABL, /**/ ABR,  A20, A21, /**/ A22 );

            if( A00.Height() > 0 )
            {
                APan_MC_Star.AlignWith( A00 );
                APan_MR_Star.AlignWith( A00 );
                APan_MC_Star.ResizeTo( ATL.Height(), A11.Width() );
                APan_MR_Star.ResizeTo( ATL.Height(), A11.Width() );
                WPan.AlignWith( A01 );
                WPan_MC_Star.AlignWith( A00 );
                WPan_MR_Star.AlignWith( A00 );
                WPan.ResizeTo( ATL.Height(), A11.Width() );
                WPan_MC_Star.ResizeTo( ATL.Height(), A11.Width() );
                WPan_MR_Star.ResizeTo( ATL.Height(), A11.Width() );
                PartitionUp
                ( APan_MC_Star, A01_MC_Star,
                                A11_MC_Star, A11.Height() );
                PartitionUp
                ( APan_MR_Star, A01_MR_Star,
                                A11_MR_Star, A11.Height() );
                PartitionUp
                ( WPan_MC_Star, W01_MC_Star,
                                W11_MC_Star, A11.Height() );
                PartitionUp
                ( WPan_MR_Star, W01_MR_Star,
                                W11_MR_Star, A11.Height() );
                //------------------------------------------------------------//
                // Accumulate the Householder vectors into A01 and form W01 
                // such that subtracting (A01 W01' + W01 A01') is equal to 
                // successively applying the similarity transformations 
                // (I-tau h h')A00(I-tau h h') for each (tau,h).
                //
                // APan[MC,* ], APan[MR,* ], WPan[MC,* ], and WPan[MR,* ] are 
                // formed during the panel factorization.
                advanced::internal::PanelTridiagU
                ( ATL, WPan, 
                  APan_MC_Star, APan_MR_Star, WPan_MC_Star, WPan_MR_Star );
                basic::internal::LocalTriangularRank2K
                ( Upper, Transpose, Transpose,
                  (R)-1, A01_MC_Star, W01_MC_Star, A01_MR_Star, W01_MR_Star,
                  (R)1, A00 );
                //------------------------------------------------------------//
                APan_MC_Star.FreeAlignments();
                APan_MR_Star.FreeAlignments();
                WPan.FreeAlignments();
                WPan_MC_Star.FreeAlignments();
                WPan_MR_Star.FreeAlignments();
            }
            else
            {
                A11_Star_Star = A11;
                advanced::Tridiag( Upper, A11_Star_Star.LocalMatrix() );
                A11 = A11_Star_Star;
            }

            SlidePartitionUpDiagonal
            ( ATL, /**/ ATR,  A00, /**/ A01, A02,
             /*************/ /******************/
                   /**/       A10, /**/ A11, A12,
              ABL, /**/ ABR,  A20, /**/ A21, A22 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
void
elemental::advanced::internal::TridiagU
( DistMatrix<complex<R>,MC,  MR  >& A,
  DistMatrix<complex<R>,Star,Star>& t )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::TridiagU");
    if( A.Grid() != t.Grid() )
        throw logic_error( "A and t must be distributed over the same grid." );
#endif
    const Grid& g = A.Grid();
#ifndef RELEASE
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( t.Viewing() )
        throw logic_error( "t must not be a view." );
#endif
    typedef complex<R> C;

    DistMatrix<C,MD,Star> tDiag(g);
    tDiag.AlignWithDiag( A, 1 );
    tDiag.ResizeTo( A.Height()-1, 1 );

    if( g.InGrid() )
    {
        // Matrix views 
        DistMatrix<C,MC,MR> 
            ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
            ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                             A20(g), A21(g), A22(g);
        DistMatrix<C,MD,Star> tT(g),  t0(g), 
                              tB(g),  t1(g),
                                      t2(g);

        // Temporary distributions
        DistMatrix<C,Star,Star> A11_Star_Star(g);
        DistMatrix<C,MC,  Star> APan_MC_Star(g);
        DistMatrix<C,MR,  Star> APan_MR_Star(g);
        DistMatrix<C,MC,  Star> A11_MC_Star(g);
        DistMatrix<C,MR,  Star> A11_MR_Star(g);
        DistMatrix<C,MC,  Star> A01_MC_Star(g);
        DistMatrix<C,MR,  Star> A01_MR_Star(g);
        DistMatrix<C,MC,  MR  > WPan(g);
        DistMatrix<C,MC,  Star> WPan_MC_Star(g);
        DistMatrix<C,MR,  Star> WPan_MR_Star(g);
        DistMatrix<C,MC,  Star> W11_MC_Star(g);
        DistMatrix<C,MR,  Star> W11_MR_Star(g);
        DistMatrix<C,MC,  Star> W01_MC_Star(g);
        DistMatrix<C,MR,  Star> W01_MR_Star(g);
        DistMatrix<C,Star,Star> t1_Star_Star(g);

        PartitionUpDiagonal
        ( A, ATL, ATR,
             ABL, ABR, 0 );
        PartitionUp
        ( tDiag, tT,
                 tB, 0 );
        while( ABR.Height() < A.Height() )
        {
            RepartitionUpDiagonal
            ( ATL, /**/ ATR,  A00, A01, /**/ A02,
                   /**/       A10, A11, /**/ A12,
             /*************/ /******************/
              ABL, /**/ ABR,  A20, A21, /**/ A22 );

            RepartitionUp
            ( tT,  t0,
                   t1,
             /**/ /**/
              tB,  t2 );
            
            if( A00.Height() > 0 )
            {
                APan_MC_Star.AlignWith( A00 );
                APan_MR_Star.AlignWith( A00 );
                APan_MC_Star.ResizeTo( ATL.Height(), A11.Width() );
                APan_MR_Star.ResizeTo( ATL.Height(), A11.Width() );
                WPan.AlignWith( A01 );
                WPan_MC_Star.AlignWith( A00 );
                WPan_MR_Star.AlignWith( A00 );
                WPan.ResizeTo( ATL.Height(), A11.Width() );
                WPan_MC_Star.ResizeTo( ATL.Height(), A11.Width() );
                WPan_MR_Star.ResizeTo( ATL.Height(), A11.Width() );
                PartitionUp
                ( APan_MC_Star, A01_MC_Star,
                                A11_MC_Star, A11.Height() );
                PartitionUp
                ( APan_MR_Star, A01_MR_Star,
                                A11_MR_Star, A11.Height() );
                PartitionUp
                ( WPan_MC_Star, W01_MC_Star,
                                W11_MC_Star, A11.Height() );
                PartitionUp
                ( WPan_MR_Star, W01_MR_Star,
                                W11_MR_Star, A11.Height() );
                //------------------------------------------------------------//
                // Accumulate the Householder vectors into A01 and form W01 such
                // that subtracting (A01 W01' + W01 A01') is equal to 
                // successively applying the similarity transformations 
                // (I-conj(tau) h h')A00(I-tau h h') for each (tau,h).
                //
                // APan[MC,* ], APan[MR,* ], WPan[MC,* ], and WPan[MR,* ] are 
                // formed during the panel factorization.
                advanced::internal::PanelTridiagU
                ( ATL, WPan, t1,
                  APan_MC_Star, APan_MR_Star, WPan_MC_Star, WPan_MR_Star );
                basic::internal::LocalTriangularRank2K
                ( Upper, ConjugateTranspose, ConjugateTranspose,
                  (C)-1, A01_MC_Star, W01_MC_Star, A01_MR_Star, W01_MR_Star,
                  (C)1, A00 );
                //------------------------------------------------------------//
                APan_MC_Star.FreeAlignments();
                APan_MR_Star.FreeAlignments();
                WPan.FreeAlignments();
                WPan_MC_Star.FreeAlignments();
                WPan_MR_Star.FreeAlignments();
            }
            else
            {
                A11_Star_Star = A11;
                t1_Star_Star.ResizeTo( t1.Height(), 1 );

                advanced::Tridiag
                ( Upper, A11_Star_Star.LocalMatrix(), 
                  t1_Star_Star.LocalMatrix() );

                A11 = A11_Star_Star;
                t1 = t1_Star_Star;
            }

            SlidePartitionUpDiagonal
            ( ATL, /**/ ATR,  A00, /**/ A01, A02,
             /*************/ /******************/
                   /**/       A10, /**/ A11, A12,
              ABL, /**/ ABR,  A20, /**/ A21, A22 );

            SlidePartitionUp
            ( tT,  t0,
             /**/ /**/
                   t1,
              tB,  t2 );
        }
    }
    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::advanced::internal::TridiagU
( DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::TridiagU
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::internal::TridiagU
( DistMatrix<scomplex,MC,  MR  >& A, 
  DistMatrix<scomplex,Star,Star>& t );

template void elemental::advanced::internal::TridiagU
( DistMatrix<dcomplex,MC,  MR  >& A, 
  DistMatrix<dcomplex,Star,Star>& t );
#endif

