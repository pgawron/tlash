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
elemental::basic::internal::TrsvUN
( Diagonal diagonal, 
  const DistMatrix<F,MC,MR>& U, 
        DistMatrix<F,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrsvUN");
    if( U.Grid() != x.Grid() )
        throw logic_error( "U and x must be distributed over the same grid." );
    if( U.Height() != U.Width() )
        throw logic_error( "U must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw logic_error( "x must be a vector." );
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( U.Width() != xLength )
        throw logic_error( "Nonconformal TrsvUN." );
#endif
    const Grid& g = U.Grid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<F,MC,MR> 
            UTL(g), UTR(g),  U00(g), U01(g), U02(g),
            UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                             U20(g), U21(g), U22(g);

        DistMatrix<F,MC,MR> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<F,Star,Star> U11_Star_Star(g);
        DistMatrix<F,Star,Star> x1_Star_Star(g);
        DistMatrix<F,MR,  Star> x1_MR_Star(g);
        DistMatrix<F,MC,  Star> z0_MC_Star(g);

        // Start the algorithm
        LockedPartitionUpDiagonal
        ( U, UTL, UTR,
             UBL, UBR, 0 );
        PartitionUp
        ( x, xT,
             xB, 0 );
        while( xT.Height() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, U01, /**/ U02,
                   /**/       U10, U11, /**/ U12,
             /*************/ /******************/
              UBL, /**/ UBR,  U20, U21, /**/ U22 );

            RepartitionUp
            ( xT,  x0,
                   x1,
             /**/ /**/
              xB,  x2 );

            x1_MR_Star.AlignWith( U01 );
            z0_MC_Star.AlignWith( U01 );
            z0_MC_Star.ResizeTo( x0.Height(), 1 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            U11_Star_Star = U11;
            basic::Trsv
            ( Upper, Normal, diagonal,
              U11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_MR_Star = x1_Star_Star;
            basic::Gemv
            ( Normal, (F)-1, 
              U01.LockedLocalMatrix(), 
              x1_MR_Star.LockedLocalMatrix(),
              (F)0, z0_MC_Star.LocalMatrix() );
            x0.SumScatterUpdate( (F)1, z0_MC_Star );
            //----------------------------------------------------------------//
            x1_MR_Star.FreeAlignments();
            z0_MC_Star.FreeAlignments();

            SlideLockedPartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, /**/ U01, U02,
             /*************/ /******************/
                   /**/       U10, /**/ U11, U12,
              UBL, /**/ UBR,  U20, /**/ U21, U22 );

            SlidePartitionUp
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );
        }
    }
    else
    {
        // Matrix views 
        DistMatrix<F,MC,MR> 
            UTL(g), UTR(g),  U00(g), U01(g), U02(g),
            UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                             U20(g), U21(g), U22(g);

        DistMatrix<F,MC,MR> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<F,Star,Star> U11_Star_Star(g);
        DistMatrix<F,Star,Star> x1_Star_Star(g);
        DistMatrix<F,Star,MR  > x1_Star_MR(g);
        DistMatrix<F,Star,MC  > z0_Star_MC(g);
        DistMatrix<F,MR,  MC  > z0_MR_MC(g);
        DistMatrix<F,MC,  MR  > z0(g);

        // Start the algorithm
        LockedPartitionUpDiagonal
        ( U, UTL, UTR,
             UBL, UBR, 0 );
        PartitionLeft( x,  xL, xR, 0 );
        while( xL.Width() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, U01, /**/ U02,
                   /**/       U10, U11, /**/ U12,
             /*************/ /******************/
              UBL, /**/ UBR,  U20, U21, /**/ U22 );

            RepartitionLeft
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );

            x1_Star_MR.AlignWith( U01 );
            z0_Star_MC.AlignWith( U01 );
            z0.AlignWith( x0 );
            z0_Star_MC.ResizeTo( 1, x0.Width() );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            U11_Star_Star = U11;
            basic::Trsv
            ( Upper, Normal, diagonal,
              U11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_Star_MR = x1_Star_Star;
            basic::Gemv
            ( Normal, (F)-1, 
              U01.LockedLocalMatrix(), 
              x1_Star_MR.LockedLocalMatrix(),
              (F)0, z0_Star_MC.LocalMatrix() );
            z0_MR_MC.SumScatterFrom( z0_Star_MC );
            z0 = z0_MR_MC;
            basic::Axpy( (F)1, z0, x0 );
            //----------------------------------------------------------------//
            x1_Star_MR.FreeAlignments();
            z0_Star_MC.FreeAlignments();
            z0.FreeAlignments(); 

            SlideLockedPartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, /**/ U01, U02,
             /*************/ /******************/
                   /**/       U10, /**/ U11, U12,
              UBL, /**/ UBR,  U20, /**/ U21, U22 );

            SlidePartitionLeft
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& x );

template void elemental::basic::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& x );

template void elemental::basic::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

