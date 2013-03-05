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
elemental::basic::internal::TrsvLN
( Diagonal diagonal, 
  const DistMatrix<F,MC,MR>& L, 
        DistMatrix<F,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrsvLN");
    if( L.Grid() != x.Grid() )
        throw logic_error( "L and x must be distributed over the same grid." );
    if( L.Height() != L.Width() )
        throw logic_error( "L must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw logic_error( "x must be a vector." );
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( L.Width() != xLength )
        throw logic_error( "Nonconformal TrsvLN." );
#endif
    const Grid& g = L.Grid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<F,MC,MR> 
            LTL(g), LTR(g),  L00(g), L01(g), L02(g),
            LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                             L20(g), L21(g), L22(g);

        DistMatrix<F,MC,MR> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<F,Star,Star> L11_Star_Star(g);
        DistMatrix<F,Star,Star> x1_Star_Star(g);
        DistMatrix<F,MR,  Star> x1_MR_Star(g);
        DistMatrix<F,MC,  Star> z2_MC_Star(g);

        // Start the algorithm
        LockedPartitionDownDiagonal
        ( L, LTL, LTR,
             LBL, LBR, 0 );
        PartitionDown
        ( x, xT,
             xB, 0 );
        while( xB.Height() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionDown
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );

            x1_MR_Star.AlignWith( L21 );
            z2_MC_Star.AlignWith( L21 );
            z2_MC_Star.ResizeTo( x2.Height(), 1 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            basic::Trsv
            ( Lower, Normal, diagonal,
              L11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_MR_Star = x1_Star_Star;
            basic::Gemv
            ( Normal, (F)-1, 
              L21.LockedLocalMatrix(), 
              x1_MR_Star.LockedLocalMatrix(),
              (F)0, z2_MC_Star.LocalMatrix() );
            x2.SumScatterUpdate( (F)1, z2_MC_Star );
            //----------------------------------------------------------------//
            x1_MR_Star.FreeAlignments();
            z2_MC_Star.FreeAlignments();

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            SlidePartitionDown
            ( xT,  x0,
                   x1,
             /**/ /**/
              xB,  x2 );
        }
    }
    else
    {
        // Matrix views 
        DistMatrix<F,MC,MR> 
            LTL(g), LTR(g),  L00(g), L01(g), L02(g),
            LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                             L20(g), L21(g), L22(g);

        DistMatrix<F,MC,MR> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<F,Star,Star> L11_Star_Star(g);
        DistMatrix<F,Star,Star> x1_Star_Star(g);
        DistMatrix<F,Star,MR  > x1_Star_MR(g);
        DistMatrix<F,Star,MC  > z2_Star_MC(g);
        DistMatrix<F,MR,  MC  > z2_MR_MC(g);
        DistMatrix<F,MC,  MR  > z2(g);

        // Start the algorithm
        LockedPartitionDownDiagonal
        ( L, LTL, LTR,
             LBL, LBR, 0 );
        PartitionRight( x,  xL, xR, 0 );
        while( xR.Width() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionRight
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );

            x1_Star_MR.AlignWith( L21 );
            z2_Star_MC.AlignWith( L21 );
            z2.AlignWith( x2 );
            z2_Star_MC.ResizeTo( 1, x2.Width() );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            basic::Trsv
            ( Lower, Normal, diagonal,
              L11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_Star_MR = x1_Star_Star;
            basic::Gemv
            ( Normal, (F)-1, 
              L21.LockedLocalMatrix(), 
              x1_Star_MR.LockedLocalMatrix(),
              (F)0, z2_Star_MC.LocalMatrix() );
            z2_MR_MC.SumScatterFrom( z2_Star_MC );
            z2 = z2_MR_MC;
            basic::Axpy( (F)1, z2, x2 );
            //----------------------------------------------------------------//
            x1_Star_MR.FreeAlignments();
            z2_Star_MC.FreeAlignments();
            z2.FreeAlignments(); 

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            SlidePartitionRight
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& x );

template void elemental::basic::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& x );

template void elemental::basic::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

