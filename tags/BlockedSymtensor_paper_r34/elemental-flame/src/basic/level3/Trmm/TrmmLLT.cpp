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

#include "./TrmmUtil.hpp"

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

// Left Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := tril(L)^T,
//   X := tril(L)^H,
//   X := trilu(L)^T, or
//   X := trilu(L)^H
template<typename T>
void
elemental::basic::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& L,
        DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrmmLLT");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        basic::internal::TrmmLLTA( orientation, diagonal, alpha, L, X );
    else
        basic::internal::TrmmLLTC( orientation, diagonal, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}
 
template<typename T>
void
elemental::basic::internal::TrmmLLTA
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& L,
        DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrmmLLTA");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrmmLLTA expects a (Conjugate)Transpose option." );
    if( L.Height() != L.Width() || L.Height() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLLTA: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR>
        XL(g), XR(g),
        X0(g), X1(g), X2(g);

    DistMatrix<T,MC,Star> X1_MC_Star(g);
    DistMatrix<T,MR,Star> Z1_MR_Star(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    PartitionRight
    ( X, XL, XR, 0 );
    while( XL.Width() < X.Width() )
    {
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        X1_MC_Star.AlignWith( L );
        Z1_MR_Star.AlignWith( L );
        Z1_MR_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_MC_Star = X1;
        Z1_MR_Star.SetToZero();
        basic::internal::LocalTrmmAccumulateLLT
        ( orientation, diagonal, alpha, L, X1_MC_Star, Z1_MR_Star );

        Z1_MR_MC.SumScatterFrom( Z1_MR_Star );
        X1 = Z1_MR_MC;
        //--------------------------------------------------------------------//
        X1_MC_Star.FreeAlignments();
        Z1_MR_Star.FreeAlignments();

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
   
template<typename T>
void
elemental::basic::internal::TrmmLLTC
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& L,
        DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrmmLLTC");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrmmLLT expects a (Conjugate)Transpose option." );
    if( L.Height() != L.Width() || L.Height() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLLTC: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR> XT(g),  X0(g),
                        XB(g),  X1(g),
                                X2(g);

    // Temporary distributions
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,MC,  Star> L21_MC_Star(g);
    DistMatrix<T,Star,VR  > X1_Star_VR(g);
    DistMatrix<T,Star,MR  > D1_Star_MR(g);

    // Start the algorithm
    basic::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 ); 

        L21_MC_Star.AlignWith( X2 );
        D1_Star_MR.AlignWith( X1 );
        D1_Star_MR.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_Star_VR = X1;
        L11_Star_Star = L11;
        basic::internal::LocalTrmm
        ( Left, Lower, orientation, diagonal, (T)1, L11_Star_Star, X1_Star_VR );
        X1 = X1_Star_VR;
 
        L21_MC_Star = L21;
        basic::internal::LocalGemm
        ( orientation, Normal, (T)1, L21_MC_Star, X2, (T)0, D1_Star_MR );
        X1.SumScatterUpdate( (T)1, D1_Star_MR );
       //--------------------------------------------------------------------//
        L21_MC_Star.FreeAlignments();
        D1_Star_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::LocalTrmmAccumulateLLT
( Orientation orientation, Diagonal diagonal, T alpha,
  const DistMatrix<T,MC,MR  >& L,
  const DistMatrix<T,MC,Star>& X_MC_Star,
        DistMatrix<T,MR,Star>& Z_MR_Star )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrmmAccumulateLLT");
    if( L.Grid() != X_MC_Star.Grid() ||
        X_MC_Star.Grid() != Z_MR_Star.Grid() )
        throw logic_error( "{L,X,Z} must be distributed over the same grid." );
    if( L.Height() != L.Width() ||
        L.Height() != X_MC_Star.Height() ||
        L.Height() != Z_MR_Star.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateLLT: " << "\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X[MC,* ] ~ " << X_MC_Star.Height() << " x "
                               << X_MC_Star.Width() << "\n"
            << "  Z[MR,* ] ` " << Z_MR_Star.Height() << " x "
                               << Z_MR_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( X_MC_Star.ColAlignment() != L.ColAlignment() ||
        Z_MR_Star.ColAlignment() != L.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = L.Grid();
    
    // Matrix views
    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,MC,Star>
        XT_MC_Star(g),  X0_MC_Star(g),
        XB_MC_Star(g),  X1_MC_Star(g),
                        X2_MC_Star(g);

    DistMatrix<T,MR,Star>
        ZT_MR_Star(g),  Z0_MR_Star(g),
        ZB_MR_Star(g),  Z1_MR_Star(g),
                        Z2_MR_Star(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionDown
    ( X_MC_Star, XT_MC_Star,
                 XB_MC_Star, 0 );
    PartitionDown
    ( Z_MR_Star, ZT_MR_Star,
                 ZB_MR_Star, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionDown
        ( XT_MC_Star,  X0_MC_Star,
         /**********/ /**********/
                       X1_MC_Star,
          XB_MC_Star,  X2_MC_Star );

        RepartitionDown
        ( ZT_MR_Star,  Z0_MR_Star,
         /**********/ /**********/
                       Z1_MR_Star,
          ZB_MR_Star,  Z2_MR_Star );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        D11.MakeTrapezoidal( Left, Lower );
        if( diagonal == Unit )
            SetDiagonalToOne( D11 );
        basic::internal::LocalGemm
        ( orientation, Normal,
          alpha, D11, X1_MC_Star, (T)1, Z1_MR_Star );

        basic::internal::LocalGemm
        ( orientation, Normal,
          alpha, L21, X2_MC_Star, (T)1, Z1_MR_Star );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlideLockedPartitionDown
        ( XT_MC_Star,  X0_MC_Star,
                       X1_MC_Star,
         /**********/ /**********/
          XB_MC_Star,  X2_MC_Star );

        SlidePartitionDown
        ( ZT_MR_Star,  Z0_MR_Star,
                       Z1_MR_Star,
         /**********/ /**********/
          ZB_MR_Star,  Z2_MR_Star );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void elemental::basic::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::basic::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

