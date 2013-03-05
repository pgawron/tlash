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

/*
   Parallelization of Variant 2 Upper Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A01^H A01
   A11 := Chol(A11)
   A12 := A12 - A01^H A02
   A12 := triu(A11)^-H A12
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A01[MC,* ] <- A01[MC,MR]
   X11^H[MR,* ] := (A01[MC,MR])^H A01[MC,* ]
   A11[MC,MR] := A11[MC,MR] - ((SumCol(X11^H[MR,* ]))[* ,MC])^H

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X12^H[MR,* ] := (A02[MC,MR])^H A01[MC,* ]
   A12[MC,MR] := A12[MC,MR] - ((SumCol(X12^H[MR,* ]))[MC,* ])^H

   A12[* ,VR] <- A12[MC,MR]
   A12[* ,VR] := triu(A11[* ,* ])^-H A12[* ,VR]
   A12[MC,MR] <- A12[* ,VR]
   -----------------------------------------------------
*/
template<typename F> // representation of real or complex number
void
elemental::advanced::internal::CholUVar2
( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CholUVar2");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,MC,  Star> A01_MC_Star(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,Star,VR  > A12_Star_VR(g);
    DistMatrix<F,MR,  Star> X11Herm_MR_Star(g);
    DistMatrix<F,MR,  MC  > X11Herm_MR_MC(g);
    DistMatrix<F,MC,  MR  > X11(g);
    DistMatrix<F,MR,  Star> X12Herm_MR_Star(g);
    DistMatrix<F,MR,  MC  > X12Herm_MR_MC(g);
    DistMatrix<F,MC,  MR  > X12(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A01_MC_Star.AlignWith( A01 );
        X11Herm_MR_Star.AlignWith( A01 );
        X11Herm_MR_MC.AlignWith( A11 );
        X11.AlignWith( A11 );
        X12Herm_MR_Star.AlignWith( A02 );
        X12Herm_MR_MC.AlignWith( A12 );
        X12.AlignWith( A12 );
        X11Herm_MR_Star.ResizeTo( A11.Width(), A11.Height() );
        X12Herm_MR_Star.ResizeTo( A12.Width(), A12.Height() );
        //--------------------------------------------------------------------//
        A01_MC_Star = A01;
        basic::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (F)1, A01, A01_MC_Star, (F)0, X11Herm_MR_Star );
        X11Herm_MR_MC.SumScatterFrom( X11Herm_MR_Star );
        basic::ConjTrans( X11Herm_MR_MC, X11 );
        basic::Axpy( (F)-1, X11, A11 );

        A11_Star_Star = A11;
        advanced::internal::LocalChol( Upper, A11_Star_Star );
        A11 = A11_Star_Star;

        basic::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (F)1, A02, A01_MC_Star, (F)0, X12Herm_MR_Star );
        X12Herm_MR_MC.SumScatterFrom( X12Herm_MR_Star );
        basic::ConjTrans( X12Herm_MR_MC, X12 );
        basic::Axpy( (F)-1, X12, A12 );

        A12_Star_VR = A12;
        basic::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (F)1, A11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_MC_Star.FreeAlignments();
        X11Herm_MR_Star.FreeAlignments();
        X11Herm_MR_MC.FreeAlignments();
        X11.FreeAlignments();
        X12Herm_MR_Star.FreeAlignments();
        X12Herm_MR_MC.FreeAlignments();
        X12.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

/*
   Naive parallelization of Variant 2 Upper Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A01^H A01
   A11 := Chol(A11)
   A12 := A12 - A01^H A02
   A12 := triu(A11)^-H A12
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A01[MC,* ] <- A01[MC,MR]
   X11[* ,MR] := (A01[MC,* ])^H A01[MC,MR]
   A11[MC,MR] := A11[MC,MR] - (SumCol(X11[* ,MR]))[MC,* ]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X12[* ,MR] := (A01[MC,* ])^H A02[MC,MR]
   A12[MC,MR] := A12[MC,MR] - (SumCol(X12[* ,MR]))[MC,* ]

   A12[* ,VR] <- A12[MC,MR]
   A12[* ,VR] := triu(A11[* ,* ])^-H A12[* ,VR]
   A12[MC,MR] <- A12[* ,VR]
   -----------------------------------------------------
*/
template<typename F> // representation of real or complex number
void
elemental::advanced::internal::CholUVar2Naive
( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CholUVar2Naive");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
    if( A.Grid().VCRank() == 0 )
    {
        cout << "CholUVar2Naive exists solely for academic purposes. Please "
                "use CholUVar2 in real applications." << endl;
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,MC,  Star> A01_MC_Star(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,Star,VR  > A12_Star_VR(g);
    DistMatrix<F,Star,MR  > X11_Star_MR(g);
    DistMatrix<F,Star,MR  > X12_Star_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A01_MC_Star.AlignWith( A01 );
        X11_Star_MR.AlignWith( A01 );
        X12_Star_MR.AlignWith( A02 );
        X11_Star_MR.ResizeTo( A11.Height(), A11.Width() );
        X12_Star_MR.ResizeTo( A12.Height(), A12.Width() );
        //--------------------------------------------------------------------//
        A01_MC_Star = A01;
        basic::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (F)1, A01_MC_Star, A01, (F)0, X11_Star_MR );
        A11.SumScatterUpdate( (F)-1, X11_Star_MR );

        A11_Star_Star = A11;
        advanced::internal::LocalChol( Upper, A11_Star_Star );
        A11 = A11_Star_Star;

        basic::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (F)1, A01_MC_Star, A02, (F)0, X12_Star_MR );
        A12.SumScatterUpdate( (F)-1, X12_Star_MR );

        A12_Star_VR = A12;
        basic::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (F)1, A11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_MC_Star.FreeAlignments();
        X11_Star_MR.FreeAlignments();
        X12_Star_MR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::advanced::internal::CholUVar2
( DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::CholUVar2Naive
( DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::CholUVar2
( DistMatrix<double,MC,MR>& A );

template void elemental::advanced::internal::CholUVar2Naive
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::internal::CholUVar2
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::advanced::internal::CholUVar2Naive
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::advanced::internal::CholUVar2
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::advanced::internal::CholUVar2Naive
( DistMatrix<dcomplex,MC,MR>& A );
#endif

