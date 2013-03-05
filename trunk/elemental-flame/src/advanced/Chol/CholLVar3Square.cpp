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
using namespace elemental::imports;

/*
   Parallelization of Variant 3 Lower Cholesky factorization for
   square process grids.

   Original serial update:
   ------------------------
   A11 := Chol(A11) 
   A21 := A21 tril(A11)^-H
   A22 := A22 - A21 A21^H
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   
   A21^T[* ,MC] <- A21[VC,* ]
   A21^H[* ,MR] <- A21^T[* ,MC]
   A22[MC,MR] := A22[MC,MR] - (A21^T[* ,MC])^T A21^H[* ,MR]
   A21[MC,MR] <- A21^T[* ,MC]
   -----------------------------------------------------
*/
template<typename F> // representation of real or complex number
void
elemental::advanced::internal::CholLVar3Square
( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CholLVar3Square");
    if( A.Height() != A.Width() )
        throw logic_error
        ("Can only compute Cholesky factor of square matrices.");
    if( A.Grid().Height() != A.Grid().Width() )
        throw logic_error
        ("CholLVar3Square requires a square process grid.");
#endif
    const Grid& g = A.Grid();

    // Find the process holding our transposed data
    const int r = g.Height();
    int transposeRank;
    {
        int colAlignment = A.ColAlignment();
        int rowAlignment = A.RowAlignment();
        int colShift = A.ColShift();
        int rowShift = A.RowShift();

        int transposeRow = (colAlignment+rowShift) % r;
        int transposeCol = (rowAlignment+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    bool onDiagonal = ( transposeRank == g.VCRank() );

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrices
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,VC,  Star> A21_VC_Star(g);
    DistMatrix<F,Star,MC  > A21Trans_Star_MC(g);
    DistMatrix<F,Star,MR  > A21Herm_Star_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/   
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A21_VC_Star.AlignWith( A22 );
        A21Trans_Star_MC.AlignWith( A22 );
        A21Herm_Star_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        advanced::internal::LocalChol( Lower, A11_Star_Star );
        A11 = A11_Star_Star;

        A21_VC_Star = A21;
        basic::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (F)1, A11_Star_Star, A21_VC_Star );

        A21Trans_Star_MC.TransposeFrom( A21_VC_Star );
        // SendRecv to from A21^T[* ,MR] from A21^T[* ,MC], then conjugate
        // the buffer to from A21^H[* ,MR]
        A21Herm_Star_MR.ResizeTo( A21.Width(), A21.Height() ); 
        {
            if( onDiagonal )
            { 
                int size = A11.Height()*A22.LocalWidth();
                memcpy
                ( A21Herm_Star_MR.LocalBuffer(), A21Trans_Star_MC.LocalBuffer(),
                  size*sizeof(F) );
            }
            else
            {
                int sendSize = A22.LocalHeight()*A11.Width();
                int recvSize = A22.LocalWidth()*A11.Height();
                // We know that the ldim is the height since we have manually 
                // created both temporary matrices.
                mpi::SendRecv 
                ( A21Trans_Star_MC.LocalBuffer(), sendSize, transposeRank, 0,
                  A21Herm_Star_MR.LocalBuffer(),  recvSize, transposeRank, 0,
                  g.VCComm() );
            }
            basic::Conj( A21Herm_Star_MR );
        }

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        basic::internal::LocalTriangularRankK
        ( Lower, Transpose, 
          (F)-1, A21Trans_Star_MC, A21Herm_Star_MR, (F)1, A22 );

        A21.TransposeFrom( A21Trans_Star_MC );
        //--------------------------------------------------------------------//
        A21_VC_Star.FreeAlignments();
        A21Trans_Star_MC.FreeAlignments();
        A21Herm_Star_MR.FreeAlignments();

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

template void elemental::advanced::internal::CholLVar3Square
( DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::CholLVar3Square
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::internal::CholLVar3Square
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::advanced::internal::CholLVar3Square
( DistMatrix<dcomplex,MC,MR>& A );
#endif

