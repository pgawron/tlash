/*
   Copyright (c) 2009-2011, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin.
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

   Note: Written by Martin Schatz in March, 2011.
*/
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::basic::Strassen
( const DistMatrix<T,MC,MR>& A,
  const DistMatrix<T,MC,MR>& B,
  DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
  PushCallStck("basic::Strassen");
#endif
    if(A.Height() <= A.Grid().Height()){
        C.SetToZero();
        Orientation o = CharToOrientation('N');
        Gemm(o, o, T(1), A,B, T(1), C);
    }
    else{

        const Grid& g = A.Grid();
        const int n = A.Width() / 2;
        const int m = A.Height() / 2;
        T no = (T)-1;
        T o = (T)1;

        if( g.VCRank() == 0 )
        {
            std::cout << "WARNING: Strassen support is still just experimental." 
                      << std::endl;
        }

        DistMatrix<T,MC,MR> ATL(g), ATR(g), 
                            ABL(g), ABR(g);

        DistMatrix<T,MC,MR> BTL(g), BTR(g), 
                            BBL(g), BBR(g);

        DistMatrix<T,MC,MR> CTL(g), CTR(g), 
                            CBL(g), CBR(g);

        DistMatrix<T,MC,MR> M1(m, n, g), M2(m, n, g), M3(m, n, g), 
                            M4(m, n, g), M5(m, n, g), M6(m, n, g), 
                            M7(m, n, g);

        DistMatrix<T,MC,MR> tA(m, n, g), tB(m, n, g);

        ATL.LockedView(A, 0, 0, m, n);
        ATR.LockedView(A, 0, n, m, n);
        ABL.LockedView(A, m, 0, m, n);
        ABR.LockedView(A, m, n, m, n);

        BTL.LockedView(B, 0, 0, m, n);
        BTR.LockedView(B, 0, n, m, n);
        BBL.LockedView(B, m, 0, m, n);
        BBR.LockedView(B, m, n, m, n);

        CTL.View(C, 0, 0, m, n);
        CTR.View(C, 0, n, m, n);
        CBL.View(C, m, 0, m, n);
        CBR.View(C, m, n, m, n);

        //Calc M1
        tA = ABR;
        basic::Axpy(o, ATL, tA);
    
        tB = BBR;
        basic::Axpy(o, BTL, tB);
        
        Strassen(tA, tB, M1);

        //Calc M2
        tA = ABL;
        basic::Axpy(o, ABR, tA);

        tB = BTL;
        
        Strassen(tA, tB, M2);
        
        //Calc M3
        tA = ATL;

        tB = BTR;
        basic::Axpy(no, BBR, tB);
        
        Strassen(tA, tB, M3);
        
        //Calc M4
        tA = ABR;

        tB = BBL;
        basic::Axpy(no, BTL, tB);
        
        Strassen(tA, tB, M4);
        
        //Calc M5
        tA = ATL;
        basic::Axpy(o, ATR, tA);

        tB = BBR;
        
        Strassen(tA, tB, M5);
        
        //Calc M6
        tA = ABL;
        basic::Axpy(no, ATL, tA);

        tB = BTL;
        basic::Axpy(o, BTR, tB);
        
        Strassen(tA, tB, M6);
        
        //Calc M7
        tA = ATR;
        basic::Axpy(no, ABR, tA);

        tB = BBL;
        basic::Axpy(o, BBR, tB);
        
        Strassen(tA, tB, M7);

        //Calc CTL
        CTL = M7;
        basic::Axpy(no, M5, CTL);
        basic::Axpy(o, M4, CTL);
        basic::Axpy(o, M1, CTL);
        
        //Calc CTR
        CTR = M5;
        basic::Axpy(o, M3, CTR);
        
        //Calc CBL
        CBL = M4;
        basic::Axpy(o, M2, CBL);
        
        //Calc CBR
        CBR = M6;
        basic::Axpy(o, M3, CBR);
        basic::Axpy(no, M2, CBR);
        basic::Axpy(o, M1, CBR);
    }

#ifndef RELEASE
  PopCallStack();
#endif
}

template void elemental::basic::Strassen
( const DistMatrix<float,MC,MR>& A,
  const DistMatrix<float,MC,MR>& B,
  DistMatrix<float,MC,MR>& C );

template void elemental::basic::Strassen
( const DistMatrix<double,MC,MR>& A,
  const DistMatrix<double,MC,MR>& B,
  DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::Strassen
( const DistMatrix<scomplex,MC,MR>& A,
  const DistMatrix<scomplex,MC,MR>& B,
  DistMatrix<scomplex,MC,MR>& C );

template void elemental::basic::Strassen
( const DistMatrix<dcomplex,MC,MR>& A,
  const DistMatrix<dcomplex,MC,MR>& B,
  DistMatrix<dcomplex,MC,MR>& C );
#endif
