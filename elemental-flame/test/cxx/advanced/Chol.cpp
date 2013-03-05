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
#include <ctime>
#include "elemental.hpp"
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

void Usage()
{
    cout << "Generates SPD matrix then solves for its Cholesky factor.\n\n"
         << "  Chol <r> <c> <shape> <m> <nb> <rankK local nb> <correctness?> "
            "<print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  shape: {L,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  rankK local nb: local blocksize for triangular rank-k update\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename F> // represents a real or complex field
void TestCorrectness
( bool printMatrices, Shape shape,
  const DistMatrix<F,MC,MR>& A,
  const DistMatrix<F,MC,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F,MC,MR> X(m,100,g);
    DistMatrix<F,MC,MR> Y(m,100,g);
    X.SetToRandom();
    Y = X;

    if( shape == Lower )
    {
        // Test correctness by comparing the application of AOrig against a 
        // random set of 100 vectors to the application of tril(A) tril(A)^H
        basic::Trmm( Left, Lower, ConjugateTranspose, NonUnit, (F)1, A, Y );
        basic::Trmm( Left, Lower, Normal, NonUnit, (F)1, A, Y );
        basic::Hemm( Left, Lower, (F)-1, AOrig, X, (F)1, Y );
        F oneNormOfError = advanced::OneNorm( Y );
        F infNormOfError = advanced::InfinityNorm( Y );
        F frobNormOfError = advanced::FrobeniusNorm( Y );
        F infNormOfA = advanced::HermitianInfinityNorm( shape, AOrig );
        F frobNormOfA = advanced::HermitianFrobeniusNorm( shape, AOrig );
        F oneNormOfX = advanced::OneNorm( X );
        F infNormOfX = advanced::InfinityNorm( X );
        F frobNormOfX = advanced::FrobeniusNorm( X );
        if( g.VCRank() == 0 )
        {
            cout << "||A||_1 = ||A||_oo   = " << Abs(infNormOfA) << "\n"
                 << "||A||_F              = " << Abs(frobNormOfA) << "\n"
                 << "||X||_1              = " << Abs(oneNormOfX) << "\n"
                 << "||X||_oo             = " << Abs(infNormOfX) << "\n"
                 << "||X||_F              = " << Abs(frobNormOfX) << "\n"
                 << "||A X - L L^H X||_1  = " << Abs(oneNormOfError) << "\n"
                 << "||A X - L L^H X||_oo = " << Abs(infNormOfError) << "\n"
                 << "||A X - L L^H X||_F  = " << Abs(frobNormOfError) << endl;
        }
    }
    else
    {
        // Test correctness by comparing the application of AOrig against a 
        // random set of 100 vectors to the application of triu(A)^H triu(A)
        basic::Trmm( Left, Upper, Normal, NonUnit, (F)1, A, Y );
        basic::Trmm( Left, Upper, ConjugateTranspose, NonUnit, (F)1, A, Y );
        basic::Hemm( Left, Upper, (F)-1, AOrig, X, (F)1, Y );
        F oneNormOfError = advanced::OneNorm( Y );
        F infNormOfError = advanced::InfinityNorm( Y );
        F frobNormOfError = advanced::FrobeniusNorm( Y );
        F infNormOfA = advanced::HermitianInfinityNorm( shape, AOrig );
        F frobNormOfA = advanced::HermitianFrobeniusNorm( shape, AOrig );
        F oneNormOfX = advanced::OneNorm( X );
        F infNormOfX = advanced::InfinityNorm( X );
        F frobNormOfX = advanced::FrobeniusNorm( X );
        if( g.VCRank() == 0 )
        {
            cout << "||A||_1 = ||A||_oo   = " << Abs(infNormOfA) << "\n"
                 << "||A||_F              = " << Abs(frobNormOfA) << "\n"
                 << "||X||_1              = " << Abs(oneNormOfX) << "\n"
                 << "||X||_oo             = " << Abs(infNormOfX) << "\n"
                 << "||X||_F              = " << Abs(frobNormOfX) << "\n"
                 << "||A X - U^H U X||_1  = " << Abs(oneNormOfError) << "\n"
                 << "||A X - U^H U X||_oo = " << Abs(infNormOfError) << "\n"
                 << "||A X - U^H U X||_F  = " << Abs(frobNormOfError) << endl;
        }
    }
}

template<typename F> // represents a real or complex field
void TestChol
( bool testCorrectness, bool printMatrices, 
  Shape shape, int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<F,MC,MR> A(g);
    DistMatrix<F,MC,MR> AOrig(g);

    A.ResizeTo( m, m );

    A.SetToRandomHPD();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting Cholesky factorization...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    advanced::Chol( shape, A );
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = advanced::internal::CholGFlops<F>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE.\n"
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, AOrig );
}

int 
main( int argc, char* argv[] )
{
    Init( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    int rank = mpi::CommRank( comm );

    if( argc < 9 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try
    {
        int argNum = 0;
        const int r = atoi(argv[++argNum]);
        const int c = atoi(argv[++argNum]);
        const Shape shape = CharToShape(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const int nbLocal = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        const Grid g( comm, r, c );
        SetBlocksize( nb );
        basic::SetLocalTriangularRankKBlocksize<double>( nbLocal );
#ifndef WITHOUT_COMPLEX
        basic::SetLocalTriangularRankKBlocksize< std::complex<double> >
        ( nbLocal );
#endif

        if( rank == 0 )
            cout << "Will test Chol" << ShapeToChar(shape) << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestChol<double>
        ( testCorrectness, printMatrices, shape, m, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestChol<dcomplex>
        ( testCorrectness, printMatrices, shape, m, g );
#endif
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:\n"
             << e.what() << endl;
    }   

    Finalize();
    return 0;
}

