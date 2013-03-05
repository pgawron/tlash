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
    cout << "Tests UT transform application.\n\n"
         << "  UT <r> <c> <side> <shape> <orientation> <m> <offset> <nb> "
         << "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  side: {L/R}\n"
         << "  shape: {L/U}\n"
         << "  orientation: {N/C}\n"
         << "  m: height of matrix\n"
         << "  offset: diagonal transforms are stored above/below\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename R> // represents a real number
void TestCorrectness
( Side side, 
  Shape shape,
  Orientation orientation,
  int offset,
  bool printMatrices,
  const DistMatrix<R,MC,MR>& H )
{
    const Grid& g = H.Grid();
    const int m = H.Height();

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of transform..." << endl;

    // Form Z := Q^H Q or Q^H Q as an approximation to identity
    DistMatrix<R,MC,MR> Y(m,m,g);
    Y.SetToIdentity();
    advanced::UT( side, shape, orientation, offset, H, Y );
    if( printMatrices )
    {
        DistMatrix<R,MC,MR> W(m,m,g);
        W.SetToIdentity();
        if( orientation == Normal )
        {
            advanced::UT( side, shape, ConjugateTranspose, offset, H, W );
            Y.Print("Q");
            W.Print("Q^H");
        }
        else
        {
            advanced::UT( side, shape, Normal, offset, H, W );
            Y.Print("Q^H");
            W.Print("Q");
        }
    }
    DistMatrix<R,MC,MR> Z(m,m,g);
    Z.SetToZero();
    basic::Syrk( shape, Normal, 1.0, Y, 0.0, Z );

    // Form X := I - Q^H Q or Q Q^H
    DistMatrix<R,MC,MR> X(m,m,g);
    X.SetToIdentity();
    basic::Axpy( (R)-1, Z, X );
    if( printMatrices )
    {
        if( orientation == Normal )
            X.Print("I - Q Q^H");
        else
            X.Print("I - Q^H Q");
    }

    R oneNormOfError = advanced::OneNorm( X );
    R infNormOfError = advanced::InfinityNorm( X );
    R frobNormOfError = advanced::FrobeniusNorm( X );
    if( g.VCRank() == 0 )
    {
        if( orientation == Normal )
        {
            cout << "    ||Q Q^H - I||_1  = " << oneNormOfError << "\n"
                 << "    ||Q Q^H - I||_oo = " << infNormOfError << "\n"
                 << "    ||Q Q^H - I||_F  = " << frobNormOfError << endl;
        }
        else
        {
            cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
                 << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
                 << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
        }
    }
}

#ifndef WITHOUT_COMPLEX
template<typename R> // represents a real number
void TestCorrectness
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  bool printMatrices,
  const DistMatrix<complex<R>,MC,MR  >& H,
  const DistMatrix<complex<R>,MD,Star>& t )
{
    typedef complex<R> C;

    const Grid& g = H.Grid();
    const int m = H.Height();

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of transform..." << endl;

    // Form Z := Q^H Q or Q Q^H as an approximation to identity
    DistMatrix<C,MC,MR> Y(m,m,g);
    Y.SetToIdentity();
    advanced::UT( side, shape, orientation, offset, H, t, Y );
    if( printMatrices )
    {
        DistMatrix<C,MC,MR> W(m,m,g);
        W.SetToIdentity();
        if( orientation == Normal )
        {
            advanced::UT( side, shape, ConjugateTranspose, offset, H, t, W );
            Y.Print("Q");
            W.Print("Q^H");
        }
        else
        {
            advanced::UT( side, shape, Normal, offset, H, t, W );
            Y.Print("Q^H");
            W.Print("Q");
        }
    }
    DistMatrix<C,MC,MR> Z(m,m,g);
    Z.SetToZero();
    basic::Herk( shape, Normal, (C)1, Y, (C)0, Z );
    
    // Form X := I - Q^H Q or Q Q^H
    DistMatrix<C,MC,MR> X(m,m,g);
    X.SetToIdentity();
    basic::Axpy( (C)-1, Z, X );
    if( printMatrices )
    {
        if( orientation == Normal )
            X.Print("I - Q Q^H");
        else
            X.Print("I - Q^H Q");
    }

    // Compute the maximum deviance
    R oneNormOfError = advanced::OneNorm( X );
    R infNormOfError = advanced::InfinityNorm( X );
    R frobNormOfError = advanced::FrobeniusNorm( X );
    if( g.VCRank() == 0 )
    {
        if( orientation == Normal )
        {
            cout << "    ||Q Q^H - I||_1  = " << oneNormOfError << "\n"
                 << "    ||Q Q^H - I||_oo = " << infNormOfError << "\n"
                 << "    ||Q Q^H - I||_F  = " << frobNormOfError << endl;
        }
        else
        {
            cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
                 << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
                 << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
        }
    }
}
#endif // WITHOUT_COMPLEX

template<typename F> // represents a real or complex number
void TestUT
( Side side, Shape shape, Orientation orientation, 
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g );

template<>
void TestUT<double>
( Side side, Shape shape, Orientation orientation,
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> H(g);
    DistMatrix<R,MC,MR> A(g);

    H.ResizeTo( m, m );
    A.ResizeTo( m, m );

    H.SetToRandom();
    A.SetToRandom();
    if( printMatrices )
    {
        H.Print("H");
        A.Print("A");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting UT transform...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    advanced::UT( side, shape, orientation, offset, H, A );
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = advanced::internal::UTGFlops<R>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( side, shape, orientation, offset, printMatrices, H );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestUT< complex<double> >
( Side side, Shape shape, Orientation orientation,
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g )
{
    typedef complex<double> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,MR  > H(g);
    DistMatrix<C,MD,Star> t(g);
    DistMatrix<C,MC,MR  > A(g);

    H.ResizeTo( m, m );
    A.ResizeTo( m, m );

    t.AlignWithDiag( H, offset );
    t.ResizeTo( H.DiagonalLength( offset ), 1 );

    H.SetToRandom();
    A.SetToRandom();
    DistMatrix<C,MC,MR> HCol(g);
    if( shape == Lower )
    {
        for( int i=0; i<t.Height(); ++i )
        {
            // View below the diagonal containing the implicit 1
            HCol.View( H, i-offset+1, i, m-(i-offset+1), 1 );
            C norm = basic::Nrm2( HCol );
            C alpha = 2./(norm*norm+1.);
            t.Set( i, 0, alpha );
        }
    }
    else
    {
        for( int i=0; i<t.Height(); ++i ) 
        {
            // View above the diagonal containing the implicit 1
            HCol.View( H, 0, i+offset, i, 1 );
            C norm = basic::Nrm2( HCol );
            C alpha = 2./(norm*norm+1.);
            t.Set( i, 0, alpha );
        }
    }

    if( printMatrices )
    {
        H.Print("H");
        A.Print("A");
        t.Print("t");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting UT transform...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    advanced::UT( side, shape, orientation, offset, H, t, A );
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = advanced::internal::UTGFlops<C>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
    {
        TestCorrectness
        ( side, shape, orientation, offset, printMatrices, H, t );
    }
}
#endif

int 
main( int argc, char* argv[] )
{
    Init( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 11 )
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
        const Side side = CharToSide(*argv[++argNum]);
        const Shape shape = CharToShape(*argv[++argNum]);
        const Orientation orientation = CharToOrientation(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int offset = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);
        if( shape == Lower && offset > 0 )
        {
            throw std::runtime_error
                  ("The offset cannot be positive if the transforms are in "
                   "the lower triangle.");
        }
        else if( shape == Upper && offset < 0 )
        {
            throw std::runtime_error
                  ("The offset cannot be negative if the transforms are in "
                   "the upper triangle.");
        }
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

        if( rank == 0 )
            cout << "Will test UT transform" << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestUT<double>
        ( side, shape, orientation, m, offset, 
          testCorrectness, printMatrices, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestUT<dcomplex>
        ( side, shape, orientation, m, offset, 
          testCorrectness, printMatrices, g );
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

