/*
   Copyright (c) 2009-2011, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
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

   Written by Martin Schatz in March, 2011.
*/
#include "elemental.hpp"
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

void Usage()
{
    cout << "Strassen <n> <rProc> <cProc>\n"
         << "\n"
         << "  <n>: size of matrix to test.\n"
         << "  <rProc>: processor grid height.\n"
         << "  <cProc>: processor grid width.\n"
         << "\n"
         << "For now <n> must be a power of 2\n"
         << "and <rProc>, <cProc> must be\n"
         << "equal powers of 2."
         << std::endl;
}

int
main( int argc, char* argv[] )
{
    Init( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    int rank = mpi::CommRank( comm );

    if( argc < 4 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try 
    {
        const int n = atoi(argv[1]);
        const int r = atoi(argv[2]);
        const int c = atoi(argv[3]);

        bool test_square = false;
        for(int i = 2; i <= n; i *= 2)
            if(i == n){
                test_square = true;
                break;
            }
            
        if(!test_square){
            if( rank == 0)
                Usage();
            Finalize();
            return 0;
        }

        test_square = false;
        for(int i = 2; i <= n; i *= 2)
            if(i == n){
                test_square = true;
                break;
            }

        if(!(test_square && r == c)){
            if(rank == 0)
                Usage();
            Finalize();
            return 0;
        }

        Grid g( comm, r, c );

        DistMatrix<double,MC,MR> A( n, n, g ), B( n, n, g ), C( n, n, g);
        A.SetToRandom();
        B.SetToRandom();
        C.SetToRandom();

        basic::Strassen(A,B,C);    

        A.Print("A");
        B.Print("B");
        C.Print("C");

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

