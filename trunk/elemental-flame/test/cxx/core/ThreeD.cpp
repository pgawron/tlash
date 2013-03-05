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
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include "elemental.hpp"

using namespace std;
using namespace elemental;
using namespace elemental::imports;
using namespace elemental::basic;

void Usage()
{
    cout << "3-D Elemental Matrix multiply.\n\n"
         << "  ThreeD <r> <c> <h> <Am> <An> <Bm> <Bn>\n\n"
         << "  r: height of grid mesh\n"
         << "  d: width of grid mesh\n"
         << "  h: depth of grid mesh\n"
         << "  Am: height of matrix A\n"
         << "  An: width of matrix A\n"
         << "  Bm: height of matrix B\n"
         << "  Bn: width of matrix B\n" << endl;
}

int 
main( int argc, char* argv[] )
{
    Init( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    int rank = mpi::CommRank( comm );
    
    if( argc < 8 )
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
        const int h = atoi(argv[++argNum]);
        const int Am = atoi(argv[++argNum]);
        const int An = atoi(argv[++argNum]);
        const int Bm = atoi(argv[++argNum]);
        const int Bn = atoi(argv[++argNum]);

		if( An != Bm)
		{
			if( rank == 0 ){
				printf("Error: Dimension mis-match (%dx%d) * (%dx%d)\n", Am, An, Bm, Bn);
				Usage();
			}
			Finalize();
			return 0;
		}

#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        int p = mpi::CommSize( comm );
		mpi::Group group, depthGroup, meshGroup;
		mpi::Comm MDComm;

		int meshSize = r*c;

		int MDRank = rank / meshSize;
		int depthColor = rank % meshSize;

		mpi::CommSplit(comm, depthColor, MDRank, MDComm);
//		printf("Global Rank: %d has depthRank: %d\n", rank, mpi::CommRank(MDComm));
/*
		int sbuf[1];
		int rbuf[1];
		sbuf[0] = rank;
		mpi::AllReduce(&(sbuf[0]), &(rbuf[0]), 1, MPI_SUM, MDComm);

		printf("Rank: %d has reduce=%d\n", rank, rbuf[0]);		
*/
        std::vector<int> depthRanks(h);
        for( int i=0; i<h; ++i )
            depthRanks[i] = depthColor+(i*meshSize);

		std::vector<int> meshRanks(meshSize);
		for( int i = 0; i < meshSize; i++)
			meshRanks[i] = i+(MDRank*meshSize);
/*
		string pMe;
		for(int i = 0; i < depthRanks.size(); i++){
			char buf[30];
			sprintf(buf, "%d ", depthRanks[i]);
			pMe.append(buf);
		}
		printf("Rank: %d has r%m: %d MDComm %s\n", rank, rank % meshSize, pMe.c_str());	
*/
		mpi::CommGroup( comm, group );
		mpi::GroupIncl( group, meshRanks.size(), &meshRanks[0], meshGroup);

		const Grid meshGrid( comm, meshGroup, r, c );

		double *sABuf;
		int sAC = utilities::LocalLength(Am, meshGrid.MCRank(), 0, meshGrid.Height()) * 
				  utilities::LocalLength(An, meshGrid.MRRank(), 0, meshGrid.Width());
		int rAC = sAC / h;

		double *rABuf = new double[rAC];

		double *sBBuf;
		int sBC = utilities::LocalLength(Bm, meshGrid.MCRank(), 0, meshGrid.Height()) * 
				  utilities::LocalLength(Bn, meshGrid.MRRank(), 0, meshGrid.Width());
		int rBC = sBC / h;
		double *rBBuf = new double[rBC];

		if(rank % meshSize == 0){
			printf("Rank: %d sAC: %d rAC: %d sBC: %d rBC: %d\n", rank, sAC, rAC, sBC, rBC);
		}

		if(rank < meshSize){
			printf("meshSize %d x %d\n", meshGrid.Height(), meshGrid.Width());
			DistMatrix<double, MC, MR> A(Am, An, meshGrid);
			A.SetToIdentity();
			basic::Gemm(Normal, Normal, (double)2, A, A, (double)2, A);

			sABuf = new double[sAC];

			const double *tmp = A.LocalBuffer();
			for(int i = 0; i < sAC; i++)
				sABuf[i] = tmp[i];

			double *tmpData = new double[sBC];
			for(int i = 0; i < sBC; i++)
				tmpData[i] = i*meshSize + rank;
			DistMatrix<double, MC, MR> B(Bm, Bn, 0, 0, tmpData, Bm/meshGrid.Height(), meshGrid);
			B.Print("Initial B");

			B.SetToIdentity();
			basic::Gemm(Normal, Normal, (double)1, B, B, (double)1, B);
			DistMatrix<double, MR, Star> B_MR_Star(meshGrid);
			B_MR_Star.TransposeFrom(B);
			B = B_MR_Star;
			B.Print("Initial BT");
			tmp = B.LocalBuffer();

			sBBuf = new double[sBC];

			for(int i = 0; i < sBC; i++)
				sBBuf[i] = tmp[i];

/*
			sBBuf = new double[Bm*Bn/meshSize];
			for(int i = 0; i < Bm*Bn/meshSize; i++)
				sBBuf[i] = (rank % meshSize)+i*meshSize;
*/
//			A.Print("A");
//			B.Print("B");

			string pMe;
			for(int i = 0; i < sBC; i++){
				char buf[10];
				sprintf(buf, "%.0f ", sBBuf[i]);
				pMe.append(buf);
			}
			printf("Rank: %d has sBBuf: %s\n", rank, pMe.c_str());	
			printf("Rank: testing\n");
			DistMatrix<double, MC, MR> testing(Bn, Bm/h, 0, 0, sBBuf, Bn/meshGrid.Height(), meshGrid);
			testing.Print("testing");
		}

//		printf("Rank: %d bling\n", rank);
		if(rank < meshSize){
			printf("Rank: %d sBC: %d rBC: %d sBC/h: %d \n", rank, sBC, rBC, sBC/h);

		}

		mpi::Scatter(sABuf, sAC/h,
					 rABuf, rAC, 0, MDComm);

		mpi::Scatter(sBBuf, sBC/h,
					 rBBuf, rBC, 0, MDComm);
/*
		if(rank < meshSize){
			delete [] sABuf;
			delete [] sBBuf;
		}
*/
		if(rank < meshSize){
			string printMe;
			for(int i = 0; i < rBC; i++){
				char buf[30];
				sprintf(buf, "%.0f ", rBBuf[i]);
				printMe.append(buf);
				printMe.append(" ");
			}
			printf("Rank: %d has values %s\n", rank, printMe.c_str());
		}

		DistMatrix<double, MC, MR> A(Am, An/h, 0, 0, rABuf, Am/meshGrid.Height(), meshGrid);
//		DistMatrix<double, MC, Star> A_MC_Star(meshGrid);
//		A_MC_Star = A;

		DistMatrix<double, MC, MR> B(Bn, Bm/h, 0, 0, rBBuf, Bm/meshGrid.Height(), meshGrid);
		if(rank < meshSize)
			B.Print("B");
		DistMatrix<double, MC, MR> B_MC_MR(meshGrid);
		B_MC_MR = B;
		DistMatrix<double, MR, Star> B_MR_Star(meshGrid);
		B_MR_Star.TransposeFrom(B);
//		B_MC_MR = B_MR_Star;

		DistMatrix<double, MC, MR> newB(meshGrid);
		newB = B_MR_Star;

		if(rank < meshSize){
		A.Print("A");
		newB.Print("newB");
		}

		DistMatrix<double, MC, MR> C(Am, Bn, meshGrid);
		C.SetToZero();

//		C.Print("C");

		basic::Gemm(Normal, Normal, (double)1, A, newB, (double)1, C);
//		C.Print("C after");

		int CLocalBufSize = C.LocalHeight() * C.LocalWidth();
		double *sCBuf = C.LocalBuffer();
		double *rCBuf = new double[CLocalBufSize];

		mpi::AllReduce(sCBuf, rCBuf, CLocalBufSize, MPI_SUM, MDComm);

		DistMatrix<double, MC, MR> resC(Am, Bn, 0, 0, rCBuf, Am/meshGrid.Height(), meshGrid);
		if(rank < meshSize)
			resC.Print("result");

/*
        // Drop down to a square grid, change the matrix, and redistribute back
        int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

        std::vector<int> sqrtRanks(pSqrt*pSqrt);
        for( int i=0; i<pSqrt*pSqrt; ++i )
            sqrtRanks[i] = i;

        mpi::Group group, sqrtGroup;
        
        mpi::CommGroup( comm, group );
        mpi::GroupIncl( group, sqrtRanks.size(), &sqrtRanks[0], sqrtGroup );

        const Grid grid( comm );
        const Grid sqrtGrid( comm, sqrtGroup );

        DistMatrix<double,MC,MR> A( m, n, grid );
        DistMatrix<double,MC,MR> ASqrt( m, n, sqrtGrid );

        A.SetToIdentity();
        A.Print("A");

        ASqrt = A;
        ASqrt.Print("ASqrt := A");

        basic::Scal( 2.0, ASqrt );
        ASqrt.Print("ASqrt := 2 ASqrt");

        A = ASqrt;
        A.Print("A := ASqrt");
*/
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

