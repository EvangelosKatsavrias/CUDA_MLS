//   CUDA_MLS Framework
//
//   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
//
//   This file is part of the CUDA_MLS Framework.
//
//   CUDA_MLS Framework is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License version 3 as published by
//   the Free Software Foundation.
//
//   CUDA_MLS Framework is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
//
//   Contact Info:
//   Evangelos D. Katsavrias
//   email/skype: vageng@gmail.com
// -----------------------------------------------------------------------

#include<cmath>
#include"myBLAS.h"
#include<iostream>
#include"gaussElimination.h"
#include"conjugateGradient.h"
#include"LUDecomposers.h"
#include"LUSolvers.h"
#include"matrix_memoryStorage.h"
#include"LinearSolver.h"
#include<random>

void myBLAStest()
{

typedef double type_t;


type_t vec1[3]={1, 2, 3};
type_t out1[3];

//initializeVector(out1, out1+3, type_t(10));
//multiplicationOfVectorbyScalar(vec1, vec1+3, out1, type_t(10));
//tensorProductOfVectors(vec1, vec1+3, vec1, vec1+3, out1);
//for (int index=0; index<9; index++) std::cout << out1[index] << std::endl;
//sumOfVectors(vec1, vec1+3, vec1, out1);
type_t res = dotProductOfVectors(vec1, vec1+3, vec1);
std::cout << "Dot product result " << res << std::endl;


//for (int index=0; index<3; index++) std::cout << out1[index] << std::endl;



type_t A[] = {2, -1, 0, -1, 2, -1, 0, -1, 2};
type_t b[] = { 5, 5, 5 };

gaussElimination_backSubst(A, b, b, 3);

//for ( int i=0; i<3; i++) std::cout << b[i] << std::endl;


//type_t A1[] = { 0.256181, -0.0770212, 0.032368, 0.152549, -0.0412452, 0.0931575, -0.0770212, 0.032368, -0.0152552, -0.0412452, 0.0165393, -0.0224688, 0.032368, -0.0152552, 0.00793668, 0.0165393, -0.0074382, 0.00862191, 0.152549, -0.0412452, 0.0165393, 0.0931575, -0.0224688, 0.0582618, -0.0412452, 0.0165393, -0.0074382, -0.0224688, 0.00862191, -0.0124272, 0.0931575, -0.0224688, 0.00862191, 0.0582618, -0.0124272, 0.0372629 };
type_t A1[] = { 2, -1, 0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, -1, 2, -1, 0, 0, 0, 0, -1, 2 };


type_t A4[] = { 15,    -3,    -2,    -1,     0,     0,     0,     0,    -3,    15,    -3,    -2,    -1,     0,     0,     0,    -2,    -3,    15,    -3,    -2,    -1,     0,     0,    -1,    -2,    -3,    15,    -3,    -2,    -1,     0,     0,    -1,    -2,    -3,    15,    -3,    -2,    -1,     0,     0,    -1,    -2,    -3,    15,    -3,    -2,     0,     0,     0,    -1,    -2,    -3,    15,    -3,     0,     0,     0,     0,    -1,    -2,    -3,    15 };

type_t A5[] = {21.9298,    0.8363,    1.0050,    1.0762,    0.7429,    1.7175,    0.4790,    0.5741,    1.7158,    1.1940,    1.0707,    0.6988,    0.1168,    1.1999,    1.3697,    0.8363,   21.5155,    1.1819,    0.8906,    1.4962,    0.4570 ,   0.8682,    0.8576,    0.7083,    0.1221 ,   0.8774 ,   0.9924   , 1.6242    ,1.2561    ,1.2613    ,1.0050    ,1.1819   ,20.7631    ,1.7253    ,1.0495    ,0.9441    ,1.2840    ,0.9839    ,1.5570    ,0.9493    ,1.1444    ,0.9251    ,1.4554    ,1.5310    ,1.0860    ,1.0762    ,0.8906    ,1.7253   ,20.6808    ,1.3996    ,0.9775    ,1.0625    ,1.2512    ,0.6878    ,0.8224    ,1.1326    ,1.6910    ,1.0340    ,0.6888    ,0.5795    ,0.7429    ,1.4962    ,1.0495    ,1.3996   ,20.4870    ,1.3097    ,0.8785    ,0.2748    ,0.5149    ,0.8000    ,0.8772    ,0.9994    ,1.2667    ,1.2077    ,0.7345    ,1.7175    ,0.4570    ,0.9441    ,0.9775    ,1.3097   ,21.1356    ,0.2415    ,0.4966    ,0.6763    ,1.0191    ,1.0305    ,0.5808    ,1.0156    ,1.0564    ,0.9168    ,0.4790    ,0.8682    ,1.2840    ,1.0625    ,0.8785    ,0.2415   ,21.2040    ,0.3696    ,0.7901    ,1.1065    ,0.8801    ,0.9976    ,0.3902    ,0.4367    ,1.8931    ,0.5741    ,0.8576    ,0.9839    ,1.2512    ,0.2748    ,0.4966    ,0.3696   ,21.9238    ,0.8739    ,0.8246    ,1.7594    ,1.1650    ,0.5929    ,0.7010    ,0.6987    ,1.7158    ,0.7083    ,1.5570    ,0.6878    ,0.5149    ,0.6763    ,0.7901    ,0.8739   ,21.1594    ,1.4526    ,1.1011    ,1.5977    ,1.1328    ,0.8219    ,0.6244    ,1.1940    ,0.1221    ,0.9493    ,0.8224    ,0.8000    ,1.0191    ,1.1065    ,0.8246    ,1.4526   ,21.8896    ,1.0661    ,0.6782    ,1.1553    ,1.1305    ,0.6273    ,1.0707    ,0.8774    ,1.1444    ,1.1326    ,0.8772    ,1.0305    ,0.8801    ,1.7594    ,1.1011    ,1.0661   ,20.1196    ,0.9216    ,1.1480    ,1.6655    ,0.4241    ,0.6988    ,0.9924    ,0.9251    ,1.6910    ,0.9994    ,0.5808    ,0.9976    ,1.1650    ,1.5977    ,0.6782    ,0.9216   ,20.3670    ,1.0128    ,0.8204    ,1.3751    ,0.1168    ,1.6242    ,1.4554    ,1.0340    ,1.2667    ,1.0156    ,0.3902    ,0.5929    ,1.1328    ,1.1553    ,1.1480    ,1.0128   ,20.7572    ,1.0375    ,0.7950    ,1.1999    ,1.2561    ,1.5310    ,0.6888    ,1.2077    ,1.0564    ,0.4367    ,0.7010    ,0.8219    ,1.1305    ,1.6655    ,0.8204    ,1.0375   ,20.3414    ,0.8305    ,1.3697    ,1.2613    ,1.0860    ,0.5795    ,0.7345    ,0.9168    ,1.8931    ,0.6987    ,0.6244    ,0.6273    ,0.4241    ,1.3751    ,0.7950    ,0.8305   ,21.4224 };


std::default_random_engine defRandEng;
std::random_device randDevice;


int n = 15; 
type_t randMatrix[n*n], randMatrix_T[n*n], L[n*n], U[n*n], D[n];
type_t* AA = randMatrix;
std::fill(L, L+n*n, 0);


double randInterval = randDevice.max() - randDevice.min();
for (int i=0; i<n*n;i++) { randMatrix[i] = randDevice()/randInterval; randMatrix_T[i] = randMatrix[i]; }

for (int i=0; i<n;i++) for (int j = 0; j<n; j++) randMatrix[ i*n +j ] = randMatrix_T[ i*n +j ] +randMatrix_T[ j*n +i ];
for (int i=0; i<n;i++) randMatrix[ i*n +i ] += 10;
//for (int i = 0; i < n*n; i++) {std::cout << randMatrix[i] << "\t"; if ( (i+1)%n==0 ) std::cout << std::endl;}


//LU_Crout_decomposer (AA, L, U, n);
//LU_Doolittle_decomposer (AA, L, U, n);
//LDU_decomposer (AA, L, U, D, n);
LDL_decomposer (AA, L, D, n);
//Cholesky_decomposer (AA, L, n);


int numOfBands = find_numOfBands(AA, n);
type_t A_band[n*numOfBands], L_band[n*numOfBands];
std::fill(L_band, L_band+n*numOfBands, 0);
convert_matrixStorage_2Banded(AA, n, numOfBands, A_band);
Cholesky_decomposer_banded (A_band, A_band, n, numOfBands );


int n_max[n+1];
find_skylineColumns(AA, n, n_max);
type_t A_skyline[ n_max[n] ];
convert_matrixStorage_2Skyline(AA, n, n_max, A_skyline);
Cholesky_decomposer_skylineInPlace(A_skyline, n, n_max);


//std::cout << "A: \n";
//for ( int i=0; i<n*n; i++) std::cout << AA[i] << std::endl;
//std::cout << "\nL: \n";
//for ( int i=0; i<n*n; i++) std::cout << L[i] << std::endl;
//std::cout << "\nU: \n";
//for ( int i=0; i<n*n; i++) std::cout << U[i] << std::endl;
//std::cout << "\nD: \n";
//for ( int i=0; i<n; i++) std::cout << D[i] << std::endl;

//for ( int i=0; i<n*numOfBands; i++) std::cout << A_band[i] << std::endl;
//std::cout << "\nL_band: \n";
//for ( int i=0; i<n*numOfBands; i++) std::cout << L_band[i] << std::endl;

//for ( int i = 0; i < n+1; i++) std::cout << n_max[i] << std::endl;
//for ( int i = 0; i < n_max[n]; i++) std::cout << A_skyline[i] << std::endl;


type_t b1[n]; std::fill( b1, b1+n, 1 );
type_t x1[n], x2[n]; std::fill( x1, x1+n, 0 );
//LU_Crout_solver(L, U, b1, x1, n);
//LU_Doolittle_solver( L, U, b1, x1, n );
//LDU_solver(L, U, D, b1, x1, n);
//LDL_solver(L, D, b1, x1, n);
//Cholesky_solver(L, b1, x1, n);
//Cholesky_solver_banded(A_band, b1, x1, n, numOfBands );
//Cholesky_solver_skyline(A_skyline, b1, x1, n, n_max);
//CG(AA, b1, x1, n);
//gaussElimination_backSubst(AA, b1, x2, n);


//BandedMatrix<type_t> bmat(AA, n);
//Cholesky_decomposer_banded (bmat.data(), bmat.data(), n, bmat.get_numOfBands() );
//SkylineMatrix<type_t> skmat(AA, n);



GaussEliminationSolver<type_t> 	gaussEl(n);
ConjugateGradientSolver<type_t> conjGrad(n);
LUSolver<type_t> 		lusol(n);
LUCroutSolver<type_t> 		luCrout(n);
LDUSolver<type_t> 		lduSol(n);
LDLSolver<type_t> 		linSolver(n);
CholeskySolver<type_t> 		cholSolver(n);
CholeskyBandedSolver<type_t> 	bandedSolver(n);
//CholeskySkylineSolver<type_t> 	skylineSolver(n);
LinearSolver<type_t> 	skylineSolver = LinearSolver<type_t>(CholeskySkyline, n);



type_t err(0);
gaussEl		(AA, b1, x1);
conjGrad	(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
lusol		(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
luCrout		(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
lduSol		(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
linSolver	(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
//linSolver	(L, D, b1, x2);for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
cholSolver	(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
bandedSolver	(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
skylineSolver	(AA, b1, x2);  for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);


type_t det = determinant (lusol.get_upperMatrix(), n, 1);

//for (int i = 0; i < n*n; i++) {std::cout << randMatrix[i] << "\t"; if ( (i+1)%n==0 ) std::cout << std::endl;}


//for (int i=0; i<n; i++) err += pow((x1[i] - x2[i]), 2);
//for (int i=0; i<n; i++) std::cout << x2[i] << std::endl;
std::cout << "Error norm: " << sqrt(err) << std::endl;
std::cout << "Determinant: " << det << std::endl;

}
