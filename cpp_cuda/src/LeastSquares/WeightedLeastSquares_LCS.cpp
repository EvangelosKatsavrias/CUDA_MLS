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

#include"LeastSquares.h"
#include"weightedLeastSquares_cuda.h"
#include"myBLAS.h"
#include<cmath>
#include<map>
#include<fstream>
#include<thread>
#include<functional>
#include<cuda_runtime.h>
#include<algorithm>
#include<numeric>


template<class T>
std::string WeightedLeastSquares_LCS<T>::getLSType()
{ return "Weighted least squares in LCS"; }

std::ofstream ff, ff2;


template<class T>
void buildSystemAtSamplePoint_LCS3D(T distance, T samplePoint_x, T samplePoint_y, T samplePoint_z, T sampleValue_x, T sampleValue_y, T sampleValue_z, T* Sum_bbT, T* Sum_bf, T stationaryPoint_x, T stationaryPoint_y, T stationaryPoint_z, int numOfMonomials, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction)
{

	T bbT[numOfMonomials*numOfMonomials];
	T b[numOfMonomials];
	T theta;
	T theta_b[numOfMonomials];
	T theta_bx[numOfMonomials], theta_by[numOfMonomials], theta_bz[numOfMonomials];

	theta = weightingFunction( distance );


	basisFunctions			(b, samplePoint_x -stationaryPoint_x, samplePoint_y -stationaryPoint_y, samplePoint_z -stationaryPoint_z );
//	for (int i = 0; i < numOfMonomials; i++ ) std::cout << b[i] << std::endl;
	multiplicationOfVectorbyScalar	(b, b+numOfMonomials, theta_b, theta);


	tensorProductOfVectors		(theta_b, theta_b+numOfMonomials, b, b+numOfMonomials, bbT);
	sumOfVectors			(bbT, bbT+numOfMonomials*numOfMonomials, Sum_bbT, Sum_bbT);


	multiplicationOfVectorbyScalar	(theta_b, theta_b+numOfMonomials, theta_bx, sampleValue_x);
	multiplicationOfVectorbyScalar	(theta_b, theta_b+numOfMonomials, theta_by, sampleValue_y);
	multiplicationOfVectorbyScalar	(theta_b, theta_b+numOfMonomials, theta_bz, sampleValue_z);
	sumOfVectors			(theta_bx, theta_bx+numOfMonomials, Sum_bf, Sum_bf);
	sumOfVectors			(theta_by, theta_by+numOfMonomials, Sum_bf+numOfMonomials, Sum_bf+numOfMonomials);
	sumOfVectors			(theta_bz, theta_bz+numOfMonomials, Sum_bf+2*numOfMonomials, Sum_bf+2*numOfMonomials);

}



template<class T>
void findCoefficientsAtStationaryPoint_LCS3D(T stationaryPoint_x, T stationaryPoint_y, T stationaryPoint_z, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, int numOfMonomials, LinearSolver<T>& linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, std::vector<T>& samplePointsDistance, T* Sum_bf)
{

	T Sum_bbT[numOfMonomials*numOfMonomials]; std::fill(Sum_bbT, Sum_bbT+numOfMonomials*numOfMonomials, T(0));
	T Sum_bfx[3*numOfMonomials]; std::fill(Sum_bfx, Sum_bfx+3*numOfMonomials, T(0) );

	for ( int samplePointIndex=0; samplePointIndex < numOfSamplePoints; samplePointIndex++ )
	{
		if ( samplePointsDistance[samplePointIndex] )
		buildSystemAtSamplePoint_LCS3D<T>( samplePointsDistance[samplePointIndex], samplePoints_x[samplePointIndex], samplePoints_y[samplePointIndex], samplePoints_z[samplePointIndex], sampleValues_x[samplePointIndex], sampleValues_y[samplePointIndex], sampleValues_z[samplePointIndex], Sum_bbT, Sum_bfx, stationaryPoint_x, stationaryPoint_y, stationaryPoint_z, numOfMonomials, basisFunctions, weightingFunction);
	}

//	for ( int i=0; i<numOfMonomials*numOfMonomials; i++ ) std::cout << Sum_bbT[i] << std::endl;
//	for ( int i=0; i<2*numOfMonomials; i++ ) std::cout << Sum_bfx[i] << std::endl;

	linearSolver.set_givenMatrix(0);
	linearSolver(Sum_bbT, Sum_bfx, Sum_bf );
	linearSolver(Sum_bbT, Sum_bfx+numOfMonomials, Sum_bf+numOfMonomials);
	linearSolver(Sum_bbT, Sum_bfx+2*numOfMonomials, Sum_bf+2*numOfMonomials);

//	for ( int i=0; i<2*numOfMonomials; i++ ) ff << Sum_bfx[i] << std::endl;
//	for ( int i=0; i<2*numOfMonomials; i++ ) std::cout << Sum_bf[i] << std::endl;

}



template<class T>
void cpuVersion_findCoefficients_LCS3D(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, int numOfMonomials, T* coefficients, LinearSolver<T>& linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, bool progressIndexFlag )
{	
	T Sum_bf[3*numOfMonomials], distance; T span( weightingFunction.get_span() ); 
	std::vector<T> samplePointsDistance; samplePointsDistance.resize( numOfSamplePoints );
	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0);


	if ( progressIndexFlag ) std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); 
//	clock_t calculationClocks = -clock();
//	this->m_calculationClocks = -clock();
//	ff.open( "bf_cpu.dat" ); ff2.open( "bbT_cpu.dat" );


	for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPoints; stationaryPointIndex++)
	{
		int counter(0); std::fill(samplePointsDistance.begin(), samplePointsDistance.end(), 0 );

		for ( size_t samplePointIndex = 0; samplePointIndex<numOfSamplePoints; samplePointIndex++ ) {
			distance = distance3D(stationaryPoints_x[stationaryPointIndex], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex], samplePoints_y[samplePointIndex], stationaryPoints_z[stationaryPointIndex], samplePoints_z[samplePointIndex]);
			if ( distance < span ) { samplePointsDistance[samplePointIndex]=distance; counter++; }
		}
		if ( counter < numOfMonomials ) {
			std::fill(coefficients +stationaryPointIndex    *3*numOfMonomials,
				  coefficients +(stationaryPointIndex+1)*3*numOfMonomials,
				  T(0) ); continue; }

		numOfEffectiveStationaryPoints++;
		std::fill(Sum_bf, Sum_bf+3*numOfMonomials, T(0));

		findCoefficientsAtStationaryPoint_LCS3D ( stationaryPoints_x[stationaryPointIndex], stationaryPoints_y[stationaryPointIndex], stationaryPoints_z[stationaryPointIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, samplePoints_z, sampleValues_x, sampleValues_y, sampleValues_z, numOfMonomials, linearSolver, basisFunctions, weightingFunction, samplePointsDistance, Sum_bf );

		for (int index=0; index<3*numOfMonomials; index++) coefficients[ stationaryPointIndex*3*numOfMonomials +index ] = Sum_bf[index];

		totalDisplacement_L2norm += Sum_bf[0]*Sum_bf[0] +Sum_bf[numOfMonomials]*Sum_bf[numOfMonomials] +Sum_bf[2*numOfMonomials]*Sum_bf[2*numOfMonomials];

		if ( progressIndexFlag ) std::cout << "\r" << stationaryPointIndex*nodesFraction << "% "; std::cout.flush();
	}

//	ff.close( ); ff2.close( );

	if ( progressIndexFlag ) std::cout << "\r" << "100.00%\n\n";
	std::cout << "CPU process" << std::endl;
	std::cout << "Number of effective stationary points: \t\t" << numOfEffectiveStationaryPoints << std::endl;
	std::cout << "Number of outlier stationary points: \t\t" << numOfStationaryPoints -numOfEffectiveStationaryPoints << std::endl;
	std::cout << "L2 norm of the evaluated 2D displacement field: " << sqrt(totalDisplacement_L2norm) << std::endl << std::endl;

}




template<class T>
void WeightedLeastSquares_LCS<T>::findCoefficients_LCS3D( )
{

	int numOfSamplePoints 			= this->m_scatteredData->get_numOfNodes();
	T* samplePoints_x 			= this->m_scatteredData->get_domainComponent(0);
	T* samplePoints_y 			= this->m_scatteredData->get_domainComponent(1);
	T* samplePoints_z 			= this->m_scatteredData->get_domainComponent(2);
	T* sampleValues_x 			= this->m_scatteredData->get_fieldComponent(0);
	T* sampleValues_y 			= this->m_scatteredData->get_fieldComponent(1);
	T* sampleValues_z 			= this->m_scatteredData->get_fieldComponent(2);

	T* stationaryPoints_x			= this->m_stationaryPoints->get_component(0);
	T* stationaryPoints_y			= this->m_stationaryPoints->get_component(1);
	T* stationaryPoints_z			= this->m_stationaryPoints->get_component(2);
	int numOfStationaryPoints		= this->get_numOfStationaryPoints();
	int numOfMonomials 			= this->m_basisFunctions->get_numOfMonomials();
	int polDegree				= this->m_basisFunctions->get_degree();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	LinearSolver<T>& linearSolver		= *(this->m_linearSolver);
	linearSolver.set_numberOfColumns(numOfMonomials);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);

	this->m_coefficients.resize(3*numOfStationaryPoints*numOfMonomials);

	this->m_calculationClocks = -clock();

	switch (m_computingMode) {
		case CPU:
			cpuVersion_findCoefficients_LCS3D<T>(numOfStationaryPoints, stationaryPoints_x, stationaryPoints_y, stationaryPoints_z, numOfSamplePoints, samplePoints_x, samplePoints_y, samplePoints_z, sampleValues_x, sampleValues_y, sampleValues_z, numOfMonomials, this->m_coefficients.data(), linearSolver, basisFunctions, weightingFunction, 1 ); 
			break;

		case CPU_MT: {
		/*		size_t numOfThreads = std::thread::hardware_concurrency();
				std::cout << "Number of available cpu threads: " << numOfThreads << std::endl;

				size_t numOfPoints[ numOfThreads ];
				size_t typicalSize = numOfStationaryPoints/numOfThreads;
				size_t boundSize = numOfStationaryPoints -(numOfThreads-1)*typicalSize;
				for ( size_t tIndex = 0; tIndex < numOfThreads-1; tIndex++ ) numOfPoints[tIndex] = typicalSize;
				numOfPoints[numOfThreads-1] = boundSize;

				size_t pointShifts[numOfThreads]; pointShifts[0] = 0;
				std::partial_sum( numOfPoints, numOfPoints+numOfThreads-1, pointShifts+1 );

				std::vector<std::thread> threads;
				for ( size_t tIndex = 0; tIndex < numOfThreads; tIndex++ )
				{
					threads.push_back( std::thread( cpuVersion_findCoefficients_LCS2D<T>, numOfPoints[tIndex], stationaryPoints_x +pointShifts[tIndex], stationaryPoints_y +pointShifts[tIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, this->m_coefficients.data() +pointShifts[tIndex]*2*numOfMonomials, std::ref(linearSolver), std::ref(basisFunctions), std::ref(weightingFunction), 0 )  );
				}
				for ( auto& thread : threads ) thread.join();

		*/		break; 
		}

		case CUDA:
			cudaVersion_findCoefficientsLCS3D2<T>(numOfStationaryPoints, stationaryPoints_x, stationaryPoints_y, stationaryPoints_z, numOfSamplePoints, samplePoints_x, samplePoints_y, samplePoints_z, sampleValues_x, sampleValues_y, sampleValues_z, numOfMonomials, polDegree, this->m_coefficients.data(), linearSolver, weightingFunction, 0); 
			break;
		case CPU_CUDA:{ 
		/*		size_t numOfThreads = std::thread::hardware_concurrency() -1;
				std::cout << "Number of available cpu threads: " << numOfThreads << std::endl;

				size_t numOfPoints[ numOfThreads ];

				size_t gpuSize = 0.5*numOfStationaryPoints;
				size_t typicalSize = 0.5*numOfStationaryPoints/(numOfThreads-1);
				size_t boundSize = numOfStationaryPoints -gpuSize -(numOfThreads-2)*typicalSize;
				numOfPoints[0] = gpuSize;
				for ( size_t tIndex = 1; tIndex < numOfThreads-1; tIndex++ ) numOfPoints[tIndex] = typicalSize;
				numOfPoints[numOfThreads-1] = boundSize;

				size_t pointShifts[numOfThreads]; pointShifts[0] = 0;
				std::partial_sum( numOfPoints, numOfPoints+numOfThreads-1, pointShifts+1 );

				std::vector<std::thread> threads;

				threads.push_back( std::thread( cudaVersion_findCoefficientsLCS2D<T>, gpuSize, stationaryPoints_x, stationaryPoints_y, numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, polDegree, this->m_coefficients.data(), std::ref(linearSolver), std::ref(weightingFunction), 0 ) );

				for ( size_t tIndex = 1; tIndex < numOfThreads; tIndex++ )
				{
					threads.push_back( std::thread( cpuVersion_findCoefficients_LCS2D<T>, numOfPoints[tIndex], stationaryPoints_x +pointShifts[tIndex], stationaryPoints_y +pointShifts[tIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, this->m_coefficients.data() +pointShifts[tIndex]*2*numOfMonomials, std::ref(linearSolver), std::ref(basisFunctions), std::ref(weightingFunction), 0 ));
				}

				for ( auto& thread : threads ) thread.join();

		*/		break;
		}
		case MultiGPU_CUDA: { 
		/*		int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
				std::cout << "\nNumber of cuda enabled GPU's: " << dev_count << std::endl;

				size_t numOfPoints[ dev_count ];
				size_t typicalSize = numOfStationaryPoints/dev_count;
				size_t boundSize = numOfStationaryPoints -(dev_count-1)*typicalSize;
				for ( size_t dIndex = 0; dIndex < dev_count-1; dIndex++ ) numOfPoints[dIndex] = typicalSize;
				numOfPoints[dev_count-1] = boundSize;

				size_t pointShifts[dev_count]; pointShifts[0] = 0;
				std::partial_sum( numOfPoints, numOfPoints+dev_count-1, pointShifts+1 );

				std::vector<std::thread> threads;
				for ( size_t dIndex = 0; dIndex < dev_count; dIndex++ )
				{
					threads.push_back( std::thread( cudaVersion_findCoefficientsLCS2D<T>, numOfPoints[dIndex], stationaryPoints_x +pointShifts[dIndex], stationaryPoints_y +pointShifts[dIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, polDegree, this->m_coefficients.data() +2*numOfMonomials*pointShifts[dIndex], std::ref(linearSolver), std::ref(weightingFunction), dIndex ) );
				}
				for ( auto& thread : threads ) thread.join();
		*/		break;
		}
	}
}


float mytime(0);

template<class T>
void buildSystemAtSamplePoint_LCS2D(T distance, T samplePoint_x, T samplePoint_y, T sampleValue_x, T sampleValue_y, T* Sum_bbT, T* Sum_bf, T stationaryPoint_x, T stationaryPoint_y, int numOfMonomials, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction)
{
	T bbT[numOfMonomials*numOfMonomials];
	T b[numOfMonomials];
	T theta;
	T theta_b[numOfMonomials];
	T theta_bx[numOfMonomials], theta_by[numOfMonomials];

	theta = weightingFunction( distance );

	basisFunctions			(b, samplePoint_x -stationaryPoint_x, samplePoint_y -stationaryPoint_y );
	multiplicationOfVectorbyScalar	(b, b+numOfMonomials, theta_b, theta);

//	clock_t calculationClocks = -clock();

	tensorProductOfVectors		(theta_b, theta_b+numOfMonomials, b, b+numOfMonomials, bbT);
	sumOfVectors			(bbT, bbT+numOfMonomials*numOfMonomials, Sum_bbT, Sum_bbT);

//	calculationClocks += clock();
//	mytime += ((float)calculationClocks*1000)/CLOCKS_PER_SEC;

	multiplicationOfVectorbyScalar	(theta_b, theta_b+numOfMonomials, theta_bx, sampleValue_x);
	multiplicationOfVectorbyScalar	(theta_b, theta_b+numOfMonomials, theta_by, sampleValue_y);
	sumOfVectors			(theta_bx, theta_bx+numOfMonomials, Sum_bf, Sum_bf);
	sumOfVectors			(theta_by, theta_by+numOfMonomials, Sum_bf+numOfMonomials, Sum_bf+numOfMonomials);

}



template<class T>
void findCoefficientsAtStationaryPoint_LCS2D(T stationaryPoint_x, T stationaryPoint_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, LinearSolver<T>& linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, std::vector<T>& samplePointsDistance, T* Sum_bf)
{

//	clock_t calculationClocks = -clock();

	T Sum_bbT[numOfMonomials*numOfMonomials]; std::fill(Sum_bbT, Sum_bbT+numOfMonomials*numOfMonomials, T(0));
	T Sum_bfx[2*numOfMonomials]; std::fill(Sum_bfx, Sum_bfx+2*numOfMonomials, T(0) );


	for ( int samplePointIndex=0; samplePointIndex < numOfSamplePoints; samplePointIndex++ )
	{
		if ( samplePointsDistance[samplePointIndex] )
		buildSystemAtSamplePoint_LCS2D<T>( samplePointsDistance[samplePointIndex], samplePoints_x[samplePointIndex], samplePoints_y[samplePointIndex], sampleValues_x[samplePointIndex], sampleValues_y[samplePointIndex], Sum_bbT, Sum_bfx, stationaryPoint_x, stationaryPoint_y, numOfMonomials, basisFunctions, weightingFunction);
		else weightingFunction(1);
	}

//	for ( int i=0; i<numOfMonomials*numOfMonomials; i++ ) std::cout << Sum_bbT[i] << std::endl;
//	for ( int i=0; i<2*numOfMonomials; i++ ) std::cout << Sum_bfx[i] << std::endl;

//	calculationClocks += clock();
//	mytime += (float)calculationClocks/CLOCKS_PER_SEC;


	linearSolver.set_givenMatrix(0);
	linearSolver(Sum_bbT, Sum_bfx, Sum_bf );
	linearSolver(Sum_bbT, Sum_bfx+numOfMonomials, Sum_bf+numOfMonomials);


//	for ( int i=0; i<2*numOfMonomials; i++ ) ff << Sum_bfx[i] << std::endl;
//	for ( int i=0; i<2*numOfMonomials; i++ ) std::cout << Sum_bf[i] << std::endl;

}


template<class T>
void cpuVersion_findCoefficients_LCS2D(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, T* coefficients, LinearSolver<T>& linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, bool progressIndexFlag )
{
	T Sum_bf[2*numOfMonomials], distance;// T span( weightingFunction.get_span() ); 
	std::vector<T> samplePointsDistance; samplePointsDistance.resize( numOfSamplePoints );
	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0);

	if ( progressIndexFlag ) std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); 

//	this->m_calculationClocks = -clock();
//	ff.open( "bf_cpu.dat" ); ff2.open( "bbT_cpu.dat" );

//	float time(0);
//	clock_t calculationClocks = -clock();
	for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPoints; stationaryPointIndex++)
	{

//		clock_t calculationClocks = -clock();
		int counter(0); std::fill(samplePointsDistance.begin(), samplePointsDistance.end(), 0 );

		for ( size_t samplePointIndex = 0; samplePointIndex<numOfSamplePoints; samplePointIndex++ ) {
			distance = distance2D(stationaryPoints_x[stationaryPointIndex], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex], samplePoints_y[samplePointIndex]);
			//if ( distance < weightingFunction.get_span( samplePointIndex ) ) { 
				samplePointsDistance[samplePointIndex]=distance; //counter++; }
		}
//		if ( counter < numOfMonomials ) {
//			std::fill(coefficients +stationaryPointIndex    *2*numOfMonomials,
//				  coefficients +(stationaryPointIndex+1)*2*numOfMonomials,
//				  T(0) ); continue; }


		numOfEffectiveStationaryPoints++;
		std::fill(Sum_bf, Sum_bf+2*numOfMonomials, T(0));

		findCoefficientsAtStationaryPoint_LCS2D ( stationaryPoints_x[stationaryPointIndex], stationaryPoints_y[stationaryPointIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, linearSolver, basisFunctions, weightingFunction, samplePointsDistance, Sum_bf );

		for (int index=0; index<2*numOfMonomials; index++) coefficients[ stationaryPointIndex*2*numOfMonomials +index ] = Sum_bf[index];

		totalDisplacement_L2norm += Sum_bf[0]*Sum_bf[0] +Sum_bf[numOfMonomials]*Sum_bf[numOfMonomials];

//		calculationClocks += clock();
//		time += (float)calculationClocks/CLOCKS_PER_SEC;

		if ( progressIndexFlag ) std::cout << "\r" << stationaryPointIndex*nodesFraction << "% "; std::cout.flush();
	}

//	calculationClocks += clock();
//	time += (float)calculationClocks/CLOCKS_PER_SEC;
	std::cout << " ============================> CPU Calculation time: " << mytime << "msec" << std::endl;

//	ff.close( ); ff2.close( );

	if ( progressIndexFlag ) std::cout << "\r" << "100.00%\n\n";
	std::cout << "CPU process" << std::endl;
	std::cout << "Number of effective stationary points: \t\t" << numOfEffectiveStationaryPoints << std::endl;
	std::cout << "Number of outlier stationary points: \t\t" << numOfStationaryPoints -numOfEffectiveStationaryPoints << std::endl;
	std::cout << "L2 norm of the evaluated 2D displacement field: " << sqrt(totalDisplacement_L2norm) << std::endl << std::endl;
}



template<class T>
void WeightedLeastSquares_LCS<T>::findCoefficients_LCS2D( )
{
	int numOfSamplePoints 			= this->m_scatteredData->get_numOfNodes();
	T* samplePoints_x 			= this->m_scatteredData->get_domainComponent(0);
	T* samplePoints_y 			= this->m_scatteredData->get_domainComponent(1);
	T* sampleValues_x 			= this->m_scatteredData->get_fieldComponent(0);
	T* sampleValues_y 			= this->m_scatteredData->get_fieldComponent(1);

	T* stationaryPoints_x			= this->m_stationaryPoints->get_component(0);
	T* stationaryPoints_y			= this->m_stationaryPoints->get_component(1);
	int numOfStationaryPoints		= this->get_numOfStationaryPoints();
	int numOfMonomials 			= this->m_basisFunctions->get_numOfMonomials();
	int polDegree				= this->m_basisFunctions->get_degree();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	LinearSolver<T>& linearSolver		= *(this->m_linearSolver);
	linearSolver.set_numberOfColumns(numOfMonomials);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);

	this->m_coefficients.resize(2*numOfStationaryPoints*numOfMonomials);

	this->m_calculationClocks = -clock();

	switch (m_computingMode) {
		case CPU: { cpuVersion_findCoefficients_LCS2D<T>(numOfStationaryPoints, stationaryPoints_x, stationaryPoints_y, numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, this->m_coefficients.data(), linearSolver, basisFunctions, weightingFunction, 1 ); break; }

		case CPU_MT: {

				weightingFunction.set_numOfVariableData(0);
				size_t numOfThreads = std::thread::hardware_concurrency();
				std::cout << "Number of available cpu threads: " << numOfThreads << std::endl;

				size_t numOfPoints[ numOfThreads ];
				size_t typicalSize = numOfStationaryPoints/numOfThreads;
				size_t boundSize = numOfStationaryPoints -(numOfThreads-1)*typicalSize;
				for ( size_t tIndex = 0; tIndex < numOfThreads-1; tIndex++ ) numOfPoints[tIndex] = typicalSize;
				numOfPoints[numOfThreads-1] = boundSize;

				size_t pointShifts[numOfThreads]; pointShifts[0] = 0;
				std::partial_sum( numOfPoints, numOfPoints+numOfThreads-1, pointShifts+1 );

				std::vector<std::thread> threads;
				for ( size_t tIndex = 0; tIndex < numOfThreads; tIndex++ )
				{
					threads.push_back( std::thread( cpuVersion_findCoefficients_LCS2D<T>, numOfPoints[tIndex], stationaryPoints_x +pointShifts[tIndex], stationaryPoints_y +pointShifts[tIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, this->m_coefficients.data() +pointShifts[tIndex]*2*numOfMonomials, std::ref(linearSolver), std::ref(basisFunctions), std::ref(weightingFunction), 0 )  );
				}
				for ( auto& thread : threads ) thread.join();

			break; }


		case CUDA: { cudaVersion_findCoefficientsLCS2D2<T>(numOfStationaryPoints, stationaryPoints_x, stationaryPoints_y, numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, polDegree, this->m_coefficients.data(), linearSolver, weightingFunction, 0); break; }


		case CPU_CUDA:{
				weightingFunction.set_numOfVariableData(0);
				size_t numOfThreads = std::thread::hardware_concurrency() -1;
				std::cout << "Number of available cpu threads: " << numOfThreads << std::endl;

				size_t numOfPoints[ numOfThreads ];

				float gpuLoadProportion = 0.5;
				size_t gpuSize = gpuLoadProportion*numOfStationaryPoints;
				size_t typicalSize = gpuLoadProportion*numOfStationaryPoints/(numOfThreads-1);
				size_t boundSize = numOfStationaryPoints -gpuSize -(numOfThreads-2)*typicalSize;
				numOfPoints[0] = gpuSize;
				for ( size_t tIndex = 1; tIndex < numOfThreads-1; tIndex++ ) numOfPoints[tIndex] = typicalSize;
				numOfPoints[numOfThreads-1] = boundSize;

				size_t pointShifts[numOfThreads]; pointShifts[0] = 0;
				std::partial_sum( numOfPoints, numOfPoints+numOfThreads-1, pointShifts+1 );

				std::vector<std::thread> threads;

				threads.push_back( std::thread( cudaVersion_findCoefficientsLCS2D<T>, gpuSize, stationaryPoints_x, stationaryPoints_y, numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, polDegree, this->m_coefficients.data(), std::ref(linearSolver), std::ref(weightingFunction), 0 ) );

				for ( size_t tIndex = 1; tIndex < numOfThreads; tIndex++ )
				{
					threads.push_back( std::thread( cpuVersion_findCoefficients_LCS2D<T>, numOfPoints[tIndex], stationaryPoints_x +pointShifts[tIndex], stationaryPoints_y +pointShifts[tIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, this->m_coefficients.data() +pointShifts[tIndex]*2*numOfMonomials, std::ref(linearSolver), std::ref(basisFunctions), std::ref(weightingFunction), 0 ));
				}

				for ( auto& thread : threads ) thread.join();

			break;}
		case MultiGPU_CUDA: {

				weightingFunction.set_numOfVariableData(0);
				int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
				std::cout << "\nNumber of cuda enabled GPU's: " << dev_count << std::endl;

				size_t numOfPoints[ dev_count ];
				size_t typicalSize = numOfStationaryPoints/dev_count;
				size_t boundSize = numOfStationaryPoints -(dev_count-1)*typicalSize;
				for ( size_t dIndex = 0; dIndex < dev_count-1; dIndex++ ) numOfPoints[dIndex] = typicalSize;
				numOfPoints[dev_count-1] = boundSize;

				size_t pointShifts[dev_count]; pointShifts[0] = 0;
				std::partial_sum( numOfPoints, numOfPoints+dev_count-1, pointShifts+1 );

				std::vector<std::thread> threads;
				for ( size_t dIndex = 0; dIndex < dev_count; dIndex++ )
				{
					threads.push_back( std::thread( cudaVersion_findCoefficientsLCS2D2<T>, numOfPoints[dIndex], stationaryPoints_x +pointShifts[dIndex], stationaryPoints_y +pointShifts[dIndex], numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, polDegree, this->m_coefficients.data() +2*numOfMonomials*pointShifts[dIndex], std::ref(linearSolver), std::ref(weightingFunction), dIndex ) );
				}
				for ( auto& thread : threads ) thread.join();

			break;}
	}
}


template<class T>
void buildSystemAtSamplePoint_LCS(T samplePoint, T sampleValue, T* Sum_bbT, T* Sum_bf, T stationaryPoint, int numOfMonomials, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction)
{
	T bbT[numOfMonomials*numOfMonomials];
	T b[numOfMonomials];
	T theta;
	T theta_b[numOfMonomials];

	theta = weightingFunction(fabs(samplePoint - stationaryPoint));
	basisFunctions(b, samplePoint -stationaryPoint);
	multiplicationOfVectorbyScalar(b, b+numOfMonomials, theta_b, theta);
	tensorProductOfVectors(theta_b, theta_b+numOfMonomials, b, b+numOfMonomials, bbT);
	sumOfVectors(bbT, bbT+numOfMonomials*numOfMonomials, Sum_bbT, Sum_bbT);
	multiplicationOfVectorbyScalar(theta_b, theta_b+numOfMonomials, theta_b, sampleValue);
	sumOfVectors(theta_b, theta_b+numOfMonomials, Sum_bf, Sum_bf);
}


template<class T>
void findCoefficientsAtStationaryPoint_LCS(T stationaryPoint, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, LinearSolver<T>& linearSolver, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction, T* Sum_bf)
{
	T Sum_bbT[numOfMonomials*numOfMonomials]; std::fill(Sum_bbT, Sum_bbT+numOfMonomials*numOfMonomials, T(0));

	for (int samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++)
	{
		buildSystemAtSamplePoint_LCS<T>(samplePoints[samplePointIndex], sampleValues[samplePointIndex], Sum_bbT, Sum_bf, stationaryPoint, numOfMonomials, basisFunctions, weightingFunction);
	}
	linearSolver(Sum_bbT, Sum_bf);
}


template<class T>
void WeightedLeastSquares_LCS<T>::findCoefficients_LCS1D( )
{
	int numOfSamplePoints 			= this->m_scatteredData->get_numOfNodes();
	T* samplePoints 			= this->m_scatteredData->get_domainComponent(0);
	T* sampleValues 			= this->m_scatteredData->get_fieldComponent(0);

	T* stationaryPoints			= this->m_stationaryPoints->get_component(0);
	int numOfStationaryPoints		= this->get_numOfStationaryPoints();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	LinearSolver<T>& linearSolver		= *(this->m_linearSolver);
	linearSolver.set_numberOfColumns(numOfMonomials);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);


	if (m_computingMode)
	{
		this->m_coefficients.resize(numOfStationaryPoints*numOfMonomials);
		cudaVersion_findCoefficientsLCS<T>(weightingFunction.get_span(), numOfStationaryPoints, stationaryPoints, numOfSamplePoints, samplePoints, sampleValues, numOfMonomials, this->m_coefficients.data());
	}
	else
	{
		this->m_coefficients.reserve(numOfStationaryPoints*numOfMonomials);
		T Sum_bf[numOfMonomials];
		for (int stationaryPointIndex = 0; stationaryPointIndex < get_numOfStationaryPoints(); stationaryPointIndex++)
		{
			std::fill(Sum_bf, Sum_bf+numOfMonomials, T(0));
			findCoefficientsAtStationaryPoint_LCS<T>(stationaryPoints[stationaryPointIndex], numOfSamplePoints, samplePoints, sampleValues, numOfMonomials, linearSolver, basisFunctions, weightingFunction, Sum_bf);

			for (int index=0; index<numOfMonomials;index++) this->m_coefficients.push_back(Sum_bf[index]);
		}
	}
}



template<class T>
void WeightedLeastSquares_LCS<T>::find_coefficients()
{
	switch ( this->m_scatteredData->get_numOfDomainDimensions() )
	{
		case 1:	{ findCoefficients_LCS1D( ); break; }
		case 2:	{ findCoefficients_LCS2D( ); break; }
		case 3:	{ findCoefficients_LCS3D( ); break; }
	}
}


template<class T>
T evaluateSinglePoint_LCS(T* stationaryPoints, T evaluationPoint, T* coefficients_local, int numOfStationaryPoints, int numOfMonomials, BasisFunctions<T>& basisFunctions, WeightingFunction<T>& weightingFunction)
{
	T weightPerStationaryPoint[numOfStationaryPoints], weight_sum(0);
	for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPoints; stationaryPointIndex++)
	{
		weightPerStationaryPoint[stationaryPointIndex] = weightingFunction(fabs(evaluationPoint - stationaryPoints[stationaryPointIndex]));
		weight_sum += weightPerStationaryPoint[stationaryPointIndex];
	}

	T evaluation(0), b[numOfMonomials];
	for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPoints; stationaryPointIndex++)
	{
		basisFunctions(b, evaluationPoint -stationaryPoints[stationaryPointIndex]);
		evaluation += weightPerStationaryPoint[stationaryPointIndex]/weight_sum*dotProductOfVectors(b, b+numOfMonomials, coefficients_local+stationaryPointIndex*numOfMonomials);
	}

	return evaluation;
}


template<class T>
void WeightedLeastSquares_LCS<T>::evaluationMethod()
{
	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	T* stationaryPoints			= this->m_stationaryPoints->get_component(0);
	T* coefficients_local			= this->m_coefficients.data();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);
	int numOfStationaryPoints		= get_numOfStationaryPoints();


	if (m_computingMode) cudaVersion_evaluationMethodLCS(numOfEvaluationPoints, evaluationPoints, weightingFunction.get_span(), numOfStationaryPoints, stationaryPoints, numOfMonomials, coefficients_local, evaluations);
	else
		for (int evaluationPointIndex = 0; evaluationPointIndex < numOfEvaluationPoints; evaluationPointIndex++) evaluations[evaluationPointIndex] = evaluateSinglePoint_LCS<T>(stationaryPoints, evaluationPoints[evaluationPointIndex], coefficients_local, numOfStationaryPoints, numOfMonomials, basisFunctions, weightingFunction);

}




template<class T>
void WeightedLeastSquares_LCS<T>::derivativesEvaluationMethod()
{
/*
	int numOfSamplePoints 			= this->m_scatteredData->get_numOfNodes();
	T* samplePoints_x 			= this->m_scatteredData->get_domainComponent(0);
	T* samplePoints_y 			= this->m_scatteredData->get_domainComponent(1);
	T* sampleValues_x 			= this->m_scatteredData->get_fieldComponent(0);
	T* sampleValues_y 			= this->m_scatteredData->get_fieldComponent(1);

	T* stationaryPoints_x			= this->m_stationaryPoints->get_component(0);
	T* stationaryPoints_y			= this->m_stationaryPoints->get_component(1);
	int numOfStationaryPoints		= this->get_numOfStationaryPoints();
	int numOfMonomials 			= this->m_basisFunctions->get_numOfMonomials();
	int polDegree				= this->m_basisFunctions->get_degree();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	LinearSolver<T>& linearSolver		= *(this->m_linearSolver);
	linearSolver.set_numberOfColumns(numOfMonomials);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);

	int numOfEvaluationPoints		= this->m_evaluationData->get_numOfNodes();
	int numOfMonomials 			= this->m_basisFunctions->get_degree()+1;
	this->m_evaluationData->set_fieldComponent(0);
	T* evaluations				= this->m_evaluationData->get_fieldComponent(0);
	T* evaluationPoints 			= this->m_evaluationData->get_domainComponent(0);
	T* stationaryPoints			= this->m_stationaryPoints->get_component(0);
	T* coefficients_local			= this->m_coefficients.data();
	BasisFunctions<T>& basisFunctions	= *(this->m_basisFunctions);
	WeightingFunction<T>& weightingFunction = *(m_weightingFunction);
	int numOfStationaryPoints		= get_numOfStationaryPoints();
*/


//	this->m_coefficients.resize(2*numOfStationaryPoints*numOfMonomials);

//	this->m_calculationClocks = -clock();

//	cpuVersion_findCoefficients_LCS2D<T>(numOfStationaryPoints, stationaryPoints_x, stationaryPoints_y, numOfSamplePoints, samplePoints_x, samplePoints_y, sampleValues_x, sampleValues_y, numOfMonomials, this->m_coefficients.data(), linearSolver, basisFunctions, weightingFunction, 1 ); break; }


}
