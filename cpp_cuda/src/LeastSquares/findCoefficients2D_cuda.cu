#include<iostream>
#include<fstream>
#include"weightedLeastSquares_cuda.h"
#include"gpuDeviceProperties.h"
#include"WeightingFunction.h"
#include<typeinfo>
#include"vectorPlotOutStreamer.h"
#include"conjugateGradient.h"
#include"fileFunctions.h"
#include"myBLAS.h"
#include<cuda.h>
#include<cuda_runtime.h>
#include<cuda_runtime_api.h>
#include<cudaProfiler.h>
#include<cublas_v2.h>
#include"findCoefficients2D_cuda_kernels.h"
#include<math.h>
#include<cmath>
#include"cuBLAS_wrappers.h"


constexpr cudaMemcpyKind H2D = cudaMemcpyKind::cudaMemcpyHostToDevice;
constexpr cudaMemcpyKind D2H = cudaMemcpyKind::cudaMemcpyDeviceToHost;


// !!!! results to be reassured !!!! ===================================== optimal design, a given size of thread-block is processing iteratively in tiles the sample points, the block contains the data for all the monomials in order to be able to calculate the partial bbT which corresponds to the current sample points, minimized data transferring by calculating bbT and bf in the same kernel launch, several stationary points launched in a single kernel which can be controlled as a parameter =====================================================
template<class T>
void cudaVersion_findCoefficientsLCS2D(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID )
{

	weightingFunction.set_numOfVariableData(0);
	cudaSetDevice( cudaDeviceID );

	//========================== size of arrays to be processed
	size_t sizeOfCoefficientsArray 		= 2*sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;
	size_t sizeOfPowersArray 		= sizeof(int)*numOfMonomials;
	size_t sizeOfbbTArray			= sizeof(T)*numOfMonomials*numOfMonomials;
	size_t sizeOfbfArray			= sizeof(T)*numOfMonomials;

	cudaError_t cudaStatus; 
	memset( coefficients, T(0), sizeOfCoefficientsArray );


	//========================= find monomials degree
	int *powers_x, *powers_y; int count(0);
	cudaStatus = cudaHostAlloc( (void**) &powers_x, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &powers_y, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	for (int i = 0; i<=polDegree; i++) for (int j = 0; j<=polDegree-i; j++) { powers_x[count] = j; powers_y[count] = i; count++; }


//	for (int i = 256+128; i<numOfSamplePoints; i++) { std::cout << sampleValues_y[i] << std::endl; }

	//========================== page lock host arrays
	cudaStatus = cudaHostRegister (stationaryPoints_x,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_y,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


// ==============================  Set the computational resources ===================================


	//========================= analyse hardware and available computing resources
//	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
//	dev_count = get_AllDevicesProperties(dev_prop);

//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);
//	plotGPUProperties(dev_prop[0]);

//	int totalAvailableThreads(0);
//	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
//		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;

//	std::cout << "totanl num of threads " << totalAvailableThreads << std::endl;
//	std::cout << "num of multiprocessors " << dev_prop[0].multiProcessorCount << std::endl;

	size_t maxSizeOfThreadBlock(1024);


	//========================== computational parameters per stationary point
	size_t numOfOperations	( numOfMonomials );
	size_t numOfDataUnits	( numOfSamplePoints );
	size_t numOfStreams	( 3 );


	if ( typeid(T) == typeid( double ) ) cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);



	size_t numOfSamplePointsPerBlock = highestPowerof2(1024/numOfMonomials);
//	std::min( pow_i(2, 8 -numOfMonomials/3 -1), numOfSamplePoints );
	switch ( polDegree ) {
		case 0: { numOfSamplePointsPerBlock = 128; break;}
		case 1: { numOfSamplePointsPerBlock = 64; break;}
		case 2: { numOfSamplePointsPerBlock = 32; break;}
		case 3: { numOfSamplePointsPerBlock = 16; break;}
		default: numOfSamplePointsPerBlock = std::min( numOfSamplePointsPerBlock, size_t(16));
	}

	// problem data analysis per stationary point
	size_t 	numOfStationaryPointsPerKernel 	= std::min( 1*1024, numOfStationaryPoints );

	size_t 	numOfBlocks 			= numOfSamplePoints /numOfSamplePointsPerBlock;
	dim3 	gridOfBlocks 			( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	blockSize 			( numOfSamplePointsPerBlock, numOfMonomials, 1); 


	size_t	dataArraySizePerKernel		= numOfSamplePoints*numOfMonomials*numOfStationaryPointsPerKernel;
	size_t	sizeOfBlockDataTile		= numOfSamplePointsPerBlock*sizeof(T);
	// powers_x +powers_y +max(b +sampleValues_x +sampleValues_y +bfx +bfy,   b +bbT)
	size_t 	sharedMemorySize 		= 2*numOfMonomials*sizeOfBlockDataTile +sizeOfbbTArray;

	std::cout << "\nCUDA execution settings:"
		<< "\nNumber of stationary points in kernel \t\t\t" 		<< numOfStationaryPointsPerKernel
		<< "\nNumber of sample points per block \t\t\t" 		<< numOfSamplePointsPerBlock
		<< "\nNumber of sample points tiles in kernel \t\t"		<< numOfBlocks
		<< "\nShared memory per block \t\t\t\t"				<< sharedMemorySize << "b"
		<< std::endl;


	//========================== create processing streams
//	cudaStream_t stream0, stream1, stream2;
//	cudaStreamCreate(&stream0);// cudaStreamCreate(&stream1); cudaStreamCreate(&stream2);


	//========================= gpu memory pointers
	T *gpu_coefficients, *gpu_stationaryPoints_x, *gpu_stationaryPoints_y, *gpu_samplePoints_x, *gpu_samplePoints_y, *gpu_sampleValues_x, *gpu_sampleValues_y, *gpu_bbT, *gpu_bfx, *gpu_bfy;
	int *gpu_powers_x, *gpu_powers_y;



	//========================== allocate gpu memory
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_x, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_y, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_x, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_y, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_x, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_y, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_x, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_y, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bbT, 			sizeOfbbTArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bfx,			sizeOfbfArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bfy,			sizeOfbfArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;



//================================== Launch the process ========================================================


	T *Sum_bbT, *Sum_bfx, *Sum_bfy;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bbT, sizeOfbbTArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &Sum_bfx, sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bfy, sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0); clock_t calculationClocks; T distance;

	size_t totalGlobalMemory = 2*sizeOfPowersArray +4*sizeOfSamplePointsArray +2*sizeOfStationaryPointsArray +sizeOfbbTArray*numOfStationaryPointsPerKernel +2*sizeOfbfArray*numOfStationaryPointsPerKernel +2*dataArraySizePerKernel*sizeof(T);
	std::cout << "Global memory required: \t\t\t\t" << totalGlobalMemory/1024 << "kb" << std::endl;


//	std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); calculationClocks = -clock();



	cudaStatus = cudaMemcpy( gpu_powers_x, 			powers_x, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_y, 			powers_y, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_x,		samplePoints_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_y, 		samplePoints_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_x, 		sampleValues_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_y, 		sampleValues_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_x, 	stationaryPoints_x, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_y, 	stationaryPoints_y, 	sizeOfStationaryPointsArray, 	H2D );


	typedef void (*kernelFun2D) ( int*, int*, T*, T*, size_t, T*, T*, T*, T*, WeightingFunction<T>, T*, T*, T*);

	kernelFun2D launchedFun( findCoefficientsLCS2D_cudaKernel );
	if (numOfMonomials == 1) launchedFun = findCoefficientsLCS2D_cudaKernel_0deg;
 


	size_t stationaryPointShift(0), coefficientsShift(0);
	for (int stationaryPoint = 0; stationaryPoint < numOfStationaryPoints; stationaryPoint += numOfStationaryPointsPerKernel )
	{

		if ( stationaryPoint +numOfStationaryPointsPerKernel > numOfStationaryPoints ) {
			numOfStationaryPointsPerKernel 	= numOfStationaryPoints -stationaryPoint;
			gridOfBlocks 			= dim3( numOfBlocks, numOfStationaryPointsPerKernel, 1 ); }


		//cuProfilerStart();

		// ===================================== gpu processing part ==================================
		launchedFun <<< dim3( numOfStationaryPointsPerKernel, 1, 1 ), dim3( numOfSamplePointsPerBlock, numOfMonomials, 1 ), sharedMemorySize >>>
		( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, gpu_sampleValues_x, gpu_sampleValues_y, weightingFunction, gpu_bbT, gpu_bfx, gpu_bfy);


		cudaDeviceSynchronize();
		cudaMemcpy( Sum_bbT, gpu_bbT, sizeOfbbTArray*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bfx, gpu_bfx, sizeOfbfArray*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bfy, gpu_bfy, sizeOfbfArray*numOfStationaryPointsPerKernel, D2H);
		cudaDeviceSynchronize();

		// ========================================= cpu post processing part ================================
		int counter(0);
		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			for ( size_t samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++ ) {
				distance = distance2D(stationaryPoints_x[stationaryPointIndex +stationaryPointShift ], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex +stationaryPointShift ], samplePoints_y[samplePointIndex]);
				if ( distance < weightingFunction.get_span() ) counter++; if ( counter >= numOfMonomials ) break;
			}
			if ( counter < numOfMonomials ) {
				std::fill(coefficients +stationaryPointIndex*2*numOfMonomials 	  +coefficientsShift,
					  coefficients +(stationaryPointIndex+1)*2*numOfMonomials +coefficientsShift,
					  T(0) );
				continue; }

			numOfEffectiveStationaryPoints++;

			linearSolver.set_givenMatrix(0);
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bfx +stationaryPointIndex*numOfMonomials, 
					coefficients +stationaryPointIndex*2*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bfy +stationaryPointIndex*numOfMonomials,
					coefficients +(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift );

			totalDisplacement_L2norm += 	 coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ] 
							*coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ];
		}

		stationaryPointShift += numOfStationaryPointsPerKernel;
		coefficientsShift = stationaryPointShift*2*numOfMonomials;

//		std::cout << "\r" << (stationaryPointShift)*nodesFraction << "% "; std::cout.flush();

	}


//	std::cout << "\r" << "100.00%\n\n";
	std::cout << "CUDA GPU process" << std::endl;
	std::cout << "Number of effective stationary points: \t\t" 	<< numOfEffectiveStationaryPoints << std::endl;
	std::cout << "Number of outlier stationary points: \t\t" 	<< numOfStationaryPoints -numOfEffectiveStationaryPoints << std::endl;
	std::cout << "L2 norm of the evaluated 2D displacement field: " << sqrt(totalDisplacement_L2norm) << std::endl << std::endl;


// ============================ Free Allocated resources  ==============================================

//	cudaStreamDestroy ( stream0 );// cudaStreamDestroy ( stream1 ); cudaStreamDestroy ( stream2 );

	cudaFreeHost ( powers_x );
	cudaFreeHost ( powers_y );
	cudaFreeHost ( Sum_bbT );
	cudaFreeHost ( Sum_bfx );
	cudaFreeHost ( Sum_bfy );

	cudaHostUnregister (stationaryPoints_x);
	cudaHostUnregister (stationaryPoints_y);
	cudaHostUnregister (samplePoints_x);
	cudaHostUnregister (samplePoints_y);
	cudaHostUnregister (sampleValues_x);
	cudaHostUnregister (sampleValues_y);

	cudaFree ( gpu_powers_x );		cudaFree ( gpu_powers_y );
	cudaFree ( gpu_stationaryPoints_x );	cudaFree ( gpu_stationaryPoints_y );
	cudaFree ( gpu_samplePoints_x ); 	cudaFree ( gpu_samplePoints_y );
	cudaFree ( gpu_sampleValues_x );	cudaFree ( gpu_sampleValues_y );
	cudaFree ( gpu_bbT );
	cudaFree ( gpu_bfx );
	cudaFree ( gpu_bfy );

//	cudaDeviceReset();

}

template void cudaVersion_findCoefficientsLCS2D(int numOfStationaryPoints, float* stationaryPoints_x, float* stationaryPoints_y, int numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, int numOfMonomials, int polDegree, float* coefficients, LinearSolver<float>& solver, WeightingFunction<float>& wFun, size_t deviceId);
template void cudaVersion_findCoefficientsLCS2D(int numOfStationaryPoints, double* stationaryPoints_x, double* stationaryPoints_y, int numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, int numOfMonomials, int polDegree, double* coefficients, LinearSolver<double>& solver, WeightingFunction<double>& wFun, size_t deviceId);





// !! assured correct resutls !!============================================== algorithm design for massive kernel launches, each processing a parametrized tile of data  ===============================================
template<class T>
void cudaVersion_findCoefficientsLCS2D2(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID )
{
	weightingFunction.set_numOfVariableData(0);
	cudaSetDevice( cudaDeviceID );

	//========================== size of arrays to be processed
	size_t sizeOfCoefficientsArray 		= 2*sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;
	size_t sizeOfPowersArray 		= sizeof(int)*numOfMonomials;
	size_t sizeOfbbTArray			= sizeof(T)*numOfMonomials*numOfMonomials;
	size_t sizeOfbfArray			= sizeof(T)*numOfMonomials;

	cudaError_t cudaStatus; 
	memset( coefficients, T(0), sizeOfCoefficientsArray );


	//========================= find monomials degree
	int *powers_x, *powers_y; int count(0);
	cudaStatus = cudaHostAlloc( (void**) &powers_x, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &powers_y, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	for (int i = 0; i<=polDegree; i++) for (int j = 0; j<=polDegree-i; j++) { powers_x[count] = j; powers_y[count] = i; count++; }


//	for (int i = 256+128; i<numOfSamplePoints; i++) { std::cout << sampleValues_y[i] << std::endl; }

	//========================== page lock host arrays
	cudaStatus = cudaHostRegister (stationaryPoints_x,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_y,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


// ==============================  Set the computational resources ===================================


	//========================= analyse hardware and available computing resources
//	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
//	dev_count = get_AllDevicesProperties(dev_prop);

//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);
//	plotGPUProperties(dev_prop[0]);

//	int totalAvailableThreads(0);
//	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
//		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;

//	std::cout << "totanl num of threads " << totalAvailableThreads << std::endl;
//	std::cout << "num of multiprocessors " << dev_prop[0].multiProcessorCount << std::endl;

	size_t maxSizeOfThreadBlock(1024);


	//========================== computational parameters per stationary point
	size_t numOfOperations	( numOfMonomials );
	size_t numOfDataUnits	( numOfSamplePoints );
	size_t numOfStreams	( 3 );


	if ( typeid(T) == typeid( double ) ) cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);


	// problem data analysis per stationary point
	size_t 	numOfStationaryPointsPerKernel 	= min( 1024, numOfStationaryPoints );
	size_t	numOfSamplePointsPerKernel	= max(16, min( pow_i(2, 8 -numOfMonomials/3 -1), numOfSamplePoints ));
	size_t	numOfSamplePointsPerBlock	= min( numOfSamplePointsPerKernel, min(maxSizeOfThreadBlock, size_t(highestPowerof2( maxSizeOfThreadBlock/numOfMonomials )) ) );
	size_t	numOfIters			= numOfSamplePoints /numOfSamplePointsPerKernel;
	size_t 	numOfBlocks 			= numOfSamplePointsPerKernel /numOfSamplePointsPerBlock;
	dim3 	gridOfBlocks 			( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	blockSize 			( numOfSamplePointsPerBlock, numOfMonomials, 1); 


	size_t	numOfResidualBlocks		= (numOfSamplePoints -numOfSamplePointsPerBlock*numOfBlocks*numOfIters) /numOfSamplePointsPerBlock;
	size_t	numOfResidualSamplePoints	= numOfSamplePoints -numOfSamplePointsPerBlock*(numOfBlocks*numOfIters +numOfResidualBlocks);
	bool	divGridFlag			= numOfResidualBlocks;
	bool	divBlockFlag			= numOfResidualSamplePoints;
	dim3 	divGridOfBlocks 		( numOfResidualBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	divBlockSize			( numOfResidualSamplePoints, numOfMonomials, 1);


	size_t	sizeOfBlockDataTile		= numOfSamplePointsPerBlock*sizeof(T);
	// powers_x +powers_y +max(b +sampleValues_x +sampleValues_y +bfx +bfy,   b +bbT)
//	size_t 	sharedMemorySize 		= 2*sizeOfPowersArray +sizeOfBlockDataTile +numOfMonomials*sizeOfBlockDataTile +max( 2*sizeOfBlockDataTile +2*4*sizeOfbfArray, sizeOfbbTArray);
	size_t 	sharedMemorySize 		= 2*numOfMonomials*sizeOfBlockDataTile +sizeOfbbTArray;
	size_t	sizeOfResidualDataTile		= numOfResidualSamplePoints*sizeof(T);
//	size_t 	sharedMemorySize_divBlc		= 2*sizeOfPowersArray +sizeOfBlockDataTile +numOfMonomials*sizeOfResidualDataTile +max( 2*sizeOfResidualDataTile +2*4*sizeOfbfArray, sizeOfbbTArray);
	size_t 	sharedMemorySize_divBlc		= 2*numOfMonomials*sizeOfResidualDataTile +sizeOfbbTArray;

/*
	std::cout << "\nCUDA execution settings:"
		<< "\nNumber of stationary points in kernel \t\t\t" 		<< numOfStationaryPointsPerKernel
		<< "\nNumber of sample points in kernel \t\t\t" 		<< numOfSamplePointsPerKernel
		<< "\nNumber of sample points per block \t\t\t" 		<< numOfSamplePointsPerBlock
		<< "\nNumber of kernel iterations per stationary point \t" 	<< numOfIters
		<< "\nNumber of sample point blocks \t\t\t\t"			<< numOfBlocks
		<< "\nNumber of residual sample point blocks \t\t\t"		<< numOfResidualBlocks
		<< "\nNumber of residual sample points \t\t\t" 			<< numOfResidualSamplePoints
		<< "\nShared memory per block \t\t\t\t"				<< sharedMemorySize << "b"
		<< "\nShared memory in residual sample points block \t\t"	<< sharedMemorySize_divBlc << "b"
		<< std::endl;
*/

	//========================== create processing streams
//	cudaStream_t stream0, stream1, stream2;
//	cudaStreamCreate(&stream0);// cudaStreamCreate(&stream1); cudaStreamCreate(&stream2);


	//========================= gpu memory pointers
	T *gpu_coefficients, *gpu_stationaryPoints_x, *gpu_stationaryPoints_y, 
	  *gpu_samplePoints_x, *gpu_samplePoints_y, *gpu_sampleValues_x, *gpu_sampleValues_y,
	  *gpu_bbT, *gpu_bf;
	int *gpu_powers_x, *gpu_powers_y;


	//========================== allocate gpu memory
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_x, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_y, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_x, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_y, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bbT, 			sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bf,			2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


//================================== Launch the process ========================================================


	T *Sum_bbT, *Sum_bf;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bbT, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &Sum_bf, 2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0); clock_t calculationClocks; T distance;

	size_t totalGlobalMemory = 2*sizeOfPowersArray +4*sizeOfSamplePointsArray +2*sizeOfStationaryPointsArray +sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel +2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel;

//	std::cout << "Global memory required: \t\t\t\t" << totalGlobalMemory/1024 << "kb" << std::endl;

//	std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); calculationClocks = -clock();

	cudaStatus = cudaMemcpy( gpu_powers_x, 			powers_x, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_y, 			powers_y, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_x,		samplePoints_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_y, 		samplePoints_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_x, 		sampleValues_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_y, 		sampleValues_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_x, 	stationaryPoints_x, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_y, 	stationaryPoints_y, 	sizeOfStationaryPointsArray, 	H2D );


	typedef void (*kernelFun2D) ( int*, int*, T*, T*, T*, T*, T*, T*, size_t, WeightingFunction<T>, T*, T*);

	kernelFun2D launchedFun(findCoefficientsLCS2D_cudaKernel2);
	if (numOfMonomials == 1) launchedFun = findCoefficientsLCS2D_cudaKernel_0deg2;



	size_t stationaryPointShift(0), coefficientsShift(0);
	for (int stationaryPoint = 0; stationaryPoint < numOfStationaryPoints; stationaryPoint += numOfStationaryPointsPerKernel )
	{

		if ( stationaryPoint +numOfStationaryPointsPerKernel > numOfStationaryPoints ) {
 			numOfStationaryPointsPerKernel 	= numOfStationaryPoints -stationaryPoint;
			gridOfBlocks 			= dim3( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
			divGridOfBlocks 		= dim3( numOfResidualBlocks, numOfStationaryPointsPerKernel, 1 ); }


		// ====================================== gpu processing part =========================================
		cudaMemset ( gpu_bf, 0, 2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel);
		cudaMemset ( gpu_bbT, 0, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel);
		cudaDeviceSynchronize();



		for ( int iteration = 0; iteration < numOfIters; iteration++ )
		{

			if (stationaryPoint == 0 && iteration == 0 ) cuProfilerStart();
			size_t dataShift = iteration*numOfBlocks*numOfSamplePointsPerBlock;
			launchedFun <<< gridOfBlocks, blockSize, sharedMemorySize>>>
			( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, 
			gpu_samplePoints_x +dataShift, gpu_samplePoints_y +dataShift, gpu_sampleValues_x +dataShift, gpu_sampleValues_y +dataShift, numOfBlocks, weightingFunction, gpu_bbT, gpu_bf);
			cudaDeviceSynchronize();

			if (stationaryPoint == 0 && iteration == 0) cuProfilerStop();
		}



		if ( divGridFlag ) {
			size_t dataShift = numOfIters*numOfBlocks*numOfSamplePointsPerBlock;
			launchedFun <<< divGridOfBlocks, blockSize, sharedMemorySize>>>
			( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint,
			gpu_samplePoints_x +dataShift, gpu_samplePoints_y +dataShift, gpu_sampleValues_x +dataShift, gpu_sampleValues_y +dataShift, numOfBlocks, weightingFunction, gpu_bbT, gpu_bf);}


		if ( divBlockFlag ) {
			size_t dataShift = (numOfIters*numOfBlocks +numOfResidualBlocks)*numOfSamplePointsPerBlock;
			launchedFun <<< dim3(1, numOfStationaryPointsPerKernel, 1), divBlockSize, sharedMemorySize_divBlc>>>
			( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, 
			gpu_samplePoints_x +dataShift, gpu_samplePoints_y +dataShift, gpu_sampleValues_x +dataShift, gpu_sampleValues_y +dataShift, numOfBlocks, weightingFunction, gpu_bbT, gpu_bf);}

		cudaDeviceSynchronize();

		cudaMemcpy( Sum_bbT, gpu_bbT, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bf, gpu_bf, 2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel, D2H);
		cudaDeviceSynchronize();



		// ====================================== cpu post processing part ====================================
		int counter(0);
		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			counter = 0;
			for ( size_t samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++ ) {
				distance = distance2D(stationaryPoints_x[stationaryPointIndex +stationaryPointShift ], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex +stationaryPointShift ], samplePoints_y[samplePointIndex]);
				if ( distance < weightingFunction.get_span() ) counter++; if ( counter >= numOfMonomials ) break;
			}
			if ( counter < numOfMonomials ) {
				std::fill(coefficients +stationaryPointIndex*2*numOfMonomials 	  +coefficientsShift,
					  coefficients +(stationaryPointIndex+1)*2*numOfMonomials +coefficientsShift,
					  T(0) );
				continue; }

			numOfEffectiveStationaryPoints++;

			for ( int block = 1; block < numOfBlocks; block++ ) {
				for (int monomial = 0; monomial < numOfMonomials*numOfMonomials; monomial++)
					Sum_bbT[ monomial +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks ] += Sum_bbT[ monomial +block*numOfMonomials*numOfMonomials +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks ];
				for (int monomial = 0; monomial < 2*numOfMonomials; monomial++)
					Sum_bf[ monomial +stationaryPointIndex*numOfMonomials*2*numOfBlocks ] += Sum_bf[ monomial +block*2*numOfMonomials +stationaryPointIndex*numOfMonomials*2*numOfBlocks ];
			}


			linearSolver.set_givenMatrix(0);
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks,
					Sum_bf  +stationaryPointIndex*numOfMonomials*2*numOfBlocks, 
					coefficients +stationaryPointIndex*2*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT+stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks, 
					Sum_bf +numOfMonomials+stationaryPointIndex*numOfMonomials*2*numOfBlocks,
					coefficients +(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift );

			totalDisplacement_L2norm += 	 coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ] 
							*coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ];

		}

		stationaryPointShift += numOfStationaryPointsPerKernel;
		coefficientsShift = stationaryPointShift*2*numOfMonomials;

//		std::cout << "\r" << (stationaryPointShift)*nodesFraction << "% "; std::cout.flush();

	}

/*
//	std::cout << "\r" << "100.00%\n\n";
	std::cout << "CUDA GPU process" << std::endl;
	std::cout << "Number of effective stationary points: \t\t" 	<< numOfEffectiveStationaryPoints << std::endl;
	std::cout << "Number of outlier stationary points: \t\t" 	<< numOfStationaryPoints -numOfEffectiveStationaryPoints << std::endl;
	std::cout << "L2 norm of the evaluated 2D displacement field: " << sqrt(totalDisplacement_L2norm) << std::endl << std::endl;
*/

// ============================ Free Allocated resources  ==============================================

//	cudaStreamDestroy ( stream0 );// cudaStreamDestroy ( stream1 ); cudaStreamDestroy ( stream2 );

	cudaFreeHost ( powers_x );
	cudaFreeHost ( powers_y );
	cudaFreeHost ( Sum_bbT );
	cudaFreeHost ( Sum_bf );

	cudaHostUnregister (stationaryPoints_x);
	cudaHostUnregister (stationaryPoints_y);
	cudaHostUnregister (samplePoints_x);
	cudaHostUnregister (samplePoints_y);
	cudaHostUnregister (sampleValues_x);
	cudaHostUnregister (sampleValues_y);

	cudaFree ( gpu_powers_x );		cudaFree ( gpu_powers_y );
	cudaFree ( gpu_stationaryPoints_x );	cudaFree ( gpu_stationaryPoints_y );
	cudaFree ( gpu_samplePoints_x ); 	cudaFree ( gpu_samplePoints_y );
	cudaFree ( gpu_sampleValues_x );	cudaFree ( gpu_sampleValues_y );
	cudaFree ( gpu_bbT );
	cudaFree ( gpu_bf );

//	cudaDeviceReset();
}

template void cudaVersion_findCoefficientsLCS2D2(int numOfStationaryPoints, float* stationaryPoints_x, float* stationaryPoints_y, int numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, int numOfMonomials, int polDegree, float* coefficients, LinearSolver<float>& solver, WeightingFunction<float>& wFun, size_t deviceId);
template void cudaVersion_findCoefficientsLCS2D2(int numOfStationaryPoints, double* stationaryPoints_x, double* stationaryPoints_y, int numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, int numOfMonomials, int polDegree, double* coefficients, LinearSolver<double>& solver, WeightingFunction<double>& wFun, size_t deviceId);


/*
// !! assured correct resutls !!============================================== algorithm design for massive kernel launches, each processing a parametrized tile of data  ===============================================
template<class T>
void cudaVersion_findCoefficientsLCS2DConstantMemory2(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID )
{
	weightingFunction.set_numOfVariableData(0);
	cudaSetDevice( cudaDeviceID );

	//========================== size of arrays to be processed
	size_t sizeOfCoefficientsArray 		= 2*sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;
	size_t sizeOfPowersArray 		= sizeof(int)*numOfMonomials;
	size_t sizeOfbbTArray			= sizeof(T)*numOfMonomials*numOfMonomials;
	size_t sizeOfbfArray			= sizeof(T)*numOfMonomials;

	cudaError_t cudaStatus; 
	memset( coefficients, T(0), sizeOfCoefficientsArray );


	//========================= find monomials degree
	int *powers_x, *powers_y; int count(0);
	cudaStatus = cudaHostAlloc( (void**) &powers_x, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &powers_y, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	for (int i = 0; i<=polDegree; i++) for (int j = 0; j<=polDegree-i; j++) { powers_x[count] = j; powers_y[count] = i; count++; }


//	for (int i = 256+128; i<numOfSamplePoints; i++) { std::cout << sampleValues_y[i] << std::endl; }

	//========================== page lock host arrays
	cudaStatus = cudaHostRegister (stationaryPoints_x,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_y,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


// ==============================  Set the computational resources ===================================


	//========================= analyse hardware and available computing resources
//	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
//	dev_count = get_AllDevicesProperties(dev_prop);

//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);
//	plotGPUProperties(dev_prop[0]);

//	int totalAvailableThreads(0);
//	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
//		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;

//	std::cout << "totanl num of threads " << totalAvailableThreads << std::endl;
//	std::cout << "num of multiprocessors " << dev_prop[0].multiProcessorCount << std::endl;

	size_t maxSizeOfThreadBlock(1024);


	//========================== computational parameters per stationary point
	size_t numOfOperations	( numOfMonomials );
	size_t numOfDataUnits	( numOfSamplePoints );
	size_t numOfStreams	( 3 );


	if ( typeid(T) == typeid( double ) ) cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);


	// problem data analysis per stationary point
	size_t 	numOfStationaryPointsPerKernel 	= min( 1024, numOfStationaryPoints );
	size_t	numOfSamplePointsPerKernel	= max(16, min( pow_i(2, 8 -numOfMonomials/3 -1), numOfSamplePoints ));
	size_t	numOfSamplePointsPerBlock	= min( numOfSamplePointsPerKernel, min(maxSizeOfThreadBlock, size_t(highestPowerof2( maxSizeOfThreadBlock/numOfMonomials )) ) );
	size_t	numOfIters			= numOfSamplePoints /numOfSamplePointsPerKernel;
	size_t 	numOfBlocks 			= numOfSamplePointsPerKernel /numOfSamplePointsPerBlock;
	dim3 	gridOfBlocks 			( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	blockSize 			( numOfSamplePointsPerBlock, numOfMonomials, 1); 


	size_t	numOfResidualBlocks		= (numOfSamplePoints -numOfSamplePointsPerBlock*numOfBlocks*numOfIters) /numOfSamplePointsPerBlock;
	size_t	numOfResidualSamplePoints	= numOfSamplePoints -numOfSamplePointsPerBlock*(numOfBlocks*numOfIters +numOfResidualBlocks);
	bool	divGridFlag			= numOfResidualBlocks;
	bool	divBlockFlag			= numOfResidualSamplePoints;
	dim3 	divGridOfBlocks 		( numOfResidualBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	divBlockSize			( numOfResidualSamplePoints, numOfMonomials, 1);


	size_t	sizeOfBlockDataTile		= numOfSamplePointsPerBlock*sizeof(T);
	// powers_x +powers_y +max(b +sampleValues_x +sampleValues_y +bfx +bfy,   b +bbT)
//	size_t 	sharedMemorySize 		= 2*sizeOfPowersArray +sizeOfBlockDataTile +numOfMonomials*sizeOfBlockDataTile +max( 2*sizeOfBlockDataTile +2*4*sizeOfbfArray, sizeOfbbTArray);
	size_t 	sharedMemorySize 		= 2*numOfMonomials*sizeOfBlockDataTile +sizeOfbbTArray;
	size_t	sizeOfResidualDataTile		= numOfResidualSamplePoints*sizeof(T);
//	size_t 	sharedMemorySize_divBlc		= 2*sizeOfPowersArray +sizeOfBlockDataTile +numOfMonomials*sizeOfResidualDataTile +max( 2*sizeOfResidualDataTile +2*4*sizeOfbfArray, sizeOfbbTArray);
	size_t 	sharedMemorySize_divBlc		= 2*numOfMonomials*sizeOfResidualDataTile +sizeOfbbTArray;


	std::cout << "\nCUDA execution settings:"
		<< "\nNumber of stationary points in kernel \t\t\t" 		<< numOfStationaryPointsPerKernel
		<< "\nNumber of sample points in kernel \t\t\t" 		<< numOfSamplePointsPerKernel
		<< "\nNumber of sample points per block \t\t\t" 		<< numOfSamplePointsPerBlock
		<< "\nNumber of kernel iterations per stationary point \t" 	<< numOfIters
		<< "\nNumber of sample point blocks \t\t\t\t"			<< numOfBlocks
		<< "\nNumber of residual sample point blocks \t\t\t"		<< numOfResidualBlocks
		<< "\nNumber of residual sample points \t\t\t" 			<< numOfResidualSamplePoints
		<< "\nShared memory per block \t\t\t\t"				<< sharedMemorySize << "b"
		<< "\nShared memory in residual sample points block \t\t"	<< sharedMemorySize_divBlc << "b"
		<< std::endl;


	//========================== create processing streams
//	cudaStream_t stream0, stream1, stream2;
//	cudaStreamCreate(&stream0);// cudaStreamCreate(&stream1); cudaStreamCreate(&stream2);


	//========================= gpu memory pointers
	T *gpu_coefficients, *gpu_stationaryPoints_x, *gpu_stationaryPoints_y, *gpu_bbT, *gpu_bf;
//	int *gpu_powers_x, *gpu_powers_y;

	//*gpu_samplePoints_x, *gpu_samplePoints_y, *gpu_sampleValues_x, *gpu_sampleValues_y,


	//========================== allocate gpu memory
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_x, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_y, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
//	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
//	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
//	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
//	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_x, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_y, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bbT, 			sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bf,			2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


//================================== Launch the process ========================================================


	T *Sum_bbT, *Sum_bf;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bbT, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &Sum_bf, 2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0); clock_t calculationClocks; T distance;

	size_t totalGlobalMemory = 2*sizeOfPowersArray +4*sizeOfSamplePointsArray +2*sizeOfStationaryPointsArray +sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel +2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel;

	std::cout << "Global memory required: \t\t\t\t" << totalGlobalMemory/1024 << "kb" << std::endl;

//	std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); calculationClocks = -clock();

	cudaMemcpyToSymbol( gpu_samplePoints_x,	samplePoints_x, sizeOfSamplePointsArray);
	cudaMemcpyToSymbol( gpu_samplePoints_y, samplePoints_y, sizeOfSamplePointsArray);
	cudaMemcpyToSymbol( gpu_sampleValues_x, sampleValues_x, sizeOfSamplePointsArray);
	cudaMemcpyToSymbol( gpu_sampleValues_y, sampleValues_y, sizeOfSamplePointsArray);
	
	cudaMemcpyToSymbol( gpu_powers_x, powers_x, sizeOfPowersArray );
	cudaMemcpyToSymbol( gpu_powers_y, powers_y, sizeOfPowersArray );


//	cudaStatus = cudaMemcpy( gpu_powers_x, 			powers_x, 		sizeOfPowersArray, 	 	H2D );
//	cudaStatus = cudaMemcpy( gpu_powers_y, 			powers_y, 		sizeOfPowersArray, 	 	H2D );
//	cudaStatus = cudaMemcpy( gpu_samplePoints_x,		samplePoints_x, 	sizeOfSamplePointsArray, 	H2D );
//	cudaStatus = cudaMemcpy( gpu_samplePoints_y, 		samplePoints_y, 	sizeOfSamplePointsArray, 	H2D );
//	cudaStatus = cudaMemcpy( gpu_sampleValues_x, 		sampleValues_x, 	sizeOfSamplePointsArray, 	H2D );
//	cudaStatus = cudaMemcpy( gpu_sampleValues_y, 		sampleValues_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_x, 	stationaryPoints_x, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_y, 	stationaryPoints_y, 	sizeOfStationaryPointsArray, 	H2D );


	typedef void (*kernelFun2D) ( T*, T*, size_t, size_t, WeightingFunction<T>, T*, T*);

	kernelFun2D launchedFun(findCoefficientsLCS2D_cudaKernelConstantMemory2);
	if (numOfMonomials == 1) launchedFun = findCoefficientsLCS2D_cudaKernel_0degConstantMemory2;



	size_t stationaryPointShift(0), coefficientsShift(0);
	for (int stationaryPoint = 0; stationaryPoint < numOfStationaryPoints; stationaryPoint += numOfStationaryPointsPerKernel )
	{

		if ( stationaryPoint +numOfStationaryPointsPerKernel > numOfStationaryPoints ) {
 			numOfStationaryPointsPerKernel 	= numOfStationaryPoints -stationaryPoint;
			gridOfBlocks 			= dim3( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
			divGridOfBlocks 		= dim3( numOfResidualBlocks, numOfStationaryPointsPerKernel, 1 ); }


		// ====================================== gpu processing part =========================================
		cudaMemset ( gpu_bf, 0, 2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel);
		cudaMemset ( gpu_bbT, 0, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel);
		cudaDeviceSynchronize();



		for ( int iteration = 0; iteration < numOfIters; iteration++ )
		{

			if (stationaryPoint == 0 && iteration == 0 ) cuProfilerStart();
			size_t dataShift = iteration*numOfBlocks*numOfSamplePointsPerBlock;
			launchedFun <<< gridOfBlocks, blockSize, sharedMemorySize>>>
			( gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, numOfBlocks, dataShift, weightingFunction, gpu_bbT, gpu_bf);
			cudaDeviceSynchronize();

			if (stationaryPoint == 0 && iteration == 0) cuProfilerStop();
		}



		if ( divGridFlag ) {
			size_t dataShift = numOfIters*numOfBlocks*numOfSamplePointsPerBlock;
			launchedFun <<< divGridOfBlocks, blockSize, sharedMemorySize>>>
			( gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, numOfBlocks, dataShift, weightingFunction, gpu_bbT, gpu_bf);}


		if ( divBlockFlag ) {
			size_t dataShift = (numOfIters*numOfBlocks +numOfResidualBlocks)*numOfSamplePointsPerBlock;
			launchedFun <<< dim3(1, numOfStationaryPointsPerKernel, 1), divBlockSize, sharedMemorySize_divBlc>>>
			( gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, numOfBlocks, dataShift, weightingFunction, gpu_bbT, gpu_bf);}

		cudaDeviceSynchronize();

		cudaMemcpy( Sum_bbT, gpu_bbT, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bf, gpu_bf, 2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel, D2H);
		cudaDeviceSynchronize();



		// ====================================== cpu post processing part ====================================
		int counter(0);
		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			counter = 0;
			for ( size_t samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++ ) {
				distance = distance2D(stationaryPoints_x[stationaryPointIndex +stationaryPointShift ], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex +stationaryPointShift ], samplePoints_y[samplePointIndex]);
				if ( distance < weightingFunction.get_span() ) counter++; if ( counter >= numOfMonomials ) break;
			}
			if ( counter < numOfMonomials ) {
				std::fill(coefficients +stationaryPointIndex*2*numOfMonomials 	  +coefficientsShift,
					  coefficients +(stationaryPointIndex+1)*2*numOfMonomials +coefficientsShift,
					  T(0) );
				continue; }

			numOfEffectiveStationaryPoints++;

			for ( int block = 1; block < numOfBlocks; block++ ) {
				for (int monomial = 0; monomial < numOfMonomials*numOfMonomials; monomial++)
					Sum_bbT[ monomial +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks ] += Sum_bbT[ monomial +block*numOfMonomials*numOfMonomials +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks ];
				for (int monomial = 0; monomial < 2*numOfMonomials; monomial++)
					Sum_bf[ monomial +stationaryPointIndex*numOfMonomials*2*numOfBlocks ] += Sum_bf[ monomial +block*2*numOfMonomials +stationaryPointIndex*numOfMonomials*2*numOfBlocks ];
			}


			linearSolver.set_givenMatrix(0);
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks,
					Sum_bf  +stationaryPointIndex*numOfMonomials*2*numOfBlocks, 
					coefficients +stationaryPointIndex*2*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT+stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks, 
					Sum_bf +numOfMonomials+stationaryPointIndex*numOfMonomials*2*numOfBlocks,
					coefficients +(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift );

			totalDisplacement_L2norm += 	 coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ] 
							*coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ];

		}

		stationaryPointShift += numOfStationaryPointsPerKernel;
		coefficientsShift = stationaryPointShift*2*numOfMonomials;

//		std::cout << "\r" << (stationaryPointShift)*nodesFraction << "% "; std::cout.flush();

	}


//	std::cout << "\r" << "100.00%\n\n";
	std::cout << "CUDA GPU process" << std::endl;
	std::cout << "Number of effective stationary points: \t\t" 	<< numOfEffectiveStationaryPoints << std::endl;
	std::cout << "Number of outlier stationary points: \t\t" 	<< numOfStationaryPoints -numOfEffectiveStationaryPoints << std::endl;
	std::cout << "L2 norm of the evaluated 2D displacement field: " << sqrt(totalDisplacement_L2norm) << std::endl << std::endl;


// ============================ Free Allocated resources  ==============================================

//	cudaStreamDestroy ( stream0 );// cudaStreamDestroy ( stream1 ); cudaStreamDestroy ( stream2 );

	cudaFreeHost ( powers_x );
	cudaFreeHost ( powers_y );
	cudaFreeHost ( Sum_bbT );
	cudaFreeHost ( Sum_bf );

	cudaHostUnregister (stationaryPoints_x);
	cudaHostUnregister (stationaryPoints_y);
	cudaHostUnregister (samplePoints_x);
	cudaHostUnregister (samplePoints_y);
	cudaHostUnregister (sampleValues_x);
	cudaHostUnregister (sampleValues_y);

	cudaFree ( gpu_powers_x );		cudaFree ( gpu_powers_y );
	cudaFree ( gpu_stationaryPoints_x );	cudaFree ( gpu_stationaryPoints_y );
	cudaFree ( gpu_samplePoints_x ); 	cudaFree ( gpu_samplePoints_y );
	cudaFree ( gpu_sampleValues_x );	cudaFree ( gpu_sampleValues_y );
	cudaFree ( gpu_bbT );
	cudaFree ( gpu_bf );

//	cudaDeviceReset();
}

template void cudaVersion_findCoefficientsLCS2DConstantMemory2(int numOfStationaryPoints, float* stationaryPoints_x, float* stationaryPoints_y, int numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, int numOfMonomials, int polDegree, float* coefficients, LinearSolver<float>& solver, WeightingFunction<float>& wFun, size_t deviceId);
template void cudaVersion_findCoefficientsLCS2DConstantMemory2(int numOfStationaryPoints, double* stationaryPoints_x, double* stationaryPoints_y, int numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, int numOfMonomials, int polDegree, double* coefficients, LinearSolver<double>& solver, WeightingFunction<double>& wFun, size_t deviceId);
*/



// ============================================================= algorithm design using the cuBLAS library ===================================================
template<class T>
void cudaVersion_findCoefficientsLCS2D3(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID )
{

	weightingFunction.set_numOfVariableData(0);
	cudaSetDevice( cudaDeviceID );

	//========================== size of arrays to be processed
	size_t sizeOfCoefficientsArray 		= 2*sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;
	size_t sizeOfPowersArray 		= sizeof(int)*numOfMonomials;
	size_t sizeOfbbTArray			= sizeof(T)*numOfMonomials*numOfMonomials;
	size_t sizeOfbfArray			= sizeof(T)*numOfMonomials;

	cudaError_t cudaStatus; 
	memset( coefficients, T(0), sizeOfCoefficientsArray );


	//========================= find monomials degree
	int *powers_x, *powers_y; int count(0);
	cudaStatus = cudaHostAlloc( (void**) &powers_x, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &powers_y, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	for (int i = 0; i<=polDegree; i++) for (int j = 0; j<=polDegree-i; j++) { powers_x[count] = j; powers_y[count] = i; count++; }


//	for (int i = 256+128; i<numOfSamplePoints; i++) { std::cout << sampleValues_y[i] << std::endl; }

	//========================== page lock host arrays
	cudaStatus = cudaHostRegister (stationaryPoints_x,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_y,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


// ==============================  Set the computational resources ===================================


	//========================= analyse hardware and available computing resources
//	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
//	dev_count = get_AllDevicesProperties(dev_prop);

//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);
//	plotGPUProperties(dev_prop[0]);

//	int totalAvailableThreads(0);
//	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
//		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;

//	std::cout << "totanl num of threads " << totalAvailableThreads << std::endl;
//	std::cout << "num of multiprocessors " << dev_prop[0].multiProcessorCount << std::endl;

	size_t maxSizeOfThreadBlock(1024);


	//========================== computational parameters per stationary point
	size_t numOfOperations	( numOfMonomials );
	size_t numOfDataUnits	( numOfSamplePoints );
	size_t numOfStreams	( 16 );


	if ( typeid(T) == typeid( double ) ) cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);


	// problem data analysis per stationary point
	size_t 	numOfStationaryPointsPerKernel 	= min( 512, numOfStationaryPoints );
	size_t	numOfSamplePointsPerKernel	= max(16, min( pow_i(2, 8 -numOfMonomials/3 -1), numOfSamplePoints ));
	size_t	numOfSamplePointsPerBlock	= min( numOfSamplePointsPerKernel, min(maxSizeOfThreadBlock, size_t(highestPowerof2( maxSizeOfThreadBlock/numOfMonomials )) ) );
	size_t	numOfIters			= numOfSamplePoints /numOfSamplePointsPerKernel;
	size_t 	numOfBlocks 			= numOfSamplePointsPerKernel /numOfSamplePointsPerBlock;
	dim3 	gridOfBlocks 			( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	blockSize 			( numOfSamplePointsPerBlock, numOfMonomials, 1); 


//	size_t	sizeOfBlockDataTile		= numOfSamplePointsPerBlock*sizeof(T);
	// powers_x +powers_y +max(b +sampleValues_x +sampleValues_y +bfx +bfy,   b +bbT)
//	size_t 	sharedMemorySize 		= 2*sizeOfPowersArray +sizeOfBlockDataTile +numOfMonomials*sizeOfBlockDataTile +max( 2*sizeOfBlockDataTile +2*4*sizeOfbfArray, sizeOfbbTArray);
//	size_t 	sharedMemorySize 		= 2*numOfMonomials*sizeOfBlockDataTile +sizeOfbbTArray;


//	std::cout << "\nCUDA execution settings:"
//		<< "\nNumber of stationary points in kernel \t\t\t" 		<< numOfStationaryPointsPerKernel
//		<< "\nNumber of sample points in kernel \t\t\t" 		<< numOfSamplePointsPerKernel
//		<< "\nNumber of sample points per block \t\t\t" 		<< numOfSamplePointsPerBlock
//		<< "\nNumber of kernel iterations per stationary point \t" 	<< numOfIters
//		<< "\nNumber of sample point blocks \t\t\t\t"			<< numOfBlocks
//		<< "\nNumber of residual sample point blocks \t\t\t"		<< numOfResidualBlocks
//		<< "\nNumber of residual sample points \t\t\t" 			<< numOfResidualSamplePoints
//		<< "\nShared memory per block \t\t\t\t"				<< sharedMemorySize << "b"
//		<< "\nShared memory in residual sample points block \t\t"	<< sharedMemorySize_divBlc << "b"
//		<< std::endl;



	//========================= gpu memory pointers
	T *gpu_coefficients, *gpu_stationaryPoints_x, *gpu_stationaryPoints_y, *gpu_samplePoints_x, *gpu_samplePoints_y, *gpu_sampleValues_x, *gpu_sampleValues_y, *gpu_bTheta, *gpu_b, *gpu_bbT, *gpu_bf;
	int *gpu_powers_x, *gpu_powers_y;



	//========================== allocate gpu memory
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_x, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_y, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_x, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_y, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bbT, 			sizeOfbbTArray*numOfStreams ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bf,			2*sizeOfbfArray*numOfStreams ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bTheta,			numOfSamplePoints*numOfMonomials*numOfStationaryPointsPerKernel*sizeof(T) ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_b,			numOfSamplePoints*numOfMonomials*numOfStationaryPointsPerKernel*sizeof(T) ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;



//================================== Launch the process ========================================================


	T *Sum_bbT, *Sum_bf;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bbT, sizeOfbbTArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &Sum_bf, 2*sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0); clock_t calculationClocks; T distance;

//	size_t totalGlobalMemory = 2*sizeOfPowersArray +4*sizeOfSamplePointsArray +2*sizeOfStationaryPointsArray +sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel +2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel;
//	std::cout << "Global memory required: \t\t\t\t" << totalGlobalMemory/1024 << "kb" << std::endl;

//	std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); calculationClocks = -clock();



	cudaStatus = cudaMemcpy( gpu_powers_x, 			powers_x, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_y, 			powers_y, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_x,		samplePoints_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_y, 		samplePoints_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_x, 		sampleValues_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_y, 		sampleValues_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_x, 	stationaryPoints_x, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_y, 	stationaryPoints_y, 	sizeOfStationaryPointsArray, 	H2D );


	cudaStream_t streams[ numOfStreams ];
	for ( size_t i = 0; i < numOfStreams; i++ ) cudaStreamCreate(streams+i);

	cublasHandle_t handles[ numOfStreams ];
	for ( size_t i = 0; i < numOfStreams; i++ ) cublasCreate( handles+i );
	for ( size_t i = 0; i < numOfStreams; i++ ) cublasSetStream( handles[i], streams[i] );
	cublasStatus_t cublasStatus;

	size_t streamIndex(0);

	const T alpha(1), beta(0);


	size_t stationaryPointShift(0), coefficientsShift(0);
	for (int stationaryPoint = 0; stationaryPoint < numOfStationaryPoints; stationaryPoint += numOfStationaryPointsPerKernel )
	{

		if ( stationaryPoint +numOfStationaryPointsPerKernel > numOfStationaryPoints ) {
			numOfStationaryPointsPerKernel 	= numOfStationaryPoints -stationaryPoint;
			gridOfBlocks 			= dim3( numOfBlocks, numOfStationaryPointsPerKernel, 1 ); }



		// ================================================ gpu processing part ========================================
		find_bTheta_LCS2D_cudaKernel <<< dim3(2, numOfStationaryPointsPerKernel ), dim3(32, numOfMonomials)>>>
		( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, 
		numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, weightingFunction, gpu_bTheta, gpu_b);
		cudaDeviceSynchronize();

		cudaMemset ( gpu_bf, 0, 2*sizeOfbfArray*numOfStreams );
		cudaMemset ( gpu_bbT, 0, sizeOfbbTArray*numOfStreams );
		cudaDeviceSynchronize();

		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			cublasStatus = cublasGgemm(handles[streamIndex], CUBLAS_OP_T, CUBLAS_OP_N, numOfMonomials, numOfMonomials, numOfSamplePoints, &alpha, gpu_b +numOfMonomials*numOfSamplePoints*stationaryPointIndex, numOfSamplePoints, gpu_bTheta +numOfMonomials*numOfSamplePoints*stationaryPointIndex, numOfSamplePoints, &beta, gpu_bbT+streamIndex*numOfMonomials*numOfMonomials, numOfMonomials);
			if (cublasStatus) std::cout << "ERROR!!" << std::endl;

			cublasStatus = cublasGgemv(handles[streamIndex], CUBLAS_OP_T, numOfSamplePoints, numOfMonomials, &alpha, gpu_bTheta+numOfMonomials*numOfSamplePoints*stationaryPointIndex, numOfSamplePoints, gpu_sampleValues_x, 1, &beta, gpu_bf +streamIndex*2*numOfMonomials, 1);
			if (cublasStatus) std::cout << "ERROR!!" << std::endl;

			cublasStatus = cublasGgemv(handles[streamIndex], CUBLAS_OP_T, numOfSamplePoints, numOfMonomials, &alpha, gpu_bTheta+numOfMonomials*numOfSamplePoints*stationaryPointIndex, numOfSamplePoints, gpu_sampleValues_y, 1, &beta, gpu_bf+numOfMonomials+streamIndex*2*numOfMonomials, 1);
			if (cublasStatus) std::cout << "ERROR!!" << std::endl;

			cudaMemcpyAsync( Sum_bbT+numOfMonomials*numOfMonomials*stationaryPointIndex, gpu_bbT+streamIndex*numOfMonomials*numOfMonomials, sizeOfbbTArray, D2H, streams[streamIndex]);
			cudaMemcpyAsync( Sum_bf+2*numOfMonomials*stationaryPointIndex, gpu_bf+streamIndex*numOfMonomials*2, 2*sizeOfbfArray, D2H, streams[streamIndex]);

			streamIndex++;
			if ( streamIndex == numOfStreams ) streamIndex = 0;
		}

		cudaDeviceSynchronize();

		//  =========================================== cpu post processing part  =========================================
		int counter(0);
		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			for ( size_t samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++ ) {
				distance = distance2D(stationaryPoints_x[stationaryPointIndex +stationaryPointShift ], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex +stationaryPointShift ], samplePoints_y[samplePointIndex]);
				if ( distance < weightingFunction.get_span() ) counter++; if ( counter >= numOfMonomials ) break;
			}
			if ( counter < numOfMonomials ) {
				std::fill(coefficients +stationaryPointIndex*2*numOfMonomials 	  +coefficientsShift,
					  coefficients +(stationaryPointIndex+1)*2*numOfMonomials +coefficientsShift,
					  T(0) );
				continue; }

			numOfEffectiveStationaryPoints++;

			linearSolver.set_givenMatrix(0);
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bf +2*stationaryPointIndex*numOfMonomials, 
					coefficients +stationaryPointIndex*2*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bf +(stationaryPointIndex*2+1)*numOfMonomials,
					coefficients +(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift );

			totalDisplacement_L2norm += 	 coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ] 
							*coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ];

//			for (int monomial = 0; monomial < 2*numOfMonomials; monomial++) std::cout << coefficients[monomial +stationaryPointIndex*2*numOfMonomials +coefficientsShift]  << std::endl;
		}



		stationaryPointShift += numOfStationaryPointsPerKernel;
		coefficientsShift = stationaryPointShift*2*numOfMonomials;

//		std::cout << "\r" << (stationaryPointShift)*nodesFraction << "% "; std::cout.flush();

	}

//	std::cout << "\r" << "100.00%\n\n";
	std::cout << "CUDA GPU process" << std::endl;
	std::cout << "Number of effective stationary points: \t\t" 	<< numOfEffectiveStationaryPoints << std::endl;
	std::cout << "Number of outlier stationary points: \t\t" 	<< numOfStationaryPoints -numOfEffectiveStationaryPoints << std::endl;
	std::cout << "L2 norm of the evaluated 2D displacement field: " << sqrt(totalDisplacement_L2norm) << std::endl << std::endl;


// ============================ Free Allocated resources  ==============================================

	for ( size_t i = 0; i < numOfStreams; i++ ) cublasDestroy( handles[i] );
	for ( size_t i = 0; i < numOfStreams; i++ ) cudaStreamDestroy( streams[i] );

	cudaFreeHost ( powers_x );
	cudaFreeHost ( powers_y );
	cudaFreeHost ( Sum_bbT );
	cudaFreeHost ( Sum_bf );


	cudaHostUnregister (stationaryPoints_x);
	cudaHostUnregister (stationaryPoints_y);
	cudaHostUnregister (samplePoints_x);
	cudaHostUnregister (samplePoints_y);
	cudaHostUnregister (sampleValues_x);
	cudaHostUnregister (sampleValues_y);

	cudaFree ( gpu_powers_x );		cudaFree ( gpu_powers_y );
	cudaFree ( gpu_stationaryPoints_x );	cudaFree ( gpu_stationaryPoints_y );
	cudaFree ( gpu_samplePoints_x ); 	cudaFree ( gpu_samplePoints_y );
	cudaFree ( gpu_sampleValues_x );	cudaFree ( gpu_sampleValues_y );
	cudaFree ( gpu_bTheta );
	cudaFree ( gpu_b );
	cudaFree ( gpu_bbT );
	cudaFree ( gpu_bf );



}

template void cudaVersion_findCoefficientsLCS2D3(int numOfStationaryPoints, float* stationaryPoints_x, float* stationaryPoints_y, int numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, int numOfMonomials, int polDegree, float* coefficients, LinearSolver<float>& solver, WeightingFunction<float>& wFun, size_t deviceId);
template void cudaVersion_findCoefficientsLCS2D3(int numOfStationaryPoints, double* stationaryPoints_x, double* stationaryPoints_y, int numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, int numOfMonomials, int polDegree, double* coefficients, LinearSolver<double>& solver, WeightingFunction<double>& wFun, size_t deviceId);






//  ===================================================  various algorithm designs for the order of the computations and the data transfers =================================
template<class T>
void cudaVersion_findCoefficientsLCS2D4(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID )
{

	weightingFunction.set_numOfVariableData(0);
	cudaSetDevice( cudaDeviceID );

	//========================== size of arrays to be processed
	size_t sizeOfCoefficientsArray 		= 2*sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;
	size_t sizeOfPowersArray 		= sizeof(int)*numOfMonomials;
	size_t sizeOfbbTArray			= sizeof(T)*numOfMonomials*numOfMonomials;
	size_t sizeOfbfArray			= sizeof(T)*numOfMonomials;

	cudaError_t cudaStatus; 
	memset( coefficients, T(0), sizeOfCoefficientsArray );


	//========================= find monomials degree
	int *powers_x, *powers_y; int count(0);
	cudaStatus = cudaHostAlloc( (void**) &powers_x, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &powers_y, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	for (int i = 0; i<=polDegree; i++) for (int j = 0; j<=polDegree-i; j++) { powers_x[count] = j; powers_y[count] = i; count++; }


//	for (int i = 256+128; i<numOfSamplePoints; i++) { std::cout << sampleValues_y[i] << std::endl; }

	//========================== page lock host arrays
	cudaStatus = cudaHostRegister (stationaryPoints_x,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_y,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


// ==============================  Set the computational resources ===================================


	//========================= analyse hardware and available computing resources
//	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
//	dev_count = get_AllDevicesProperties(dev_prop);

//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);
//	plotGPUProperties(dev_prop[0]);

//	int totalAvailableThreads(0);
//	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
//		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;

//	std::cout << "totanl num of threads " << totalAvailableThreads << std::endl;
//	std::cout << "num of multiprocessors " << dev_prop[0].multiProcessorCount << std::endl;

	size_t maxSizeOfThreadBlock(1024);


	//========================== computational parameters per stationary point
	size_t numOfOperations	( numOfMonomials );
	size_t numOfDataUnits	( numOfSamplePoints );
	size_t numOfStreams	( 3 );


	if ( typeid(T) == typeid( double ) ) cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);


	// problem data analysis per stationary point
	size_t 	numOfStationaryPointsPerKernel 	= min( 1*2048, numOfStationaryPoints );
	size_t	numOfSamplePointsPerKernel	= min( pow_i(2, 8 -numOfMonomials/3 -1), numOfSamplePoints );
	size_t	numOfSamplePointsPerBlock	= min( numOfSamplePointsPerKernel, min(maxSizeOfThreadBlock, size_t(highestPowerof2( maxSizeOfThreadBlock/numOfMonomials )) ) );
	size_t	numOfIters			= numOfSamplePoints /numOfSamplePointsPerKernel;
	size_t 	numOfBlocks 			= numOfSamplePointsPerKernel /numOfSamplePointsPerBlock;
	dim3 	gridOfBlocks 			( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	blockSize 			( numOfSamplePointsPerBlock, numOfMonomials, 1); 


//	size_t	sizeOfBlockDataTile		= numOfSamplePointsPerBlock*sizeof(T);
//	// powers_x +powers_y +max(b +sampleValues_x +sampleValues_y +bfx +bfy,   b +bbT)
//	size_t 	sharedMemorySize 		= 2*numOfMonomials*sizeOfBlockDataTile +sizeOfbbTArray;


//	std::cout << "\nCUDA execution settings:"
//		<< "\nNumber of stationary points in kernel \t\t\t" 		<< numOfStationaryPointsPerKernel
//		<< "\nNumber of sample points in kernel \t\t\t" 		<< numOfSamplePointsPerKernel
//		<< "\nNumber of sample points per block \t\t\t" 		<< numOfSamplePointsPerBlock
//		<< "\nNumber of kernel iterations per stationary point \t" 	<< numOfIters
//		<< "\nNumber of sample point blocks \t\t\t\t"			<< numOfBlocks
//		<< "\nShared memory per block \t\t\t\t"				<< sharedMemorySize << "b"
//		<< std::endl;


	//========================== create processing streams
//	cudaStream_t stream0, stream1, stream2;
//	cudaStreamCreate(&stream0);// cudaStreamCreate(&stream1); cudaStreamCreate(&stream2);


	//========================= gpu memory pointers
	T *gpu_coefficients, *gpu_stationaryPoints_x, *gpu_stationaryPoints_y, *gpu_samplePoints_x, *gpu_samplePoints_y, *gpu_sampleValues_x, *gpu_sampleValues_y, *gpu_bTheta, *gpu_b, *gpu_bbT, *gpu_bfx, *gpu_bfy;
	int *gpu_powers_x, *gpu_powers_y;


	//========================== allocate gpu memory
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_x, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_y, 	sizeOfStationaryPointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_x, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_y, 		sizeOfSamplePointsArray );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_x, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_y, 		sizeOfPowersArray );		if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bbT, 			sizeOfbbTArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bfx,			sizeOfbfArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bfy,			sizeOfbfArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bTheta,			numOfSamplePoints*numOfMonomials*numOfStationaryPointsPerKernel*sizeof(T) ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_b,			numOfSamplePoints*numOfMonomials*numOfStationaryPointsPerKernel*sizeof(T) ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;



//================================== Launch the process ========================================================


	T *Sum_bbT, *Sum_bfx, *Sum_bfy;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bbT, sizeOfbbTArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &Sum_bfx, sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bfy, sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0); clock_t calculationClocks; T distance;

//	size_t totalGlobalMemory = 2*sizeOfPowersArray +4*sizeOfSamplePointsArray +2*sizeOfStationaryPointsArray +sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel +2*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel;
//	std::cout << "Global memory required: \t\t\t\t" << totalGlobalMemory/1024 << "kb" << std::endl;


//	std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); calculationClocks = -clock();



	cudaStatus = cudaMemcpy( gpu_powers_x, 			powers_x, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_y, 			powers_y, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_x,		samplePoints_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_y, 		samplePoints_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_x, 		sampleValues_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_y, 		sampleValues_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_x, 	stationaryPoints_x, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_y, 	stationaryPoints_y, 	sizeOfStationaryPointsArray, 	H2D );


	typedef void (*kernelFun2D) ( int*, int*, T*, T*, size_t, T*, T*, T*, T*, WeightingFunction<T>, T*, T*, T*);

	kernelFun2D launchedFun( findCoefficientsLCS2D_cudaKernel );
	if (numOfMonomials == 1) launchedFun = findCoefficientsLCS2D_cudaKernel_0deg;
 


	size_t stationaryPointShift(0), coefficientsShift(0);
	for (int stationaryPoint = 0; stationaryPoint < numOfStationaryPoints; stationaryPoint += numOfStationaryPointsPerKernel )
	{

		if ( stationaryPoint +numOfStationaryPointsPerKernel > numOfStationaryPoints ) {
			numOfStationaryPointsPerKernel 	= numOfStationaryPoints -stationaryPoint;
			gridOfBlocks 			= dim3( numOfBlocks, numOfStationaryPointsPerKernel, 1 ); }


		// =================================== gpu processing part ========================================
		find_bTheta_LCS2D_cudaKernel <<< dim3(2, numOfStationaryPointsPerKernel ), dim3(32, numOfMonomials)>>>
		( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, 
		numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, weightingFunction, gpu_bTheta, gpu_b);
		cudaDeviceSynchronize();

//		find_bf_LCS2D_cudaKernel <<< dim3( 1, numOfMonomials, numOfStationaryPointsPerKernel ), dim3( numOfSamplePoints, 1, 1 ), sizeof(T)*2*numOfSamplePoints >>>
//		( numOfSamplePoints, gpu_sampleValues_x, gpu_sampleValues_y, numOfMonomials, gpu_bTheta, gpu_bfx, gpu_bfy );

		find_bf_LCS2D_cudaKernel2 <<< dim3( numOfMonomials, numOfStationaryPointsPerKernel, 1 ), dim3( 128, 1, 1 ), sizeof(T)*2*128 >>>
		( numOfSamplePoints, gpu_sampleValues_x, gpu_sampleValues_y, numOfMonomials, gpu_bTheta, gpu_bfx, gpu_bfy );


//		find_bTheta_bf_LCS2D_cudaKernel <<< dim3( 1, numOfMonomials, numOfStationaryPointsPerKernel ), dim3( numOfSamplePoints, 1, 1), sizeof(T)*2*numOfSamplePoints >>>
//		( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, numOfMonomials, numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, gpu_sampleValues_x, gpu_sampleValues_y, weightingFunction, gpu_bTheta, gpu_b, gpu_bfx, gpu_bfy);



//		find_bTheta_bf_LCS2D_cudaKernel2 <<< dim3( numOfMonomials, numOfStationaryPointsPerKernel, 1 ), dim3( 128, 1, 1), sizeof(T)*2*128 >>>
//		( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, numOfMonomials, numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, gpu_sampleValues_x, gpu_sampleValues_y, weightingFunction, gpu_bTheta, gpu_b, gpu_bfx, gpu_bfy);



		size_t numOfTiles = (numOfMonomials +16 -1) /16;
		find_bbT_LCS2D_cudaKernel <<< dim3( numOfTiles, numOfTiles, numOfStationaryPointsPerKernel ), dim3(16, 16, 1) >>>
		( numOfSamplePoints, numOfMonomials, gpu_bTheta, gpu_b, gpu_bbT );



//		launchedFun <<< dim3( numOfStationaryPointsPerKernel, 1, 1 ), dim3( 128, numOfMonomials, 1 ), sizeof(T)*2*128*numOfMonomials +sizeof(T)*numOfMonomials*numOfMonomials >>>
//		( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, gpu_sampleValues_x, gpu_sampleValues_y, weightingFunction, gpu_bbT, gpu_bfx, gpu_bfy);



//		findCoefficientsLCS2D_cudaKernel3 <<< dim3( numOfStationaryPointsPerKernel, 1, 1 ), dim3( 16, numOfMonomials, 1 ), sizeof(T)*2*16*numOfMonomials +4*sizeof(T)*numOfMonomials*numOfMonomials >>>
//		( gpu_powers_x, gpu_powers_y, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, gpu_sampleValues_x, gpu_sampleValues_y, weightingFunction, gpu_bbT, gpu_bfx, gpu_bfy);


		cudaDeviceSynchronize();
		cudaMemcpy( Sum_bbT, gpu_bbT, sizeOfbbTArray*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bfx, gpu_bfx, sizeOfbfArray*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bfy, gpu_bfy, sizeOfbfArray*numOfStationaryPointsPerKernel, D2H);
		cudaDeviceSynchronize();
//		for (int monomial = 0; monomial < numOfMonomials*numOfMonomials; monomial++) std::cout << Sum_bbT[monomial] << std::endl;
//		for (int monomial = 0; monomial < numOfMonomials; monomial++) std::cout << Sum_bfx[monomial] << std::endl;



		//  =========================================== cpu post processing part  =========================================
		int counter(0);
		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			for ( size_t samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++ ) {
				distance = distance2D(stationaryPoints_x[stationaryPointIndex +stationaryPointShift ], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex +stationaryPointShift ], samplePoints_y[samplePointIndex]);
				if ( distance < weightingFunction.get_span() ) counter++; if ( counter >= numOfMonomials ) break;
			}
			if ( counter < numOfMonomials ) {
				std::fill(coefficients +stationaryPointIndex*2*numOfMonomials 	  +coefficientsShift,
					  coefficients +(stationaryPointIndex+1)*2*numOfMonomials +coefficientsShift,
					  T(0) );
				continue; }

			numOfEffectiveStationaryPoints++;

			linearSolver.set_givenMatrix(0);
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bfx +stationaryPointIndex*numOfMonomials, 
					coefficients +stationaryPointIndex*2*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bfy +stationaryPointIndex*numOfMonomials,
					coefficients +(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift );

			totalDisplacement_L2norm += 	 coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ] 
							*coefficients[stationaryPointIndex*2*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*2+1)*numOfMonomials +coefficientsShift ];

//			for (int monomial = 0; monomial < 2*numOfMonomials; monomial++) std::cout << coefficients[monomial +stationaryPointIndex*2*numOfMonomials +coefficientsShift]  << std::endl;
		}

		stationaryPointShift += numOfStationaryPointsPerKernel;
		coefficientsShift = stationaryPointShift*2*numOfMonomials;

//		std::cout << "\r" << (stationaryPointShift)*nodesFraction << "% "; std::cout.flush();

	}


//	std::cout << "\r" << "100.00%\n\n";
	std::cout << "CUDA GPU process" << std::endl;
	std::cout << "Number of effective stationary points: \t\t" 	<< numOfEffectiveStationaryPoints << std::endl;
	std::cout << "Number of outlier stationary points: \t\t" 	<< numOfStationaryPoints -numOfEffectiveStationaryPoints << std::endl;
	std::cout << "L2 norm of the evaluated 2D displacement field: " << sqrt(totalDisplacement_L2norm) << std::endl << std::endl;


// ============================ Free Allocated resources  ==============================================

//	cudaStreamDestroy ( stream0 );// cudaStreamDestroy ( stream1 ); cudaStreamDestroy ( stream2 );

	cudaFreeHost ( powers_x );
	cudaFreeHost ( powers_y );
	cudaFreeHost ( Sum_bbT );
	cudaFreeHost ( Sum_bfx );
	cudaFreeHost ( Sum_bfy );


	cudaHostUnregister (stationaryPoints_x);
	cudaHostUnregister (stationaryPoints_y);
	cudaHostUnregister (samplePoints_x);
	cudaHostUnregister (samplePoints_y);
	cudaHostUnregister (sampleValues_x);
	cudaHostUnregister (sampleValues_y);

	cudaFree ( gpu_powers_x );		cudaFree ( gpu_powers_y );
	cudaFree ( gpu_stationaryPoints_x );	cudaFree ( gpu_stationaryPoints_y );
	cudaFree ( gpu_samplePoints_x ); 	cudaFree ( gpu_samplePoints_y );
	cudaFree ( gpu_sampleValues_x );	cudaFree ( gpu_sampleValues_y );
	cudaFree ( gpu_bTheta );
	cudaFree ( gpu_b );
	cudaFree ( gpu_bbT );
	cudaFree ( gpu_bfx );
	cudaFree ( gpu_bfy );

}

template void cudaVersion_findCoefficientsLCS2D4(int numOfStationaryPoints, float* stationaryPoints_x, float* stationaryPoints_y, int numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, int numOfMonomials, int polDegree, float* coefficients, LinearSolver<float>& solver, WeightingFunction<float>& wFun, size_t deviceId);
template void cudaVersion_findCoefficientsLCS2D4(int numOfStationaryPoints, double* stationaryPoints_x, double* stationaryPoints_y, int numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, int numOfMonomials, int polDegree, double* coefficients, LinearSolver<double>& solver, WeightingFunction<double>& wFun, size_t deviceId);



//	cudaMemcpyAsync( d_A0, h_A_pinned, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream0 );
//	cudaMemcpyAsync( d_B0, h_B_pinned, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream0 );
//	findCoefficientsLCS2D_cudaKernel <<< gridOfBlocksSize, blockSize, sharedMemorySize >>> 
//	(span, numOfStationaryPoints, gpu_stationaryPoints_x, gpu_stationaryPoints_y, numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, gpu_sampleValues_x, gpu_sampleValues_y, numOfMonomials, gpu_powers_x, gpu_powers_y, gpu_coefficients);


//	kernel <<< tileSize/256, 256, 0, stream0 >>> (d_A0, d_B0, d_C0 );
//	for ( size_t iteration = 0; iteration < numOfIterations; iteration += dataTileSize*numOfStreams )
//	{
//		cudaMemcpyAsync( d_A1, h_A_pinned +iteration +tileSize, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream1 );
//		cudaMemcpyAsync( d_B1, h_B_pinned +iteration +tileSize, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream1 );
//		kernel <<< tileSize/256, 256, 0, stream1 >>> ( d_A1, d_B1, d_C1 );
//		cudaMemcpyAsync( d_C0, h_C_pinned +iteration, sizeof(T)*tileSize, cudaMemcpyDeviceToHost, stream0 );
//		cudaMemcpy(coefficients, gpu_coefficients, sizeOfCoefficientsArray, cudaMemcpyDeviceToHost);


//		cudaMemcpyAsync( d_A2, h_A_pinned +iteration +2*tileSize, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream2 );
//		cudaMemcpyAsync( d_B2, h_B_pinned +iteration +2*tileSize, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream2 );
//		kernel <<< tileSize/256, 256, 0, stream2 >>> (d_A2, d_B2, d_C2 );
//		cudaMemcpyAsync( d_C1, h_C_pinned +iteration +tileSize, sizeof(T)*tileSize, cudaMemcpyDeviceToHost, stream1 );


//		cudaMemcpyAsync( d_A0, h_A_pinned +iteration +tileSize*3, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream0 );
//		cudaMemcpyAsync( d_B0, h_B_pinned +iteration +tileSize*3, sizeof(T)*tileSize, cudaMemcpyHostToDevice, stream0 );
//		kernel <<< tileSize/256, 256, 0, stream0 >>> (d_A0, d_B0, d_C0 );
//		cudaMemcpyAsync( d_C2, h_C_pinned +iteration +2*tileSize, sizeof(T)*tileSize, cudaMemcpyDeviceToHost, stream2 );
//	}



