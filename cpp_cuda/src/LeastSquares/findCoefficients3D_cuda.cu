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
#include"findCoefficients3D_cuda_kernels.h"
#include<math.h>
#include<cmath>
#include"cuBLAS_wrappers.h"


constexpr cudaMemcpyKind H2D = cudaMemcpyKind::cudaMemcpyHostToDevice;
constexpr cudaMemcpyKind D2H = cudaMemcpyKind::cudaMemcpyDeviceToHost;


// ===================================== optimal design, a given size of thread-block is processing iteratively in tiles the sample points, the block contains the data for all the monomials in order to be able to calculate the partial bbT which corresponds to the current sample points, minimized data transferring by calculating bbT and bf in the same kernel launch, several stationary points launched in a single kernel which can be controlled as a parameter =====================================================
template<class T>
void cudaVersion_findCoefficientsLCS3D(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID )
{

	cudaSetDevice( cudaDeviceID );

	//========================== size of arrays to be processed
	size_t sizeOfCoefficientsArray 		= 3*sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;
	size_t sizeOfPowersArray 		= sizeof(int)*numOfMonomials;
	size_t sizeOfbbTArray			= sizeof(T)*numOfMonomials*numOfMonomials;
	size_t sizeOfbfArray			= sizeof(T)*numOfMonomials;

	cudaError_t cudaStatus;
	memset( coefficients, T(0), sizeOfCoefficientsArray );


	//========================= find monomials degree
	int *powers_x, *powers_y, *powers_z; int count(0);
	cudaStatus = cudaHostAlloc( (void**) &powers_x, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &powers_y, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostAlloc( (void**) &powers_z, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	for (int i = 0; i<=polDegree; i++)
		for (int j = 0; j<=polDegree-i; j++)
			for (int k = 0; k<=polDegree-i-j; k++) 
			{ powers_x[count] = k; powers_y[count] = j; powers_z[count] = i; count++; }



//	for (int i = 256+128; i<numOfSamplePoints; i++) { std::cout << sampleValues_y[i] << std::endl; }

	//========================== page lock host arrays
	cudaStatus = cudaHostRegister (stationaryPoints_x,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_y,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_z,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_z,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_z,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


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
	size_t 	numOfStationaryPointsPerKernel 	= std::min( 1*2048, numOfStationaryPoints );

	size_t 	numOfBlocks 			= numOfSamplePoints /numOfSamplePointsPerBlock;
	dim3 	gridOfBlocks 			( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
	dim3	blockSize 			( numOfSamplePointsPerBlock, numOfMonomials, 1); 


	size_t	dataArraySizePerKernel		= numOfSamplePoints*numOfMonomials*numOfStationaryPointsPerKernel;
	size_t	sizeOfBlockDataTile		= numOfSamplePointsPerBlock*sizeof(T);
	// powers_x +powers_y +max(b +sampleValues_x +sampleValues_y +bfx +bfy,   b +bbT)
	size_t 	sharedMemorySize 		= 3*numOfMonomials*sizeOfBlockDataTile +sizeOfbbTArray;

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
	T *gpu_coefficients, *gpu_stationaryPoints_x, *gpu_stationaryPoints_y, *gpu_stationaryPoints_z, *gpu_samplePoints_x, *gpu_samplePoints_y, *gpu_samplePoints_z, *gpu_sampleValues_x, *gpu_sampleValues_y, *gpu_sampleValues_z, *gpu_bbT, *gpu_bfx, *gpu_bfy, *gpu_bfz;
	int *gpu_powers_x, *gpu_powers_y, *gpu_powers_z;



	//========================== allocate gpu memory
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_x, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_y, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_z, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_x, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_y, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_z, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_x, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_y, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_z, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_x, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_y, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_z, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bbT, 			sizeOfbbTArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bfx,			sizeOfbfArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bfy,			sizeOfbfArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bfz,			sizeOfbfArray*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;



//================================== Launch the process ========================================================


	T *Sum_bbT, *Sum_bfx, *Sum_bfy, *Sum_bfz;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bbT, sizeOfbbTArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &Sum_bfx, sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bfy, sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bfz, sizeOfbfArray*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0); clock_t calculationClocks; T distance;

	size_t totalGlobalMemory = 3*sizeOfPowersArray +6*sizeOfSamplePointsArray +3*sizeOfStationaryPointsArray +sizeOfbbTArray*numOfStationaryPointsPerKernel +3*sizeOfbfArray*numOfStationaryPointsPerKernel +3*dataArraySizePerKernel*sizeof(T);
	std::cout << "Global memory required: \t\t\t\t" << totalGlobalMemory/1024 << "kb" << std::endl;


//	std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); calculationClocks = -clock();



	cudaStatus = cudaMemcpy( gpu_powers_x, 			powers_x, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_y, 			powers_y, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_z, 			powers_z, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_x,		samplePoints_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_y, 		samplePoints_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_z, 		samplePoints_z, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_x, 		sampleValues_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_y, 		sampleValues_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_z, 		sampleValues_z, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_x, 	stationaryPoints_x, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_y, 	stationaryPoints_y, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_z, 	stationaryPoints_z, 	sizeOfStationaryPointsArray, 	H2D );


	typedef void (*kernelFun3D) ( int*, int*, int*, T*, T*, T*, size_t, T*, T*, T*, T*, T*, T*, WeightingFunction<T>, T*, T*, T*, T*);

	kernelFun3D launchedFun( findCoefficientsLCS3D_cudaKernel );
	if ( polDegree == 0 ) launchedFun = findCoefficientsLCS3D_cudaKernel_0deg;
 


	size_t stationaryPointShift(0), coefficientsShift(0);
	for (int stationaryPoint = 0; stationaryPoint < numOfStationaryPoints; stationaryPoint += numOfStationaryPointsPerKernel )
	{

		if ( stationaryPoint +numOfStationaryPointsPerKernel > numOfStationaryPoints ) {
			numOfStationaryPointsPerKernel 	= numOfStationaryPoints -stationaryPoint;
			gridOfBlocks 			= dim3( numOfBlocks, numOfStationaryPointsPerKernel, 1 ); }


		cuProfilerStart();

		// ===================================== gpu processing part ==================================
		launchedFun <<< dim3( numOfStationaryPointsPerKernel, 1, 1 ), dim3( numOfSamplePointsPerBlock, numOfMonomials, 1 ), sharedMemorySize >>>
		( gpu_powers_x, gpu_powers_y, gpu_powers_z, gpu_stationaryPoints_x +stationaryPointShift, gpu_stationaryPoints_y +stationaryPointShift, gpu_stationaryPoints_z +stationaryPointShift, numOfSamplePoints, gpu_samplePoints_x, gpu_samplePoints_y, gpu_samplePoints_z, gpu_sampleValues_x, gpu_sampleValues_y, gpu_sampleValues_z, weightingFunction, gpu_bbT, gpu_bfx, gpu_bfy, gpu_bfz );


		cudaDeviceSynchronize();
		cudaMemcpy( Sum_bbT, gpu_bbT, sizeOfbbTArray*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bfx, gpu_bfx, sizeOfbfArray*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bfy, gpu_bfy, sizeOfbfArray*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bfz, gpu_bfz, sizeOfbfArray*numOfStationaryPointsPerKernel, D2H);
		cudaDeviceSynchronize();

		cuProfilerStop();

		// ========================================= cpu post processing part ================================
		int counter(0);
		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			for ( size_t samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++ ) {
				distance = distance3D(stationaryPoints_x[stationaryPointIndex +stationaryPointShift ], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex +stationaryPointShift ], samplePoints_y[samplePointIndex], stationaryPoints_z[stationaryPointIndex +stationaryPointShift ], samplePoints_z[samplePointIndex] );
				if ( distance < weightingFunction.get_span() ) counter++; if ( counter >= numOfMonomials ) break;
			}
			if ( counter < numOfMonomials ) {
				std::fill(coefficients +stationaryPointIndex*3*numOfMonomials 	  +coefficientsShift,
					  coefficients +(stationaryPointIndex+1)*3*numOfMonomials +coefficientsShift,
					  T(0) );
				continue; }

			numOfEffectiveStationaryPoints++;

			linearSolver.set_givenMatrix(0);
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bfx +stationaryPointIndex*numOfMonomials, 
					coefficients +stationaryPointIndex*3*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bfy +stationaryPointIndex*numOfMonomials,
					coefficients +(stationaryPointIndex*3+1)*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials, Sum_bfz +stationaryPointIndex*numOfMonomials,
					coefficients +(stationaryPointIndex*3+2)*numOfMonomials +coefficientsShift );


			totalDisplacement_L2norm += 	 coefficients[stationaryPointIndex*3*numOfMonomials +coefficientsShift ] 
							*coefficients[stationaryPointIndex*3*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*3+1)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*3+1)*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*3+2)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*3+2)*numOfMonomials +coefficientsShift ];
		}

		stationaryPointShift += numOfStationaryPointsPerKernel;
		coefficientsShift = stationaryPointShift*3*numOfMonomials;

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
	cudaFreeHost ( powers_z );
	cudaFreeHost ( Sum_bbT );
	cudaFreeHost ( Sum_bfx );
	cudaFreeHost ( Sum_bfy );
	cudaFreeHost ( Sum_bfz );

	cudaHostUnregister (stationaryPoints_x);
	cudaHostUnregister (stationaryPoints_y);
	cudaHostUnregister (stationaryPoints_z);
	cudaHostUnregister (samplePoints_x);
	cudaHostUnregister (samplePoints_y);
	cudaHostUnregister (samplePoints_z);
	cudaHostUnregister (sampleValues_x);
	cudaHostUnregister (sampleValues_y);
	cudaHostUnregister (sampleValues_z);

	cudaFree ( gpu_powers_x );		cudaFree ( gpu_powers_y );		cudaFree ( gpu_powers_z );
	cudaFree ( gpu_stationaryPoints_x );	cudaFree ( gpu_stationaryPoints_y );	cudaFree ( gpu_stationaryPoints_z );
	cudaFree ( gpu_samplePoints_x ); 	cudaFree ( gpu_samplePoints_y );        cudaFree ( gpu_samplePoints_z );
	cudaFree ( gpu_sampleValues_x );	cudaFree ( gpu_sampleValues_y );        cudaFree ( gpu_sampleValues_z );
	cudaFree ( gpu_bbT );
	cudaFree ( gpu_bfx );
	cudaFree ( gpu_bfy );
	cudaFree ( gpu_bfz );

}

template void cudaVersion_findCoefficientsLCS3D(int numOfStationaryPoints, float* stationaryPoints_x, float* stationaryPoints_y, float* stationaryPoints_z, int numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* samplePoints_z, float* sampleValues_x, float* sampleValues_y, float* sampleValues_z, int numOfMonomials, int polDegree, float* coefficients, LinearSolver<float>& solver, WeightingFunction<float>& wFun, size_t deviceId);
template void cudaVersion_findCoefficientsLCS3D(int numOfStationaryPoints, double* stationaryPoints_x, double* stationaryPoints_y, double* stationaryPoints_z, int numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* samplePoints_z, double* sampleValues_x, double* sampleValues_y, double* sampleValues_z, int numOfMonomials, int polDegree, double* coefficients, LinearSolver<double>& solver, WeightingFunction<double>& wFun, size_t deviceId);




// ============================================== algorithm design for massive kernel launches, each processing a parametrized tile of data  ===============================================
template<class T>
void cudaVersion_findCoefficientsLCS3D2(int numOfStationaryPoints, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, int numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, int numOfMonomials, int polDegree, T* coefficients, LinearSolver<T>& linearSolver, WeightingFunction<T>& weightingFunction, size_t cudaDeviceID )
{
	cudaSetDevice( cudaDeviceID );

	//========================== size of arrays to be processed
	size_t sizeOfCoefficientsArray 		= 3*sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;
	size_t sizeOfPowersArray 		= sizeof(int)*numOfMonomials;
	size_t sizeOfbbTArray			= sizeof(T)*numOfMonomials*numOfMonomials;
	size_t sizeOfbfArray			= sizeof(T)*numOfMonomials;

	cudaError_t cudaStatus;
	memset( coefficients, T(0), sizeOfCoefficientsArray );


	//========================= find monomials degree
	int *powers_x, *powers_y, *powers_z; int count(0);
	cudaStatus = cudaHostAlloc( (void**) &powers_x, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &powers_y, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostAlloc( (void**) &powers_z, sizeOfPowersArray, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	for (int i = 0; i<=polDegree; i++)
		for (int j = 0; j<=polDegree-i; j++)
			for (int k = 0; k<=polDegree-i-j; k++) 
			{ powers_x[count] = k; powers_y[count] = j; powers_z[count] = i; count++; }



//	for (int i = 256+128; i<numOfSamplePoints; i++) { std::cout << sampleValues_y[i] << std::endl; }

	//========================== page lock host arrays
	cudaStatus = cudaHostRegister (stationaryPoints_x,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_y,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (stationaryPoints_z,sizeOfStationaryPointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (samplePoints_z,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_x,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_y,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaHostRegister (sampleValues_z,	sizeOfSamplePointsArray,	cudaHostRegisterDefault );	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


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
	size_t 	numOfStationaryPointsPerKernel 	= min( 2048, numOfStationaryPoints );
	size_t	numOfSamplePointsPerKernel	= min( pow_i(2, 8 -numOfMonomials/3 -1), numOfSamplePoints );
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
	size_t 	sharedMemorySize 		= 3*numOfMonomials*sizeOfBlockDataTile +sizeOfbbTArray;
	size_t	sizeOfResidualDataTile		= numOfResidualSamplePoints*sizeof(T);
//	size_t 	sharedMemorySize_divBlc		= 2*sizeOfPowersArray +sizeOfBlockDataTile +numOfMonomials*sizeOfResidualDataTile +max( 2*sizeOfResidualDataTile +2*4*sizeOfbfArray, sizeOfbbTArray);
	size_t 	sharedMemorySize_divBlc		= 3*numOfMonomials*sizeOfResidualDataTile +sizeOfbbTArray;


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
	T *gpu_coefficients, *gpu_stationaryPoints_x, *gpu_stationaryPoints_y, *gpu_stationaryPoints_z, *gpu_samplePoints_x, *gpu_samplePoints_y, *gpu_samplePoints_z, *gpu_sampleValues_x, *gpu_sampleValues_y, *gpu_sampleValues_z, *gpu_bbT, *gpu_bf;
	int *gpu_powers_x, *gpu_powers_y, *gpu_powers_z;


	//========================== allocate gpu memory
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_x, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_y, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_stationaryPoints_z, 	sizeOfStationaryPointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_x, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_y, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_samplePoints_z, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_x, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_y, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_sampleValues_z, 		sizeOfSamplePointsArray );				if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_x, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_y, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_powers_z, 		sizeOfPowersArray );					if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bbT, 			sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
	cudaStatus = cudaMalloc ((void**) &gpu_bf,			3*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel ); 	if (cudaStatus) std::cout << "ERROR!!" << std::endl;



//================================== Launch the process ========================================================


	T *Sum_bbT, *Sum_bf;
	cudaStatus = cudaHostAlloc( (void**) &Sum_bbT, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl; 
	cudaStatus = cudaHostAlloc( (void**) &Sum_bf, 3*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel, cudaHostAllocDefault ); if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	size_t numOfEffectiveStationaryPoints(0); T totalDisplacement_L2norm(0); clock_t calculationClocks; T distance;

	size_t totalGlobalMemory = 3*sizeOfPowersArray +6*sizeOfSamplePointsArray +3*sizeOfStationaryPointsArray +sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel +3*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel;

	std::cout << "Global memory required: \t\t\t\t" << totalGlobalMemory/1024 << "kb" << std::endl;

//	std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfStationaryPoints); calculationClocks = -clock();



	cudaStatus = cudaMemcpy( gpu_powers_x, 			powers_x, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_y, 			powers_y, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_powers_z, 			powers_z, 		sizeOfPowersArray, 	 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_x,		samplePoints_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_y, 		samplePoints_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_samplePoints_z, 		samplePoints_z, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_x, 		sampleValues_x, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_y, 		sampleValues_y, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_sampleValues_z, 		sampleValues_z, 	sizeOfSamplePointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_x, 	stationaryPoints_x, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_y, 	stationaryPoints_y, 	sizeOfStationaryPointsArray, 	H2D );
	cudaStatus = cudaMemcpy( gpu_stationaryPoints_z, 	stationaryPoints_z, 	sizeOfStationaryPointsArray, 	H2D );


	typedef void (*kernelFun3D) ( int*, int*, int*, T*, T*, T*, T*, T*, T*, T*, T*, T*, size_t, WeightingFunction<T>, T*, T*);

	kernelFun3D launchedFun(findCoefficientsLCS3D_cudaKernel2);
	if ( polDegree == 0 ) launchedFun = findCoefficientsLCS3D_cudaKernel_0deg2;


	size_t stationaryPointShift(0), coefficientsShift(0);
	for (int stationaryPoint = 0; stationaryPoint < numOfStationaryPoints; stationaryPoint += numOfStationaryPointsPerKernel )
	{

		if ( stationaryPoint +numOfStationaryPointsPerKernel > numOfStationaryPoints ) {
 			numOfStationaryPointsPerKernel 	= numOfStationaryPoints -stationaryPoint;
			gridOfBlocks 			= dim3( numOfBlocks, numOfStationaryPointsPerKernel, 1 );
			divGridOfBlocks 		= dim3( numOfResidualBlocks, numOfStationaryPointsPerKernel, 1 ); }


		// ====================================== gpu processing part =========================================
		cudaMemset ( gpu_bf, 0, 3*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel);
		cudaMemset ( gpu_bbT, 0, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel);
		cudaDeviceSynchronize();

		cuProfilerStart();

		for ( int iteration = 0; iteration < numOfIters; iteration++ )
		{
			size_t dataShift = iteration*numOfBlocks*numOfSamplePointsPerBlock;
			launchedFun <<< gridOfBlocks, blockSize, sharedMemorySize>>>
			( gpu_powers_x, gpu_powers_y, gpu_powers_z, gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, gpu_stationaryPoints_z +stationaryPoint, 
			gpu_samplePoints_x +dataShift, gpu_samplePoints_y +dataShift, gpu_samplePoints_z +dataShift, gpu_sampleValues_x +dataShift, gpu_sampleValues_y +dataShift, gpu_sampleValues_z +dataShift, numOfBlocks, weightingFunction, gpu_bbT, gpu_bf);
			cudaDeviceSynchronize();
		}

		cuProfilerStop();

		if ( divGridFlag ) {
			size_t dataShift = numOfIters*numOfBlocks*numOfSamplePointsPerBlock;
			launchedFun <<< divGridOfBlocks, blockSize, sharedMemorySize>>>
			( gpu_powers_x, gpu_powers_y, gpu_powers_z, gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, gpu_stationaryPoints_z +stationaryPoint,
			gpu_samplePoints_x +dataShift, gpu_samplePoints_y +dataShift, gpu_samplePoints_z +dataShift, gpu_sampleValues_x +dataShift, gpu_sampleValues_y +dataShift, gpu_sampleValues_z +dataShift, numOfBlocks, weightingFunction, gpu_bbT, gpu_bf); }


		if ( divBlockFlag ) {
			size_t dataShift = (numOfIters*numOfBlocks +numOfResidualBlocks)*numOfSamplePointsPerBlock;
			launchedFun <<< dim3(1, numOfStationaryPointsPerKernel, 1), divBlockSize, sharedMemorySize_divBlc>>>
			( gpu_powers_x, gpu_powers_y, gpu_powers_z, gpu_stationaryPoints_x +stationaryPoint, gpu_stationaryPoints_y +stationaryPoint, gpu_stationaryPoints_z +stationaryPoint, 
			gpu_samplePoints_x +dataShift, gpu_samplePoints_y +dataShift, gpu_samplePoints_z +dataShift, gpu_sampleValues_x +dataShift, gpu_sampleValues_y +dataShift, gpu_sampleValues_z +dataShift, numOfBlocks, weightingFunction, gpu_bbT, gpu_bf);}

		cudaDeviceSynchronize();

		cudaMemcpy( Sum_bbT, gpu_bbT, sizeOfbbTArray*numOfBlocks*numOfStationaryPointsPerKernel, D2H);
		cudaMemcpy( Sum_bf, gpu_bf, 3*sizeOfbfArray*numOfBlocks*numOfStationaryPointsPerKernel, D2H);
		cudaDeviceSynchronize();



		// ====================================== cpu post processing part ====================================
		int counter(0);
		for (int stationaryPointIndex = 0; stationaryPointIndex < numOfStationaryPointsPerKernel; stationaryPointIndex++ ) {

			counter = 0;
			for ( size_t samplePointIndex = 0; samplePointIndex < numOfSamplePoints; samplePointIndex++ ) {
				distance = distance3D(stationaryPoints_x[stationaryPointIndex +stationaryPointShift ], samplePoints_x[samplePointIndex], stationaryPoints_y[stationaryPointIndex +stationaryPointShift ], samplePoints_y[samplePointIndex], stationaryPoints_z[stationaryPointIndex +stationaryPointShift ], samplePoints_z[samplePointIndex] );
				if ( distance < weightingFunction.get_span() ) counter++; if ( counter >= numOfMonomials ) break;
			}
			if ( counter < numOfMonomials ) {
				std::fill(coefficients +stationaryPointIndex*3*numOfMonomials 	  +coefficientsShift,
					  coefficients +(stationaryPointIndex+1)*3*numOfMonomials +coefficientsShift,
					  T(0) );
				continue; }

			numOfEffectiveStationaryPoints++;


			for ( int block = 1; block < numOfBlocks; block++ ) {
				for (int monomial = 0; monomial < numOfMonomials*numOfMonomials; monomial++)
					Sum_bbT[ monomial +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks ] += Sum_bbT[ monomial +block*numOfMonomials*numOfMonomials +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks ];
				for (int monomial = 0; monomial < 3*numOfMonomials; monomial++)
					Sum_bf[ monomial +stationaryPointIndex*numOfMonomials*3*numOfBlocks ] += Sum_bf[ monomial +block*3*numOfMonomials +stationaryPointIndex*numOfMonomials*3*numOfBlocks ];
			}


			linearSolver.set_givenMatrix(0);
			linearSolver( 	Sum_bbT +stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks,
					Sum_bf  +stationaryPointIndex*numOfMonomials*3*numOfBlocks, 
					coefficients +stationaryPointIndex*3*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT+stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks, 
					Sum_bf +numOfMonomials+stationaryPointIndex*numOfMonomials*3*numOfBlocks,
					coefficients +(stationaryPointIndex*3+1)*numOfMonomials +coefficientsShift );
			linearSolver( 	Sum_bbT+stationaryPointIndex*numOfMonomials*numOfMonomials*numOfBlocks, 
					Sum_bf +2*numOfMonomials+stationaryPointIndex*numOfMonomials*3*numOfBlocks,
					coefficients +(stationaryPointIndex*3+2)*numOfMonomials +coefficientsShift );


			totalDisplacement_L2norm += 	 coefficients[stationaryPointIndex*3*numOfMonomials +coefficientsShift ] 
							*coefficients[stationaryPointIndex*3*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*3+1)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*3+1)*numOfMonomials +coefficientsShift ]
							+coefficients[(stationaryPointIndex*3+2)*numOfMonomials +coefficientsShift ] 
							*coefficients[(stationaryPointIndex*3+2)*numOfMonomials +coefficientsShift ];
		}

		stationaryPointShift += numOfStationaryPointsPerKernel;
		coefficientsShift = stationaryPointShift*3*numOfMonomials;

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
	cudaFreeHost ( powers_z );
	cudaFreeHost ( Sum_bbT );
	cudaFreeHost ( Sum_bf );

	cudaHostUnregister (stationaryPoints_x);
	cudaHostUnregister (stationaryPoints_y);
	cudaHostUnregister (stationaryPoints_z);
	cudaHostUnregister (samplePoints_x);
	cudaHostUnregister (samplePoints_y);
	cudaHostUnregister (samplePoints_z);
	cudaHostUnregister (sampleValues_x);
	cudaHostUnregister (sampleValues_y);
	cudaHostUnregister (sampleValues_z);

	cudaFree ( gpu_powers_x );		cudaFree ( gpu_powers_y );		cudaFree ( gpu_powers_z );
	cudaFree ( gpu_stationaryPoints_x );	cudaFree ( gpu_stationaryPoints_y );    cudaFree ( gpu_stationaryPoints_z );
	cudaFree ( gpu_samplePoints_x ); 	cudaFree ( gpu_samplePoints_y );        cudaFree ( gpu_samplePoints_z );
	cudaFree ( gpu_sampleValues_x );	cudaFree ( gpu_sampleValues_y );        cudaFree ( gpu_sampleValues_z );
	cudaFree ( gpu_bbT );
	cudaFree ( gpu_bf );

}

template void cudaVersion_findCoefficientsLCS3D2(int numOfStationaryPoints, float* stationaryPoints_x, float* stationaryPoints_y, float* stationaryPoints_z, int numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* samplePoints_z, float* sampleValues_x, float* sampleValues_y, float* sampleValues_z, int numOfMonomials, int polDegree, float* coefficients, LinearSolver<float>& solver, WeightingFunction<float>& wFun, size_t deviceId);
template void cudaVersion_findCoefficientsLCS3D2(int numOfStationaryPoints, double* stationaryPoints_x, double* stationaryPoints_y, double* stationaryPoints_z, int numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* samplePoints_z, double* sampleValues_x, double* sampleValues_y, double* sampleValues_z, int numOfMonomials, int polDegree, double* coefficients, LinearSolver<double>& solver, WeightingFunction<double>& wFun, size_t deviceId);

