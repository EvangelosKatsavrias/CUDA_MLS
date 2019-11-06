#include<iostream>
#include<stdio.h>
#include"weightedLeastSquares_cuda.h"
#include"gpuDeviceProperties.h"
#include"linearizers_cuda.h"
#include"atomicOperations_doubles.h"
#include"weightingFunction_cuda.h"
#include"findStationaryPointsRanges_cuda.h"


template<class T>
__global__
void find_reciprocalWeightsSum_cudaKernel(int *ranges, T reciprocal_span, T* stationaryPoints, T* evaluationPoints, T* recip_weightsSums)
{
	__shared__ T weight_sum[1]; if ( threadIdx.x == 0 ) *weight_sum = 0;


	int lowerStationaryPointIndex = ranges[blockIdx.x*2], numOfEffectiveStationaryPoints = ranges[blockIdx.x*2+1] -lowerStationaryPointIndex +1;


	T evaluationPoint(evaluationPoints[blockIdx.x]),
	stationaryPoint(stationaryPoints[lowerStationaryPointIndex +threadIdx.x]),
	weightPerStationaryPoint( weightingFunction3(reciprocal_span, fabs(evaluationPoint - stationaryPoint) ) );


//	if ( threadIdx.x < numOfEffectiveStationaryPoints ) atomicAdd( weight_sum, weightPerStationaryPoint ); __syncthreads();
	if ( threadIdx.x == 0 )	recip_weightsSums[blockIdx.x] = 1/(*weight_sum);
}


//==========================================
template<class T>
__global__
void evaluationMethod_cudaKernel(int *ranges, T reciprocal_span, T* stationaryPoints, T* evaluationPoints, int numOfMonomials, T* coefficients, T* recip_weightsSums, T* evaluations)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *local_coefficients	= reinterpret_cast<T*>(sharedMemory),
	  *evaluation		= local_coefficients +numOfMonomials,
	  *recip_weightsSum	= evaluation +1;
	if ( threadIdx.y == 0 ) { local_coefficients[threadIdx.x] = 0; *evaluation = 0; }


	int lowerStationaryPointIndex = ranges[blockIdx.x*2], numOfEffectiveStationaryPoints = ranges[blockIdx.x*2+1] -lowerStationaryPointIndex +1;

	if ( threadIdx.x == 0 )
		for (int memoryBank = 0; memoryBank < blockDim.x*blockDim.y/32+1; memoryBank++) 
			if ( threadIdx.y == memoryBank*32 && threadIdx.y < numOfEffectiveStationaryPoints ) *recip_weightsSum = recip_weightsSums[blockIdx.x];
	__syncthreads();


	T evaluationPoint(evaluationPoints[blockIdx.x]),
	stationaryPoint(stationaryPoints[lowerStationaryPointIndex +threadIdx.y]),
	weightPerStationaryPoint(weightingFunction3(reciprocal_span, fabs(evaluationPoint - stationaryPoint) )),
	local_coefficient(coefficients[(lowerStationaryPointIndex+threadIdx.y)*numOfMonomials+threadIdx.x] *weightPerStationaryPoint *(*recip_weightsSum)),
	b(pow( evaluationPoint, threadIdx.x));

/*
	if ( threadIdx.y < numOfEffectiveStationaryPoints ) atomicAdd(local_coefficients+threadIdx.x, local_coefficient); __syncthreads();
	if ( threadIdx.y == 0 )	atomicAdd(evaluation, b *local_coefficients[threadIdx.x]);__syncthreads();
	if ( threadIdx.x == 0 && threadIdx.y == 0 ) evaluations[blockIdx.x] = *evaluation;
*/}


template<class T>
void cudaVersion_evaluationMethod(int numOfEvaluationPoints, T* evaluationPoints, T span, int numOfStationaryPoints, T* stationaryPoints, int numOfMonomials, T* host_coefficients, T* host_evaluations)
{
	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
	dev_count = get_AllDevicesProperties(dev_prop);
//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);

	int totalAvailableThreads(0);
	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;


	T *gpu_coefficients, *gpu_stationaryPoints, *gpu_evaluationPoints, *gpu_evaluations, *gpu_recip_weightsSums;
	int *gpu_ranges, maxNumOfEffectiveStationaryPoints, *gpu_maxNumOfEffectiveStationaryPoints;
	size_t sizeOfCoefficientsArray		= sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfEvaluationPointsArray	= sizeof(T)*numOfEvaluationPoints;
	size_t sizeOfRanges			= sizeof(int)*2*numOfEvaluationPoints;
	

	
	cudaError_t 	cudaMallocStatus = cudaMalloc( (void**) &gpu_coefficients,	sizeOfCoefficientsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_stationaryPoints,	sizeOfStationaryPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_evaluationPoints,	sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_evaluations,	sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_ranges,		sizeOfRanges);			if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_maxNumOfEffectiveStationaryPoints, sizeof(int));	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_recip_weightsSums, sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;



	cudaMemcpy(gpu_stationaryPoints, stationaryPoints, sizeOfStationaryPointsArray, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_evaluationPoints, evaluationPoints, sizeOfEvaluationPointsArray, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_coefficients, 	 host_coefficients,sizeOfCoefficientsArray, 	cudaMemcpyDeviceToHost);



	int maxBlockSize = min(numOfEvaluationPoints, 256);
	find_samplePointsRanges <<< dim3( ceil(numOfEvaluationPoints/maxBlockSize), 1, 1), dim3(maxBlockSize, 1, 1) >>> 
	(numOfEvaluationPoints, gpu_evaluationPoints, numOfStationaryPoints, gpu_stationaryPoints, span, gpu_ranges, gpu_maxNumOfEffectiveStationaryPoints);

	cudaMemcpy(&maxNumOfEffectiveStationaryPoints, gpu_maxNumOfEffectiveStationaryPoints, sizeof(int), cudaMemcpyDeviceToHost);
//	printf("Max number of effective stationary points %d\n", maxNumOfEffectiveStationaryPoints);


	find_reciprocalWeightsSum_cudaKernel <<< dim3(numOfEvaluationPoints, 1, 1), dim3(maxNumOfEffectiveStationaryPoints, 1, 1)  >>> 
	(gpu_ranges, T(1/span), gpu_stationaryPoints, gpu_evaluationPoints, gpu_recip_weightsSums);


	evaluationMethod_cudaKernel <<< dim3(numOfEvaluationPoints, 1, 1), dim3(numOfMonomials, maxNumOfEffectiveStationaryPoints, 1), size_t((numOfMonomials+2)*sizeof(T)) >>> 
	(gpu_ranges, T(1/span), gpu_stationaryPoints, gpu_evaluationPoints, numOfMonomials, gpu_coefficients, gpu_recip_weightsSums, gpu_evaluations);


	cudaMemcpy(host_evaluations, gpu_evaluations, sizeOfEvaluationPointsArray, cudaMemcpyDeviceToHost);
	cudaFree(gpu_coefficients); cudaFree(gpu_stationaryPoints);
	cudaFree(gpu_evaluationPoints); cudaFree(gpu_evaluations); cudaFree(gpu_ranges);
	cudaFree(gpu_maxNumOfEffectiveStationaryPoints); cudaFree(gpu_recip_weightsSums);
}


template void cudaVersion_evaluationMethod(int numOfEvaluationPoints, float* evaluationPoints, float span, int numOfStationaryPoints, float* stationaryPoints, int numOfMonomials, float* host_coefficients, float* evaluations);
template void cudaVersion_evaluationMethod(int numOfEvaluationPoints, double* evaluationPoints, double span, int numOfStationaryPoints, double* stationaryPoints, int numOfMonomials, double* host_coefficients, double* evaluations);



// ====================================================
template<class T>
__global__
void evaluationMethodLCS_cudaKernel(int *ranges, T reciprocal_span, T* stationaryPoints, T* evaluationPoints, int numOfMonomials, T* coefficients, T* recip_weightsSums, T* evaluations)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *local_coefficients	= reinterpret_cast<T*>(sharedMemory),
	  *evaluation		= local_coefficients +numOfMonomials,
	  *recip_weightsSum	= evaluation +1;
	if ( threadIdx.y == 0 ) { local_coefficients[threadIdx.x] = 0; *evaluation = 0; }


	int lowerStationaryPointIndex = ranges[blockIdx.x*2], numOfEffectiveStationaryPoints = ranges[blockIdx.x*2+1] -lowerStationaryPointIndex +1;

	if ( threadIdx.x == 0 )
		for (int memoryBank = 0; memoryBank < blockDim.x*blockDim.y/32+1; memoryBank++) 
			if ( threadIdx.y == memoryBank*32 && threadIdx.y < numOfEffectiveStationaryPoints ) *recip_weightsSum = recip_weightsSums[blockIdx.x];
	__syncthreads();


	T evaluationPoint 	(evaluationPoints[blockIdx.x]),
	stationaryPoint 	(stationaryPoints[lowerStationaryPointIndex +threadIdx.y]),
	weightPerStationaryPoint (weightingFunction3(reciprocal_span, fabs(evaluationPoint - stationaryPoint) )),
	local_coefficient 	(coefficients[(lowerStationaryPointIndex+threadIdx.y)*numOfMonomials+threadIdx.x] *weightPerStationaryPoint *(*recip_weightsSum)),
	b 			(pow( evaluationPoint -stationaryPoint, threadIdx.x));

//	if ( threadIdx.y < numOfEffectiveStationaryPoints ) atomicAdd(evaluation, b *local_coefficient); __syncthreads();
	if ( threadIdx.x == 0 && threadIdx.y == 0 ) evaluations[blockIdx.x] = *evaluation;
}


template<class T>
void cudaVersion_evaluationMethodLCS(int numOfEvaluationPoints, T* evaluationPoints, T span, int numOfStationaryPoints, T* stationaryPoints, int numOfMonomials, T* host_coefficients, T* host_evaluations)
{
	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
	dev_count = get_AllDevicesProperties(dev_prop);
//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);

	int totalAvailableThreads(0);
	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;


	T *gpu_coefficients, *gpu_stationaryPoints, *gpu_evaluationPoints, *gpu_evaluations, *gpu_recip_weightsSums;
	int *gpu_ranges, maxNumOfEffectiveStationaryPoints, *gpu_maxNumOfEffectiveStationaryPoints;
	size_t sizeOfCoefficientsArray		= sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfEvaluationPointsArray	= sizeof(T)*numOfEvaluationPoints;
	size_t sizeOfRanges			= sizeof(int)*2*numOfEvaluationPoints;
	

	
	cudaError_t 	cudaMallocStatus = cudaMalloc( (void**) &gpu_coefficients,	sizeOfCoefficientsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_stationaryPoints,	sizeOfStationaryPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_evaluationPoints,	sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_evaluations,	sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_ranges,		sizeOfRanges);			if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_maxNumOfEffectiveStationaryPoints, sizeof(int));	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_recip_weightsSums, sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;



	cudaMemcpy(gpu_stationaryPoints, stationaryPoints, sizeOfStationaryPointsArray, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_evaluationPoints, evaluationPoints, sizeOfEvaluationPointsArray, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_coefficients, 	 host_coefficients,sizeOfCoefficientsArray, 	cudaMemcpyDeviceToHost);



	int maxBlockSize = min(numOfEvaluationPoints, 256);
	find_samplePointsRanges <<< dim3( ceil(numOfEvaluationPoints/maxBlockSize), 1, 1), dim3(maxBlockSize, 1, 1) >>> 
	(numOfEvaluationPoints, gpu_evaluationPoints, numOfStationaryPoints, gpu_stationaryPoints, span, gpu_ranges, gpu_maxNumOfEffectiveStationaryPoints);

	cudaMemcpy(&maxNumOfEffectiveStationaryPoints, gpu_maxNumOfEffectiveStationaryPoints, sizeof(int), cudaMemcpyDeviceToHost);
//	printf("Max number of effective stationary points %d\n", maxNumOfEffectiveStationaryPoints);


	find_reciprocalWeightsSum_cudaKernel <<< dim3(numOfEvaluationPoints, 1, 1), dim3(maxNumOfEffectiveStationaryPoints, 1, 1)  >>> 
	(gpu_ranges, T(1/span), gpu_stationaryPoints, gpu_evaluationPoints, gpu_recip_weightsSums);


	evaluationMethodLCS_cudaKernel <<< dim3(numOfEvaluationPoints, 1, 1), dim3(numOfMonomials, maxNumOfEffectiveStationaryPoints, 1), size_t((numOfMonomials+2)*sizeof(T)) >>> 
	(gpu_ranges, T(1/span), gpu_stationaryPoints, gpu_evaluationPoints, numOfMonomials, gpu_coefficients, gpu_recip_weightsSums, gpu_evaluations);


	cudaMemcpy(host_evaluations, gpu_evaluations, sizeOfEvaluationPointsArray, cudaMemcpyDeviceToHost);
	cudaFree(gpu_coefficients); cudaFree(gpu_stationaryPoints);
	cudaFree(gpu_evaluationPoints); cudaFree(gpu_evaluations); cudaFree(gpu_ranges);
	cudaFree(gpu_maxNumOfEffectiveStationaryPoints); cudaFree(gpu_recip_weightsSums);
}


template void cudaVersion_evaluationMethodLCS(int numOfEvaluationPoints, float* evaluationPoints, float span, int numOfStationaryPoints, float* stationaryPoints, int numOfMonomials, float* host_coefficients, float* evaluations);
template void cudaVersion_evaluationMethodLCS(int numOfEvaluationPoints, double* evaluationPoints, double span, int numOfStationaryPoints, double* stationaryPoints, int numOfMonomials, double* host_coefficients, double* evaluations);



// ============================================
template<class T>
__global__
void evaluationMethodMovingLS_cudaKernel(T* evaluationPoints, int numOfMonomials, T* coefficients, T* evaluations)
{
	__shared__ T evaluationPoint[1], evaluation[1];

	if ( threadIdx.x == 0 ) { *evaluation = 0; *evaluationPoint = evaluationPoints[blockIdx.x]; }
	__syncthreads();

	T local_coefficient ( coefficients[blockIdx.x*numOfMonomials+threadIdx.x] ),
	  b ( pow( *evaluationPoint, threadIdx.x) );

//	atomicAdd(evaluation, b *local_coefficient); __syncthreads();
	if ( threadIdx.x == 0 ) evaluations[blockIdx.x] = *evaluation;
}



template<class T>
void cudaVersion_evaluationMethodMovingLS(int numOfEvaluationPoints, T* evaluationPoints, int numOfMonomials, T* host_coefficients, T* host_evaluations)
{
	int dev_count(0); cudaGetDeviceCount(&dev_count); cudaDeviceProp dev_prop[dev_count];
	dev_count = get_AllDevicesProperties(dev_prop);
//	printf("number of devices %d\n", dev_count);
//	plotGPUBasicComputeResources(dev_prop[0]);
//	plotGPUMemoryResources(dev_prop[0]);
//	plotGPUPerformanceProperties(dev_prop[0]);
//	plotGPUGridSizeProperties(dev_prop[0]);

	int totalAvailableThreads(0);
	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++)
		totalAvailableThreads += dev_prop[deviceIndex].multiProcessorCount*dev_prop[deviceIndex].maxThreadsPerMultiProcessor;


	T *gpu_coefficients, *gpu_evaluationPoints, *gpu_evaluations;
	size_t sizeOfCoefficientsArray		= sizeof(T)*numOfEvaluationPoints*numOfMonomials;
	size_t sizeOfEvaluationPointsArray	= sizeof(T)*numOfEvaluationPoints;

	
	cudaError_t 	cudaMallocStatus = cudaMalloc( (void**) &gpu_coefficients,	sizeOfCoefficientsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_evaluationPoints,	sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;
			cudaMallocStatus = cudaMalloc( (void**) &gpu_evaluations,	sizeOfEvaluationPointsArray);	if (cudaMallocStatus) std::cout << "ERROR!!" << std::endl;


	cudaMemcpy(gpu_evaluationPoints, evaluationPoints, sizeOfEvaluationPointsArray, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_coefficients, 	 host_coefficients,sizeOfCoefficientsArray, 	cudaMemcpyDeviceToHost);


	evaluationMethodMovingLS_cudaKernel <<< dim3(numOfEvaluationPoints, 1, 1), dim3(numOfMonomials, 1, 1), size_t((numOfMonomials+2)*sizeof(T)) >>> 
	(gpu_evaluationPoints, numOfMonomials, gpu_coefficients, gpu_evaluations);


	cudaMemcpy(host_evaluations, gpu_evaluations, sizeOfEvaluationPointsArray, cudaMemcpyDeviceToHost);
	cudaFree(gpu_coefficients); cudaFree(gpu_evaluationPoints); cudaFree(gpu_evaluations); 
}


template void cudaVersion_evaluationMethodMovingLS(int numOfEvaluationPoints, float* evaluationPoints, int numOfMonomials, float* host_coefficients, float* evaluations);
template void cudaVersion_evaluationMethodMovingLS(int numOfEvaluationPoints, double* evaluationPoints, int numOfMonomials, double* host_coefficients, double* evaluations);
