#include"findCoefficients_cuda_kernels.h"
#include<iostream>
#include<fstream>
#include"weightedLeastSquares_cuda.h"
#include"gpuDeviceProperties.h"
#include"linearizers_cuda.h"
#include"findStationaryPointsRanges_cuda.h"
#include<typeinfo>
#include"vectorPlotOutStreamer.h"
#include"conjugateGradient.h"
#include"fileFunctions.h"
#include"myBLAS.h"
#include<cuda.h>
//#include<cudaProfiler.h>
//#include<cublas_v2.h>



template<class T>
void cudaVersion_findCoefficientsLCS(T span, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* host_coefficients)
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
		

	T *gpu_coefficients, *gpu_stationaryPoints, *gpu_samplePoints, *gpu_sampleValues;
	size_t sizeOfCoefficientsArray 		= sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;


	cudaError_t 	cudaStatus = cudaMalloc((void**) &gpu_coefficients, 	sizeOfCoefficientsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
			cudaStatus = cudaMalloc((void**) &gpu_stationaryPoints, 	sizeOfStationaryPointsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
			cudaStatus = cudaMalloc((void**) &gpu_samplePoints, 	sizeOfSamplePointsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
			cudaStatus = cudaMalloc((void**) &gpu_sampleValues, 	sizeOfSamplePointsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	cudaMemcpy(gpu_stationaryPoints, stationaryPoints, sizeOfStationaryPointsArray, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_samplePoints, 	 samplePoints,	   sizeOfSamplePointsArray, 	cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_sampleValues, 	 sampleValues,	   sizeOfSamplePointsArray, 	cudaMemcpyHostToDevice);



	int *gpu_ranges, maxNumOfEffectiveSamplePoints, *gpu_maxNumOfEffectiveSamplePoints;
	cudaStatus = cudaMalloc( (void**) &gpu_ranges, 2*numOfStationaryPoints*sizeof(int));
	cudaStatus = cudaMalloc( (void**) &gpu_maxNumOfEffectiveSamplePoints, sizeof(int));
	int maxBlockSize = min(numOfStationaryPoints, 256);

	find_samplePointsRanges <<< dim3( ceil(numOfStationaryPoints/maxBlockSize), 1, 1), dim3(maxBlockSize, 1, 1) >>> 
	(numOfStationaryPoints, gpu_stationaryPoints, numOfSamplePoints, gpu_samplePoints, span, gpu_ranges, gpu_maxNumOfEffectiveSamplePoints);

	cudaMemcpy(&maxNumOfEffectiveSamplePoints, gpu_maxNumOfEffectiveSamplePoints, sizeof(int), cudaMemcpyDeviceToHost);
//	printf("Max number of effective sample points %d\n", maxNumOfEffectiveSamplePoints);


	size_t sharedMemorySize = ((2+numOfMonomials)*maxNumOfEffectiveSamplePoints+numOfMonomials*(numOfMonomials+1)+1)*sizeof(T);
	findCoefficientsLCS_cudaKernel <<< dim3(numOfStationaryPoints, 1, 1), dim3(numOfMonomials, maxNumOfEffectiveSamplePoints, 1), sharedMemorySize >>> 
	(gpu_ranges, span, maxNumOfEffectiveSamplePoints, numOfStationaryPoints, gpu_stationaryPoints, numOfSamplePoints, gpu_samplePoints, gpu_sampleValues, numOfMonomials, gpu_coefficients);


	cudaMemcpy(host_coefficients, gpu_coefficients, sizeOfCoefficientsArray, cudaMemcpyDeviceToHost);
	cudaFree(gpu_coefficients); cudaFree(gpu_stationaryPoints);
	cudaFree(gpu_samplePoints); cudaFree(gpu_sampleValues); cudaFree(gpu_ranges); cudaFree(gpu_maxNumOfEffectiveSamplePoints);
}

template void cudaVersion_findCoefficientsLCS(float span, int, float* stationaryPoints, int numOfSamplePoints, float* samplePoints, float* sampleValues, int numOfMonomials, float* coefficients);
template void cudaVersion_findCoefficientsLCS(double span, int, double* stationaryPoints, int numOfSamplePoints, double* samplePoints, double* sampleValues, int numOfMonomials, double* coefficients);


template<class T>
void cudaVersion_findCoefficients(T span, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* host_coefficients)
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
		

	T *gpu_coefficients, *gpu_stationaryPoints, *gpu_samplePoints, *gpu_sampleValues;
	size_t sizeOfCoefficientsArray 		= sizeof(T)*numOfStationaryPoints*numOfMonomials;
	size_t sizeOfStationaryPointsArray 	= sizeof(T)*numOfStationaryPoints;
	size_t sizeOfSamplePointsArray		= sizeof(T)*numOfSamplePoints;


	cudaError_t 	cudaStatus = cudaMalloc((void**) &gpu_coefficients, 	sizeOfCoefficientsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
			cudaStatus = cudaMalloc((void**) &gpu_stationaryPoints, 	sizeOfStationaryPointsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
			cudaStatus = cudaMalloc((void**) &gpu_samplePoints, 	sizeOfSamplePointsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;
			cudaStatus = cudaMalloc((void**) &gpu_sampleValues, 	sizeOfSamplePointsArray);	if (cudaStatus) std::cout << "ERROR!!" << std::endl;


	cudaMemcpy(gpu_stationaryPoints, stationaryPoints, sizeOfStationaryPointsArray, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_samplePoints, 	 samplePoints,	   sizeOfSamplePointsArray, 	cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_sampleValues, 	 sampleValues,	   sizeOfSamplePointsArray, 	cudaMemcpyHostToDevice);



	int *gpu_ranges, maxNumOfEffectiveSamplePoints, *gpu_maxNumOfEffectiveSamplePoints;
	cudaStatus = cudaMalloc( (void**) &gpu_ranges, 2*numOfStationaryPoints*sizeof(int));
	cudaStatus = cudaMalloc( (void**) &gpu_maxNumOfEffectiveSamplePoints, sizeof(int));
	int maxBlockSize = min(numOfStationaryPoints, 256);

	find_samplePointsRanges <<< dim3( ceil(numOfStationaryPoints/maxBlockSize), 1, 1), dim3(maxBlockSize, 1, 1) >>> 
	(numOfStationaryPoints, gpu_stationaryPoints, numOfSamplePoints, gpu_samplePoints, span, gpu_ranges, gpu_maxNumOfEffectiveSamplePoints);

	cudaMemcpy(&maxNumOfEffectiveSamplePoints, gpu_maxNumOfEffectiveSamplePoints, sizeof(int), cudaMemcpyDeviceToHost);
//	printf("Max number of effective sample points %d\n", maxNumOfEffectiveSamplePoints);


	size_t sharedMemorySize = ((2+numOfMonomials)*maxNumOfEffectiveSamplePoints+numOfMonomials*(numOfMonomials+1)+1)*sizeof(T);
	findCoefficients_cudaKernel <<< dim3(numOfStationaryPoints, 1, 1), dim3(numOfMonomials, maxNumOfEffectiveSamplePoints, 1), sharedMemorySize >>> 
	(gpu_ranges, span, maxNumOfEffectiveSamplePoints, numOfStationaryPoints, gpu_stationaryPoints, numOfSamplePoints, gpu_samplePoints, gpu_sampleValues, numOfMonomials, gpu_coefficients);


	cudaMemcpy(host_coefficients, gpu_coefficients, sizeOfCoefficientsArray, cudaMemcpyDeviceToHost);
	cudaFree(gpu_coefficients); cudaFree(gpu_stationaryPoints);
	cudaFree(gpu_samplePoints); cudaFree(gpu_sampleValues); cudaFree(gpu_ranges); cudaFree(gpu_maxNumOfEffectiveSamplePoints);
}


template void cudaVersion_findCoefficients(float span, int, float* stationaryPoints, int numOfSamplePoints, float* samplePoints, float* sampleValues, int numOfMonomials, float* coefficients);
template void cudaVersion_findCoefficients(double span, int, double* stationaryPoints, int numOfSamplePoints, double* samplePoints, double* sampleValues, int numOfMonomials, double* coefficients);

