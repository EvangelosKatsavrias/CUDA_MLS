#include"findCoefficients_cuda_kernels.h"
#include<stdio.h>
#include"atomicOperations_doubles.h"
#include"WeightingFunction.h"


template<class T>
__global__
void findCoefficientsLCS_cudaKernel(int *ranges, T span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* coefficients)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *local_samplePoints	= reinterpret_cast<T*>(sharedMemory);
	T *local_sampleValues	= local_samplePoints + maxNumOfEffectiveSamplePoints;
	T *Sum_bbT		= local_sampleValues + maxNumOfEffectiveSamplePoints;
	T *Sum_bf		= Sum_bbT + numOfMonomials*numOfMonomials;
	T *b			= Sum_bf +numOfMonomials;

	if (threadIdx.y < (numOfMonomials+1)) Sum_bbT[threadIdx.y*numOfMonomials+threadIdx.x] = 0;

	int lowerSamplePointIndex = ranges[blockIdx.x*2], numOfEffectiveSamplePoints = ranges[blockIdx.x*2+1] -lowerSamplePointIndex +1;

	if (threadIdx.x == 0 && threadIdx.y < numOfEffectiveSamplePoints )
	{
		local_samplePoints[threadIdx.y] = samplePoints[lowerSamplePointIndex +threadIdx.y];
		local_sampleValues[threadIdx.y] = sampleValues[lowerSamplePointIndex +threadIdx.y];
	}
	__syncthreads();



	T stationaryPoint(stationaryPoints[blockIdx.x]), reciprocal_span(1/span);
	if ( threadIdx.y < numOfEffectiveSamplePoints )
	{
		T theta(wendlandDistribution(reciprocal_span, fabs(local_samplePoints[threadIdx.y]- stationaryPoint))); 
		b[threadIdx.x +threadIdx.y*numOfMonomials] = pow(local_samplePoints[threadIdx.y] -stationaryPoint, threadIdx.x);
		__syncthreads();
		for (int index = 0; index < numOfMonomials; index++) atomicAdd(Sum_bbT+(index*numOfMonomials +threadIdx.x), theta*b[index+threadIdx.y*numOfMonomials]*b[threadIdx.x+threadIdx.y*numOfMonomials]);
		atomicAdd(Sum_bf+threadIdx.x, theta*b[threadIdx.x +threadIdx.y*numOfMonomials]*local_sampleValues[threadIdx.y]);
	}
	__syncthreads();


	if (threadIdx.x == 0 && threadIdx.y == 0)
	{
		for (int j = 0; j<(numOfMonomials-1); j++)
		{
			T inv = 1/Sum_bbT[j*numOfMonomials+j];
			for (int i = j+1; i<numOfMonomials; i++)
			{
				T lambda = Sum_bbT[j+i*numOfMonomials]*inv;
				for (int k = j+1; k<numOfMonomials; k++)
					Sum_bbT[k+i*numOfMonomials] = Sum_bbT[k+i*numOfMonomials] -lambda*Sum_bbT[k+j*numOfMonomials];
				Sum_bf[i] = Sum_bf[i] - lambda*Sum_bf[j];
			}
		}

		for (int i=(numOfMonomials-1); i>-1; i--)
		{
			for (int j=(numOfMonomials-1); j>i; j--)
				Sum_bf[i] = Sum_bf[i] - Sum_bbT[j+i*numOfMonomials]*Sum_bf[j];
			Sum_bf[i] = Sum_bf[i]/Sum_bbT[i+i*numOfMonomials];
		}
	}


	if (threadIdx.y == 0) coefficients[blockIdx.x*numOfMonomials +threadIdx.x] = Sum_bf[threadIdx.x];
}


template __global__ void findCoefficientsLCS_cudaKernel(int *ranges, float span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, float* stationaryPoints, int numOfSamplePoints, float* samplePoints, float* sampleValues, int numOfMonomials, float* coefficients);
template __global__ void findCoefficientsLCS_cudaKernel(int *ranges, double span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, double* stationaryPoints, int numOfSamplePoints, double* samplePoints, double* sampleValues, int numOfMonomials, double* coefficients);


template<class T>
__global__
void findCoefficients_cudaKernel(int *ranges, T span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* coefficients)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *local_samplePoints	= reinterpret_cast<T*>(sharedMemory);
	T *local_sampleValues	= local_samplePoints + maxNumOfEffectiveSamplePoints;
	T *Sum_bbT		= local_sampleValues + maxNumOfEffectiveSamplePoints;
	T *Sum_bf		= Sum_bbT + numOfMonomials*numOfMonomials;
	T *b			= Sum_bf +numOfMonomials;

	if (threadIdx.y < (numOfMonomials+1)) Sum_bbT[threadIdx.y*numOfMonomials+threadIdx.x] = 0;

	int lowerSamplePointIndex = ranges[blockIdx.x*2], numOfEffectiveSamplePoints = ranges[blockIdx.x*2+1] -lowerSamplePointIndex +1;

	if (threadIdx.x == 0 && threadIdx.y < numOfEffectiveSamplePoints )
	{
		local_samplePoints[threadIdx.y] = samplePoints[lowerSamplePointIndex +threadIdx.y];
		local_sampleValues[threadIdx.y] = sampleValues[lowerSamplePointIndex +threadIdx.y];
	}
	__syncthreads();



	T stationaryPoint(stationaryPoints[blockIdx.x]), reciprocal_span(1/span);
	if ( threadIdx.y < numOfEffectiveSamplePoints )
	{
		T theta(wendlandDistribution(reciprocal_span, fabs(local_samplePoints[threadIdx.y]- stationaryPoint))); 
		b[threadIdx.x +threadIdx.y*numOfMonomials] = pow(local_samplePoints[threadIdx.y], threadIdx.x);
		__syncthreads();
		for (int index = 0; index < numOfMonomials; index++) atomicAdd(Sum_bbT+(index*numOfMonomials +threadIdx.x), theta*b[index+threadIdx.y*numOfMonomials]*b[threadIdx.x+threadIdx.y*numOfMonomials]);
		atomicAdd(Sum_bf+threadIdx.x, theta*b[threadIdx.x +threadIdx.y*numOfMonomials]*local_sampleValues[threadIdx.y]);
	}
	__syncthreads();


	if (threadIdx.x == 0 && threadIdx.y == 0)
	{
		for (int j = 0; j<(numOfMonomials-1); j++)
		{
			T inv = 1/Sum_bbT[j*numOfMonomials+j];
			for (int i = j+1; i<numOfMonomials; i++)
			{
				T lambda = Sum_bbT[j+i*numOfMonomials]*inv;
				for (int k = j+1; k<numOfMonomials; k++)
					Sum_bbT[k+i*numOfMonomials] = Sum_bbT[k+i*numOfMonomials] -lambda*Sum_bbT[k+j*numOfMonomials];
				Sum_bf[i] = Sum_bf[i] - lambda*Sum_bf[j];
			}
		}

		for (int i=(numOfMonomials-1); i>-1; i--)
		{
			for (int j=(numOfMonomials-1); j>i; j--)
				Sum_bf[i] = Sum_bf[i] - Sum_bbT[j+i*numOfMonomials]*Sum_bf[j];
			Sum_bf[i] = Sum_bf[i]/Sum_bbT[i+i*numOfMonomials];
		}
	}


	if (threadIdx.y == 0) coefficients[blockIdx.x*numOfMonomials +threadIdx.x] = Sum_bf[threadIdx.x];
}


template __global__ void findCoefficients_cudaKernel(int *ranges, float span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, float* stationaryPoints, int numOfSamplePoints, float* samplePoints, float* sampleValues, int numOfMonomials, float* coefficients);
template __global__ void findCoefficients_cudaKernel(int *ranges, double span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, double* stationaryPoints, int numOfSamplePoints, double* samplePoints, double* sampleValues, int numOfMonomials, double* coefficients);

