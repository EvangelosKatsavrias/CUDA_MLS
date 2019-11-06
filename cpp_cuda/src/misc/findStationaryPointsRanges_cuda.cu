#include"findStationaryPointsRanges_cuda.h"

template<class T>
__global__
void find_samplePointsRanges(int numOfStationaryPoints, T *stationaryPoints, int numOfSamplePoints, T* samplePoints, T span, int *ranges, int* maxNumOfEffectiveSamplePoints)
{
	__shared__ int local_maxRange[1];

	if (threadIdx.x == 0) *local_maxRange = 0;
	T stationaryPoint;

	if (blockIdx.x*blockDim.x+threadIdx.x < numOfStationaryPoints)
	{
		stationaryPoint = stationaryPoints[blockIdx.x*blockDim.x +threadIdx.x];

		int distance = numOfSamplePoints/2, location = distance;
		T upperBound, eq_tol(span/20);

		while (distance > 1)
		{
			upperBound = samplePoints[location]; distance /= 2;
			if ( fabs(upperBound - stationaryPoint) < eq_tol ) break;
			if ( upperBound < stationaryPoint ) location += distance;
			else location -= distance;
		}

		int lowerValue_shift(1), upperValue_shift(1);
		while ( stationaryPoint - samplePoints[location -lowerValue_shift] < span && location -lowerValue_shift > 0) lowerValue_shift += 1;
		while ( samplePoints[location +upperValue_shift] - stationaryPoint < span && numOfSamplePoints -location -upperValue_shift > 1) upperValue_shift += 1;

		__syncthreads();
		atomicMax(local_maxRange, upperValue_shift +lowerValue_shift);

		int index = (blockIdx.x*blockDim.x +threadIdx.x)*2;
		ranges[index] 	= location -lowerValue_shift;
		ranges[index +1]= location +upperValue_shift;
	}

	if (threadIdx.x ==0 && blockIdx.x == 0) *maxNumOfEffectiveSamplePoints =0; __syncthreads();
	if (threadIdx.x ==0 ) atomicMax(maxNumOfEffectiveSamplePoints, *local_maxRange); __syncthreads();
	if (threadIdx.x ==0 && blockIdx.x == 0) (*maxNumOfEffectiveSamplePoints)+=1;
}

template __global__ void find_samplePointsRanges(int numOfStationaryPoints, float *stationaryPoints, int numOfSamplePoints, float *samplePoints, float span, int *ranges, int* maxNumOfEffectiveSamplePoints);
template __global__ void find_samplePointsRanges(int numOfStationaryPoints, double *stationaryPoints, int numOfSamplePoints, double *samplePoints, double span, int *ranges, int* maxNumOfEffectiveSamplePoints);
