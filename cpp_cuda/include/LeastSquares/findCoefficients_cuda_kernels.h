#include<cuda_runtime.h>

template<class T>
__global__
void findCoefficientsLCS_cudaKernel(int *ranges, T span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* coefficients);


template<class T>
__global__
void findCoefficients_cudaKernel(int *ranges, T span, int maxNumOfEffectiveSamplePoints, int numOfStationaryPoints, T* stationaryPoints, int numOfSamplePoints, T* samplePoints, T* sampleValues, int numOfMonomials, T* coefficients);


