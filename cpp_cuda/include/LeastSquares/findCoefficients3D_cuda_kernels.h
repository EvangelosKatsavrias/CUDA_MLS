#include"WeightingFunction.h"


template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel_0deg ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy, T* Sum_bfz );


template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy, T* Sum_bfz );


template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel_0deg2 ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf);


template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel2 ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf);

