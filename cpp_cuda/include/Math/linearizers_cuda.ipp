#include<cuda.h>

template __device__ int linearize_rowWiseZeroBased_2d(int, int, int);
template __device__ int linearize_rowWiseZeroBased_3d(int, int, int, int, int);


template __device__ int linearizeBlockId_in3dGrid<int>();
template __device__ int linearizeBlockId_in2dGrid<int>();
template __device__ int linearizeThreadId_in3dBlock<int>();
template __device__ int linearizeThreadId_in2dBlock<int>();
template __device__ int threadsShift_in3dGrid<int>();
template __device__ int threadsShift_in2dGrid<int>();
template __device__ int linearizeThreadId_in3dGrid<int>();
template __device__ int linearizeThreadId_in2dGrid<int>();

