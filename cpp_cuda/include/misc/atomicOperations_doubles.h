#ifndef CUDADOUBLESATOMICOPERATIONS
#define CUDADOUBLESATOMICOPERATIONS

#include<cuda.h>
#include<cuda_runtime.h>

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val);
#endif

#endif
