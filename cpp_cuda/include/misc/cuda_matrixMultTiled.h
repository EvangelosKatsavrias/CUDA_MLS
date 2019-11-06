#include<cuda.h>

#define TILE_WIDTH 16

template<class T> __global__
void tiledMatrixMulKernel_Const(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, T *d_M, T *d_N, T *d_P);

template<class T> __global__
void tiledMatrixMulKernel_Dynam(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, T *d_M, T *d_N, T *d_P);

template<class T>
void cuda_matrixMultTiled(int m, int n, int k, T *h_M, T *h_N, T *h_P);
template<class T>
void cuda_matrixMultConstTiled(int m, int n, int k, T *h_M, T *h_N, T *h_P);

