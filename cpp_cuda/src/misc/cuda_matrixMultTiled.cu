#include"cuda_matrixMultTiled.h"
#include<stdio.h>

template<class T>
__global__
void tiledMatrixMulKernel_Const(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, T *d_M, T *d_N, T *d_P)
{
	__shared__ T Mds[TILE_WIDTH][TILE_WIDTH], Nds[TILE_WIDTH][TILE_WIDTH];

	T Pvalue(0);

	int M_Row = blockIdx.y*blockDim.y +threadIdx.y, N_Col = blockIdx.x*blockDim.x +threadIdx.x;
	int numberOfTiles = (sizeOfContractedDimension%TILE_WIDTH ? sizeOfContractedDimension/TILE_WIDTH + 1:sizeOfContractedDimension/TILE_WIDTH);

	for (int tileIndex = 0; tileIndex < numberOfTiles; tileIndex++)
	{
		Mds[threadIdx.y][threadIdx.x] = d_M[M_Row *sizeOfContractedDimension +tileIndex*TILE_WIDTH +threadIdx.x];
		Nds[threadIdx.y][threadIdx.x] = d_N[(tileIndex*TILE_WIDTH +threadIdx.y) *numOfColumnsOfRightMatrix +N_Col];
		__syncthreads();

		for (int k = 0; k < TILE_WIDTH; k++) Pvalue += Mds[threadIdx.y][k] *Nds[k][threadIdx.x];
		__syncthreads();

		if (threadIdx.x < numOfColumnsOfRightMatrix && threadIdx.y < numOfRowsOfLeftMatrix)
			d_P[M_Row*numOfColumnsOfRightMatrix + N_Col] = Pvalue;
	}

//	printf("tileIndex = %d, block %d, %d, map %d --> %d, %d, %f\n", tileIndex, blockIdx.y, blockIdx.x, Row*sizeOfContractedDimension+tileIndex*TILE_WIDTH+threadIdx.x, threadIdx.y, threadIdx.x, Mds[threadIdx.y][threadIdx.x]);
//	printf("block %d, %d, map %d --> %d, %d, %f\n", blockIdx.y, blockIdx.x, (tileIndex*TILE_WIDTH+threadIdx.y)*numOfColumnsOfRightMatrix+Col, threadIdx.y, threadIdx.x, Nds[threadIdx.y][threadIdx.x]);

}

template __global__ void tiledMatrixMulKernel_Const(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, float *d_M, float *d_N, float *d_P);
template __global__ void tiledMatrixMulKernel_Const(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, double *d_M, double *d_N, double *d_P);


template<class T>
__global__
void tiledMatrixMulKernel_Dynam(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, T *d_M, T *d_N, T *d_P)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* Mds = reinterpret_cast<T*>(sharedMemory);

	T *Nds(Mds +blockDim.x*blockDim.y);
	T Pvalue(0), *M(Mds+threadIdx.y*blockDim.x), *N(Nds+threadIdx.x*blockDim.y);

	int M_Row = blockIdx.y*blockDim.y +threadIdx.y, N_Col = blockIdx.x*blockDim.x +threadIdx.x;
	int numberOfTiles = (sizeOfContractedDimension%blockDim.x ? sizeOfContractedDimension/blockDim.x + 1:sizeOfContractedDimension/blockDim.x);

	for (int tileIndex = 0; tileIndex < numberOfTiles; tileIndex++)
	{
		Mds[threadIdx.y*blockDim.x +threadIdx.x] = d_M[M_Row*sizeOfContractedDimension +tileIndex*blockDim.x +threadIdx.x];
		Nds[threadIdx.x*blockDim.y +threadIdx.y] = d_N[(tileIndex*blockDim.y +threadIdx.y)*numOfColumnsOfRightMatrix + N_Col];
		__syncthreads();

		for (int k = 0; k < blockDim.x; k++) Pvalue +=  M[k]*N[k];
		__syncthreads();
	
		if (threadIdx.x < numOfColumnsOfRightMatrix && threadIdx.y < numOfRowsOfLeftMatrix)
			d_P[M_Row*numOfColumnsOfRightMatrix + N_Col] = Pvalue;
	}

//	printf("tileIndex = %d, block %d, %d, %d --> %d, %f\n", tileIndex, blockIdx.y, blockIdx.x, M_Row*sizeOfContractedDimension +tileIndex*blockDim.x+threadIdx.x, threadIdx.y*blockDim.x+threadIdx.x, Mds[threadIdx.y*blockDim.x+threadIdx.x]);
//	printf("block %d, %d, %d --> %d, %f\n", blockIdx.y, blockIdx.x, (tileIndex*blockDim.y+threadIdx.y)*numOfColumnsOfRightMatrix +N_Col, threadIdx.x*blockDim.y+threadIdx.y, Nds[threadIdx.x*blockDim.y+threadIdx.y]);
//	printf("tileIndex %d, block %d, %d, thread %d, %d, Pvalue %f\n", tileIndex, blockIdx.y, blockIdx.x, threadIdx.y, threadIdx.x, Pvalue);
}

template __global__ void tiledMatrixMulKernel_Dynam(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, float *d_M, float *d_N, float *d_P);
template __global__ void tiledMatrixMulKernel_Dynam(int numOfRowsOfLeftMatrix, int sizeOfContractedDimension, int numOfColumnsOfRightMatrix, double *d_M, double *d_N, double *d_P);

template<class T>
void cuda_matrixMultConstTiled(int m, int n, int k, T *h_M, T *h_N, T *h_P)
{
	dim3 blockSize(TILE_WIDTH, TILE_WIDTH, 1);
	dim3 gridSize(ceil((float)n/blockSize.x), ceil((float)m/blockSize.y), 1);
	printf("block size, %d, %d,  grid size, %d, %d\n", blockSize.x, blockSize.y, gridSize.x, gridSize.y);

	size_t numOfThreads = (gridSize.x*blockSize.x)*(gridSize.y*blockSize.y);
	size_t sizeOf_M = m*k*sizeof(T),		sizeOf_N = k*n*sizeof(T),		sizeOf_P = m*n*sizeof(T);
	size_t sizeOf_dM = numOfThreads*sizeof(T),	sizeOf_dN = numOfThreads*sizeof(T),	sizeOf_dP = numOfThreads*sizeof(T);

	T *d_M, *d_N, *d_P;
	cudaError_t	cudaMallocStatus = cudaMalloc((void **) &d_M, sizeOf_dM);
  			cudaMallocStatus = cudaMalloc((void **) &d_N, sizeOf_dN);
  			cudaMallocStatus = cudaMalloc((void **) &d_P, sizeOf_dP);

	cudaMemcpy(d_M, h_M, sizeOf_M, cudaMemcpyHostToDevice);
	cudaMemcpy(d_N, h_N, sizeOf_N, cudaMemcpyHostToDevice);

	cudaMemset(d_M+m*k, 0, (numOfThreads -m*k)*sizeof(T) );	cudaMemset(d_N+n*k, 0, (numOfThreads -n*k)*sizeof(T) );
	tiledMatrixMulKernel_Const<<<gridSize, blockSize>>>(m, k, n, d_M, d_N, d_P);
	cudaMemcpy(h_P, d_P, sizeOf_P, cudaMemcpyDeviceToHost);

	cudaFree(d_N);  cudaFree(d_M); cudaFree(d_P);
}

template void cuda_matrixMultConstTiled(int m, int n, int k, float *h_M, float *h_N, float *h_P);
template void cuda_matrixMultConstTiled(int m, int n, int k, double *h_M, double *h_N, double *h_P);


template<class T>
void cuda_matrixMultTiled(int m, int n, int k, T *h_M, T *h_N, T *h_P)
{
	dim3 blockSize(TILE_WIDTH, TILE_WIDTH, 1);
	dim3 gridSize(ceil((float)n/blockSize.x), ceil((float)m/blockSize.y), 1);
	printf("block size, %d, %d,  grid size, %d, %d\n", blockSize.x, blockSize.y, gridSize.x, gridSize.y);

	size_t numOfThreads = (gridSize.x*blockSize.x)*(gridSize.y*blockSize.y);
	size_t sizeOf_M = m*k*sizeof(T),		sizeOf_N = k*n*sizeof(T),		sizeOf_P = m*n*sizeof(T);
	size_t sizeOf_dM = numOfThreads*sizeof(T),	sizeOf_dN = numOfThreads*sizeof(T),	sizeOf_dP = numOfThreads*sizeof(T);

	T *d_M, *d_N, *d_P;
	cudaError_t	cudaMallocStatus = cudaMalloc((void **) &d_M, sizeOf_dM);
  			cudaMallocStatus = cudaMalloc((void **) &d_N, sizeOf_dN);
  			cudaMallocStatus = cudaMalloc((void **) &d_P, sizeOf_dP);

	cudaMemcpy(d_M, h_M, sizeOf_M, cudaMemcpyHostToDevice);
	cudaMemcpy(d_N, h_N, sizeOf_N, cudaMemcpyHostToDevice);

	cudaMemset(d_M+m*k, 0, (numOfThreads -m*k)*sizeof(T) );	cudaMemset(d_N+n*k, 0, (numOfThreads -n*k)*sizeof(T) );
	tiledMatrixMulKernel_Dynam<<<gridSize, blockSize, 2*TILE_WIDTH*TILE_WIDTH*sizeof(T)>>>(m, k, n, d_M, d_N, d_P);
	cudaMemcpy(h_P, d_P, sizeOf_P, cudaMemcpyDeviceToHost);

	cudaFree(d_N);  cudaFree(d_M); cudaFree(d_P);
}

template void cuda_matrixMultTiled(int m, int n, int k, float *h_M, float *h_N, float *h_P);
template void cuda_matrixMultTiled(int m, int n, int k, double *h_M, double *h_N, double *h_P);

