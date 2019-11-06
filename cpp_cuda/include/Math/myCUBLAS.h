#include<cuda.h>

template <typename T = float, unsigned int blockSize = 256>
__global__ void reduce6(T *g_idata, T *g_odata, unsigned int n)  // required num of blocks = n/2/blockSize, shared memory required for block size number of data
{
	extern __shared__ char sdata[];

	unsigned int tid 	= threadIdx.x;
	unsigned int i 		= blockIdx.x*(blockSize*2) + tid; 	// global index over the data, every block will read twice its size number of data in every iteration
	unsigned int gridSize 	= blockSize*2*gridDim.x; 		// actual number of data that can be read from the running total number of threads in every iteration
	sdata[tid] = 0;

	while ( i < n ) { sdata[tid] += g_idata[i] + g_idata[i+blockSize]; i += gridSize; } //  n/gridsize number of iterations to read all the data
	__syncthreads();


	if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
	if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
	if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }

	if (tid < 32) {
		if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
		if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
		if (blockSize >= 16) sdata[tid] += sdata[tid +  8];
		if (blockSize >=  8) sdata[tid] += sdata[tid +  4];
		if (blockSize >=  4) sdata[tid] += sdata[tid +  2];
		if (blockSize >=  2) sdata[tid] += sdata[tid +  1];
	}
	if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


template<class T>
__global__
void tensorProduct_cudaKernel( T* columnVector, T* rowVector, size_t columnVectorSize, size_t rowVectorSize, T* array )
{

	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *columnVector_local 	= reinterpret_cast<T*>(sharedMemory);
	T *rowVector_local 	= columnVector +columnVectorSize;
	

	if ( threadIdx.y == 0 ) rowVector_local		[ threadIdx.x ] = rowVector	[ threadIdx.x ];
	if ( threadIdx.x == 0 ) columnVector_local	[ threadIdx.y ] = columnVector	[ threadIdx.y ];
	__syncthreads();

	array[ rowVectorSize*threadIdx.y +threadIdx.x ] = columnVector_local[ threadIdx.y ] *rowVector_local[ threadIdx.x ];

}


/*
template<class T>
__global__
void tensorProduct_cudaKernel( T* vector, size_t vectorSize, T* array )
{

	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *vector_local = reinterpret_cast<T*>(sharedMemory);
	

	if ( threadIdx.y == 0 ) vector_local [ threadIdx.x ] = vector [ threadIdx.x ];
	__syncthreads();


	array[ vectorSize*threadIdx.y +threadIdx.x ] = vector_local[ threadIdx.y ] *vector_local[ threadIdx.x ];

}
*/



