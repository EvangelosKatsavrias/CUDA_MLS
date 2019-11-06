#include"findCoefficients2D_cuda_kernels.h"
#include<stdio.h>
#include"atomicOperations_doubles.h"
#include"myBLAS.h"


template<class T>
__global__
void findCoefficientsLCS2D_cudaKernel_0deg ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.x;

	T theta;

	size_t numOfTiles = (numOfSamplePoints +blockDim.x -1)/blockDim.x;
	for (size_t tile = 0; tile < numOfTiles; tile++ ) {

		size_t pIndex = threadIdx.x +tile*blockDim.x;
		if ( pIndex < numOfSamplePoints ) {
			T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[ blockIdx.x ];
			T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[ blockIdx.x ];
			theta = weightingFunction( sqrt(dx*dx + dy*dy) );

			bTheta[ threadIdx.x ] = theta*sampleValues_x [ pIndex ];
			b_sh  [ threadIdx.x ] = theta*sampleValues_y [ pIndex ];
		}
		__syncthreads();


		size_t 	sind 	= threadIdx.x;
		int 	resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[0] += bTheta[sind]; b_sh[0] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 && tile == 0 ) {	Sum_bfx[blockIdx.x ] = bTheta[ 0 ];
							Sum_bfy[blockIdx.x ] = b_sh  [ 0 ]; }
		if ( threadIdx.x == 0 && tile  > 0 ) {	Sum_bfx[blockIdx.x ] += bTheta[ 0 ];
							Sum_bfy[blockIdx.x ] += b_sh  [ 0 ]; }
		__syncthreads();



		bTheta[ threadIdx.x ] = theta;
		__syncthreads();


		sind 	= threadIdx.x;
		resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		newReductionSize = 1024;
	
		while ( resid > 0 ) {
	
			if (resid >=1024){ 	if ( threadIdx.x < 512) bTheta[ sind ] 	+= bTheta[ sind +512]; }
			else newReductionSize = 512; __syncthreads();
			if (resid >=512) {	if ( threadIdx.x < 256) bTheta[ sind ] 	+= bTheta[ sind +256]; }
			else newReductionSize = 256; __syncthreads();
			if (resid >=256) {	if ( threadIdx.x < 128) bTheta[ sind ] 	+= bTheta[ sind +128]; }
			else newReductionSize = 128; __syncthreads();
			if (resid >=128) {	if ( threadIdx.x <  64) bTheta[ sind ]	+= bTheta[ sind +64]; }
			else newReductionSize = 64; __syncthreads();
			if (resid >= 64) {	if ( threadIdx.x <  32) bTheta[ sind ] 	+= bTheta[ sind +32]; }
			else newReductionSize = 32; __syncthreads();
			if (resid >= 32) {	if ( threadIdx.x <  16) bTheta[ sind ] 	+= bTheta[ sind +16]; }
			else newReductionSize = 16; __syncthreads();
			if (resid >= 16) {	if ( threadIdx.x <   8) bTheta[ sind ] 	+= bTheta[ sind +8]; }
			else newReductionSize = 8; __syncthreads();
			if (resid >= 8 ) {	if ( threadIdx.x <   4) bTheta[ sind ] 	+= bTheta[ sind +4]; }
			else newReductionSize = 4; __syncthreads();
			if (resid >= 4 ) {	if ( threadIdx.x <   2) bTheta[ sind ] 	+= bTheta[ sind +2]; }
			else newReductionSize = 2; __syncthreads();
			if (resid >= 2 ) {	if ( threadIdx.x == 0 ) bTheta[ sind ] 	+= bTheta[ sind +1]; }
			else newReductionSize = 1; __syncthreads();
	
			if ( resid < blockDim.x && threadIdx.x == 0 ) bTheta[0] += bTheta[sind];
			__syncthreads();
	
			resid -= newReductionSize; sind += newReductionSize;
		}
	

		if ( threadIdx.x == 0 && tile == 0 ) Sum_bbT[ blockIdx.x ]  = bTheta[ 0 ];
		if ( threadIdx.x == 0 && tile  > 0 ) Sum_bbT[ blockIdx.x ] += bTheta[ 0 ];

		__syncthreads();
	
	}

}

template __global__ void findCoefficientsLCS2D_cudaKernel_0deg ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bfx, float* Sum_bfy);
template __global__ void findCoefficientsLCS2D_cudaKernel_0deg ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bfx, double* Sum_bfy);


/*

template<class T>
__global__
void findCoefficientsLCS2D_cudaKernel_1deg ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.x;

	T theta;

	size_t numOfTiles = (numOfSamplePoints +blockDim.x -1)/blockDim.x;
	for (size_t tile = 0; tile < numOfTiles; tile++ ) {

		size_t pIndex = threadIdx.x +tile*blockDim.x;
		if ( pIndex < numOfSamplePoints ) {
			T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[ blockIdx.x ];
			T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[ blockIdx.x ];
			theta = weightingFunction( sqrt(dx*dx + dy*dy) );

			bTheta[ threadIdx.x ] = theta*sampleValues_x [ pIndex ];
			b_sh  [ threadIdx.x ] = theta*sampleValues_y [ pIndex ];
		}
		__syncthreads();


		size_t 	sind 	= threadIdx.x;
		int 	resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[0] += bTheta[sind]; b_sh[0] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 && tile == 0 ) {	Sum_bfx[blockIdx.x ] = bTheta[ 0 ];
							Sum_bfy[blockIdx.x ] = b_sh  [ 0 ]; }
		if ( threadIdx.x == 0 && tile  > 0 ) {	Sum_bfx[blockIdx.x ] += bTheta[ 0 ];
							Sum_bfy[blockIdx.x ] += b_sh  [ 0 ]; }
		__syncthreads();



		bTheta[ threadIdx.x ] = theta;
		__syncthreads();


		sind 	= threadIdx.x;
		resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		newReductionSize = 1024;
	
		while ( resid > 0 ) {
	
			if (resid >=1024){ 	if ( threadIdx.x < 512) bTheta[ sind ] 	+= bTheta[ sind +512]; }
			else newReductionSize = 512; __syncthreads();
			if (resid >=512) {	if ( threadIdx.x < 256) bTheta[ sind ] 	+= bTheta[ sind +256]; }
			else newReductionSize = 256; __syncthreads();
			if (resid >=256) {	if ( threadIdx.x < 128) bTheta[ sind ] 	+= bTheta[ sind +128]; }
			else newReductionSize = 128; __syncthreads();
			if (resid >=128) {	if ( threadIdx.x <  64) bTheta[ sind ]	+= bTheta[ sind +64]; }
			else newReductionSize = 64; __syncthreads();
			if (resid >= 64) {	if ( threadIdx.x <  32) bTheta[ sind ] 	+= bTheta[ sind +32]; }
			else newReductionSize = 32; __syncthreads();
			if (resid >= 32) {	if ( threadIdx.x <  16) bTheta[ sind ] 	+= bTheta[ sind +16]; }
			else newReductionSize = 16; __syncthreads();
			if (resid >= 16) {	if ( threadIdx.x <   8) bTheta[ sind ] 	+= bTheta[ sind +8]; }
			else newReductionSize = 8; __syncthreads();
			if (resid >= 8 ) {	if ( threadIdx.x <   4) bTheta[ sind ] 	+= bTheta[ sind +4]; }
			else newReductionSize = 4; __syncthreads();
			if (resid >= 4 ) {	if ( threadIdx.x <   2) bTheta[ sind ] 	+= bTheta[ sind +2]; }
			else newReductionSize = 2; __syncthreads();
			if (resid >= 2 ) {	if ( threadIdx.x == 0 ) bTheta[ sind ] 	+= bTheta[ sind +1]; }
			else newReductionSize = 1; __syncthreads();
	
			if ( resid < blockDim.x && threadIdx.x == 0 ) bTheta[0] += bTheta[sind];
			__syncthreads();
	
			resid -= newReductionSize; sind += newReductionSize;
		}
	

		if ( threadIdx.x == 0 && tile == 0 ) Sum_bbT[ blockIdx.x ]  = bTheta[ 0 ];
		if ( threadIdx.x == 0 && tile  > 0 ) Sum_bbT[ blockIdx.x ] += bTheta[ 0 ];

		__syncthreads();
	
	}

}

template __global__ void findCoefficientsLCS2D_cudaKernel_1deg ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bfx, float* Sum_bfy);
template __global__ void findCoefficientsLCS2D_cudaKernel_1deg ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bfx, double* Sum_bfy);

*/



template<class T>
__global__
void findCoefficientsLCS2D_cudaKernel_0deg2 ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, size_t numOfStoredBlocks, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.x;

	int pIndex = blockIdx.x*blockDim.x + threadIdx.x;
	{

		T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[blockIdx.y];
		T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[blockIdx.y];

		T theta = weightingFunction( sqrt(dx*dx + dy*dy) );

		bTheta[ threadIdx.x ] = theta*sampleValues_x [ pIndex ];
		b_sh  [ threadIdx.x ] = theta*sampleValues_y [ pIndex ];
		__syncthreads();


		size_t 	sind 	= threadIdx.x;
		int 	resid 	= blockDim.x;
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[0] += bTheta[sind]; b_sh[0] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 ) Sum_bf[blockIdx.x*2 +2*numOfStoredBlocks*blockIdx.y ] 		+= bTheta[ 0 ];
		if ( threadIdx.x == 0 )	Sum_bf[blockIdx.x*2 +2*numOfStoredBlocks*blockIdx.y +1 ] 	+= b_sh  [ 0 ]; 
		__syncthreads();
	

		bTheta[ threadIdx.x ] = theta;
		__syncthreads();
	}


	size_t 	sind 	= threadIdx.x;
	int 	resid 	= blockDim.x;
	size_t 	newReductionSize(1024);

	while ( resid > 0 ) {

		if (resid >=1024){ 	if ( threadIdx.x < 512) bTheta[ sind ] 	+= bTheta[ sind +512]; }
		else newReductionSize = 512; __syncthreads();
		if (resid >=512) {	if ( threadIdx.x < 256) bTheta[ sind ] 	+= bTheta[ sind +256]; }
		else newReductionSize = 256; __syncthreads();
		if (resid >=256) {	if ( threadIdx.x < 128) bTheta[ sind ] 	+= bTheta[ sind +128]; }
		else newReductionSize = 128; __syncthreads();
		if (resid >=128) {	if ( threadIdx.x <  64) bTheta[ sind ]	+= bTheta[ sind +64]; }
		else newReductionSize = 64; __syncthreads();
		if (resid >= 64) {	if ( threadIdx.x <  32) bTheta[ sind ] 	+= bTheta[ sind +32]; }
		else newReductionSize = 32; __syncthreads();
		if (resid >= 32) {	if ( threadIdx.x <  16) bTheta[ sind ] 	+= bTheta[ sind +16]; }
		else newReductionSize = 16; __syncthreads();
		if (resid >= 16) {	if ( threadIdx.x <   8) bTheta[ sind ] 	+= bTheta[ sind +8]; }
		else newReductionSize = 8; __syncthreads();
		if (resid >= 8 ) {	if ( threadIdx.x <   4) bTheta[ sind ] 	+= bTheta[ sind +4]; }
		else newReductionSize = 4; __syncthreads();
		if (resid >= 4 ) {	if ( threadIdx.x <   2) bTheta[ sind ] 	+= bTheta[ sind +2]; }
		else newReductionSize = 2; __syncthreads();
		if (resid >= 2 ) {	if ( threadIdx.x == 0 ) bTheta[ sind ] 	+= bTheta[ sind +1]; }
		else newReductionSize = 1; __syncthreads();

		if ( resid < blockDim.x && threadIdx.x == 0 ) bTheta[0] += bTheta[sind];
		__syncthreads();

		resid -= newReductionSize; sind += newReductionSize;
	}

	if ( threadIdx.x == 0 ) Sum_bbT[ blockIdx.x+numOfStoredBlocks*blockIdx.y ] += bTheta[ 0 ];

}

template __global__ void findCoefficientsLCS2D_cudaKernel_0deg2 ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, size_t numOfStoredBlocks, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bf);
template __global__ void findCoefficientsLCS2D_cudaKernel_0deg2 ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, size_t numOfStoredBlocks, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bf);

/*
template<class T>
__global__
void findCoefficientsLCS2D_cudaKernel_0degConstantMemory2 ( T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfStoredBlocks, size_t dataShift, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.x;

	int pIndex = blockIdx.x*blockDim.x + threadIdx.x;
	{

		T dx = gpu_samplePoints_x[ pIndex +dataShift ] -stationaryPoints_x[blockIdx.y];
		T dy = gpu_samplePoints_y[ pIndex +dataShift ] -stationaryPoints_y[blockIdx.y];

		T theta = weightingFunction( sqrt(dx*dx + dy*dy) );

		bTheta[ threadIdx.x ] = theta*gpu_sampleValues_x [ pIndex +dataShift ];
		b_sh  [ threadIdx.x ] = theta*gpu_sampleValues_y [ pIndex +dataShift ];
		__syncthreads();


		size_t 	sind 	= threadIdx.x;
		int 	resid 	= blockDim.x;
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[0] += bTheta[sind]; b_sh[0] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 ) Sum_bf[blockIdx.x*2 +2*numOfStoredBlocks*blockIdx.y ] 		+= bTheta[ 0 ];
		if ( threadIdx.x == 0 )	Sum_bf[blockIdx.x*2 +2*numOfStoredBlocks*blockIdx.y +1 ] 	+= b_sh  [ 0 ]; 
		__syncthreads();
	

		bTheta[ threadIdx.x ] = theta;
		__syncthreads();
	}


	size_t 	sind 	= threadIdx.x;
	int 	resid 	= blockDim.x;
	size_t 	newReductionSize(1024);

	while ( resid > 0 ) {

		if (resid >=1024){ 	if ( threadIdx.x < 512) bTheta[ sind ] 	+= bTheta[ sind +512]; }
		else newReductionSize = 512; __syncthreads();
		if (resid >=512) {	if ( threadIdx.x < 256) bTheta[ sind ] 	+= bTheta[ sind +256]; }
		else newReductionSize = 256; __syncthreads();
		if (resid >=256) {	if ( threadIdx.x < 128) bTheta[ sind ] 	+= bTheta[ sind +128]; }
		else newReductionSize = 128; __syncthreads();
		if (resid >=128) {	if ( threadIdx.x <  64) bTheta[ sind ]	+= bTheta[ sind +64]; }
		else newReductionSize = 64; __syncthreads();
		if (resid >= 64) {	if ( threadIdx.x <  32) bTheta[ sind ] 	+= bTheta[ sind +32]; }
		else newReductionSize = 32; __syncthreads();
		if (resid >= 32) {	if ( threadIdx.x <  16) bTheta[ sind ] 	+= bTheta[ sind +16]; }
		else newReductionSize = 16; __syncthreads();
		if (resid >= 16) {	if ( threadIdx.x <   8) bTheta[ sind ] 	+= bTheta[ sind +8]; }
		else newReductionSize = 8; __syncthreads();
		if (resid >= 8 ) {	if ( threadIdx.x <   4) bTheta[ sind ] 	+= bTheta[ sind +4]; }
		else newReductionSize = 4; __syncthreads();
		if (resid >= 4 ) {	if ( threadIdx.x <   2) bTheta[ sind ] 	+= bTheta[ sind +2]; }
		else newReductionSize = 2; __syncthreads();
		if (resid >= 2 ) {	if ( threadIdx.x == 0 ) bTheta[ sind ] 	+= bTheta[ sind +1]; }
		else newReductionSize = 1; __syncthreads();

		if ( resid < blockDim.x && threadIdx.x == 0 ) bTheta[0] += bTheta[sind];
		__syncthreads();

		resid -= newReductionSize; sind += newReductionSize;
	}

	if ( threadIdx.x == 0 ) Sum_bbT[ blockIdx.x+numOfStoredBlocks*blockIdx.y ] += bTheta[ 0 ];

}

template __global__ void findCoefficientsLCS2D_cudaKernel_0degConstantMemory2 ( float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfStoredBlocks, size_t dataShift, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bf);
template __global__ void findCoefficientsLCS2D_cudaKernel_0degConstantMemory2 ( double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfStoredBlocks, size_t dataShift, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bf);
*/



template<class T>
__global__
void findCoefficientsLCS2D_cudaKernel ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.y*blockDim.x;
	T* bbT 		= b_sh +blockDim.y*blockDim.x;


	size_t mIndex = threadIdx.y*blockDim.x;
	int serialInBlockIndex = threadIdx.x +mIndex;
	size_t bfMonomIndex = threadIdx.y +blockDim.y*blockIdx.x;

	size_t bIndex = threadIdx.y*blockDim.y +threadIdx.x +blockDim.y*blockDim.y*blockIdx.x;
	size_t rowIndex = threadIdx.y*blockDim.y;
	size_t rowIndex2 = threadIdx.y*blockDim.x;

//	T stationaryPoint_x = stationaryPoints_x[ blockIdx.x ];
//	T stationaryPoint_y = stationaryPoints_y[ blockIdx.x ];

	T theta, b;

	size_t numOfTiles = (numOfSamplePoints +blockDim.x -1)/blockDim.x;
	for (size_t tile = 0; tile < numOfTiles; tile++ )
	{

		size_t pIndex = threadIdx.x +tile*blockDim.x;

		if ( pIndex < numOfSamplePoints ) {
			int power_x = powers_x[threadIdx.y];
			int power_y = powers_y[threadIdx.y];

			T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[ blockIdx.x ];
			T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[ blockIdx.x ];
			theta = weightingFunction( sqrt(dx*dx + dy*dy) );
			b = pow_i ( dx, power_x ) *pow_i( dy, power_y);

			bTheta[ serialInBlockIndex ] = b*theta*sampleValues_x [ pIndex ];
			b_sh  [ serialInBlockIndex ] = b*theta*sampleValues_y [ pIndex ];
		}

		__syncthreads();


		size_t 	sind 	= serialInBlockIndex;
		int 	resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		size_t 	newReductionSize(1024);
		
		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads(); 
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads(); 
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[mIndex] += bTheta[sind]; b_sh[mIndex] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 && tile == 0 ) { 	Sum_bfx[ bfMonomIndex ] = bTheta[ mIndex ];
					  		Sum_bfy[ bfMonomIndex ] = b_sh  [ mIndex ]; }


		if ( threadIdx.x == 0 && tile > 0 )  { 	Sum_bfx[ bfMonomIndex ] += bTheta[ mIndex ];
					  		Sum_bfy[ bfMonomIndex ] += b_sh  [ mIndex ]; }
		__syncthreads();


		bTheta [ serialInBlockIndex ] 	= b*theta; // weighted Vandermonde matrix
		//b_sh   [ serialInBlockIndex ] 	= b; // transpose of Vandermonde matrix
		b_sh   [ threadIdx.y +threadIdx.x*blockDim.y ] 	= b; // transpose of Vandermonde matrix

		__syncthreads();



		for ( size_t colsTile = 0; colsTile < (blockDim.y+blockDim.x-1)/blockDim.x; colsTile++ ) {

			size_t colIndex = threadIdx.x +colsTile*blockDim.x;

			if ( colIndex < blockDim.y ) {
				bbT[ rowIndex +colIndex ]  = b_sh[ colIndex ] *bTheta[ rowIndex2 ];
				for ( size_t pp = 1; pp < min( size_t(blockDim.x), size_t(numOfSamplePoints -tile*blockDim.x) ); pp++ ) // iterate over the sample points!!
					bbT[ rowIndex +colIndex ] += b_sh[ colIndex +pp*blockDim.y ] *bTheta[ rowIndex2 +pp ]; //  bbT[tid.y*blkdim.y +tid.x](the Gramian matrix is stored for the current sample point) = b_sh[tid.x] (all threads will get a monomial term, the term will be the same for threads with the same tid.x, but different along dim.x, for the current sample point) *bTheta[ tid.y*blkdim.x ] (all threads will get a monomial term, the term will be the same for threads with the same tid.y, but different along dim.y, for the current sample point)


				if ( tile ==0 ) Sum_bbT[ colsTile*blockDim.x +bIndex ]  = bbT[ rowIndex +colIndex ];
				if ( tile > 0 ) Sum_bbT[ colsTile*blockDim.x +bIndex ] += bbT[ rowIndex +colIndex ];
			}
			__syncthreads();
		}

/*
		for ( int col = 0; col < blockDim.y; col+=2 ) {

			size_t 	sind 	= serialInBlockIndex;
			int 	resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
			size_t 	newReductionSize(1024);

			bbT[ threadIdx.y*blockDim.y ] = bTheta[  ]*b;

			while ( resid > 0 ) {

				if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bbT[ sind ] 	+= bbT[ sind +512];
							if ( threadIdx.x > 511 && threadIdx.x < 1024) 	bbT[ sind -512]	+= bbT[ sind +512]; }
				else newReductionSize = 512; __syncthreads(); 
		
				if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
							if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
				else newReductionSize = 256; __syncthreads(); 
		
				if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
							if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
				else newReductionSize = 128; __syncthreads();
	
				if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
							if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
				else newReductionSize = 64; __syncthreads();
	
				if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
				 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
				else newReductionSize = 32; __syncthreads();
		
				if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
							if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
				else newReductionSize = 16; __syncthreads();
		
				if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
							if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
				else newReductionSize = 8; __syncthreads();
		
				if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
							if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
				else newReductionSize = 4; __syncthreads();
		
				if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
							if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
				else newReductionSize = 2; __syncthreads();
		
				if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
							if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
				else newReductionSize = 1; __syncthreads();
	
				if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[mIndex] += bTheta[sind]; b_sh[mIndex] += b_sh[sind]; }
				__syncthreads();
	
				resid -= newReductionSize; sind += newReductionSize;
			}
	
	
			if ( threadIdx.x == 0 && tile == 0 ) { 	Sum_bfx[ bfMonomIndex ] = bTheta[ mIndex ];
						  		Sum_bfy[ bfMonomIndex ] = b_sh  [ mIndex ]; }
	
	
			if ( threadIdx.x == 0 && tile > 0 )  { 	Sum_bfx[ bfMonomIndex ] += bTheta[ mIndex ];
						  		Sum_bfy[ bfMonomIndex ] += b_sh  [ mIndex ]; }
			__syncthreads();
		
	
			if ( threadIdx.x == 0 && tile == 0 ) Sum_bbT[ blockIdx.x ]  = bTheta[ 0 ];
			if ( threadIdx.x == 0 && tile  > 0 ) Sum_bbT[ blockIdx.x ] += bTheta[ 0 ];
	
			__syncthreads();
		}
*/
	}

}

template __global__ void findCoefficientsLCS2D_cudaKernel ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bfx, float* Sum_bfy);
template __global__ void findCoefficientsLCS2D_cudaKernel ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bfx, double* Sum_bfy);


template<class T>
__global__
void findCoefficientsLCS2D_cudaKernel2 ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, size_t numOfStoredBlocks, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.y*blockDim.x;
	T* bbT 		= b_sh +blockDim.y*blockDim.x;


	size_t mIndex = threadIdx.y*blockDim.x;
	int serialInBlockIndex = threadIdx.x +mIndex;
	int pIndex = blockIdx.x*blockDim.x + threadIdx.x;

	{

		T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[blockIdx.y];
		T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[blockIdx.y];
		int power_x = powers_x[threadIdx.y];
		int power_y = powers_y[threadIdx.y];

		T theta = weightingFunction( sqrt(dx*dx + dy*dy) );
		T b = pow_i ( dx, power_x ) *pow_i( dy, power_y);

		bTheta[ serialInBlockIndex ] = b*theta*sampleValues_x [ pIndex ];
		b_sh  [ serialInBlockIndex ] = b*theta*sampleValues_y [ pIndex ];
		__syncthreads();


		size_t 	sind 	= serialInBlockIndex;
		int 	resid 	= blockDim.x; // residual sample points
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads(); 
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads(); 
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[mIndex] += bTheta[sind]; b_sh[mIndex] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 ) { Sum_bf[threadIdx.y +blockIdx.x*2*blockDim.y +2*blockDim.y*numOfStoredBlocks*blockIdx.y ] 		+= bTheta[ mIndex ];
					  Sum_bf[threadIdx.y +blockIdx.x*2*blockDim.y +2*blockDim.y*numOfStoredBlocks*blockIdx.y +blockDim.y ] 	+= b_sh  [ mIndex ]; }
		__syncthreads();

		
		bTheta [ serialInBlockIndex ] = b*theta;
		b_sh   [ threadIdx.y +threadIdx.x*blockDim.y ] = b;
	
		__syncthreads();
	}



	size_t bIndex = threadIdx.y*blockDim.y +threadIdx.x +blockIdx.x*blockDim.y*blockDim.y +blockDim.y*blockDim.y*numOfStoredBlocks*blockIdx.y;
	size_t rowIndex = threadIdx.y*blockDim.y;
	size_t rowIndex2 = threadIdx.y*blockDim.x;

	for ( size_t colsTile = 0; colsTile < (blockDim.y+blockDim.x-1)/blockDim.x; colsTile++ ) {
		size_t colIndex = threadIdx.x +colsTile*blockDim.x;

		if ( colIndex < blockDim.y ) {
			bbT[ rowIndex +colIndex ]  = b_sh[ colIndex ] *bTheta[ rowIndex2 ];
			for ( size_t pp = 1; pp < blockDim.x; pp++ )
				bbT[ rowIndex +colIndex ] += b_sh[ colIndex +pp*blockDim.y ] *bTheta[ rowIndex2 +pp ];

			Sum_bbT[ colsTile*blockDim.x +bIndex ] += bbT[ rowIndex +colIndex ]; }
		__syncthreads();
	}

}

template __global__ void findCoefficientsLCS2D_cudaKernel2 ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, size_t numOfStoredBlocks, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bf);
template __global__ void findCoefficientsLCS2D_cudaKernel2 ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, size_t numOfStoredBlocks, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bf);

/*
template<class T>
__global__
void findCoefficientsLCS2D_cudaKernelConstantMemory2 (T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfStoredBlocks, size_t dataShift, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.y*blockDim.x;
	T* bbT 		= b_sh +blockDim.y*blockDim.x;


	size_t mIndex = threadIdx.y*blockDim.x;
	int serialInBlockIndex = threadIdx.x +mIndex;
	int pIndex = blockIdx.x*blockDim.x + threadIdx.x;

	{

		T dx = gpu_samplePoints_x[ pIndex+dataShift ] -stationaryPoints_x[blockIdx.y];
		T dy = gpu_samplePoints_y[ pIndex+dataShift ] -stationaryPoints_y[blockIdx.y];
		int power_x = gpu_powers_x[threadIdx.y];
		int power_y = gpu_powers_y[threadIdx.y];

		T theta = weightingFunction( sqrt(dx*dx + dy*dy) );
		T b = pow_i ( dx, power_x ) *pow_i( dy, power_y);

		bTheta[ serialInBlockIndex ] = b*theta*gpu_sampleValues_x [ pIndex+dataShift ];
		b_sh  [ serialInBlockIndex ] = b*theta*gpu_sampleValues_y [ pIndex+dataShift ];
		__syncthreads();


		size_t 	sind 	= serialInBlockIndex;
		int 	resid 	= blockDim.x; // residual sample points
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads(); 
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads(); 
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[mIndex] += bTheta[sind]; b_sh[mIndex] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 ) { Sum_bf[threadIdx.y +blockIdx.x*2*blockDim.y +2*blockDim.y*numOfStoredBlocks*blockIdx.y ] 		+= bTheta[ mIndex ];
					  Sum_bf[threadIdx.y +blockIdx.x*2*blockDim.y +2*blockDim.y*numOfStoredBlocks*blockIdx.y +blockDim.y ] 	+= b_sh  [ mIndex ]; }
		__syncthreads();

		
		bTheta [ serialInBlockIndex ] = b*theta;
		b_sh   [ threadIdx.y +threadIdx.x*blockDim.y ] = b;
	
		__syncthreads();
	}



	size_t bIndex = threadIdx.y*blockDim.y +threadIdx.x +blockIdx.x*blockDim.y*blockDim.y +blockDim.y*blockDim.y*numOfStoredBlocks*blockIdx.y;
	size_t rowIndex = threadIdx.y*blockDim.y;
	size_t rowIndex2 = threadIdx.y*blockDim.x;

	for ( size_t colsTile = 0; colsTile < (blockDim.y+blockDim.x-1)/blockDim.x; colsTile++ ) {
		size_t colIndex = threadIdx.x +colsTile*blockDim.x;

		if ( colIndex < blockDim.y ) {
			bbT[ rowIndex +colIndex ]  = b_sh[ colIndex ] *bTheta[ rowIndex2 ];
			for ( size_t pp = 1; pp < blockDim.x; pp++ )
				bbT[ rowIndex +colIndex ] += b_sh[ colIndex +pp*blockDim.y ] *bTheta[ rowIndex2 +pp ];

			Sum_bbT[ colsTile*blockDim.x +bIndex ] += bbT[ rowIndex +colIndex ]; }
		__syncthreads();
	}

}

template __global__ void findCoefficientsLCS2D_cudaKernelConstantMemory2 ( float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfStoredBlocks, size_t dataShift, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bf);
template __global__ void findCoefficientsLCS2D_cudaKernelConstantMemory2 ( double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfStoredBlocks, size_t dataShift, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bf);
*/



template<class T>
__global__
void findCoefficientsLCS2D_cudaKernel3 ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.y*blockDim.x;
	T* bbT 		= b_sh +blockDim.y*blockDim.x;


	size_t mIndex = threadIdx.y*blockDim.x;
	int serialInBlockIndex = threadIdx.x +mIndex;
	size_t bfMonomIndex = threadIdx.y +blockDim.y*blockIdx.x;


	T theta, b;

	size_t numOfTiles = (numOfSamplePoints +blockDim.x -1)/blockDim.x;
	for (size_t tile = 0; tile < numOfTiles; tile++ )
	{

		size_t pIndex = threadIdx.x +tile*blockDim.x;

		if ( pIndex < numOfSamplePoints ) {
			int power_x = powers_x[threadIdx.y];
			int power_y = powers_y[threadIdx.y];

			T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[ blockIdx.x ];
			T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[ blockIdx.x ];
			theta = weightingFunction( sqrt(dx*dx + dy*dy) );
			b = pow_i ( dx, power_x ) *pow_i( dy, power_y);

			bTheta[ serialInBlockIndex ] = b*theta*sampleValues_x [ pIndex ];
			b_sh  [ serialInBlockIndex ] = b*theta*sampleValues_y [ pIndex ];
		}

		__syncthreads();


		size_t 	sind 	= serialInBlockIndex;
		int 	resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		size_t 	newReductionSize(1024);
		
		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind]; }
			else newReductionSize = 512; __syncthreads(); 
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind]; }
			else newReductionSize = 256; __syncthreads(); 
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; }
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; }
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; }
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[mIndex] += bTheta[sind]; b_sh[mIndex] += b_sh[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 && tile == 0 ) { 	Sum_bfx[ bfMonomIndex ] = bTheta[ mIndex ];
					  		Sum_bfy[ bfMonomIndex ] = b_sh  [ mIndex ]; }


		if ( threadIdx.x == 0 && tile > 0 )  { 	Sum_bfx[ bfMonomIndex ] += bTheta[ mIndex ];
					  		Sum_bfy[ bfMonomIndex ] += b_sh  [ mIndex ]; }
		__syncthreads();
		

		bTheta [ serialInBlockIndex ] 			= b*theta;
		b_sh   [ threadIdx.y +threadIdx.x*blockDim.y ] 	= b;
	
		__syncthreads();



		for ( size_t colsTile = 0; colsTile < (blockDim.y+blockDim.x-1)/blockDim.x; colsTile++ ) {

			size_t colIndex = threadIdx.x +colsTile*blockDim.x;
			size_t numOfPoints = min( size_t(blockDim.x), size_t(numOfSamplePoints -tile*blockDim.x) );
			size_t bIndex = threadIdx.y*blockDim.y +threadIdx.x +blockDim.y*blockDim.y*blockIdx.x;
			size_t rowIndex2 = threadIdx.y*blockDim.x;
			size_t bbTIndex = threadIdx.y*blockDim.y +colIndex;
			size_t bbTSize = blockDim.y*blockDim.y;
			
			if ( colIndex < blockDim.y ) {

							bbT[ bbTIndex 		 ] = b_sh[ colIndex 		  ] *bTheta[ rowIndex2 ];
				if (numOfPoints > 1) 	bbT[ bbTIndex +  bbTSize ] = b_sh[ colIndex +  blockDim.y ] *bTheta[ rowIndex2 +1 ];
				else			bbT[ bbTIndex +  bbTSize ] = 0;
				if (numOfPoints > 2) 	bbT[ bbTIndex +2*bbTSize ] = b_sh[ colIndex +2*blockDim.y ] *bTheta[ rowIndex2 +2 ];
				else			bbT[ bbTIndex +2*bbTSize ] = 0;
				if (numOfPoints > 3) 	bbT[ bbTIndex +3*bbTSize ] = b_sh[ colIndex +3*blockDim.y ] *bTheta[ rowIndex2 +3 ];
				else			bbT[ bbTIndex +3*bbTSize ] = 0;

				for ( size_t pp = 4; pp < numOfPoints; pp +=4 ) {
									bbT[ bbTIndex		 ] += b_sh[ colIndex +	  pp*blockDim.y ] *bTheta[ rowIndex2   +pp ];
					if (numOfPoints >   pp) 	bbT[ bbTIndex +  bbTSize ] += b_sh[ colIndex +(1+pp)*blockDim.y ] *bTheta[ rowIndex2 +1+pp ];
					if (numOfPoints > 2+pp) 	bbT[ bbTIndex +2*bbTSize ] += b_sh[ colIndex +(2+pp)*blockDim.y ] *bTheta[ rowIndex2 +2+pp ];
					if (numOfPoints > 3+pp) 	bbT[ bbTIndex +3*bbTSize ] += b_sh[ colIndex +(3+pp)*blockDim.y ] *bTheta[ rowIndex2 +3+pp ];
				}

				bbT[ bbTIndex		] += bbT[ bbTIndex +2*bbTSize];
				bbT[ bbTIndex +bbTSize  ] += bbT[ bbTIndex +3*bbTSize];
				bbT[ bbTIndex		] += bbT[ bbTIndex +bbTSize];

				if ( tile ==0 ) Sum_bbT[ colsTile*blockDim.x +bIndex ]  = bbT[ bbTIndex ];
				if ( tile > 0 ) Sum_bbT[ colsTile*blockDim.x +bIndex ] += bbT[ bbTIndex ];
			}
			__syncthreads();
		}

	}
}

template __global__ void findCoefficientsLCS2D_cudaKernel3 ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bfx, float* Sum_bfy);
template __global__ void findCoefficientsLCS2D_cudaKernel3 ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bfx, double* Sum_bfy);



template<class T>
__global__
void find_bTheta_LCS2D_cudaKernel ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, WeightingFunction<T> weightingFunction, T* bTheta, T* b_out)
{
	size_t serialInBlockIndex 		= threadIdx.x +threadIdx.y*blockDim.x;
	size_t serialInGridIndex  		= serialInBlockIndex +blockIdx.x*blockDim.x*blockDim.y;
	size_t totalGridSize	  		= blockDim.x*blockDim.y*gridDim.x;
	size_t shiftOfDataPerMonomial		= numOfSamplePoints*threadIdx.y;
	size_t shiftOfDataPerStationaryPoint 	= blockIdx.y*numOfSamplePoints*blockDim.y;
	size_t xGridIndex			= blockIdx.x*blockDim.x + threadIdx.x;
	size_t xGridSize			= blockDim.x*gridDim.x;


	int power_x = powers_x[threadIdx.y];
	int power_y = powers_y[threadIdx.y];
	T stationaryPoint_x = stationaryPoints_x[blockIdx.y];
       	T stationaryPoint_y = stationaryPoints_y[blockIdx.y];
	T dx, dy, theta, b;
	size_t numOfTiles = (numOfSamplePoints +xGridSize -1)/(xGridSize);


	for ( size_t tile = 0; tile < numOfTiles; tile++ ) {

		size_t point = xGridIndex +tile*xGridSize;
		if ( point < numOfSamplePoints ) {
			dx = samplePoints_x[ point ] -stationaryPoint_x;
			dy = samplePoints_y[ point ] -stationaryPoint_y;
			theta = weightingFunction( sqrt(dx*dx + dy*dy) );
			b = pow_i ( dx, power_x ) *pow_i( dy, power_y);

			bTheta[ point +shiftOfDataPerMonomial +shiftOfDataPerStationaryPoint ] = b*theta;
			b_out [ point +shiftOfDataPerMonomial +shiftOfDataPerStationaryPoint ] = b; }
		__syncthreads(); 
	}
}

template __global__ void find_bTheta_LCS2D_cudaKernel ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, WeightingFunction<float> weightingFunction, float* bTheta, float* b_out);
template __global__ void find_bTheta_LCS2D_cudaKernel ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, WeightingFunction<double> weightingFunction, double* bTheta, double* b_out);


template<class T>
__global__
void find_bTheta_bf_LCS2D_cudaKernel ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfMonomials, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, WeightingFunction<T> weightingFunction, T* bTheta, T* b_out, T* bfx_out, T* bfy_out)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *bfx = reinterpret_cast<T*>(sharedMemory);
	T *bfy = bfx +blockDim.x;


	size_t bf_shift	= blockIdx.z*numOfMonomials;
	size_t b_shift	= blockIdx.z*numOfMonomials*numOfSamplePoints;
	size_t mono_shift = blockIdx.y*numOfSamplePoints;
	size_t xGridIndex = blockIdx.x*blockDim.x + threadIdx.x;
	size_t xGridSize  = blockDim.x*gridDim.x;


	int power_x = powers_x[blockIdx.y];
	int power_y = powers_y[blockIdx.y];
	T stationaryPoint_x = stationaryPoints_x[blockIdx.z];
	T stationaryPoint_y = stationaryPoints_y[blockIdx.z];
	T dx, dy, theta, b;


	if ( xGridIndex < numOfSamplePoints ) {
		dx = samplePoints_x[ xGridIndex ] -stationaryPoint_x;
		dy = samplePoints_y[ xGridIndex ] -stationaryPoint_y;
		theta = weightingFunction( sqrt(dx*dx + dy*dy) );
		b = pow_i ( dx, power_x ) *pow_i( dy, power_y);

		bfx [ threadIdx.x ] = b*theta*sampleValues_x [ xGridIndex ];
		bfy [ threadIdx.x ] = b*theta*sampleValues_y [ xGridIndex ]; }
	__syncthreads();


	size_t 	sind 	= threadIdx.x;
	int 	resid 	= min(int(blockDim.x), int(numOfSamplePoints -blockDim.x*blockIdx.x) ); // num of residual sample points
	size_t 	newReductionSize(1024);

	while ( resid > 0 ) {

		if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bfx[ sind ] 	+= bfx[ sind +512];
					if ( threadIdx.x > 511 && threadIdx.x < 1024) 	bfy[ sind -512] += bfy[ sind]; }
		else newReductionSize = 512; __syncthreads(); 

		if (resid >=512) {	if ( threadIdx.x < 256  ) 			bfx[ sind ] 	+= bfx[ sind +256];
					if ( threadIdx.x > 255 && threadIdx.x < 512) 	bfy[ sind -256] += bfy[ sind]; }
		else newReductionSize = 256; __syncthreads(); 

		if (resid >=256) {	if ( threadIdx.x < 128) 			bfx[ sind ] 	+= bfx[ sind +128];
					if ( threadIdx.x > 127 && threadIdx.x < 256) 	bfy[ sind -128] += bfy[ sind]; }
		else newReductionSize = 128; __syncthreads();

		if (resid >=128) {	if ( threadIdx.x <  64) 			bfx[ sind ]	+= bfx[ sind +64];
					if ( threadIdx.x >  63 && threadIdx.x < 128) 	bfy[ sind - 64] += bfy[ sind]; }
		else newReductionSize = 64; __syncthreads();

		if (resid >= 64) {	if ( threadIdx.x <  32) 			bfx[ sind ] 	+= bfx[ sind +32];
		 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	bfy[ sind - 32] += bfy[ sind]; }
		else newReductionSize = 32; __syncthreads();

		if (resid >= 32) {	if ( threadIdx.x <  16) 			bfx[ sind ] 	+= bfx[ sind +16];
					if ( threadIdx.x >  15 && threadIdx.x <  32) 	bfy[ sind - 16] += bfy[ sind]; }
		else newReductionSize = 16; __syncthreads();

		if (resid >= 16) {	if ( threadIdx.x <   8) 			bfx[ sind ] 	+= bfx[ sind +8];
					if ( threadIdx.x >   7 && threadIdx.x <  16) 	bfy[ sind -  8] += bfy[ sind]; }
		else newReductionSize = 8; __syncthreads();

		if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bfx[ sind ] 	+= bfx[ sind +4];
					if ( threadIdx.x >   3 && threadIdx.x <   8) 	bfy[ sind -  4] += bfy[ sind]; }
		else newReductionSize = 4; __syncthreads();

		if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bfx[ sind ] 	+= bfx[ sind +2];
					if ( threadIdx.x >   1 && threadIdx.x <   4) 	bfy[ sind -  2] += bfy[ sind]; }
		else newReductionSize = 2; __syncthreads();

		if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bfx[ sind ] 	+= bfx[ sind +1];
					if ( threadIdx.x == 1 )				bfy[ sind -  1] += bfy[ sind]; }
		else newReductionSize = 1; __syncthreads();

		if ( resid < blockDim.x && threadIdx.x == 0 )  { bfx[0] += bfx[sind]; bfy[0] += bfy[sind]; }
		__syncthreads();

		resid -= newReductionSize; sind += newReductionSize;
	}


	if ( threadIdx.x == 0 ) {
		bfx_out [blockIdx.y +bf_shift] = bfx[0];
		bfy_out [blockIdx.y +bf_shift] = bfy[0]; }
	__syncthreads();


	if ( threadIdx.x == 0 && blockIdx.x > 0 ) {
		atomicAdd( bfx_out +blockIdx.y +bf_shift, bfx[0] );
		atomicAdd( bfy_out +blockIdx.y +bf_shift, bfy[0] ); }
	__syncthreads();


	bTheta [ xGridIndex +mono_shift +b_shift ] = b*theta;
	b_out  [ xGridIndex +mono_shift +b_shift ] = b;
	__syncthreads();

}

template __global__ void find_bTheta_bf_LCS2D_cudaKernel ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfMonomials, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, WeightingFunction<float> weightingFunction, float* bTheta, float* b_out, float* bfx_out, float* bfy_out);
template __global__ void find_bTheta_bf_LCS2D_cudaKernel ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfMonomials, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, WeightingFunction<double> weightingFunction, double* bTheta, double* b_out, double* bfx_out, double* bfy_out);


template<class T>
__global__
void find_bTheta_bf_LCS2D_cudaKernel2 ( int* powers_x, int* powers_y, T* stationaryPoints_x, T* stationaryPoints_y, size_t numOfMonomials, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* sampleValues_x, T* sampleValues_y, WeightingFunction<T> weightingFunction, T* bTheta, T* b_out, T* bfx_out, T* bfy_out)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *bfx = reinterpret_cast<T*>(sharedMemory);
	T *bfy = bfx +blockDim.x;


	size_t bf_shift		= blockIdx.y*numOfMonomials;
	size_t b_shift		= blockIdx.y*numOfMonomials*numOfSamplePoints;
	size_t mono_shift 	= blockIdx.x*numOfSamplePoints;

	size_t numOfTiles = (numOfSamplePoints +blockDim.x -1)/blockDim.x;
	for ( size_t tile = 0; tile < numOfTiles; tile++ ) {


		if ( blockDim.x*tile +threadIdx.x < numOfSamplePoints ) {

			int power_x 		= powers_x[blockIdx.x];
			int power_y 		= powers_y[blockIdx.x];
	
			T dx 	= samplePoints_x[ threadIdx.x +tile*blockDim.x ] -stationaryPoints_x[blockIdx.y];
			T dy 	= samplePoints_y[ threadIdx.x +tile*blockDim.x ] -stationaryPoints_y[blockIdx.y];
			T theta = weightingFunction( sqrt(dx*dx + dy*dy) );
			T b 	= pow_i ( dx, power_x ) *pow_i( dy, power_y);


			bfx [ threadIdx.x ] = b*theta*sampleValues_x [ threadIdx.x +tile*blockDim.x ];
			bfy [ threadIdx.x ] = b*theta*sampleValues_y [ threadIdx.x +tile*blockDim.x ];


			bTheta [ threadIdx.x +tile*blockDim.x +mono_shift +b_shift ] = b*theta;
			b_out  [ threadIdx.x +tile*blockDim.x +mono_shift +b_shift ] = b;
		}
		__syncthreads();
	
	
		size_t 	sind 	= threadIdx.x;
		int 	resid 	= min( int(blockDim.x), int( numOfSamplePoints -blockDim.x*tile ) ); // num of residual sample points
		size_t 	newReductionSize(1024);
	
		while ( resid > 0 ) {
	
			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bfx[ sind ] 	+= bfx[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	bfy[ sind -512] += bfy[ sind]; }
			else newReductionSize = 512; __syncthreads(); 
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bfx[ sind ] 	+= bfx[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	bfy[ sind -256] += bfy[ sind]; }
			else newReductionSize = 256; __syncthreads(); 
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bfx[ sind ] 	+= bfx[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	bfy[ sind -128] += bfy[ sind]; }
			else newReductionSize = 128; __syncthreads();
	
			if (resid >=128) {	if ( threadIdx.x <  64) 			bfx[ sind ]	+= bfx[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	bfy[ sind - 64] += bfy[ sind]; }
			else newReductionSize = 64; __syncthreads();
	
			if (resid >= 64) {	if ( threadIdx.x <  32) 			bfx[ sind ] 	+= bfx[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	bfy[ sind - 32] += bfy[ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bfx[ sind ] 	+= bfx[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	bfy[ sind - 16] += bfy[ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bfx[ sind ] 	+= bfx[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	bfy[ sind -  8] += bfy[ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bfx[ sind ] 	+= bfx[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	bfy[ sind -  4] += bfy[ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bfx[ sind ] 	+= bfx[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	bfy[ sind -  2] += bfy[ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bfx[ sind ] 	+= bfx[ sind +1];
						if ( threadIdx.x == 1 )				bfy[ sind -  1] += bfy[ sind]; }
			else newReductionSize = 1; __syncthreads();
	
			if ( resid < blockDim.x && threadIdx.x == 0 )  { bfx[0] += bfx[sind]; bfy[0] += bfy[sind]; }
			__syncthreads();
	
			resid -= newReductionSize; sind += newReductionSize;
		}
	
	
		if ( threadIdx.x == 0 && tile == 0 ) {
			bfx_out [blockIdx.x +bf_shift] = bfx[0];
			bfy_out [blockIdx.x +bf_shift] = bfy[0]; }

		if ( threadIdx.x == 0 && tile > 0 ) {
			bfx_out [blockIdx.x +bf_shift] += bfx[0];
			bfy_out [blockIdx.x +bf_shift] += bfy[0]; }
		__syncthreads();

	}

}

template __global__ void find_bTheta_bf_LCS2D_cudaKernel2 ( int* powers_x, int* powers_y, float* stationaryPoints_x, float* stationaryPoints_y, size_t numOfMonomials, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* sampleValues_x, float* sampleValues_y, WeightingFunction<float> weightingFunction, float* bTheta, float* b_out, float* bfx_out, float* bfy_out);
template __global__ void find_bTheta_bf_LCS2D_cudaKernel2 ( int* powers_x, int* powers_y, double* stationaryPoints_x, double* stationaryPoints_y, size_t numOfMonomials, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* sampleValues_x, double* sampleValues_y, WeightingFunction<double> weightingFunction, double* bTheta, double* b_out, double* bfx_out, double* bfy_out);



template<class T, int tileWidth>
__global__
void find_bbT_LCS2D_cudaKernel ( size_t numOfSamplePoints, size_t numOfMonomials, T* bTheta, T* b, T* bbT )
{
	size_t bbT_shift	= blockIdx.z*numOfMonomials*numOfMonomials;
	size_t b_shift		= blockIdx.z*numOfMonomials*numOfSamplePoints;

	size_t numOfInnerTiles 	= (numOfSamplePoints +tileWidth -1) /tileWidth;

	__shared__ float Mds[tileWidth][tileWidth], Nds[tileWidth][tileWidth];

	T Pvalue(0);
	size_t M_Row = blockIdx.y*blockDim.y +threadIdx.y, N_Col = blockIdx.x*blockDim.y +threadIdx.y;

	for (size_t tileIndex = 0; tileIndex < numOfInnerTiles; tileIndex++)
	{
		if ( tileIndex*tileWidth +threadIdx.x < numOfSamplePoints && M_Row < numOfMonomials && N_Col < numOfMonomials ) {
			Mds [threadIdx.y][threadIdx.x] = bTheta [ tileIndex*tileWidth +threadIdx.x +M_Row*numOfSamplePoints +b_shift ];
			Nds [threadIdx.y][threadIdx.x] = b      [ tileIndex*tileWidth +threadIdx.x +N_Col*numOfSamplePoints +b_shift ]; }
		__syncthreads();

		if ( M_Row < numOfMonomials && blockIdx.x*blockDim.x +threadIdx.x < numOfMonomials ) {
			for (size_t k = 0; k < tileWidth; k++) Pvalue += Mds[threadIdx.y][k] *Nds[threadIdx.x][k];
			bbT [ M_Row*numOfMonomials +blockIdx.x*blockDim.x +threadIdx.x +bbT_shift ] = Pvalue; }
		 __syncthreads();
	}

}

template __global__ void find_bbT_LCS2D_cudaKernel ( size_t numOfSamplePoints, size_t numOfMonomials, float* bTheta, float* b, float* bbT );
template __global__ void find_bbT_LCS2D_cudaKernel ( size_t numOfSamplePoints, size_t numOfMonomials, double* bTheta, double* b, double* bbT );



template<class T>
__global__
void find_bf_LCS2D_cudaKernel ( size_t numOfSamplePoints, T* sampleValues_x, T* sampleValues_y, size_t numOfMonomials, T* b, T* bfx_out, T* bfy_out )
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *bfx = reinterpret_cast<T*>(sharedMemory);
	T *bfy = bfx +blockDim.x;


	size_t bf_shift	= blockIdx.z*numOfMonomials;
	size_t b_shift	= blockIdx.z*numOfMonomials*numOfSamplePoints;
	size_t mono_shift = blockIdx.y*numOfSamplePoints;
	size_t xGridIndex = blockIdx.x*blockDim.x + threadIdx.x;
	size_t xGridSize  = blockDim.x*gridDim.x;


	if ( xGridIndex < numOfSamplePoints ) {
		bfx [ threadIdx.x ] = b [ xGridIndex +mono_shift +b_shift ]*sampleValues_x [ xGridIndex ];
		bfy [ threadIdx.x ] = b [ xGridIndex +mono_shift +b_shift ]*sampleValues_y [ xGridIndex ]; }
	__syncthreads();


	size_t 	sind 	= threadIdx.x;
	int 	resid 	= min(int(blockDim.x), int(numOfSamplePoints -blockDim.x*blockIdx.x) ); // num of residual sample points
	size_t 	newReductionSize(1024);


	while ( resid > 0 ) {

		if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bfx[ sind ] 	+= bfx[ sind +512];
					if ( threadIdx.x > 511 && threadIdx.x < 1024) 	bfy[ sind -512] += bfy[ sind]; }
		else newReductionSize = 512; __syncthreads(); 

		if (resid >=512) {	if ( threadIdx.x < 256  ) 			bfx[ sind ] 	+= bfx[ sind +256];
					if ( threadIdx.x > 255 && threadIdx.x < 512) 	bfy[ sind -256] += bfy[ sind]; }
		else newReductionSize = 256; __syncthreads(); 

		if (resid >=256) {	if ( threadIdx.x < 128) 			bfx[ sind ] 	+= bfx[ sind +128];
					if ( threadIdx.x > 127 && threadIdx.x < 256) 	bfy[ sind -128] += bfy[ sind]; }
		else newReductionSize = 128; __syncthreads();

		if (resid >=128) {	if ( threadIdx.x <  64) 			bfx[ sind ]	+= bfx[ sind +64];
					if ( threadIdx.x >  63 && threadIdx.x < 128) 	bfy[ sind - 64] += bfy[ sind]; }
		else newReductionSize = 64; __syncthreads();

		if (resid >= 64) {	if ( threadIdx.x <  32) 			bfx[ sind ] 	+= bfx[ sind +32];
		 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	bfy[ sind - 32] += bfy[ sind]; }
		else newReductionSize = 32; __syncthreads();

		if (resid >= 32) {	if ( threadIdx.x <  16) 			bfx[ sind ] 	+= bfx[ sind +16];
					if ( threadIdx.x >  15 && threadIdx.x <  32) 	bfy[ sind - 16] += bfy[ sind]; }
		else newReductionSize = 16; __syncthreads();

		if (resid >= 16) {	if ( threadIdx.x <   8) 			bfx[ sind ] 	+= bfx[ sind +8];
					if ( threadIdx.x >   7 && threadIdx.x <  16) 	bfy[ sind -  8] += bfy[ sind]; }
		else newReductionSize = 8; __syncthreads();

		if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bfx[ sind ] 	+= bfx[ sind +4];
					if ( threadIdx.x >   3 && threadIdx.x <   8) 	bfy[ sind -  4] += bfy[ sind]; }
		else newReductionSize = 4; __syncthreads();

		if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bfx[ sind ] 	+= bfx[ sind +2];
					if ( threadIdx.x >   1 && threadIdx.x <   4) 	bfy[ sind -  2] += bfy[ sind]; }
		else newReductionSize = 2; __syncthreads();

		if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bfx[ sind ] 	+= bfx[ sind +1];
					if ( threadIdx.x == 1 )				bfy[ sind -  1] += bfy[ sind]; }
		else newReductionSize = 1; __syncthreads();

		if ( resid < blockDim.x && threadIdx.x == 0 )  { bfx[0] += bfx[sind]; bfy[0] += bfy[sind]; }
		__syncthreads();

		resid -= newReductionSize; sind += newReductionSize;
	}


	if ( threadIdx.x == 0 ) {
		bfx_out [blockIdx.y +bf_shift] = bfx[0];
		bfy_out [blockIdx.y +bf_shift] = bfy[0]; }
	__syncthreads();


	if ( threadIdx.x == 0 && blockIdx.x > 0 ) {
		atomicAdd( bfx_out +blockIdx.y +bf_shift, bfx[0] );
		atomicAdd( bfy_out +blockIdx.y +bf_shift, bfy[0] ); }
	__syncthreads();

}
template __global__ void find_bf_LCS2D_cudaKernel ( size_t numOfSamplePoints, float* sampleValues_x, float* sampleValues_y, size_t numOfMonomials, float* b, float* bfx_out, float* bfy_out );
template __global__ void find_bf_LCS2D_cudaKernel ( size_t numOfSamplePoints, double* sampleValues_x, double* sampleValues_y, size_t numOfMonomials, double* b, double* bfx_out, double* bfy_out );


template<class T>
__global__
void find_bf_LCS2D_cudaKernel2 ( size_t numOfSamplePoints, T* sampleValues_x, T* sampleValues_y, size_t numOfMonomials, T* b, T* bfx_out, T* bfy_out )
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T *bfx = reinterpret_cast<T*>(sharedMemory);
	T *bfy = bfx +blockDim.x;


	size_t bf_shift	= blockIdx.y*numOfMonomials;
	size_t b_shift	= blockIdx.y*numOfMonomials*numOfSamplePoints;
	size_t mono_shift = blockIdx.x*numOfSamplePoints;


	size_t numOfTiles = ( numOfSamplePoints +blockDim.x -1)/blockDim.x;
	for ( size_t tile = 0; tile < numOfTiles; tile++ ) {

		if ( blockDim.x*tile +threadIdx.x < numOfSamplePoints  ) {
			bfx [ threadIdx.x ] = b [ threadIdx.x +tile*blockDim.x +mono_shift +b_shift ]*sampleValues_x [ threadIdx.x +tile*blockDim.x ];
			bfy [ threadIdx.x ] = b [ threadIdx.x +tile*blockDim.x +mono_shift +b_shift ]*sampleValues_y [ threadIdx.x +tile*blockDim.x ]; }
		__syncthreads();


		size_t 	sind 	= threadIdx.x;
		int 	resid 	= min(int(blockDim.x), int(numOfSamplePoints -blockDim.x*tile) ); // num of residual sample points
		size_t 	newReductionSize(1024);
	
	
		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bfx[ sind ] 	+= bfx[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	bfy[ sind -512] += bfy[ sind]; }
			else newReductionSize = 512; __syncthreads(); 
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bfx[ sind ] 	+= bfx[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	bfy[ sind -256] += bfy[ sind]; }
			else newReductionSize = 256; __syncthreads(); 
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bfx[ sind ] 	+= bfx[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	bfy[ sind -128] += bfy[ sind]; }
			else newReductionSize = 128; __syncthreads();
	
			if (resid >=128) {	if ( threadIdx.x <  64) 			bfx[ sind ]	+= bfx[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	bfy[ sind - 64] += bfy[ sind]; }
			else newReductionSize = 64; __syncthreads();
	
			if (resid >= 64) {	if ( threadIdx.x <  32) 			bfx[ sind ] 	+= bfx[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	bfy[ sind - 32] += bfy[ sind]; }
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bfx[ sind ] 	+= bfx[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	bfy[ sind - 16] += bfy[ sind]; }
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bfx[ sind ] 	+= bfx[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	bfy[ sind -  8] += bfy[ sind]; }
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bfx[ sind ] 	+= bfx[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	bfy[ sind -  4] += bfy[ sind]; }
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bfx[ sind ] 	+= bfx[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	bfy[ sind -  2] += bfy[ sind]; }
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bfx[ sind ] 	+= bfx[ sind +1];
						if ( threadIdx.x == 1 )				bfy[ sind -  1] += bfy[ sind]; }
			else newReductionSize = 1; __syncthreads();
	
			if ( resid < blockDim.x && threadIdx.x == 0 )  { bfx[0] += bfx[sind]; bfy[0] += bfy[sind]; }
			__syncthreads();
	
			resid -= newReductionSize; sind += newReductionSize;
		}

		if ( threadIdx.x == 0 && tile == 0 ) {
			bfx_out [blockIdx.x +bf_shift] = bfx[0];
			bfy_out [blockIdx.x +bf_shift] = bfy[0]; }

		if ( threadIdx.x == 0 && tile > 0 ) {
			bfx_out [blockIdx.x +bf_shift] += bfx[0];
			bfy_out [blockIdx.x +bf_shift] += bfy[0]; }
		__syncthreads();

	}

}

template __global__ void find_bf_LCS2D_cudaKernel2 ( size_t numOfSamplePoints, float* sampleValues_x, float* sampleValues_y, size_t numOfMonomials, float* b, float* bfx_out, float* bfy_out );
template __global__ void find_bf_LCS2D_cudaKernel2 ( size_t numOfSamplePoints, double* sampleValues_x, double* sampleValues_y, size_t numOfMonomials, double* b, double* bfx_out, double* bfy_out );
