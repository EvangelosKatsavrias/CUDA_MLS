#include"findCoefficients3D_cuda_kernels.h"
#include<stdio.h>
#include"atomicOperations_doubles.h"
#include"myBLAS.h"


template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel_0deg ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy, T* Sum_bfz)
{

	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.x;
	T* b_sh2	= b_sh +blockDim.x;

	T theta;

	size_t numOfTiles = (numOfSamplePoints +blockDim.x -1)/blockDim.x;
	for (size_t tile = 0; tile < numOfTiles; tile++ ) {

		size_t pIndex = threadIdx.x +tile*blockDim.x;
		if ( pIndex < numOfSamplePoints ) {
			T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[ blockIdx.x ];
			T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[ blockIdx.x ];
			T dz = samplePoints_z[ pIndex ] -stationaryPoints_z[ blockIdx.x ];
			theta = weightingFunction( sqrt(dx*dx + dy*dy + dz*dz) );

			bTheta[ threadIdx.x ] = theta*sampleValues_x [ pIndex ];
			b_sh  [ threadIdx.x ] = theta*sampleValues_y [ pIndex ];
			b_sh2 [ threadIdx.x ] = theta*sampleValues_z [ pIndex ];
		}
		__syncthreads();


		size_t 	sind 	= threadIdx.x;
		int 	resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind];
						if ( threadIdx.x < 512	) 			b_sh2 [ sind ]	+= b_sh2 [ sind +512]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind];
						if ( threadIdx.x < 256  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +256];       
			}
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; 
						if ( threadIdx.x < 128  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +128];       
			}
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; 
						if ( threadIdx.x <  64) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +64];       
			}
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; 
						if ( threadIdx.x <  32) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +32];       
			}
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; 
						if ( threadIdx.x <  16) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +16];       
			}
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; 
						if ( threadIdx.x <   8) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +8];       
			}
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; 
						if ( threadIdx.x <   4) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +4];       
			}
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; 
						if ( threadIdx.x <   2) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +2];       
			}
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; 
						if ( threadIdx.x == 0 ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +1];       
			}
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[0] += bTheta[sind]; b_sh[0] += b_sh[sind]; b_sh2[0] += b_sh2[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 && tile == 0 ) {	Sum_bfx[blockIdx.x ] = bTheta[ 0 ];
							Sum_bfy[blockIdx.x ] = b_sh  [ 0 ];
							Sum_bfz[blockIdx.x ] = b_sh2  [ 0 ]; }
		if ( threadIdx.x == 0 && tile  > 0 ) {	Sum_bfx[blockIdx.x ] += bTheta[ 0 ];
							Sum_bfy[blockIdx.x ] += b_sh  [ 0 ];
							Sum_bfz[blockIdx.x ] += b_sh2  [ 0 ]; }
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

template __global__ void findCoefficientsLCS3D_cudaKernel_0deg ( int* powers_x, int* powers_y, int* powers_z, float* stationaryPoints_x, float* stationaryPoints_y, float* stationaryPoints_z, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* samplePoints_z, float* sampleValues_x, float* sampleValues_y, float* sampleValues_z, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bfx, float* Sum_bfy, float* Sum_bfz);
template __global__ void findCoefficientsLCS3D_cudaKernel_0deg ( int* powers_x, int* powers_y, int* powers_z, double* stationaryPoints_x, double* stationaryPoints_y, double* stationaryPoints_z, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* samplePoints_z, double* sampleValues_x, double* sampleValues_y, double* sampleValues_z, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bfx, double* Sum_bfy, double* Sum_bfz);



template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, size_t numOfSamplePoints, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bfx, T* Sum_bfy, T* Sum_bfz )
{

	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.y*blockDim.x;
	T* b_sh2	= b_sh   +blockDim.y*blockDim.x;
	T* bbT 		= b_sh2  +blockDim.y*blockDim.x;


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
			int power_z = powers_z[threadIdx.y];

			T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[ blockIdx.x ];
			T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[ blockIdx.x ];
			T dz = samplePoints_z[ pIndex ] -stationaryPoints_z[ blockIdx.x ];
			theta = weightingFunction( sqrt(dx*dx + dy*dy + dz*dz) );
			b = pow_i ( dx, power_x ) *pow_i( dy, power_y) *pow_i( dz, power_z);

			bTheta[ serialInBlockIndex ] = b*theta*sampleValues_x [ pIndex ];
			b_sh  [ serialInBlockIndex ] = b*theta*sampleValues_y [ pIndex ];
			b_sh2 [ serialInBlockIndex ] = b*theta*sampleValues_z [ pIndex ];
		}

		__syncthreads();


		size_t 	sind 	= serialInBlockIndex;
		int 	resid 	= min ( int(blockDim.x), int(numOfSamplePoints -tile*blockDim.x) ); // residual sample points
		size_t 	newReductionSize(1024);
		
		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind];
						if ( threadIdx.x < 512	) 			b_sh2 [ sind ]	+= b_sh2 [ sind +512]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind];
						if ( threadIdx.x < 256  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +256];       
			}
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; 
						if ( threadIdx.x < 128  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +128];       
			}
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; 
						if ( threadIdx.x <  64) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +64];       
			}
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; 
						if ( threadIdx.x <  32) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +32];       
			}
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; 
						if ( threadIdx.x <  16) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +16];       
			}
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; 
						if ( threadIdx.x <   8) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +8];       
			}
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; 
						if ( threadIdx.x <   4) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +4];       
			}
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; 
						if ( threadIdx.x <   2) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +2];       
			}
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; 
						if ( threadIdx.x == 0 ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +1];       
			}
			else newReductionSize = 1; __syncthreads();


			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[mIndex] += bTheta[sind]; b_sh[mIndex] += b_sh[sind]; b_sh2[mIndex] += b_sh2[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 && tile == 0 ) { 	Sum_bfx[ bfMonomIndex ] = bTheta[ mIndex ];
					  		Sum_bfy[ bfMonomIndex ] = b_sh  [ mIndex ]; 
					  		Sum_bfz[ bfMonomIndex ] = b_sh2 [ mIndex ]; }


		if ( threadIdx.x == 0 && tile > 0 )  { 	Sum_bfx[ bfMonomIndex ] += bTheta[ mIndex ];
					  		Sum_bfy[ bfMonomIndex ] += b_sh  [ mIndex ];
					  		Sum_bfz[ bfMonomIndex ] += b_sh2 [ mIndex ]; }
		__syncthreads();


		bTheta [ serialInBlockIndex ] 			= b*theta;
		b_sh   [ threadIdx.y +threadIdx.x*blockDim.y ] 	= b;

		__syncthreads();



		for ( size_t colsTile = 0; colsTile < (blockDim.y+blockDim.x-1)/blockDim.x; colsTile++ ) {

			size_t colIndex = threadIdx.x +colsTile*blockDim.x;

			if ( colIndex < blockDim.y ) {
				bbT[ rowIndex +colIndex ]  = b_sh[ colIndex ] *bTheta[ rowIndex2 ];
				for ( size_t pp = 1; pp < min( size_t(blockDim.x), size_t(numOfSamplePoints -tile*blockDim.x) ); pp++ )
					bbT[ rowIndex +colIndex ] += b_sh[ colIndex +pp*blockDim.y ] *bTheta[ rowIndex2 +pp ];


				if ( tile ==0 ) Sum_bbT[ colsTile*blockDim.x +bIndex ]  = bbT[ rowIndex +colIndex ];
				if ( tile > 0 ) Sum_bbT[ colsTile*blockDim.x +bIndex ] += bbT[ rowIndex +colIndex ];
			}
			__syncthreads();
		}

	}
}

template __global__ void findCoefficientsLCS3D_cudaKernel ( int* powers_x, int* powers_y, int* powers_z, float* stationaryPoints_x, float* stationaryPoints_y, float* stationaryPoints_z, size_t numOfSamplePoints, float* samplePoints_x, float* samplePoints_y, float* samplePoints_z, float* sampleValues_x, float* sampleValues_y, float* sampleValues_z, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bfx, float* Sum_bfy, float* Sum_bfz );
template __global__ void findCoefficientsLCS3D_cudaKernel ( int* powers_x, int* powers_y, int* powers_z, double* stationaryPoints_x, double* stationaryPoints_y, double* stationaryPoints_z, size_t numOfSamplePoints, double* samplePoints_x, double* samplePoints_y, double* samplePoints_z, double* sampleValues_x, double* sampleValues_y, double* sampleValues_z, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bfx, double* Sum_bfy, double* Sum_bfz );



template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel_0deg2 ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.x;
	T* b_sh2	= b_sh +blockDim.x;


	int pIndex = blockIdx.x*blockDim.x + threadIdx.x;
	{

		T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[blockIdx.y];
		T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[blockIdx.y];
		T dz = samplePoints_z[ pIndex ] -stationaryPoints_z[blockIdx.y];

		T theta = weightingFunction( sqrt(dx*dx + dy*dy + dz*dz) );

		bTheta[ threadIdx.x ] = theta*sampleValues_x [ pIndex ];
		b_sh  [ threadIdx.x ] = theta*sampleValues_y [ pIndex ];
		b_sh2 [ threadIdx.x ] = theta*sampleValues_z [ pIndex ];
		__syncthreads();


		size_t 	sind 	= threadIdx.x;
		int 	resid 	= blockDim.x;
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind];
						if ( threadIdx.x < 512	) 			b_sh2 [ sind ]	+= b_sh2 [ sind +512]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind];
						if ( threadIdx.x < 256  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +256];       
			}
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; 
						if ( threadIdx.x < 128  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +128];       
			}
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; 
						if ( threadIdx.x <  64) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +64];       
			}
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; 
						if ( threadIdx.x <  32) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +32];       
			}
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; 
						if ( threadIdx.x <  16) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +16];       
			}
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; 
						if ( threadIdx.x <   8) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +8];       
			}
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; 
						if ( threadIdx.x <   4) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +4];       
			}
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; 
						if ( threadIdx.x <   2) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +2];       
			}
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; 
						if ( threadIdx.x == 0 ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +1];       
			}
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[0] += bTheta[sind]; b_sh[0] += b_sh[sind]; b_sh2[0] += b_sh2[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 ) Sum_bf[blockIdx.x*3 +3*numOfStoredBlocks*blockIdx.y ] 	 += bTheta[ 0 ];
		if ( threadIdx.x == 0 )	Sum_bf[blockIdx.x*3 +3*numOfStoredBlocks*blockIdx.y +1 ] += b_sh  [ 0 ]; 
		if ( threadIdx.x == 0 )	Sum_bf[blockIdx.x*3 +3*numOfStoredBlocks*blockIdx.y +2 ] += b_sh2 [ 0 ]; 
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

template __global__ void findCoefficientsLCS3D_cudaKernel_0deg2 ( int* powers_x, int* powers_y, int* powers_z, float* stationaryPoints_x, float* stationaryPoints_y, float* stationaryPoints_z, float* samplePoints_x, float* samplePoints_y, float* samplePoints_z, float* sampleValues_x, float* sampleValues_y, float* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bf);
template __global__ void findCoefficientsLCS3D_cudaKernel_0deg2 ( int* powers_x, int* powers_y, int* powers_z, double* stationaryPoints_x, double* stationaryPoints_y, double* stationaryPoints_z, double* samplePoints_x, double* samplePoints_y, double* samplePoints_z, double* sampleValues_x, double* sampleValues_y, double* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bf);



template<class T>
__global__
void findCoefficientsLCS3D_cudaKernel2 ( int* powers_x, int* powers_y, int* powers_z, T* stationaryPoints_x, T* stationaryPoints_y, T* stationaryPoints_z, T* samplePoints_x, T* samplePoints_y, T* samplePoints_z, T* sampleValues_x, T* sampleValues_y, T* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<T> weightingFunction, T* Sum_bbT, T* Sum_bf)
{
	extern __shared__ __align__(sizeof(T)) unsigned char sharedMemory[];
	T* bTheta 	= reinterpret_cast<T*>(sharedMemory);
	T* b_sh		= bTheta +blockDim.y*blockDim.x;
	T* b_sh2	= b_sh +blockDim.y*blockDim.x;
	T* bbT 		= b_sh2 +blockDim.y*blockDim.x;


	size_t mIndex = threadIdx.y*blockDim.x;
	int serialInBlockIndex = threadIdx.x +mIndex;
	int pIndex = blockIdx.x*blockDim.x + threadIdx.x;

	{

		T dx = samplePoints_x[ pIndex ] -stationaryPoints_x[blockIdx.y];
		T dy = samplePoints_y[ pIndex ] -stationaryPoints_y[blockIdx.y];
		T dz = samplePoints_z[ pIndex ] -stationaryPoints_z[blockIdx.y];
		int power_x = powers_x[threadIdx.y];
		int power_y = powers_y[threadIdx.y];
		int power_z = powers_z[threadIdx.y];

		T theta = weightingFunction( sqrt(dx*dx + dy*dy + dz*dz) );
		T b = pow_i ( dx, power_x ) *pow_i( dy, power_y)* pow_i( dz, power_z);

		bTheta[ serialInBlockIndex ] = b*theta*sampleValues_x [ pIndex ];
		b_sh  [ serialInBlockIndex ] = b*theta*sampleValues_y [ pIndex ];
		b_sh2 [ serialInBlockIndex ] = b*theta*sampleValues_z [ pIndex ];
		__syncthreads();


		size_t 	sind 	= serialInBlockIndex;
		int 	resid 	= blockDim.x; // residual sample points
		size_t 	newReductionSize(1024);

		while ( resid > 0 ) {

			if (resid >=1024) {	if ( threadIdx.x < 512  ) 			bTheta[ sind ] 	+= bTheta[ sind +512];
						if ( threadIdx.x > 511 && threadIdx.x < 1024) 	b_sh  [ sind -512]+= b_sh  [ sind];
						if ( threadIdx.x < 512	) 			b_sh2 [ sind ]	+= b_sh2 [ sind +512]; }
			else newReductionSize = 512; __syncthreads();
	
			if (resid >=512) {	if ( threadIdx.x < 256  ) 			bTheta[ sind ] 	+= bTheta[ sind +256];
						if ( threadIdx.x > 255 && threadIdx.x < 512) 	b_sh  [ sind -256]+= b_sh  [ sind];
						if ( threadIdx.x < 256  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +256];       
			}
			else newReductionSize = 256; __syncthreads();
	
			if (resid >=256) {	if ( threadIdx.x < 128) 			bTheta[ sind ] 	+= bTheta[ sind +128];
						if ( threadIdx.x > 127 && threadIdx.x < 256) 	b_sh  [ sind -128]+= b_sh  [ sind]; 
						if ( threadIdx.x < 128  ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +128];       
			}
			else newReductionSize = 128; __syncthreads();

			if (resid >=128) {	if ( threadIdx.x <  64) 			bTheta[ sind ]	+= bTheta[ sind +64];
						if ( threadIdx.x >  63 && threadIdx.x < 128) 	b_sh  [ sind - 64]+= b_sh  [ sind]; 
						if ( threadIdx.x <  64) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +64];       
			}
			else newReductionSize = 64; __syncthreads();

			if (resid >= 64) {	if ( threadIdx.x <  32) 			bTheta[ sind ] 	+= bTheta[ sind +32];
			 			if ( threadIdx.x >  31 && threadIdx.x <  64) 	b_sh  [ sind - 32]+= b_sh  [ sind]; 
						if ( threadIdx.x <  32) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +32];       
			}
			else newReductionSize = 32; __syncthreads();
	
			if (resid >= 32) {	if ( threadIdx.x <  16) 			bTheta[ sind ] 	+= bTheta[ sind +16];
						if ( threadIdx.x >  15 && threadIdx.x <  32) 	b_sh  [ sind - 16]+= b_sh  [ sind]; 
						if ( threadIdx.x <  16) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +16];       
			}
			else newReductionSize = 16; __syncthreads();
	
			if (resid >= 16) {	if ( threadIdx.x <   8) 			bTheta[ sind ] 	+= bTheta[ sind +8];
						if ( threadIdx.x >   7 && threadIdx.x <  16) 	b_sh  [ sind -  8]+= b_sh  [ sind]; 
						if ( threadIdx.x <   8) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +8];       
			}
			else newReductionSize = 8; __syncthreads();
	
			if (resid >= 8 ) {	if ( threadIdx.x <   4) 			bTheta[ sind ] 	+= bTheta[ sind +4];
						if ( threadIdx.x >   3 && threadIdx.x <   8) 	b_sh  [ sind -  4]+= b_sh  [ sind]; 
						if ( threadIdx.x <   4) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +4];       
			}
			else newReductionSize = 4; __syncthreads();
	
			if (resid >= 4 ) {	if ( threadIdx.x <   2) 			bTheta[ sind ] 	+= bTheta[ sind +2];
						if ( threadIdx.x >   1 && threadIdx.x <   4) 	b_sh  [ sind -  2]+= b_sh  [ sind]; 
						if ( threadIdx.x <   2) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +2];       
			}
			else newReductionSize = 2; __syncthreads();
	
			if (resid >= 2 ) {	if ( threadIdx.x == 0 )				bTheta[ sind ] 	+= bTheta[ sind +1];
						if ( threadIdx.x == 1 )				b_sh  [ sind -  1]+= b_sh  [ sind]; 
						if ( threadIdx.x == 0 ) 			b_sh2 [ sind ] 	+= b_sh2 [ sind +1];       
			}
			else newReductionSize = 1; __syncthreads();

			if ( resid < blockDim.x && threadIdx.x == 0 )  { bTheta[mIndex] += bTheta[sind]; b_sh[mIndex] += b_sh[sind]; b_sh2[mIndex] += b_sh2[sind]; }
			__syncthreads();

			resid -= newReductionSize; sind += newReductionSize;
		}


		if ( threadIdx.x == 0 ) { Sum_bf[threadIdx.y +blockIdx.x*3*blockDim.y +3*blockDim.y*numOfStoredBlocks*blockIdx.y ] 		+= bTheta[ mIndex ];
					  Sum_bf[threadIdx.y +blockIdx.x*3*blockDim.y +3*blockDim.y*numOfStoredBlocks*blockIdx.y +blockDim.y ] 	+= b_sh  [ mIndex ]; 
					  Sum_bf[threadIdx.y +blockIdx.x*3*blockDim.y +3*blockDim.y*numOfStoredBlocks*blockIdx.y +blockDim.y ] 	+= b_sh2 [ mIndex ]; 
		}
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

template __global__ void findCoefficientsLCS3D_cudaKernel2 ( int* powers_x, int* powers_y, int* powers_z, float* stationaryPoints_x, float* stationaryPoints_y, float* stationaryPoints_z, float* samplePoints_x, float* samplePoints_y, float* samplePoints_z, float* sampleValues_x, float* sampleValues_y, float* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<float> weightingFunction, float* Sum_bbT, float* Sum_bf);
template __global__ void findCoefficientsLCS3D_cudaKernel2 ( int* powers_x, int* powers_y, int* powers_z, double* stationaryPoints_x, double* stationaryPoints_y, double* stationaryPoints_z, double* samplePoints_x, double* samplePoints_y, double* samplePoints_z, double* sampleValues_x, double* sampleValues_y, double* sampleValues_z, size_t numOfStoredBlocks, WeightingFunction<double> weightingFunction, double* Sum_bbT, double* Sum_bf);



