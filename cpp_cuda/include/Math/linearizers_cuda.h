//   CUDA_MLS Framework
//
//   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
//
//   This file is part of the CUDA_MLS Framework.
//
//   CUDA_MLS Framework is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License version 3 as published by
//   the Free Software Foundation.
//
//   CUDA_MLS Framework is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
//
//   Contact Info:
//   Evangelos D. Katsavrias
//   email/skype: vageng@gmail.com
// -----------------------------------------------------------------------

#include<cuda.h>
/*
template<class T>
__device__ T linearize_rowWiseZeroBased_2d(T i, T j, T size_x);

template<class T>
__device__ T linearize_rowWiseZeroBased_3d(T i, T j, T k, T size_x, T size_y);

template<class T>
__device__ T linearizeBlockId_in3dGrid();

template<class T>
__device__ T linearizeBlockId_in2dGrid();

template<class T>
__device__ T linearizeThreadId_in3dBlock();

template<class T>
__device__ T linearizeThreadId_in2dBlock();

template<class T>
__device__ T threadsShift_in3dGrid();

template<class T>
__device__ T threadsShift_in2dGrid();

template<class T>
__device__ T linearizeThreadId_in3dGrid();

template<class T>
__device__ T linearizeThreadId_in2dGrid();
*/

template<class T>
__device__ T linearize_rowWiseZeroBased_2d(T i, T j, T size_x)
{
	return i + j*size_x;
}

template<class T>
__device__ T linearize_rowWiseZeroBased_3d(T i, T j, T k, T size_x, T size_y)
{
	return i + j*size_x + k*size_x*size_y;
}


// ============ CUDA linearizers ================================================

template<class T>
__device__ T linearizeBlockId_in3dGrid()
{
	return linearize_rowWiseZeroBased_3d(blockIdx.x, blockIdx.y, blockIdx.z, gridDim.x, gridDim.y);
}


template<class T>
__device__ T linearizeBlockId_in2dGrid()
{
	return linearize_rowWiseZeroBased_2d<T>(blockIdx.x, blockIdx.y, gridDim.x);
}


template<class T>
__device__ T linearizeThreadId_in3dBlock()
{
	return linearize_rowWiseZeroBased_3d<T>(threadIdx.x, threadIdx.y, threadIdx.z, blockDim.x, blockDim.y);
}


template<class T>
__device__ T linearizeThreadId_in2dBlock()
{
	return linearize_rowWiseZeroBased_2d<T>(threadIdx.x, threadIdx.y, blockDim.x);
}

template<class T>
__device__ T threadsShift_in3dGrid()
{
	return linearizeBlockId_in3dGrid<T>()*blockDim.x*blockDim.y*blockDim.z;
}

template<class T>
__device__ T threadsShift_in2dGrid()
{
	return linearizeBlockId_in2dGrid<T>()*blockDim.x*blockDim.y;
}


template<class T>
__device__ T linearizeThreadId_in3dGrid()
{
	return threadsShift_in3dGrid<T>() + linearizeThreadId_in3dBlock<T>();
}

template<class T>
__device__ T linearizeThreadId_in2dGrid()
{
	return threadsShift_in2dGrid<T>() + linearizeThreadId_in2dBlock<T>();
}


