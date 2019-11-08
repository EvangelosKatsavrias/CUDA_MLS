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

