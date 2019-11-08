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

#ifndef GPUDEVICEPROPERTIES
#define GPUDEVICEPROPERTIES

#include<cuda.h>
#include<iostream>
#include<string>
#include<vector>
#include<math.h>

//#include"helper_cuda.h"
/*
cudaSetDevice(dev);
cudaDeviceProp deviceProp;
cudaGetDeviceProperties(&deviceProp, dev);
_ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount);
*/

void createDevPropertiesContainers(cudaDeviceProp&, std::vector<std::string>&, std::vector<int>&);
int get_AllDevicesProperties(cudaDeviceProp* dev_prop);

template<class Type>
void printProperty(std::string propertyDescription, Type propertyValue);
void plotGPUProperties(cudaDeviceProp& dev_prop);
void plotAllGPUsProperties();
void plotGPUBasicComputeResources(cudaDeviceProp& dev_prop);
void plotGPUMemoryResources(cudaDeviceProp& dev_prop);
void plotGPUPerformanceProperties(cudaDeviceProp& dev_prop);
void plotGPUGridSizeProperties(cudaDeviceProp& dev_prop);


#endif
