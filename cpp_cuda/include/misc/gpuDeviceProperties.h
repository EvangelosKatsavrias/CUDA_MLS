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
