#include"gpuDeviceProperties.h"


int get_AllDevicesProperties(cudaDeviceProp* dev_prop)
{
	int dev_count; cudaGetDeviceCount(&dev_count);
	for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++) cudaGetDeviceProperties(&dev_prop[deviceIndex], deviceIndex);
	
	return dev_count;
}


void plotAllGPUsProperties()
{
  int dev_count; cudaGetDeviceCount(&dev_count);
  cudaDeviceProp dev_prop[dev_count];
  for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++) cudaGetDeviceProperties(&dev_prop[deviceIndex], deviceIndex);


  std::cout << "\n----------------------------\n| Number of GPU devices: " << dev_count << " |\n----------------------------\n" << std::endl;


  for (int deviceIndex = 0; deviceIndex < dev_count; deviceIndex++) 
  {
     std::cout << "\nProperties of GPU device nr. " << deviceIndex << "\n------------------------------" << std::endl;
     plotGPUProperties(dev_prop[deviceIndex]);
  }

}
