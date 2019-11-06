#include"gpuDeviceProperties.h"

void createDevPropertiesContainers(cudaDeviceProp& dev_prop, std::vector<std::string>& propertyDescriptions, std::vector<int>& propertyValue)
{
// ====== Basic compute resources =====
	propertyDescriptions.push_back("Major compute capability");			propertyValue.push_back(dev_prop.major);
	propertyDescriptions.push_back("Minor compute capability");			propertyValue.push_back(dev_prop.minor);

 	propertyDescriptions.push_back("Number of streaming multiprocessors (SM)");	propertyValue.push_back(dev_prop.multiProcessorCount);// decisive for the computational performance
	propertyDescriptions.push_back("Max threads per multiprocessor");		propertyValue.push_back(dev_prop.maxThreadsPerMultiProcessor);
	propertyDescriptions.push_back("Maximum permissible threads per block");	propertyValue.push_back(dev_prop.maxThreadsPerBlock);
	propertyDescriptions.push_back("Warp size");					propertyValue.push_back(dev_prop.warpSize);


// ====== Memory resources ========
	propertyDescriptions.push_back("Total global memory (Mbytes)");			propertyValue.push_back(dev_prop.totalGlobalMem/1000000);
	propertyDescriptions.push_back("Total constant memory (kbytes)");		propertyValue.push_back(dev_prop.totalConstMem/1000);
	propertyDescriptions.push_back("L2 cache memory size (kbytes)");		propertyValue.push_back(dev_prop.l2CacheSize/1000);
	propertyDescriptions.push_back("Size of shared memory per block (kbytes)");	propertyValue.push_back(dev_prop.sharedMemPerBlock/1000);
	propertyDescriptions.push_back("Shared memory per multiprocessor (kbytes)");	propertyValue.push_back(dev_prop.sharedMemPerMultiprocessor/1000);
	propertyDescriptions.push_back("Number of registers per block");		propertyValue.push_back(dev_prop.regsPerBlock);
	propertyDescriptions.push_back("Number of registers per multiprocessor");	propertyValue.push_back(dev_prop.regsPerMultiprocessor);


// ====== Performance properties ========
	propertyDescriptions.push_back("Clock frequency (MHz)");			propertyValue.push_back(dev_prop.clockRate/1000);// decisive for the computational performance
	propertyDescriptions.push_back("Memory clock rate (MHz)");			propertyValue.push_back(dev_prop.memoryClockRate/1000);
	propertyDescriptions.push_back("Memory bus width (bit)");			propertyValue.push_back(dev_prop.memoryBusWidth);
	propertyDescriptions.push_back("Memory pitch");					propertyValue.push_back(dev_prop.memPitch);
	propertyDescriptions.push_back("Global L1 cache memory, support");		propertyValue.push_back(dev_prop.globalL1CacheSupported);
	propertyDescriptions.push_back("Local L1 cache memory, support");		propertyValue.push_back(dev_prop.localL1CacheSupported);
//	propertyDescriptions.push_back("Single to double precision performance ratio");	propertyValue.push_back(dev_prop.singleToDoublePrecisionPerfRatio);


// ====== Grid size properties ======
	propertyDescriptions.push_back("Maximum permissible gridded threads per block, in dimension x");propertyValue.push_back(dev_prop.maxThreadsDim[0]);
	if(sizeof(dev_prop.maxThreadsDim) > sizeof(int)){propertyDescriptions.push_back("Maximum permissible gridded threads per block, in dimension y"); propertyValue.push_back(dev_prop.maxThreadsDim[1]);}
	if(sizeof(dev_prop.maxThreadsDim) > 2*sizeof(int)){propertyDescriptions.push_back("Maximum permissible gridded threads per block, in dimension z"); propertyValue.push_back(dev_prop.maxThreadsDim[2]);}
	propertyDescriptions.push_back("Maximum permissible grid size (thread blocks), in dimension x");	propertyValue.push_back(dev_prop.maxGridSize[0]);
	propertyDescriptions.push_back("Maximum permissible grid size (thread blocks), in dimension y");	propertyValue.push_back(dev_prop.maxGridSize[1]);
	propertyDescriptions.push_back("Maximum permissible grid size (thread blocks), in dimension z");	propertyValue.push_back(dev_prop.maxGridSize[2]);



// ====== Further device abilities ========
	propertyDescriptions.push_back("Integrated device");				propertyValue.push_back(dev_prop.integrated);
	propertyDescriptions.push_back("The device is multi-GPU board");		propertyValue.push_back(dev_prop.isMultiGpuBoard);
	propertyDescriptions.push_back("Multi-GPU board group-ID");			propertyValue.push_back(dev_prop.multiGpuBoardGroupID);
	propertyDescriptions.push_back("The device can map host memory");		propertyValue.push_back(dev_prop.canMapHostMemory);
	propertyDescriptions.push_back("Compute mode");					propertyValue.push_back(dev_prop.computeMode);
	propertyDescriptions.push_back("Concurrent kernels");				propertyValue.push_back(dev_prop.concurrentKernels);
	propertyDescriptions.push_back("Kernel execution timeout");			propertyValue.push_back(dev_prop.kernelExecTimeoutEnabled);
	propertyDescriptions.push_back("Data transfer overlap ability");		propertyValue.push_back(dev_prop.deviceOverlap);
	propertyDescriptions.push_back("Async engine count");				propertyValue.push_back(dev_prop.asyncEngineCount);
	propertyDescriptions.push_back("Unified addressing");				propertyValue.push_back(dev_prop.unifiedAddressing);
	propertyDescriptions.push_back("Stream priorities supported");			propertyValue.push_back(dev_prop.streamPrioritiesSupported);
	propertyDescriptions.push_back("Managed memory");				propertyValue.push_back(dev_prop.managedMemory);
	propertyDescriptions.push_back("ECC enabled");					propertyValue.push_back(dev_prop.ECCEnabled);
	propertyDescriptions.push_back("PCI bus ID");					propertyValue.push_back(dev_prop.pciBusID);
	propertyDescriptions.push_back("PCI device ID");				propertyValue.push_back(dev_prop.pciDeviceID);
	propertyDescriptions.push_back("PCI domain ID");				propertyValue.push_back(dev_prop.pciDomainID);
	propertyDescriptions.push_back("TCC driver");					propertyValue.push_back(dev_prop.tccDriver);
//	propertyDescriptions.push_back("Pageable memory access");			propertyValue.push_back(dev_prop.pageableMemoryAccess);
//	propertyDescriptions.push_back("Concurrent managed access");			propertyValue.push_back(dev_prop.concurrentManagedAccess);



// ===== graphics shader processor abilities ========
	propertyDescriptions.push_back("Max 1D texture size");						propertyValue.push_back(dev_prop.maxTexture1D);
	propertyDescriptions.push_back("Max 1D mipmapped texture size");				propertyValue.push_back(dev_prop.maxTexture1DMipmap);
	propertyDescriptions.push_back("Max 1D bound to linear memory texture size");			propertyValue.push_back(dev_prop.maxTexture1DLinear);
	propertyDescriptions.push_back("Max 1D layered texture dimensions, x");				propertyValue.push_back(dev_prop.maxTexture1DLayered[0]);
	propertyDescriptions.push_back("Max 1D layered texture dimensions, y");				propertyValue.push_back(dev_prop.maxTexture1DLayered[1]);
	propertyDescriptions.push_back("Max 2D texture dimensions, x");					propertyValue.push_back(dev_prop.maxTexture2D[0]);
	propertyDescriptions.push_back("Max 2D texture dimensions, y");					propertyValue.push_back(dev_prop.maxTexture2D[1]);
	propertyDescriptions.push_back("Max 2D mipmapped texture dimensions, x");			propertyValue.push_back(dev_prop.maxTexture2DMipmap[0]);
	propertyDescriptions.push_back("Max 2D mipmapped texture dimensions, y");			propertyValue.push_back(dev_prop.maxTexture2DMipmap[1]);
	propertyDescriptions.push_back("Max 2D texture dimensions bound to pitched memory, x");		propertyValue.push_back(dev_prop.maxTexture2DLinear[0]);
	propertyDescriptions.push_back("Max 2D texture dimensions bound to pitched memory, y");		propertyValue.push_back(dev_prop.maxTexture2DLinear[1]);
	propertyDescriptions.push_back("Max 2D texture dimensions bound to pitched memory, z");		propertyValue.push_back(dev_prop.maxTexture2DLinear[2]);
	propertyDescriptions.push_back("Max 2D texture dimensions under gather operations, x");		propertyValue.push_back(dev_prop.maxTexture2DGather[0]);
	propertyDescriptions.push_back("Max 2D texture dimensions under gather operations, y");		propertyValue.push_back(dev_prop.maxTexture2DGather[1]);
	propertyDescriptions.push_back("Max 2D layered texture dimensions, x");				propertyValue.push_back(dev_prop.maxTexture2DLayered[0]);
	propertyDescriptions.push_back("Max 2D layered texture dimensions, y");				propertyValue.push_back(dev_prop.maxTexture2DLayered[1]);
	propertyDescriptions.push_back("Max 2D layered texture dimensions, z");				propertyValue.push_back(dev_prop.maxTexture2DLayered[2]);
	propertyDescriptions.push_back("Max 3D texture dimensions, x");					propertyValue.push_back(dev_prop.maxTexture3D[0]);
	propertyDescriptions.push_back("Max 3D texture dimensions, y");					propertyValue.push_back(dev_prop.maxTexture3D[1]);
	propertyDescriptions.push_back("Max 3D texture dimensions, z");					propertyValue.push_back(dev_prop.maxTexture3D[2]);
	propertyDescriptions.push_back("Max 3D alternate texture dimensions, x");			propertyValue.push_back(dev_prop.maxTexture3DAlt[0]);
	propertyDescriptions.push_back("Max 3D alternate texture dimensions, y");			propertyValue.push_back(dev_prop.maxTexture3DAlt[1]);
	propertyDescriptions.push_back("Max 3D alternate texture dimensions, z");			propertyValue.push_back(dev_prop.maxTexture3DAlt[2]);
	propertyDescriptions.push_back("Max cube map texture dimensions");				propertyValue.push_back(dev_prop.maxTextureCubemap);
	propertyDescriptions.push_back("Max cube map layered texture dimensions, x");			propertyValue.push_back(dev_prop.maxTextureCubemapLayered[0]);
	propertyDescriptions.push_back("Max cube map layered texture dimensions, y");			propertyValue.push_back(dev_prop.maxTextureCubemapLayered[1]);
	propertyDescriptions.push_back("Alignment requirement for textures");				propertyValue.push_back(dev_prop.textureAlignment);
	propertyDescriptions.push_back("Pitch alignment requirement for texture references bound to pitched memory");		propertyValue.push_back(dev_prop.texturePitchAlignment);

 
	propertyDescriptions.push_back("Max 1D surface size");				propertyValue.push_back(dev_prop.maxSurface1D);
	propertyDescriptions.push_back("Max 1D layered surface dimensions, x");		propertyValue.push_back(dev_prop.maxSurface1DLayered[0]);
	propertyDescriptions.push_back("Max 1D layered surface dimensions, y");		propertyValue.push_back(dev_prop.maxSurface1DLayered[1]);
	propertyDescriptions.push_back("Max 2D surface dimensions, x");			propertyValue.push_back(dev_prop.maxSurface2D[0]);
	propertyDescriptions.push_back("Max 2D surface dimensions, y");			propertyValue.push_back(dev_prop.maxSurface2D[1]);
	propertyDescriptions.push_back("Max 3D surface dimensions, x");			propertyValue.push_back(dev_prop.maxSurface3D[0]);
	propertyDescriptions.push_back("Max 3D surface dimensions, y");			propertyValue.push_back(dev_prop.maxSurface3D[1]);
	propertyDescriptions.push_back("Max 3D surface dimensions, z");			propertyValue.push_back(dev_prop.maxSurface3D[2]);
	propertyDescriptions.push_back("Max 2D layered surface dimensions, x");		propertyValue.push_back(dev_prop.maxSurface2DLayered[0]);
	propertyDescriptions.push_back("Max 2D layered surface dimensions, y");		propertyValue.push_back(dev_prop.maxSurface2DLayered[1]);
	propertyDescriptions.push_back("Max 2D layered surface dimensions, z");		propertyValue.push_back(dev_prop.maxSurface2DLayered[2]);
	propertyDescriptions.push_back("Max cube map surface size");			propertyValue.push_back(dev_prop.maxSurfaceCubemap);
	propertyDescriptions.push_back("Max cube map layered surface dimensions, x");	propertyValue.push_back(dev_prop.maxSurfaceCubemapLayered[0]);
	propertyDescriptions.push_back("Max cube map layered surface dimensions, y");	propertyValue.push_back(dev_prop.maxSurfaceCubemapLayered[1]);
	propertyDescriptions.push_back("Surface alignment");				propertyValue.push_back(dev_prop.surfaceAlignment);

}


