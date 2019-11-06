#include"gpuDeviceProperties.h"


template<class Type>
void printProperty(std::string propertyDescription, Type propertyValue)
{
	int maxStringLength = 88;
        int numOfTabs = ceil((float)(maxStringLength - propertyDescription.length())/8); std::string tabs = "";
	for (int tabIndex = 0; tabIndex < numOfTabs; tabIndex++) tabs += "\t";

	std::cout << propertyDescription << tabs << propertyValue << std::endl;
}


void plotGPUProperties(cudaDeviceProp& dev_prop)
{
	std::cout << "\nDevice name: " << dev_prop.name << std::endl;
	std::vector<std::string> propertiesDescriptions; std::vector<int> propertiesValues; createDevPropertiesContainers(dev_prop, propertiesDescriptions, propertiesValues);

	std::vector<int>::iterator iter2 = propertiesValues.begin();
	for (std::vector<std::string>::iterator iter = propertiesDescriptions.begin(); iter != propertiesDescriptions.end(); iter++, iter2++) printProperty(*iter, *iter2);
}


void plotGPUBasicComputeResources(cudaDeviceProp& dev_prop)
{
	std::cout << "\nBasic compute resources: " << std::endl;
	std::vector<std::string> propertiesDescriptions; std::vector<int> propertiesValues; createDevPropertiesContainers(dev_prop, propertiesDescriptions, propertiesValues);

	std::vector<int> resourcesProps = {0, 1, 2, 3, 4, 5};
	for ( int index:resourcesProps ) printProperty(propertiesDescriptions[index], propertiesValues[index]);
}

void plotGPUMemoryResources(cudaDeviceProp& dev_prop)
{
	std::cout << "\nMemory resources: " << std::endl;
	std::vector<std::string> propertiesDescriptions; std::vector<int> propertiesValues; createDevPropertiesContainers(dev_prop, propertiesDescriptions, propertiesValues);

	std::vector<int> resourcesProps = {6, 7, 8, 9, 10, 11, 12};
	for ( int index:resourcesProps ) printProperty(propertiesDescriptions[index], propertiesValues[index]);
}

void plotGPUPerformanceProperties(cudaDeviceProp& dev_prop)
{
	std::cout << "\nPerformance properties: " << std::endl;
	std::vector<std::string> propertiesDescriptions; std::vector<int> propertiesValues; createDevPropertiesContainers(dev_prop, propertiesDescriptions, propertiesValues);

	std::vector<int> resourcesProps = {13, 14, 15, 16, 17, 18};
	for ( int index:resourcesProps ) printProperty(propertiesDescriptions[index], propertiesValues[index]);
}

void plotGPUGridSizeProperties(cudaDeviceProp& dev_prop)
{
	std::cout << "\nGrid size properties: " << std::endl;
	std::vector<std::string> propertiesDescriptions; std::vector<int> propertiesValues; createDevPropertiesContainers(dev_prop, propertiesDescriptions, propertiesValues);

	std::vector<int> resourcesProps = {19, 20, 21, 22, 23, 24};
	for ( int index:resourcesProps ) printProperty(propertiesDescriptions[index], propertiesValues[index]);
}
