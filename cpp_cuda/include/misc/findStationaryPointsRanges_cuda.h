#include<cuda.h>

template<class T>
__global__
void find_samplePointsRanges(int numOfStationaryPoints, T *stationaryPoints, int numOfSamplePoints, T* samplePoints, T span, int *ranges, int* maxNumOfEffectiveSamplePoints);
