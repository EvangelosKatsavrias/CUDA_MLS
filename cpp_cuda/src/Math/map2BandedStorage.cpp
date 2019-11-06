#include"matrix_memoryStorage.h"


int map2BandedStorage(int i, int j, int n)
{ return (j-i)*n+i; }


int mapFromBandedStorage(int i_band, int n, int* i, int *j)
{
	int bandIndex = i_band /n; int inBandShift = i_band %n;
	*i = inBandShift; *j = inBandShift +bandIndex;
	return *j;
}
