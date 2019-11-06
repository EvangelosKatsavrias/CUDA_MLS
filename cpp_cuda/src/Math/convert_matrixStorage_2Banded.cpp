#include"matrix_memoryStorage.h"
#include<algorithm>
#include<iostream>

template<class T>
void convert_matrixStorage_2Banded(T* A, int n, int numOfBands, T* A_band)
{
	for (int i =0; i<numOfBands*n; i++) A_band[ i ] = 0;
	for (int i =0; i<n; i++) for (int j= i; j < std::min(i+numOfBands, n); j++) A_band[ map2BandedStorage(i,j,n) ] = A[i*n+j];
}


template void convert_matrixStorage_2Banded(float* A, int n, int numOfBands, float* A_band);
template void convert_matrixStorage_2Banded(double* A, int n, int numOfBands, double* A_band);
