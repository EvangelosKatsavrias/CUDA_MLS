#include"matrix_memoryStorage.h"
#include<algorithm>
#include<iostream>

template<class T>
void convert_matrixStorage_2Skyline(T* A, int n, int* n_max, T* A_skyline)
{
	int k(0);
	for (int i =0; i<n; i++)
		for (int j = i; j > i -n_max[i+1] +n_max[i]; j--)
		{ A_skyline[ k ] = A[j +i*n]; k++; }
}


template void convert_matrixStorage_2Skyline(float* A, int n, int* n_max, float* A_skyline);
template void convert_matrixStorage_2Skyline(double* A, int n, int* n_max, double* A_skyline);
