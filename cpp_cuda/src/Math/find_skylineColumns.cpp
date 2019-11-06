#include"matrix_memoryStorage.h"
#include<algorithm>
#include<iostream>

template<class T>
int find_skylineColumns(T* A, int n, int* n_max)
{

	n_max[0] = 0; n_max[1] = 1;

	for (int i = 1; i<n; i++)
		for (int j = 0; j <i+1; j++)
			if ( A[i*n +j] != 0 ) { n_max[i+1] = n_max[i] +i -j +1; break; }

}


template int find_skylineColumns(float* A, int n, int* n_max);
template int find_skylineColumns(double* A, int n, int* n_max);
