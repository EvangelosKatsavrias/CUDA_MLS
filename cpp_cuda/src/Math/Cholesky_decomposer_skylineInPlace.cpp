#include"LUDecomposers.h"
#include<math.h>
#include<cmath>
#include<algorithm>
#include<iostream>

template<typename T>
void Cholesky_decomposer_skylineInPlace(T* A, int n, int* n_max)
{

// Calculation of the rest L, D and U matrices elements
for (int i=0; i<n; i++)
{
	// Calculation of the L lower triangular matrix elements
	int mi = i+1 -(n_max[i+1] -n_max[i]);
	for (int j = mi; j < i; j++)
	{
		int mj = j+1 - (n_max[j+1] -n_max[j]);
		int m = std::max(mi, mj);
		for (int k = m; k < j; k++) A[ n_max[i]+i-j] -= A[ n_max[i]+i-k] *A[ n_max[j]+j-k];
		A[ n_max[i]+i-j] /= A[ n_max[j] ];
	}

	// Calculation of the diagonal L matrix elements
	for (int k=mi; k< i; k++) A[ n_max[i] ] -= pow(A[ n_max[i]+i-k],2);
	A[n_max[i]] = sqrt( A[n_max[i]] );
}
}

template void Cholesky_decomposer_skylineInPlace(float* A, int n, int* n_max);
template void Cholesky_decomposer_skylineInPlace(double* A, int n, int* n_max);
