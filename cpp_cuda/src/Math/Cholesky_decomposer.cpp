#include"LUDecomposers.h"
#include<math.h>
#include<cmath>


template<typename T>
void Cholesky_decomposer(T* A, T* L, int n)
{

for (int i=0; i<n; i++)
{
	// Calculation of the diagonal L matrix elements
	L[i*n+i] = A[i*n+i];
	for (int k = 0; k < i; k++) L[i*n+i] -= pow(L[i*n+k],2);
	L[i*n+i] = sqrt( L[i*n+i] );

	// Calculation of the L lower triangular matrix elements
	T invd = 1/L[i*n+i];
	for (int j = i+1; j<n; j++)
	{
	    L[j*n+i] = A[j*n+i];
	    for (int k = 0; k < i; k++) L[j*n+i] -= L[j*n+k] *L[i*n+k];
	    L[j*n+i] *= invd;
	}
}
}

template void Cholesky_decomposer(float* A, float* L, int n);
template void Cholesky_decomposer(double* A, double* L, int n);


