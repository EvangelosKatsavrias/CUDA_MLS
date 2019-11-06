#include"LUDecomposers.h"
#include<math.h>
#include<cmath>

template<typename T>
void LDL_decomposer(T* A, T* L, T* D, int n)
{

// Calculation of the rest L, D and U matrices elements
for (int i=0; i<n; i++)
{
	// Calculation of the D diagonal matrix elements
	D[i] = A[i*n+i];
	
	for (int k=0;k<=i-1;k++) D[i] -= D[k] *pow(L[i*n+k],2);
	
	// Calculation of the L lower triangular matrix elements
	T invd = 1/D[i];
	for (int j = i+1; j<n; j++)
	{
	    L[j*n+i] = A[j*n+i];
	    for (int k= 0; k<=i-1; k++) L[j*n+i] -= D[k]*L[j*n+k] *L[i*n+k];
	    L[j*n+i] *= invd;
	}
}

// Filling the main diagonal of the lower triangular matrix L with ones
for (int i = 0; i<n; i++) L[i*n+i] = 1;
}

template void LDL_decomposer(float* A, float* L, float* D, int n);
template void LDL_decomposer(double* A, double* L, double* D, int n);


