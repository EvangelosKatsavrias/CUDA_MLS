#include"LUDecomposers.h"

template<typename T>
void LDU_decomposer(T* A, T* L, T* U, T* D, int n)
{

for (int i=0;i<n*n;i++) {L[i] = 0; U[i] = 0;}

for (int i=0;i<n;i++)
{
	// Calculation of the L lower triangular matrix elements
	for (int j=0;j<=i-1;j++)
	{
		L[i*n+j] = A[i*n+j];
		for (int k=0;k<=j-1;k++) L[i*n+j] -= D[k]*L[i*n+k]*U[k*n+j];
		L[i*n+j] /= D[j];
	}
	
	// Calculation of the D diagonal matrix elements
	D[i] = A[i*n+i];
	for (int k=0;k<=i-1;k++) D[i] -= D[k]*L[i*n+k]*U[k*n+i];
	
	// Calculation of the U lower triangular matrix elements
	T invd = 1/D[i];
	for (int j=i+1;j<n;j++)
	{
		U[i*n+j] = A[i*n+j];
		for (int k=0;k<=i-1;k++) U[i*n+j] -= D[k]*L[i*n+k]*U[k*n+j];
		U[i*n+j] *= invd;
	}
	
	// Filling the main diagonal of L and U with ones
	U[i*n+i]=1; L[i*n+i] = 1;
}
}


template void LDU_decomposer(float* A, float* L, float* U, float* D, int n);
template void LDU_decomposer(double* A, double* L, double* U, double* D, int n);
