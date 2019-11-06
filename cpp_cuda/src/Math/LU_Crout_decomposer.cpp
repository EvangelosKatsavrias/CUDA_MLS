#include"LUDecomposers.h"

template<typename T>
void LU_Crout_decomposer(T* A, T* L, T* U, int n)
{

for (int i=0;i<n;i++)
{
	// Calculation of the L lower triangular matrix elements
	for (int j=0;j<=i;j++)
	{
		L[i*n+j] = A[i*n+j];
		for (int k=0;k<=j-1;k++) L[i*n+j] -= L[i*n+k]*U[k*n+j];
	}
	
	
	// Calculation of the U upper triangular matrix elements
	T invd = 1/L[i*n+i];
	for (int j=i+1;j<n;j++)
	{
		U[i*n+j] = A[i*n+j];
		for (int k=0;k<=i-1;k++) U[i*n+j] -= L[i*n+k]*U[k*n+j];
		U[i*n+j] *= invd;
	}
	
	// Filling the main diagonal of U with ones
	U[i*n+i]=1;

}

}

template void LU_Crout_decomposer(float* A, float* L, float* U, int n);
template void LU_Crout_decomposer(double* A, double* L, double* U, int n);

