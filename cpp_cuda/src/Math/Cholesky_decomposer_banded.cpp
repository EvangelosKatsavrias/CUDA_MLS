#include"LUDecomposers.h"
#include<math.h>
#include<cmath>
#include<algorithm>


template<typename T>
void Cholesky_decomposer_banded(T* A, T* L, int n, int n_b)
{

for (int i=0; i<n; i++)
{
	// Calculation of the diagonal L matrix elements
	L[i] = A[i]; 

	for (int k = 1; k < std::min(i+1, n_b); k++) L[i] -= pow(L[i-k +k*n],2); 
	L[i] = sqrt( L[i] );

	// Calculation of the L lower triangular matrix elements
	T invd = 1/L[i];
	for (int j = 1; j < std::min(n-i+1, n_b); j++)
	{
		L[ i +j*n ] = A[ i +j*n ];
		for (int k = 1; k < std::min(i+1, n_b); k++)
			if (k+j+1<=n_b) L[ i +j*n ] -= L[ i-k +k*n ] *L[ i-k +(k+j)*n ];
		L[ i +j*n ] *= invd;
	}
}
}

template void Cholesky_decomposer_banded(float* A, float* L, int n, int n_b);
template void Cholesky_decomposer_banded(double* A, double* L, int n, int n_b);
