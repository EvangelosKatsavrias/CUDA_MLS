#include"LUDecomposers.h"
#include<math.h>
#include<cmath>

template<typename T>
void LLT_Cholesky_inplace(T* A, int n)
{
// Calculation of the rest L, D and U matrices elements
for (int i=0; i<n; i++)
{
	// Calculation of the L lower triangular matrix elements
	for (int j = 0; j<=i-1; j++)
	{
		for (int k = 0; k<=j-1; k++) A[i*n+j] -= A[i*n+k] *A[j*n+k];
		A[i*n+j] /= A[j*n+j];
	}

	// Calculation of the diagonal L matrix elements
	for (int k=0;k<=i-1;k++) A[i*n+i] -= pow(A[i*n+k],2);
	A[i*n+i] = sqrt( A[i*n+i] );
}
}

template void LLT_Cholesky_inplace(float* A, int n);
template void LLT_Cholesky_inplace(double* A, int n);


