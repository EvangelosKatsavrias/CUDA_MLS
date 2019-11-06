#include"LUSolvers.h"

template<typename T>
void LDU_solver(T* L, T* U, T* D, T* b, T* x, int n)
{
	if ( b != x ) for (int i =0; i<n; i++) x[i] = b[i];

	for (int i=0; i<n; i++)	{
		for (int j=0; j<i; j++) x[i] -= L[i*n+j]*x[j];
		x[i] /= L[i*n+i];
	}

	for (int i=n-1; i>=0; i--) {
		x[i] /= D[i];
		for ( int j=n-1; j>=i+1; j--) x[i] -= U[i*n+j]*x[j];
	}
}

template void LDU_solver(float* L, float* U, float* D, float* b, float* x, int n);
template void LDU_solver(double* L, double* U, double* D, double* b, double* x, int n);

