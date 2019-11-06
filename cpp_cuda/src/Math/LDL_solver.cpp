#include"LUSolvers.h"

template<typename T>
void LDL_solver(T* L, T* D, T* b, T* x, int n)
{
	if ( b != x ) for (int i =0; i<n; i++) x[i] = b[i];

	for (int i=0; i<n; i++) for (int j=0; j<i; j++) x[i] -= L[i*n+j]*x[j];

	for (int i=n-1; i>=0; i--) {
		x[i] /= D[i];
		for ( int j=n-1; j>=i+1; j--) x[i] -= L[j*n+i]*x[j];
	}
}


template void LDL_solver(float* L, float* D, float* b, float* x, int n);
template void LDL_solver(double* L, double* D, double* b, double* x, int n);

