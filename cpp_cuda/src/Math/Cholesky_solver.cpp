#include"LUSolvers.h"

template<typename T>
void Cholesky_solver(T* L, T* b, T* x, int n)
{
	if ( b != x ) for (int i =0; i<n; i++) x[i] = b[i];

	for ( int i = 0; i < n; i++ ) {
		for ( int j = 0; j < i; j++ ) x[i] -= L[i*n+j] *x[j];
		x[i] /= L[i*n+i];
	}

	for ( int i=n-1; i>=0; i-- ) {
		for ( int j=n-1; j>=i+1; j--) x[i] -= L[j*n+i]*x[j];
		x[i] /= L[i*n+i];
	}
}

template void Cholesky_solver(float* L, float* b, float * x, int n);
template void Cholesky_solver(double* L, double* b, double* x, int n);


