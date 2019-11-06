#include"LUSolvers.h"
#include<algorithm>

template<typename T>
void Cholesky_solver_banded(T* L, T* b, T* x, int n, int n_b)
{
	if ( b != x ) for (int i =0; i<n; i++) x[i] = b[i];

	for (int i=0; i<n; i++)	{
		int m = std::max(0, i-n_b );
		for (int j=m; j<i; j++) x[i] -= L[j +(i-j)*n ]*x[j];
		x[i] /= L[i];
	}

	for (int i=n-1; i>=0; i--) {
		for ( int j=std::min(n-1, i+n_b); j >= i+1; j--) x[i] -= L[i +(j-i)*n ]*x[j];
		x[i] /= L[i];
	}
}

template void Cholesky_solver_banded(float* L, float* b, float* x, int n, int n_b);
template void Cholesky_solver_banded(double* L, double* b, double* x, int n, int n_b);


