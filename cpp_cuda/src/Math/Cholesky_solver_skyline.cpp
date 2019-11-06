#include"LUSolvers.h"
#include<algorithm>

template<typename T>
void Cholesky_solver_skyline(T* L, T* b, T* x, int n, int* n_max )
{
	if ( b != x ) for (int i =0; i<n; i++) x[i] = b[i];

	for (int i=0; i<n; i++)	{
		for (int j=i-(n_max[i+1]-n_max[i]-1); j<i; j++) x[i] -= L[n_max[i]+i-j]*x[j];
		x[i] /= L[n_max[i]];
	}

	for (int i=n-1; i>=0; i--) {
		for ( int j=n-1; j>=i+1; j--) if (j-(n_max[j+1]-n_max[j]-1)<=i) x[i] -= L[n_max[j]+j-i]*x[j];
		x[i] /= L[n_max[i]];
	}
}


template void Cholesky_solver_skyline(float* L, float* b, float* x, int n, int* n_max);
template void Cholesky_solver_skyline(double* L, double* b, double* x, int n, int* n_max);


