#include"LUDecomposers.h"

template<typename T>
void LU_decomposer_inPlace(T* A, int n)
{
	for (int i=0;i<n;i++)
	{
		for (int j=i+1;j<n;j++)
		{
			A[j*n+i] /= A[i*n+i];
			for (int k=i+1; k<n; k++) A[j*n+k] -= A[j*n+i]*A[i*n+k];
		}
	}
}

template void LU_decomposer_inPlace(float* A, int n);
template void LU_decomposer_inPlace(double* A, int n);
