//   CUDA_MLS Framework
//
//   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
//
//   This file is part of the CUDA_MLS Framework.
//
//   CUDA_MLS Framework is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License version 3 as published by
//   the Free Software Foundation.
//
//   CUDA_MLS Framework is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
//
//   Contact Info:
//   Evangelos D. Katsavrias
//   email/skype: vageng@gmail.com
// -----------------------------------------------------------------------

#include"LUDecomposers.h"

template<typename T>
void LU_Doolittle_decomposer(T* A, T* L, T* U, int n)
{

for (int i=0;i<n;i++)
{
	// Calculation of the L lower triangular matrix elements
	for (int j=0;j<=i-1;j++)
	{
		L[i*n+j] = A[i*n+j];
		for (int k=0;k<=j-1;k++) L[i*n+j] -= L[i*n+k]*U[k*n+j];
		L[i*n+j] /= U[j*n+j];
	}
	
	// Calculation of the U upper triangular matrix elements
	for (int j=i;j<n;j++)
	{
		U[i*n+j] = A[i*n+j];
		for (int k=0;k<=i-1;k++) U[i*n+j] -= L[i*n+k]*U[k*n+j];
	}
	
	// Filling the main diagonal of L with ones
	L[i*n+i]=1;
}
}


template void LU_Doolittle_decomposer(float* A, float* L, float* U, int n);
template void LU_Doolittle_decomposer(double* A, double* L, double* U, int n);
