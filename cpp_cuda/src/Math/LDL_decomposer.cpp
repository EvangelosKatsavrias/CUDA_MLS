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
#include<math.h>
#include<cmath>

template<typename T>
void LDL_decomposer(T* A, T* L, T* D, int n)
{

// Calculation of the rest L, D and U matrices elements
for (int i=0; i<n; i++)
{
	// Calculation of the D diagonal matrix elements
	D[i] = A[i*n+i];
	
	for (int k=0;k<=i-1;k++) D[i] -= D[k] *pow(L[i*n+k],2);
	
	// Calculation of the L lower triangular matrix elements
	T invd = 1/D[i];
	for (int j = i+1; j<n; j++)
	{
	    L[j*n+i] = A[j*n+i];
	    for (int k= 0; k<=i-1; k++) L[j*n+i] -= D[k]*L[j*n+k] *L[i*n+k];
	    L[j*n+i] *= invd;
	}
}

// Filling the main diagonal of the lower triangular matrix L with ones
for (int i = 0; i<n; i++) L[i*n+i] = 1;
}

template void LDL_decomposer(float* A, float* L, float* D, int n);
template void LDL_decomposer(double* A, double* L, double* D, int n);


