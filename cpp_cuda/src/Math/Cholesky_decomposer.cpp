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
void Cholesky_decomposer(T* A, T* L, int n)
{

for (int i=0; i<n; i++)
{
	// Calculation of the diagonal L matrix elements
	L[i*n+i] = A[i*n+i];
	for (int k = 0; k < i; k++) L[i*n+i] -= pow(L[i*n+k],2);
	L[i*n+i] = sqrt( L[i*n+i] );

	// Calculation of the L lower triangular matrix elements
	T invd = 1/L[i*n+i];
	for (int j = i+1; j<n; j++)
	{
	    L[j*n+i] = A[j*n+i];
	    for (int k = 0; k < i; k++) L[j*n+i] -= L[j*n+k] *L[i*n+k];
	    L[j*n+i] *= invd;
	}
}
}

template void Cholesky_decomposer(float* A, float* L, int n);
template void Cholesky_decomposer(double* A, double* L, int n);


