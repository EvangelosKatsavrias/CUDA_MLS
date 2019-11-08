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
#include<algorithm>
#include<iostream>

template<typename T>
void Cholesky_decomposer_skylineInPlace(T* A, int n, int* n_max)
{

// Calculation of the rest L, D and U matrices elements
for (int i=0; i<n; i++)
{
	// Calculation of the L lower triangular matrix elements
	int mi = i+1 -(n_max[i+1] -n_max[i]);
	for (int j = mi; j < i; j++)
	{
		int mj = j+1 - (n_max[j+1] -n_max[j]);
		int m = std::max(mi, mj);
		for (int k = m; k < j; k++) A[ n_max[i]+i-j] -= A[ n_max[i]+i-k] *A[ n_max[j]+j-k];
		A[ n_max[i]+i-j] /= A[ n_max[j] ];
	}

	// Calculation of the diagonal L matrix elements
	for (int k=mi; k< i; k++) A[ n_max[i] ] -= pow(A[ n_max[i]+i-k],2);
	A[n_max[i]] = sqrt( A[n_max[i]] );
}
}

template void Cholesky_decomposer_skylineInPlace(float* A, int n, int* n_max);
template void Cholesky_decomposer_skylineInPlace(double* A, int n, int* n_max);
