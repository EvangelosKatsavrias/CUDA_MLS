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


template<typename T>
void LLT_Cholesky_bandedInplace(T* A, int n, int n_b)
{
// Calculation of the rest L, D and U matrices elements
for (int i=0; i<n; i++)
{
	// Calculation of the L lower triangular matrix elements
	int m = std::max(1, i+1-n_b);
	for (int j = m; j<i-1; j++)
	{
		for (int k = m; k<j-1; k++) A[j*n+i-j+1] -= A[k*n+i-k+1] *A[k*n+j-k+1];
		A[j*n+i-j+1] /= A[j*n];
	}

	// Calculation of the diagonal L matrix elements
	for (int k=m;k<i-1;k++) A[i*n] -= pow(A[k*n+i-k+1],2);
	A[i*n] = sqrt( A[i*n] );
}
}


template void LLT_Cholesky_bandedInplace(float* A, int n, int n_b);
template void LLT_Cholesky_bandedInplace(double* A, int n, int n_b);


