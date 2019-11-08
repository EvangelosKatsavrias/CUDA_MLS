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
void LDL_Cholesky_inplace(T* A, T* D, int n)
{

// Calculation of the rest L, D and U matrices elements
for (int i=0; i<n; i++)
{
    // Calculation of the L lower triangular matrix elements
    for (int j = 0; j<=i-1; j++)
    {
        for (int k= 0; k<=j-1; k++) A[i*n+j] -= D[k]*A[i*n+k] *A[j*n+k];
        A[i*n+j] = A[i*n+j]/D[j];
    }


    // Calculation of the D diagonal matrix elements
    D[i] = A[i*n+i];
    for (int k=0;k<=i-1;k++) D[i] -= D[k] *pow(A[i*n+k],2);

}

}

template void LDL_Cholesky_inplace(float* A, float* D, int n);
template void LDL_Cholesky_inplace(double* A, double* D, int n);
