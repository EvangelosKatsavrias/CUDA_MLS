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
