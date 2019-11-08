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

#include"LUSolvers.h"

template<typename T>
void LDL_solver(T* L, T* D, T* b, T* x, int n)
{
	if ( b != x ) for (int i =0; i<n; i++) x[i] = b[i];

	for (int i=0; i<n; i++) for (int j=0; j<i; j++) x[i] -= L[i*n+j]*x[j];

	for (int i=n-1; i>=0; i--) {
		x[i] /= D[i];
		for ( int j=n-1; j>=i+1; j--) x[i] -= L[j*n+i]*x[j];
	}
}


template void LDL_solver(float* L, float* D, float* b, float* x, int n);
template void LDL_solver(double* L, double* D, double* b, double* x, int n);

