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
#include<algorithm>

template<typename T>
void Cholesky_solver_banded(T* L, T* b, T* x, int n, int n_b)
{
	if ( b != x ) for (int i =0; i<n; i++) x[i] = b[i];

	for (int i=0; i<n; i++)	{
		int m = std::max(0, i-n_b );
		for (int j=m; j<i; j++) x[i] -= L[j +(i-j)*n ]*x[j];
		x[i] /= L[i];
	}

	for (int i=n-1; i>=0; i--) {
		for ( int j=std::min(n-1, i+n_b); j >= i+1; j--) x[i] -= L[i +(j-i)*n ]*x[j];
		x[i] /= L[i];
	}
}

template void Cholesky_solver_banded(float* L, float* b, float* x, int n, int n_b);
template void Cholesky_solver_banded(double* L, double* b, double* x, int n, int n_b);


