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

#include"matrix_memoryStorage.h"
#include<algorithm>
#include<iostream>

template<class T>
int find_skylineColumns(T* A, int n, int* n_max)
{

	n_max[0] = 0; n_max[1] = 1;

	for (int i = 1; i<n; i++)
		for (int j = 0; j <i+1; j++)
			if ( A[i*n +j] != 0 ) { n_max[i+1] = n_max[i] +i -j +1; break; }

}


template int find_skylineColumns(float* A, int n, int* n_max);
template int find_skylineColumns(double* A, int n, int* n_max);
