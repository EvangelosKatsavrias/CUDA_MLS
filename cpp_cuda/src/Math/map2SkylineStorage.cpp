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

#include<limits>
#include<numeric>

template<typename T>
void find_skyline_mapVector( T* A, int* maxa, int n)
{
	maxa[0] = 1; maxa[1] = 2;
	for ( int j=2; j<n; j++ ) for ( int i=1; i<j; i++ ) if ( A[i*n+j] != 0 ) { maxa[j] = maxa[j-1] +j-i+1; break; }
}

template void find_skyline_mapVector( double* A, int* maxa, int n);
template void find_skyline_mapVector( float* A, int* maxa, int n);


int map2SkylineStorage(int i, int j, int* maxa)
{
	if (i > j) { int temp(i); i=j; j=temp; }; // if the queried element belongs to the lower triangular part of the matrix, consider the corresponding reciprocal element
	int firstElement = j -( maxa[j+1] -maxa[j] ); // the row of the most top non zero element of the queried matrix column
	if ( firstElement > i ) return std::numeric_limits<int>::quiet_NaN();
	return maxa[j] +(j-i);
}
