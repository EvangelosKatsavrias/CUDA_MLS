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

#include"gaussElimination.h"
#include"linearizers.h"
#include<iostream>

int(*linearizer)(int row, int col, int numCols) = linearize_rowWiseZeroBased_2d;


template<class T>
void gaussElimination_backSubst(T* A, T* b, T* x, int n)
{

T* At;
 
if ( b != x ) { // then the solution will not be in place of the data memory place

	for (int i =0; i<n; i++) x[i] = b[i];

	At = new T[ n*n ];
	for (int i =0; i<n*n; i++) At[i] = A[i];
	A = At;
}


//-----------------Elimination phase-------------

for (int j = 0; j<(n-1); j++)
{
	T inv = 1/A[linearizer(j, j, n)];

	for (int i = j+1; i<n; i++)
	{
		T lambda = A[linearizer(j, i, n)]*inv;
		for (int k = j+1; k<n; k++) A[linearizer(k, i, n)] = A[linearizer(k, i, n)] - lambda*A[linearizer(k, j, n)];

		x[i] = x[i] - lambda*x[j];
	}
}


//--------------Back substitution phase-----------

for (int i=(n-1); i>-1; i--)
{
	for (int j=(n-1); j>i; j--)
	{
		x[i] = x[i] - A[linearizer(j, i, n)]*x[j];
	}

	x[i] = x[i]/A[linearizer(i, i, n)];
}

if ( b != x ) { delete[] At; }

}



template<class T>
void gaussElimination_frontSubst(T* A, T* b, T* x, int n)
{

T* At;

if ( b != x ) { // then the solution will not be in place of the data memory place

	for (int i =0; i<n; i++) x[i] = b[i];

	At = new T[ n*n ];
	for (int i =0; i<n*n; i++) At[i] = A[i];
	A = At;
}



//----------------- Elimination phase -------------

for (int j = (n-1); j>0; j--)
{
	float inv = 1/A[linearizer(j, j, n)];

	for (int i = j-1; i>-1; i--)
	{
		float lambda = A[linearizer(j, i, n)]*inv;

		for (int k = j-1; k>-1; k--) A[linearizer(k, i, n)] = A[linearizer(k, i, n)] - lambda*A[linearizer(k, j, n)];

		x[i] = x[i] - lambda*x[j];
	}
}


//-------------- Front substitution phase -----------

for (int i=0; i<n; i++)
{
	for (int j=0; j<i; j++)
	{
		x[i] = x[i] - A[linearizer(j, i, n)]*x[j];
	}

	x[i] = x[i]/A[linearizer(i, i, n)];
}


if ( b != x ) { delete[] At; }


}



template void gaussElimination_backSubst(float*, float*, float*, int);
template void gaussElimination_backSubst(double*, double*, double*, int);


template void gaussElimination_frontSubst(float*, float*, float*, int);
template void gaussElimination_frontSubst(double*, double*, double*, int);


