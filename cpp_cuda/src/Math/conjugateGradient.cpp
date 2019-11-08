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

#include"conjugateGradient.h"
#include<iostream>
#include"myBLAS.h"
#include<math.h>
#include<cmath>


template void CG<double>( double* A, double* b, double* x, int n);
template void CG<float>	( float* A, float* b, float* x, int n);


template<class T>
void CG(T* A, T* b, T* x, int n)
{

	T epsilon = T(1e-4), epsilon2(1e-16), rnew[n], rold[n], p[n], bid[n], sum, sum1, sum2, residual, alpha, vita;

	for (int i = 0; i<n; i++) x[i] = 0;
	
	for (int i = 0; i<n; i++) {
		sum1 = 0;
		for (int j = 0; j<n; j++) sum1 += A[i*n+j] *x[j];
		rold[i] = b[i] -sum1;
		p[i] = rold[i];
	}

	
	for (int k = 0; k<n*n; k++) {  //     !   NCP iterations (at most)
	
		sum = 0; sum1 = 0;
	
		for (int i = 0; i<n; i++ ) {
			sum2 = 0;
			sum1 += rold[i]*rold[i];
			for (int j = 0; j<n; j++) sum2 += p[j]*A[j*n + i ];
			bid[i] = sum2;
			sum += sum2*p[i];
		}
		
		
//		std::cout << sum << "\t" << std::abs(sum) << "\t" << pow(epsilon, 2) << std::endl;

		if ( std::abs(sum) < epsilon2 ) {/* std::cout << "first exit at " << k << std::endl;*/ return;}  //  homogeneous system
		
		alpha = sum1/sum;    //   this is eta
	         
		for (int i = 0; i<n; i++) {
			x[i] += alpha*p[i];
			rnew[i] = rold[i] -alpha*bid[i];
		}
	
		if ( k%60 == 0 ) {//  !   regularly refresh residual 
			for (int i = 0; i<n; i++) {
				sum1 = 0;
				for (int j = 0; j<n; j++) sum1 += A[i*n+j]*x[j];
				rnew[i] = b[i] -sum1;
			}
		}
	
		residual = 0;
		for (int i = 0; i<n; i++) residual += rnew[i]*rnew[i];
		residual = sqrt(residual);
		if (residual < epsilon) {/* std::cout << "second exit at " << k << std::endl;*/ return;}
	
		sum1 = 0; sum2 = 0;
		for (int i = 0; i<n; i++) {
			sum1 += rold[i]*rold[i];
			sum2 += rnew[i]*rnew[i];
		}
	
		vita = sum2/sum1;
		for (int j = 0; j<n; j++) {
			p[j] = rnew[j] +vita*p[j];
			rold[j] = rnew[j];
		}

		//for ( int i=0;i<n;i++ ) std::cout << rold[i] << "\t"; std::cout << std::endl;
	}

	std::cout << "Residual = " << residual << std::endl;
	std::cout << "CG: Unable to Converge !!! \n";
}
