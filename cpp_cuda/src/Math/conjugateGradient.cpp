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
