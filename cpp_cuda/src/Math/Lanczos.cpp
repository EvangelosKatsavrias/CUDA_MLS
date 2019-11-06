/* Easily adaptable code for calculating the
Lanczos vectors of a given matrix. Extensions
can be made to improve numerical stability by
any
orthogonlization
procedure
due
to
the
inherent instabilities of finite precision. -
Travis Hoppe, 3.7.05 */

/*
#include <iostream>
#include <iomanip>
#include <cmath>
#include "nrutil.h"
using namespace std;
int N=3;

// Matrix operations
void MdV( double ** A, double * b, double * out) 
{
	int i,j;double value=0;
	for(j=0;j<N;j++){
		for(i=0,value=0;i<N;i++)
			value+=A[j][i]*b[i];
			out[j] = value; 
	}
}


void VTdM( double * bT, double ** A, double *out) 
{
	int i,j;double value=0;
	for(j=0;j<N;j++){
		for(i=0,value=0;i<N;i++)
			value+=A[i][j]*bT[i];
			out[j] = value; }
}

void VTdV( double *bT, double * c, double & out)
{
	out = 0;
	for(int i=0;i<N;i++) out += bT[i]*c[i]; 
}

void VmV( double *b, double & s, double *c) 
{
	for(int i=0;i<N;i++)
		b[i] -= s * c[i]; 
}

void VdS( double *b, double s) 
{
	for(int i=0;i<N;i++)
	b[i] *= s; 
}


// Various outputs
void printV( double * b ) 
{
	cout << " [";
	for(int i=0;i<N;i++) cout << setw(8) << b[i]; cout << "] "; 
}

void printV( double * b, int offset ) 
{
	cout << " [";
	for(int i=0+offset;i<N+offset;i++) cout << setw(8) << b[i]; cout << "] "; 
}
*/
/*
int main() 
{
	// allocate space for the variables
	int i,j,k;
	double ** A = dmatrix(0,N,0,N);
	double ** T = dmatrix(0,N,0,N);
	double *x[N+1],*y[N+1],*z[N+1],*w[N+1];
	double *alpha = dvector(0,N+1);
	double *beta = dvector(0,N+1);
	double *gamma = dvector(0,N+1);
	double delta;

	for(i=0;i<N+1;i++) {
		x[i] = dvector(0,N);
		y[i] = dvector(0,N);
		z[i] = dvector(0,N);
		w[i] = dvector(0,N); 
	}

	// set the sample matrix up
	A[0][0]=.5; A[0][1]=.5; A[0][2]=-.5;
	A[1][0]= 0; A[1][1]= 0; A[1][2]=- 2;
	A[2][0]=1.5; A[2][1]=-.5; A[2][2]=4.5;
	x[1][0]=1;
	x[1][1]=0;
	x[1][2]=-1;
	y[1][0]=.5; y[1][1]=-.5; y[1][2]=-.5;

	// calculate the Lanczos vectors
	for(k=1;k<N;k++) {
		MdV (A,x[k],z[k]);
		VTdM(y[k],A,w[k]);
		VTdV(y[k],z[k],alpha[k]);

		VmV (z[k],alpha[k],x[k]);
		VmV (z[k],gamma[k-1],x[k-1]);

		VmV (w[k],alpha[k],y[k]);
		VmV (w[k],beta[k-1] ,y[k-1]);

		VTdV(w[k],z[k],delta);
		beta[k] = sqrt(abs(delta));

		gamma[k] = beta[k];
		if(gamma[k]<0) gamma[k] *= -1;

		x[k+1] = z[k];
		VdS(x[k+1], 1.0/beta[k] );

		y[k+1] = w[k];
		VdS(y[k+1], 1.0/gamma[k] );

		cout << k << " -> " <<"x "; printV(x[k]);
		cout << "y "; printV(y[k]);
		cout << alpha[k] << ' ' << beta[k] << ' ' << gamma[k] << endl;

		cout << "     " << "v "; printV(z[k]);
		cout << "w "; printV(w[k]);
		cout << endl;

	}

	MdV (A,x[N],z[N]);
	VmV (z[N],alpha[N],x[N]);
	VTdV(y[N],z[N],alpha[N]);
	// display the output
	cout << k << " -> ";
	cout << "x "; printV(x[N]);
	cout << "y "; printV(y[N]);
	cout << alpha[N] << endl;
	cout << " -=-=-=-=-=-=-= " << endl;
	cout << " diagonal: "; printV(alpha,1); 
	cout << " super-di: "; printV(gamma,0); 
	cout << endl;
	return 0;
}
*/
