#include<math.h>
#include<cmath>
#include"myBLAS.h"
#include"LUDecomposers.h"
#include<iostream>

template<class T>
void initializeVector(T* start, T* end, T initialValue)
{
	while (start != end) *(start++) = initialValue;
}


template<class T>
void multiplicationOfVectorbyScalar(T* start, T* end, T* result, T scalar)
{
	while (start != end) *(result++) = *(start++)*scalar;
}


template<class T>
void tensorProductOfVectors(T* start, T* end, T* start_2, T* end_2, T* result)
{
	int vectorSize = end_2-start_2; int iter(0);
	while (start != end)
	{
		multiplicationOfVectorbyScalar(start_2, end_2, result+vectorSize*(iter++), *(start++));
	}
}


template<class T>
void sumOfVectors(T* start, T* end, T* start_2, T* result)
{
	while (start != end) *(result++) = *(start++) + *(start_2++);
}


template<class T>
T dotProductOfVectors(T* start, T* end, T* start_2)
{
	T result(0);
	while (start != end) result += *(start++) * *(start_2++);

	return result;
}


template<class T>
void linspace(T start, T end, int numOfElements, T* output)
{
	T step = (end - start) /(numOfElements-1);
	for (int index = 0; index < numOfElements; index++) {output[index] = start; start += step;}

}

/*
template<class T>
T pow_i( T base, int exp ) {
	T result = 1;
	while (1) {
		if (exp & 1) result *= base; exp >>= 1;
		if (!exp) break; base *= base; }
	return result;
}
*/

int highestPowerof2(unsigned int n) {
	// Invalid input
	if (n < 1)
	    return 0;
	
	int res = 1;
	
	// Try all powers starting from 2^1
	for (int i=0; i<8*sizeof(unsigned int); i++)
	{
	    int curr = 1 << i;
	
	    // If current power is more than n, break
	    if (curr > n)
	       break;
	
	    res = curr;
	}
	
	return res;
}




int factorial(int n)
{ return (n==0 || n==1) ? 1:n*factorial(n-1); }

int choose_k_outOf_n(int k, int n)
{ return factorial(n)/(factorial(n-k)*factorial(k)); }

int numOfMonomials_ofDegree(int dim, int deg )
{ return choose_k_outOf_n( deg, dim+deg-1 ); }

int numOfMonomials_upToDegree(int dim, int deg)
{ return choose_k_outOf_n( deg, dim+deg ); }

template<typename T>
void evaluateMonomials2D(int deg, T x, T y, T* eval)
{
	int count(0);
	for (int i = 0; i<=deg; i++) for (int j = 0; j<=deg-i; j++) { eval[count] = pow_i(x, j)*pow_i(y, i); count++; }
}


template<typename T>
void evaluateMonomials3D(int deg, T x, T y, T z, T* eval)
{
	int count(0);
	for (int i = 0; i<=deg; i++)
		for (int j = 0; j<=deg-i; j++)
			for (int k = 0; k<=deg-i-j; k++) { 
				eval[count] = pow_i(x, k)*pow_i(y, j)*pow_i(z, i); count++; }
}


template<typename T>
T distance2D(T x1, T x2, T y1, T y2) {
	return sqrt( (x1 -x2)*(x1 -x2) +(y1 -y2)*(y1 -y2) );
}


template<typename T>
T distance3D(T x1, T x2, T y1, T y2, T z1, T z2) {
	return sqrt( (x1 -x2)*(x1 -x2) +(y1 -y2)*(y1 -y2) +(z1 -z2)*(z1 -z2) );
}


template<typename T>
T norm(T* x, int dim, int normType )
{
	T res(0);
	for (int comp=0; comp<dim; comp++ ) res += pow_i( fabs( x[comp] ), normType );

	return pow( res, 1/normType );
}


template<typename T>
T determinant(T* A, int numCols, bool luFlag)
{
	if ( !luFlag )
	{
		T U[numCols*numCols];
		for (int i=0;i<numCols;i++)
			for (int j=0;j<numCols;j++) U[i*numCols+j] = A[i*numCols+j];

		for (int i=0;i<numCols;i++)
		{
			for (int j=i+1;j<numCols;j++)
			{
				U[j*numCols+i] /= U[i*numCols+i];
				for (int k=i+1; k<numCols; k++) U[j*numCols+k] -= U[j*numCols+i]*U[i*numCols+k];
			}
		}
	
		A = U;
	}

	T det(1);
	for (int i=0; i<numCols; i++) det *= A[i*numCols+i];
}


template<typename T>
T condNumber_Hadamard( T* A, int numCols )
{
	T res(1);
	for (int row=0; row<numCols; row++ ) res *= norm(A+row*numCols, numCols, 2);

	return determinant(A, numCols)/res;
}


template<typename T>
T condNumber(T* A, int numCols)
{


}



template<typename T>
void crossProduct3D( const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, T* res)
{	res[0] = y1*z2 -z1*y2;
	res[1] = z1*x2 -x1*z2;
	res[2] = x1*y2 -y1*x2;
}


template<typename T>
void planeUnitVector3D( const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, T* res)
{
	crossProduct3D( x1, y1, z1, x2, y2, z2, res);
	T norm = T(1)/sqrt( res[0]*res[0] +res[1]*res[1] +res[2]*res[2] );
	res[0] *= norm; res[1] *= norm; res[2] *= norm;
}

template <class T>
std::pair<size_t, T> absoluteMax(size_t size, T* values)
{
	T max(values[0]); size_t pos(0);
	for ( size_t i = 0; i < size; i++ ) {
		if ( values[i] < 0 ) { if ( -values[i] > max ) { max = -values[i]; pos = i;} }
		else { if( values[i] > max ) { max = values[i]; pos = i;} }
	}

	return std::pair<size_t, T>(pos, max);
}


template <class T>
std::pair<size_t, T> absoluteMin(size_t size, T* values)
{
	T min(values[0]); size_t pos(0);
	for ( size_t i = 0; i < size; i++ ) {
		if ( values[i] < 0 ) { if ( -values[i] < min ) { min = -values[i]; pos = i; } }
		else { if ( values[i] < min ) {min = values[i]; pos = i;} }
	}

	return std::pair<size_t, T>(pos, min);
}




// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// iSurfer.org makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
 

// a Point (or vector) is defined by its coordinates
typedef struct {int x, y, z;} Point;  // set z=0 for a 2D Point

// a Triangle is given by three points: Point V0, V1, V2 

// a Polygon is given by:
//       int n = number of vertex points
//       Point* V[] = an array of n+1 vertex points with V[n]=V[0]
 

// Note: for efficiency low-level short functions are declared to be inline.


// isLeft(): test if a point is Left|On|Right of an infinite 2D line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 to P1
//          =0 for P2 on the line
//          <0 for P2 right of the line
inline int isLeft( Point P0, Point P1, Point P2 )
{
    return ( (P1.x - P0.x) * (P2.y - P0.y)
           - (P2.x - P0.x) * (P1.y - P0.y) );
}
//===================================================================


// orientation2D_Triangle(): test the orientation of a 2D triangle
//  Input:  three vertex points V0, V1, V2
//  Return: >0 for counterclockwise 
//          =0 for none (degenerate)
//          <0 for clockwise
inline int orientation2D_Triangle( Point V0, Point V1, Point V2 )
{
    return isLeft(V0, V1, V2);
}
//===================================================================


// area2D_Triangle(): compute the area of a 2D triangle
//  Input:  three vertex points V0, V1, V2
//  Return: the (float) area of triangle T
inline float area2D_Triangle( Point V0, Point V1, Point V2 )
{
    return (float)isLeft(V0, V1, V2) / 2.0;
}
//===================================================================


// orientation2D_Polygon(): test the orientation of a simple 2D polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+1 vertex points with V[n]=V[0]
//  Return: >0 for counterclockwise 
//          =0 for none (degenerate)
//          <0 for clockwise
//  Note: this algorithm is faster than computing the signed area.
int orientation2D_Polygon( int n, Point* V )
{
    // first find rightmost lowest vertex of the polygon
    int rmin = 0;
    int xmin = V[0].x;
    int ymin = V[0].y;

    for (int i=1; i<n; i++) {
        if (V[i].y > ymin)
            continue;
        if (V[i].y == ymin) {   // just as low
            if (V[i].x < xmin)  // and to left
                continue;
        }
        rmin = i;      // a new rightmost lowest vertex
        xmin = V[i].x;
        ymin = V[i].y;
    }

    // test orientation at the rmin vertex
    // ccw <=> the edge leaving V[rmin] is left of the entering edge
    if (rmin == 0)
        return isLeft( V[n-1], V[0], V[1] );
    else
        return isLeft( V[rmin-1], V[rmin], V[rmin+1] );
}
 //===================================================================


// area2D_Polygon(): compute the area of a 2D polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+1 vertex points with V[n]=V[0]
//  Return: the (float) area of the polygon
float area2D_Polygon( int n, Point* V )
{
    float area = 0;
    int  i, j, k;   // indices

    if (n < 3) return 0;  // a degenerate polygon

    for (i=1, j=2, k=0; i<n; i++, j++, k++) {
        area += V[i].x * (V[j].y - V[k].y);
    }
    area += V[n].x * (V[1].y - V[n-1].y);  // wrap-around term
    return area / 2.0;
}
//===================================================================


// area3D_Polygon(): compute the area of a 3D planar polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+1 points in a 2D plane with V[n]=V[0]
//          Point N = a normal vector of the polygon's plane
//  Return: the (float) area of the polygon
float area3D_Polygon( int n, Point* V, Point N )
{
    float area = 0;
    float an, ax, ay, az; // abs value of normal and its coords
    int  coord;           // coord to ignore: 1=x, 2=y, 3=z
    int  i, j, k;         // loop indices

    if (n < 3) return 0;  // a degenerate polygon

    // select largest abs coordinate to ignore for projection
    ax = (N.x>0 ? N.x : -N.x);    // abs x-coord
    ay = (N.y>0 ? N.y : -N.y);    // abs y-coord
    az = (N.z>0 ? N.z : -N.z);    // abs z-coord

    coord = 3;                    // ignore z-coord
    if (ax > ay) {
        if (ax > az) coord = 1;   // ignore x-coord
    }
    else if (ay > az) coord = 2;  // ignore y-coord


    // compute area of the 2D projection
    switch (coord) {
      case 1:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].y * (V[j].z - V[k].z));
        break;
      case 2:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].z * (V[j].x - V[k].x));
        break;
      case 3:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].x * (V[j].y - V[k].y));
        break;
    }
    switch (coord) {    // wrap-around term
      case 1:
        area += (V[n].y * (V[1].z - V[n-1].z));
        break;
      case 2:
        area += (V[n].z * (V[1].x - V[n-1].x));
        break;
      case 3:
        area += (V[n].x * (V[1].y - V[n-1].y));
        break;
    }

    // scale to get area before projection
    an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
    switch (coord) {
      case 1:
        area *= (an / (2 * N.x));
        break;
      case 2:
        area *= (an / (2 * N.y));
        break;
      case 3:
        area *= (an / (2 * N.z));
    }
    return area;
}
//=================================================================== 


#include"myBLAS.ipp"
