#ifndef MYBLASHEADER
#define MYBLASHEADER

#include<utility>
#include<cuda_runtime.h>
#include<iostream>

template<class T>
void initializeVector(T* start, T* end, T initialValue);

template<class T>
void multiplicationOfVectorbyScalar(T* start, T* end, T* result, T scalar);

template<class T>
void tensorProductOfVectors(T* start, T* end, T* start_2, T* end_2, T* result);

template<class T>
void sumOfVectors(T* start, T* end, T* start_2, T* result);

template<class T>
T dotProductOfVectors(T* start, T* end, T* start_2);

template<class T>
void linspace(T start, T end, int numOfElements, T* output);

//template<class T>
//T pow_i( T base, int exp );

template<class T>
__host__ __device__
T pow_i( T base, int exp ) {
	T result = 1;
	while (1) {
		if (exp & 1) result *= base; exp >>= 1;
		if (!exp) break; base *= base; }
	return result;
}


int highestPowerof2(unsigned int n);

int factorial(int n);

int choose_k_outOf_n(int k, int n);

int numOfMonomials_ofDegree(int dim, int deg );

int numOfMonomials_upToDegree(int dim, int deg);

template<typename T>
void evaluateMonomials2D(int deg, T x, T y, T* eval);

template<typename T>
void evaluateMonomials3D(int deg, T x, T y, T z, T* eval);


template<typename T>
T distance2D(T x1, T x2, T y1, T y2);


template<typename T>
T distance3D(T x1, T x2, T y1, T y2, T z1, T z2);


template<typename T>
T norm(T* P1, T* P2, int dim, int normType );

template<typename T>
T determinant(T* A, int numCols, bool luFlag=0);

template<typename T>
T condNumber_Hadamard( T* A, int numCols );

template<typename T>
T condNumber(T* A, int numCols);

template<typename T>
void crossProduct3D( const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& y3, T* res);

template<typename T>
void planeUnitVector3D( const T& x1, const T& y1, const T& z1, const T& x2, const T& y2, const T& z2, T* res);



template<typename T, size_t numOfVertices = 3>
class PolygonArea
{

public:
	T operator() (T** x, T** y)
	{
		T area(0);
		for (int vertexIndex = 0; vertexIndex < numOfVertices-1; vertexIndex++ )
		{ area += det(*x[vertexIndex], *x[vertexIndex+1], *y[vertexIndex], *y[vertexIndex+1]); }
		  area += det(*x[numOfVertices-1], *x[0], *y[numOfVertices-1], *y[0]);
		return 0.5*area;
	}

T static det(T x1, T x2, T y1, T y2) { return x1*y2 -x2*y1; }

};


template <class T>
std::pair<size_t, T> absoluteMax(size_t size, T* values);

template <class T>
std::pair<size_t, T> absoluteMin(size_t size, T* values);



template<typename T, size_t numOfVertices = 3>
class PolygonArea3D
{

public:
	T operator() (T** x, T** y, T** z, T* planeUnitVector = 0)
	{
		T planeVectorTemp[3];
		if ( !planeUnitVector ) {
			T x1, x2, y1, y2, z1, z2;
			planeUnitVector = planeVectorTemp;
			x1 = *x[ 1 ] -*x[ 0 ];
			x2 = *x[ 2 ] -*x[ 0 ];
			y1 = *y[ 1 ] -*y[ 0 ];
			y2 = *y[ 2 ] -*y[ 0 ];
			z1 = *z[ 1 ] -*z[ 0 ];
			z2 = *z[ 2 ] -*z[ 0 ];
			planeUnitVector3D ( x1, y1, z1, x2, y2, z2, planeUnitVector );
		}

		T area;
		std::pair<size_t, T> maxComp = absoluteMax( 3, planeUnitVector );
		switch (maxComp.first) {
			case 0:	
				area = m_area2D(y, z); area = (area*planeUnitVector[0] > 0) ? area : -area; 
				return area/planeUnitVector[0];
			case 1: 
				area = m_area2D(x, z); area = (area*planeUnitVector[1] > 0) ? area : -area; 
				return area/planeUnitVector[1];
			case 2: 
				area = m_area2D(x, y); area = (area*planeUnitVector[2] > 0) ? area : -area; 
				return area/planeUnitVector[2];
		}
	}

protected:
	PolygonArea<T, numOfVertices> m_area2D;
};



template<typename T, size_t numOfVerticesPerFace = 4, size_t numOfFaces = 4>
class PolyhedronVolume
{

public:
	T operator() (T** x, T** y, T** z)
	{
		T volume(0); T faceDiv; T unitVector[3];
		size_t faceVertexShift;
		
		for (int faceIndex = 0; faceIndex < numOfFaces-1; faceIndex++ )
		{
			faceVertexShift = faceIndex*numOfVerticesPerFace; 
			T x1, x2, y1, y2, z1, z2;
			x1 = *x[ faceVertexShift ];
			x2 = *x[ faceVertexShift +1 ];
			y1 = *y[ faceVertexShift ];
			y2 = *y[ faceVertexShift +1 ];
			z1 = *z[ faceVertexShift ];
			z2 = *z[ faceVertexShift +1 ];

			planeUnitVector3D ( x1, y1, z1, x2, y2, z2, unitVector );
			volume += (unitVector[0]*x1 +unitVector[1]*y1 +unitVector[2]*z1) *m_faceAreaFunctor( x +faceVertexShift, y +faceVertexShift, z +faceVertexShift, unitVector );
		}

		return volume/T(3);
	}

protected:
	PolygonArea3D<T, numOfVerticesPerFace> m_faceAreaFunctor;

};



#endif
