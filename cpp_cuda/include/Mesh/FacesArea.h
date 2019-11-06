#ifndef FACESAREAHEADERFILE
#define FACESAREAHEADERFILE

#include<stdexcept>
#include<exception>
#include"ContainerBase.h"
#include"Topology.h"
#include"myBLAS.h"
#include<math.h>
#include<algorithm>

template< typename nodalComponentType, typename connectivityType, size_t numOfConnectivityNodes = 3 >
class ElementsArea: public ContainerBase<nodalComponentType, std::vector<nodalComponentType> >
{

typedef ContainerBase<nodalComponentType, std::vector<nodalComponentType> > containerBase_t;

protected:
	Nodes<nodalComponentType>*					m_nodeCoords;
	Connectivities<connectivityType, numOfConnectivityNodes>*	m_connectivities;
	PolygonArea<nodalComponentType, numOfConnectivityNodes>*	m_areaEstimator;
	std::vector<bool>						m_internalMemberFlags;


public:
	ElementsArea(Nodes<nodalComponentType>* nodes=0, Connectivities<connectivityType, numOfConnectivityNodes>* connectivities=0, PolygonArea<nodalComponentType, numOfConnectivityNodes>* estimator=0 );
	ElementsArea( ElementsArea<nodalComponentType, connectivityType, numOfConnectivityNodes>* referenceObject);
	ElementsArea( ElementsArea<nodalComponentType, connectivityType, numOfConnectivityNodes>& referenceObject);
	~ElementsArea();


	Nodes<nodalComponentType>*					get_nodeCoords() 	{ return m_nodeCoords; }
	Connectivities<connectivityType, numOfConnectivityNodes>*	get_connectivities()	{ return m_connectivities; }
	PolygonArea<nodalComponentType, numOfConnectivityNodes>*	get_areaEstimator()	{ return m_areaEstimator; }
	int								set_nodeCoords(Nodes<nodalComponentType>* nodeCoords) 	{ m_nodeCoords=nodeCoords; }
	int								set_connectivities(Connectivities<connectivityType, numOfConnectivityNodes>* connectivities)	{ m_connectivities=connectivities; }
	int								set_areaEstimator(PolygonArea<nodalComponentType, numOfConnectivityNodes>* areaEstimator)	{ m_areaEstimator=areaEstimator; }


	int find_areas();

};


template<typename T, typename T1, size_t T2>
ElementsArea<T, T1, T2>::ElementsArea(Nodes<T>* nodes, Connectivities<T1, T2>* connectivities, PolygonArea<T, T2>* estimator ):
m_nodeCoords(nodes), m_connectivities(connectivities), m_areaEstimator(estimator), m_internalMemberFlags(std::vector<bool>(3,0))
{
	if (!nodes) 		{ m_nodeCoords 		= new Nodes<T>; 		m_internalMemberFlags[0]=1; }
	if (!connectivities) 	{ m_connectivities 	= new Connectivities<T1, T2>; 	m_internalMemberFlags[1]=1; }
	if (!estimator) 	{ m_areaEstimator 	= new PolygonArea<T, T2>; 	m_internalMemberFlags[2]=1; }
}


template<typename nodalComponentType, typename connectivityType, size_t numOfConnectivityNodes>
ElementsArea<nodalComponentType, connectivityType, numOfConnectivityNodes>::ElementsArea( ElementsArea<nodalComponentType, connectivityType, numOfConnectivityNodes>* referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType>>(referenceObject), m_internalMemberFlags(std::vector<bool>(3,0))
{
	m_nodeCoords 	= referenceObject->get_nodeCoords();
	m_connectivities= referenceObject->get_connectivities();
	m_areaEstimator = referenceObject->get_areaEstimator();
}


template<typename nodalComponentType, typename connectivityType, size_t numOfConnectivityNodes>
ElementsArea<nodalComponentType, connectivityType, numOfConnectivityNodes>::ElementsArea( ElementsArea<nodalComponentType, connectivityType, numOfConnectivityNodes>& referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType> >(referenceObject), m_internalMemberFlags(std::vector<bool>(3,1))
{
	m_nodeCoords 		= new Nodes<nodalComponentType>( *(referenceObject.get_nodeCoords() ) );
	m_connectivities 	= new Connectivities<connectivityType, numOfConnectivityNodes>( *(referenceObject.get_connectivities() ) );
	m_areaEstimator 	= new PolygonArea<nodalComponentType, numOfConnectivityNodes>;
}


template<typename T, typename T1, size_t T2>
ElementsArea<T, T1, T2>::~ElementsArea()
{
	if (m_internalMemberFlags[0]) delete m_nodeCoords;
	if (m_internalMemberFlags[1]) delete m_connectivities;
	if (m_internalMemberFlags[2]) delete m_areaEstimator;
}


template<typename T, typename T1, size_t T2>
int ElementsArea<T, T1, T2>::find_areas()
{
	T *x[T2], *y[T2];

	this->resize(m_connectivities->get_numOfConnectivities());

	for ( size_t elementIndex = 0; elementIndex < m_connectivities->get_numOfConnectivities(); elementIndex++ ) {

		T1* elementVertices = (*m_connectivities)[elementIndex];

		for (int vertexIndex = 0; vertexIndex < T2; vertexIndex++ ) {
			x[vertexIndex] = m_nodeCoords[0][0]+elementVertices[vertexIndex];
			y[vertexIndex] = m_nodeCoords[0][1]+elementVertices[vertexIndex];

		}

		this->m_data[elementIndex] = (*m_areaEstimator)(x, y);

	}

}


template< class nodalComponentType, class connectivityType, size_t numOfVerticesPerFace = 3>
class FacesArea: public ContainerBase<nodalComponentType, std::vector<nodalComponentType> >
{

protected:
	Nodes<nodalComponentType>*					m_nodeCoords;
	Connectivities<connectivityType, numOfVerticesPerFace>*		m_connectivities;
	PolygonArea3D<nodalComponentType, numOfVerticesPerFace>*	m_areaEstimator;
	std::vector<bool>						m_allocatedMemberFlags;


public:
	FacesArea ( Nodes<nodalComponentType>* nodes=0, Connectivities<connectivityType, numOfVerticesPerFace>* connectivities=0, PolygonArea3D<nodalComponentType, numOfVerticesPerFace>* estimator=0 );
	FacesArea ( FacesArea <nodalComponentType, connectivityType, numOfVerticesPerFace>* referenceObject);
	FacesArea ( FacesArea <nodalComponentType, connectivityType, numOfVerticesPerFace>& referenceObject);
	~FacesArea ( );


	Nodes <nodalComponentType>*					get_nodeCoords() 	{ return m_nodeCoords; }
	Connectivities <connectivityType, numOfVerticesPerFace>*	get_connectivities()	{ return m_connectivities; }
	PolygonArea3D <nodalComponentType, numOfVerticesPerFace>*	get_areaEstimator()	{ return m_areaEstimator; }


	int set_nodeCoords (Nodes<nodalComponentType>* nodeCoords) 						{ m_nodeCoords 		= nodeCoords; }
	int set_vertexConnectivities (Connectivities<connectivityType, numOfVerticesPerFace>* connectivities)	{ m_connectivities 	= connectivities; }
	int set_areaEstimator (PolygonArea3D<nodalComponentType, numOfVerticesPerFace>* areaEstimator)		{ m_areaEstimator 	= areaEstimator; }


	int virtual find_areas();

};


template<typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace>
FacesArea<nodalComponentType, connectivityType, numOfVerticesPerFace>::FacesArea(Nodes<nodalComponentType>* nodes, Connectivities<connectivityType, numOfVerticesPerFace>* connectivities, PolygonArea3D<nodalComponentType, numOfVerticesPerFace>* estimator ):
m_nodeCoords(nodes), m_connectivities(connectivities), m_areaEstimator(estimator), m_allocatedMemberFlags(std::vector<bool>(3,0))
{
	if (!nodes) 		{ m_nodeCoords 		= new Nodes<nodalComponentType>; 		m_allocatedMemberFlags[0]=1; }
	if (!connectivities) 	{ m_connectivities 	= new Connectivities<connectivityType, numOfVerticesPerFace>; 	m_allocatedMemberFlags[1]=1; }
	if (!estimator) 	{ m_areaEstimator 	= new PolygonArea3D<nodalComponentType, numOfVerticesPerFace>; 	m_allocatedMemberFlags[2]=1; }
}


template<typename nodalComponentType, typename connectivityType, size_t connectivityNodes>
FacesArea<nodalComponentType, connectivityType, connectivityNodes>::FacesArea( FacesArea<nodalComponentType, connectivityType, connectivityNodes>* referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType>>(referenceObject), m_allocatedMemberFlags(std::vector<bool>(3,0))
{
	m_nodeCoords 	= referenceObject->get_nodeCoords();
	m_connectivities= referenceObject->get_connectivities();
	m_areaEstimator = referenceObject->get_areaEstimator();
}


template<typename nodalComponentType, typename connectivityType, size_t connectivityNodes>
FacesArea<nodalComponentType, connectivityType, connectivityNodes>::FacesArea( FacesArea<nodalComponentType, connectivityType, connectivityNodes>& referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType> >(referenceObject), m_allocatedMemberFlags(std::vector<bool>(3,1))
{
	m_nodeCoords 		= new Nodes<nodalComponentType>( *(referenceObject.get_nodeCoords() ) );
	m_connectivities 	= new Connectivities<connectivityType, connectivityNodes>( *(referenceObject.get_connectivities() ) );
	m_areaEstimator 	= new PolygonArea3D<nodalComponentType, connectivityNodes>;
}


template<typename T, typename connectivityType, size_t numOfVerticesPerFace>
FacesArea<T, connectivityType, numOfVerticesPerFace>::~FacesArea()
{
	if (m_allocatedMemberFlags[0]) delete m_nodeCoords;
	if (m_allocatedMemberFlags[1]) delete m_connectivities;
	if (m_allocatedMemberFlags[2]) delete m_areaEstimator;
}


template<typename nodalComponentType, typename connectivityType, size_t numOfVerticesPerFace>
int FacesArea<nodalComponentType, connectivityType, numOfVerticesPerFace>::find_areas()
{
	nodalComponentType *x[numOfVerticesPerFace], *y[numOfVerticesPerFace], *z[numOfVerticesPerFace];

	this->resize(m_connectivities->get_numOfConnectivities());

	for ( size_t elementIndex = 0; elementIndex < m_connectivities->get_numOfConnectivities(); elementIndex++ ) {

		connectivityType* elementVertices = m_connectivities->at(elementIndex);

		for (int vertexIndex = 0; vertexIndex < numOfVerticesPerFace; vertexIndex++ ) {
			x[vertexIndex] = m_nodeCoords[0][0]+elementVertices[vertexIndex];
			y[vertexIndex] = m_nodeCoords[0][1]+elementVertices[vertexIndex];
			z[vertexIndex] = m_nodeCoords[0][2]+elementVertices[vertexIndex];
		}

		this->m_data[elementIndex] = (*m_areaEstimator)(x, y, z);
	}
}




template< typename nodalComponentType, typename connectivityType, size_t numOfConnectivityNodes = 3 >
class ElementsQuality: public ContainerBase<nodalComponentType, std::vector<nodalComponentType> >
{

public:
typedef ContainerBase<nodalComponentType, std::vector<nodalComponentType> > containerBase_t;
typedef ElementsArea<nodalComponentType, size_t, numOfConnectivityNodes> 	elementsArea_t;

protected:
	Nodes<nodalComponentType>*					m_nodeCoords;
	Connectivities<connectivityType, numOfConnectivityNodes>*	m_connectivities;
//	PolygonArea<nodalComponentType, numOfConnectivityNodes>*	m_areaEstimator;
	std::vector<bool>						m_internalMemberFlags;
	elementsArea_t*							m_elementsArea;

public:
	ElementsQuality(Nodes<nodalComponentType>* nodes=0, Connectivities<connectivityType, numOfConnectivityNodes>* connectivities=0, elementsArea_t* elementsArea=0 );
	ElementsQuality( ElementsQuality<nodalComponentType, connectivityType, numOfConnectivityNodes>* referenceObject);
	ElementsQuality( ElementsQuality<nodalComponentType, connectivityType, numOfConnectivityNodes>& referenceObject);
	~ElementsQuality();


	elementsArea_t*							get_elementsArea()	{ return m_elementsArea; }
	Nodes<nodalComponentType>*					get_nodeCoords() 	{ return m_nodeCoords; }
	Connectivities<connectivityType, numOfConnectivityNodes>*	get_connectivities()	{ return m_connectivities; }
//	PolygonArea<nodalComponentType, numOfConnectivityNodes>*	get_areaEstimator()	{ return m_areaEstimator; }
	int								set_nodeCoords(Nodes<nodalComponentType>* nodeCoords) 	{ m_nodeCoords=nodeCoords; }
	int								set_connectivities(Connectivities<connectivityType, numOfConnectivityNodes>* connectivities)	{ m_connectivities=connectivities; }
//	int								set_areaEstimator(PolygonArea<nodalComponentType, numOfConnectivityNodes>* areaEstimator)	{ m_areaEstimator=areaEstimator; }
	int								set_elementsArea(elementsArea_t* elementsArea) { m_elementsArea = elementsArea; }

	int find_qualities();

};


template<typename T, typename T1, size_t T2>
ElementsQuality<T, T1, T2>::ElementsQuality(Nodes<T>* nodes, Connectivities<T1, T2>* connectivities, elementsArea_t* elementsArea):
m_nodeCoords(nodes), m_connectivities(connectivities), m_elementsArea(elementsArea), m_internalMemberFlags(std::vector<bool>(10,0))
{
	if (!nodes) 		{ m_nodeCoords 		= new Nodes<T>; 		m_internalMemberFlags[0]=1; }
	if (!connectivities) 	{ m_connectivities 	= new Connectivities<T1, T2>; 	m_internalMemberFlags[1]=1; }
//	if (!estimator) 	{ m_areaEstimator 	= new PolygonArea<T, T2>; 	m_internalMemberFlags[2]=1; }
	if (!elementsArea) 	{ m_elementsArea 	= new ElementsArea<T, T1, T2>; 	m_internalMemberFlags[3]=1; }
}


template<typename nodalComponentType, typename connectivityType, size_t numOfConnectivityNodes>
ElementsQuality<nodalComponentType, connectivityType, numOfConnectivityNodes>::ElementsQuality( ElementsQuality<nodalComponentType, connectivityType, numOfConnectivityNodes>* referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType>>(referenceObject), m_internalMemberFlags(std::vector<bool>(10,0))
{
	m_nodeCoords 	= referenceObject->get_nodeCoords();
	m_connectivities= referenceObject->get_connectivities();
//	m_areaEstimator = referenceObject->get_areaEstimator();
	m_elementsArea = referenceObject->get_elementsArea();
}


template<typename nodalComponentType, typename connectivityType, size_t numOfConnectivityNodes>
ElementsQuality<nodalComponentType, connectivityType, numOfConnectivityNodes>::ElementsQuality( ElementsQuality<nodalComponentType, connectivityType, numOfConnectivityNodes>& referenceObject):
ContainerBase<nodalComponentType, std::vector<nodalComponentType> >(referenceObject), m_internalMemberFlags(std::vector<bool>(10,1))
{
	m_nodeCoords 		= new Nodes<nodalComponentType>( *(referenceObject.get_nodeCoords() ) );
	m_connectivities 	= new Connectivities<connectivityType, numOfConnectivityNodes>( *(referenceObject.get_connectivities() ) );
//	m_areaEstimator 	= new PolygonArea<nodalComponentType, numOfConnectivityNodes>;
	m_elementsArea 		= new ElementsArea<nodalComponentType, connectivityType, numOfConnectivityNodes>;
}


template<typename T, typename T1, size_t T2>
ElementsQuality<T, T1, T2>::~ElementsQuality()
{
	if (m_internalMemberFlags[0]) delete m_nodeCoords;
	if (m_internalMemberFlags[1]) delete m_connectivities;
//	if (m_internalMemberFlags[2]) delete m_areaEstimator;
	if (m_internalMemberFlags[3]) delete m_elementsArea;
}


template<typename T, typename T1, size_t T2>
int ElementsQuality<T, T1, T2>::find_qualities()
{
	T *x[T2], *y[T2];

	this->resize(m_connectivities->get_numOfConnectivities());

	for ( size_t elementIndex = 0; elementIndex < m_connectivities->get_numOfConnectivities(); elementIndex++ ) {

		T1* elementVertices = (*m_connectivities)[elementIndex];

		for (int vertexIndex = 0; vertexIndex < T2; vertexIndex++ ) {
			x[vertexIndex] = m_nodeCoords[0][0]+elementVertices[vertexIndex];
			y[vertexIndex] = m_nodeCoords[0][1]+elementVertices[vertexIndex];

		}

		T sideLengths[ T2 ];
		for (int vertexIndex = 0; vertexIndex < T2-1; vertexIndex++ ) {
			sideLengths[ vertexIndex ] = distance2D(*x[vertexIndex], *x[vertexIndex+1], 
								*y[vertexIndex], *y[vertexIndex+1]);
		}
		sideLengths[ T2-1 ] = distance2D(*x[T2-1], *x[0], *y[T2-1], *y[0]);

		const T pi = 3.141592654;

		T angles[T2];
		angles[ 0 ] = ( *x[T2-1] -*x[0] )*( *x[1] -*x[0] )
				+( *y[T2-1] -*y[0] )*( *y[1] -*y[0] );
		angles[ 0 ] /= sideLengths[ T2-1 ]*sideLengths[ 0 ];
		if ( angles[ 0 ] >= T(1) ) angles[ 0 ] = 0;
		else if ( angles[ 0 ] <= -T(1) ) angles[ 0 ] = pi;
		else angles[ 0 ] = acos( angles[0] );
		for (int vertexIndex = 1; vertexIndex < T2-1; vertexIndex++ ) {
			angles[ vertexIndex ] = ( *x[vertexIndex-1] -*x[vertexIndex] )*( *x[vertexIndex+1] -*x[vertexIndex] ) 
						+( *y[vertexIndex-1] -*y[vertexIndex] )*( *y[vertexIndex+1] -*y[vertexIndex] );
			angles[ vertexIndex ] /= sideLengths[ vertexIndex-1 ]*sideLengths[ vertexIndex ];
			if ( angles[ vertexIndex ] >= T(1) ) angles[ vertexIndex ] = 0;
			else if ( angles[ vertexIndex ] <= -T(1) ) angles[ vertexIndex ] = pi;
			else angles[ vertexIndex ] = acos( angles[vertexIndex] );
		}
		angles[ T2-1 ] = ( *x[T2-2] -*x[T2-1] )*( *x[0] -*x[T2-1] ) 
					+( *y[T2-2] -*y[T2-1] )*( *y[0] -*y[T2-1] );
		angles[ T2-1 ] /= sideLengths[ T2-2 ]*sideLengths[ T2-1 ];
		if ( angles[ T2-1 ] >= T(1) ) angles[ T2-1 ] = 0;
		else if ( angles[ T2-1 ] <= -T(1) ) angles[ T2-1 ] = pi;
		else angles[ T2-1 ] = acos( angles[T2-1] );

		auto minmaxAngles = std::minmax_element( angles, angles+T2 );

		T m1, m2; 
		switch ( T2 ) {
			case 3: { 
					const T equiAngle = 1.047197551; 
					m1 = (*minmaxAngles.second-equiAngle)/( pi-equiAngle );
					m2 = (equiAngle -*minmaxAngles.first)/( equiAngle );
				}
		}
		
		this->m_data[elementIndex] = std::max( m1, m2 );
	}
}

/*
template<typename T, typename T1, size_t T2>
int ElementsQuality<T, T1, T2>::find_qualities()
{
	T *x[T2], *y[T2];

	this->resize(m_connectivities->get_numOfConnectivities());

	for ( size_t elementIndex = 0; elementIndex < m_connectivities->get_numOfConnectivities(); elementIndex++ ) {

		T1* elementVertices = (*m_connectivities)[elementIndex];

		for (int vertexIndex = 0; vertexIndex < T2; vertexIndex++ ) {
			x[vertexIndex] = m_nodeCoords[0][0]+elementVertices[vertexIndex];
			y[vertexIndex] = m_nodeCoords[0][1]+elementVertices[vertexIndex];

		}

		T sides[ T2 ];
		for (int vertexIndex = 0; vertexIndex < T2-1; vertexIndex++ ) {
			sides[ vertexIndex ] = distance2D(*x[vertexIndex], 
							*x[vertexIndex+1], 
							*y[vertexIndex], 
							*y[vertexIndex+1]);
		}
		sides[ T2-1 ] = distance2D(*x[T2-1], *x[0], *y[T2-1], *y[0]);

		T sumOfSides(0); T prodOfSides(1);
		for (int vertexIndex = 0; vertexIndex < T2; vertexIndex++ ) {
			sumOfSides += sides[vertexIndex]; prodOfSides *= sides[ vertexIndex ];
		}

		T area = m_elementsArea[0][elementIndex];
		T innerCircleRadius = 2*area/sumOfSides; 
		T outerCircleRadius = prodOfSides/(4*area);

		this->m_data[elementIndex] = outerCircleRadius/innerCircleRadius;
	}
}

*/



#endif
