#include"BasisFunctions.h"

template<class T>
BasisFunctions<T>::BasisFunctions(int dimensions, int degree)
{ set_dimensions(dimensions); set_degree(degree); }


template<class T>
int BasisFunctions<T>::get_degree() {return m_degree;}
template<class T>
int BasisFunctions<T>::get_dimensions() {return m_dimensions;}
template<class T>
void BasisFunctions<T>::set_degree(int newDegree) {m_degree = newDegree;}
template<class T>
void BasisFunctions<T>::set_dimensions(int newDimensions) {m_dimensions = newDimensions;}

template<class T>
void BasisFunctions<T>::operator () (T* monomials, T evalPoint)
{ univariateMonomials(monomials, evalPoint, m_degree); }

template<class T>
void BasisFunctions<T>::operator () (T* monomials, T evalPoint_x, T evalPoint_y)
{ bivariateMonomials(monomials, evalPoint_x, evalPoint_y, m_degree );  }

template<class T>
void BasisFunctions<T>::operator () (T* monomials, T evalPoint_x, T evalPoint_y, T evalPoint_z)
{ trivariateMonomials(monomials, evalPoint_x, evalPoint_y, evalPoint_z, m_degree );  }
