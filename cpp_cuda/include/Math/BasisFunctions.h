#ifndef BASISFUNCTIONSHEADER
#define BASISFUNCTIONSHEADER

#include<functional>
#include"myBLAS.h"


using namespace std::placeholders;

template<class T>
class BasisFunctions
{

private:
	int m_degree;
	int m_dimensions;

public:
	BasisFunctions(int dimensions=1, int degree=2);

	int get_degree();
	int get_dimensions();
	void set_degree(int newDegree);
	void set_dimensions(int newDimensions);
	int  get_numOfMonomials()
	{ return numOfMonomials_upToDegree(m_dimensions, m_degree); }


	static void univariateMonomials(T* monomials, T evalPoint, int degree)
	{
		monomials[0] = 1;
		for (int index = 0; index < degree; index++) monomials[index+1] = monomials[index]*evalPoint;
	}

	static void bivariateMonomials(T* monomials, T evalPoint_x, T evalPoint_y, int degree)
	{  evaluateMonomials2D(degree, evalPoint_x, evalPoint_y, monomials); }


	static void trivariateMonomials(T* monomials, T evalPoint_x, T evalPoint_y, T evalPoint_z, int degree)
	{  evaluateMonomials3D(degree, evalPoint_x, evalPoint_y, evalPoint_z, monomials); }


	void operator () (T* monomials, T evalPoint);
	void operator () (T* monomials, T evalPoint_x, T evalPoint_y);
	void operator () (T* monomials, T evalPoint_x, T evalPoint_y, T evalPoint_z);

};


template class BasisFunctions<float>;
template class BasisFunctions<double>;

#endif
