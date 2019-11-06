#ifndef LEASTSQUARESEXCEPTIONS
#define LEASTSQUARESEXCEPTIONS

#include<exception>
#include<iostream>
#include"LeastSquares.h"

class leastSquaresException: public std::exception
{
	protected:
		LeastSquares* m_referenceLeastSquares;

	public:
		leastSquaresException(LeastSquares* refLeastSquares=0);
};

class leastSquaresException_numOfValues: public leastSquaresException
{
	public:
		leastSquaresException_numOfValues();
};

#endif
