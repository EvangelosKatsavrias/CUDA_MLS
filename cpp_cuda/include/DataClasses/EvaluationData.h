#ifndef EVALUATIONDATAHEADER
#define EVALUATIONDATAHEADER

#include"ScatteredData.h"

template<class T=float>
class EvaluationData : public ScatteredData<T>
{
public:
	using ScatteredData<T>::ScatteredData;
};

template class EvaluationData<float>;
template class EvaluationData<double>;


#endif
