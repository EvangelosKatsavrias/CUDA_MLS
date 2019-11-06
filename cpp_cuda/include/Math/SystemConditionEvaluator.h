#ifndef SYSTEMCONDITIONEVALUATORHEADER
#define SYSTEMCONDITIONEVALUATORHEADER 

template<class T> 
class SystemConditionEvaluator
{

	T*	m_data;
	int	m_numOfColumns;
	T	m_determinant_preestimate;
	T	m_HadamardCondNumber;
	T	m_spectralCondNumber;

public:
	SystemConditionEvaluator(T* data, int numOfColumns): m_data(data), m_numOfColumns(numOfColumns) {  }
	~SystemConditionEvaluator() {}


	static T preestimate_spectralCondNumber(T* data, int numOfColumns) { 
		T maxVal(0);
		for (int i=0; i<numOfColumns*numOfColumns; i+numOfColumns)
			if ( maxVal < data[i] ) maxVal = data[i];

		
	}

};


template class SystemConditionEvaluator<float>;
template class SystemConditionEvaluator<double>;

#endif
