#include"LeastSquares.h"
#include"myBLAS.h"
#include<cmath>


template<class T>
LeastSquares<T>::~LeastSquares()
{
	if (m_flags_createdData[0]) delete(this->m_scatteredData);
	if (m_flags_createdData[1]) delete(this->m_evaluationData);
	if (m_flags_createdData[2]) delete(this->m_basisFunctions);
	if (m_flags_createdData[3]) delete(this->m_linearSolver);
	if (m_flags_createdData[4]) delete(this->m_derivativesData);
}


template<class T>
LeastSquares<T>::LeastSquares(ScatteredData<T>* scatteredData, EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData)
{
	m_flags_createdData = std::vector<bool>(10, 0);

	if (scatteredData) 	m_scatteredData 	= scatteredData;
	else 			{m_scatteredData 	= new ScatteredData<T>; m_flags_createdData[0]=1;}

	if (evaluationData) 	m_evaluationData 	= evaluationData;
	else			{m_evaluationData 	= new EvaluationData<T>; m_flags_createdData[1]=1;}

	if (basisFunctions) 	m_basisFunctions 	= basisFunctions;
	else			{m_basisFunctions 	= new BasisFunctions<T>; m_flags_createdData[2]=1;}

	if (linearSolver) 	m_linearSolver 		= linearSolver; 
	else 			{m_linearSolver 	= new LinearSolver<T>; m_flags_createdData[3]=1;}

	if (derivativesData) 	m_derivativesData 	= derivativesData;
	else			{m_derivativesData 	= new EvaluationData<T>; m_flags_createdData[4]=1;}

}


template<class T>
LeastSquares<T>::LeastSquares(LeastSquares<T>* referenceLeastSquaresObject, ScatteredData<T>* scatteredData,  EvaluationData<T>* evaluationData, BasisFunctions<T>* basisFunctions, LinearSolver<T>* linearSolver, EvaluationData<T>* derivativesData )
{
	m_flags_createdData = std::vector<bool>(10, 0);
	copy_properties(referenceLeastSquaresObject);
	if (scatteredData) 	m_scatteredData 	= scatteredData;
	if (evaluationData) 	m_evaluationData 	= evaluationData;
	if (derivativesData) 	m_derivativesData 	= derivativesData;
	if (basisFunctions) 	m_basisFunctions 	= basisFunctions;
	if (linearSolver) 	m_linearSolver 		= linearSolver;
}


template<class T>
int LeastSquares<T>::set_scatteredData(ScatteredData<T>* newScatteredData)
{
	if (m_flags_createdData[0]) delete(this->m_scatteredData); m_flags_createdData[0]=0;
	m_scatteredData = newScatteredData; m_evaluationNeedFlag = 1;
	return 1;
}


template<class T>
int LeastSquares<T>::set_evaluationData(EvaluationData<T>* newEvaluationData)
{
	if (m_flags_createdData[1]) delete(this->m_evaluationData); m_flags_createdData[1]=0;
	m_evaluationData = newEvaluationData; m_evaluationNeedFlag = 1;
	return 1;
}

template<class T>
int LeastSquares<T>::set_derivativesData(EvaluationData<T>* newDerivativesData)
{
	if (m_flags_createdData[4]) delete(this->m_derivativesData); m_flags_createdData[4]=0;
	m_derivativesData = newDerivativesData;
	return 1;
}


template<class T>
int LeastSquares<T>::set_basisFunctions(BasisFunctions<T>* newBasisFunctions)
{
	if (m_flags_createdData[2]) delete(this->m_basisFunctions); m_flags_createdData[2]=0;
	m_basisFunctions = newBasisFunctions; m_evaluationNeedFlag = 1;
	return 1;
}


template<class T>
int LeastSquares<T>::set_linearSolver(LinearSolver<T>* newLinearSolver)
{
	if (m_flags_createdData[3]) delete(this->m_linearSolver); m_flags_createdData[3]=0;
	m_linearSolver = newLinearSolver; m_evaluationNeedFlag = 1;
	return 1;
}


template<class T>
std::vector<T>& LeastSquares<T>::get_coefficients() { return m_coefficients; }


template<class T>
void LeastSquares<T>::copy_properties(LeastSquares<T>* referenceLeastSquaresObject)
{
	m_scatteredData 	= referenceLeastSquaresObject->m_scatteredData;
	m_evaluationData 	= referenceLeastSquaresObject->m_evaluationData;
	m_derivativesData 	= referenceLeastSquaresObject->m_derivativesData;
	m_basisFunctions 	= referenceLeastSquaresObject->m_basisFunctions;
	m_linearSolver 		= referenceLeastSquaresObject->m_linearSolver;
}


template<class T>
void LeastSquares<T>::evaluatePoints()
{
	if (this->m_evaluationNeedFlag||this->m_updateNeedFlag)
	{
		cudaEvent_t start, stop; cudaEventCreate(&start); cudaEventCreate(&stop); cudaEventRecord(start);
		find_coefficients(); 
		cudaEventRecord(stop); cudaEventSynchronize(stop);
		cudaEventElapsedTime(&m_evalTime_coefficients, start, stop);

		this->m_evaluationNeedFlag = 0; 
		this->m_updateNeedFlag = 0;
	}
	cudaEvent_t start, stop; cudaEventCreate(&start); cudaEventCreate(&stop); cudaEventRecord(start);
	evaluationMethod();
	cudaEventRecord(stop); cudaEventSynchronize(stop);
	cudaEventElapsedTime(&m_evalTime_evalPoints, start, stop);
}


template<class T>
void LeastSquares<T>::evaluateDerivatives()
{
	cudaEvent_t start, stop; cudaEventCreate(&start); cudaEventCreate(&stop); cudaEventRecord(start);
	derivativesEvaluationMethod ();
	cudaEventRecord(stop); cudaEventSynchronize(stop);
	cudaEventElapsedTime(&m_evalTime_evalPoints, start, stop);
}




template <typename T>
void LeastSquares<T>::streamAllProperties(std::ostream& out)
{
	out << "Least squares fitting type: " << this->getLSType() << std::endl;
	out << "Number of spatial dimensions: " << this->get_scatteredData()->get_numOfDomainDimensions() << ", Number of fields:  " << this->get_scatteredData()->get_numOfFieldsPerPoint() << ", Number of data points: " << this->get_scatteredData()->get_numOfNodes() << std::endl;
	out << "Number of evaluation points: " << this->get_evaluationData()->get_numOfNodes() << std::endl;
	out << "Polynomial degree: " << this->get_basisFunctions()->get_degree();
	out << ", Polynomial terms: " << this->get_basisFunctions()->get_numOfMonomials() << std::endl;
	out << "Linear solver info: \n"; this->get_linearSolver()->write_solverInfo(out);
	std::cout << "Execution time of coefficients evaluation: " << m_evalTime_coefficients << "ms" << std::endl;
	std::cout << "Execution time to evaluate points: " << m_evalTime_evalPoints << "ms" << std::endl;
	std::cout << "Total execution time: " << m_evalTime_coefficients +m_evalTime_evalPoints << "ms" << std::endl;
}

std::istream& operator >> (std::istream& stream, computingMode& var )
{
	int temp; stream >> temp; 
	var = static_cast<computingMode> ( temp );
	return stream;
}


std::string get_computingMode_name( computingMode& mode )
{
	std::string compMode;
	switch (mode) {
		case CPU: { compMode="CPU"; break;  }
		case CPU_MT: { compMode="CPU multithreading"; break;  }
		case CUDA: { compMode="CUDA"; break;  }
		case CPU_CUDA: { compMode="CPU+CUDA"; break;  }
		case MultiGPU_CUDA: { compMode="CUDA MultiGPU"; break;  }
		case CPU_MultiGPU_CUDA: { compMode="CPU+CUDA MultiGPU"; break;  }
		case OPENCL: { compMode="OPENCL"; break;  }
	}
	return compMode;
}

std::ostream& operator << ( std::ostream& out, computingMode& mode )
{
	out << get_computingMode_name( mode );
	return out;
}

