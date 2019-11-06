#include<string>
#include"LeastSquares.h"
#include"fileFunctions.h"

void weightedLSTest()
{
	using chosenType = float;


//========== Load data ===================================================
	size_t numOfDataPoints;
	std::string dataSetPath 	= "dataset_1";
	std::string sampleDataPath 	= dataSetPath + "/data";
	std::string evalPointsPath 	= dataSetPath + "/data_evalPoints";
	std::string stationPointsPath 	= dataSetPath + "/data_stationaryPoints2";


	getFromTextFile_numOfColumnsAndLines(sampleDataPath, 0, &numOfDataPoints);
	ScatteredData<chosenType> samplePoints(1, 1, numOfDataPoints);
	samplePoints.set_AllComponents(sampleDataPath);

	getFromTextFile_numOfColumnsAndLines(evalPointsPath, 0, &numOfDataPoints);
	EvaluationData<chosenType> wLS_EvalData(1, 1, numOfDataPoints);
	wLS_EvalData.set_AllDomainComponents(evalPointsPath);
	
	getFromTextFile_numOfColumnsAndLines(stationPointsPath, 0, &numOfDataPoints);
	StationaryPoints<chosenType> stationaryPoints(1, numOfDataPoints);
	stationaryPoints.set_components(stationPointsPath);


//========== Create LS objects ===========================================
	WeightingFunction<chosenType> weightFunction(Wendland, 0.2);
	BasisFunctions<chosenType> basisFunctions(1, 2);
	LinearSolver<chosenType> linearSolver(gaussEliminationWithBackSubst);

	WeightedLeastSquares<chosenType> 	wLS (&weightFunction, &samplePoints, &wLS_EvalData, &stationaryPoints, &basisFunctions, &linearSolver); wLS.set_computingMode(CUDA);
	WeightedLeastSquares_LCS<chosenType> 	wLS_LCS (&wLS);
	MovingLeastSquares<chosenType> 		movingLS (&wLS);
	MovingLeastSquares_LCS<chosenType> 	movingLS_LCS (&wLS);


//========== Execute and report ==========================================
	wLS.evaluatePoints();		std::cout << wLS;
	wLS_LCS.evaluatePoints();	std::cout << wLS_LCS;
	movingLS.evaluatePoints();	std::cout << movingLS;
	movingLS_LCS.evaluatePoints();	std::cout << movingLS_LCS;

}
