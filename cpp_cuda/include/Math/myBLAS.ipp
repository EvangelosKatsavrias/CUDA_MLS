
template void initializeVector(float* start, float* end, float initialValue);
template void multiplicationOfVectorbyScalar(float* start, float* end, float* result, float scalar);
template void tensorProductOfVectors(float* start, float* end, float* start_2, float* end_2, float* result);
template void sumOfVectors(float* start, float* end, float* start_2, float* result);
template float dotProductOfVectors(float* start, float* end, float* start_2);
template void linspace(float start, float end, int numOfElements, float* output);

template float pow_i( float base, int exp );
template double pow_i( double base, int exp );

template void initializeVector(double* start, double* end, double initialValue);
template void multiplicationOfVectorbyScalar(double* start, double* end, double* result, double scalar);
template void tensorProductOfVectors(double* start, double* end, double* start_2, double* end_2, double* result);
template void sumOfVectors(double* start, double* end, double* start_2, double* result);
template double dotProductOfVectors(double* start, double* end, double* start_2);
template void linspace(double start, double end, int numOfElements, double* output);

template void evaluateMonomials2D(int deg, float x, float y, float* eval);
template void evaluateMonomials2D(int deg, double x, double y, double* eval);

template void evaluateMonomials3D(int deg, float x, float y, float z, float* eval);
template void evaluateMonomials3D(int deg, double x, double y, double z, double* eval);

template float distance2D(float x1, float x2, float y1, float y2);
template double distance2D(double x1, double x2, double y1, double y2);

template float distance3D(float x1, float x2, float y1, float y2, float z1, float z2);
template double distance3D(double x1, double x2, double y1, double y2, double z1, double z2);


template float norm(float* x, int dim, int normType );
template double norm(double* x, int dim, int normType );

template float determinant(float* A, int numCols, bool luFlag=0);
template double determinant(double* A, int numCols, bool luFlag=0);

template float condNumber_Hadamard(float* A, int numCols );
template double condNumber_Hadamard(double* A, int numCols );

template float condNumber(float* A, int numCols);
template double condNumber(double* A, int numCols);

template void crossProduct3D( const float& x1, const float& y1, const float& z1, const float& x2, const float& y2, const float& z2, float* res);
template void crossProduct3D( const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2, double* res);

template void planeUnitVector3D( const float& x1, const float& y1, const float& z1, const float& x2, const float& y2, const float& z2, float* res);
template void planeUnitVector3D( const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2, double* res);

template std::pair<size_t, float> absoluteMax(size_t size, float* values);
template std::pair<size_t, double> absoluteMax(size_t size, double* values);

template std::pair<size_t, float> absoluteMin(size_t size, float* values);
template std::pair<size_t, double> absoluteMin(size_t size, double* values);


