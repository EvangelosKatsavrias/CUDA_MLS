template<typename T>
void LU_Crout_decomposer(T* A, T* L, T* U, int n);

template<typename T>
void LU_Doolittle_decomposer(T* A, T* L, T* U, int n);

template<typename T>
void LU_decomposer_inPlace(T* A, int n);

template<typename T>
void LDU_decomposer(T* A, T* L, T* U, T* D, int n);

template<typename T>
void LDL_decomposer(T* A, T* L, T* D, int n);

template<typename T>
void Cholesky_decomposer(T* A, T* L, int n);

template<typename T>
void Cholesky_decomposer_banded(T* A, T* L, int n, int n_b);

template<typename T>
void Cholesky_decomposer_skylineInPlace(T* A, int n, int* n_max);
