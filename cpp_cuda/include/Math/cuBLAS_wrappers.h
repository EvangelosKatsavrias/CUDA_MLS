__host__ __device__
inline cublasStatus_t cublasGgemv(cublasHandle_t handle, cublasOperation_t trans, int m, int n, const float *alpha, const float *A, int lda, const float *x, int incx, const float *beta, float *y, int incy) {
	return cublasSgemv(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy); }
__host__ __device__
inline cublasStatus_t cublasGgemv(cublasHandle_t handle, cublasOperation_t trans, int m, int n, const double *alpha, const double *A, int lda, const double *x, int incx, const double *beta, double *y, int incy) { 
	return cublasDgemv(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy); }

__host__ __device__
inline cublasStatus_t cublasGgemm(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, int m, int n, int k, const float *alpha, const float *A, int lda, const float *B, int ldb, const float *beta, float *C, int ldc) {
	return cublasSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); }
__host__ __device__
inline cublasStatus_t cublasGgemm(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb, int m, int n, int k, const double *alpha, const double *A, int lda, const double *B, int ldb, const double *beta, double *C, int ldc) { 
	return cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); }


