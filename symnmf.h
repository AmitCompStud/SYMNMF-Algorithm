#ifndef SYMNMF_H_ 
#define SYMNMF_H_
# define HANDLE_ERROR() do { \
    printf("An Error Has Occurred\n"); \
    exit(1); \
} while (0)

/**
 * Initializes a 2D matrix of doubles with the specified dimensions.
 * Allocates memory for a length x width matrix and stores the data pointer for later freeing.
 * 
 * @param length Number of rows in the matrix
 * @param width Number of columns in the matrix  
 * @param pPointer Pointer to store the address of the allocated data array
 * @return Pointer to the initialized matrix, or NULL if allocation failed
 */
double** initializeMatrix(int length, int width, double** pPointer);

/**
 * Computes the similarity matrix A from input data points.
 * A[i][j] = exp(-||x_i - x_j||^2 / 2) for i != j, A[i][i] = 0.
 * The matrix is symmetric, so only upper triangle is computed.
 * 
 * @param A Pointer to the similarity matrix to be filled (N x N)
 * @param N Number of data points
 * @param D Dimensionality of each data point
 * @param inputs Input data points (N x D matrix)
 */
void similarityMatrix(double**A, int N, int D, double** inputs);

/**
 * Computes the diagonal degree matrix D from similarity matrix A.
 * D[i][i] = sum of row i in A, D[i][j] = 0 for i != j.
 * 
 * @param D Pointer to the diagonal matrix to be filled (N x N)
 * @param A Similarity matrix (N x N)
 * @param N Size of the matrices
 */
void diagonalMatrix(double** D, double** A, int N);

/**
 * Computes the inverse square root of diagonal elements.
 * DTag[i][i] = 1/sqrt(D[i][i]), DTag[i][j] = 0 for i != j.
 * 
 * @param DTag Pointer to the result matrix (N x N)
 * @param matrix Input diagonal matrix D (N x N)
 * @param N Size of the matrices
 */
void raiseNegativeHalfDiag(double** DTag, double** matrix, int N);

/**
 * Computes the normalized similarity matrix W.
 * W[i][j] = DTag[i][i] * A[i][j] * DTag[j][j] for i != j, W[i][i] = 0.
 * The matrix is symmetric.
 * 
 * @param W Pointer to the normalized matrix to be filled (N x N)
 * @param A Similarity matrix (N x N)
 * @param DTag Diagonal matrix with inverse square roots (N x N)
 * @param N Size of the matrices
 */
void normalizedSimilarityMatrix(double** W, double** A, double** DTag, int N);

/**
 * Performs the SymNMF algorithm to factorize matrix W into H*H^T.
 * Iteratively updates H until convergence or max iterations reached.
 * 
 * @param W Normalized similarity matrix (N x N)
 * @param H Initial H matrix to be updated (N x k)
 * @param N Number of data points
 * @param k Number of clusters
 * @return 0 on success, 1 on memory allocation failure
 */
int symnmfAlgorithm(double** W, double** H, int N, int k);

/**
 * Frees memory allocated for up to 5 matrices and their data arrays.
 * Safely handles NULL pointers to avoid double-free errors.
 * 
 * @param M1 First matrix pointer
 * @param M2 Second matrix pointer
 * @param M3 Third matrix pointer
 * @param M4 Fourth matrix pointer
 * @param M5 Fifth matrix pointer
 * @param M1DataPoints Data array for M1
 * @param M2DataPoints Data array for M2
 * @param M3DataPoints Data array for M3
 * @param M4DataPoints Data array for M4
 * @param M5DataPoints Data array for M5
 */
void freeData(double **M1, double **M2, double **M3, double **M4, double **M5, double *M1DataPoints, double *M2DataPoints, double *M3DataPoints, double *M4DataPoints, double *M5DataPoints);

/**
 * Checks for memory allocation failures in 3 matrices.
 * Frees any allocated memory if failure detected.
 * 
 * @param M1 First matrix pointer
 * @param M2 Second matrix pointer
 * @param M3 Third matrix pointer
 * @param M1DataPoints Data array for M1
 * @param M2DataPoints Data array for M2
 * @param M3DataPoints Data array for M3
 * @return 0 if all allocations successful, 1 if any failed
 */
int checkMallocFailure3(double** M1, double** M2, double** M3, double* M1DataPoints, double* M2DataPoints, double* M3DataPoints);

#endif