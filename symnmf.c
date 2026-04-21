# define _GNU_SOURCE
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include "symnmf.h"


/**
 * Initializes a 2D matrix of doubles with the specified dimensions.
 * Allocates memory for a length x width matrix and stores the data pointer for later freeing.
 * 
 * @param length Number of rows in the matrix
 * @param width Number of columns in the matrix  
 * @param pPointer Pointer to store the address of the allocated data array
 * @return Pointer to the initialized matrix, or NULL if allocation failed
 */
double** initializeMatrix(int length, int width, double** pPointer){
    double **A,*p;
    int i;

    p = calloc(length*width , sizeof(double));
    if (p == NULL)
        return NULL; /*if calloc failed*/
    *pPointer = p; /*saving the address of p for freeing later*/
    A = malloc(length*sizeof(double*));
    if (A == NULL)
       return NULL; /*if malloc failed*/
    for( i=0 ; i<length ; i++ )
        A[i] = p+i*width;
    return A;
}

/**
 * Prints a matrix to stdout in CSV format with 4 decimal places.
 * Each row on a new line, values separated by commas.
 * 
 * @param matrix The matrix to print
 * @param rows Number of rows in the matrix
 * @param cols Number of columns in the matrix
 */
void printMatrix(double** matrix, int rows, int cols){
    int i, j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            printf("%.4f", matrix[i][j]);
            if (j != cols-1) /*if we didn't get to end of line*/
                printf(",");
        }
        printf("\n"); /*space between lines*/
    }
}


/**
 * Reads input data from a file and creates a 2D matrix.
 * Parses CSV format, determines N (rows) and D (columns) automatically.
 * 
 * @param filePointer Open file pointer to read from
 * @param N Pointer to store number of data points (rows)
 * @param D Pointer to store dimensionality (columns)
 * @param inputMatrixDataPoints Pointer to store data array address for freeing
 * @return 2D matrix containing the input data, or NULL on error
 */
double** readInput(FILE* filePointer, int* N, int* D, double** inputMatrixDataPoints){
    double currInput, **inputMatrix = NULL;
    int counterN = 0, counterD = 0, total_elements = 0, i, j;
    char c = '0';
    /*first pass for finding N and D*/
    while (fscanf(filePointer, "%lf%c", &currInput, &c) == 2){
        total_elements++;
        if (counterN == 0) counterD++; /*count D by the first vector*/
        if (c == '\n') counterN++;
    }
    if (total_elements != counterN * counterD){ /*if N*D != all elements, then the input file is wrong*/
        HANDLE_ERROR(); }
    *N = counterN;
    *D = counterD;
    rewind(filePointer);
    inputMatrix = initializeMatrix(*N,*D,inputMatrixDataPoints);
    if (inputMatrix == NULL){
        free(*inputMatrixDataPoints);
        free(inputMatrix);
        HANDLE_ERROR();}
    /*second loop: updating the input matrix*/
    for (i = 0; i < *N; i++) {
        for (j = 0; j < *D; j++) {
            if (fscanf(filePointer, "%lf%c", &inputMatrix[i][j], &c) != 2) {
                /*handles the case where the file ends unexpectedly*/
                HANDLE_ERROR();}}}
    return inputMatrix;
}


/**
 * Calculates the squared Euclidean distance between two vectors.
 * Returns squared distance to avoid unnecessary sqrt computation.
 * 
 * @param vector1 First vector
 * @param vector2 Second vector
 * @param dim Dimensionality of the vectors
 * @return Squared Euclidean distance between the vectors
 */
double euclideanDistance(double* vector1, double* vector2, int dim){
    int i;
    double result = 0;/*final distance*/
    for(i = 0; i < dim; i++){
        result += (vector1[i]-vector2[i])*(vector1[i]-vector2[i]); /*distance between two points squared */
    }
    return result;
}


/**
 * Calculates the sum of all elements in a specified row of a matrix.
 * 
 * @param matrix The input matrix
 * @param cols Number of columns in the matrix
 * @param row The row index to sum
 * @return Sum of elements in the specified row
 */
double sumRowMatrix(double** matrix, int cols,int row ){
    int j;
    double sum = 0;
    for (j = 0; j < cols; j++){
        sum += matrix[row][j];
    }
    return sum;
}


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
void similarityMatrix(double**A, int N, int D, double** inputs){
    double dist_sq;
    int i, j;
    
    for (i = 0; i < N; i++){
        for (j = i+1; j < N; j++){
            dist_sq = euclideanDistance(inputs[i], inputs[j], D);
            A[i][j] = exp(-dist_sq / 2.0);
            A[j][i] = A[i][j]; /*the matrix is symmteric*/
        }
    }
}


/**
 * Computes the diagonal degree matrix D from similarity matrix A.
 * D[i][i] = sum of row i in A, D[i][j] = 0 for i != j.
 * 
 * @param D Pointer to the diagonal matrix to be filled (N x N)
 * @param A Similarity matrix (N x N)
 * @param N Size of the matrices
 */
void diagonalMatrix(double** D, double** A, int N){
    int i;
    for(i=0;i<N;i++){
        D[i][i]= sumRowMatrix(A,N,i);
    }
}


/**
 * Computes the inverse square root of diagonal elements.
 * DTag[i][i] = 1/sqrt(D[i][i]), DTag[i][j] = 0 for i != j.
 * 
 * @param DTag Pointer to the result matrix (N x N)
 * @param matrix Input diagonal matrix D (N x N)
 * @param N Size of the matrices
 */
void raiseNegativeHalfDiag(double** DTag, double** matrix, int N){
    int i;
    for (i = 0; i < N; i++){
        DTag[i][i] = 1.0 / sqrt(matrix[i][i]);
    }
}


/**
 * Computes the normalized similarity matrix W.
 * W = DTag * A * DTag
 * W[i][j] = DTag[i][i] * A[i][j] * DTag[j][j] for i != j, W[i][i] = 0.
 * The matrix is symmetric.
 * 
 * @param W Pointer to the normalized matrix to be filled (N x N)
 * @param A Similarity matrix (N x N)
 * @param DTag Diagonal matrix with inverse square roots (N x N)
 * @param N Size of the matrices
 */
void normalizedSimilarityMatrix(double** W, double** A, double** DTag, int N){
    int i,j;
    for (i = 0;i<N;i++){
        for(j = i + 1;j<N;j++){ /*w_ii = 0 for all i, initialized before*/
            W[i][j] = DTag[i][i]*A[i][j]*DTag[j][j]; /*w_ij = d'_ii * a_ij * d'_jj*/
            W[j][i] = W[i][j]; /*W is symmetric*/
        }
    }
}


/**
 * Calculates the squared Frobenius norm of the difference between two matrices.
 * Returns ||A - B||_F^2 without taking square root.
 * 
 * @param A First matrix
 * @param B Second matrix
 * @param rows Number of rows in matrices
 * @param cols Number of columns in matrices
 * @return Squared Frobenius norm of (A - B)
 */
double frobeniusNormSquared(double** A, double** B, int rows, int cols){
    int i,j;
    double sum = 0;
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            sum += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return sum;
}


/**
 * Computes a single element of matrix multiplication A * B at position (i,j).
 * 
 * @param A First matrix
 * @param B Second matrix
 * @param i Row index in result matrix
 * @param j Column index in result matrix
 * @param ACols Number of columns in A (rows in B)
 * @return Value at position (i,j) in A * B
 */
double matrixMultSingleValue(double** A, double** B, int i, int j, int ACols){
    int k;
    double res = 0;
    for (k = 0; k < ACols; k++){
        res += A[i][k] * B[k][j];
    }
    return res;
}


/**
 * Computes matrix multiplication C = A * B.
 * 
 * @param C Result matrix (ARows x BCols)
 * @param A First matrix (ARows x ACols)
 * @param B Second matrix (ACols x BCols)
 * @param ARows Number of rows in A
 * @param ACols Number of columns in A (rows in B)
 * @param BCols Number of columns in B
 */
void matrixMult(double** C, double** A, double** B, int ARows, int ACols,int BCols){
    int i, j;
    for (i = 0; i < ARows; i++){
        for (j = 0; j < BCols; j++){
           C[i][j] = matrixMultSingleValue(A, B, i, j, ACols);
        }
    }
}


/**
 * Computes the transpose of a matrix.
 * 
 * @param transposeMatrix Result matrix (cols x rows)
 * @param matrix Input matrix (rows x cols)
 * @param rows Number of rows in input matrix
 * @param cols Number of columns in input matrix
 * @return Pointer to the transposed matrix (same as transposeMatrix)
 */
double** matrixTranspose(double** transposeMatrix, double** matrix, int rows, int cols){
    int i,j;
    for (i = 0; i < rows;i++){
        for(j = 0; j < cols;j++){
            transposeMatrix[j][i] = matrix[i][j];
        }
    }
    return transposeMatrix;
}


/**
 * Calculates the updated value for H[i][j] in the SymNMF algorithm.
 * Uses the multiplicative update rule.
 * 
 * @param H Current H matrix
 * @param WMultH Matrix W * H
 * @param HMultHTmultH Matrix H * (H^T * H)
 * @param beta Parameter
 * @param i Row index
 * @param j Column index
 * @return Updated value for H[i][j]
 */
double calcNewHValue(double** H, double** WMultH, double** HMultHTmultH, double beta, int i, int j){
        return H[i][j] * (1 - beta + beta * (WMultH[i][j] / (HMultHTmultH[i][j] + 1e-6))); /*adding 1e-6 to the denominator to not divide by 0*/
}


/**
 * Copies all values from matrix B to matrix A.
 * 
 * @param A Destination matrix
 * @param B Source matrix
 * @param rows Number of rows
 * @param cols Number of columns
 */
void deepCopyMatrix(double** A, double** B, int rows, int cols){
    int i, j;
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            A[i][j] = B[i][j];
        }
    }
}


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
void freeData(double **M1, double **M2, double **M3, double **M4, double **M5, double *M1DataPoints, double *M2DataPoints, double *M3DataPoints, double *M4DataPoints, double *M5DataPoints){
    free(M1);
    free(M2);
    free(M3);
    free(M4);
    free(M5);
    free(M1DataPoints);
    free(M2DataPoints);
    free(M3DataPoints);
    free(M4DataPoints);
    free(M5DataPoints);
}


/**
 * Checks for memory allocation failures in 5 matrices.
 * Frees any allocated memory if failure detected.
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
 * @return 0 if all allocations successful, 1 if any failed
 */
int checkMallocFailure5(double** M1, double** M2, double** M3,double** M4, double** M5, double* M1DataPoints, double* M2DataPoints, double* M3DataPoints, double* M4DataPoints, double* M5DataPoints){
    if (M1 == NULL || M2 == NULL || M3 == NULL || M4 == NULL || M5 == NULL || M1DataPoints == NULL || M2DataPoints == NULL || M3DataPoints == NULL || M4DataPoints == NULL || M5DataPoints == NULL){  /*Matrix being null means malloc failed*/
        freeData(M1,M2,M3,M4,M5,M1DataPoints,M2DataPoints,M3DataPoints,M4DataPoints,M5DataPoints);
        return 1;
    }
    return 0;
}


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
int symnmfAlgorithm(double** W, double** H, int N, int k){
    double *newHDataPoints = NULL, *HTDataPoints = NULL, *WMultHDataPoints = NULL, *HTMultHDataPoints = NULL, *HMultHTmultHDataPoints = NULL;
    double **newH = initializeMatrix(N, k, &newHDataPoints), **WMultH = initializeMatrix(N, k, &WMultHDataPoints), **HT = initializeMatrix(k, N, &HTDataPoints), **HTMultH=initializeMatrix(k, k, &HTMultHDataPoints),**HMultHTmultH=initializeMatrix(N, k, &HMultHTmultHDataPoints);
    double eps = 1e-4, beta = 0.5;
    int maxIter = 300, currIter = 0, i, j;
    /*checking if any malloc failed when initializing matrices*/
    if (checkMallocFailure5(newH, WMultH, HT, HTMultH, HMultHTmultH, newHDataPoints, WMultHDataPoints, HTDataPoints, HTMultHDataPoints, HMultHTmultHDataPoints) == 1){
        return 1;
    }
    do {
        if(currIter!=0){ /*if not first iteration then update H^(i)=newH^(i-1)*/
            deepCopyMatrix(H,newH,N,k); /* H gets newH values*/
        }
        matrixMult(WMultH, W, H, N, N, k); /*WMultH = W*H*/
        matrixTranspose(HT, H, N, k); /*HT=H^T*/
        matrixMult(HTMultH, HT, H, k, N, k); /*HTMultH = (H^T)*H*/
        matrixMult(HMultHTmultH, H, HTMultH, N, k, k); /*HMultHTmultH = H*((H^T)*H)*/
        
        for (i = 0; i < N; i++){
            for(j = 0; j < k; j++){
                newH[i][j] = calcNewHValue(H,WMultH, HMultHTmultH, beta, i, j);
            }
        }
        currIter++;
    } while (currIter < maxIter && (frobeniusNormSquared(H, newH, N, k) >= eps));
    deepCopyMatrix(H,newH,N,k); /* H gets newH values*/
    /*Freeing memory*/
    freeData(HT, WMultH, HTMultH,HMultHTmultH, newH, HTDataPoints,WMultHDataPoints,HTMultHDataPoints,HMultHTmultHDataPoints,newHDataPoints);
    return 0;
}


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
int checkMallocFailure3(double** M1, double** M2, double** M3, double* M1DataPoints, double* M2DataPoints, double* M3DataPoints){
    if (M1 == NULL || M2 == NULL || M3 == NULL || M1DataPoints == NULL || M2DataPoints == NULL || M3DataPoints == NULL){  /*Matrix being null means malloc failed*/
        freeData(M1,M2,M3,NULL,NULL,M1DataPoints,M2DataPoints,M3DataPoints,NULL,NULL);
        return 1;
    }
    return 0;
}


/**
 * Processes the specified goal operation on input data.
 * Computes and prints the appropriate matrix based on goal.
 * 
 * @param goal Operation to perform: "sym", "ddg", "norm", or "symnmf"
 * @param inputMatrix Input data points (N x d)
 * @param N Number of data points
 * @param d Dimensionality of data
 * @return 0 on success, 1 on error or invalid goal
 */
int handleGoal(char* goal, double** inputMatrix, int N, int d) {
    double **A = NULL, *ADataPoints = NULL, **D = NULL, *DDataPoints = NULL, **W = NULL, *WDataPoints = NULL;
    /*initalizing the matrices*/
    A = initializeMatrix(N, N, &ADataPoints);
    D = initializeMatrix(N, N, &DDataPoints);
    W = initializeMatrix(N, N, &WDataPoints);
    if (checkMallocFailure3(A, D, W, ADataPoints, DDataPoints, WDataPoints) == 1){
        return 1;}
    similarityMatrix(A, N, d, inputMatrix); /*creating matrix A*/
    if (strcmp(goal, "sym") == 0){ /*goal == sym*/
        printMatrix(A, N, N);
        freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
        return 0;}
    diagonalMatrix(D,A, N);/*creating matrix D*/
    if (strcmp(goal, "ddg") == 0){ /*goal == ddg*/
        printMatrix(D, N, N);
        freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
        return 0;}
    raiseNegativeHalfDiag(D, D, N); /*updates D = D^(-1/2)*/
    normalizedSimilarityMatrix(W,A, D, N); /*creating matrix W*/
    if (strcmp(goal, "norm") == 0){ /*goal == norm*/
        printMatrix(W, N, N);
        freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
        return 0;}
    /*if goal is not sym, ddg or norm*/
    freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
    return 1;
}


/**
 * Main entry point for the SymNMF command-line program.
 * Expects 2 arguments: goal and input file path.
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line arguments
 * @return 0 on success, 1 on error
 */
int main(int argc, char **argv){
    /*intializing variables*/
    FILE *ifp = NULL;
    char *goal;
    double **inputMatrix=NULL, *inputMatrixDataPoints=NULL;
    int N, d, result;
    /*checking the arguments*/
    if (argc != 3){
       HANDLE_ERROR();
    }
    goal = argv[1];
    /*opening the file*/
    ifp = fopen(argv[2], "r");
    if (ifp == NULL){
        HANDLE_ERROR();
    }
    inputMatrix = readInput(ifp, &N, &d, &inputMatrixDataPoints); /*reading the input from the file*/
    if (fclose(ifp) == EOF){ /*closing the file*/
        HANDLE_ERROR();
    }
    /*handling goal*/
    result = handleGoal(goal, inputMatrix, N, d);
    /*freeing memory*/
    free(inputMatrixDataPoints);
    free(inputMatrix);
    if (result == 1){
        HANDLE_ERROR();
    }
    return 0;
}