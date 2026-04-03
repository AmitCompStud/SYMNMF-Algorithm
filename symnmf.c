# define _GNU_SOURCE
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include "symnmf.h"


double** initializeMatrix(int length, int width, double** pPointer){
    double **A,*p;
    int i;

    p = calloc(length*width , sizeof(double));
    if (p == NULL)
        return NULL; /*if calloc failed*/
    A = malloc(length*sizeof(double*));
    if (A == NULL)
       return NULL; /*if malloc failed*/
    for( i=0 ; i<length ; i++ )
        A[i] = p+i*width;
    *pPointer = p; /*saving the address of p for freeing later*/
    return A;
}

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


double** readInput(FILE* filePointer, int* N, int* D, double** inputMatrixDataPoints){
    double currInput, **inputMatrix = NULL;
    int counterN = 0, counterD = 0, total_elements = 0, i, j;
    char c = '0';
    /*first pass for finding N and D*/
    while (fscanf(filePointer, "%lf%c", &currInput, &c) == 2){
        total_elements++;
        if (counterN == 0){ 
            counterD++; /*count D by the first vector*/
        }
        if (c == '\n'){
            counterN++;
        }
    }
    if (total_elements != counterN * counterD){ /*if N*D != all elements, then the input file is wrong*/
        HANDLE_ERROR();
    }
    *N = counterN;
    *D = counterD;
    rewind(filePointer);
    inputMatrix = initializeMatrix(*N,*D,inputMatrixDataPoints);
    if (inputMatrix == NULL){
        free(*inputMatrixDataPoints);
        free(inputMatrix);
        HANDLE_ERROR();
    }
    /*second loop: updating the input matrix*/
    for (i = 0; i < *N; i++) {
        for (j = 0; j < *D; j++) {
            if (fscanf(filePointer, "%lf%c", &inputMatrix[i][j], &c) != 2) {
                /*handles the case where the file ends unexpectedly*/
                HANDLE_ERROR();
            }
        }
    }
    return inputMatrix;
}


double euclideanDistance(double* vector1, double* vector2, int dim){
    int i;
    double result = 0;/*final distance*/
    for(i = 0; i < dim; i++){
        result += (vector1[i]-vector2[i])*(vector1[i]-vector2[i]); /*distance between two points squared */
    }
    return result;
}


double sumRowMatrix(double** matrix, int cols,int row ){
    int j;
    double sum = 0;
    for (j = 0; j < cols; j++){
        sum += matrix[row][j];
    }
    return sum;
}


void similarityMatrix(double**A, int N, int D, double** inputs){ /*updates matrix A*/
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


void diagonalMatrix(double** D, double** A, int N){ /*updates matrix D*/
    int i;
    for(i=0;i<N;i++){
        D[i][i]= sumRowMatrix(A,N,i);
    }
}


void raiseNegativeHalfDiag(double** DTag, double** matrix, int N){ /*updates DTag*/
    int i;
    for (i = 0; i < N; i++){
        DTag[i][i] = 1.0 / sqrt(matrix[i][i]);
    }
}


void normalizedSimilarityMatrix(double** W, double** A, double** DTag, int N){ /*updates matrix W*/
    int i,j;
    for (i = 0;i<N;i++){
        for(j = i + 1;j<N;j++){ /*w_ii = 0 for all i, initialized before*/
            W[i][j] = DTag[i][i]*A[i][j]*DTag[j][j]; /*w_ij = d'_ii * a_ij * d'_jj*/
            W[j][i] = W[i][j]; /*W is symmetric*/
        }
    }
}


double frobeniusNormSquared(double** A, double** B, int rows, int cols){
    int i,j;
    double sum;
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            sum += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return sum;
}


double matrixMultSingleValue(double** A, double** B, int i, int j, int ACols){
    int k;
    double res = 0;
    for (k = 0; k < ACols; k++){
        res += A[i][k] * B[k][j];
    }
    return res;
}


void matrixMult(double** C, double** A, double** B, int ARows, int ACols,int BCols){/*C=A*B*/
    int i, j;
    for (i = 0; i < ARows; i++){
        for (j = 0; i < BCols; j++){
           C[i][j] = matrixMultSingleValue(A, B, i, j, ACols);
        }
    }
}


double** matrixTranspose(double** transposeMatrix, double** matrix, int rows, int cols){
    int i,j;
    for (i = 0; i < rows;i++){
        for(j = 0; j < cols;j++){
            transposeMatrix[j][i] = matrix[i][j];
        }
    }
    return transposeMatrix;
}


double calcNewHValue(double** H, double** WMultH, double** HMultHTmultH, double beta, int i, int j){
        return H[i][j] * (1 - beta + beta * (WMultH[i][j] / (HMultHTmultH[i][j] + 1e-6))); /*adding 1e-6 to the denominator to not divide by 0*/
}


void deepCopyMatrix(double** A, double** B, int rows, int cols){ /*A[i][j] = B[i][j] for all i,j*/
    int i, j;
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            A[i][j] = B[i][j];
        }
    }
}


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


int checkMallocFailure5(double** M1, double** M2, double** M3,double** M4, double* M5, double* M1DataPoints, double* M2DataPoints, double* M3DataPoints, double* M4DataPoints, double* M5DataPoints){
    if (M1 == NULL || M2 == NULL || M3 == NULL || M4 == NULL || M5 == NULL || M1DataPoints == NULL || M2DataPoints == NULL || M3DataPoints == NULL || M4DataPoints == NULL || M5DataPoints == NULL){  /*Matrix being null means malloc failed*/
        freeData(M1,M2,M3,M4,M5,M1DataPoints,M2DataPoints,M3DataPoints,M4DataPoints,M5DataPoints);
        return 1;
    }
    return 0;
}


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
        matrixMult(HMultHTmultH, H, HT, N, k, k); /*HMultHTmultH = H*((H^T)*H)*/
        
        for (i = 0; i < N; j++){
            for(j=0;j<k;j++){
                newH[i][j] = calcNewHValue(H,WMultH, HMultHTmultH,beta, i,j);
            }
        }
        currIter++;
    } while (currIter < maxIter && (frobeniusNormSquared(H, newH, N, k) >= eps));
    
    deepCopyMatrix(H,newH,N,k); /* H gets newH values*/
     
    /*Freeing memory*/
    freeData(HT, WMultH, HTMultH,HMultHTmultH, newH, HTDataPoints,WMultHDataPoints,HTMultHDataPoints,HMultHTmultHDataPoints,newHDataPoints);
    return 0;
}


int checkMallocFailure3(double** M1, double** M2, double** M3, double* M1DataPoints, double* M2DataPoints, double* M3DataPoints){
    if (M1 == NULL || M2 == NULL || M3 == NULL || M1DataPoints == NULL || M2DataPoints == NULL || M3DataPoints == NULL){  /*Matrix being null means malloc failed*/
        freeData(M1,M2,M3,NULL,NULL,M1DataPoints,M2DataPoints,M3DataPoints,NULL,NULL);
        return 1;
    }
    return 0;
}


int handleGoal(char* goal, double** inputMatrix, int N, int d) {
    double **A = NULL, *ADataPoints = NULL, **D = NULL, *DDataPoints = NULL, **W = NULL, *WDataPoints = NULL;
    /*initalizing the matrices*/
    A = initializeMatrix(N, N, &ADataPoints);
    D = initializeMatrix(N, N, &DDataPoints);
    W = initializeMatrix(N, N, &WDataPoints);
    if (checkMallocFailure3(A, D, W, ADataPoints, DDataPoints, WDataPoints) == 1){
        return 1;
    }
    
    similarityMatrix(A, N, d, inputMatrix); /*creating matrix A*/

    if (strcmp(goal, "sym") == 0){ /*goal == sym*/
        printMatrix(A, N, N);
        freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
        return 0;
    }

    diagonalMatrix(D,A, N);/*creating matrix D*/

    if (strcmp(goal, "ddg") == 0){ /*goal == ddg*/
        printMatrix(D, N, N);
        freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
        return 0;
    }

    raiseNegativeHalfDiag(D, D, N); /*updates D = D^(-1/2)*/
    normalizedSimilarityMatrix(W,A, D, N); /*creating matrix W*/

    if (strcmp(goal, "norm") == 0){ /*goal == norm*/
        printMatrix(W, N, N);
        freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
        return 0;
    }
    /*if goal is not sym, ddg or norm*/
    freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
    return 1;
}


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