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
        HANDLE_ERROR();
    A = malloc(length*sizeof(double*));
    if (A == NULL)
        HANDLE_ERROR();
    for( i=0 ; i<length ; i++ )
        A[i] = p+i*width;
    *pPointer = p; /*saving the adress of p for freeing later*/
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
    double currInput, **inputMatrix;
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


double** similarityMatrix(int N, int D, double** inputs, double** dataPoints){
    double** A = initializeMatrix(N,N,dataPoints), dist_sq;
    int i, j;
    
    for (i = 0; i < N; i++){
        for (j = i+1; j < N; j++){
            dist_sq = euclideanDistance(inputs[i], inputs[j], D);
            A[i][j] = exp(-dist_sq / 2.0);
            A[j][i] = A[i][j]; /*the matrix is symmteric*/
        }
    }
    return A;

}


double** diagonalMatrix(double** A, int N, double** dataPoints){
    double** D = initializeMatrix(N,N, dataPoints);
    int i;
    for(i=0;i<N;i++){
        D[i][i]= sumRowMatrix(A,N,i);
    }
    return D;
}


double** raiseNegativeHalfDiag(double** matrix, int N, double** dataPoints){
    double** DTag = initializeMatrix(N, N, dataPoints);
    int i;
    for (i = 0; i < N; i++){
        DTag[i][i] = 1.0 / sqrt(matrix[i][i]);
    }
    return DTag;
}


double** normalizedSimilarityMatrix(double** A, double** DTag, int N, double ** dataPoints){
    double** W = initializeMatrix(N,N, dataPoints);
    int i,j;
    for (i = 0;i<N;i++){
        for(j = i + 1;j<N;j++){ /*w_ii = 0 for all i, initialized before*/
            W[i][j] = DTag[i][i]*A[i][j]*DTag[j][j]; /*w_ij = d'_ii * a_ij * d'_jj*/
            W[j][i] = W[i][j]; /*W is symmetric*/
        }
    }
    return W;
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


void symnmfAlgorithm(double** W, double** H, int N, int k){
    double *newHDataPoints, *HTDataPoints, *WMultHDataPoints, *MTMultHDataPoints, *HMultHTmultHDataPoints;
    double **newH = initializeMatrix(N, k, &newHDataPoints), **WMultH = initializeMatrix(N, k, &WMultHDataPoints), **HT = initializeMatrix(k, N, &HTDataPoints), **HTMultH=initializeMatrix(k, k, &MTMultHDataPoints),**HMultHTmultH=initializeMatrix(N, k, &HMultHTmultHDataPoints);
    double eps = 1e-4, beta = 0.5;
    int maxIter = 300, currIter = 0, i, j;
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
    free(HT);
    free(HTDataPoints);
    free(WMultH);
    free(WMultHDataPoints);
    free(HMultHTmultH);
    free(HMultHTmultHDataPoints);
    free(newH);
    free(newHDataPoints);
}


void freeData(double **A, double *ADataPoints, double **D, double *DDataPoints, double **DTag, double *DTagDataPoints, double **W, double *WDataPoints){
        free(ADataPoints);
        free(A);
        free(DDataPoints);
        free(D);
        free(DTagDataPoints);
        free(DTag);
        free(WDataPoints);
        free(W);
}


int handleGoal(char* goal, double** inputMatrix, int N, int d) {
    double **A = NULL, *ADataPoints = NULL, **D = NULL, *DDataPoints = NULL, **DTag = NULL, *DTagDataPoints = NULL, **W = NULL, *WDataPoints = NULL;
    /*creating matrix A*/
    A = similarityMatrix(N, d, inputMatrix, &ADataPoints);

    if (strcmp(goal, "sym") == 0){ /*goal == sym*/
        printMatrix(A, N, N);
        freeData(A, ADataPoints, D, DDataPoints, DTag, DTagDataPoints, W, WDataPoints);
        return 0;
    }

    /*creating matrix D*/
    D = diagonalMatrix(A, N, &DDataPoints);

    if (strcmp(goal, "ddg") == 0){ /*goal == ddg*/
        printMatrix(D, N, N);
        freeData(A, ADataPoints, D, DDataPoints, DTag, DTagDataPoints, W, WDataPoints);
        return 0;
    }

    /*creating matrix W*/
    /*D'= D^(-1/2)*/
    DTag = raiseNegativeHalfDiag(D,N, &DTagDataPoints);
    W = normalizedSimilarityMatrix(A, DTag, N, &WDataPoints);

    if (strcmp(goal, "norm") == 0){ /*goal == norm*/
        printMatrix(W, N, N);
        freeData(A, ADataPoints, D, DDataPoints, DTag, DTagDataPoints, W, WDataPoints);
        return 0;
    }
    /*if goal is not sym, ddg or norm*/
    freeData(A, ADataPoints, D, DDataPoints, DTag, DTagDataPoints, W, WDataPoints);
    return 1;
}


int main(int argc, char **argv){
    /*intializing variables*/
    FILE *ifp = NULL;
    char *goal;
    double **inputMatrix, *inputMatrixDataPoints;
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
    result = handleGoal(goal, inputMatrix, N, d) == 1;
    /*freeing memory*/
    free(inputMatrixDataPoints);
    free(inputMatrix);

    if (result == 1){
        HANDLE_ERROR();
    }
    
    return 0;


}