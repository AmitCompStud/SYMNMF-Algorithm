# define _GNU_SOURCE
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# define HANDLE_ERROR() do { \
    printf("An Error Has Occurred\n"); \
    exit(1); \
} while (0)


double** initalizeMatrix(int length, int width, double** pPointer){
    double **A,*p;
    int i;

    p = calloc(length*width , sizeof(double));
    A = malloc(length*sizeof(double*));
    for( i=0 ; i<length ; i++ )
        A[i] = p+i*width;
    *pPointer = p; //saving the adress of p for freeing later
    return A;
}

void printMatrix(double** matrix, int rows, int cols){
    int i, j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            printf("%.4f", matrix[i][j]);
            if (j != cols-1) //if we didn't get to end of line
                printf(",");
        }
        printf("\n"); //space between lines
    }
}


double** readInput(FILE* filePointer, int* N, int* D, double** inputMatrixDataPoints){
    double currInput, **inputMatrix;  
    int counterN = 0, counterD = 0, total_elements = 0, i, j;
    char c = '0';
    //first pass for finding N and D
    while (fscanf(filePointer, "%lf%c", &currInput, &c) == 2){
        total_elements++;
        if (counterN == 0){ 
            counterD++; //count D by the first vector
        }
        if (c == '\n'){
            counterN++;
        }
    }
    if (total_elements != counterN * counterD){ //if N*D != all elements, then the input file is wrong
        HANDLE_ERROR();
    }
    *N = counterN;
    *D = counterD;
    rewind(filePointer);
    inputMatrix = initalizeMatrix(*N,*D,inputMatrixDataPoints);
    if (inputMatrix == NULL){
        HANDLE_ERROR();
    }
    //second loop: updating the input matrix
    for (i = 0; i < *N; i++) {
        for (j = 0; j < *D; j++) {
            if (fscanf(filePointer, "%lf%c", &inputMatrix[i][j], &c) != 2) {
                //handles the case where the file ends unexpectedly
                HANDLE_ERROR();
            }
        }
    }
    return inputMatrix;
}


double euclideanDistance(double* vector1, double* vector2, int dim){
    int i;
    double result = 0;//final distance
    for(i = 0; i < dim; i++){
        result += (vector1[i]-vector2[i])*(vector1[i]-vector2[i]); //distance between two points squared 
    }
    return result;
}


double** similarityMatrix(int N, int D, double** inputs, double** dataPoints){
    double** A=initalizeMatrix(N,N,dataPoints), dist_sq;
    int i, j;
    
    for (i = 0; i < N; i++){
        for (j = i+1; j < N; j++){
            dist_sq = euclideanDistance(inputs[i], inputs[j], D);
            A[i][j] = exp(-dist_sq / 2.0);
            A[j][i] = A[i][j]; //the matrix is symmteric
        }
    }
    return A;

}


int main(int argc, char **argv){
    //intializing variables
    FILE *ifp = NULL;
    char *goal;
    double **inputMatrix, *inputMatrixDataPoints, **A, *ADataPoints;
    int N, d;
    
    //checking the arguments
    if (argc != 3){
       HANDLE_ERROR();
    }
    goal = argv[1];
    //opening the file
    ifp = fopen(argv[2], "r");
    if (ifp == NULL){
        HANDLE_ERROR();
    }
    
    inputMatrix = readInput(ifp, &N, &d, &inputMatrixDataPoints); //reading the input from the file
    if (fclose(ifp) == EOF){ //closing the file
        HANDLE_ERROR();
    }

    //creating matrix A
    A=similarityMatrix(N,d, inputMatrix, &ADataPoints);

    //if goal is sym
    printMatrix(A);


        
}