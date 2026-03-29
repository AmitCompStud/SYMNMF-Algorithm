# define _GNU_SOURCE
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# define HANDLE_ERROR() do { \
    printf("An Error Has Occurred\n"); \
    exit(1); \
} while (0)

double* readInput(FILE* filePointer, int* N){
    //TODO: fix.
    double currInput, *inputVector;  
    int counter = 0;
    //first loop: counting how many vectors there are in the file
    while (fscanf(filePointer, "%lf", &currInput) == 1)
    {
        counter++; 
    }
    rewind(filePointer);
    inputVector = (double*)malloc(counter * sizeof(double));
    if (inputVector == NULL){
        HANDLE_ERROR();
    }

    //second loop: putting the numbers in the input vector
    counter = 0;//Reseting the counter
    while (fscanf(filePointer, "%lf", &currInput) == 1)
    {
        inputVector[counter++] = currInput;
    }
    *N = counter;
    return inputVector;
}


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



int main(int argc, char **argv){
    //intializing variables
    FILE *ifp = NULL;
    char *goal;
    double **inputMatrix, **A, *ADataPoints;
    int i, N;
    
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
    
    inputMatrix = readInput(ifp, &N); //intialzing the input vector

    //creating matrix A
    A=initalizeMatrix(N,N,&ADataPoints);

        
}