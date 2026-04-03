#ifndef SYMNMF_H_ 
#define SYMNMF_H_
# define HANDLE_ERROR() do { \
    printf("An Error Has Occurred\n"); \
    exit(1); \
} while (0)

double** initializeMatrix(int length, int width, double** pPointer);

void similarityMatrix(double**A, int N, int D, double** inputs);

void diagonalMatrix(double** D, double** A, int N);

void raiseNegativeHalfDiag(double** DTag, double** matrix, int N);

void normalizedSimilarityMatrix(double** W, double** A, double** DTag, int N);

int symnmfAlgorithm(double** W, double** H, int N, int k);

void freeData(double **M1, double **M2, double **M3, double **M4, double **M5, double *M1DataPoints, double *M2DataPoints, double *M3DataPoints, double *M4DataPoints, double *M5DataPoints);

int checkMallocFailure3(double** M1, double** M2, double** M3, double* M1DataPoints, double* M2DataPoints, double* M3DataPoints);

#endif