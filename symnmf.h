#ifndef SYMNMF_H_ 
#define SYMNMF_H_
# define HANDLE_ERROR() do { \
    printf("An Error Has Occurred\n"); \
    exit(1); \
} while (0)

double** similarityMatrix(int N, int D, double** inputs, double** dataPoints);

double** diagonalMatrix(double** A, int N, double** dataPoints);

double** normalizedSimilarityMatrix(double** A, double** DTag, int N, double ** dataPoints);

void symnmfAlgorithm(double** W, double** H, int N, int k);

double** raiseNegativeHalfDiag(double** matrix, int N, double** dataPoints);

#endif