#define PY_SSIZE_T_CLEAN
#include <Python.h>


static double* PyListToRow(PyObject* lst){
    int i, n = PyObject_Length(lst);
    double* row = malloc(n * sizeof(double));
    //TODO: error
    for (i = 0; i < n; i++){
        row[i] = PyList_GetItem(lst, i);
    }
    return row;
}


static double** PyMatrixToMatrix(PyObject* matrix){
    int i, n = PyObject_Length(lst);
    double** returnMatrix = malloc(n *  sizeof(double*));
    //TODO: error
    for (i = 0; i < n; i++){
        returnMatrix[i] = PyListToRow(PyList_GetItem(matrix, i));
    }
    return returnMatrix;
}


static PyObject* arrToPyList(double* arr, int length){
    PyObject* pyList = PyList_New(length), *pyDouble;
    int i;
    for (i = 0; i < length; i++){
        pyDouble = Py_BuildValue("d", arr[i]);
        PyList_SetItem(pyList, i, pyDouble);
    }
    return PyList;
}


static PyObject* matrixToPyMatrix(double** matrix, int rows, int cols){
    PyObject* pyMatrix = PyList_New(rows);
    int i;
    for (i = 0; i < rows; i++){
        PyList_SetItem(arrToPyList(matrix[i], cols), i);
    }
    return pyMatrix;
}


static **double convertPyMatrixToCMatrix(PyObject* args, int* N, int* d){
    double** inputMatrix=NULL;
    PyObject* pyMatrix;
    PyArg_ParseTuple(args, "O", &pyMatrix);
    *N = PyObject_Length(pyMatrix); //amount of rows in the matrix
    *d = PyObject_Length(PyList_GetItem(pyMatrix, 0)); //length of the first row in matrix
    inputMatrix = PyMatrixToMatrix(pyMatrix);
    return inputMatrix;
}


void freeInputMatrix(double** inputMatrix, int rows){
    int i;
    for(i=0;i<rows;i++){
        free(inputMatrix[i]);
    }
    free(inputMatrix);
}


static PyObject* sym(PyObject* self, PyObject *args){
    double** inputMatrix, **A, *ADataPoints;
    PyObject* pyMatrix;
    int N, d;

    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);

    A = similarityMatrix(N, d, inputMatrix, &ADataPoints); //creating A
    pyMatrix = matrixToPyMatrix(A, N, N); //turning A into python matrix
    //freeing memory
    free(A);
    free(ADataPoints);
    freeInputMatrix(inputMatrix, N);
    return pyMatrix;
}


static PyObject* ddg(PyObject* self, PyObject *args){
    double** inputMatrix, **A, *ADataPoints, **D, *DDataPoints;
    PyObject* pyMatrix;
    int N, d;
    
    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);

    A = similarityMatrix(N, d, inputMatrix, &ADataPoints); //creating A
    D = diagonalMatrix(A, N, &DDataPoints); //creating D
    pyMatrix = matrixToPyMatrix(D, N, N); //turning D into python matrix
    //freeing memory
    free(A);
    free(ADataPoints);
    free(D);
    free(DDataPoints);
    freeInputMatrix(inputMatrix, N);
    return pyMatrix;
}


static PyObject* norm(PyObject* self, PyObject *args){
    double** inputMatrix, **A, *ADataPoints, **D, *DDataPoints, **DTag, *DTagDataPoints, **W, *WDataPoints;
    PyObject* pyMatrix;
    int N, d;
    
    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);

    A = similarityMatrix(N, d, inputMatrix, &ADataPoints); //creating A
    D = diagonalMatrix(A, N, &DDataPoints); //creating D
    DTag = raiseNegativeHalfDiag(D, N, &DTagDataPoints); //D' = D^(-0.5)
    W = normalizedSimilarityMatrix(A, DTag, N, &WDataPoints); //creating W
    pyMatrix = matrixToPyMatrix(W, N, N); //turning W into python matrix
    //freeing memory
    free(A);
    free(ADataPoints);
    free(D);
    free(DDataPoints);
    free(DTag);
    free(DTagDataPoints);
    free(W);
    free(WDataPoints);
    freeInputMatrix(inputMatrix, N);
    return pyMatrix;
}

static PyObject* symnmf(PyObject* self, PyObject *args){
    double **W, **H;
    PyObject* pyMatrixW, *pyMatrixH;
    int N, k;
    
    //converting python W and H to C matrices
    PyArg_ParseTuple(args, "OO", &pyMatrixW, &pyMatrixH);
    N = PyObject_Length(pyMatrixW); //amount of rows in the matrix
    k = PyObject_Length(PyList_GetItem(pyMatrixH, 0)); //length of the first row in matrix H
    W = PyMatrixToMatrix(pyMatrixW);
    H = PyMatrixToMatrix(pyMatrixH);
    
    //freeing memory

    free(WDataPoints);
    return pyMatrix;
}


