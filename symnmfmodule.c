# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "symnmf.h"


static void freeMatrix(double** matrix, int rows){
    /**
     * Frees a 2D matrix and its rows.
     * 
     * @param matrix The matrix to free
     * @param rows Number of rows in the matrix
     */
    int i;
    for(i=0;i<rows;i++){
        free(matrix[i]);
    }
    free(matrix);
}


static double* PyListToRow(PyObject* lst){
    /**
     * Converts a Python list to a C double array.
     * 
     * @param lst Python list of floats
     * @return C array of doubles, or NULL on failure
     */
    int i, n = PyObject_Length(lst);
    double* row = malloc(n * sizeof(double));
    if (row == NULL){ /*if malloc failed*/
        return NULL;
    }
    for (i = 0; i < n; i++){
        row[i] = PyFloat_AsDouble(PyList_GetItem(lst, i));
    }
    return row;
}


static double** PyMatrixToMatrix(PyObject* matrix){
    /**
     * Converts a Python list of lists to a C 2D matrix.
     * 
     * @param matrix Python matrix (list of lists)
     * @return C 2D matrix, or NULL on failure
     */
    int i, n = PyObject_Length(matrix);
    double** returnMatrix = malloc(n * sizeof(double*));
    if (returnMatrix == NULL){
        return NULL;
    }
    for (i = 0; i < n; i++){
        returnMatrix[i] = NULL;
    }
    for (i = 0; i < n; i++){
        returnMatrix[i] = PyListToRow(PyList_GetItem(matrix, i));
        if (returnMatrix[i] == NULL){ /*if malloc failed on row i, free data and return NULL*/
            freeMatrix(returnMatrix, n);
            return NULL;
        }
    }
    return returnMatrix;
}


static PyObject* arrToPyList(double* arr, int length){
    /**
     * Converts a C double array to a Python list.
     * 
     * @param arr C array of doubles
     * @param length Length of the array
     * @return Python list of floats
     */
    PyObject* pyList = PyList_New(length), *pyDouble;
    int i;
    for (i = 0; i < length; i++){
        pyDouble = Py_BuildValue("d", arr[i]);
        PyList_SetItem(pyList, i, pyDouble);
    }
    return pyList;
}


static PyObject* matrixToPyMatrix(double** matrix, int rows, int cols){
    /**
     * Converts a C 2D matrix to a Python list of lists.
     * 
     * @param matrix C 2D matrix
     * @param rows Number of rows
     * @param cols Number of columns
     * @return Python matrix (list of lists)
     */
    PyObject* pyMatrix = PyList_New(rows);
    int i;
    for (i = 0; i < rows; i++){
        PyList_SetItem(pyMatrix, i, arrToPyList(matrix[i], cols));
    }
    return pyMatrix;
}


static double** convertPyMatrixToCMatrix(PyObject* args, int* N, int* d){
    /**
     * Converts Python matrix argument to C matrix and extracts dimensions.
     * 
     * @param args Python arguments tuple
     * @param N Pointer to store number of rows
     * @param d Pointer to store number of columns
     * @return C matrix, or NULL on failure
     */
    double** inputMatrix=NULL;
    PyObject* pyMatrix;
    PyArg_ParseTuple(args, "O", &pyMatrix);
    *N = PyObject_Length(pyMatrix); /*amount of rows in the matrix*/
    *d = PyObject_Length(PyList_GetItem(pyMatrix, 0)); /*length of the first row in matrix*/
    inputMatrix = PyMatrixToMatrix(pyMatrix);
    return inputMatrix;
}


static PyObject* sym(PyObject* self, PyObject *args){
    /**
     * Computes the similarity matrix A from input data.
     * 
     * Python signature: sym(matrix)
     * @param self Module object
     * @param args Python tuple containing one matrix argument
     * @return Similarity matrix as Python list of lists
     */
    double** inputMatrix, **A = NULL, *ADataPoints = NULL;
    PyObject* pyMatrix;
    int N, d;

    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);
    if (inputMatrix == NULL){
        /*send error to python*/
        PyErr_NoMemory(); 
        return NULL;
    }
    /*initalizing matrix A*/
    A = initializeMatrix(N, N, &ADataPoints);
     /*if malloc failed*/
    if (A == NULL || ADataPoints == NULL){
        /*free input and send error to python*/
        free(A);
        free(ADataPoints);
        freeMatrix(inputMatrix, N);
        PyErr_NoMemory(); 
        return NULL;
    }
    similarityMatrix(A, N, d, inputMatrix); /*creating A*/
    pyMatrix = matrixToPyMatrix(A, N, N); /*turning A into python matrix*/
    /*freeing memory*/
    free(A);
    free(ADataPoints);
    freeMatrix(inputMatrix, N);
    return pyMatrix;
}


static PyObject* ddg(PyObject* self, PyObject *args){
    /**
     * Computes the diagonal degree matrix D from input data.
     * 
     * Python signature: ddg(matrix)
     * @param self Module object
     * @param args Python tuple containing one matrix argument
     * @return Diagonal degree matrix as Python list of lists
     */
    double** inputMatrix, **A, *ADataPoints, **D, *DDataPoints;
    PyObject* pyMatrix;
    int N, d;
    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);
    if (inputMatrix == NULL){
        /*send error to python*/
        PyErr_NoMemory(); 
        return NULL;}
    /*initalizing matrices*/
    A = initializeMatrix(N, N, &ADataPoints);
    D = initializeMatrix(N, N, &DDataPoints);
    /*if malloc failed*/
    if(A==NULL||ADataPoints==NULL||D==NULL||DDataPoints==NULL){
        free(A);
        free(ADataPoints);
        free(D);
        free(DDataPoints);
        freeMatrix(inputMatrix, N);
        PyErr_NoMemory(); 
        return NULL;}
    similarityMatrix(A, N, d, inputMatrix); /*creating A*/
    diagonalMatrix(D, A, N); /*creating D*/
    pyMatrix = matrixToPyMatrix(D, N, N); /*turning D into python matrix*/
    /*freeing memory*/
    free(A);
    free(ADataPoints);
    free(D);
    free(DDataPoints);
    freeMatrix(inputMatrix, N);
    return pyMatrix;
}


static PyObject* norm(PyObject* self, PyObject *args){
    /**
     * Computes the normalized similarity matrix W from input data.
     * 
     * Python signature: norm(matrix)
     * @param self Module object
     * @param args Python tuple containing one matrix argument
     * @return Normalized similarity matrix as Python list of lists
     */
    double** inputMatrix = NULL, **A = NULL, *ADataPoints = NULL, **D = NULL, *DDataPoints = NULL, **W = NULL, *WDataPoints = NULL;
    PyObject* pyMatrix;
    int N, d;
    
    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);
    if (inputMatrix == NULL){
        /*send error to python*/
        PyErr_NoMemory(); 
        return NULL;
    }
    /*initalizing matrices*/
    A = initializeMatrix(N, N, &ADataPoints);
    D = initializeMatrix(N, N, &DDataPoints);
    W = initializeMatrix(N, N, &WDataPoints);
    /*if malloc failed*/
    if (checkMallocFailure3(A, D, W, ADataPoints, DDataPoints, WDataPoints) == 1){
        freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
        freeMatrix(inputMatrix, N);
        PyErr_NoMemory(); 
        return NULL;
    }
    similarityMatrix(A, N, d, inputMatrix); /*creating A*/
    diagonalMatrix(D, A, N); /*creating D*/
    raiseNegativeHalfDiag(D, D, N); /*D = D^(-0.5)*/
    normalizedSimilarityMatrix(W, A, D, N); /*creating W*/
    pyMatrix = matrixToPyMatrix(W, N, N); /*turning W into python matrix*/
    /*freeing memory*/
    freeData(A, D, W, NULL, NULL, ADataPoints, DDataPoints, WDataPoints, NULL, NULL);
    freeMatrix(inputMatrix, N);
    return pyMatrix;
}

static PyObject* symnmf(PyObject* self, PyObject *args){
    /**
     * Performs SymNMF factorization to update H matrix.
     * 
     * Python signature: symnmf(W, H)
     * @param self Module object
     * @param args Python tuple containing W and H matrices
     * @return Updated H matrix as Python list of lists
     */
    double **W = NULL, **H = NULL;
    PyObject *pyMatrixW, *pyMatrixH, *resMatrix;
    int N, k, res;
    
    /*converting python W and H to C matrices*/
    PyArg_ParseTuple(args, "OO", &pyMatrixW, &pyMatrixH);
    N = PyObject_Length(pyMatrixW); /*amount of rows in the matrix*/
    k = PyObject_Length(PyList_GetItem(pyMatrixH, 0)); /*length of the first row in matrix H*/
    W = PyMatrixToMatrix(pyMatrixW);
    H = PyMatrixToMatrix(pyMatrixH);
    if (H == NULL || W == NULL){
        freeMatrix(W, N);
        freeMatrix(H, N);
        PyErr_NoMemory(); 
        return NULL;
    }
    res = symnmfAlgorithm(W, H, N, k); /*running the symnmf algo, updates H*/
    if (res == 1){
        freeMatrix(W, N);
        freeMatrix(H, N);
        PyErr_NoMemory(); 
        return NULL;
    }
    resMatrix = matrixToPyMatrix(H, N, k); /*turning H into python matrix*/
    /*freeing memory*/
    freeMatrix(W, N);
    freeMatrix(H, N);
    return resMatrix;
}


static PyMethodDef symnmfMethods[] = {
    {"sym",                 
      (PyCFunction) sym,
      METH_VARARGS,           
      PyDoc_STR("input: N x d matrix of floats. output: returns the similarity matrix A")},
    {"ddg",                 
      (PyCFunction) ddg,
      METH_VARARGS,           
      PyDoc_STR("input: N x d matrix of floats. output: returns the diagonal degree matrix D")},
    {"norm",                 
      (PyCFunction) norm,
      METH_VARARGS,           
      PyDoc_STR("input: N x d matrix of floats. output: returns the normalized similarity matrix W")},
    {"symnmf",                 
      (PyCFunction) symnmf,
      METH_VARARGS,           
      PyDoc_STR("input: N x N normalized similarity matrix, N x k initial decomposition matrix H. output: returns the updated decomposition matrix H")},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    symnmfMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}