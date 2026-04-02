# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "symnmf.h"


static double* PyListToRow(PyObject* lst){
    int i, n = PyObject_Length(lst);
    double* row = malloc(n * sizeof(double));
    if (row == NULL){
        HANDLE_ERROR();
    }
    for (i = 0; i < n; i++){
        row[i] = PyList_GetItem(lst, i);
    }
    return row;
}


static double** PyMatrixToMatrix(PyObject* matrix){
    int i, n = PyObject_Length(lst);
    double** returnMatrix = malloc(n * sizeof(double*));
    if (returnMatrix == NULL){
        HANDLE_ERROR();
    }
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
    *N = PyObject_Length(pyMatrix); /*amount of rows in the matrix*/
    *d = PyObject_Length(PyList_GetItem(pyMatrix, 0)); /*length of the first row in matrix*/
    inputMatrix = PyMatrixToMatrix(pyMatrix);
    return inputMatrix;
}


void freeMatrix(double** matrix, int rows){
    int i;
    for(i=0;i<rows;i++){
        free(matrix[i]);
    }
    free(matrix);
}


static PyObject* sym(PyObject* self, PyObject *args){
    double** inputMatrix, **A, *ADataPoints;
    PyObject* pyMatrix;
    int N, d;

    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);

    A = similarityMatrix(N, d, inputMatrix, &ADataPoints); /*creating A*/
    pyMatrix = matrixToPyMatrix(A, N, N); /*turning A into python matrix*/
    /*freeing memory*/
    free(A);
    free(ADataPoints);
    freeMatrix(inputMatrix, N);
    return pyMatrix;
}


static PyObject* ddg(PyObject* self, PyObject *args){
    double** inputMatrix, **A, *ADataPoints, **D, *DDataPoints;
    PyObject* pyMatrix;
    int N, d;
    
    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);

    A = similarityMatrix(N, d, inputMatrix, &ADataPoints); /*creating A*/
    D = diagonalMatrix(A, N, &DDataPoints); /*creating D*/
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
    double** inputMatrix, **A, *ADataPoints, **D, *DDataPoints, **DTag, *DTagDataPoints, **W, *WDataPoints;
    PyObject* pyMatrix;
    int N, d;
    
    inputMatrix = convertPyMatrixToCMatrix(args, &N, &d);

    A = similarityMatrix(N, d, inputMatrix, &ADataPoints); /*creating A*/
    D = diagonalMatrix(A, N, &DDataPoints); /*creating D*/
    DTag = raiseNegativeHalfDiag(D, N, &DTagDataPoints); /*D' = D^(-0.5)*/
    W = normalizedSimilarityMatrix(A, DTag, N, &WDataPoints); /*creating W*/
    pyMatrix = matrixToPyMatrix(W, N, N); /*turning W into python matrix*/
    /*freeing memory*/
    free(A);
    free(ADataPoints);
    free(D);
    free(DDataPoints);
    free(DTag);
    free(DTagDataPoints);
    free(W);
    free(WDataPoints);
    freeMatrix(inputMatrix, N);
    return pyMatrix;
}

static PyObject* symnmf(PyObject* self, PyObject *args){
    double **W, **H;
    PyObject* pyMatrixW, *pyMatrixH;
    int N, k;
    
    /*converting python W and H to C matrices*/
    PyArg_ParseTuple(args, "OO", &pyMatrixW, &pyMatrixH);
    N = PyObject_Length(pyMatrixW); /*amount of rows in the matrix*/
    k = PyObject_Length(PyList_GetItem(pyMatrixH, 0)); /*length of the first row in matrix H*/
    W = PyMatrixToMatrix(pyMatrixW);
    H = PyMatrixToMatrix(pyMatrixH);
    symnmfAlgorithm(W,H,N,k); /*running the symnmf algo, updates H*/
    
    pyMatrix = matrixToPyMatrix(H, N, N); /*turning H into python matrix*/
    /*freeing memory*/
    freeMatrix(W);
    freeMatrix(H);

    return pyMatrix;
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