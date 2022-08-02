#include "useful.h"
#include "math.h"
#include "stdlib.h"
#include "matrix_calculations.h"

double **createWeightedAdjacencyMatrix(int n, int d, double** data_points) {
    double **matrix = allocate_memory_array_of_points(n, n);
    double val;
    int i, j;
    for (i=0; i < n; i++) {
        for (j=i+1; j < n; j++) {
            val = exp(-euclidean_norm(d, data_points[i], data_points[j])/2);
            matrix[i][j] = val;
            matrix[j][i] = val;
        }
    }
    return matrix;
}

double **createDiagonalDegreeMatrix(int n, double ** weightedMatrix) {
    double **matrix = allocate_memory_array_of_points(n, n);
    double val;
    int i, j;
    for (i=0; i < n; i++) {
        val = 0;
        for (j=0; j < n; j++) {
            val += weightedMatrix[i][j];
        }
        matrix[i][i] = val;
    }
    return matrix;
}

double ** getDiagonalMatrixPoweredByMinusHalf(int n, double ** matrix) {
    double **resultMatrix = allocate_memory_array_of_points(n, n);
    int i;
    for (i=0; i < n; i++) {
        resultMatrix[i][i] = 1 / sqrt(matrix[i][i]);
    }
    return resultMatrix;
}

double **multiplyMatrix(int n, double **matrixA, double **matrixB) {
    // matrixA and matrixB of size nxn
    double **matrix = allocate_memory_array_of_points(n, n);
    int i, j, k;
    double val;
    for (i=0; i < n ; i++) {
        for (j=0; j < n ; j++) {
            val = 0;
            for (k=0; k < n; k++) {
                val += matrixA[i][k] * matrixB[k][j];
            }
            matrix[i][j] = val;
        }
    }
    return matrix;
}


double **multiply3Matrices(int n, double **matrixA, double **matrixB, double **matrixC) {
    // all matrices of size nxn
    double ** tmp = multiplyMatrix(n, matrixA, matrixB);
    double ** result = multiplyMatrix(n, tmp, matrixC);
    free(tmp);
    return result;
}

double **subtractIbyMatrix(int n, double ** matrix) {
    double **result = allocate_memory_array_of_points(n, n);
    int i, j;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            if (i==j) {
                result[i][j] = 1;
            }
            result[i][j] -= matrix[i][j];
        }
    }
    return result;
}

double **getUnitMatrix(int n) {
    double **result = allocate_memory_array_of_points(n, n);
    int i, j;
    for (i=0; i < n; i++) {
        result[i][j] = 1;
    }
    return result;
}

double **createNormalizedGraphLaplacian(int n, double** diagonalDegreeMatrix, double ** weightedAdjacencyMatrix) {
    double ** d = getDiagonalMatrixPoweredByMinusHalf(n, diagonalDegreeMatrix);
    double ** tmp = multiply3Matrices(n, d, weightedAdjacencyMatrix, d);
    double ** result = subtractIbyMatrix(n, tmp);
    free(tmp);
    return result;
}

double **transposeMatrix(int n, double **matrix) {
    double **result = allocate_memory_array_of_points(n, n);
    int i, j;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
                result[i][j] = matrix[j][i];
        }
    }
    return result;
}

//return matrix P^tAP
double **multipleFromBothSides(int n, double** matrixA, double ** matrixP) {
    double **matrixPTransposed = transposeMatrix(n, matrixP);
    double **result = multiply3Matrices(n, matrixPTransposed, matrixA, matrixP);
    free(matrixPTransposed);
    return result;
}

int *getIndicesOfLargestAbsoluteValue(int n, double **matrix) {
    int indices[] = {0,0};
    int i, j;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            if (fabs(matrix[indices[0]][indices[1]]) <= fabs(matrix[i][j])) {
                indices[0] = i;
                indices[1] = j;
            }
        }
    }
    return indices;
}
