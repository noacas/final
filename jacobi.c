#include "useful.h"
#include "math.h"
#include "stdlib.h"
#include "jacobi.h"
#include "matrix_calculations.h"

const int MAX_ROTATION = 100;
const double EPSILON = 0.00001;

struct eigenValueAndVector {
    double eigenvalues;
    double * eigenvectors;
} eigenValuesAndVectors;


struct eigenValueAndVector* jacobiAlg(double ** matrix, int n) {
    struct eigenValuesAndVectors* result = calloc(n, sizeof(eigenValuesAndVectors));
    int i, j;
    double ** p, **nextMatrix, **rotationMatrix, **eigenVectorsMatrix=NULL, **tmp, off, prevOff = calcOff(matrix, n), *vector;
    for (i=0; i < MAX_ROTATION; i++) {
        rotationMatrix = getRotationMatrix(n, matrix);
        if (eigenVectorsMatrix == NULL) {
            eigenVectorsMatrix = rotationMatrix;
        }
        else {
            tmp = eigenVectorsMatrix;
            eigenVectorsMatrix = multipleFromBothSides(n, eigenVectorsMatrix, rotationMatrix);
            free(tmp);
        }
        nextMatrix = multipleFromBothSides(n, matrix, rotationMatrix);
        free(matrix);
        free(rotationMatrix);
        off = calcOff(nextMatrix, n);
        matrix = nextMatrix;
        if (checkConvergence(off, prevOff)) {
            break;
        }
        prevOff = off;
    }
    for (i=0; i<n; i++) {
        result[i].eigenvalues = matrix[i][i];
        vector = calloc(n, sizeof(double));
        for (j=0; j<n; j++) {
            vector[j] = eigenVectorsMatrix[j][i];
        }
        result[i].eigenvector = vector;
    }
    return result;
}

int checkConvergence(double prevOff, double off) {
    return (prevOff-off <= EPSILON);
}

double calcOff(double ** matrix, int n) {
    double sum = 0;
    int i,j;
    for (i=0; i < n; i++) {
        for (j=0; j < n; j++) {
            if (i!=j) {
                sum += pow(matrix[i][j], 2);
            }
        }
    }
    return sum;
}

double **getRotationMatrix(int n, double **matrix) {
    int i, j, * indices = getIndicesOfLargestAbsoluteValue(n, matrix);
    i = indices[0];
    j = indices[1];
    free(indices);
    double c, t, s, ** result, *cAndT = obtainCAndT(n, matrix, i, j);
    c = cAndT[0];
    t = cAndT[1];
    free(cAndT);
    s = t * c;
    result = getUnitMatrix(n);
    result[i][i] = c;
    result[j][j] = c;
    result[i][j] = sign(j-i) * s;
    result[j][i] = sign(i-j) * s;
    return result;
}

double *obtainCAndT(int n, double **matrix, int pivotI, int pivotJ) {
    double teta, t, c, s, result[2];
    teta = (matrix[pivotJ][pivotJ] - matrix[pivotI][pivotI]) / (2*matrix[pivotI][pivotJ]);
    t = sign(teta) / (fabs(teta) + sqrt(pow(teta, 2) + 1));
    c = 1 / sqrt(pow(t, 2) + 1);
    result[0] = c;
    result[1] = t;
    return result;
}
