#include "useful.h"
#include "math.h"
#include "stdlib.h"

const int MAX_ROTATION = 100;
const double EPSILON = 0.00001;

struct eigenValuesAndVectors {
    double * eigenvalues;
    double ** eigenvectors;
} eigenValuesAndVectors;


struct eigenValuesAndVectors jacobiAlg(double ** matrix, int n) {
    int i;
    double ** p, **nextMatrix, **rotationMatrix, **eigenVectorsMatrix=NULL;
    double off, prevOff = calcOff(matrix, n);
    for (i=0; i < MAX_ROTATION; i++) {
        rotationMatrix = getRotaionMatrix();
        if (eigenVectorsMatrix == NULL) {
            eigenVectorsMatrix = rotationMatrix;
        }
        nextMatrix = multipleFromBothSides();
        free(matrix);
        free(rotationMatrix);
        off = calcOff(nextMatrix);
        matrix = nextMatrix;
        if (checkConvergence(off, prevOff)) {
            break;
        }
        prevOff = off;
    }
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