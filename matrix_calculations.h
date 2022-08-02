//
// Created by Noa Caspi on 02/08/2022.
//

#ifndef FINAL_MATRIX_CALCULATIONS_H
#define FINAL_MATRIX_CALCULATIONS_H

double **createWeightedAdjacencyMatrix(int n, int d, double** data_points);
double **createDiagonalDegreeMatrix(int n, double ** weightedMatrix);
double ** getDiagonalMatrixPoweredByMinusHalf(int n, double ** matrix);
double **multiplyMatrix(int n, double **matrixA, double **matrixB);
double **multiply3Matrices(int n, double **matrixA, double **matrixB, double **matrixC);
double **subtractIbyMatrix(int n, double ** matrix);
double **createNormalizedGraphLaplacian(int n, double** diagonalDegreeMatrix, double ** weightedAdjacencyMatrix);
double **transposeMatrix(int n, double **matrix);
double **multipleFromBothSides(int n, double** matrixA, double ** matrixP);
int *getIndicesOfLargestAbsoluteValue(int n, double **matrix);
double *obtainCAndT(int n, double **matrix, int pivotI, int pivotJ);
double **getUnitMatrix(int n);
#endif //FINAL_MATRIX_CALCULATIONS_H
