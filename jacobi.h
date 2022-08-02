#ifndef FINAL_JACOBI_H
#define FINAL_JACOBI_H

struct eigenValueAndVector *jacobiAlg(double **, int);
int checkConvergence(double, double);
double calcOff(double **, int);
double **getRotationMatrix(int n, double **matrix);

#endif //FINAL_JACOBI_H
