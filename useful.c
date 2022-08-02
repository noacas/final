#include "math.h"
#include "stdlib.h"
#include <ctype.h>
#include <assert.h>
#define PY_SSIZE_T_CLEAN

double euclidean_norm_powered(int d, double *p1, double *p2) {
    int i;
    double result = 0;
    for (i=0; i < d ; i++) {
        result += pow(p1[i] - p2[i], 2);
    }
    return result;
}

double euclidean_norm(int d, double *p1, double *p2) {
    return sqrt(euclidean_norm_powered(d, p1, p2));
}

int* allocate_memory_array_of_size(int k) {
    int* a = calloc(k, sizeof(int));
    assert(a != NULL);
    return a;
}

double** allocate_memory_array_of_points(int d, int array_size) {
    double *p;
    double **a;
    int i;
    p = calloc(d * array_size, sizeof(double));
    a = calloc(array_size, sizeof(double *));
    for(i=0 ; i < array_size ; i++ )
        a[i] = p+ i * d;
    assert(a != NULL);
    assert(p != NULL);
    return a;
}

int sign(double num) {
    if (num >= 0) {
        return 1;
    }
    return 0;
}
