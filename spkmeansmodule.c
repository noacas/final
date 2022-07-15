#include <Python.h>
#include <ctype.h>
#include "useful.h"
#include "kmeans.h"
#include "stdlib.h"

void fill_data_list(int d, int n, double **data, PyObject * data_python);
PyObject* save_to_output(int d, int k, double **centroids);
static PyObject* fit(PyObject *self, PyObject *args);


static PyObject* fit(PyObject *self, PyObject *args)
{
    int k, max_iter, n, d, *cluster_size;
    double epsilon;
    PyObject * init_centroids_python, * data_points_python, *result;
    double **data_points, **centroids, **new_centroids, **centroids_result;
    if (!PyArg_ParseTuple(args, "iidiOiO", &k, &max_iter, &epsilon, &d, &init_centroids_python, &n, &data_points_python))
        return NULL;
    data_points = allocate_memory_array_of_points(d, n);
    centroids = allocate_memory_array_of_points(d, k);
    new_centroids = allocate_memory_array_of_points(d, k);
    cluster_size = allocate_memory_array_of_size(k);
    fill_data_list(d, n, data_points, data_points_python);
    fill_data_list(d, k, centroids, init_centroids_python);
    save_to_output(d, k, centroids);
    centroids_result = runAlg(d, k, n, max_iter, epsilon, data_points, centroids, new_centroids, cluster_size);
    free(data_points);
    free(cluster_size);
    result = save_to_output(d, k, centroids_result);
    free(centroids);
    free(new_centroids);
    return result;
}

void fill_data_list(int d, int n, double **data, PyObject * data_python) {
    int i, j;
    PyObject* inner_list, * py_float;
    for (i = 0 ; i < n ; i++) {
        inner_list = PyList_GetItem(data_python, i);
        for (j = 0; j < d; j++) {
            py_float = PyList_GetItem(inner_list, j);
            data[i][j] = PyFloat_AsDouble(py_float);
        }
    }
}

PyObject* save_to_output(int d, int k, double **centroids) {
    int i, j;
    PyObject *val, *inner_list, *result;
    result = PyList_New(k);
    for (i = 0 ; i < k ; i++) {
        inner_list = PyList_New(d);
        for (j = 0; j < d; j++) {
            val = PyFloat_FromDouble(centroids[i][j]);
            PyList_SetItem(inner_list, j, val);
        }
        PyList_SetItem(result, i, inner_list);
    }
    return result;
}

static PyMethodDef kmeansMethods[] = {
        {"fit",
                (PyCFunction) fit,
                     METH_VARARGS,
                        PyDoc_STR("kmeans alg")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        kmeansMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}


