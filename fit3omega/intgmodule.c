/*
(Nov. 16, 2019)

Python wrapper for Borca eq. (1) implemented in C.

This is probably the simplest possible working solution. I would strongly advise
taking a look at the site linked below if, like me, you are not too familiar
with C.

<https://stackabuse.com/enhancing-python-with-custom-c-extensions/>

There are many tutorials like the above out there, but this is the only one
(that I could find) that implements a simple way of digesting a Python list into
a C array.

*/

#include <Python.h>
#include "intglib.h"


/* WRAPPED INTEGRAND FUNCTIONS */
// semi-infinite case
static PyObject *F_s(PyObject *self, PyObject *args) {
  double L;
  double b;
  double omega;
  PyObject *kxs;
  PyObject *kys;
  PyObject *Cvs;
  PyObject *ds;

  if (!PyArg_ParseTuple(args, "dddOOOO",
    &L, &b, &omega, &kxs, &kys, &Cvs, &ds)) {
    return NULL;
    }

  int n = PyObject_Length(ds);
  if (n < 0) {
    return NULL;
  }

  double c_ds[n];
  double c_kxs[n];
  double c_kys[n];
  double c_Cvs[n];

  for (int i = 0; i < n; i++) {
    PyObject *kx_item = PyList_GetItem(kxs, i);
    double kx = PyFloat_AsDouble(kx_item);
    c_kxs[i] = kx;

    PyObject *ky_item = PyList_GetItem(kys, i);
    double ky = PyFloat_AsDouble(ky_item);
    c_kys[i] = ky;

    PyObject *Cv_item = PyList_GetItem(Cvs, i);
    double Cv = PyFloat_AsDouble(Cv_item);
    c_Cvs[i] = Cv;

    PyObject *d_item = PyList_GetItem(ds, i);
    double d = PyFloat_AsDouble(d_item);
    c_ds[i] = d;
  }

  double result = f_s(n, L, b, omega, c_kxs, c_kys, c_Cvs, c_ds);
  return Py_BuildValue("d", result);
}


// adiabatic case
static PyObject *F_a(PyObject *self, PyObject *args) {
  double L;
  double b;
  double omega;
  PyObject *kxs;
  PyObject *kys;
  PyObject *Cvs;
  PyObject *ds;

  if (!PyArg_ParseTuple(args, "dddOOOO",
    &L, &b, &omega, &kxs, &kys, &Cvs, &ds)) {
    return NULL;
    }

  int n = PyObject_Length(ds);
  if (n < 0) {
    return NULL;
  }

  double c_ds[n];
  double c_kxs[n];
  double c_kys[n];
  double c_Cvs[n];

  for (int i = 0; i < n; i++) {
    PyObject *kx_item = PyList_GetItem(kxs, i);
    double kx = PyFloat_AsDouble(kx_item);
    c_kxs[i] = kx;

    PyObject *ky_item = PyList_GetItem(kys, i);
    double ky = PyFloat_AsDouble(ky_item);
    c_kys[i] = ky;

    PyObject *Cv_item = PyList_GetItem(Cvs, i);
    double Cv = PyFloat_AsDouble(Cv_item);
    c_Cvs[i] = Cv;

    PyObject *d_item = PyList_GetItem(ds, i);
    double d = PyFloat_AsDouble(d_item);
    c_ds[i] = d;
  }

  double result = f_a(n, L, b, omega, c_kxs, c_kys, c_Cvs, c_ds);
  return Py_BuildValue("d", result);
}


// isothermal case
static PyObject *F_o(PyObject *self, PyObject *args) {
  double L;
  double b;
  double omega;
  PyObject *kxs;
  PyObject *kys;
  PyObject *Cvs;
  PyObject *ds;

  if (!PyArg_ParseTuple(args, "dddOOOO",
    &L, &b, &omega, &kxs, &kys, &Cvs, &ds)) {
    return NULL;
    }

  int n = PyObject_Length(ds);
  if (n < 0) {
    return NULL;
  }

  double c_ds[n];
  double c_kxs[n];
  double c_kys[n];
  double c_Cvs[n];

  for (int i = 0; i < n; i++) {
    PyObject *kx_item = PyList_GetItem(kxs, i);
    double kx = PyFloat_AsDouble(kx_item);
    c_kxs[i] = kx;

    PyObject *ky_item = PyList_GetItem(kys, i);
    double ky = PyFloat_AsDouble(ky_item);
    c_kys[i] = ky;

    PyObject *Cv_item = PyList_GetItem(Cvs, i);
    double Cv = PyFloat_AsDouble(Cv_item);
    c_Cvs[i] = Cv;

    PyObject *d_item = PyList_GetItem(ds, i);
    double d = PyFloat_AsDouble(d_item);
    c_ds[i] = d;
  }

  double result = f_o(n, L, b, omega, c_kxs, c_kys, c_Cvs, c_ds);
  return Py_BuildValue("d", result);
}



/* WRAPPED INTEGRAL FUNCTIONS */
// semi-infinite case
static PyObject *Intg_a(PyObject *self, PyObject *args) {
  double xi;
  double xf;
  int N;
  double b;
  double omega;
  PyObject *ds;
  PyObject *kxs;
  PyObject *kys;
  PyObject *Cvs;

  if (!PyArg_ParseTuple(args, "ddiddOOOO",
    &xi, &xf, &N, &b, &omega, &ds, &kxs, &kys, &Cvs)) {
    return NULL;
  }

  int n = PyObject_Length(ds);
  if (n < 0) {
    return NULL;
  }

  double c_ds[n];
  double c_kxs[n];
  double c_kys[n];
  double c_Cvs[n];

  for (int i = 0; i < n; i++) {
    PyObject *d_item = PyList_GetItem(ds, i);
    double d = PyFloat_AsDouble(d_item);
    c_ds[i] = d;

    PyObject *kx_item = PyList_GetItem(kxs, i);
    double kx = PyFloat_AsDouble(kx_item);
    c_kxs[i] = kx;

    PyObject *ky_item = PyList_GetItem(kys, i);
    double ky = PyFloat_AsDouble(ky_item);
    c_kys[i] = ky;

    PyObject *Cv_item = PyList_GetItem(Cvs, i);
    double Cv = PyFloat_AsDouble(Cv_item);
    c_Cvs[i] = Cv;
  }

  double result;
  result = integrate_f_a(n, xi, xf, N, b, omega, c_ds, c_kxs, c_kys, c_Cvs);

  return Py_BuildValue("d", result);
}


// adiabatic case
static PyObject *Intg_s(PyObject *self, PyObject *args) {
  double xi;
  double xf;
  int N;
  double b;
  double omega;
  PyObject *ds;
  PyObject *kxs;
  PyObject *kys;
  PyObject *Cvs;

  if (!PyArg_ParseTuple(args, "ddiddOOOO",
    &xi, &xf, &N, &b, &omega, &ds, &kxs, &kys, &Cvs)) {
    return NULL;
  }

  int n = PyObject_Length(ds);
  if (n < 0) {
    return NULL;
  }

  double c_ds[n];
  double c_kxs[n];
  double c_kys[n];
  double c_Cvs[n];

  for (int i = 0; i < n; i++) {
    PyObject *d_item = PyList_GetItem(ds, i);
    double d = PyFloat_AsDouble(d_item);
    c_ds[i] = d;

    PyObject *kx_item = PyList_GetItem(kxs, i);
    double kx = PyFloat_AsDouble(kx_item);
    c_kxs[i] = kx;

    PyObject *ky_item = PyList_GetItem(kys, i);
    double ky = PyFloat_AsDouble(ky_item);
    c_kys[i] = ky;

    PyObject *Cv_item = PyList_GetItem(Cvs, i);
    double Cv = PyFloat_AsDouble(Cv_item);
    c_Cvs[i] = Cv;
  }

  double result;
  result = integrate_f_s(n, xi, xf, N, b, omega, c_ds, c_kxs, c_kys, c_Cvs);

  return Py_BuildValue("d", result);
}


// isothermal case
static PyObject *Intg_o(PyObject *self, PyObject *args) {
  double xi;
  double xf;
  int N;
  double b;
  double omega;
  PyObject *ds;
  PyObject *kxs;
  PyObject *kys;
  PyObject *Cvs;

  if (!PyArg_ParseTuple(args, "ddiddOOOO",
    &xi, &xf, &N, &b, &omega, &ds, &kxs, &kys, &Cvs)) {
    return NULL;
  }

  int n = PyObject_Length(ds);
  if (n < 0) {
    return NULL;
  }

  double c_ds[n];
  double c_kxs[n];
  double c_kys[n];
  double c_Cvs[n];

  for (int i = 0; i < n; i++) {
    PyObject *d_item = PyList_GetItem(ds, i);
    double d = PyFloat_AsDouble(d_item);
    c_ds[i] = d;

    PyObject *kx_item = PyList_GetItem(kxs, i);
    double kx = PyFloat_AsDouble(kx_item);
    c_kxs[i] = kx;

    PyObject *ky_item = PyList_GetItem(kys, i);
    double ky = PyFloat_AsDouble(ky_item);
    c_kys[i] = ky;

    PyObject *Cv_item = PyList_GetItem(Cvs, i);
    double Cv = PyFloat_AsDouble(Cv_item);
    c_Cvs[i] = Cv;
  }

  double result;
  result = integrate_f_o(n, xi, xf, N, b, omega, c_ds, c_kxs, c_kys, c_Cvs);

  return Py_BuildValue("d", result);
}


static PyMethodDef IntgLib_FunctionsTable[] = {
  {"intg_a", Intg_a, METH_VARARGS, "adiabatic boundary condition (integral)"},
  {"intg_s", Intg_s, METH_VARARGS, "semi-infinite boundary condition (integral)"},
  {"intg_o", Intg_o, METH_VARARGS, "isothermal boundary condition (integral)"},
  {"f_a", F_a, METH_VARARGS, "adiabatic boudnary condition (integrand)"},
  {"f_s", F_s, METH_VARARGS, "semi-infinite boundary condition (integrand)"},
  {"f_o", F_o, METH_VARARGS, "isothermal boundary condition (integrand)"},
  {NULL, NULL, 0, NULL}
};


static PyModuleDef IntgLib_Module = {
  PyModuleDef_HEAD_INIT,
  "intglib",
  "C library for fast implementation of Borca eq. (1) for general anisotropic DeltaT",
  -1,
  IntgLib_FunctionsTable
};

PyMODINIT_FUNC PyInit_intglib(void) {
  return PyModule_Create(&IntgLib_Module);
}
