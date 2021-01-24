#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "integrate.h"


PyObject* convert_complex_array(std::complex<double>* arr, int nf)
{
  PyObject* reals = PyList_New(nf);
  PyObject* imags = PyList_New(nf);
  PyObject* output = PyList_New(2);
  for (int i = 0; i < nf; i++) {
    PyObject* real = PyFloat_FromDouble(arr[i].real());
    PyObject* imag = PyFloat_FromDouble(arr[i].imag());
    PyList_SetItem(reals, i, real);
    PyList_SetItem(imags, i, imag);
  }
  PyList_SetItem(output, 0, reals);
  PyList_SetItem(output, 1, imags);
  return output;
}


static PyObject *Intg(PyObject *self, PyObject *args)
{
  int nf;
  int nl;
  PyObject* omegas;
  PyObject* ds;
  PyObject* kxs;
  PyObject* kys;
  PyObject* Cvs;
  double b;
  double lambda_i;
  double lambda_f;
  int N;
  char b_type;

  if (!PyArg_ParseTuple(args,
                        "iiOOOOOdddic",
                        &nf,
                        &nl,
                        &omegas,
                        &ds,
                        &kxs,
                        &kys,
                        &Cvs,
                        &b,
                        &lambda_i,
                        &lambda_f,
                        &N,
                        &b_type))
    return NULL;

  int n = PyObject_Length(ds);

  if (n < 0)
    return NULL;

  double omegas_[n];
  double ds_[n];
  double kxs_[n];
  double kys_[n];
  double Cvs_[n];

  for (int i = 0; i < n; i++) {
    PyObject* omega = PyList_GetItem(omegas, i);
    omegas_[i] = PyFloat_AsDouble(omega);
    
    PyObject* d = PyList_GetItem(ds, i);
    ds_[i] = PyFloat_AsDouble(d);
    
    PyObject* kx = PyList_GetItem(kxs, i);
    kxs_[i] = PyFloat_AsDouble(kx);

    PyObject* ky = PyList_GetItem(kys, i);
    kys_[i] = PyFloat_AsDouble(ky);

    PyObject* Cv = PyList_GetItem(Cvs, i);
    Cvs_[i] = PyFloat_AsDouble(Cv);
  }

  Sample s {nf, nl, omegas_, ds_, kxs_, kys_, Cvs_, b};
  Integral intg {lambda_i, lambda_f, N, s, b_type};

  std::complex<double>* result = intg.integrate();
  
  return convert_complex_array(result, nf);
}


static PyMethodDef IntgLib_FunctionsTable[] = {
  {"BTeq1", Intg, METH_VARARGS, "the integral term from Borca-Tasciuc Eq.(1)"},
  {NULL, NULL, 0, NULL}
};


static PyModuleDef IntgLib_Module = {
  PyModuleDef_HEAD_INIT,
  "BTeq1",
  "C++ function for fast implementation of Borca-Tasciuc Eq. (1)",
  -1,
  IntgLib_FunctionsTable
};


PyMODINIT_FUNC PyInit_BTeq1(void) {
  return PyModule_Create(&IntgLib_Module);
}
