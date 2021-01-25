#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include "intg.h"


IntegralTermBT_EQ1* INTG = NULL;
IntegralTermBT_EQ1* get_integrator(vector<double>& o_,
                                   double b,
                                   double lambda_i,
                                   double lambda_f,
                                   int N)
{

  if (INTG == NULL)
    INTG = new IntegralTermBT_EQ1{o_,b,lambda_i,lambda_f,N};

  return INTG;
}


PyObject* convert_complex_array(complex<double>* arr, int nf)
{
  // returns list of two lists, real parts and complex parts
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


static PyObject* Integrate(PyObject* self, PyObject *args)
{
  // for args from Python side
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

  // assign above variables
  if (!PyArg_ParseTuple(args,
                        "OOOOOdddic",
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

  int nl = PyObject_Length(ds);
  int nf = PyObject_Length(omegas);

  if (nl < 0 || nf < 0)
    return NULL;

  vector<double> o_(nf);
  vector<double> d_(nl);
  vector<double> kx_(nl);
  vector<double> ky_(nl);
  vector<double> C_(nl);

  for (int i = 0; i < nf; i++) {
    PyObject* omega = PyList_GetItem(omegas, i);
    o_[i] = PyFloat_AsDouble(omega);
  }

  for (int i = 0; i < nl; i++) {
    PyObject* d = PyList_GetItem(ds, i);
    d_[i] = PyFloat_AsDouble(d);

    PyObject* kx = PyList_GetItem(kxs, i);
    kx_[i] = PyFloat_AsDouble(kx);

    PyObject* ky = PyList_GetItem(kys, i);
    ky_[i] = PyFloat_AsDouble(ky);

    PyObject* Cv = PyList_GetItem(Cvs, i);
    C_[i] = PyFloat_AsDouble(Cv);
  }

  // new or same integrator object
  IntegralTermBT_EQ1* intg = get_integrator(o_,b,lambda_i,lambda_f,N);
  complex<double>* result = intg->integral(d_, kx_, ky_, C_, b_type);
  return convert_complex_array(result, nf);
}


static PyMethodDef IntgLib_FunctionsTable[] = {
  {"integrate", Integrate, METH_VARARGS, "evaluate the integral term from Borca Eq. 1"},
  {NULL, NULL, 0, NULL}
};


static PyModuleDef IntgLib_Module = {
  PyModuleDef_HEAD_INIT,
  "intglib",
  "C++ function for fast implementation of Borca-Tasciuc Eq. (1)",
  -1,
  IntgLib_FunctionsTable
};


PyMODINIT_FUNC PyInit_intglib(void) {
  return PyModule_Create(&IntgLib_Module);
}
