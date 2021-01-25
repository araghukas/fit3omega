#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include "intg.h"


PyObject* convert_complex_array(complex<double>* arr, int nf)
{
  // returns list of two lists, real parts and complex parts
  PyObject* reals = PyList_New(nf);
  PyObject* imags = PyList_New(nf);
  PyObject* output = PyList_New(2);
  for (int i = 0; i < nf; i++) {
    PyObject* re = PyFloat_FromDouble(arr[i].real());
    PyObject* im = PyFloat_FromDouble(arr[i].imag());
    PyList_SetItem(reals, i, re);
    PyList_SetItem(imags, i, im);
  }
  PyList_SetItem(output, 0, reals);
  PyList_SetItem(output, 1, imags);
  return output;
}


IntegralTermBT_EQ1* INTG = NULL;
PyObject* Integrator(PyObject* self, PyObject* args)
{
  if (INTG == NULL)
    return NULL;

  PyObject* ds;
  PyObject* kxs;
  PyObject* kys;
  PyObject* Cvs;
  char b_type;

  if (!PyArg_ParseTuple(args,"OOOOc",&ds,&kxs,&kys,&Cvs,&b_type))
    return NULL;

  int nl = PyObject_Length(ds);

  if (nl < 0)
    return NULL;

  vector<double> d_(nl);
  vector<double> kx_(nl);
  vector<double> ky_(nl);
  vector<double> Cv_(nl);

  for (int i = 0; i < nl; i++) {
    PyObject* d = PyList_GetItem(ds, i);
    d_[i] = PyFloat_AsDouble(d);

    PyObject* kx = PyList_GetItem(kxs, i);
    kx_[i] = PyFloat_AsDouble(kx);

    PyObject* ky = PyList_GetItem(kys, i);
    ky_[i] = PyFloat_AsDouble(ky);

    PyObject* Cv = PyList_GetItem(Cvs, i);
    Cv_[i] = PyFloat_AsDouble(Cv);
  }

  complex<double>* result = INTG->integral(d_, kx_, ky_, Cv_, b_type);
  return convert_complex_array(result, INTG->get_nf());
}

static PyMethodDef integrator = {"integrator", Integrator, METH_VARARGS,
                                 "from-Python callable integrator function"};


PyObject* GetIntegrator(PyObject* self, PyObject* args)
// parametrize and return a function that Python can call
{
  PyObject* omegas_;
  double b;
  double lambda_i;
  double lambda_f;
  int N;

  if (!PyArg_ParseTuple(args,"Odddi",&omegas_,&b,&lambda_i,&lambda_f,&N))
    return NULL;

  int nf = PyObject_Length(omegas_);

  if (nf < 0)
    return NULL;

  vector<double> omegas(nf);
  for (int i = 0; i < nf; i++) {
    PyObject* omega = PyList_GetItem(omegas_, i);
    omegas[i] = PyFloat_AsDouble(omega);
  }

  if (INTG != NULL)
    delete INTG;

  INTG = new IntegralTermBT_EQ1{omegas, b, lambda_i, lambda_f, N};
  return PyCFunction_New(&integrator, NULL);
}


static PyMethodDef IntgLib_FunctionsTable[] = {
  {"get_integrator", GetIntegrator, METH_VARARGS,
   "get the from-Python callable integrator function"},
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
