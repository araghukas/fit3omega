#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "integrate.h"


// flags == "Python-side `set` method called successfully; parameters received."
int BT_PARAMS_SET;
int OGC_PARAMS_SET;


// =================================================================================================
//
//
// PYTHON BOILERPLATE BELOW
//
//
// =================================================================================================


static PyObject *as_complex_nparray(double complex *arr, int size)
{
  const npy_intp dims = size;
  return PyArray_SimpleNewFromData(1, &dims, NPY_COMPLEX128, arr);
}


/*
----------------------------------------------------------------------------------------------------
THE WRAPPED INTEGRAL FUNCTIONS
----------------------------------------------------------------------------------------------------
*/


static PyObject *BT_Integral(PyObject *self, PyObject *args)
{
  if (!BT_PARAMS_SET)
    return NULL;

  PyObject *Py_ds;
  PyObject *Py_kys;
  PyObject *Py_ratio_xys;
  PyObject *Py_Cvs;

  if (!PyArg_ParseTuple(args,"OOOO",&Py_ds,&Py_kys,&Py_ratio_xys,&Py_Cvs))
    return NULL;

  if (PyObject_Length(Py_ds) != n_layers_)
    return NULL;

  for (int j = 0; j < n_layers_; j++) {
    PyObject *d = PyList_GetItem(Py_ds, j);
    ds_[j] = PyFloat_AsDouble(d);

    PyObject *ky = PyList_GetItem(Py_kys, j);
    kys_[j] = PyFloat_AsDouble(ky);

    PyObject *rat = PyList_GetItem(Py_ratio_xys, j);
    kxs_[j] = PyFloat_AsDouble(rat) * kys_[j];

    PyObject *Cv = PyList_GetItem(Py_Cvs, j);
    Cvs_[j] = PyFloat_AsDouble(Cv);
  }

  return as_complex_nparray(bt_integral(ds_,kxs_,kys_,Cvs_), n_OMEGAS);
}


/*
----------------------------------------------------------------------------------------------------
SET STATIC VALUES
----------------------------------------------------------------------------------------------------
*/
static PyObject *BT_Set(PyObject *self, PyObject *args)
{
  PyArrayObject *omegas_Py;
  if (!PyArg_ParseTuple(args,"O!dddic",&PyArray_Type,&omegas_Py,
                                       &b_,&lambda_i_,&lambda_f_,&n_layers_,&boundary_type_)) // globals
    return NULL;

   if (n_layers_ <= 0 || n_layers_ > MAX_N_LAYERS)
    return NULL;

  n_OMEGAS = (int) PyArray_Size((PyObject *) omegas_Py);
  if (n_OMEGAS <= 0 || n_OMEGAS > MAX_N_OMEGAS)
    return NULL;

  npy_intp start_index = 0;
  OMEGAS = PyArray_GetPtr(omegas_Py, &start_index); // global
  if (OMEGAS == NULL)
    return NULL;
  Py_INCREF(omegas_Py);

  make_logspace(LAMBDAS, lambda_i_, lambda_f_, N_LAMBDAS);
  BT_PARAMS_SET = 1;
  Py_RETURN_NONE;
}




/*
----------------------------------------------------------------------------------------------------
MODULE DEFINITIONS
----------------------------------------------------------------------------------------------------
*/
static PyMethodDef Integrate_FunctionsTable[] = {
  {"bt_set", BT_Set, METH_VARARGS, "mandatory set static variable values"},
  {"bt_integral", BT_Integral, METH_VARARGS, "computes the integral term in Borca-Tascuic Eq. (1)"},
  {NULL, NULL, 0, NULL}
};


static PyModuleDef Integrate_Module = {
  PyModuleDef_HEAD_INIT,
  "integrate",
  "C library for a fast implementation of 3ω data fitting methods",
  -1,
  Integrate_FunctionsTable
};


PyMODINIT_FUNC PyInit_integrate(void) {
  import_array();
  BT_PARAMS_SET = 0;
  OGC_PARAMS_SET = 0;
  return PyModule_Create(&Integrate_Module);
}
