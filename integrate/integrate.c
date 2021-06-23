#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "integrate.h"
#include "borca_tasciuc.h"
#include "olson_graham_chen.h"
#include "exceptions.h"


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

  if (!PyArg_ParseTuple(args,"OOOO",&Py_ds,&Py_kys,&Py_ratio_xys,&Py_Cvs)) {
    PyErr_SetString(BT_IntegralError, ARGS_ERROR_MSG);
    return NULL;
   }

  if (PyObject_Length(Py_ds) != N_LAYERS) {
    PyErr_SetString(BT_IntegralError, LENGTH_ERROR_MSG);
    return NULL;
  }

  for (int j = 0; j < N_LAYERS; j++) {
    PyObject *d = PyList_GetItem(Py_ds, j);
    ds_[j] = PyFloat_AsDouble(d);

    PyObject *ky = PyList_GetItem(Py_kys, j);
    kys_[j] = PyFloat_AsDouble(ky);

    PyObject *rat = PyList_GetItem(Py_ratio_xys, j);
    psis_[j] = PyFloat_AsDouble(rat);

    PyObject *Cv = PyList_GetItem(Py_Cvs, j);
    Cvs_[j] = PyFloat_AsDouble(Cv);
  }

  return as_complex_nparray(bt_integral(), n_OMEGAS);
}


static PyObject *OGC_Integral(PyObject *self, PyObject * args)
{
  if (!OGC_PARAMS_SET)
    return NULL;

  PyObject *Py_ds;
  PyObject *Py_kys;
  PyObject *Py_ratio_xys;
  PyObject *Py_Cvs;
  PyObject *Py_Rcs;

  if (!PyArg_ParseTuple(args,"OOOOO",&Py_ds,&Py_kys,&Py_ratio_xys,&Py_Cvs,&Py_Rcs)) {
    PyErr_SetString(OGC_IntegralError, ARGS_ERROR_MSG);
    return NULL;
  }

  if (PyObject_Length(Py_ds) != N_LAYERS) {
    PyErr_SetString(OGC_IntegralError, LENGTH_ERROR_MSG);
    return NULL;
  }

  for (int j = 0; j < N_LAYERS; j++) {
    PyObject *d = PyList_GetItem(Py_ds, j);
    ds_[j] = PyFloat_AsDouble(d);

    PyObject *ky = PyList_GetItem(Py_kys, j);
    kys_[j] = PyFloat_AsDouble(ky);

    PyObject *rat = PyList_GetItem(Py_ratio_xys, j);
    psis_[j] = PyFloat_AsDouble(rat);

    PyObject *Cv = PyList_GetItem(Py_Cvs, j);
    Cvs_[j] = PyFloat_AsDouble(Cv);

    PyObject *Rc = PyList_GetItem(Py_Rcs, j);
    Rcs_[j] = PyFloat_AsDouble(Rc);
  }

  return as_complex_nparray(ogc_integral(), n_OMEGAS);
}


/*
----------------------------------------------------------------------------------------------------
INITIALIZER FUNCTIONS
----------------------------------------------------------------------------------------------------
*/
static PyObject *BT_Set(PyObject *self, PyObject *args)
{
  PyArrayObject *omegas_Py;
  double lambda_i_, lambda_f_;
  if (!PyArg_ParseTuple(args,"O!dddic",
      &PyArray_Type,&omegas_Py,&HALF_WIDTH,&lambda_i_,&lambda_f_,&N_LAYERS,&boundary_type_)) {
    PyErr_SetString(BT_SetError, ARGS_ERROR_MSG);
    return NULL;
  }

  if (N_LAYERS <= 0 || N_LAYERS > MAX_N_LAYERS) {
    PyErr_SetString(N_LAYERS_Error, N_LAYERS_Error_MSG);
    return NULL;
  }

  n_OMEGAS = (int) PyArray_Size((PyObject *) omegas_Py);
  if (n_OMEGAS <= 0 || n_OMEGAS > MAX_N_OMEGAS) {
    PyErr_SetString(n_OMEGAS_Error, n_OMEGAS_Error_MSG);
    return NULL;
  }

  npy_intp start_index = 0;
  OMEGAS = PyArray_GetPtr(omegas_Py, &start_index); // global
  if (OMEGAS == NULL) {
    PyErr_SetString(NullOmegasError, NullOmegasError_MSG);
    return NULL;
  }
  Py_INCREF(omegas_Py);

  make_logspace(LAMBDAS, lambda_i_, lambda_f_, N_XPTS);
  BT_PARAMS_SET = 1;
  Py_RETURN_NONE;
}


static PyObject *OGC_Set(PyObject *self, PyObject *args)
{
  PyArrayObject *omegas_Py;
  double chi_i_, chi_f_;
  if (!PyArg_ParseTuple(args,"O!dddi",
      &PyArray_Type,&omegas_Py,&HALF_WIDTH,&chi_i_,&chi_f_,&N_LAYERS)) {
    PyErr_SetString(OGC_SetError, ARGS_ERROR_MSG);
    return NULL;
  }

  if (N_LAYERS <= 0 || N_LAYERS > MAX_N_LAYERS) {
    PyErr_SetString(N_LAYERS_Error, N_LAYERS_Error_MSG);
    return NULL;
  }
  n_OMEGAS = PyArray_Size((PyObject *) omegas_Py);
  if (n_OMEGAS <= 0 || n_OMEGAS > MAX_N_OMEGAS) {
    PyErr_SetString(n_OMEGAS_Error, n_OMEGAS_Error_MSG);
    return NULL;
  }

  npy_intp start_index = 0;
  OMEGAS = PyArray_GetPtr(omegas_Py, &start_index); // global
  if (OMEGAS == NULL) {
    PyErr_SetString(NullOmegasError, NullOmegasError_MSG);
    return NULL;
  }
  Py_INCREF(omegas_Py);

  make_logspace(CHIS, chi_i_, chi_f_, N_XPTS);
  OGC_PARAMS_SET = 1;
  Py_RETURN_NONE;
}


/*
----------------------------------------------------------------------------------------------------
MODULE DEFINITIONS
----------------------------------------------------------------------------------------------------
*/
static PyMethodDef Integrate_FunctionsTable[] = {
  {"bt_set", BT_Set, METH_VARARGS, "mandatory initializer method"},
  {"ogc_set", OGC_Set, METH_VARARGS, "mandatory initializer method"},
  {"bt_integral", BT_Integral, METH_VARARGS, "computes the integral term in Borca-Tascuic Eq. (1)"},
  {"ogc_integral", OGC_Integral, METH_VARARGS, "computes the entire integral in OGC Eq. (4)"},
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
  init_exceptions();

  BT_PARAMS_SET = 0;
  OGC_PARAMS_SET = 0;
  return PyModule_Create(&Integrate_Module);
};
