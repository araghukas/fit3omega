#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define MAX_N_LAYERS 10
#define MAX_N_OMEGAS 150
#define N_XPTS       200


// SAMPLE PARAMETERS
const double* omegas;
double b;
double lambda_i;
double lambda_f;
int n_layers;
char b_type;

npy_intp n_omegas; // length of `omegas` array

int PARAMS_SET; // flag == "Python-side `set` method called successfully"

// INTERMEDIATE VALUES
double complex fB_[MAX_N_OMEGAS];
double complex fA_ [MAX_N_OMEGAS];
double complex phi_[MAX_N_OMEGAS];
double complex AB_next_[MAX_N_OMEGAS];
double complex kkB_[MAX_N_OMEGAS];
double complex tanh_term_[MAX_N_OMEGAS];
double complex f_prev_[MAX_N_OMEGAS];
double complex integrand_result_[MAX_N_OMEGAS];
double complex integral_result_[MAX_N_OMEGAS];
double lambdas_[N_XPTS];


// INTEGRATION PARAMETERS
double ds_ [MAX_N_LAYERS];
double kxs_ [MAX_N_LAYERS];
double kys_ [MAX_N_LAYERS];
double Cvs_ [MAX_N_LAYERS];


void make_logspace(double *arr, double min, double max, int size)
{
  for (int k = 0; k < size; k++)
    arr[k] = min * pow(max / min, ((double) k) / (size - 1));
}


double complex sinc_sq(double x)
{
  double sinc = sin(x) / x;
  return sinc*sinc + I*0.0;
}


double complex *fB(int i_layer, double lambda, double *kxs, double *kys, double *Cvs)
{
  /* Borca-Tasciuc Eq.(3) */
  for (int i = 0; i < n_omegas; i++) {
    double re = kxs[i_layer] / kys[i_layer] * lambda * lambda;
    double im = 2.0 * omegas[i] * Cvs[i_layer] / kys[i_layer];
    fB_[i] = csqrt(re + I*im);
  }
  return fB_;
}


double complex *phi(int i_layer, double lambda, double* ds, double* kxs, double* kys, double* Cvs)
{
  /* Borca-Tasciuc Eq.(4) */
  double complex *fB__ = fB(i_layer,lambda,kxs,kys,Cvs);
  for (int i = 0; i < n_omegas; i++)
    phi_[i] = ds[i_layer] * fB__[i];

  return phi_;
}


double complex *fA(int i_layer, double lambda, double* ds, double* kxs, double* kys, double* Cvs)
{
  /* Borca-Tasciuc Eq.(2) */
  if (i_layer == n_layers - 1) {
    switch (b_type)
    {
      case 's':
        for (int i = 0; i < n_omegas; i++)
          fA_[i] = -1.0 + I*0.0;
        break;
      case 'a':
      {
        double complex *phi__ = phi(i_layer,lambda,ds,kxs,kys,Cvs);
        for (int i = 0; i < n_omegas; i++)
          fA_[i] = -1.0 * ctanh(phi__[i]);
        break;
      }
      case 'i':
      {
        double complex *phi__ = phi(i_layer,lambda,ds,kxs,kys,Cvs);
        for (int i = 0; i < n_omegas; i++)
          fA_[i] = -1.0 / ctanh(phi__[i]);
        break;
      }
      default:
        return NULL;
    }
    return fA_;
  }

  // note: adding 1 to the index moves DOWNWARD in layers
  double complex *A_next = fA(i_layer+1,lambda,ds,kxs,kys,Cvs);
  double complex *B_next = fB(i_layer+1,lambda,kxs,kys,Cvs);

  for (int i = 0; i < n_omegas; i++)
    AB_next_[i] = A_next[i] * B_next[i];

  double complex *B_i = fB(i_layer,lambda,kxs,kys,Cvs);
  double k_i = kys[i_layer];
  double k_next = kys[i_layer+1];

  for (int i = 0; i < n_omegas; i++)
    kkB_[i] = k_next / (k_i * B_i[i]);

  double complex *phi__ = phi(i_layer,lambda,ds,kxs,kys,Cvs);

  for (int i = 0; i < n_omegas; i++)
    tanh_term_[i] = ctanh(phi__[i]);

  for (int i = 0; i < n_omegas; i++) {
    double complex ABkkB = AB_next_[i] * kkB_[i];
    double complex tanh_term = tanh_term_[i];
    fA_[i] = (ABkkB - tanh_term) / (1.0 - ABkkB * tanh_term);
  }
  return fA_;
}


double complex *integrand(double lambda, double *ds, double *kxs, double *kys, double *Cvs)
{
  /* Borca-Tasciuc Eq.(1) integrand */
  double complex *A_top = fA(0,lambda,ds,kxs,kys,Cvs);
  double complex *B_top = fB(0,lambda,kxs,kys,Cvs);
  double complex sinc_sq_ = sinc_sq(b * lambda);
  for (int i = 0; i < n_omegas; i++)
    integrand_result_[i] = sinc_sq_ / (A_top[i] * B_top[i]);
  return integrand_result_;
}


double complex *integral(double *ds, double *kxs, double *kys, double *Cvs)
{
  /* Borca-Tasciuc Eq.(1) integral */
  double complex *f0 = integrand(lambdas_[0],ds,kxs,kys,Cvs);
  for (int i = 0; i < n_omegas; i++) {
    integral_result_[i] = 0.0 + I*0.0;
    f_prev_[i] = f0[i];
  }

  for (int k = 1; k < N_XPTS; k++) {
    double complex *fk = integrand(lambdas_[k],ds,kxs,kys,Cvs);
    double dx = lambdas_[k] - lambdas_[k - 1];
    for (int i = 0; i < n_omegas; i++) {
      integral_result_[i] += (dx / 2.0 + 0.0*I) * (fk[i] + f_prev_[i]);
      f_prev_[i] = fk[i];
    }
  }
  return integral_result_;
}


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
THE WRAPPED INTEGRAL FUNCTION
----------------------------------------------------------------------------------------------------
*/


static PyObject *Integral(PyObject *self, PyObject *args)
{
  if (!PARAMS_SET)
    return NULL;

  PyObject *Py_ds;
  PyObject *Py_kys;
  PyObject *Py_ratio_xys;
  PyObject *Py_Cvs;

  if (!PyArg_ParseTuple(args,"OOOO",&Py_ds,&Py_kys,&Py_ratio_xys,&Py_Cvs))
    return NULL;

  if (PyObject_Length(Py_ds) != n_layers)
    return NULL;

  for (int j = 0; j < n_layers; j++) {
    PyObject *d = PyList_GetItem(Py_ds, j);
    ds_[j] = PyFloat_AsDouble(d);

    PyObject *ky = PyList_GetItem(Py_kys, j);
    kys_[j] = PyFloat_AsDouble(ky);

    PyObject *rat = PyList_GetItem(Py_ratio_xys, j);
    kxs_[j] = PyFloat_AsDouble(rat) * kys_[j];

    PyObject *Cv = PyList_GetItem(Py_Cvs, j);
    Cvs_[j] = PyFloat_AsDouble(Cv);
  }

  return as_complex_nparray(integral(ds_,kxs_,kys_,Cvs_), n_omegas);
}


/*
----------------------------------------------------------------------------------------------------
THE WRAPPED INTEGRAND FUNCTION
----------------------------------------------------------------------------------------------------
*/


static PyObject *Integrand(PyObject *self, PyObject *args)
{
  if (!PARAMS_SET)
    return NULL;

  double lambda;
  PyObject *Py_ds;
  PyObject *Py_kys;
  PyObject *Py_ratio_xys;
  PyObject *Py_Cvs;


  if (!PyArg_ParseTuple(args,"dOOOO",&lambda,&Py_ds,&Py_kys,&Py_ratio_xys,&Py_Cvs))
    return NULL;

  if ((int) PyObject_Length(Py_ds) != n_layers)
    return NULL;

  for (int j = 0; j < n_layers; j++) {
    PyObject *d = PyList_GetItem(Py_ds, j);
    ds_[j] = PyFloat_AsDouble(d);

    PyObject *ky = PyList_GetItem(Py_kys, j);
    kys_[j] = PyFloat_AsDouble(ky);

    PyObject *rat = PyList_GetItem(Py_ratio_xys, j);
    kxs_[j] = PyFloat_AsDouble(rat) * kys_[j];

    PyObject *Cv = PyList_GetItem(Py_Cvs, j);
    Cvs_[j] = PyFloat_AsDouble(Cv);
  }

  return as_complex_nparray(integrand(lambda,ds_,kxs_,kys_,Cvs_), n_omegas);
}


/*
----------------------------------------------------------------------------------------------------
SET STATIC VALUES
----------------------------------------------------------------------------------------------------
*/
static PyObject *Set(PyObject *self, PyObject *args)
{
  PyArrayObject *omegas_Py;
  if (!PyArg_ParseTuple(args,"O!dddic",&PyArray_Type,&omegas_Py,
                                       &b,&lambda_i,&lambda_f,&n_layers,&b_type)) // globals
    return NULL;

   if (n_layers <= 0 || n_layers > MAX_N_LAYERS)
    return NULL;

  n_omegas = PyArray_Size((PyObject *) omegas_Py);
  if (n_omegas <= 0 || n_omegas > MAX_N_OMEGAS)
    return NULL;

  npy_intp start_index = 0;
  omegas = PyArray_GetPtr(omegas_Py, &start_index); // global
  if (omegas == NULL)
    return NULL;
  Py_INCREF(omegas_Py);

  make_logspace(lambdas_, lambda_i, lambda_f, N_XPTS);
  PARAMS_SET = 1;
  Py_RETURN_NONE;
}




/*
----------------------------------------------------------------------------------------------------
MODULE DEFINITIONS
----------------------------------------------------------------------------------------------------
*/
static PyMethodDef Integrate_FunctionsTable[] = {
  {"set", Set, METH_VARARGS, "mandatory set static variable values"},
  {"integral", Integral, METH_VARARGS, "computes the integral term in Borca-Tascuic Eq. (1)"},
  {"integrand", Integrand, METH_VARARGS, "computes the integrand in Bora-Tascuic Eq. (1)"},
  {NULL, NULL, 0, NULL}
};


static PyModuleDef Integrate_Module = {
  PyModuleDef_HEAD_INIT,
  "intg",
  "C library for a fast implementation of Borca-Tascuic Eq. (1)",
  -1,
  Integrate_FunctionsTable
};


PyMODINIT_FUNC PyInit_intg(void) {
  import_array();
  PARAMS_SET = 0;
  return PyModule_Create(&Integrate_Module);
}
