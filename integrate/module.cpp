#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "intg.h"


struct IntegralArgs {
  const vector<double>& o_;
  const vector<double>& d_;
  const vector<double>& kx_;
  const vector<double>& ky_;
  const vector<double>& C_;
};


bool essEqual(double a, double b)
// returns true if doubles a and b are "essenetially" equal
{
  static double EPS = 1e-7;
  return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPS );
}


bool operator==(IntegralArgs& a1, IntegralArgs& a2)
{
  int nf1 = static_cast<int>(a1.o_.size());
  int nf2 = static_cast<int>(a2.o_.size());
  int nl1 = static_cast<int>(a1.d_.size());
  int nl2 = static_cast<int>(a2.d_.size());
  
  if (!(nf1 == nf2 && nl1 == nl2))
    return false;

  // NOTE: not checking omegas to equality

  for (int i = 0; i < nl1; i++)
    if (!essEqual(a1.d_[i],a2.d_[i]))
      return false;

  for (int i = 0; i < nl1; i++)
    if (!essEqual(a1.kx_[i],a2.kx_[i]))
      return false;

  for (int i = 0; i < nl1; i++)
    if (!essEqual(a1.ky_[i],a2.ky_[i]))
      return false;

  for (int i = 0; i < nl1; i++)
    if (!essEqual(a1.C_[i],a2.C_[i]))
      return false;

  return true;
}


IntegralTermBT_EQ1* INTG = NULL;
IntegralArgs* ARGS = NULL;
IntegralTermBT_EQ1* get_integrator(vector<double>& o_,
                                   vector<double>& d_,
                                   vector<double>& kx_,
                                   vector<double>& ky_,
                                   vector<double>& C_,
                                   double b,
                                   double lambda_i,
                                   double lambda_f,
                                   int N,
                                   char b_type)
{
  IntegralArgs* args = new IntegralArgs{o_, d_, kx_, ky_, C_};
  if (INTG == NULL || args != ARGS) {
    INTG = new IntegralTermBT_EQ1{o_,d_,kx_,ky_,C_,b,lambda_i,lambda_f,N,b_type};
    ARGS = args;
  }
  else
    delete args;

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

  IntegralTermBT_EQ1* intg = get_integrator(
    o_,d_,kx_,ky_,C_,b,lambda_i,lambda_f,N,b_type
  );
  
  return convert_complex_array(intg->integral(), nf);
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
