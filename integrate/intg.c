
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

const int MAX_n_layers = 10;


// SAMPLE PARAMETERS
double *omegas;
double b;
double lambda_i;
double lambda_f;
int n_pts;
int n_layers;
char b_type;
int n_omegas;


// INTERMEDIATE VALUES
double *lambdas_;
double complex *fB_;
double complex *fA_;
double complex *phi_;
double complex *AB_next_;
double complex *kkB_;
double complex *tanh_term_;
double complex *integrand_result_;
double complex *integral_result_;
double complex *f_prev_;


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
  for (int i = 0; i < n_omegas; i++) {
    double re = kxs[i_layer] / kys[i_layer] * lambda * lambda;
    double im = 2.0 * omegas[i] * Cvs[i_layer] / kys[i_layer];
    fB_[i] = csqrt(re + I*im);
  }
  return fB_;
}


double complex *phi(int i_layer, double lambda, double* ds, double* kxs, double* kys, double* Cvs)
{
  double complex *fB__ = fB(i_layer,lambda,kxs,kys,Cvs);
  for (int i = 0; i < n_omegas; i++)
    phi_[i] = ((double complex) ds[i_layer]) * fB__[i];

  return phi_;
}


double complex *fA(int i_layer, double lambda, double* ds, double* kxs, double* kys, double* Cvs)
{
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
          fA_[i] = (-1.0 + I*0.0) * ctanh(phi__[i]);
        break;
      }
      case 'i':
      {
        double complex *phi__ = phi(i_layer,lambda,ds,kxs,kys,Cvs);
        for (int i = 0; i < n_omegas; i++)
          fA_[i] = (-1.0 + I*0.0) / ctanh(phi__[i]);
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
    AB_next_[i] = A_next[i] + B_next[i];

  double complex *B_i = fB(i_layer,lambda,kxs,kys,Cvs);
  double k_i = kys[i_layer];
  double k_next = kys[i_layer+1];

  for (int i = 0; i < n_omegas; i++)
    kkB_[i] = (k_next / (k_i * B_i[i]) + I*0.0);

  double complex *phi__ = phi(i_layer,lambda,ds,kxs,kys,Cvs);

  for (int i = 0; i < n_omegas; i++)
    tanh_term_[i] = tanh(phi__[i]);

  for (int i = 0; i < n_omegas; i++) {
    double complex ABkkB = AB_next_[i] * kkB_[i];
    double complex tanh_term = tanh_term_[i];
    fA_[i] = (ABkkB - tanh_term) / ((1.0 + I*0.0) - ABkkB * tanh_term);
  }
  return fA_;
}


double complex *integrand(double lambda, double *ds, double *kxs, double *kys, double *Cvs)
{
  double complex *A_top = fA(0,lambda,ds,kxs,kys,Cvs);
  double complex *B_top = fB(0,lambda,kxs,kys,Cvs);
  double complex sinc_sq_ = sinc_sq(b * lambda);
  for (int i = 0; i < n_omegas; i++)
    integrand_result_[i] = sinc_sq_ / (A_top[i] * B_top[i]);
  return integrand_result_;
}


double complex *integral(double *ds, double *kxs, double *kys, double *Cvs)
{
  double complex *f0 = integrand(lambdas_[0],ds,kxs,kys,Cvs);
  for (int i = 0; i < n_omegas; i++) {
    integral_result_[i] = 0.0 + I*0.0;
    f_prev_[i] = f0[i];
  }

  for (int k = 1; k < n_pts; k++) {
    double complex *fk = integrand(lambdas_[k],ds,kxs,kys,Cvs);
    double dx = lambdas_[k] - lambdas_[k - 1];
    for (int i = 0; i < n_omegas; i++) {
      integral_result_[i] += (dx / 2.0 + 0.0*I) * (fk[i] + f_prev_[i]);
      f_prev_[i] = fk[i];
    }
  }
  return integral_result_;
}


void start_module(double b_,double lambda_i_,double lambda_f_,int n_pts_,int n_layers_,char b_type_,int n_omegas_)
{
  // TODO: check allocations!!

  b = b_;
  lambda_i = lambda_i_;
  lambda_f = lambda_f_;
  n_pts = n_pts_;
  n_layers = n_layers_;
  b_type = b_type_;
  n_omegas = n_omegas_;

  // omegas = (double *) malloc(n_omegas * sizeof(double));
  lambdas_ = (double *) PyMem_RawMalloc(n_omegas * sizeof(double));
  fB_ = (double complex *) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  fA_ = (double complex *) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  phi_ = (double complex *) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  AB_next_ = (double complex*) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  kkB_ = (double complex*) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  tanh_term_ = (double complex*) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  integrand_result_ = (double complex*) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  integral_result_ = (double complex*) PyMem_RawMalloc(n_omegas * sizeof(double complex));
  f_prev_ = (double complex*) PyMem_RawMalloc(n_omegas * sizeof(double complex));

  make_logspace(lambdas_,lambda_i,lambda_f,n_pts);
  // populate omegas and init sample inside wrapped Py func
}


void stop_module()
{
  PyMem_RawFree(lambdas_);
  PyMem_RawFree(fB_);
  PyMem_RawFree(fA_);
  PyMem_RawFree(phi_);
  PyMem_RawFree(AB_next_);
  PyMem_RawFree(kkB_);
  PyMem_RawFree(tanh_term_);
  PyMem_RawFree(integrand_result_);
  PyMem_RawFree(integral_result_);
  PyMem_RawFree(f_prev_);
  PyMem_RawFree(omegas);
}


// =================================================================================================
//
//
// PYTHON BOILERPLACE BELOW
//
//
// =================================================================================================

// ARGUMENT ARRAYS
double ds_[MAX_n_layers];
double kxs_[MAX_n_layers];
double kys_[MAX_n_layers];
double Cvs_[MAX_n_layers];


static PyObject *convert_complex_array(double complex *arr, int size)
{
  PyObject *reals = PyList_New(n_omegas);
  PyObject *imags = PyList_New(n_omegas);
  for (int i = 0; i < n_omegas; i++) {
    double complex z = arr[i];
    PyObject *re = PyFloat_FromDouble(creal(z));
    PyObject *im = PyFloat_FromDouble(cimag(z));
    PyList_SetItem(reals, i, re);
    PyList_SetItem(imags, i, im);
  }
  PyObject *output = PyList_New(2);
  PyList_SetItem(output, 0, reals);
  PyList_SetItem(output, 1, imags);
  return output;
}


/*
----------------------------------------------------------------------------------------------------
MODULE INITIALIZING
----------------------------------------------------------------------------------------------------
*/
static PyObject *StartModule(PyObject *self, PyObject *args)
{
  PyObject *Py_omegas;
  double b;
  double lambda_i;
  double lambda_f;
  int n_pts;
  int n_layers;
  char b_type;

  if (!PyArg_ParseTuple(args,"Odddiic",&Py_omegas,&b,&lambda_i,&lambda_f,&n_pts,&n_layers,&b_type))
    return NULL;

  // load provided omegas into global `omegas` array
  int n_omegas = PyObject_Length(Py_omegas);
  if (n_omegas < 0)
    return NULL;

  omegas = (double *) PyMem_RawMalloc(n_omegas * sizeof(double));

  for (int i = 0; i < n_omegas; i++) {
    PyObject *Py_omega = PyList_GetItem(Py_omegas, i);
    omegas[i] = PyFloat_AsDouble(Py_omega);
  }


  start_module(b,lambda_i,lambda_f,n_pts,n_layers,b_type,n_omegas);
  Py_RETURN_NONE;
}


static PyObject *StopModule(PyObject *self, PyObject *args)
{
  stop_module();
  Py_RETURN_NONE;
}


/*
----------------------------------------------------------------------------------------------------
THE WRAPPED INTEGRAL FUNCTION
----------------------------------------------------------------------------------------------------
*/


static PyObject *Integral(PyObject *self, PyObject *args)
{
  PyObject *Py_ds;
  PyObject *Py_kxs;
  PyObject *Py_kys;
  PyObject *Py_Cvs;

  if (!PyArg_ParseTuple(args,"OOOO",&Py_ds,&Py_kxs,&Py_kys,&Py_Cvs))
    return NULL;

  int n_layers = PyObject_Length(Py_ds);
  if (n_layers < 0 || n_layers > MAX_n_layers)
    return NULL;

  for (int j = 0; j < n_layers; j++) {
    PyObject *d = PyList_GetItem(Py_ds, j);
    ds_[j] = PyFloat_AsDouble(d);

    PyObject *kx = PyList_GetItem(Py_kxs, j);
    kxs_[j] = PyFloat_AsDouble(kx);

    PyObject *ky = PyList_GetItem(Py_kys, j);
    kys_[j] = PyFloat_AsDouble(ky);

    PyObject *Cv = PyList_GetItem(Py_Cvs, j);
    Cvs_[j] = PyFloat_AsDouble(Cv);
  }

  return convert_complex_array(integral(ds_,kxs_,kys_,Cvs_), n_omegas);
}


/*
----------------------------------------------------------------------------------------------------
THE WRAPPED INTEGRAND FUNCTION
----------------------------------------------------------------------------------------------------
*/


static PyObject *Integrand(PyObject *self, PyObject *args)
{
  double lambda;
  PyObject *Py_ds;
  PyObject *Py_kxs;
  PyObject *Py_kys;
  PyObject *Py_Cvs;

  if (!PyArg_ParseTuple(args,"dOOOO",&lambda,&Py_ds,&Py_kxs,&Py_kys,&Py_Cvs))
    return NULL;

  int n_layers = PyObject_Length(Py_ds);
  if (n_layers < 0 || n_layers > MAX_n_layers)
    return NULL;

  for (int j = 0; j < n_layers; j++) {
    PyObject *d = PyList_GetItem(Py_ds, j);
    ds_[j] = PyFloat_AsDouble(d);

    PyObject *kx = PyList_GetItem(Py_kxs, j);
    kxs_[j] = PyFloat_AsDouble(kx);

    PyObject *ky = PyList_GetItem(Py_kys, j);
    kys_[j] = PyFloat_AsDouble(ky);

    PyObject *Cv = PyList_GetItem(Py_Cvs, j);
    Cvs_[j] = PyFloat_AsDouble(Cv);
  }

  return convert_complex_array(integrand(lambda,ds_,kxs_,kys_,Cvs_), n_omegas);
}


/*
----------------------------------------------------------------------------------------------------
MODULE DEFINITIONS
----------------------------------------------------------------------------------------------------
*/
static PyMethodDef Integrate_FunctionsTable[] = {
  {"start", StartModule, METH_VARARGS, "allocate memory for the integrate module"},
  {"stop", StopModule, METH_VARARGS, "deallocate memory for the integrate module"},
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
  return PyModule_Create(&Integrate_Module);
}






/*
---------------------------------------------------------------------------------------




                          MAIN FOR TESTING




----------------------------------------------------------------------------------------
*/

/*
int main(void)
{
  // pretend these were parsed from Py input
  int n_omegas0 = 39;
  double omegas0[39] = {
    3141.592653589793, 3453.2524929747874, 3795.8303622241842,
    4172.393466187293, 4586.313289954719, 5041.295784799573,
    5541.414548697101, 6091.147298498122, 6695.415960302848,
    7359.630736976494, 8089.73854735669, 8892.2762708438,
    9774.429274089282, 10744.095743789525, 11809.957401576448,
    12981.5572341369, 14269.384934502563, 15684.970819492475,
    17240.989064179787, 18951.371177672692, 22898.00049022626,
    25169.583073304595, 27666.51666176646, 30411.157068693716,
    33428.07790237048, 36744.29058134312, 40389.4861759443,
    44396.30124156193, 48800.61002374043, 53641.84565131603,
    58963.3531933268, 64812.77774072464, 71242.49098750357,
    78310.06013055627, 86078.7632864636, 94618.15603990141,
    104004.69419614648, 114322.41831337663, 125663.70614359173};
  double b0 = 3.4285e-05;
  double lambda_i0 = 0.001;
  double lambda_f0 = 1e7;
  int n_pts0 = 200;
  int n_layers0 = 3;
  char b_type0 = 's';

  // this block would happen in the wrapped function
  start_module(b0,lambda_i0,lambda_f0,n_pts0,n_layers0,b_type0,n_omegas0);
  for (int i = 0; i < n_omegas; i++)
    omegas[i] = omegas0[i];  // "parsing" PyList_GetItem

  for (int i = 0; i < n_pts; i++)
    printf("%.10f\n", lambdas_[i]);
  printf("lambdas-------------------------------------\n");

  // fake iteration variables
  double ds0[3] = {50e-9, 3e-6, 300e-6};
  double kxs0[3] = {1.0, 10.0, 100.0};
  double kys0[3] = {1.0, 10.0, 100.0};
  double Cvs0[3] = {1706100.0, 2124260.0, 1630300.0};
  double complex *test_fB = fB(0,4500.0,kxs0,kys0,Cvs0);
  for (int i = 0; i < n_omegas; i++) {
    printf("%.10f ", creal(test_fB[i]));
    printf("%.10f\n", cimag(test_fB[i]));
  }
  printf("fB-------------------------------------\n");

 double complex *test_phi = phi(0,4500.0,ds0,kxs0,kys0,Cvs0);
 for (int i = 0; i < n_omegas; i++) {
   printf("%.10f ", creal(test_phi[i]));
   printf("%.10f\n", cimag(test_phi[i]));
 }
 printf("fA-------------------------------------\n");

double complex *test_fA = fA(0,4500.0,ds0,kxs0,kys0,Cvs0);
for (int i = 0; i < n_omegas; i++) {
  printf("%.10f ", creal(test_fA[i]));
  printf("%.10f\n", cimag(test_fA[i]));
}
printf("integrand-------------------------------------\n");

double complex *test_integrand = integrand(4500.0,ds0,kxs0,kys0,Cvs0);
for (int i = 0; i < n_omegas; i++) {
 printf("%.10f ", creal(test_integrand[i]));
 printf("%.10f\n", cimag(test_integrand[i]));
}
printf("integral-------------------------------------\n");

double complex *test_integral = integral(ds0,kxs0,kys0,Cvs0);
for (int i = 0; i < n_omegas; i++) {
 printf("%.10f ", creal(test_integral[i]));
 printf("%.10f\n", cimag(test_integral[i]));
}

  stop_module();

  return 0;
}
*/
