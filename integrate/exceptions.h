// compile-time constants for exception names
#ifndef _PYTHON_H
#include <Python.h>
#endif


const char *ARGS_ERROR_MSG = "invalid (missing or extra) arguments";
const char *LENGTH_ERROR_MSG = "array length incompatible with sample configuration";


PyObject *BT_IntegralError;
const char *BT_IntegralError_NAME = "integrate.BT_IntegralError";

PyObject *BT_SetError;
const char *BT_SetError_NAME = "integrate.BT_SetError";

PyObject *OGC_IntegralError;
const char *OGC_IntegralError_NAME = "integrate.OGC_IntegralError";

PyObject *OGC_SetError;
const char *OGC_SetError_NAME = "integrate.OGC_SetError";

PyObject *N_LAYERS_Error;
const char *N_LAYERS_Error_NAME = "integrate.NumberOfLayersError";
const char *N_LAYERS_Error_MSG = "invalid number of sample layers (exceeds max?)";

PyObject *n_OMEGAS_Error;
const char *n_OMEGAS_Error_NAME = "integrate.NumberOfOmegasError";
const char *n_OMEGAS_Error_MSG = "invalid length of omegas array (exceeds max?)";

PyObject *NullOmegasError;
const char *NullOmegasError_NAME = "integrate.NullOmegasError";
const char *NullOmegasError_MSG = "pointer to omegas array is null";


void init_exceptions(void) {
    BT_IntegralError = PyErr_NewException(BT_IntegralError_NAME, NULL, NULL);
    BT_SetError = PyErr_NewException(BT_SetError_NAME, NULL, NULL);

    OGC_IntegralError = PyErr_NewException(OGC_IntegralError_NAME, NULL, NULL);
    OGC_SetError = PyErr_NewException(OGC_SetError_NAME, NULL, NULL);

    N_LAYERS_Error = PyErr_NewException(N_LAYERS_Error_NAME, NULL, NULL);

    NullOmegasError = PyErr_NewException(NullOmegasError_NAME, NULL, NULL);
    n_OMEGAS_Error = PyErr_NewException(n_OMEGAS_Error_NAME, NULL, NULL);
}