// compile-time constants for exception names
#ifndef _PYTHON_H
#include <Python.h>
#endif


const char *ARGS_ERROR_MSG = "invalid (missing or extra) arguments";
const char *LENGTH_ERROR_MSG = "array length incompatible with sample configuration";


PyObject *BT_IntegralError;
const char *BT_IntegralError_NAME = "integrate.BT_IntegralError";

PyObject *BT_SetArgsError;
const char *BT_SetArgsError_NAME = "integrate.BT_SetArgsError";

PyObject *BT_NotSetError;
const char *BT_NotSetError_NAME = "integral.BT_NotSetError";
const char *BT_NotSetError_MSG = "BT_Set method has not been called yet";

PyObject *OGC_IntegralError;
const char *OGC_IntegralError_NAME = "integrate.OGC_IntegralError";

PyObject *OGC_IntegralDerError;
const char *OGC_IntegralDerError_NAME = "integrate.OGC_IntegralDerError";

PyObject *OGC_SetArgsError;
const char *OGC_SetArgsError_NAME = "integrate.OGC_SetArgsError";

PyObject *OGC_NotSetError;
const char *OGC_NotSetError_NAME = "integral.OGC_NotSetError";
const char *OGC_NotSetError_MSG = "OGC_Set method has not been called yet";

PyObject *n_LAYERS_Error;
const char *n_LAYERS_Error_NAME = "integrate.NumberOfLayersError";
const char *n_LAYERS_Error_MSG = "invalid number of sample layers (exceeds max?)";

PyObject *n_OMEGAS_Error;
const char *n_OMEGAS_Error_NAME = "integrate.NumberOfOmegasError";
const char *n_OMEGAS_Error_MSG = "invalid length of omegas array (exceeds max?)";

PyObject *NullOmegasError;
const char *NullOmegasError_NAME = "integrate.NullOmegasError";
const char *NullOmegasError_MSG = "pointer to omegas array is null";

PyObject *ParameterIDError;
const char *ParameterIDError_NAME = "integrate.ParameterIDError";
const char *ParameterIDError_MSG = "encountered invalid parameter ID";


void init_exceptions(void) {
    BT_IntegralError = PyErr_NewException(BT_IntegralError_NAME, NULL, NULL);
    BT_NotSetError = PyErr_NewException(BT_NotSetError_NAME, NULL, NULL);

    OGC_IntegralError = PyErr_NewException(OGC_IntegralError_NAME, NULL, NULL);
    OGC_NotSetError = PyErr_NewException(OGC_NotSetError_NAME, NULL, NULL);
    OGC_IntegralDerError = PyErr_NewException(OGC_IntegralError_NAME, NULL, NULL);

    n_LAYERS_Error = PyErr_NewException(n_LAYERS_Error_NAME, NULL, NULL);

    NullOmegasError = PyErr_NewException(NullOmegasError_NAME, NULL, NULL);
    n_OMEGAS_Error = PyErr_NewException(n_OMEGAS_Error_NAME, NULL, NULL);

    ParameterIDError = PyErr_NewException(ParameterIDError_NAME, NULL, NULL);
}