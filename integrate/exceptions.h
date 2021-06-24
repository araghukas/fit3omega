#include <Python.h>


static const char *ARGS_ERROR_MSG = "invalid (missing or extra) arguments";
static const char *LENGTH_ERROR_MSG = "array length incompatible with sample configuration";


PyObject *BT_IntegralError;
static const char *BT_IntegralError_NAME = "integrate.BT_IntegralError";

PyObject *BT_SetArgsError;
static const char *BT_SetArgsError_NAME = "integrate.BT_SetArgsError";

PyObject *BT_NotSetError;
static const char *BT_NotSetError_NAME = "integral.BT_NotSetError";
static const char *BT_NotSetError_MSG = "BT_Set method has not been called yet";

PyObject *OGC_IntegralError;
static const char *OGC_IntegralError_NAME = "integrate.OGC_IntegralError";

PyObject *OGC_IntegralDerError;
static const char *OGC_IntegralDerError_NAME = "integrate.OGC_IntegralDerError";

PyObject *OGC_SetArgsError;
static const char *OGC_SetArgsError_NAME = "integrate.OGC_SetArgsError";

PyObject *OGC_NotSetError;
static const char *OGC_NotSetError_NAME = "integral.OGC_NotSetError";
static const char *OGC_NotSetError_MSG = "OGC_Set method has not been called yet";

PyObject *n_LAYERS_Error;
static const char *n_LAYERS_Error_NAME = "integrate.NumberOfLayersError";
static const char *n_LAYERS_Error_MSG = "invalid number of sample layers (exceeds max?)";

PyObject *n_OMEGAS_Error;
static const char *n_OMEGAS_Error_NAME = "integrate.NumberOfOmegasError";
static const char *n_OMEGAS_Error_MSG = "invalid length of omegas array (exceeds max?)";

PyObject *NullOmegasError;
static const char *NullOmegasError_NAME = "integrate.NullOmegasError";
static const char *NullOmegasError_MSG = "pointer to omegas array is null";

PyObject *ParameterIDError;
static const char *ParameterIDError_NAME = "integrate.ParameterIDError";
static const char *ParameterIDError_MSG = "encountered invalid parameter ID";


void init_exceptions(void);
