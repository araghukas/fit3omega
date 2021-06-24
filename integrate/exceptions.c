#include "exceptions.h"


void init_exceptions(void)
{
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