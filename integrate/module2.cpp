#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>


typedef struct {
  PyObject_HEAD
  double ob_fval; // <--- would be a standard python float
  // contents of object
} Integrator;

static PyTypeObject IntegratorType = {
  // flags and function pointers that the interpreter inspects
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "intglib.Integrator",
  .tp_doc = "Fast integrator of Borca-Tascuic Eq (1)",
  .tp_basicsize = sizeof(Integrator),
  .tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_new = PyType_GenericNew, // <--- default __new__() equivalent
};

static PyModuleDef intglib = {
  PyModuleDef_HEAD_INIT,
  .m_name = "intglib",
  .m_doc = "C++ module for Borca-Tascuic Eq (1) integration.",
  .m_size = -1,
};

PyMODINIT_FUNC PyInit_intglib(void) {
  if (PyType_Ready(&IntegratorType) < 0)
    return NULL;

  PyObject* m = PyModule_Create(&intglib);
  if (m == NULL)
    return NULL;

  Py_INCREF(&Integrator);
  if (PyModule_AddObject(m, "Integrator", (PyObject*) &IntegratorType) < 0) {
    Py_DECREF(&IntegratorType);
    Py_DECREF(m);
    return NULL;
  }

  return m;
}
