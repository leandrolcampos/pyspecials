// Copyright 2023 The PySpecials Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif /* PY_SSIZE_T_CLEAN */

#include <Python.h>
#include <math.h>

#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"

#ifndef PYSPECIALS_DOC_GENERATED
#define STRINGIZE(x) #x
#define UFUNC_DOC(NAME) STRINGIZE(NAME##_docstring)
#else
#include "toms708_doc_generated.h"
#define UFUNC_DOC(NAME) DOC_##NAME
#endif

/* ---------------------- Fortran function declaration ---------------------- */

extern void
bratio_(double* a,
        double* b,
        double* x,
        double* y,
        double* w,
        double* w1,
        int* ierr);

/* --------------------------- wrapper functions ---------------------------- */

double
betainc(double a, double b, double x)
{
  if (isnan(a) || isnan(b) || isnan(x)) {
    return NAN;
  }

  double y = 1.0 - x;
  double w = 0.0;
  double w1 = 0.0;
  int ierr = 0;

  bratio_(&a, &b, &x, &y, &w, &w1, &ierr);

  if (ierr != 0) {
    return NAN;
  }

  return w;
}

double
betaincc(double a, double b, double x)
{
  if (isnan(a) || isnan(b) || isnan(x)) {
    return NAN;
  }

  double y = 1.0 - x;
  double w = 0.0;
  double w1 = 0.0;
  int ierr = 0;

  bratio_(&a, &b, &x, &y, &w, &w1, &ierr);

  if (ierr != 0) {
    return NAN;
  }

  return w1;
}

/* ---------------------------- ufunc definition ---------------------------- */

static void
double_betainc(char** args,
               const npy_intp* dimensions,
               const npy_intp* steps,
               void* data)
{
  npy_intp i;
  npy_intp n = dimensions[0];
  char *in1 = args[0], *in2 = args[1], *in3 = args[2];
  char* out = args[3];
  npy_intp in1_step = steps[0], in2_step = steps[1], in3_step = steps[2];
  npy_intp out_step = steps[3];

  double a = 0.0;
  double b = 0.0;
  double x = 0.0;

  for (i = 0; i < n; i++) {
    /* BEGIN main ufunc computation */
    a = *(double*)in1;
    b = *(double*)in2;
    x = *(double*)in3;
    *((double*)out) = betainc(a, b, x);
    /* END main ufunc computation */

    in1 += in1_step;
    in2 += in2_step;
    in3 += in3_step;
    out += out_step;
  }
}

static void
double_betaincc(char** args,
                const npy_intp* dimensions,
                const npy_intp* steps,
                void* data)
{
  npy_intp i;
  npy_intp n = dimensions[0];
  char *in1 = args[0], *in2 = args[1], *in3 = args[2];
  char* out = args[3];
  npy_intp in1_step = steps[0], in2_step = steps[1], in3_step = steps[2];
  npy_intp out_step = steps[3];

  double a = 0.0;
  double b = 0.0;
  double x = 0.0;

  for (i = 0; i < n; i++) {
    /* BEGIN main ufunc computation */
    a = *(double*)in1;
    b = *(double*)in2;
    x = *(double*)in3;
    *((double*)out) = betaincc(a, b, x);
    /* END main ufunc computation */

    in1 += in1_step;
    in2 += in2_step;
    in3 += in3_step;
    out += out_step;
  }
}

/* --------------------------- ufunc registration --------------------------- */

#define FUNC_ARRAY_NAME(NAME) NAME##_funcs

#define UFUNC_FUNC_ARRAY_DOUBLE(NAME)                                          \
  static PyUFuncGenericFunction FUNC_ARRAY_NAME(NAME)[] = { &double_##NAME }

UFUNC_FUNC_ARRAY_DOUBLE(betainc);
UFUNC_FUNC_ARRAY_DOUBLE(betaincc);

static const char double_4_times[] = { NPY_DOUBLE,
                                       NPY_DOUBLE,
                                       NPY_DOUBLE,
                                       NPY_DOUBLE };

typedef struct ufunc_descriptor_struct
{
  PyUFuncGenericFunction* funcs;
  char* types;
  int ntypes;
  int nin;
  int nout;
  const char* name;
  const char* doc;
} UFUNC_DESCRIPTOR_t;

UFUNC_DESCRIPTOR_t ufunc_descriptors[] = { { FUNC_ARRAY_NAME(betainc),
                                             double_4_times,
                                             1,
                                             3,
                                             1,
                                             "betainc",
                                             UFUNC_DOC(betainc) },
                                           { FUNC_ARRAY_NAME(betaincc),
                                             double_4_times,
                                             1,
                                             3,
                                             1,
                                             "betaincc",
                                             UFUNC_DOC(betaincc) } };

static int
addUfuncs(PyObject* dictionary)
{
  PyObject* fun;
  int i;
  const int ufunc_count =
    sizeof(ufunc_descriptors) / sizeof(ufunc_descriptors[0]);
  for (i = 0; i < ufunc_count; i++) {
    UFUNC_DESCRIPTOR_t* d = &ufunc_descriptors[i];
    fun = PyUFunc_FromFuncAndData(d->funcs,
                                  NULL,
                                  d->types,
                                  d->ntypes,
                                  d->nin,
                                  d->nout,
                                  PyUFunc_None,
                                  d->name,
                                  d->doc,
                                  0);
    if (fun == NULL) {
      return -1;
    }

    int ret = PyDict_SetItemString(dictionary, d->name, fun);
    Py_DECREF(fun);
    if (ret < 0) {
      return -1;
    }
  }
  return 0;
}

/* ------------------------- module initialization -------------------------- */

#define TOMS708_MODULE_NAME "toms708"

static PyMethodDef Toms708Methods[] = { { NULL, NULL, 0, NULL } };

static struct PyModuleDef moduledef = { PyModuleDef_HEAD_INIT,
                                        TOMS708_MODULE_NAME,
                                        NULL,
                                        -1,
                                        Toms708Methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL };

PyMODINIT_FUNC
PyInit_toms708(void)
{
  PyObject* m;
  PyObject* d;
  PyObject* version;

  m = PyModule_Create(&moduledef);
  if (m == NULL) {
    return NULL;
  }

  import_array();
  import_ufunc();

  d = PyModule_GetDict(m);
  if (d == NULL) {
    return NULL;
  }

  /* Load the ufunc operators into the module's namespace */
  if (addUfuncs(d) < 0) {
    return NULL;
  }

  return m;
}