/*===========================================================================

 Copyright (C) 2004-2020 Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <Python.h>
#include "numpy/arrayobject.h"
#include "structmember.h"
#include <string.h>

#include "gfi_array.h"
#include "getfem_interface.h"
#include "getfem_arch_config.h"
#include <assert.h>

#if PY_MAJOR_VERSION >= 3
#define PyString_AsString(o) PyUnicode_AsUTF8(o)
#define PyString_FromFormat(a,b,c) PyUnicode_FromFormat(a,b,c)
#define PyString_Check(o) PyUnicode_Check(o)
#define PyInt_Check(o) PyLong_Check(o)
#define PyInt_AsLong(o) PyLong_AsLong(o)
#define PyString_FromString(o) PyUnicode_FromString(o)
#define PyString_FromStringAndSize(o,l) PyUnicode_FromStringAndSize(o,l)
#define PyInt_FromLong(o) PyLong_FromLong(o)
#endif

static PyObject *call_getfem(PyObject *self, PyObject *args);
static PyObject *getfem_env(PyObject *self, PyObject *args);
static PyObject *call_getfem_from_constructor(PyObject *self, PyObject *args);
//static PyObject *register_types(PyObject *self, PyObject *args);
static PyObject *register_python_factory(PyObject *dummy, PyObject *args);

//static PyObject *PyDerivedTypes[GETFEMINT_NB_CLASS] = {NULL};
static PyObject *python_factory = NULL; // PyCallable

/* definition of the root GetfemObject */

typedef struct PyGetfemObject {
  PyObject_HEAD
  unsigned classid, objid;
} PyGetfemObject;

static PyObject *
GetfemObject_name(PyGetfemObject *self) {
  return PyString_FromFormat("getfem.GetfemObject(classid=%d,objid=%d)",
                             self->classid, self->objid);
}

static int
GetfemObject_hash(PyGetfemObject *key) {
  return key->objid + (key->classid << 14);
}

static int
GetfemObject_compare(PyGetfemObject *self, PyGetfemObject *other) {
  if (self->classid < other->classid) return -1;
  if (self->classid > other->classid) return +1;
  if (self->objid < other->objid) return -1;
  if (self->objid > other->objid) return +1;
  return 0;
}

static PyObject *
GfObject_richcompare(PyGetfemObject *self, PyGetfemObject *other, int op) {
  int bc = GetfemObject_compare(self, other);
  switch(op) {
  case Py_LT : if (bc <  0) Py_RETURN_TRUE; else Py_RETURN_FALSE;
  case Py_LE : if (bc <= 0) Py_RETURN_TRUE; else Py_RETURN_FALSE;
  case Py_EQ : if (bc == 0) Py_RETURN_TRUE; else Py_RETURN_FALSE;
  case Py_NE : if (bc != 0) Py_RETURN_TRUE; else Py_RETURN_FALSE;
  case Py_GT : if (bc == 1) Py_RETURN_TRUE; else Py_RETURN_FALSE;
  case Py_GE : if (bc >= 0) Py_RETURN_TRUE; else Py_RETURN_FALSE;
  }
  return NULL;
}

static PyMethodDef module_methods[] = {
    {"getfem", call_getfem, METH_VARARGS, "Execute a getfem command."},
    {"getfem_env", getfem_env, METH_VARARGS,
     "Builder variables for documentation"},
    {"getfem_from_constructor",  call_getfem_from_constructor, METH_VARARGS,
     "internal -- Execute a getfem command for building a new object."},
    //{"register_types", register_types, METH_VARARGS,
    // "register the derived types (internal function)"},
    {"register_python_factory", register_python_factory, METH_VARARGS,
     "register (on initialization) the python function which is used to "
     "build objects from a GetfemObject type (internal function)"},
    {NULL}        /* Sentinel */
};

static PyMethodDef GetfemObject_methods[] = {
    {"name", (PyCFunction)GetfemObject_name, METH_NOARGS,
     "Return a string combining classID and objID."
    },
    {NULL}
};

static PyMemberDef GetfemObject_members[] = {
    {"classid", T_INT, offsetof(PyGetfemObject, classid), 0,
     "Class ID"},
    {"objid", T_INT, offsetof(PyGetfemObject, objid), 0,
     "Object ID in the GetFEM workspace"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyGetfemObject_Type = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                                 /* ob_size (deprecated) */
#endif
    "_getfem.GetfemObject",            /* tp_name */
    sizeof(PyGetfemObject),            /* tp_basicsize */
    0,                                 /* tp_itemsize */
    0,                                 /* tp_dealloc */
    0,                                 /* tp_print */
    0,                                 /* tp_getattr */
    0,                                 /* tp_setattr */
#if PY_MAJOR_VERSION >= 3
    GetfemObject_compare,              /* tp_compare, necessary for dictionary*/
#else
    (cmpfunc)GetfemObject_compare,     /* tp_compare, necessary for dictionary*/
#endif
    0,                                 /* tp_repr */
    0,                                 /* tp_as_number */
    0,                                 /* tp_as_sequence */
    0,                                 /* tp_as_mapping */
    (hashfunc)GetfemObject_hash,       /* tp_hash */
    0,                                 /* tp_call */
    0,                                 /* tp_str */
    0,                                 /* tp_getattro */
    0,                                 /* tp_setattro */
    0,                                 /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,/* tp_flags */
    "Generic GetFEM objects",          /* tp_doc */
    0,                                 /* tp_traverse */
    0,                                 /* tp_clear */
    (richcmpfunc)GfObject_richcompare, /* tp_richcompare */
    0,                                 /* tp_weaklistoffset */
    0,                                 /* tp_iter */
    0,                                 /* tp_iternext */
    GetfemObject_methods,              /* tp_methods */
    GetfemObject_members,              /* tp_members */
    0,                                 /* tp_getset */
    0,                                 /* tp_base */
    0,                                 /* tp_dict */
    0,                                 /* tp_descr_get */
    0,                                 /* tp_descr_set */
    0,                                 /* tp_dictoffset */
    0,/*(initproc)GetfemObject_init*/  /* tp_init */
    0,                                 /* tp_alloc */
    0,/*GetfemObject_new*/             /* tp_new */
    0,                                 /* tp_free */
    0,                                 /* tp_is_gc */
    0,                                 /* tp_bases */
    0,                                 /* tp_mro */
    0,                                 /* tp_cache */
    0,                                 /* tp_subclasses */
    0,                                 /* tp_weaklist */
    0,/*(destructor)GetfemObject_del*/ /* tp_del */
};

#ifndef PyMODINIT_FUNC        /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_getfem",     /* m_name */
  "getfem-python3 interface module.",  /* m_doc */
  -1,                  /* m_size */
  module_methods,    /* m_methods */
  NULL,                /* m_reload */
  NULL,                /* m_traverse */
  NULL,                /* m_clear */
  NULL,                /* m_free */
};

PyMODINIT_FUNC
PyInit__getfem(void)
{
  PyObject *m;
  PyGetfemObject_Type.tp_new = PyType_GenericNew;
  if (PyType_Ready(&PyGetfemObject_Type) < 0)
    return NULL;
  m = PyModule_Create(&moduledef);
  import_array(); /* init Numpy */
  Py_INCREF(&PyGetfemObject_Type);
  PyModule_AddObject(m, "GetfemObject", (PyObject *)&PyGetfemObject_Type);
  return m;
}

#else

PyMODINIT_FUNC
init_getfem(void)
{
  PyObject *m;
  PyGetfemObject_Type.tp_new = PyType_GenericNew;
  if (PyType_Ready(&PyGetfemObject_Type) < 0)
    return;
  m = Py_InitModule3("_getfem", module_methods,
                     "python-getfem interface module.");
  import_array(); /* init Numpy */
  Py_INCREF(&PyGetfemObject_Type);
  PyModule_AddObject(m, "GetfemObject", (PyObject *)&PyGetfemObject_Type);
}

#endif

#define COLLECTCHUNK 2
typedef struct ptr_collect {
  void *p[COLLECTCHUNK];
  int n;
  struct ptr_collect *next;
  int owned[COLLECTCHUNK];
} ptr_collect;

typedef struct gcollect {
  ptr_collect *allocated;
  ptr_collect *pyobjects;
} gcollect;

static ptr_collect *
ptr_collect_push_front(ptr_collect *col, void *p, int owned) {
  if (col == NULL || col->n == COLLECTCHUNK) {
    ptr_collect *pcol = col;
    col = malloc(sizeof(ptr_collect));
    col->next = pcol;
    col->n = 1;
    col->p[0] = p;
    col->owned[0] = owned;
  } else {
    col->p[col->n] = p;
    col->owned[col->n++] = owned;
  }
  return col;
}

/* mark a pyobject as referenced */
static void
gc_ref(gcollect *gc, PyObject *o, int owned) {
  gc->pyobjects = ptr_collect_push_front(gc->pyobjects, o, owned);
}

/* allocate a collectable chunk of memory */
static void *
gc_alloc(gcollect *gc, size_t sz) {
  //printf("gc_alloc(%lu)\n", sz);
  void *p = malloc(sz == 0 ? 1 : sz);
  if (p) {
    gc->allocated = ptr_collect_push_front(gc->allocated, p, 1);
  } else {
    PyErr_Format(PyExc_RuntimeError,
                 "could not allocate %d bytes: memory exhausted", (int)sz);
  }
  return p;
}

/* release allocated memory and decrement refcount of
   referenced objects */
static void
gc_release(gcollect *gc) {
  ptr_collect *p, *np;
  int i;
  if (!PyErr_Occurred())
    for (p = gc->pyobjects; p; p = np) {
      for (i=0; i < p->n; ++i)
        if (p->owned[i])
          Py_DECREF((PyObject*)p->p[i]);
      np = p->next; free(p);
    }
  gc->pyobjects = NULL;
  //fprintf(stderr, "gc->allocated=%p\n", gc->allocated);
  for (p = gc->allocated; p; p = np) {
    //fprintf(stderr, "release bloc: n=%d, next=%p\n", p->n, p->next);
    for (i=0; i < p->n; ++i) {
      if (p->owned[i])
      //fprintf(stderr, " i=%d release %p\n", i, p->p[i]);
        free(p->p[i]);
    }
    np = p->next; free(p);
  }
  gc->allocated = NULL;
}



int
PyObject_is_GetfemObject(PyObject *o, gfi_object_id *pid)
{
  PyGetfemObject *go = NULL;
  PyObject *attr_id = NULL;
  if (PyObject_TypeCheck(o, &PyGetfemObject_Type)) {
    go = (PyGetfemObject*)o;
  } else if ((attr_id = PyObject_GetAttrString(o,"id"))) {
    if (PyObject_TypeCheck(attr_id, &PyGetfemObject_Type))
      go = (PyGetfemObject*)attr_id;
  }
  PyErr_Clear(); /* clear error flag if attribute "id" was not found */
  if (go && pid) {
    pid->cid = go->classid; pid->id = go->objid;
  }
  if (attr_id) { Py_DECREF(attr_id); }
  return go != NULL;
}

static gfi_array *
PyObject_to_gfi_array(gcollect *gc, PyObject *o)
{
  gfi_object_id id;
  gfi_array *t = gc_alloc(gc, sizeof(gfi_array)); if (!t) return NULL;
#define TGFISTORE0(T) t->storage.gfi_storage_u.data_##T
#define TGFISTORE1(T,field) data_##T##_##field
#define TGFISTORE(T,field) TGFISTORE0(T).TGFISTORE1(T,field)
  PyErr_Clear();
  if (PyString_Check(o)) {
    //printf("String\n");
    /* for strings, the pointer is shared, no copy */
    int L = (int)(strlen(PyString_AsString(o)));
    char *s = PyString_AsString(o);
    gc_ref(gc, o, 0);

    t->storage.type = GFI_CHAR;
    t->dim.dim_len = 1; t->dim.dim_val = &TGFISTORE(char,len);
    TGFISTORE(char,len)=L;
    TGFISTORE(char,val)=s;
  } else if (PyInt_Check(o) || PyLong_Check(o)) {
    //printf("Int or Long\n");
    /* usual python integer */
    int d = (int)PyInt_AsLong(o);
    if (PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_OverflowError))
      return (gfi_array *)PyErr_Format(PyExc_OverflowError,
                                       "in getfem interface.");

    t->storage.type = GFI_INT32;
    t->dim.dim_len = 0; t->dim.dim_val = &TGFISTORE(int32,len);
    TGFISTORE(int32,len)=1;
    if (!(TGFISTORE(int32,val)=gc_alloc(gc,sizeof(int)))) return NULL;
    TGFISTORE(int32,val)[0] = d;
  } else if (PyFloat_Check(o)) {
    //printf("Float\n");
    /* usual python float */
    double df = PyFloat_AsDouble(o);
    t->storage.type = GFI_DOUBLE;
    t->dim.dim_len = 0; t->dim.dim_val = &TGFISTORE(double,len);
    TGFISTORE(double,len)=1;
    t->storage.gfi_storage_u.data_double.is_complex = 0;
    if (!(TGFISTORE(double,val)=gc_alloc(gc,sizeof(double)))) return NULL;
    TGFISTORE(double,val)[0] = df;
  } else if (PyComplex_Check(o)) {
    //printf("Complex\n");
    /* usual python complex number */
    double real = PyComplex_RealAsDouble(o);
    double imag = PyComplex_ImagAsDouble(o);
    t->storage.type = GFI_DOUBLE;
    t->dim.dim_len = 0; t->dim.dim_val = &TGFISTORE(double,len);
    TGFISTORE(double,len)=2;
    t->storage.gfi_storage_u.data_double.is_complex = 1;
    if (!(TGFISTORE(double,val)=gc_alloc(gc,sizeof(double)*2))) return NULL;
    TGFISTORE(double,val)[0] = real;
    TGFISTORE(double,val)[1] = imag;
  } else if (PyTypeNum_ISNUMBER(PyArray_ObjectType(o,0))) {
    //printf("Numerical Array");
    /* python numeric sequences are stored in numerical array */
    int dtype = PyArray_ObjectType(o,0);

    PyObject *po = NULL;
    switch (dtype) {
      case NPY_BOOL:
      case NPY_BYTE:
      case NPY_UBYTE:
      case NPY_SHORT:
      case NPY_USHORT:
      case NPY_INT:
      case NPY_UINT:
      case NPY_LONG:
      case NPY_ULONG:
      case NPY_LONGLONG:
      case NPY_ULONGLONG:
        t->storage.type = GFI_INT32;

        if (PyArray_NDIM((PyArrayObject *)o) == 1) /* Is there a bug in
                                                      PyArray_CheckFromAny ? */
          po = PyArray_CheckFromAny(o,PyArray_DescrFromType(NPY_INT),0,0,
                                    NPY_ARRAY_FORCECAST | NPY_ARRAY_OUT_ARRAY
                                    | NPY_ARRAY_ELEMENTSTRIDES, NULL);
        else
          po = PyArray_CheckFromAny(o,PyArray_DescrFromType(NPY_INT),0,0,
                                    NPY_ARRAY_FORCECAST | NPY_ARRAY_OUT_FARRAY
                                    | NPY_ARRAY_ELEMENTSTRIDES, NULL);
        if (!po) { PyErr_NoMemory(); return NULL;}

        gc_ref(gc, po, 1);
        /* No new copy. */
        TGFISTORE(int32,val) = (int *)PyArray_DATA((PyArrayObject *)po);
        break;
      case NPY_FLOAT:
      case NPY_DOUBLE:
      case NPY_LONGDOUBLE:
        t->storage.type = GFI_DOUBLE;
        t->storage.gfi_storage_u.data_double.is_complex = 0;

        if (PyArray_NDIM((PyArrayObject *)o) == 1) /* Is there a bug in
                                                      PyArray_CheckFromAny ? */
          po = PyArray_CheckFromAny(o,PyArray_DescrFromType(NPY_DOUBLE),0,0,
                                    NPY_ARRAY_FORCECAST | NPY_ARRAY_OUT_ARRAY
                                    | NPY_ARRAY_ELEMENTSTRIDES, NULL);
        else
          po = PyArray_CheckFromAny(o,PyArray_DescrFromType(NPY_DOUBLE),0,0,
                                    NPY_ARRAY_FORCECAST | NPY_ARRAY_OUT_FARRAY
                                    | NPY_ARRAY_ELEMENTSTRIDES, NULL);
        if (!po) { PyErr_NoMemory(); return NULL;}

        gc_ref(gc, po, 1);
        /* No new copy. */
        TGFISTORE(double,val) = (double *)PyArray_DATA((PyArrayObject *)po);
        break;
      case NPY_CFLOAT:
      case NPY_CDOUBLE:
      case NPY_CLONGDOUBLE:
        t->storage.type = GFI_DOUBLE;
        t->storage.gfi_storage_u.data_double.is_complex = 1;

        if (PyArray_NDIM((PyArrayObject *)o) == 1) /* is there a bug in
                                                      PyArray_CheckFromAny ? */
          po = PyArray_CheckFromAny(o,PyArray_DescrFromType(NPY_CDOUBLE),0,0,
                                    NPY_ARRAY_FORCECAST | NPY_ARRAY_OUT_ARRAY
                                    | NPY_ARRAY_ELEMENTSTRIDES, NULL);
        else
          po = PyArray_CheckFromAny(o,PyArray_DescrFromType(NPY_CDOUBLE),0,0,
                                    NPY_ARRAY_FORCECAST | NPY_ARRAY_OUT_FARRAY
                                    | NPY_ARRAY_ELEMENTSTRIDES, NULL);
        if (!po) { PyErr_NoMemory(); return NULL;}

        gc_ref(gc, po, 1);
        /* No new copy. */
        TGFISTORE(double,val) = (double *)PyArray_DATA((PyArrayObject *)po);
        break;
      default: {
        PyObject *sdtype =PyObject_Str((PyObject*)PyArray_DescrFromType(dtype));
        PyErr_Format(PyExc_RuntimeError, "invalid numeric dtype: %s",
          PyString_AsString(sdtype));
        Py_DECREF(sdtype);
        return NULL;
      }
    }
    t->dim.dim_len = PyArray_NDIM((PyArrayObject *)po);
    t->dim.dim_val = (u_int *)gc_alloc(gc, t->dim.dim_len * sizeof(u_int));

    int i;
    for (i=0; i < t->dim.dim_len; ++i)
      t->dim.dim_val[i] = (u_int)PyArray_DIM((PyArrayObject *)po,i);
  } else if (PyTuple_Check(o) || PyList_Check(o)) {
    //printf("Tuple or List\n");
    /* python tuples and lists are stored in 'cell arrays'
       (i.e. matlab's lists of inhomogeneous elements) */
    int i;
    t->storage.type = GFI_CELL;
    t->dim.dim_len = 1; t->dim.dim_val = &TGFISTORE(cell,len);

    if (PyTuple_Check(o)) TGFISTORE(cell,len) = (unsigned)(PyTuple_GET_SIZE(o));
    else TGFISTORE(cell,len) = (unsigned)(PyList_GET_SIZE(o));

    if (!(TGFISTORE(cell,val)
          = gc_alloc(gc,sizeof(gfi_array*)*TGFISTORE(cell,len)))) return NULL;
    gfi_array **p = TGFISTORE(cell,val);

    for (i=0; i < TGFISTORE(cell,len); ++i) {
      if (PyTuple_Check(o))
        p[i] = PyObject_to_gfi_array(gc, PyTuple_GET_ITEM(o,i));
      else p[i] = PyObject_to_gfi_array(gc, PyList_GET_ITEM(o,i));
      if (!p[i]) return NULL;
    }
  } else if (PyObject_is_GetfemObject(o, &id)) {
    //printf("Object\n");
    /* getfem objects are refered to with a couple (classid, objectid) */
    t->storage.type = GFI_OBJID;
    t->dim.dim_len = 1; t->dim.dim_val = &TGFISTORE(cell,len);
    t->storage.gfi_storage_u.objid.objid_len = 1;
    if (!(t->storage.gfi_storage_u.objid.objid_val =
          gc_alloc(gc,sizeof(gfi_object_id)))) return NULL;
    t->storage.gfi_storage_u.objid.objid_val[0] = id;
  } else {
    int dtype = PyArray_ObjectType(o,0);
    PyObject *stype = PyObject_Str((PyObject*)o->ob_type);
    PyObject *sdtype = PyObject_Str((PyObject*)PyArray_DescrFromType(dtype));
    PyErr_Format(PyExc_RuntimeError,
                 "unhandled argument (type, dtype): (%s, %s)",
                 PyString_AsString(stype), PyString_AsString(sdtype));
    Py_DECREF(stype);
    Py_DECREF(sdtype);
    return NULL;
  }
  return t;
}

static PyObject*
PyGetfemObject_FromObjId(gfi_object_id id, int in__init__) {
  PyObject *o;
  PyGetfemObject *go = PyObject_New(PyGetfemObject, &PyGetfemObject_Type);
  Py_INCREF(go);
  //printf("PyGetfemObject_FromObjId(cid=%d, oid=%d,in__init__=%d)\n",
  //       id.cid,id.id,in__init__);
  if (!go) return NULL;
  go->classid = id.cid; go->objid = id.id;
  if (!in__init__) {
    PyObject *arg;
    if (!(arg = Py_BuildValue("(O)", go))) return NULL;
    //printf("  -> arg= "); PyObject_Print(arg,stdout,0); printf("\n");
    o = PyEval_CallObject(python_factory, arg);
    Py_DECREF(arg);
  } else o = (PyObject*)go;
  //printf("  -> return "); PyObject_Print(o,stdout,0); printf("\n");
  return o;
}

static const gfi_array **
build_gfi_array_list(gcollect *gc, PyObject *tuple, char **pfunction_name,
                     int *nb) {
  const gfi_array **l;
  int i, j;
  if (PyTuple_GET_SIZE(tuple) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "missing function name"); return NULL;
  }
  if (!PyString_Check(PyTuple_GET_ITEM(tuple,0))) {
    PyErr_SetString(PyExc_RuntimeError, "expecting function name as a string");
    return NULL;
  }
  *pfunction_name = PyString_AsString(PyTuple_GET_ITEM(tuple,0));
  *nb = (int)(PyTuple_GET_SIZE(tuple) - 1);
  if (!(l = gc_alloc(gc, sizeof(gfi_array*) * *nb))) return NULL;
  for (i=0, j = 0; i < *nb; ++i) {
    PyObject *o = PyTuple_GET_ITEM(tuple,i+1);
    if (o != Py_None) {
      gfi_array *g = PyObject_to_gfi_array(gc, o);
      if (g) { l[j++] = g; } else return NULL;
    }
  }
  *nb = j;
  return l;
}

PyObject*
gfi_array_to_PyObject(gfi_array *t, int in__init__) {
  PyObject *o = NULL;
  // assert(t);
  switch (t->storage.type) {
  case GFI_UINT32:
  case GFI_INT32: {
    //printf("GFI_INT32\n");
    if (t->dim.dim_len == 0) return PyInt_FromLong(TGFISTORE(int32,val)[0]);
    else {
      npy_intp *dim = PyDimMem_NEW(t->dim.dim_len);
      int i;
      for(i=0; i < t->dim.dim_len; i++)
        dim[i] = (npy_intp)t->dim.dim_val[i];
      if (!(o = PyArray_EMPTY(t->dim.dim_len, dim, NPY_INT, 1))) return NULL;
      PyDimMem_FREE(dim);

      npy_intp itemsize = PyArray_ITEMSIZE((PyArrayObject *)o);
      npy_intp size = PyArray_Size(o); /* Number of elements. */
      memcpy(PyArray_DATA((PyArrayObject *)o), TGFISTORE(int32,val),
             size*itemsize); // new copy
    }
  } break;
  case GFI_DOUBLE: {
    // printf("GFI_DOUBLE\n");
    if (!gfi_array_is_complex(t)) {
      if (t->dim.dim_len == 0)
        return PyFloat_FromDouble(TGFISTORE(double,val)[0]);
      else {
        npy_intp *dim = PyDimMem_NEW(t->dim.dim_len);
        int i;
        for(i=0; i< t->dim.dim_len; i++)
          dim[i] = (npy_intp)t->dim.dim_val[i];
        if (!(o = PyArray_EMPTY(t->dim.dim_len, dim, NPY_DOUBLE, 1)))
          return NULL;
        PyDimMem_FREE(dim);
      }
    } else {
      if (t->dim.dim_len == 0)
        return PyComplex_FromDoubles(TGFISTORE(double,val)[0],
                                     TGFISTORE(double,val)[1]);
      else {
        npy_intp *dim = PyDimMem_NEW(t->dim.dim_len);
        int i;
        for(i=0; i< t->dim.dim_len; i++)
          dim[i] = (npy_intp)t->dim.dim_val[i];
        if (!(o = PyArray_EMPTY(t->dim.dim_len, dim, NPY_CDOUBLE, 1)))
          return NULL;
        PyDimMem_FREE(dim);
      }
    }
    npy_intp itemsize = PyArray_ITEMSIZE((PyArrayObject *)o);
    npy_intp size = PyArray_Size(o); /* Number of elements. */
    memcpy(PyArray_DATA((PyArrayObject *)o), TGFISTORE(double,val),
           size*itemsize); // new copy
  } break;
  case GFI_CHAR: {
    //printf("GFI_CHAR\n");
    o = PyString_FromStringAndSize(TGFISTORE(char,val),TGFISTORE(char,len));
  } break;
  case GFI_CELL: {
    //printf("GFI_CELL\n");
    unsigned i;
    if (!(o = PyTuple_New(TGFISTORE(cell,len)))) return NULL;
    for (i=0; i < TGFISTORE(cell,len); ++i) {
      PyObject *to = gfi_array_to_PyObject(TGFISTORE(cell,val)[i], in__init__);
      if (!to) return NULL;
      PyTuple_SET_ITEM(o,i,to);
    }
  } break;
  case GFI_OBJID: {
    //printf("GFI_OBJID\n");
    int nb = t->storage.gfi_storage_u.objid.objid_len;
    if (nb != 1) {
#if 0
      /* PyArray_OBJECT is not supported in numarray ... */
      npy_intp *dim = PyDimMem_NEW(t->dim.dim_len);
      int i;
      for(i=0; i< t->dim.dim_len; i++)
        dim[i] = (npy_intp)t->dim.dim_val[i];
      if (!(o = PyArray_EMPTY(t->dim.dim_len, dim, NPY_OBJECT,1))) return NULL;

      if (!PyArray_ISFARRAY(PyArray_DATA((PyArrayObject *)o))) {
        // I'm just too lazy to transpose matrices
        PyErr_Format(PyExc_RuntimeError,
                     "cannot return %d-D array of %d getfem objects",
                     t->dim.dim_len, nb);
        return NULL;
      }
      for (i = 0; i<nb; ++i) {
        (PyObject*)PyArray_GETPTR1((PyArrayObject*)o,i) = // not compiling
          PyGetfemObject_FromObjId(t->storage.gfi_storage_u.objid.objid_val[i],
                                   in__init__);
      }
#else
      /* return a python list to be on the safe side */
      if (t->dim.dim_len != 1) {
        PyErr_Format(PyExc_RuntimeError,
                     "cannot return %d-D array of %d getfem objects",
                     t->dim.dim_len, nb);
      }
      if (!(o = PyList_New(nb))) return NULL;

      int i;
      for (i=0; i<nb; ++i) {
        PyList_SetItem(o, i,
                       PyGetfemObject_FromObjId
                       (t->storage.gfi_storage_u.objid.objid_val[i],
                        in__init__));
      }
#endif
    } else {
      o = PyGetfemObject_FromObjId(t->storage.gfi_storage_u.objid.objid_val[0],
                                   in__init__);
    }
  } break;
  case GFI_SPARSE: {
    //printf("GFI_SPARSE\n");
    PyErr_SetString(PyExc_RuntimeError,
                    "Numpy does not have Native sparse matrices. "
                    "Use getfem sparse objects instead.");
  } break;
  default:  {
    assert(0);
  } break;
  }
  return o;
}

static PyObject *
call_getfem_(PyObject *self, PyObject *args, int in__init__)
{
  const gfi_array **in = 0;
  gfi_array **out = 0;
  int in_cnt = 0, out_cnt = -1;
  char *function_name, *infomsg, *errmsg;
  gcollect gc;
  PyObject *result = NULL;
  gc.allocated = gc.pyobjects = NULL;
  assert(PyTuple_Check(args));
  //printf("calling getfem: args=%d\n", PyTuple_GET_SIZE(args));fflush(stdout);
  in = build_gfi_array_list(&gc, args, &function_name, &in_cnt);
  if (in) {
    //fprintf(stdout,"  -> function = %s\n", function_name);
    Py_BEGIN_ALLOW_THREADS;
    errmsg = getfem_interface_main(PYTHON_INTERFACE, function_name, in_cnt,
                                   in, &out_cnt, &out, &infomsg,0);
    Py_END_ALLOW_THREADS;
    if (infomsg) {
      printf("message from gf_%s follow:\n%s\n", function_name, infomsg);
      fflush(stdout);
    }
    if (errmsg) {
      if (strstr(errmsg, "Internal error:"))
        PyErr_Format(PyExc_AssertionError, "(Getfem::InternalError) -- %s",
                     errmsg);
      else
        PyErr_Format(PyExc_RuntimeError, "(Getfem::InterfaceError) -- %s",
                     errmsg);
    } else {
      //fprintf(stderr, "%s : success, nb_out = %d\n", function_name, out_cnt);
      if (out_cnt == 0) {
        result = Py_None; Py_INCREF(Py_None);
      } else if (out) {
        int i, err = 0;
        PyObject *d[out_cnt];
        for (i = 0; i < out_cnt; ++i) {
          if (!err && !(d[i] = gfi_array_to_PyObject(out[i], in__init__)))
            err = 1;
          gfi_array_destroy(out[i]);
        }

        free(out);
        if (!err) {
          if (out_cnt > 1) {
            result = PyTuple_New(out_cnt);
            for (i = 0; i < out_cnt; ++i) PyTuple_SET_ITEM(result,i,d[i]);
          } else result = d[0];
        }
      }
    }
  }
  gc_release(&gc);
  return PyErr_Occurred() ? NULL : result;
}

static PyObject*
call_getfem(PyObject *self, PyObject *args)
{ return call_getfem_(self,args, 0); }
static PyObject*
call_getfem_from_constructor(PyObject *self, PyObject *args)
{ return call_getfem_(self,args, 1); }

static PyObject *
getfem_env(PyObject *self, PyObject *args) {
  char* word_in;

  size_t size = PyTuple_GET_SIZE(args);
  if (size != 1) {
    PyErr_Format(PyExc_TypeError,
                 "getfem_env() takes exactly 1 argument (%d given)",
                 (int)size);
    return NULL;
  } else if (!PyArg_ParseTuple(args,"s",&word_in)) {
    return NULL;
  }

  PyObject* word_out;

  if (strcmp(word_in,"project") == 0) {
    word_out = PyString_FromString("GetFEM");
  } else if (strcmp(word_in,"copyright") == 0) {
    word_out = PyString_FromString
    ("2004-2020 GetFEM project");
  } else if (strcmp(word_in,"authors") == 0) {
    word_out = PyString_FromString
    ("Yves Renard, Julien Pommier, Konstantinos Poulios");
  } else if (strcmp(word_in,"url") == 0) {
    word_out = PyString_FromString("http://home.gna.org/getfem/");
  } else if (strcmp(word_in,"license") == 0) {
    word_out = PyString_FromString("GNU LGPL v3");
  } else if (strcmp(word_in,"package") == 0) {
    word_out = PyString_FromString(GETFEM_PACKAGE);
  } else if (strcmp(word_in,"package_name") == 0) {
    word_out = PyString_FromString(GETFEM_PACKAGE_NAME);
  } else if (strcmp(word_in,"package_string") == 0) {
    word_out = PyString_FromString(GETFEM_PACKAGE_STRING);
  } else if(strcmp(word_in,"package_tarname") == 0) {
    word_out = PyString_FromString(GETFEM_PACKAGE_TARNAME);
  } else if(strcmp(word_in,"package_version") == 0 ||
            strcmp(word_in,"release") == 0) {
    word_out = PyString_FromString(GETFEM_PACKAGE_VERSION);
  } else if(strcmp(word_in,"version") == 0) {
    word_out = PyString_FromString(GETFEM_VERSION);
  } else {
    word_out = PyString_FromString("");
  }

  Py_INCREF(word_out);
  return word_out;
}

/* Copied verbatim from the "Extending and Embedding the Python Interpreter"
   tutorial */
static PyObject *
register_python_factory(PyObject *dummy, PyObject *args)
{
  PyObject *result = NULL;
  PyObject *temp;

  if (PyArg_ParseTuple(args, "O:register_python_factory", &temp)) {
    if (!PyCallable_Check(temp)) {
      PyErr_SetString(PyExc_TypeError, "parameter must be callable");
      return NULL;
    }
    Py_XINCREF(temp);         /* Add a reference to new callback */
    Py_XDECREF(python_factory);  /* Dispose of previous callback */
    python_factory = temp;       /* Remember new callback */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    result = Py_None;
  }
  return result;
}
