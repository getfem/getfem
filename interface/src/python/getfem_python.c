// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2004-2008 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License,
// or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================
#include <Python.h>
//#include "numarray/libnumarray.h"
//#include "numarray/arrayobject.h"
////#include "Numeric/arrayobject.h"
#include "numpy/arrayobject.h"
#include "structmember.h"
#include "gfi_array.h"
#include "getfem_interface.h"
static PyObject *call_getfem(PyObject *self, PyObject *args);
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
GetfemObject_name(PyGetfemObject* self)
{
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
  else if (self->objid < other->objid) return +1;
  else return 0;
}

static PyMethodDef module_methods[] = {
    {"getfem",  call_getfem, METH_VARARGS,
     "Execute a getfem command."},
    {"getfem_from_constructor",  call_getfem_from_constructor, METH_VARARGS,
     "internal -- Execute a getfem command for building a new object."},
    //{"register_types", register_types, METH_VARARGS, "register the derived types (internal function)"},
    {"register_python_factory", register_python_factory, METH_VARARGS, "register (on initialization) the python function which is used to build objects from a GetfemObject type (internal function)"},
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
     "Object ID in the Getfem++ workspace"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyGetfemObject_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size (deprecated) */
    "_getfem.GetfemObject",     /*tp_name*/
    sizeof(PyGetfemObject),    /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    (cmpfunc)GetfemObject_compare, /*tp_compare -- necessary for dictionary */
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    (hashfunc)GetfemObject_hash,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "Generic Getfem++ objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    GetfemObject_methods,             /* tp_methods */
    GetfemObject_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,/*(initproc)GetfemObject_init,*/      /* tp_init */
    0,                       /* tp_alloc */
    0/*GetfemObject_new*/,                 /* tp_new */
    0, /* tp_free */
    0, /* tp_is_gc */
    0, /* tp_bases */
    0, /* tp_mro */
    0, /* tp_cache */
    0, /* tp_subclasses */
    0, /* tp_weaklist */
    0/*(destructor)GetfemObject_del*/, /* tp_del */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_getfem(void)
{
  PyObject *m;
  PyGetfemObject_Type.tp_new = PyType_GenericNew;
  if (PyType_Ready(&PyGetfemObject_Type) < 0)
    return;
  m = Py_InitModule3("_getfem", module_methods, "The Getfem-python interface module.");
  import_array(); /* init Numeric */
  //import_libnumarray(); /* init Numarray */
  Py_INCREF(&PyGetfemObject_Type);
  PyModule_AddObject(m, "GetfemObject", (PyObject *)&PyGetfemObject_Type);
}





#define COLLECTCHUNK 2
typedef struct ptr_collect {
  void *p[COLLECTCHUNK];
  int n;
  struct ptr_collect *next;
} ptr_collect;

typedef struct gcollect {
  ptr_collect *allocated;
  ptr_collect *pyobjects;
} gcollect;

static ptr_collect *
ptr_collect_push_front(ptr_collect *col, void *p) {
  if (col == NULL || col->n == COLLECTCHUNK) {
    ptr_collect *pcol = col;
    col = malloc(sizeof(ptr_collect));
    col->next = pcol;
    col->n = 1; col->p[0] = p;
  } else {
    col->p[col->n++] = p;
  }
  return col;
}

/* mark a pyobject as referenced */
static void
gc_ref(gcollect *gc, PyObject *o) {
  gc->pyobjects = ptr_collect_push_front(gc->pyobjects, o);
}

/* allocate a collectable chunk of memory */
static void *
gc_alloc(gcollect *gc, size_t sz) {
  //printf("gc_alloc(%lu)\n", sz);
  void *p = malloc(sz == 0 ? 1 : sz);
  if (p) {
    gc->allocated = ptr_collect_push_front(gc->allocated, p);
  } else {
    PyErr_Format(PyExc_RuntimeError, "could not allocate %d bytes: memory exhausted", (int)sz);
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
      /*for (i=0; i < p->n; ++i)
	Py_DECREF((PyObject*)p->p[i]);*/
      np = p->next; free(p);
    }
  gc->pyobjects = NULL;
  //fprintf(stderr, "gc->allocated=%p\n", gc->allocated);
  for (p = gc->allocated; p; p = np) {
    //fprintf(stderr, "release bloc: n=%d, next=%p\n", p->n, p->next);
    for (i=0; i < p->n; ++i) {
      //fprintf(stderr, " i=%d release %p\n", i, p->p[i]);
      free(p->p[i]);
    }
    np = p->next; free(p);
  }
  gc->allocated = NULL;
}

/* returns 0 for unsupported types */
static inline unsigned
itemsize_of_PyArray(PyArrayObject *pa) {
  switch (pa->descr->type_num) {
    case PyArray_DOUBLE : return sizeof(double);
    case PyArray_INT:     return sizeof(int);
    case PyArray_CDOUBLE: return sizeof(double)*2;
    default: PyErr_SetString(PyExc_RuntimeError, "invalid numeric array type"); return 0;
  }
}

/* Check if the strides are compatible with fortran ..
   always return true for 1D arrays, or 1x1xN arrays etc.
*/
static int PyArray_IsFortranCompatible(PyArrayObject *pa) {
  unsigned itemsize=itemsize_of_PyArray(pa);
  if (itemsize == 0) return 0;
  int i, s = itemsize;
  for (i=0; i < pa->nd; ++i) {
    if (pa->dimensions[i] > 1) {
      if (pa->strides[i] != s) return 0;
    }
    s *= pa->dimensions[i];
  }
  return 1;
}

/* transfert the content of pa into a fortran_array (all gfi_arrays follow
 fortran ordering convention ... why did the numpy people choose c-order ?)
 if revert is non null, the fortran data is copied to the PyArray
 returns -1 on error, 0 on success
*/
static int
copy_PyArray_data(PyArrayObject *pa, void *fortran_array, int revert) {
  unsigned itemsize, size;
  if ((itemsize = itemsize_of_PyArray(pa)) == 0) return -1;
  size = PyArray_Size((PyObject*)pa);
  /*printf("copying array of size %d (x%d) dir=%d\n", size, itemsize, revert);
  { int i; printf("dim/strides = "); for (i = 0; i < pa->nd; ++i) printf("%d/%d ", pa->dimensions[i], pa->strides[i]); printf("\n"); }
  */
  if (pa->nd == 0) return 0;
  else if (PyArray_IsFortranCompatible(pa)) {
    /* fine, we can just copy raw memory */
    if (revert == 0) memcpy(fortran_array,pa->data,size*itemsize);
    else memcpy(pa->data,fortran_array,size*itemsize);
  } else {
    /* arg, we need to convert between the pyarray and fortran ordering */
    int i,j;
    int cnt[pa->nd], pos = 0;
    //printf("not fortran compatible..\n");
    for (i = 0; i < pa->nd; ++i) cnt[i] = 0;
    for (i = 0; i < size; ++i, fortran_array += itemsize) {
      //printf("i=%d, cnt= ",i); for (j = 0; j < pa->nd; ++j) printf("%d ",cnt[j]); printf("pos = %d",pos);
      switch (pa->descr->type_num) {
        case PyArray_CDOUBLE: {
          if (revert == 0) {
            ((double*)fortran_array)[0] = ((double*)(pa->data+pos))[0];
            ((double*)fortran_array)[1] = ((double*)(pa->data+pos))[1];
          } else {
            ((double*)(pa->data+pos))[0] = ((double*)fortran_array)[0];
            ((double*)(pa->data+pos))[1] = ((double*)fortran_array)[1];
          }
        } break;
        case PyArray_DOUBLE: {
          if (revert == 0) *((double*)fortran_array) = *((double*)(pa->data+pos));
          else *((double*)(pa->data+pos)) = *((double*)fortran_array);
          /*if (revert == 0) printf(" -> %g\n", *((double*)fortran_array));
            else printf(" -> %g\n", *((double*)(pa->data+pos)));*/
        } break;
        case PyArray_INT: {
          if (revert == 0) *((int*)fortran_array) = *((int*)(pa->data+pos));
          else *((int*)(pa->data+pos)) = *((int*)fortran_array);
          //if (revert == 0) printf(" -> %d\n", *((int*)fortran_array));
          //else printf(" -> %d\n", *((int*)(pa->data+pos)));
        } break;
      }
      j = 0;
      while (1) {
	cnt[j]++; pos += pa->strides[j];
	if (cnt[j] >= pa->dimensions[j]) {
	  pos -= pa->strides[j]*pa->dimensions[j];
	  cnt[j] = 0; ++j;
	  if (j == pa->nd) break;
	} else break;
      }
    }
  }
  return 0;
}

int
PyObject_is_GetfemObject(PyObject *o, gfi_object_id *pid) {
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
  if (PyString_Check(o)) {
    /* for strings, the pointer is shared, no copy */
    char *s = PyString_AsString(o);
    gc_ref(gc, o);
    t->storage.type = GFI_CHAR;
    t->dim.dim_len = 1; t->dim.dim_val = &TGFISTORE(char,len);
    TGFISTORE(char,len)=strlen(s);
    TGFISTORE(char,val)=s;
  } else if (PyInt_Check(o)) {
    /* usual python integer */
    t->storage.type = GFI_INT32;
    t->dim.dim_len = 0; t->dim.dim_val = &TGFISTORE(int32,len);
    TGFISTORE(int32,len)=1;
    if (!(TGFISTORE(int32,val)=gc_alloc(gc,sizeof(int)))) return NULL;
    TGFISTORE(int32,val)[0] = (int)PyInt_AsLong(o); /* should check overflow */
  } else if (PyFloat_Check(o)) {
    /* usual python float */
    t->storage.type = GFI_DOUBLE;
    t->dim.dim_len = 0; t->dim.dim_val = &TGFISTORE(double,len);
    TGFISTORE(double,len)=1;
    t->storage.gfi_storage_u.data_double.is_complex = 0;
    if (!(TGFISTORE(double,val)=gc_alloc(gc,sizeof(double)))) return NULL;
    TGFISTORE(double,val)[0] = PyFloat_AsDouble(o);
  } else if (PyComplex_Check(o)) {
    /* usual python complex number */
    t->storage.type = GFI_DOUBLE;
    t->dim.dim_len = 0; t->dim.dim_val = &TGFISTORE(double,len);
    TGFISTORE(double,len)=2;
    t->storage.gfi_storage_u.data_double.is_complex = 1;
    if (!(TGFISTORE(double,val)=gc_alloc(gc,sizeof(double)*2))) return NULL;
    TGFISTORE(double,val)[0] = PyComplex_RealAsDouble(o);
    TGFISTORE(double,val)[1] = PyComplex_ImagAsDouble(o);
  } else if (PyTuple_Check(o)) {
    /* python tuples are stored in 'cell arrays' (i.e. matlab's lists of inhomogeneous elements) */
    int i;
    t->storage.type = GFI_CELL;
    t->dim.dim_len = 1; t->dim.dim_val = &TGFISTORE(cell,len);
    TGFISTORE(cell,len) = PyTuple_GET_SIZE(o);
    if (!(TGFISTORE(cell,val) = gc_alloc(gc,sizeof(gfi_array*)*TGFISTORE(cell,len)))) return NULL;
    gfi_array **p = TGFISTORE(cell,val);
    for (i=0; i < PyTuple_GET_SIZE(o); ++i) {
      p[i] = PyObject_to_gfi_array(gc, PyTuple_GET_ITEM(o,i));
      if (!p[i]) return NULL;
    }
  } else if (PyObject_is_GetfemObject(o, &id)) {
    /* getfem objects are refered to with a couple (classid, objectid) */
    t->storage.type = GFI_OBJID;
    t->dim.dim_len = 1; t->dim.dim_val = &TGFISTORE(cell,len);
    t->storage.gfi_storage_u.objid.objid_len = 1;
    if (!(t->storage.gfi_storage_u.objid.objid_val =
	  gc_alloc(gc,sizeof(gfi_object_id)))) return NULL;
    t->storage.gfi_storage_u.objid.objid_val[0] = id;
  } else {
    /* and finally, numpy arrays (and anything convertible to a numpy array) */
    int type = -1;
    if (!PyArray_Check(o)) {
      switch (PyArray_ObjectType(o, PyArray_INT)) {
        case PyArray_INT:
#ifdef PyArray_UNSIGNED_TYPES
        case PyArray_UINT:
#endif
          //        case PyArray_LONG:
          type = PyArray_INT; break;
        case PyArray_FLOAT:
        case PyArray_DOUBLE:
          type = PyArray_DOUBLE; break;
        case PyArray_CDOUBLE:
        case PyArray_CFLOAT:
          type = PyArray_CDOUBLE; break;
      }
    } else if (PyArray_CanCastSafely(((PyArrayObject*)o)->descr->type_num, PyArray_INT)) {
      //printf("got a pyarray(INT)\n");
      type = PyArray_INT;
    } else if (PyArray_CanCastSafely(((PyArrayObject*)o)->descr->type_num, PyArray_DOUBLE)) {
      type = PyArray_DOUBLE;
    } else if (PyArray_CanCastSafely(((PyArrayObject*)o)->descr->type_num, PyArray_CDOUBLE)) {
      type = PyArray_CDOUBLE;
    }
    PyArrayObject *ao = NULL;
    if (type != -1) ao = (PyArrayObject*) PyArray_FromObject(o, type,0,0);

    if (!ao) {
      PyObject *stype = PyObject_Str((PyObject*)o->ob_type);
      PyErr_Format(PyExc_RuntimeError, "unhandled argument type: %s", PyString_AsString(stype));
      Py_DECREF(stype);
      return NULL;
    }
    gc_ref(gc, (PyObject*)ao);

    switch (type) {
      case PyArray_CDOUBLE:
      case PyArray_DOUBLE: {
        t->storage.type = GFI_DOUBLE;
        t->storage.gfi_storage_u.data_double.is_complex = (type == PyArray_CDOUBLE) ? 1 : 0;
        unsigned ndouble = PyArray_Size((PyObject*)ao) * ((type == PyArray_CDOUBLE) ? 2 : 1);
        if (PyArray_IsFortranCompatible(ao)) {
          TGFISTORE(double,val) = (double*)ao->data; // no copy
        } else {
          if (!(TGFISTORE(double,val) = gc_alloc(gc, ndouble * sizeof(double))))
            return NULL;
          if (copy_PyArray_data(ao, TGFISTORE(double,val), 0) != 0) return NULL;
        }
        TGFISTORE(double,len) = ndouble;
      } break;
      case PyArray_INT: {
        t->storage.type = GFI_INT32;
        unsigned nint = PyArray_Size((PyObject*)ao);
        if (PyArray_IsFortranCompatible(ao)) {
          TGFISTORE(int32,val) = (int*)ao->data; // no copy
        } else {
          if (!(TGFISTORE(int32,val) = gc_alloc(gc, nint * sizeof(int))))
            return NULL;
          if (copy_PyArray_data(ao, TGFISTORE(int32,val), 0) != 0) return NULL;
        }
      } break;
    }

    t->dim.dim_len = ao->nd;
    int *d = ao->dimensions;
    t->dim.dim_val = (u_int*)d;
  }
  return t;
}

static PyObject*
PyGetfemObject_FromObjId(gfi_object_id id, int in__init__) {
  PyObject *o;
  PyGetfemObject *go = PyObject_New(PyGetfemObject, &PyGetfemObject_Type); Py_INCREF(go);
  //printf("PyGetfemObject_FromObjId(cid=%d, oid=%d,in__init__=%d)\n", id.cid,id.id,in__init__);
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
build_gfi_array_list(gcollect *gc, PyObject *tuple, char **pfunction_name, int *nb) {
  const gfi_array **l;
  int i, j;
  if (PyTuple_GET_SIZE(tuple) == 0) {
    PyErr_SetString(PyExc_RuntimeError, "missing function name"); return NULL;
  }
  if (!PyString_Check(PyTuple_GET_ITEM(tuple,0))) {
    PyErr_SetString(PyExc_RuntimeError, "expecting function name as a string"); return NULL;
  }
  *pfunction_name = PyString_AsString(PyTuple_GET_ITEM(tuple,0));
  *nb = PyTuple_GET_SIZE(tuple) - 1;
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
  assert(t);
  switch (t->storage.type) {
  case GFI_UINT32:
  case GFI_INT32: {
    if (t->dim.dim_len == 0) return PyInt_FromLong(TGFISTORE(int32,val)[0]);
    else {
      if (!(o = PyArray_FromDims(t->dim.dim_len, (int*)t->dim.dim_val, PyArray_INT))) return NULL;
      if (copy_PyArray_data((PyArrayObject*)o, TGFISTORE(int32,val), 1) != 0) return NULL;
    }
  } break;
  case GFI_DOUBLE: {
    if (!gfi_array_is_complex(t)) {
      if (t->dim.dim_len == 0) return PyFloat_FromDouble(TGFISTORE(double,val)[0]);
      else {
        //int i;  printf("received an array of dimension %d, val=", t->dim.dim_len); for (i=0; i < t->dim.dim_len; ++i) printf("%d ",t->dim.dim_val[i]); printf("\n");
        if (!(o = PyArray_FromDims(t->dim.dim_len, (int*)t->dim.dim_val, PyArray_DOUBLE))) return NULL;
      }
    } else {
      if (t->dim.dim_len == 0) return PyComplex_FromDoubles(TGFISTORE(double,val)[0], TGFISTORE(double,val)[1]);
      else {
        if (!(o = PyArray_FromDims(t->dim.dim_len, (int*)t->dim.dim_val, PyArray_CDOUBLE))) return NULL;
      }
    }
    if (copy_PyArray_data((PyArrayObject*)o, TGFISTORE(double,val), 1) != 0) return NULL;
  } break;
  case GFI_CHAR: {
    o = PyString_FromStringAndSize(TGFISTORE(char,val),TGFISTORE(char,len));
  } break;
  case GFI_CELL: {
    unsigned i;
    if (!(o = PyTuple_New(TGFISTORE(cell,len)))) return NULL;
    for (i=0; i < TGFISTORE(cell,len); ++i) {
      PyObject *to = gfi_array_to_PyObject(TGFISTORE(cell,val)[i], in__init__);
      if (!to) return NULL;
      PyTuple_SET_ITEM(o,i,to);
    }
  } break;
  case GFI_OBJID: {
    if (t->storage.gfi_storage_u.objid.objid_len != 1) {
#if 0
      /* PyArray_OBJECT is not supported in numarray ... */
      int i;
      if (!(o = PyArray_FromDims(t->dim.dim_len, (int*)t->dim.dim_val, PyArray_OBJECT))) return NULL;
      if (!PyArray_IsFortranCompatible((PyArrayObject*)o)) { // I'm just too lazy to transpose matrices
	PyErr_Format(PyExc_RuntimeError, "cannot return %d-D array of %d getfem objects",
		     t->dim.dim_len, t->storage.gfi_storage_u.objid.objid_len);
	return NULL;
      }
      for (i = 0; i < t->storage.gfi_storage_u.objid.objid_len; ++i) {
	((PyObject**)(((PyArrayObject*)o)->data))[i] =
	  PyGetfemObject_FromObjId(t->storage.gfi_storage_u.objid.objid_val[i], in__init__);
      }
#else
      /* return a python list to be on the safe side */
      int nb = t->storage.gfi_storage_u.objid.objid_len, i;
      if (t->dim.dim_len != 1) {
	PyErr_Format(PyExc_RuntimeError, "cannot return %d-D array of %d getfem objects",
		     t->dim.dim_len, nb);
      }
      if (!(o = PyList_New(nb))) return NULL;
      for (i=0; i < nb; ++i) {
	PyList_SetItem(o, i, PyGetfemObject_FromObjId(t->storage.gfi_storage_u.objid.objid_val[i], in__init__));
      }
#endif
    } else {
      o = PyGetfemObject_FromObjId(t->storage.gfi_storage_u.objid.objid_val[0], in__init__);
    }
  } break;
  case GFI_SPARSE: {
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
    errmsg = getfem_interface_main(PYTHON_INTERFACE, function_name, in_cnt, in, &out_cnt, &out, &infomsg);
    if (infomsg) {
      printf("message from gf_%s follow:\n%s\n", function_name, infomsg); fflush(stdout);
    }
    if (errmsg) {
      if (strstr(errmsg, "Internal error:"))
        PyErr_Format(PyExc_AssertionError, "(Getfem::InternalError) -- %s", errmsg);
      else
        PyErr_Format(PyExc_RuntimeError, "(Getfem::InterfaceError) -- %s", errmsg);
    } else {
      //fprintf(stderr, "%s : success, nb_out = %d\n", function_name, out_cnt);
      if (out_cnt == 0) {
	result = Py_None; Py_INCREF(Py_None);
      } else if (out) {
	int i, err = 0;
	PyObject *d[out_cnt];
	for (i = 0; i < out_cnt; ++i) {
	  if (!err && !(d[i] = gfi_array_to_PyObject(out[i], in__init__))) err = 1;
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
call_getfem(PyObject *self, PyObject *args) { return call_getfem_(self,args, 0); }
static PyObject*
call_getfem_from_constructor(PyObject *self, PyObject *args) { return call_getfem_(self,args, 1); }

/*static PyObject *
register_types(PyObject *self, PyObject *args)
{
  printf("registering types..\n");
  if (PyArg_ParseTuple(args,"OOOOOOOO",
		       &PyDerivedTypes[MESH_CLASS_ID],
		       &PyDerivedTypes[MESHFEM_CLASS_ID],
		       &PyDerivedTypes[GEOTRANS_CLASS_ID],
		       &PyDerivedTypes[FEM_CLASS_ID],
		       &PyDerivedTypes[INTEG_CLASS_ID],
		       &PyDerivedTypes[ELTM_CLASS_ID],
		       &PyDerivedTypes[CVSTRUCT_CLASS_ID],
		       &PyDerivedTypes[POLY_CLASS_ID],
		       &PyDerivedTypes[SLICE_CLASS_ID])) return NULL;
  //Py_INCREF(PyDerivedTypes[MESH_CLASS_ID]);
  PyObject_Print(PyDerivedTypes[MESH_CLASS_ID],stderr,0);
  if (!PyClass_Check(PyDerivedTypes[MESH_CLASS_ID])) {
    PyErr_Format(PyExc_RuntimeError, "Not a class..");
    return NULL;
  }
  return Py_None;
  }*/

/* copied verbatim from the "Extending and Embedding the Python Interpreter" tutorial */
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
