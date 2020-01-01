/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gfi_array.h"


/* -------------------- creation ------------------------*/

void* gfi_malloc(size_t sz) {
  if (sz == 0) return malloc(1);
  else return malloc(sz);
}

void* gfi_calloc(size_t n, size_t m) {
  void *p = (n*m == 0) ? malloc(1) : calloc(n,m);
  /*printf("gfi_calloc(%d,%d) -> %p\n", n, m, p);*/
  return p;
}

void gfi_free(void *p) {
  /*printf("gfi_free  (%p)\n", p);*/
  if (p) free(p); 
}

#define FREE(p) { /*printf("%s@%d ", __FILE__, __LINE__); fflush(stdout); */gfi_free(p); p = NULL; }

gfi_array*
gfi_array_create(int ndim, int *dims, 
                 gfi_type_id type, gfi_complex_flag is_complex) {
  int i,sz;
  gfi_array *t = gfi_calloc(1, sizeof(gfi_array)); if (t == NULL) return NULL;
  t->dim.dim_len = ndim;
  t->dim.dim_val = gfi_calloc(ndim, sizeof(int)); if ( t->dim.dim_val == NULL) { gfi_free(t); return NULL; }

  /*printf("gfi_array_create(ndim = %d, type = %d, @%p %p\n", ndim, type, t, t->dim.dim_val);*/

  for (i=0,sz=1; i < ndim; ++i) { t->dim.dim_val[i] = dims[i]; sz *= dims[i]; }
  t->storage.type = type;
  switch (t->storage.type) {
  case GFI_CHAR: {
    t->storage.gfi_storage_u.data_char.data_char_len = sz;
    t->storage.gfi_storage_u.data_char.data_char_val = gfi_malloc(sz*sizeof(char));
    if (t->storage.gfi_storage_u.data_char.data_char_val == NULL) goto not_enough_mem;
  } break;
  case GFI_INT32: {
    t->storage.gfi_storage_u.data_int32.data_int32_len = sz;
    t->storage.gfi_storage_u.data_int32.data_int32_val = gfi_malloc(sz*sizeof(int));
    if (t->storage.gfi_storage_u.data_int32.data_int32_val == NULL) goto not_enough_mem;
  } break;
  case GFI_UINT32: {
    t->storage.gfi_storage_u.data_uint32.data_uint32_len = sz;
    t->storage.gfi_storage_u.data_uint32.data_uint32_val = gfi_malloc(sz*sizeof(int));
    if (t->storage.gfi_storage_u.data_uint32.data_uint32_val == NULL) goto not_enough_mem;
  } break;
  case GFI_CELL: {
    t->storage.gfi_storage_u.data_cell.data_cell_len = sz;
    t->storage.gfi_storage_u.data_cell.data_cell_val = (gfi_array**)gfi_calloc(sz,sizeof(gfi_array*));
    if (t->storage.gfi_storage_u.data_cell.data_cell_val == NULL) goto not_enough_mem;
  } break;
  case GFI_DOUBLE: {
    t->storage.gfi_storage_u.data_double.is_complex = is_complex;
    t->storage.gfi_storage_u.data_double.data_double_len = sz * (is_complex ? 2 : 1);
    t->storage.gfi_storage_u.data_double.data_double_val = gfi_calloc(sz,sizeof(double) * (is_complex ? 2 : 1));
    if (t->storage.gfi_storage_u.data_double.data_double_val == NULL) goto not_enough_mem;
  } break;
  case GFI_OBJID: {
    t->storage.gfi_storage_u.objid.objid_len = sz;
    t->storage.gfi_storage_u.objid.objid_val = gfi_calloc(sz, sizeof(gfi_object_id));  
    if (t->storage.gfi_storage_u.objid.objid_val == NULL) goto not_enough_mem;
  } break;
  default: {
    printf("internal error"); return NULL; 
  } break;
  }  
  return t;
 not_enough_mem:
  gfi_array_destroy(t); gfi_free(t); return NULL;
}

/* 0-D arrays are used for non-array scalars */
gfi_array*
gfi_array_create_0(gfi_type_id type, gfi_complex_flag is_complex) {
  return gfi_array_create(0, NULL, type, is_complex);
}

gfi_array*
gfi_array_create_1(int M, gfi_type_id type, gfi_complex_flag is_complex) {
  return gfi_array_create(1, &M, type, is_complex);
}

gfi_array*
gfi_array_create_2(int M, int N, gfi_type_id type, gfi_complex_flag is_complex) {
  int d[2]; d[0]=M; d[1]=N;
  return gfi_array_create(2, d, type, is_complex);
}

gfi_array*
gfi_array_from_string(const char *s) {
  gfi_array*t;
  int n = (int)strlen(s);
  t = gfi_array_create_1(n,GFI_CHAR, GFI_REAL);
  if (t) memcpy(gfi_char_get_data(t), s, n);
  return t;
}

gfi_array*
gfi_create_sparse(int m, int n, int nzmax, gfi_complex_flag is_complex) {
  gfi_array *t = gfi_calloc(1, sizeof(gfi_array));
  t->dim.dim_len = 2; 
  t->dim.dim_val = gfi_calloc(2, sizeof(int)); 
  t->dim.dim_val[0] = m; t->dim.dim_val[1] = n;
  t->storage.type = GFI_SPARSE;
  t->storage.gfi_storage_u.sp.is_complex = is_complex;
  t->storage.gfi_storage_u.sp.ir.ir_len = nzmax; t->storage.gfi_storage_u.sp.ir.ir_val = gfi_calloc(nzmax,sizeof(int));
  t->storage.gfi_storage_u.sp.jc.jc_len = n+1;   t->storage.gfi_storage_u.sp.jc.jc_val = gfi_calloc(n+1,sizeof(int));
  t->storage.gfi_storage_u.sp.pr.pr_len = nzmax * (is_complex ? 2 : 1); 
  t->storage.gfi_storage_u.sp.pr.pr_val = gfi_calloc(nzmax,sizeof(double) * (is_complex ? 2 : 1));
  if ((nzmax && (t->storage.gfi_storage_u.sp.ir.ir_val == NULL ||
		 t->storage.gfi_storage_u.sp.pr.pr_val == NULL)) ||
      t->storage.gfi_storage_u.sp.jc.jc_val == NULL) { gfi_array_destroy(t); t = NULL; }
  return t;
}
/*
gfi_array*
gfi_create_objid(int nid, unsigned *ids, unsigned cid) {
  int i;
  gfi_array *t = gfi_calloc(1, sizeof(gfi_array));
  t->dim.dim_len = nid; 
  t->dim.dim_val = gfi_calloc(1, sizeof(int)); t->dim.dim_val[0]=1;
  t->storage.type = GFI_OBJID;
  t->storage.gfi_storage_u.objid.objid_len = nid;
  t->storage.gfi_storage_u.objid.objid_val = gfi_calloc(nid, sizeof(gfi_object_id));
  for (i=0; i < nid; ++i) {
    t->storage.gfi_storage_u.objid.objid_val[i].id = ids[i];
    t->storage.gfi_storage_u.objid.objid_val[i].cid = cid;
  }
  return t;
}
*/
/* ----------------- destruction ------------ */

void 
gfi_array_destroy(gfi_array *t) {
  if (t == NULL) return;
  FREE(t->dim.dim_val);
  switch (t->storage.type) {
  case GFI_CHAR: {
    FREE(t->storage.gfi_storage_u.data_char.data_char_val);
  } break;
  case GFI_INT32: {
    FREE(t->storage.gfi_storage_u.data_int32.data_int32_val);
  } break;
  case GFI_UINT32: {
    FREE(t->storage.gfi_storage_u.data_uint32.data_uint32_val);
  } break;
  case GFI_CELL: {
    int i;
    if (t->storage.gfi_storage_u.data_cell.data_cell_len)
      assert(t->storage.gfi_storage_u.data_cell.data_cell_val);
    for (i=0; i < t->storage.gfi_storage_u.data_cell.data_cell_len; ++i)
      gfi_array_destroy(t->storage.gfi_storage_u.data_cell.data_cell_val[i]);
    FREE(t->storage.gfi_storage_u.data_cell.data_cell_val);
  } break;
  case GFI_DOUBLE: {
    FREE(t->storage.gfi_storage_u.data_double.data_double_val);
  } break;
  case GFI_SPARSE: {
    FREE(t->storage.gfi_storage_u.sp.ir.ir_val);
    FREE(t->storage.gfi_storage_u.sp.jc.jc_val);
    FREE(t->storage.gfi_storage_u.sp.pr.pr_val);
  }
  case GFI_OBJID: {
    FREE(t->storage.gfi_storage_u.objid.objid_val);
  } break;
  default: {
    assert(0);
  }
  }
}

/* ----------------- inquiry ----------------- */

int
gfi_array_get_ndim(const gfi_array*t) {
  assert(t);
  return t->dim.dim_len;
}

const int*
gfi_array_get_dim(const gfi_array*t) {
  assert(t);
  return (const int*)t->dim.dim_val;
}

unsigned 
gfi_array_nb_of_elements(const gfi_array *t) {
  unsigned i,sz=1;
  assert(t);
  if (t->storage.type != GFI_SPARSE) {
    for (i=0,sz=1; i < t->dim.dim_len; ++i) { sz *= t->dim.dim_val[i]; }
    return sz;
  } else return t->storage.gfi_storage_u.sp.pr.pr_len;
}

unsigned int*
gfi_sparse_get_ir(const gfi_array *t) {
  assert(t); 
  assert(t->storage.type == GFI_SPARSE);
  return (unsigned int*)t->storage.gfi_storage_u.sp.ir.ir_val;
}

unsigned int*
gfi_sparse_get_jc(const gfi_array *t) {
  assert(t); 
  assert(t->storage.type == GFI_SPARSE);
  return (unsigned int*)t->storage.gfi_storage_u.sp.jc.jc_val;
}

double*
gfi_sparse_get_pr(const gfi_array *t) {
  assert(t); 
  assert(t->storage.type == GFI_SPARSE);
  return t->storage.gfi_storage_u.sp.pr.pr_val;
}

char*
gfi_char_get_data(const gfi_array* t) {
  assert(t);
  assert(t->storage.type == GFI_CHAR);
  return t->storage.gfi_storage_u.data_char.data_char_val;
}

int*
gfi_int32_get_data(const gfi_array *t) {
  assert(t);
  assert(t->storage.type == GFI_INT32);
  return t->storage.gfi_storage_u.data_int32.data_int32_val;
}

unsigned *
gfi_uint32_get_data(const gfi_array *t) {
  assert(t);
  assert(t->storage.type == GFI_UINT32);
  return (unsigned*)t->storage.gfi_storage_u.data_uint32.data_uint32_val;
}

double*
gfi_double_get_data(const gfi_array* t) {
  assert(t);
  assert(t->storage.type == GFI_DOUBLE);
  return t->storage.gfi_storage_u.data_double.data_double_val;
}

int
gfi_array_is_complex(const gfi_array* t) {
  assert(t);
  if (t->storage.type == GFI_DOUBLE) 
    return t->storage.gfi_storage_u.data_double.is_complex;
  else if (t->storage.type == GFI_SPARSE)
    return t->storage.gfi_storage_u.sp.is_complex;
  else return GFI_REAL;
}

gfi_array**
gfi_cell_get_data(const gfi_array *t) {
  assert(t);
  assert(t->storage.type == GFI_CELL);
  return t->storage.gfi_storage_u.data_cell.data_cell_val;
}

gfi_object_id*
gfi_objid_get_data(const gfi_array *t) {
  assert(t);
  assert(t->storage.type == GFI_OBJID);
  return t->storage.gfi_storage_u.objid.objid_val;
}

gfi_type_id
gfi_array_get_class(const gfi_array *t) {
  assert(t);
  return (t->storage.type);
}

/* ----------------- debugging ------------------- */

const char*
gfi_array_get_class_name(const gfi_array *t) {
  assert(t);
  return gfi_type_id_name(gfi_array_get_class(t), gfi_array_is_complex(t));
}

const char *
gfi_type_id_name(gfi_type_id id, gfi_complex_flag is_complex) {
  switch (id) {
    case GFI_CHAR: return "CHAR";
    case GFI_INT32: return "INT32";
    case GFI_UINT32: return "UINT32";
    case GFI_CELL: return "CELL";
    case GFI_DOUBLE: return is_complex == GFI_REAL ? "DOUBLE" : "DOUBLE COMPLEX";
    case GFI_SPARSE: return is_complex == GFI_REAL ? "SPARSE" : "SPARSE COMPLEX";
    case GFI_OBJID: return "GETFEM OBJECT ID";
    default: return "UNKNOWN..";
  }
}

#define PRINT_ARR(format,separ,len,maxprintlen,nbperline,p) { \
int i; for (i=0; i < len && i < maxprintlen; ++i) { \
if ((i+1) % nbperline == 0) printf("\n"); else if (i>0) printf("%s",separ); \
printf(format,p[i]); } if (i<len) printf("...");}

void gfi_array_print_(gfi_array *t, int lev) {
  unsigned int i;
  if (t == NULL) { printf("NULL array ...\n"); return; }
  for (i=0; i < lev; ++i) printf("  ");
  printf("dim : "); 
  for (i=0; i < t->dim.dim_len; ++i) printf("%s%d",i>0?"x":"", t->dim.dim_val[i]);
  printf(" of %s, content={", gfi_array_get_class_name(t));
  switch (t->storage.type) {
  case GFI_CHAR: {    
    PRINT_ARR("%c","",t->storage.gfi_storage_u.data_char.data_char_len,400,80,t->storage.gfi_storage_u.data_char.data_char_val);
  } break;
  case GFI_INT32: {
    PRINT_ARR("%4d",", ",t->storage.gfi_storage_u.data_int32.data_int32_len,60,15,t->storage.gfi_storage_u.data_int32.data_int32_val);
  } break;
  case GFI_UINT32: {
    PRINT_ARR("%4d",", ",t->storage.gfi_storage_u.data_uint32.data_uint32_len,60,15,t->storage.gfi_storage_u.data_uint32.data_uint32_val);
  } break;
  case GFI_CELL: {
    int i;
    printf("\n");
    for (i=0; i < t->storage.gfi_storage_u.data_cell.data_cell_len; ++i) {
      gfi_array_print_(t->storage.gfi_storage_u.data_cell.data_cell_val[i],lev+1);
    }
    printf("\n"); for (i=0; i < lev; ++i) printf("  ");
  } break;
  case GFI_DOUBLE: {
    PRINT_ARR("%8g",", ",t->storage.gfi_storage_u.data_double.data_double_len,40,10,t->storage.gfi_storage_u.data_double.data_double_val);    
  } break;
  case GFI_SPARSE: {
    printf("\n"); for (i=0; i < lev+1; ++i) printf("  "); printf("ir="); 
    PRINT_ARR("%4d",", ",t->storage.gfi_storage_u.sp.ir.ir_len,15,16,t->storage.gfi_storage_u.sp.ir.ir_val);
    printf("\n"); for (i=0; i < lev+1; ++i) printf("  "); printf("jc="); 
    PRINT_ARR("%4d",", ",t->storage.gfi_storage_u.sp.jc.jc_len,15,16,t->storage.gfi_storage_u.sp.jc.jc_val);
    printf("\n"); for (i=0; i < lev+1; ++i) printf("  "); printf("pr="); 
    PRINT_ARR("%8g",", ",t->storage.gfi_storage_u.sp.pr.pr_len,15,8,t->storage.gfi_storage_u.sp.pr.pr_val);
    printf("\n"); for (i=0; i < lev; ++i) printf("  ");
  } break;
  case GFI_OBJID: {
    int i;
    printf("cid,id=[");
    for (i=0; i < gfi_array_nb_of_elements(t); ++i) printf("%s{%d,%d}", i?", ":"", gfi_objid_get_data(t)[i].cid, gfi_objid_get_data(t)[i].id);
    printf("]\n");
  } break;
  default: {
    printf("internal error"); return; 
  } break;
  }
  printf("}\n");
}

void 
gfi_array_print(gfi_array *t) {
  gfi_array_print_(t,0);
}

static int cancel_flag = 0;

int is_cancel_flag_set() {
  return cancel_flag;
}

/* set by ctrl-c handlers when a long computation should be interrupted
   (typically mdbrick.solve) */
void set_cancel_flag(int v) {
  cancel_flag = v;
}
