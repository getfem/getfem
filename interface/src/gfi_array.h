/* -*- c++ -*- (enables emacs c++ mode)
===========================================================================

 Copyright (C) 2006-2008 Yves Renard, Julien Pommier.

 This file is a part of GETFEM++

 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you may use this file as part of a free software
 library without restriction.  Specifically, if other files instantiate
 templates or use macros or inline functions from this file, or you compile
 this file and link it with other files to produce an executable, this
 file does not by itself cause the resulting executable to be covered by
 the GNU General Public License.  This exception does not however
 invalidate any other reasons why the executable file might be covered by
 the GNU General Public License.

 ===========================================================================*/

#ifndef GFI_ARRAY
#define GFI_ARRAY

#include <sys/types.h>
#ifdef USE_RPC
# include <rpc/types.h>
#else
# ifndef __u_char_defined
  typedef unsigned int u_int;
#  endif
#endif

/* remember to update also the list in getfem.base.py !! */
typedef enum { MESH_CLASS_ID, MESHFEM_CLASS_ID, MESHIM_CLASS_ID, 
	       MDBRICK_CLASS_ID, MDSTATE_CLASS_ID, 
	       GEOTRANS_CLASS_ID, 
	       FEM_CLASS_ID, INTEG_CLASS_ID, 
	       ELTM_CLASS_ID, CVSTRUCT_CLASS_ID,
	       POLY_CLASS_ID, SLICE_CLASS_ID, 
	       GSPARSE_CLASS_ID, PRECOND_CLASS_ID, 
               LEVELSET_CLASS_ID, MESH_LEVELSET_CLASS_ID, 
	       GETFEMINT_NB_CLASS } getfemint_class_id;

#ifdef __cplusplus
extern "C" {
#endif

  /* extrait du gfm_rpc.h genere par rpcgen */
enum gfi_type_id {
	GFI_INT32 = 0,
	GFI_UINT32 = 1,
	GFI_DOUBLE = 2,
	GFI_CHAR = 4,
	GFI_CELL = 5,
	GFI_OBJID = 6,
	GFI_SPARSE = 7
};
typedef enum gfi_type_id gfi_type_id;

struct gfi_object_id {
	int id;
	int cid;
};
typedef struct gfi_object_id gfi_object_id;

struct gfi_sparse {
	struct {
		u_int ir_len;
		int *ir_val;
	} ir;
	struct {
		u_int jc_len;
		int *jc_val;
	} jc;
	struct {
		u_int pr_len; /* == nnz*2 when is_complex == 1 */
		double *pr_val;
	} pr;
        int is_complex;
};
typedef struct gfi_sparse gfi_sparse;

typedef struct gfi_array *pgfi_array;

struct gfi_storage {
	gfi_type_id type;
	union {
		struct {
			u_int data_int32_len;
			int *data_int32_val;
		} data_int32;
		struct {
			u_int data_uint32_len;
			u_int *data_uint32_val;
		} data_uint32;
		struct {
			u_int data_double_len; /* twice the real size of the array when is_complex == 1 */
			double *data_double_val;
                        int is_complex;
		} data_double;
		struct {
			u_int data_char_len;
			char *data_char_val;
		} data_char;
		struct {
			u_int data_cell_len;
			pgfi_array *data_cell_val;
		} data_cell;
		struct {
			u_int objid_len;
			struct gfi_object_id *objid_val;
		} objid;
		struct gfi_sparse sp;
	} gfi_storage_u;
};
typedef struct gfi_storage gfi_storage;

struct gfi_array {
	struct {
		u_int dim_len;
		u_int *dim_val;
	} dim;
	gfi_storage storage;
};
typedef struct gfi_array gfi_array;

struct gfi_array_list {
	struct {
		u_int arg_len;
		gfi_array *arg_val;
	} arg;
};
typedef struct gfi_array_list gfi_array_list;

enum gfi_status {
	GFI_STATUS_OK = 0,
	GFI_STATUS_ERROR = 1
};
typedef enum gfi_status gfi_status;

struct gfi_output {
	gfi_status status;
	union {
		gfi_array_list output;
		char *errmsg;
	} gfi_output_u;
  char *infomsg;
};
typedef struct gfi_output gfi_output;

typedef enum {GFI_REAL=0, GFI_COMPLEX=1} gfi_complex_flag;

void gfi_free(void *p);
void* gfi_calloc(size_t n, size_t m);

gfi_array*
gfi_array_create(int ndim, int *dims, gfi_type_id type, gfi_complex_flag);
gfi_array*
gfi_array_create_1(int M, gfi_type_id type, gfi_complex_flag);
gfi_array*
gfi_array_create_2(int M, int N, gfi_type_id type, gfi_complex_flag);
gfi_array*
gfi_array_from_string(const char *s);
gfi_array*
gfi_create_sparse(int m, int n, int nzmax, gfi_complex_flag);
  /*gfi_array*
    gfi_create_objid(int nid, unsigned *ids, unsigned cid);*/
void 
gfi_array_destroy(gfi_array *t);
int
gfi_array_get_ndim(const gfi_array*t);
const int*
gfi_array_get_dim(const gfi_array*t);
unsigned 
gfi_array_nb_of_elements(const gfi_array *t);
unsigned int*
gfi_sparse_get_ir(const gfi_array *t);
unsigned int*
gfi_sparse_get_jc(const gfi_array *t);
double*
gfi_sparse_get_pr(const gfi_array *t);
char*
gfi_char_get_data(const gfi_array* t);
int*
gfi_int32_get_data(const gfi_array *t);
unsigned *
gfi_uint32_get_data(const gfi_array *t);
double*
gfi_double_get_data(const gfi_array* t);
int
gfi_array_is_complex(const gfi_array* t);
gfi_array**
gfi_cell_get_data(const gfi_array *t);
gfi_object_id*
gfi_objid_get_data(const gfi_array *t);
gfi_type_id
gfi_array_get_class(const gfi_array *t);
const char*
gfi_array_get_class_name(const gfi_array *t);
const char *
gfi_type_id_name(gfi_type_id id, gfi_complex_flag is_complex);
void gfi_array_print(gfi_array *t);

typedef enum { MATLAB_INTERFACE, PYTHON_INTERFACE } gfi_interface_type;

  int is_cancel_flag_set(void);
  void set_cancel_flag(int v);
#ifdef __cplusplus
}
#endif

#endif
