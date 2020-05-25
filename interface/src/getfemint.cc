/*===========================================================================

 Copyright (C) 2001-2020 Y. Renard, J. Pommier.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

#include <getfem/getfem_mesh_slice.h>
#include <getfem/getfem_global_function.h>
#include <getfem/getfem_mat_elem_type.h>
#include <getfem/getfem_mesher.h>
#include <getfem/getfem_continuation.h>
#include <getfem/getfem_mesh_level_set.h>
#include <getfem/getfem_im_data.h>
#include <getfem/getfem_models.h>
#include <getfemint_misc.h>
#include <getfemint_workspace.h>
#include <getfemint_precond.h>
#include <getfemint_gsparse.h>

//#ifdef MAINTAINER_MODE
# ifdef HAVE_CSIGNAL
#  include <csignal>
# else
#  include <signal.h>
# endif
//#endif
#include <exception>

namespace getfemint {

  void array_dimensions::assign_dimensions(const gfi_array *mx) {
    sz = gfi_array_nb_of_elements(mx);
    ndim_ = gfi_array_get_ndim(mx);
    const int *d = gfi_array_get_dim(mx);
    for (unsigned i=0; i < ndim_; ++i)
      if (i < ARRAY_DIMENSIONS_MAXDIM) sizes_[i] = d[i];
      else sizes_[ARRAY_DIMENSIONS_MAXDIM-1] *= d[i];
  }

  /* append the dimensions [d0,d0+n[ of other, but ignore the first dimension
     if the array is a row matrix in matlab */
  size_type array_dimensions::push_back(const array_dimensions& other,
                                        unsigned d0, unsigned n,
                                        bool matlab_row_matrix_is_a_vector) {
    size_type qqdim = 1;
    for (unsigned i=d0; i < d0+n; ++i) {
      if (i || !matlab_row_matrix_is_a_vector ||
          !(!config::has_1D_arrays() && other.ndim() == 2 && other.dim(0) == 1))
        push_back(other.dim(i));
      qqdim *= other.dim(i);
    }
    return qqdim;
  }

  void array_dimensions::reshape(unsigned n_, unsigned m_, unsigned p_) {
    if (sz != n_*m_*p_) THROW_INTERNAL_ERROR;
    ndim_ = 3; sizes_[0] = n_; sizes_[1] = m_; sizes_[2] = p_;
  }

  void array_dimensions::opt_transform_col_vect_into_row_vect() {
    if (ndim() == 1 && !config::has_1D_arrays()) {
      ndim_ = 2;
      sizes_[1] = sizes_[0];
      sizes_[0] = 1;
    }
  }

  static std::string dim_of_gfi_array(const gfi_array *t) {
    std::stringstream ss;
    for (size_type i=0; i<static_cast<size_type>(gfi_array_get_ndim(t)); ++i) {
      if (i) ss << "x";
      ss << gfi_array_get_dim(t)[i];
    }
    return ss.str();
  }

  const sub_index &sub_index::check_range(size_type n) const {
    /*for (gmm::unsorted_sub_index::const_iterator it = begin(); it != end(); ++it) {
      if (*it >= n) {
        THROW_BADARG("wrong matrix sub index: " << *it + config::base_index() << " not in range [" <<
                     config::base_index() << ".." << n-1 + config::base_index() << "]");
      }
      }*/
    if (last() >= n)
      THROW_BADARG("wrong matrix sub index: " << last() + config::base_index() << " not in range [" <<
                   config::base_index() << ".." << n-1 + config::base_index() << "]");
    return *this;
  }

  double
  mexarg_in::to_scalar_(bool isint) {
    double dv;
    if (gfi_array_nb_of_elements(arg) != 1) {
      THROW_BADARG("Argument " << argnum <<
                   " has dimensions " << dim_of_gfi_array(arg) <<
                   " but a [1x1] " << std::string(isint ? "integer" : "scalar") <<
                   " was expected");
    }
    switch (gfi_array_get_class(arg)) {
    case GFI_DOUBLE: {
      if (gfi_array_is_complex(arg))
        THROW_BADARG("Argument " << argnum <<
                     " was expected to be a REAL number and we got a COMPLEX number!");
      dv = *gfi_double_get_data(arg);
    } break;
    case GFI_INT32: {
      dv = (double)( *((dal::int32_type*)gfi_int32_get_data(arg)));
    } break;
    case GFI_UINT32: {
      dv = (double)( *((dal::int32_type*)gfi_uint32_get_data(arg)));
    } break;
    default: {
      THROW_BADARG("Argument " << argnum << " of class " << gfi_array_get_class_name(arg) <<
                   " is not a scalar value");
    } break;
    }
    return dv;
  }

  double
  mexarg_in::to_scalar(double minval, double maxval) {
    double dv = to_scalar_(false);
    if (dv < minval || dv > maxval) {
      THROW_BADARG("Argument " << argnum <<
                   " is out of bounds : " << dv << " not in " <<
                   "[" << minval << "..." << maxval << "]");
    }
    return dv;
  }

  complex_type
  mexarg_in::to_scalar(complex_type) {
    if (gfi_array_nb_of_elements(arg) != 1) {
      THROW_BADARG("Argument " << argnum <<
                   " has dimensions " << dim_of_gfi_array(arg) <<
                   " but a [1x1] complex number was expected");
    }
    garray<complex_type> g = to_garray(complex_type());
    return g[0];
  }

  bool mexarg_in::is_integer() {
    if (gfi_array_nb_of_elements(arg) != 1) return false;
    if (is_complex()) return false;
    switch (gfi_array_get_class(arg)) {
      case GFI_DOUBLE: {
        double dv = *gfi_double_get_data(arg);
        if (dv != int(dv)) return false;
        return true;
      } break;
      case GFI_INT32:
      case GFI_UINT32: return true; break;
      default: return false;
    }
    return false;
  }

  bool mexarg_in::is_bool() {
    if (gfi_array_nb_of_elements(arg) != 1) return false;
    if (is_complex()) return false;

    switch (gfi_array_get_class(arg)) {
      case GFI_UINT32: {
        unsigned uv = *gfi_uint32_get_data(arg);
        if (uv > 1.0) return false;
        return true;
      } break;
      case GFI_INT32: {
        int iv = *gfi_int32_get_data(arg);
        if (iv < 0.0 || iv > 1.0) return false;
        return true;
      } break;
      case GFI_DOUBLE: {
        double dv = *gfi_double_get_data(arg);
        if (dv != int(dv) || dv < 0.0 || dv > 1.0) return false;
        return true;
      } break;
      default: return false;
    }
    return false;
  }

  bool
  mexarg_in::to_bool() {
    double dv = to_scalar_(true);
    if (dv != floor(dv) || dv < 0.0 || dv > 1.0) {
      THROW_BADARG("Argument " << argnum <<
                   " is not an bool value");
    }
    return (bool)dv;
  }

  int
  mexarg_in::to_integer(int minval, int maxval) {
    double dv = to_scalar_(true);
    if (dv != floor(dv)) {
      THROW_BADARG("Argument " << argnum <<
                   " is not an integer value");
    }
    if (dv < minval || dv > maxval) {
      THROW_BADARG("Argument " << argnum <<
                   " is out of bounds : " << dv << " not in " <<
                   "[" << minval << "..." << maxval << "]");
    }
    return (int)dv;
  }

  size_type
  mexarg_in::to_convex_number(const getfem::mesh &m) {
    int cv = to_integer(config::base_index(), INT_MAX) - config::base_index();
    if (!m.convex_index().is_in(cv))
      THROW_BADARG( "Convex " << cv << " is not part of the mesh");
    return cv;
  }

  short_type
  mexarg_in::to_face_number(short_type nbf) {
    int f = to_integer(config::base_index(),
                       int(config::base_index()+nbf-1));
    f -= config::base_index();
    return short_type(f);
  }

  bool
  mexarg_in::is_object_id(id_type *pid, id_type *pcid) const {
    if (gfi_array_get_class(arg) == GFI_OBJID &&
        gfi_array_nb_of_elements(arg) == 1) {
      if (pid)  *pid  = gfi_objid_get_data(arg)->id;
      if (pcid)  *pcid  = gfi_objid_get_data(arg)->cid;
      return true;
    }
    return false;
  }

  bool mexarg_in::is_sparse()
  { return (gfi_array_get_class(arg) == GFI_SPARSE || is_spmat_object(*this)); }

  bool
  mexarg_in::is_complex() {
    if (!is_spmat_object(*this)) return gfi_array_is_complex(arg);
    else return to_spmat_object(*this)->is_complex();
  }

  std::string
  mexarg_in::to_string() {
    /* string => row vector. */
    if (!is_string())
      THROW_BADARG("Argument " << argnum << " must be a string.");
    return std::string(gfi_char_get_data(arg), gfi_array_nb_of_elements(arg));
  }

  id_type
  mexarg_in::to_object_id(id_type *pid, id_type *pcid) {
    id_type id,cid;
    if (!is_object_id(&id, &cid)) {
      THROW_BADARG("wrong type for argument " << argnum << ": expecting a "
		   "getfem object, got a " << gfi_array_get_class_name(arg));
    }
    if (pid) *pid = id;
    if (pcid) *pcid = cid;
    return id;
  }


  // const getfem::mesh *
  // mexarg_in::to_const_mesh(id_type& mid) {
  //   id_type id, cid;
  //   to_object_id(&id,&cid);
  //   if (cid != MESHFEM_CLASS_ID && cid != MESH_CLASS_ID && cid != MESHIM_CLASS_ID
  //       && cid != MESHIMDATA_CLASS_ID) {
  //     THROW_BADARG("argument " << argnum << " should be a mesh or mesh_fem "
  //                  "or mesh_im or mesh_im_data descriptor, its class is "
  //                  << name_of_getfemint_class_id(cid));
  //   }
  //   const getfem::mesh *mesh = NULL;
  //   getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
  //   if (object_is_mesh(o)) {
  //     mid = id;
  //     mesh = &object_to_mesh(o)->mesh();
  //   } else if (dynamic_cast<getfem::mesh_fem *>(o)) {
  //     mesh = &dynamic_cast<getfem::mesh_fem *>(o)->linked_mesh();
  //     mid = workspace().object(mesh); // A arranger avec ce qui précède...
  //     if (mid == id_type(-1)) THROW_INTERNAL_ERROR;
  //   } else if (dynamic_cast<getfem::mesh_im *>(o)) {
  //     // o n'est pas bon pour le moment. Le sera quand il viendra du workspace
  //     mesh = &dynamic_cast<getfem::mesh_im *>(o)->linked_mesh();
  //     mid = workspace().object(mesh); // A arranger avec ce qui précède...
  //     if (mid == id_type(-1)) THROW_INTERNAL_ERROR;
  //   } else if (dynamic_cast<getfem::im_data *>(o)) {
  //     // o n'est pas bon pour le moment. Le sera quand il viendra du workspace
  //     mesh =&dynamic_cast<getfem::im_data *>(o)->linked_mesh_im().linked_mesh();
  //     mid = workspace().object(mesh); // A arranger avec ce qui précède...
  //     if (mid == id_type(-1)) THROW_INTERNAL_ERROR;
  //   } else THROW_INTERNAL_ERROR;
  //   return mesh;
  // }

  getfem::mesh_region
  mexarg_in::to_mesh_region() {
    if (gfi_array_get_class(arg) != GFI_INT32 &&
        gfi_array_get_class(arg) != GFI_UINT32 &&
        gfi_array_get_class(arg) != GFI_DOUBLE) {
      THROW_BADARG("expected a mesh region!");
    }
    iarray v = to_iarray();
    return getfemint::to_mesh_region(v);
  }

  void mexarg_in::check_dimensions(const array_dimensions &v, int expected_dim) {
    if (v.getn() != 1 && v.getm() != 1 && v.size() != 0) {
      THROW_BADARG("Argument " << argnum <<
                   " should be a vector, not a matrix");
    }
    if (expected_dim != -1 && (int)v.size() != expected_dim) {
      THROW_BADARG("Argument " << argnum <<
                   " has wrong dimensions: expected " << expected_dim <<
                   ", found " << v.size());
    }
  }

  void mexarg_in::check_dimensions(array_dimensions &v, int expected_m, int expected_n, int expected_p, int expected_q) {
    /* check for preference for 'row' vectors */
    if (expected_m == -2 && expected_n == -1 && v.size() == v.getm()) {
      v.reshape(1,v.getm());
    }
    if (expected_m >= 0 && (int)v.getm() != expected_m)
      THROW_BADARG("Argument " << argnum <<
                   " has a wrong number of rows (" << v.getm() <<
                   ") , " << expected_m << " rows were expected");
    if (expected_n >= 0 && (int)v.getn() != expected_n)
      THROW_BADARG("Argument " << argnum <<
                   " has a wrong number of columns (" << v.getn() <<
                   ") , " << expected_n << " columns were expected");
    if (expected_p >= 0 && (int)v.getp() != expected_p)
      THROW_BADARG("Argument " << argnum <<
                   " was expected to be a three-dimensional array, with " <<
                   expected_p << " elements in its third dimension (got " <<
                   v.getp() << ")");
    if (expected_q >= 0 && (int)v.getq() != expected_q)
      THROW_BADARG("Argument " << argnum <<
                   " was expected to be a four-dimensional array, with " <<
                   expected_q << " elements in its fourth dimension (got " <<
                   v.getq() << ")");
  }

  mexarg_in &mexarg_in::check_trailing_dimension(int expected_dim) {
    int nd = gfi_array_get_ndim(arg);
    int d;
    if (nd == 0) d = 1; else d = gfi_array_get_dim(arg)[nd-1];
    if (d != expected_dim) {
      array_dimensions ad(arg);
      std::string tip_of_the_day;
      if (nd == 2 && int(ad.dim(0)) == expected_dim)
        tip_of_the_day = "\n You should probably transpose your array..";
      THROW_BADARG("The trailing dimension of argument " << argnum
                   << " (an array of size " << ad << ")" <<
                   " has " << d << " elements, " << expected_dim << " were expected"
                   << tip_of_the_day);
    }
    return *this;
  }

  darray
  mexarg_in::to_darray() {
    if(gfi_array_is_complex(arg) ||
       (gfi_array_get_class(arg) != GFI_DOUBLE &&
        gfi_array_get_class(arg) != GFI_INT32 &&
        gfi_array_get_class(arg) != GFI_UINT32)) {
      THROW_BADARG("Argument " << argnum <<
                   " should be a DOUBLE REAL data array");
    }
    return darray(arg);
  }

  /* check that the supplied argument IS a vector,
     with the right number of elements */
  darray
  mexarg_in::to_darray(int expected_dim) {
    darray v = to_darray();
    check_dimensions(v, expected_dim);
    return v;
  }

  /* check that the supplied array has a good number of rows
     and/or a good number of columns (and or a good number of third dimension) */
  darray
  mexarg_in::to_darray(int expected_m, int expected_n,
                       int expected_p, int expected_q) {
    darray v = to_darray();
    check_dimensions(v, expected_m, expected_n, expected_p, expected_q);
    return v;
  }

  rcarray
  mexarg_in::to_rcarray() {
    if((gfi_array_get_class(arg) != GFI_DOUBLE &&
        gfi_array_get_class(arg) != GFI_INT32 &&
        gfi_array_get_class(arg) != GFI_UINT32)) {
      THROW_BADARG("Argument " << argnum <<
                   " should be a DOUBLE REAL or COMPLEX data array");
    }
    return rcarray(arg);
  }

  /* check that the supplied argument IS a vector,
     with the right number of elements */
  rcarray
  mexarg_in::to_rcarray(int expected_dim) {
    rcarray v = to_rcarray();
    check_dimensions(v.sizes(), expected_dim);
    return v;
  }

  /* check that the supplied array has a good number of rows
     and/or a good number of columns (and or a good number of third dimension) */
  rcarray
  mexarg_in::to_rcarray(int expected_m, int expected_n,
                        int expected_p, int expected_q) {
    rcarray v = to_rcarray();
    check_dimensions(v.sizes(), expected_m, expected_n,
                     expected_p, expected_q);
    return v;
  }

  carray
  mexarg_in::to_carray() {
    if(gfi_array_get_class(arg) != GFI_DOUBLE &&
       gfi_array_get_class(arg) != GFI_INT32 &&
       gfi_array_get_class(arg) != GFI_UINT32) {
      THROW_BADARG("Argument " << argnum <<
                   " should be a DOUBLE COMPLEX data array");
    }
    return carray(arg);
  }

  /* check that the supplied argument IS a vector,
     with the right number of elements */
  carray
  mexarg_in::to_carray(int expected_dim) {
    carray v = to_carray();
    check_dimensions(v, expected_dim);
    return v;
  }

  /* check that the supplied array has a good number of rows
     and/or a good number of columns (and or a good number of third dimension) */
  carray
  mexarg_in::to_carray(int expected_m, int expected_n,
                       int expected_p, int expected_q) {
    carray v = to_carray();
    check_dimensions(v, expected_m, expected_n, expected_p, expected_q);
    return v;
  }

  iarray
  mexarg_in::to_iarray() {
    //check_int_values(INT_MIN,INT_MAX);
    if (gfi_array_get_class(arg) == GFI_INT32 ||
        gfi_array_get_class(arg) == GFI_UINT32) {
      return iarray(arg);
    } else if (gfi_array_get_class(arg) == GFI_DOUBLE) {
      darray v(arg);
      //cout << "conversion DOUBLE -> INT for array of size " << v.getm() << "x" << v.getn() << "x" << v.getp() << "\n";
      iarray ia(new int[v.size()], v.size());  ia.assign_dimensions(arg);
      for (unsigned i=0; i < v.size(); i++) {
        ia[i] = int(v[i]);
        if (ia[i] != v[i]) {
          THROW_BADARG("Argument " << argnum <<
                       " should be a DOUBLE REAL data array containing only "
                       "INTEGER values --- at index " << i+config::base_index() <<
                       " the scalar value " << v[i] << " was found");
        }
      }
      return ia;
    } else
      THROW_BADARG("Argument " << argnum <<
                   " should be an INTEGER data array");

  }

  iarray
  mexarg_in::to_iarray(int expected_dim) {
    //check_int_values();
    iarray v = to_iarray();
    check_dimensions(v, expected_dim);
    return v;
  }

  iarray
  mexarg_in::to_iarray(int expected_m, int expected_n,
                       int expected_p, int expected_q) {
    //check_int_values();
    iarray v = to_iarray(); //(arg);
    check_dimensions(v, expected_m, expected_n,
                     expected_p, expected_q);
    return v;
  }

  getfem::base_node
  mexarg_in::to_base_node(int expected_dim) {
    darray w = to_darray(expected_dim,1);
    getfem::base_node bn(w.size()); std::copy(w.begin(), w.end(), bn.begin());
    return bn;
  }

  /* get a (native only) sparse matrix */
  void mexarg_in::to_sparse(gf_real_sparse_csc_const_ref& M) {
    if (gfi_array_get_class(arg) != GFI_SPARSE) {
      THROW_BADARG("Argument " << argnum <<
		   " was expected to be a sparse matrix");
    }
    if (is_complex()) {
      THROW_BADARG("Argument " << argnum <<
		   " cannot be a complex sparse matrix");
    }
    assert(gfi_array_get_ndim(arg)==2);
    M = gf_real_sparse_csc_const_ref(gfi_sparse_get_pr(arg),
				     gfi_sparse_get_ir(arg),
				     gfi_sparse_get_jc(arg),
                                     gfi_array_get_dim(arg)[0],
				     gfi_array_get_dim(arg)[1]);
  }

  void mexarg_in::to_sparse(gf_cplx_sparse_csc_const_ref& M) {
    if (gfi_array_get_class(arg) != GFI_SPARSE) {
      THROW_BADARG("Argument " << argnum <<
		   " was expected to be a sparse matrix");
    }
    if (!is_complex()) {
      THROW_BADARG("Argument " << argnum << " cannot be a real sparse matrix");
    }
    assert(gfi_array_get_ndim(arg)==2);
    M = gf_cplx_sparse_csc_const_ref((complex_type*)gfi_sparse_get_pr(arg),
				     gfi_sparse_get_ir(arg),
				     gfi_sparse_get_jc(arg),
                                     gfi_array_get_dim(arg)[0],
				     gfi_array_get_dim(arg)[1]);
  }

  /* get a (native or getfem) sparse matrix */
  std::shared_ptr<gsparse> mexarg_in::to_sparse() {
    if (gfi_array_get_class(arg) == GFI_SPARSE) {
      return std::make_shared<gsparse>(arg);
    } else {
      id_type id,cid;
      to_object_id(&id,&cid);
      if (cid != SPMAT_CLASS_ID)
        THROW_BADARG("Argument " << argnum <<
		     " was expected to be a sparse matrix");
      auto gsp=workspace().shared_pointer(id,name_of_getfemint_class_id(cid));
      auto gsp2= std::const_pointer_cast<gsparse>
	(std::dynamic_pointer_cast<const gsparse>(gsp));
      GMM_ASSERT1(gsp2.get(), "Internal error");
      return gsp2;
    }
  }


  /* converts the gfi_array into a bit vector , shift all its values by
     'shift' and checking that they are a subset of 'subsetof' ( if the pointer
     is non-nul ) */
  dal::bit_vector
  mexarg_in::to_bit_vector(const dal::bit_vector *subsetof, int shift) {
    dal::bit_vector bv;
    iarray v = to_iarray();
    for (unsigned i = 0; i < v.size(); i++) {
      int idx = int(v[i]) + shift;
      if (idx < 0 || idx > 1000000000) {
        THROW_BADARG("Argument " << argnum <<
                     " should only contain values greater or equal to "
                     << -shift << " ([found " << v[i] << ")");
      } else if (subsetof && !subsetof->is_in(idx)) {
        THROW_BADARG("Argument " << argnum <<
                     " is not a valid set (contains values not allowed, such as " << (int)v[i] << ")");
      }
      bv.add(idx);
    }
    return bv;
  }


  /* maybe we should check that there is no dup values in the sub_index */
  sub_index
  mexarg_in::to_sub_index() {
    iarray ia = to_iarray(-1);
    std::vector<size_type> va(ia.size());
    for (unsigned i=0; i < ia.size(); ++i) {
      va[i] = ia[i] - config::base_index();
    }
    return sub_index(va);
  }

  gfi_array *
  create_object_id(int nid, id_type *ids, id_type cid, bool not_as_a_vector) {
    gfi_array *arg;
    if (not_as_a_vector) {
      assert(nid==1);
      arg = checked_gfi_array_create_0(GFI_OBJID);
    } else {
      arg = checked_gfi_array_create_1(nid, GFI_OBJID);
    }
    for (size_type i=0; i < size_type(nid); ++i) {
      gfi_objid_get_data(arg)[i].id = ids[i];
      gfi_objid_get_data(arg)[i].cid = cid;
    }
    return arg;
  }

  void
  mexarg_out::from_object_id(id_type id, id_type cid)
  { arg = create_object_id(id, cid); }

  void
  mexarg_out::from_object_id(std::vector<id_type> ids, id_type cid)
  { arg = create_object_id(int(ids.size()), &ids[0], cid); }

  void
  mexarg_out::from_integer(int i) {
    if (config::can_return_integer()) {
      arg = checked_gfi_array_create_0(GFI_INT32);
      *gfi_int32_get_data(arg) = i;
    } else from_scalar(i);
  }

  void
  mexarg_out::from_scalar(double v) {
    arg = checked_gfi_array_create_0(GFI_DOUBLE);
    *gfi_double_get_data(arg) = v;
  }

  void
  mexarg_out::from_string(const char *s) {
    arg = checked_gfi_array_from_string(s);
  }

  void
  mexarg_out::from_bit_vector(const dal::bit_vector& bv, int shift) {
    iarray w = create_iarray_h(unsigned(bv.card()));
    size_type j = 0;
    for (dal::bv_visitor i(bv); !i.finished(); ++i) {
      w[unsigned(j++)] = int(i + shift);
    }
    if (j != bv.card()) THROW_INTERNAL_ERROR;
  }

  void
  mexarg_out::from_mesh_region(const getfem::mesh_region &region) {
    iarray w = create_iarray(2, unsigned(region.size()));
    unsigned j=0;
    for (getfem::mr_visitor i(region); !i.finished(); ++i, ++j) {
      w(0,j) = int(i.cv() + config::base_index());
      w(1,j) = i.f()  + config::base_index();
    }
  }

  /* remember that M will be erased by these functions */
  void mexarg_out::from_sparse(gf_real_sparse_by_col& M,
			       output_sparse_fmt fmt) {
     gsparse gsp;
     from_sparse(gsp.destructive_assign(M), fmt);
  }

  void mexarg_out::from_sparse(gf_cplx_sparse_by_col& M,
			       output_sparse_fmt fmt) {
    gsparse gsp;
    from_sparse(gsp.destructive_assign(M), fmt);
  }

  void mexarg_out::from_sparse(gsparse& M, output_sparse_fmt fmt) {
    if (fmt == USE_DEFAULT_SPARSE) {
      fmt = (config::prefer_native_sparse() ? USE_NATIVE_SPARSE : USE_GSPARSE);
    }
    if (fmt == USE_GSPARSE) {
      auto gsp = std::make_shared<gsparse>();
      gsp->swap(M);
      id_type id = store_spmat_object(gsp);
      from_object_id(id, SPMAT_CLASS_ID);
    } else {
      M.to_csc();
      size_type nnz = M.nnz();
      size_type ni = M.nrows(), nj = M.ncols();
      arg = checked_gfi_create_sparse(int(ni), int(nj), int(nnz),
                                      M.is_complex() ? GFI_COMPLEX : GFI_REAL);
      assert(arg != NULL);
      double *pr;
      unsigned *ir, *jc;
      pr = gfi_sparse_get_pr(arg); assert(pr != NULL);
      ir = gfi_sparse_get_ir(arg); assert(ir != NULL);
      jc = gfi_sparse_get_jc(arg); assert(jc != NULL); /* dim == nj+1 */
      if (!M.is_complex()) {
        memcpy(pr, M.real_csc().pr, sizeof(double)*nnz);
        memcpy(ir, M.real_csc().ir, sizeof(int)*nnz);
        memcpy(jc, M.real_csc().jc, sizeof(int)*(nj+1));
      } else {
        memcpy(pr, M.cplx_csc().pr, sizeof(complex_type)*nnz);
        memcpy(ir, M.cplx_csc().ir, sizeof(int)*nnz);
        memcpy(jc, M.cplx_csc().jc, sizeof(int)*(nj+1));
      }
      M.deallocate(); // avoid nasty leak
    }
  }

  void
  mexarg_out::from_tensor(const getfem::base_tensor& t) {
    std::vector<int> tab(t.sizes().begin(), t.sizes().end());
    arg = checked_gfi_array_create(int(t.order()), &(tab.begin()[0]),
                                   GFI_DOUBLE);
    double *q = (double *)(gfi_double_get_data(arg));
    std::copy(t.begin(), t.end(), q);
  }

  carray
  mexarg_out::create_carray_h(unsigned dim) {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE, GFI_COMPLEX);
    else arg = checked_gfi_array_create_2(1,dim, GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  carray
  mexarg_out::create_carray_v(unsigned dim) {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE, GFI_COMPLEX);
    else arg = checked_gfi_array_create_2(dim, 1, GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  carray
  mexarg_out::create_carray(unsigned n,unsigned m,unsigned p) {
    int sz[3];
    sz[0] = n; sz[1] = m; sz[2] = p;
    arg = checked_gfi_array_create(3,sz,GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  carray
  mexarg_out::create_carray(unsigned dim, unsigned nbdof) {
    arg = checked_gfi_array_create_2(dim, nbdof, GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  darray
  mexarg_out::create_darray_h(unsigned dim) {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE);
    else
      arg = checked_gfi_array_create_2(1,dim, GFI_DOUBLE);
    return darray(arg);
  }

  darray
  mexarg_out::create_darray_v(unsigned dim) {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE);
    else
      arg = checked_gfi_array_create_2(dim, 1, GFI_DOUBLE);
    return darray(arg);
  }

  darray
  mexarg_out::create_darray(unsigned n,unsigned m,unsigned p) {
    int sz[3];
    sz[0] = n; sz[1] = m; sz[2] = p;
    arg = checked_gfi_array_create(3,sz,GFI_DOUBLE);
    return darray(arg);
  }

  /* creates a 'matrix' (from the matlab point of view.. for getfem
     it is still a vector, stored in fortran style) */
  darray
  mexarg_out::create_darray(unsigned dim, unsigned nbdof) {
    arg = checked_gfi_array_create_2(dim, nbdof, GFI_DOUBLE);
    return darray(arg);
  }

  iarray
  mexarg_out::create_iarray_h(unsigned dim) {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_INT32);
    else arg = checked_gfi_array_create_2(1,dim, GFI_INT32);
    return iarray(arg);
  }

  iarray
  mexarg_out::create_iarray_v(unsigned dim) {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_INT32);
    else arg = checked_gfi_array_create_2(dim, 1, GFI_INT32);
    return iarray(arg);
  }

  iarray
  mexarg_out::create_iarray(unsigned n,unsigned m,unsigned p) {
    int sz[3];
    sz[0] = n; sz[1] = m; sz[2] = p;
    arg = checked_gfi_array_create(3,sz,GFI_INT32);
    return iarray(arg);
  }

  /* creates a 'matrix' (from the matlab point of view.. for getfem
     it is still a vector, stored in fortran style) */
  iarray
  mexarg_out::create_iarray(unsigned dim, unsigned nbdof) {
    arg = checked_gfi_array_create_2(dim, nbdof, GFI_INT32);
    return iarray(arg);
  }

  darray mexarg_out::create_array(const array_dimensions &dims, double) {
    arg = checked_gfi_array_create(dims.ndim(),(const int*)dims.sizes(),GFI_DOUBLE);
    return darray(arg);
  }

  carray mexarg_out::create_array(const array_dimensions &dims, complex_type) {
    arg = checked_gfi_array_create(dims.ndim(),(const int*)dims.sizes(),GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }


  std::string cmd_normalize(const std::string& a) {
    std::string b = a;
    for (size_type i = 0; i < b.size(); ++i) {
      b[i] = char(toupper(b[i]));
      if (b[i] == '_' || b[i] == '-') b[i] = ' ';
    }
    return b;
  }

  /* very tolerant case-insensitive string comparison:
     spaces are matched with underscores.*/
  bool cmd_strmatch(const std::string& a, const char *s) {
    return
      cmd_strmatchn(a,s,unsigned(-1));
  }

  bool cmd_strmatchn(const std::string& a, const char *s, unsigned n) {
    size_type i;
    for (i=0; s[i] && i < n && i < a.size(); ++i) {
      if ((a[i] == ' ' || a[i] == '_') && (s[i] == ' ' || s[i] == '_' || s[i] == '-')) continue;
      if (toupper(a[i]) != toupper(s[i])) return false;
    }
    if (i == n || (s[i] == 0 && i == a.size()))
      return true;
    else return false;
  }

  bool check_cmd(const std::string& cmdname, const char *s,
                 const mexargs_in& in,
                 int min_argin, int max_argin) {
    if (cmd_strmatch(cmdname,s)) {
      if ((int)in.remaining() < min_argin) {
        THROW_BADARG("Not enough input arguments for command '"<<
                     cmdname << "' (got " << in.narg() <<
                     ", expected at least " << min_argin + in.narg()- in.remaining() << ")");
      }
      if ((int)in.remaining() > max_argin && max_argin != -1) {
        THROW_BADARG("Too much input arguments for command '"<<
                     cmdname << "' (got " << in.narg() <<
                     ", expected at most " << max_argin + in.narg()- in.remaining()<< ")");
      }
      return true;
    }
    return false;
  }

  bool check_cmd(const std::string& cmdname, const char *s,
                 const mexargs_out& out,
                 int min_argout, int max_argout) {
    int Index = 0;
    if (cmd_strmatch(cmdname,s)) {
      if (out.get_scilab()) {
          // For Scilab. Without any parameters, myinterface(x) is equivalent to
          // ans = myinterface(x). So, we have a minima one parameter.
          if (min_argout==0 && max_argout==0) {
              min_argout = -1;
              max_argout = -1;
          }
          Index = 1;
      }

      if (min_argout > 0 && out.narg_known() &&
          out.narg_in_range(Index, min_argout-1)) {
        THROW_BADARG("Not enough output arguments for command '"<<
                     cmdname << "' (expected at least " << min_argout << ")");
      }
      if (out.narg_known() && out.narg_known() &&
          out.narg_in_range(max_argout+1,-1) && max_argout != -1) {
        THROW_BADARG("Too much output arguments for command '"<<
                     cmdname << "' (expected at most " << max_argout << ")");
      }
      return true;
    }
    return false;
  }

  bool check_cmd(const std::string& cmdname, const char *s,
                 const mexargs_in& in, const mexargs_out& out,
                 int min_argin, int max_argin,
                 int min_argout, int max_argout) {
    return
      check_cmd(cmdname, s,  in, min_argin, max_argin) &&
      check_cmd(cmdname, s, out, min_argout, max_argout);
  }

  mexargs_in::mexargs_in(int n, const gfi_array *p[], bool use_cell_) {
    nb_arg = n; use_cell = use_cell_;
    if (!use_cell) {
      in = p; idx.add(0, n);
    } else {
      assert(n == 1);
      assert(p[0]!=0);
      if (gfi_array_get_class(p[0])!=GFI_CELL)
	THROW_BADARG("Need a argument of type list");

      // assert(gfi_array_get_class(p[0])==GFI_CELL);
      nb_arg = gfi_array_nb_of_elements(p[0]);
      in = new const gfi_array*[nb_arg];
      for (int i = 0; i < nb_arg; i++) { in[i] = gfi_cell_get_data(p[0])[i]; idx.add(i); }
    }
  }

  mexargs_in::~mexargs_in() {
    if (in && use_cell) delete[] in;
  }

  mexargs_out::mexargs_out(int n) {
    idx = 0;
    okay = 0;
    nb_arg = n;
    scilab_flag = false;
  }

  mexargs_out::~mexargs_out() {
    if (!okay) {
      for (size_type i=0; i < out.size(); ++i) {
        if (out[i]) { gfi_array_destroy(out[i]); free(out[i]); }
      }
      out.clear();
      workspace().destroy_newly_created_objects();
    } else {
      workspace().commit_newly_created_objects();
    }
  }

  /* ensure that out[idx] is valid */
  void mexargs_out::check() const {
    if (nb_arg != -1) {
      if ((idx >= nb_arg) && !(idx==0))
        GMM_ASSERT1(false, "Insufficient number of output arguments");
    }
    if (size_type(idx) >= out.size()) out.resize(idx+1, 0);
  }

  mexarg_out mexargs_out::pop() {
    check(); idx++;
    return mexarg_out(out[idx-1], idx);
  }

  void
  mexargs_out::return_packed_obj_ids(const std::vector<id_type>& id, id_type class_id) {
    std::vector<id_type> uid(id);
    std::sort(uid.begin(), uid.end());
    /* remove duplicates */
    uid.erase(std::unique(uid.begin(), uid.end()), uid.end());
    /* remove id_type(-1) */
    std::vector<id_type>::iterator it =
      std::find(uid.begin(), uid.end(), id_type(-1));
    if (it != uid.end())
      uid.erase(it);
    pop().from_object_id(uid, class_id);
    if (remaining()) {
      std::map<id_type,id_type> m;
      for (size_type i=0; i < uid.size(); ++i)
        m[uid[i]]= unsigned(i + config::base_index());
      iarray v = pop().create_iarray_h(unsigned(id.size()));
      for (size_type i=0; i < id.size(); ++i)
        v[unsigned(i)] = (id[i] != id_type(-1)) ? m[id[i]] : id_type(-1);
    }
  }

  /* Associate the class ID found in the matlab structures referencing
     getfem object to a class name which coincides with the class name
     given by matlab to the structure.
     
     IMPORTANT: Should correspond to the getfemint_class_id
                In particular, it should be in alphabetic order.
     
     To be completed when an object class is added.
  */
  const char *name_of_getfemint_class_id(id_type cid) {
    switch (cid) {
    case CONT_STRUCT_CLASS_ID:      return "gfContStruct";
    case CVSTRUCT_CLASS_ID:         return "gfCvStruct";
    case ELTM_CLASS_ID:             return "gfEltm";
    case FEM_CLASS_ID:              return "gfFem";
    case GEOTRANS_CLASS_ID:         return "gfGeoTrans";
    case GLOBAL_FUNCTION_CLASS_ID:  return "gfGlobalFunction";
    case INTEG_CLASS_ID:            return "gfInteg";
    case LEVELSET_CLASS_ID:         return "gfLevelSet";
    case MESH_CLASS_ID:             return "gfMesh";
    case MESHFEM_CLASS_ID:          return "gfMeshFem";
    case MESHIM_CLASS_ID:           return "gfMeshIm";
    case MESHIMDATA_CLASS_ID:       return "gfMeshImData";
    case MESH_LEVELSET_CLASS_ID:    return "gfMeshLevelSet";
    case MESHER_OBJECT_CLASS_ID:    return "gfMesherObject";
    case MODEL_CLASS_ID:            return "gfModel";
    case PRECOND_CLASS_ID:          return "gfPrecond";
    case SLICE_CLASS_ID:            return "gfSlice";
    case SPMAT_CLASS_ID:            return "gfSpmat";
    case POLY_CLASS_ID:             return "gfPoly";
    default :                       return "not_a_getfem_class";
    }
  }

  // Gives the class id of an object.
  // To be completed when an object class is added.
  id_type class_id_of_object(const dal::pstatic_stored_object &p,
			     const void **q_) {
    const void *qq; const void **q(&qq); if (q_) { q = q_; *q_ = 0; }
    if ((*q=dynamic_cast<const getfem::cont_struct_getfem_model *>(p.get())))
      return CONT_STRUCT_CLASS_ID;
    if ((*q=dynamic_cast<const bgeot::convex_structure *>(p.get())))
      return CVSTRUCT_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::mat_elem_type *>(p.get())))
      return ELTM_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::virtual_fem *>(p.get())))
      return FEM_CLASS_ID;
    if ((*q=dynamic_cast<const bgeot::geometric_trans *>(p.get())))
      return GEOTRANS_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::abstract_xy_function *>(p.get())))
      return GLOBAL_FUNCTION_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::integration_method *>(p.get())))
      return INTEG_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::level_set *>(p.get())))
      return LEVELSET_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::mesh *>(p.get())))
      return MESH_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::mesh_fem *>(p.get())))
      return MESHFEM_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::mesh_im *>(p.get())))
      return MESHIM_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::im_data *>(p.get())))
      return MESHIMDATA_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::mesh_level_set *>(p.get())))
      return MESH_LEVELSET_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::mesher_signed_distance *>(p.get())))
      return MESHER_OBJECT_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::model *>(p.get())))
      return MODEL_CLASS_ID;
    if ((*q=dynamic_cast<const gprecond_base *>(p.get())))
      return PRECOND_CLASS_ID;
    if ((*q=dynamic_cast<const getfem::stored_mesh_slice *>(p.get())))
      return SLICE_CLASS_ID;
    if ((*q=dynamic_cast<const gsparse *>(p.get())))
      return SPMAT_CLASS_ID;
    if ((*q=dynamic_cast<const getfemint_poly *>(p.get())))
      return POLY_CLASS_ID;
    return id_type(-1);
  }


  // Version of the interface functions for an object managed preferabily by
  // a raw pointer, allowing to be modified.
# define SIMPLE_RAW_POINTER_MANAGED_OBJECT(NAME, TYPE, CLASS_ID)	\
  bool is_##NAME##_object(const mexarg_in &p) {				\
    id_type id, cid;							\
    return (p.is_object_id(&id, &cid) && cid == CLASS_ID);		\
  }									\
									\
  id_type store_##NAME##_object(const std::shared_ptr<TYPE> &shp) {	\
    auto &w = workspace();						\
    id_type id = w.object((const void *)(shp.get()));			\
    if (id == id_type(-1)) {						\
      auto p =	std::dynamic_pointer_cast				\
	<const dal::static_stored_object>(shp);				\
      if (!(p.get())) THROW_INTERNAL_ERROR;				\
      id = w.push_object(p, (const void *)(shp.get()), CLASS_ID);	\
    }									\
    return id;								\
  }									\
									\
  TYPE *to_##NAME##_object(const mexarg_in &p) {			\
    id_type id, cid;							\
    if (p.is_object_id(&id, &cid) && cid == CLASS_ID) {			\
      return const_cast<TYPE *>						\
	((const TYPE *)							\
	 (workspace().object(id, name_of_getfemint_class_id(cid))));	\
    } else {								\
      THROW_BADARG("argument " << p.argnum << " should be a " <<	\
		   name_of_getfemint_class_id(CLASS_ID) <<		\
		   " descriptor, its class is "				\
		   << name_of_getfemint_class_id(cid));			\
    }									\
  }
  
  // Version of the interface functions for an object managed only by the
  // shared pointer, and assumed not to be modified.
# define SIMPLE_SHARED_POINTER_MANAGED_OBJECT(NAME, TYPE, CLASS_ID)	\
  bool is_##NAME##_object(const mexarg_in &p) {				\
    id_type id, cid;							\
    return (p.is_object_id(&id, &cid) && cid == CLASS_ID);		\
  }									\
									\
  id_type store_##NAME##_object(const std::shared_ptr<const TYPE> &shp)	\
  {									\
    auto &w = workspace();						\
    id_type id = w.object((const void *)(shp.get()));			\
    if (id == id_type(-1)) {						\
      auto p =								\
       std::dynamic_pointer_cast<const dal::static_stored_object>(shp); \
      if (!(p.get())) THROW_INTERNAL_ERROR;				\
      id = w.push_object(p, (const void *)(shp.get()), CLASS_ID);	\
    }									\
    return id;								\
  }									\
									\
  std::shared_ptr<const TYPE> to_##NAME##_object(const mexarg_in &p) {	\
    id_type id, cid;							\
    if (p.is_object_id(&id, &cid) && cid == CLASS_ID) {			\
      return std::dynamic_pointer_cast<const TYPE>			\
	(workspace().shared_pointer(id,					\
				     name_of_getfemint_class_id(cid)));	\
    } else {								\
      THROW_BADARG("argument " << p.argnum << " should be a " <<	\
		   name_of_getfemint_class_id(CLASS_ID) <<		\
		   " descriptor, its class is "				\
		   << name_of_getfemint_class_id(cid));			\
    }									\
  }

  // Functions for CONT_STRUCT_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(cont_struct,
				    getfem::cont_struct_getfem_model,
				    CONT_STRUCT_CLASS_ID)
  
  // Functions for CVSTRUCT_CLASS_ID
  SIMPLE_SHARED_POINTER_MANAGED_OBJECT(cvstruct, bgeot::convex_structure,
				       CVSTRUCT_CLASS_ID)
  
  // Functions for ELTM_CLASS_ID
  SIMPLE_SHARED_POINTER_MANAGED_OBJECT(eltm,
				       getfem::mat_elem_type,
				       ELTM_CLASS_ID)

  // Functions for FEM_CLASS_ID
  SIMPLE_SHARED_POINTER_MANAGED_OBJECT(fem, getfem::virtual_fem,
				       FEM_CLASS_ID)

  // Functions for GEOTRANS_CLASS_ID
  SIMPLE_SHARED_POINTER_MANAGED_OBJECT(geotrans, bgeot::geometric_trans,
				       GEOTRANS_CLASS_ID)

  // Functions for GLOBAL_FUNCTION_CLASS_ID
  SIMPLE_SHARED_POINTER_MANAGED_OBJECT(global_function,
				       getfem::abstract_xy_function,
				       GLOBAL_FUNCTION_CLASS_ID)

  // Functions for INTEG_CLASS_ID
  SIMPLE_SHARED_POINTER_MANAGED_OBJECT(integ, getfem::integration_method,
				       INTEG_CLASS_ID)

  // Functions for LEVELSET_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(levelset, getfem::level_set,
				    LEVELSET_CLASS_ID)

  // Functions for MESH_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(mesh, getfem::mesh, MESH_CLASS_ID)
  
  bool has_mesh_object(const mexarg_in &p) {
    return is_mesh_object(p) || is_meshfem_object(p) || is_meshim_object(p) ||
      is_meshimdata_object(p) || is_mesh_levelset_object(p);
  }

  getfem::mesh *extract_mesh_object(const mexarg_in &p) {
    id_type id, cid;
    if (p.is_object_id(&id, &cid)) {
      switch (cid) {
      case MESH_CLASS_ID: return to_mesh_object(p);
      case MESH_LEVELSET_CLASS_ID:
	return const_cast<getfem::mesh *>
	  (&to_mesh_levelset_object(p)->linked_mesh());
      case MESHIMDATA_CLASS_ID:
	return const_cast<getfem::mesh *>
	  (&to_meshimdata_object(p)->linked_mesh());
	case MESHIM_CLASS_ID:
	  return const_cast<getfem::mesh *>
	    (&to_meshim_object(p)->linked_mesh());
      case MESHFEM_CLASS_ID:
	return const_cast<getfem::mesh *>
	  (&to_meshfem_object(p)->linked_mesh());
      default:THROW_BADARG("This object do not have a mesh");
      }
    } else THROW_BADARG("Not a getfem object");
  }

  // Functions for MESHFEM_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(meshfem, getfem::mesh_fem, MESHFEM_CLASS_ID)

  // Functions for MESHIM_CLASS_ID
   SIMPLE_RAW_POINTER_MANAGED_OBJECT(meshim, getfem::mesh_im, MESHIM_CLASS_ID)

  // Functions for MESHIMDATA_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(meshimdata, getfem::im_data,
				    MESHIMDATA_CLASS_ID)

  // Functions for MESH_LEVELSET_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(mesh_levelset, getfem::mesh_level_set,
				    MESH_LEVELSET_CLASS_ID)

  // Functions for MESHER_OBJECT_CLASS_ID
  SIMPLE_SHARED_POINTER_MANAGED_OBJECT(mesher, getfem::mesher_signed_distance,
				       MESHER_OBJECT_CLASS_ID)
  
  // Functions for MODEL_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(model, getfem::model, MODEL_CLASS_ID)

  // Functions for PRECOND_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(precond, gprecond_base, PRECOND_CLASS_ID)

  // Functions for PRECOND_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(slice, getfem::stored_mesh_slice,
				    SLICE_CLASS_ID)

  // Functions for SPMAT_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(spmat, gsparse, SPMAT_CLASS_ID)

  // Functions for POLY_CLASS_ID
  SIMPLE_RAW_POINTER_MANAGED_OBJECT(poly, getfemint_poly, POLY_CLASS_ID)

} /* namespace getfemint */
