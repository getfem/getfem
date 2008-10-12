// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2001-2008 Y. Renard, J. Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include <map>
#include <getfemint.h>
#include <getfemint_misc.h>
#include <getfemint_workspace.h>
#include <getfemint_poly.h>
#include <getfemint_mesh.h>
#include <getfemint_mesh_slice.h>
#include <getfemint_mesh_fem.h>
#include <getfemint_mesh_im.h>
#include <getfemint_mdbrick.h>
#include <getfemint_mdstate.h>
#include <getfemint_matelemtype.h>
#include <getfemint_matelem.h>
#include <getfemint_pfem.h>
#include <getfemint_integ.h>
#include <getfemint_pgt.h>
#include <getfemint_convex_structure.h>
#include <getfemint_precond.h>
#include <getfemint_gsparse.h>
#include <getfemint_levelset.h>
#include <getfemint_mesh_levelset.h>

#include <getfemint_misc.h>
//#ifdef MAINTAINER_MODE
# ifdef HAVE_CSIGNAL
#  include <csignal>
# else 
#  include <signal.h>
# endif
//#endif
#include <exception> // NE PAS METTRE CE FICHIER EN PREMIER !!!! ou bien ennuis avec dec cxx 6.3 garantis

namespace getfemint {
  /* associate the class ID found in the matlab structures referencing
     getfem object to a class name which coincides with the class name
     given by matlab to the structure (hum..) */
  const char *name_of_getfemint_class_id(unsigned cid) {
    static const char *cname[GETFEMINT_NB_CLASS] = {
      "gfMesh", "gfMeshFem", "gfMeshIm", "gfMdBrick", "gfMdState", 
      "gfGeoTrans", 
      "gfFem", "gfInteg","gfEltm","gfCvStruct","gfPoly", "gfSlice", 
      "gfSpmat", "gfPrecond", "gfLevelSet", "gfMeshLevelSet"
    };

    if (cid >= GETFEMINT_NB_CLASS) return "not_a_getfem_class";
    else return cname[cid];
  }

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

  /*
  bool is_static_object(id_type id, id_type cid) {
    if (cid == FEM_CLASS_ID && id < 0x80000000)
      return false;
    else return
      (cid != MESHFEM_CLASS_ID && cid != MESHIM_CLASS_ID &&
       cid != MESH_CLASS_ID && 
       cid != MDBRICK_CLASS_ID && cid != MDSTATE_CLASS_ID && 
       cid != SLICE_CLASS_ID && cid != POLY_CLASS_ID && 
       cid != PRECOND_CLASS_ID && cid != GSPARSE_CLASS_ID && 
       cid != LEVELSET_CLASS_ID && cid != MESH_LEVELSET_CLASS_ID);
  }
  */

  static std::string dim_of_gfi_array(const gfi_array *t) {
    std::stringstream ss;
    for (size_type i=0; i < static_cast<size_type>(gfi_array_get_ndim(t)); ++i) {
      if (i) ss << "x"; ss << gfi_array_get_dim(t)[i];
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
  mexarg_in::to_scalar_(bool isint)
  {
    double dv;
    if (gfi_array_nb_of_elements(arg) != 1) {
      THROW_BADARG("Argument " << argnum << 
		   "has dimensions " << dim_of_gfi_array(arg) <<
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
                   " is not an scalar value");
    } break;
    }
    return dv;
  }

  double
  mexarg_in::to_scalar(double minval, double maxval)
  {
    double dv = to_scalar_(false);
    if (dv < minval || dv > maxval) {
      THROW_BADARG("Argument " << argnum << 
		" is out of bounds : " << dv << " not in [" << 
		minval << "..." << maxval << "]");
    }
    return dv;
  }

  complex_type
  mexarg_in::to_scalar(complex_type) {
    if (gfi_array_nb_of_elements(arg) != 1) {
      THROW_BADARG("Argument " << argnum << 
		   "has dimensions " << dim_of_gfi_array(arg) <<
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

  int
  mexarg_in::to_integer(int minval, int maxval)
  {
    double dv = to_scalar_(true);
    if (dv != floor(dv)) {
      THROW_BADARG("Argument " << argnum << 
		" is not an integer value");
    }
    if (dv < minval || dv > maxval) {
      THROW_BADARG("Argument " << argnum << 
		" is out of bounds : " << dv << " not in [" << 
		minval << "..." << maxval << "]");
    }
    return (int)dv;
  }

  size_type 
  mexarg_in::to_convex_number(const getfem::mesh &m) 
  {
    int cv = to_integer(config::base_index(), INT_MAX) - config::base_index();
    if (!m.convex_index().is_in(cv)) 
      THROW_BADARG( "Convex " << cv << " is not part of the mesh"); 
    return cv;
  }

  size_type 
  mexarg_in::to_face_number(size_type nbf) 
  {
    size_type f = to_integer(config::base_index(), 
                             int(config::base_index()+nbf-1));
    f -= config::base_index();
    return f;
  }

  bool
  mexarg_in::is_object_id(id_type *pid, id_type *pcid) 
  {    
    if (gfi_array_get_class(arg) == GFI_OBJID && 
        gfi_array_nb_of_elements(arg) == 1) {
      if (pid)  *pid  = gfi_objid_get_data(arg)->id;
      if (pcid)  *pcid  = gfi_objid_get_data(arg)->cid;
      return true;
    }
    return false;
  }

  bool
  mexarg_in::is_mesh() 
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == MESH_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_mesh(o));
    } else return false;
  }

  bool
  mexarg_in::is_mesh_fem() 
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == MESHFEM_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_mesh_fem(o));
    } else return false;
  }

  bool
  mexarg_in::is_mesh_im() 
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == MESHIM_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_mesh_im(o));
    } else return false;
  }

  bool 
  mexarg_in::is_mdbrick()
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == MDBRICK_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_mdbrick(o));
    } else return false;
  }

  bool 
  mexarg_in::is_mdstate()
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == MDSTATE_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_mdstate(o));
    } else return false;
  }

  bool
  mexarg_in::is_mesh_slice() 
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == SLICE_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_mesh_slice(o));
    } else return false;
  }

  bool
  mexarg_in::is_levelset() 
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == LEVELSET_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_levelset(o));
    } else return false;
  }

  bool
  mexarg_in::is_mesh_levelset() 
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == MESH_LEVELSET_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_mesh_levelset(o));
    } else return false;
  }

  bool 
  mexarg_in::is_gsparse()
  {
    id_type id, cid;
    if (is_object_id(&id, &cid) && cid == GSPARSE_CLASS_ID) {
      getfem_object *o = workspace().object(id, name_of_getfemint_class_id(cid));
      return (object_is_gsparse(o));
    } else return false;
  }

  bool
  mexarg_in::is_complex()
  { 
    if (!is_gsparse()) return gfi_array_is_complex(arg); 
    else return to_sparse()->is_complex(); 
  }

  std::string 
  mexarg_in::to_string()
  {
    /* string => row vector. */
    if (!is_string()) 
      THROW_BADARG("Argument " << argnum << " must be a string.");
    return std::string(gfi_char_get_data(arg), gfi_array_nb_of_elements(arg));
  }

  id_type
  mexarg_in::to_object_id(id_type *pid, id_type *pcid) 
  {
    id_type id,cid;
    if (!is_object_id(&id, &cid)) {
      THROW_BADARG("wrong type for argument " << argnum <<
		   ": expecting a getfem object, got a " << gfi_array_get_class_name(arg));
    }
    if (pid) *pid = id;
    if (pcid) *pcid = cid;
    return id;
  }


  /*
    check if the argument is a valid handle to an intset
    and returns it
  */
  /*
  dal::bit_vector * mexarg_in::to_intset()
  {
    getfem_object *o = workspace().object(to_object_id());
    if (!object_is_intset(o)) {
      THROW_BADARG("Argument " << argnum << " should be an intset descriptor");    }
    return &object_to_intset(o)->intset();
  }
  */


  /*
    check if the argument is a valid handle to a base_poly
    and returns it
  */
  bgeot::base_poly *
  mexarg_in::to_poly()
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != POLY_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a polynom descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    return &object_to_poly(o)->poly();
  }

  void mexarg_in::error_if_nonwritable(getfem_object *o, bool want_writeable) {
    if (want_writeable && o->is_const()) 
      THROW_BADARG("argument " << argnum << 
		   " should be a modifiable " << 
		   name_of_getfemint_class_id(o->class_id()) << 
		   ", this one is marked as read-only");
  }

  /*
    check if the argument is a valid handle to a mesh_fem,
    and returns it
  */
  getfemint_mesh_fem *
  mexarg_in::to_getfemint_mesh_fem(bool writeable)
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != MESHFEM_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a mesh_fem descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o, writeable);
    return object_to_mesh_fem(o);
  }
  getfem::mesh_fem *
  mexarg_in::to_mesh_fem() {
    return &to_getfemint_mesh_fem(true)->mesh_fem();
  }

  const getfem::mesh_fem *
  mexarg_in::to_const_mesh_fem() {
    return &to_getfemint_mesh_fem(false)->mesh_fem();
  }


  /*
    check if the argument is a valid handle to a mesh_im,
    and returns it
  */
  getfemint_mesh_im *
  mexarg_in::to_getfemint_mesh_im(bool writeable)
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != MESHIM_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a mesh_im descriptor, "
		   "its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o,writeable);
    return object_to_mesh_im(o);
  }

  getfem::mesh_im *
  mexarg_in::to_mesh_im() {
    return &to_getfemint_mesh_im(true)->mesh_im();
  }

  const getfem::mesh_im *
  mexarg_in::to_const_mesh_im() {
    return &to_getfemint_mesh_im(false)->mesh_im();
  }

  /*
    check if the argument is a valid mesh handle
    if a mesh_fem handle is given, its associated
    mesh is used
  */
  const getfem::mesh *
  mexarg_in::to_const_mesh() { id_type mid; return to_const_mesh(mid); }

  const getfem::mesh *
  mexarg_in::to_const_mesh(id_type& mid)
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != MESHFEM_CLASS_ID && cid != MESH_CLASS_ID && cid != MESHIM_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a mesh or mesh_fem "
		   "or mesh_im descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    const getfem::mesh *mesh = NULL;
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    if (object_is_mesh(o)) {
      mid = id;
      mesh = &object_to_mesh(o)->mesh();
    } else if (object_is_mesh_fem(o)) {
      mid = object_to_mesh_fem(o)->linked_mesh_id();
      mesh = &object_to_mesh_fem(o)->mesh_fem().linked_mesh();
    } else if (object_is_mesh_im(o)) {
      mid = object_to_mesh_im(o)->linked_mesh_id();
      mesh = &object_to_mesh_im(o)->mesh_im().linked_mesh();
    } else THROW_INTERNAL_ERROR;
    return mesh;
  }


  /*
    check if the argument is a valid mesh handle
    (the returned arg is not const, so we can't use the mesh from mesh_fem objects)
  */
  getfemint_mesh *
  mexarg_in::to_getfemint_mesh(bool writeable)
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != MESH_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a mesh descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o,writeable);
    return object_to_mesh(o);
  }

  getfem::mesh *
  mexarg_in::to_mesh() {
    return &to_getfemint_mesh()->mesh();
  }

  getfemint_mdbrick *
  mexarg_in::to_getfemint_mdbrick(bool writeable)
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != MDBRICK_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a md-brick descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o,writeable);
    return object_to_mdbrick(o);
  }

  getfemint_mdstate *
  mexarg_in::to_getfemint_mdstate(bool writeable)
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != MDSTATE_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a md-state descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o,writeable);
    return object_to_mdstate(o);
  }

  getfemint_mesh_slice *
  mexarg_in::to_getfemint_mesh_slice(bool writeable)
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != SLICE_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a mesh slice descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o,writeable);
    return object_to_mesh_slice(o);
  }

  /*getfem::stored_mesh_slice *
  mexarg_in::to_mesh_slice() {
    return &to_getfemint_mesh_slice()->mesh_slice();
  }
  */

  getfemint_levelset *
  mexarg_in::to_getfemint_levelset(bool writeable) {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != LEVELSET_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a levelset descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o,writeable);
    return object_to_levelset(o);
  }
    
  getfemint_mesh_levelset *
  mexarg_in::to_getfemint_mesh_levelset(bool writeable) {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != MESH_LEVELSET_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a mesh_levelset descriptor, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    error_if_nonwritable(o,writeable);
    return object_to_mesh_levelset(o);
  }

  /*  getfem::mesh_level_set *
  mexarg_in::to_mesh_levelset() {
    return &to_getfemint_mesh_levelset(true)->mesh_levelset();
  }
  */

  getfemint_precond *
  mexarg_in::to_precond()
  {
    id_type id, cid;
    to_object_id(&id,&cid);
    if (cid != PRECOND_CLASS_ID) {
      THROW_BADARG("argument " << argnum << " should be a preconditioner, its class is " << name_of_getfemint_class_id(cid));
    }
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    return object_to_precond(o);
  }

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
  
  getfem::pintegration_method
  mexarg_in::to_integration_method()
  {
    id_type id,cid;
    to_object_id(&id,&cid);
    if (cid != INTEG_CLASS_ID)
      THROW_BADARG("Argument " << argnum << 
		   " should be an integration method descriptor");
    if (!exists_integ(id)) {
      THROW_BADARG("Argument " << argnum << 
		   " is not a valid integration method handle");
    }
    return addr_integ(id);
  }

  getfem::pmat_elem_type
  mexarg_in::to_mat_elem_type()
  {
    id_type id,cid;
    to_object_id(&id,&cid);
    if (cid != ELTM_CLASS_ID)
      THROW_BADARG("Argument " << argnum << 
		   " should be a elementary matrix descriptor.");
    if (!exists_matelemtype(id))
      THROW_BADARG("Argument " << argnum << 
		   " is not a valid elementary matrix handle");
    return addr_matelemtype(id);
  }

  getfemint_pfem*
  mexarg_in::to_getfemint_pfem()
  {
    id_type id,cid;
    to_object_id(&id,&cid);
    if (cid != FEM_CLASS_ID)
      THROW_BADARG("Argument " << argnum << 
		   " should be a fem descriptor");
    getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
    return object_to_pfem(o);
  }

  getfem::pfem
  mexarg_in::to_fem()
  {
    getfemint_pfem *p = to_getfemint_pfem();
    return p->pfem();
  }


  bgeot::pgeometric_trans
  mexarg_in::to_pgt() {
    id_type id,cid;
    to_object_id(&id,&cid);
    if (cid != GEOTRANS_CLASS_ID)
      THROW_BADARG("Argument " << argnum << 
		   " is not a geometric transformation handle");
    if (!getfemint::exists_pgt(id))
      THROW_BADARG("Argument " << argnum << 
		   " refers to a geometric transformation that does not exists");
    return getfemint::addr_pgt(id);
  }

  bgeot::pconvex_structure
  mexarg_in::to_convex_structure() {
    id_type id,cid;
    to_object_id(&id,&cid);
    if (cid != CVSTRUCT_CLASS_ID)
      THROW_BADARG("Argument " << argnum << 
		   " is not a convex structure handle");
    if (!getfemint::exists_convex_structure(id))
      THROW_BADARG("Argument " << argnum << 
		   " refers to a convex structure that does not exists");
    return getfemint::addr_convex_structure(id);
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
  void
  mexarg_in::to_sparse(gf_real_sparse_csc_const_ref& M) {
    if (gfi_array_get_class(arg) != GFI_SPARSE) {
      THROW_BADARG("Argument " << argnum << " was expected to be a sparse matrix");
    }
    if (is_complex()) {
      THROW_BADARG("Argument " << argnum << " cannot be a complex sparse matrix");
    }
    assert(gfi_array_get_ndim(arg)==2);
    M = gf_real_sparse_csc_const_ref(gfi_sparse_get_pr(arg), gfi_sparse_get_ir(arg), gfi_sparse_get_jc(arg),
				     gfi_array_get_dim(arg)[0],gfi_array_get_dim(arg)[1]);
  }

  void
  mexarg_in::to_sparse(gf_cplx_sparse_csc_const_ref& M) {
    if (gfi_array_get_class(arg) != GFI_SPARSE) {
      THROW_BADARG("Argument " << argnum << " was expected to be a sparse matrix");
    }
    if (!is_complex()) {
      THROW_BADARG("Argument " << argnum << " cannot be a real sparse matrix");
    }
    assert(gfi_array_get_ndim(arg)==2);
    M = gf_cplx_sparse_csc_const_ref((complex_type*)gfi_sparse_get_pr(arg), gfi_sparse_get_ir(arg), gfi_sparse_get_jc(arg),
				     gfi_array_get_dim(arg)[0],gfi_array_get_dim(arg)[1]);
  }

  /* get a (native or getfem) sparse matrix */ 
  dal::shared_ptr<gsparse> 
  mexarg_in::to_sparse() {
    if (gfi_array_get_class(arg) == GFI_SPARSE) {
      dal::shared_ptr<gsparse> pgsp(new gsparse(arg));
      return pgsp;
    } else {
      id_type id,cid;
      to_object_id(&id,&cid);
      if (cid != GSPARSE_CLASS_ID) 
	THROW_BADARG("Argument " << argnum << " was expected to be a sparse matrix");
      getfem_object *o = workspace().object(id,name_of_getfemint_class_id(cid));
      //cout << "id = " << id << ", cid= "<< cid << " " << name_of_getfemint_class_id(cid) << ", ->o = " << o << "\n";
      return object_to_gsparse(o)->ref();
    }
  }
  
  getfemint_gsparse *
  mexarg_in::to_getfemint_gsparse() {
    if (gfi_array_get_class(arg) == GFI_SPARSE) {
      THROW_BADARG("Argument " << argnum << " was expected as a GETFEM sparse matrix, not a native sparse matrix");
    } else {
      id_type id,cid;
      to_object_id(&id,&cid);
      if (cid != GSPARSE_CLASS_ID) 
	THROW_BADARG("Argument " << argnum << " was expected to be a sparse matrix");
      return object_to_gsparse(workspace().object(id,name_of_getfemint_class_id(cid)));
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
  mexarg_out::from_integer(int i)
  {
    if (config::can_return_integer()) {
      arg = checked_gfi_array_create_0(GFI_INT32);
      *gfi_int32_get_data(arg) = i;
    } else from_scalar(i);
  }

  void
  mexarg_out::from_scalar(double v)
  {
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
  void
  mexarg_out::from_sparse(gf_real_sparse_by_col& M, output_sparse_fmt fmt) {
    gsparse gsp; 
    from_sparse(gsp.destructive_assign(M), fmt);
  }

  void
  mexarg_out::from_sparse(gf_cplx_sparse_by_col& M, output_sparse_fmt fmt) {
    gsparse gsp; 
    from_sparse(gsp.destructive_assign(M), fmt);
  }

  void
  mexarg_out::from_sparse(gsparse& M, output_sparse_fmt fmt) {
    if (fmt == USE_DEFAULT_SPARSE) {
      fmt = (config::prefer_native_sparse() ? USE_NATIVE_SPARSE : USE_GSPARSE);
    }
    if (fmt == USE_GSPARSE) {
      gsparse &gsp = create_gsparse();
      gsp.swap(M);
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

  gsparse & 
  mexarg_out::create_gsparse() {
    getfemint_gsparse *ggsparse = new getfemint_gsparse();
    from_object_id(workspace().push_object(ggsparse), GSPARSE_CLASS_ID);
    return ggsparse->sparse();
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
  mexarg_out::create_carray_h(unsigned dim)
  {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE, GFI_COMPLEX);
    else arg = checked_gfi_array_create_2(1,dim, GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  carray
  mexarg_out::create_carray_v(unsigned dim)
  {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE, GFI_COMPLEX);
    else arg = checked_gfi_array_create_2(dim, 1, GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  carray
  mexarg_out::create_carray(unsigned n,unsigned m,unsigned p)
  {
    int sz[3];
    sz[0] = n; sz[1] = m; sz[2] = p;
    arg = checked_gfi_array_create(3,sz,GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  carray
  mexarg_out::create_carray(unsigned dim, unsigned nbdof)
  {
    arg = checked_gfi_array_create_2(dim, nbdof, GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
  }

  darray
  mexarg_out::create_darray_h(unsigned dim)
  {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE);
    else 
      arg = checked_gfi_array_create_2(1,dim, GFI_DOUBLE);
    return darray(arg);
  }

  darray
  mexarg_out::create_darray_v(unsigned dim)
  {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_DOUBLE);
    else 
      arg = checked_gfi_array_create_2(dim, 1, GFI_DOUBLE);
    return darray(arg);
  }

  darray
  mexarg_out::create_darray(unsigned n,unsigned m,unsigned p)
  {
    int sz[3];
    sz[0] = n; sz[1] = m; sz[2] = p;
    arg = checked_gfi_array_create(3,sz,GFI_DOUBLE);
    return darray(arg);
  }

  /* creates a 'matrix' (from the matlab point of view.. for getfem 
     it is still a vector, stored in fortran style) */
  darray
  mexarg_out::create_darray(unsigned dim, unsigned nbdof)
  {
    arg = checked_gfi_array_create_2(dim, nbdof, GFI_DOUBLE);
    return darray(arg);
  }

  iarray
  mexarg_out::create_iarray_h(unsigned dim)
  {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_INT32);
    else arg = checked_gfi_array_create_2(1,dim, GFI_INT32);
    return iarray(arg);
  }

  iarray
  mexarg_out::create_iarray_v(unsigned dim)
  {
    if (config::has_1D_arrays())
      arg = checked_gfi_array_create_1(dim, GFI_INT32);
    else arg = checked_gfi_array_create_2(dim, 1, GFI_INT32);
    return iarray(arg);
  }

  iarray
  mexarg_out::create_iarray(unsigned n,unsigned m,unsigned p)
  {
    int sz[3];
    sz[0] = n; sz[1] = m; sz[2] = p;
    arg = checked_gfi_array_create(3,sz,GFI_INT32);
    return iarray(arg);
  }

  /* creates a 'matrix' (from the matlab point of view.. for getfem 
     it is still a vector, stored in fortran style) */
  iarray
  mexarg_out::create_iarray(unsigned dim, unsigned nbdof)
  {
    arg = checked_gfi_array_create_2(dim, nbdof, GFI_INT32);
    return iarray(arg);
  }

  darray mexarg_out::create_array(const array_dimensions &dims, double) 
  {
    arg = checked_gfi_array_create(dims.ndim(),(const int*)dims.sizes(),GFI_DOUBLE);
    return darray(arg);
  }

  carray mexarg_out::create_array(const array_dimensions &dims, complex_type) 
  {
    arg = checked_gfi_array_create(dims.ndim(),(const int*)dims.sizes(),GFI_DOUBLE, GFI_COMPLEX);
    return carray(arg);
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
    if (cmd_strmatch(cmdname,s)) {
      if (min_argout > 0 && out.narg_known() && 
          out.narg_in_range(0, min_argout-1)) {
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
      assert(gfi_array_get_class(p[0])==GFI_CELL);
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
  }
  
  mexargs_out::~mexargs_out() {
    if (!okay) {
      for (size_type i=0; i < out.size(); ++i) {
	if (out[i]) { gfi_array_destroy(out[i]); free(out[i]); }	
      }
      out.clear();
      workspace().destroy_newly_created_objects();
    } else workspace().commit_newly_created_objects();
  }
  
  /* ensure that out[idx] is valid */
  void mexargs_out::check() const { 
    if (nb_arg != -1) {
      if ((idx >= nb_arg) && !(idx==0)) THROW_INTERNAL_ERROR;
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

  void check_cv_fem(const getfem::mesh_fem& mf, size_type cv) {
    if (!mf.convex_index()[cv]) THROW_ERROR( "convex " << cv+config::base_index() << " has no FEM");
  }

  void check_cv_im(const getfem::mesh_im& mim, size_type cv) {
    if (!mim.convex_index()[cv]) THROW_ERROR( "convex " << cv+config::base_index() << " has no integration method!");
  }

  static double NaN = 0.;
  const double& get_NaN() {
    if (NaN == 0.) {
#ifdef NAN
    NaN = NAN;
#else
      NaN=1e308;

      /* disabled for now ... the dec alpha throws SIGFPE when i run the make check .. */

      /*double a=0.,b=0.;
      struct sigaction osa,nsa;
      sigaction(SIGFPE, 0, &nsa);
      nsa.sa_handler = SIG_IGN;
      sigaction(SIGFPE, &nsa, &osa);
      NaN = a / b;
      sigaction(SIGFPE,&osa, 0);      
      //assert(NaN != NaN);//great ... assertion fails on dec/alpha and origin2k ..
      */
#endif
    }
    return NaN;
  }

  bool is_NaN(const double& v) {
    double w = v;
    void *p = &w;
    if (memcmp(p, (void*)const_cast<double*>(&get_NaN()), sizeof(v)) == 0) return true;
    return ((v!=w) || !(v==w)); /* add a small prayer to avoid any optimization by the compiler.. */
  }

} /* namespace getfemint */
