/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2001-2020 Yves Renard, Julien Pommier.

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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

#ifndef GETFEM_MISC_H__
#define GETFEM_MISC_H__

#include <getfemint.h>
#include <getfem/getfem_mesh.h>
#include <getfem/bgeot_convex.h>
#include <gmm/gmm_iter.h>

namespace getfem {
  class abstract_hyperelastic_law;
  typedef std::shared_ptr<const abstract_hyperelastic_law> phyperelastic_law;
}

namespace getfem {
  class abstract_constraints_projection;
  typedef std::shared_ptr<const abstract_constraints_projection> pconstraints_projection;
}

namespace getfemint {
  typedef getfem::convex_face convex_face;
  typedef getfem::convex_face_ct convex_face_ct;

  gfi_array* convert_to_gfi_sparse(const gf_real_sparse_by_row & smat, double threshold=1e-12);
  gfi_array* convert_to_gfi_sparse(const gf_real_sparse_by_col & smat, double threshold=1e-12);


  void build_edge_list(const getfem::mesh &m, bgeot::edge_list &el, mexargs_in &in);

  void build_convex_face_lst(const getfem::mesh& m, std::vector<convex_face>& l, const iarray *v);

  getfem::mesh_region to_mesh_region(const iarray &v);
  getfem::mesh_region to_mesh_region(const getfem::mesh& m, const iarray *v=0);
  inline getfem::mesh_region to_mesh_region(const getfem::mesh& m, mexargs_in &in) {
    if (in.remaining()) { iarray v = in.pop().to_iarray(); return to_mesh_region(m, &v); }
    else return to_mesh_region(m);
  }

  void interpolate_on_convex_ref(const getfem::mesh_fem *mf, getfem::size_type cv, 
				 const std::vector<getfem::base_node> &pt, 
				 const darray& U,
				 getfem::base_matrix &pt_val);
  void
  eval_on_triangulated_surface(const getfem::mesh* mesh, int Nrefine,
			       const std::vector<convex_face>& cvf,
			       mexargs_out& out,
			       const getfem::mesh_fem *pmf, const darray& U);
  
  const getfem::phyperelastic_law &abstract_hyperelastic_law_from_name(const std::string &lawname, size_type N);


  const getfem::pconstraints_projection &abstract_constraints_projection_from_name(const std::string &projname);


  class interruptible_iteration : public gmm::iteration {
  public:
    interruptible_iteration(double r=1e-8);
  };
}
#endif
