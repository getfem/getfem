// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_interpolated_fem.h : definition of a finite element
//           method which interpolates a fem on a different mesh.
// Date    : October 29, 2004.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================
// To be corrected : dependencies. The mesh fem using this fem will not
//                   depend on the mesh fem arguments.



#ifndef GETFEM_INTERPOLATED_FEM_H__
#define GETFEM_INTERPOLATED_FEM_H__

#include <getfem_fem.h>
#include <getfem_mesh_fem.h>
#include <getfem_mesh_im.h>
#include <bgeot_rtree.h>
#include <bgeot_geotrans_inv.h>

namespace getfem {

  // Object representing global transformation. To be derived.

  struct virtual_interpolated_func {
    virtual void val(const base_node&, base_node &) const
    { DAL_THROW(dal::failure_error,"this interpolated_func has no value"); }
    virtual void grad(const base_node&, base_matrix &) const
    { DAL_THROW(dal::failure_error,"this interpolated_func has no gradient"); }
    virtual void hess(const base_node&, base_matrix &) const
    { DAL_THROW(dal::failure_error,"this interpolated_func has no hessian"); }
    virtual ~virtual_interpolated_func() {}
  };


  typedef const virtual_interpolated_func *pinterpolated_func;


  struct gausspt_interpolation_data {
    size_type elt;        // element of mf1 under this gauss point
    size_type flags;      // flags & 1 : there is an element or not
    // flags & 2 : base_val is stored
    // flags & 4 : grad_val is stored
    base_node ptref;      // coords on reference element of mf1 element
    base_tensor base_val; // optional storage of the base values
    base_tensor grad_val; // optional storage of the grad base values
    std::vector<size_type> local_dof; // correspondance between dof of the
    // mf1 element and dof of the interpolated element.
  };

  class interpolated_fem : public virtual_fem, public context_dependencies {
    
  protected :



    struct elt_interpolation_data {
      size_type nb_dof;
      std::vector<gausspt_interpolation_data> gausspt;
      std::vector<size_type> inddof;
    }; 

    const mesh_fem &mf;    // mf represents the original finite element method
                            // to be interpolated.
    const mesh_im &mim;    // mesh on which mf1 is interpolated. contains
                            // also the integration method.
    pinterpolated_func pif; // optional transformation

    bool store_values;
    dal::bit_vector blocked_dof;

    // auxiliary variables
    mutable std::vector<elt_interpolation_data> elements;
    mutable bgeot::rtree boxtree; // Tree containing the bounding box
                                  // of mf1 elements
    mutable std::vector<size_type> ind_dof;
    mutable size_type cv_stored;
    mutable bgeot::rtree::pbox_set boxlst;
    mutable bgeot::geotrans_inv_convex gic;
    mutable base_tensor taux;
    mutable fem_interpolation_context fictx;
    mutable size_type fictx_cv;
    mutable base_matrix G;
    mutable std::vector<base_node> node_tab_;
    mutable bgeot::multi_index mi2, mi3;
    mutable base_node ptref;
    mutable gmm::dense_matrix<scalar_type> trans;

    void build_rtree(void) const;

    bool find_a_point(base_node pt, base_node &ptr,
		      size_type &cv) const;

    void update_from_context(void) const;
    inline void actualize_fictx(pfem pf, size_type cv,
				const base_node &ptr) const;

  public :

    virtual size_type nb_dof(size_type cv) const;
    virtual size_type index_of_global_dof(size_type cv, size_type i) const;
    virtual bgeot::pconvex_ref ref_convex(size_type cv) const;
    virtual const bgeot::convex<base_node> &node_convex(size_type cv) const;
    virtual bgeot::pstored_point_tab node_tab(size_type) const;
    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;
    void real_base_value(const fem_interpolation_context& c, 
			 base_tensor &t) const;
    void real_grad_base_value(const fem_interpolation_context& c, 
			      base_tensor &t) const;
    void real_hess_base_value(const fem_interpolation_context&, 
			      base_tensor &) const;


    interpolated_fem(const mesh_fem &mef, const mesh_im &meim,
		     pinterpolated_func pif_ = 0,
		     dal::bit_vector blocked_dof = dal::bit_vector(),
		     bool store_val = true);
  };
  
  
}  /* end of namespace getfem.                                            */

#endif
  
