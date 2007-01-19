// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2007-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
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

#include <getfem_trace_mesh_fem_level_set.h>

namespace getfem {

  //======================================================================
  // Fem part
  //======================================================================

  void sub_space_fem::init() {
    cvr = org_fem->ref_convex(cv);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = 1;

    std::stringstream nm;
    nm << "FEM_SUB_SPACE(" << org_fem->debug_name() << ", " << B << ")";
    debug_name_ = nm.str();
    
    init_cvs_node();
    if (org_fem->target_dim() != 1)
      DAL_THROW(to_be_done_error, "Vectorial fems not (yet) supported");

    base_node P(dim()); gmm::fill(P, 1./4.);
    for (size_type k = 0; k < ind.size(); ++k)
	add_node(global_dof(dim()), P);
  }

  size_type sub_space_fem::index_of_global_dof(size_type, size_type j) const
  { return ind[j]; }

  void sub_space_fem::base_value(const base_node &, 
				 base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void sub_space_fem::grad_base_value(const base_node &, 
				      base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void sub_space_fem::hess_base_value(const base_node &, 
			     base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }

  void sub_space_fem::real_base_value(const fem_interpolation_context &c,
				    base_tensor &t, bool) const {
    fem_interpolation_context c0 = c;
    base_tensor val_e;
    
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.base_value(val_e);
    
    t.mat_transp_reduction(val_e, B, 0);
  }

  void sub_space_fem::real_grad_base_value(const fem_interpolation_context &c,
					   base_tensor &t, bool) const {
    fem_interpolation_context c0 = c;
    base_tensor grad_e;
    
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.grad_base_value(grad_e);

    t.mat_transp_reduction(grad_e, B, 0);
  }
  
  void sub_space_fem::real_hess_base_value(const fem_interpolation_context &c,
				  base_tensor &t, bool) const {
    fem_interpolation_context c0 = c;
    base_tensor hess_e;
   
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.hess_base_value(hess_e);
    
    t.mat_transp_reduction(hess_e, B, 0);
  }

  //======================================================================
  // Mesh fem part
  //======================================================================

  
  void trace_mesh_fem_level_set::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }
  void trace_mesh_fem_level_set::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }
  void trace_mesh_fem_level_set::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void trace_mesh_fem_level_set::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    is_adapted = false;
  }
  
  trace_mesh_fem_level_set::trace_mesh_fem_level_set(const mesh_level_set &me,
						     const mesh_fem &mef)
    : mesh_fem(mef.linked_mesh()), mls(me), mf(mef) {
    if (mf.get_qdim() != 1)
      DAL_THROW(to_be_done_error, "base mesh_fem for mesh_fem_level_set has "
		"to be of qdim one for the moment ...");
    this->add_dependency(mls);
    is_adapted = false;
  }

  DAL_SIMPLE_KEY(special_tracemf_key, pfem);
  void trace_mesh_fem_level_set::adapt(void) {
    context_check();
    clear();
    
//     dal::bit_vector kept_dof;
      
//     for (dal::bv_visitor cv(linked_mesh().convex_index());
// 	 !cv.finished(); ++cv)
//       if (!rejected_elt.is_in(cv)) {
// 	dal::bit_vector selected_dofs;
// 	for (size_type i = 0; i < mf.nb_dof_of_element(cv); ++i)
// 	  if (kept_dofs.is_in(mf.ind_dof_of_element(cv)[i])) 
// 	    selected_dofs.add(i);
// 	if (selected_dofs.card() == mf.nb_dof_of_element(cv))
// 	  set_finite_element(cv, mf.fem_of_element(cv));
// 	else if (selected_dofs.card()) {
// 	  assert(mf.fem_of_element(cv) != 0);
// 	  pfem pf = new sub_space_fem(mf.fem_of_element(cv), std::vector<unsigned>(), base_matrix(), cv);
// 	  dal::add_stored_object(new special_tracemf_key(pf), pf,
// 				 pf->ref_convex(0),
// 				 pf->node_tab(0));
// 	  build_methods.push_back(pf);
// 	  set_finite_element(cv, pf);
// 	}
//       }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                            */

