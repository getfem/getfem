// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_fem_product.cc : definition of a finite element
//           method reprensenting a "product" of two mesh_fems.
// Date    : April 8, 2005.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 1999-2005 Yves Renard
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


#include <getfem_mesh_fem_product.h>

namespace getfem {
    
  void fem_product::init() {
    
    if (pfems[0]->target_dim() > 1) 
      DAL_THROW(to_be_done_error, "To be done");
    if (pfems[1]->target_dim() > 1) 
      DAL_THROW(failure_error, "The second finite element should be scalar");

    cvr = pfems[0]->ref_convex(cv);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = 1;
    
    init_cvs_node();
    for (dal::bv_visitor i(enriched_dof1); !i.finished(); ++i) {
      for (size_type j = 0; j < pfems[1]->nb_dof(cv); ++j)
	add_node(xfem_dof(pfems[0]->dof_types()[i], xfem_index+j),
		 pfems[0]->node_of_dof(cv,i));
    }
  }

  void fem_product::base_value(const base_node &, 
				 base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_product::grad_base_value(const base_node &, 
				      base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_product::hess_base_value(const base_node &, 
			     base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }

  void fem_product::real_base_value(const fem_interpolation_context &c,
				      base_tensor &t) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;

    fem_interpolation_context c0 = c;
    std::vector<base_tensor> val_e(2);
    for (size_type k = 0; k < 2; ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab()));
      } else { c0.set_pf(pfems[k]); }
      c0.base_value(val_e[k]);
    }
    
    assert(target_dim() == 1);
    for (dal::bv_visitor i(enriched_dof1); !i.finished(); ++i) {
      itf = val_e[1].begin();
      scalar_type e = val_e[0][i];
      for (size_type j = 0; j < pfems[1]->nb_dof(cv); ++j)
	*it++ = *itf++ * e;
    }
    assert(it == t.end());
  }

  void fem_product::real_grad_base_value(const fem_interpolation_context &c,
					   base_tensor &t) const {
    bgeot::multi_index mi(3);
    mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin();
    
    fem_interpolation_context c0 = c;
    std::vector<base_tensor> grad_e(2), val_e(2);
    for (size_type k = 0; k < 2; ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab()));
      } else { c0.set_pf(pfems[k]); }
      c0.grad_base_value(grad_e[k]);
      c0.base_value(val_e[k]);
    }

    assert(target_dim() == 1);
    for (dim_type k = 0; k < c.N() ; ++k) {
      for (dal::bv_visitor i(enriched_dof1); !i.finished(); ++i) {
	size_type posg0 = k * pfems[0]->nb_base(cv);
	size_type posg1 = k * pfems[1]->nb_base(cv);
	for (size_type j = 0; j < pfems[1]->nb_base(cv); ++j)
	  *it++ = grad_e[0][i + posg0] * val_e[1][j]
	    + grad_e[1][j + posg1] * val_e[0][i];
      }
    }
    assert(it == t.end());
  }
  
  void fem_product::real_hess_base_value(const fem_interpolation_context &,
				  base_tensor &) const {
    DAL_THROW(to_be_done_error,
	      "Sorry order 2 derivatives for fem_product to be done.");
  }

  void mesh_fem_product::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }
  void mesh_fem_product::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }
  void mesh_fem_product::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void mesh_fem_product::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    is_adapted = false;
  }
  
  DAL_SIMPLE_KEY(special_mflproduct_key, pfem);
  
  void mesh_fem_product::adapt(void) {
    clear();
    for (dal::bv_visitor cv(linked_mesh().convex_index()); !cv.finished();
	 ++cv) {
      dal::bit_vector local_enriched_dof;
      for (size_type i = 0; i < mf1.nb_dof_of_element(cv); ++i)
	if (enriched_dof.is_in(mf1.ind_dof_of_element(cv)[i]))
	  local_enriched_dof.add(i);
      if (local_enriched_dof.card() > 0) {
	pfem pf = new fem_product(mf1.fem_of_element(cv),
				  mf2.fem_of_element(cv), cv,
				  xfem_index, local_enriched_dof);
	dal::add_stored_object(new special_mflproduct_key(pf), pf,
			       pf->ref_convex(0),
			       pf->node_tab(0));
	build_methods.push_back(pf);
	set_finite_element(cv, pf);
      }
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                             */

