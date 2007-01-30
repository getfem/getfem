// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2007 Yves Renard
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


#include "getfem/getfem_partial_fem.h"

namespace getfem {
    
  void partial_fem::init() {
    cvr = org_fem->ref_convex(cv);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = 1;

    std::stringstream nm;
    nm << "FEM_PARTIAL(" << org_fem->debug_name() << ", "
       << selected_dofs << ")";
    debug_name_ = nm.str();
    
    init_cvs_node();
    if (org_fem->target_dim() != 1)
      DAL_THROW(to_be_done_error, "Vectorial fems not supported");

    for (size_type k = 0; k < org_fem->nb_dof(cv); ++k) {
      if (selected_dofs.is_in(k)) {
	ind.push_back(k);
	add_node(org_fem->dof_types()[k], org_fem->node_of_dof(cv,k));
      }
    }
  }

  size_type partial_fem::index_of_global_dof(size_type, size_type j) const {
    return org_fem->index_of_global_dof(cv, ind[j]);
  }

  void partial_fem::base_value(const base_node &, 
				 base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void partial_fem::grad_base_value(const base_node &, 
				      base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void partial_fem::hess_base_value(const base_node &, 
			     base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }

  void partial_fem::real_base_value(const fem_interpolation_context &c,
				    base_tensor &t, bool) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;

    fem_interpolation_context c0 = c;
    base_tensor val_e;
    
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.base_value(val_e);
    
    for (dim_type q = 0; q < target_dim(); ++q) {
      itf = val_e.begin() + q * org_fem->nb_base(cv);
      for (size_type i = 0; i <  org_fem->nb_base(cv); ++i, ++itf)
	if (selected_dofs.is_in(i)) *it++ = *itf;
    }
    assert(it == t.end());
  }

  void partial_fem::real_grad_base_value(const fem_interpolation_context &c,
					   base_tensor &t, bool) const {
    bgeot::multi_index mi(3);
    mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;
    
    fem_interpolation_context c0 = c;
    base_tensor grad_e;
    
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.grad_base_value(grad_e);
    
    for (dim_type k = 0; k < c.N() ; ++k) {
      for (dim_type q = 0; q < target_dim(); ++q) {
	  itf = grad_e.begin()
	    + (k * target_dim() + q) * org_fem->nb_base(cv); 
	  for (size_type i = 0; i < org_fem->nb_base(cv); ++i, ++itf)
	    if (selected_dofs.is_in(i))  *it++ = *itf;
      }
    }
    assert(it == t.end());
  }
  
  void partial_fem::real_hess_base_value(const fem_interpolation_context &c,
				  base_tensor &t, bool) const {
    bgeot::multi_index mi(4);
    mi[3] = mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;
    
    fem_interpolation_context c0 = c;
    base_tensor hess_e;
   
    if (c0.have_pfp()) {
      c0.set_pfp(fem_precomp(org_fem, &c0.pfp()->get_point_tab()));
    } else { c0.set_pf(org_fem); }
    c0.hess_base_value(hess_e);
    
    for (dim_type j = 0; j < c.N() ; ++j) {
      for (dim_type k = 0; k < c.N() ; ++k) {
	for (dim_type q = 0; q < target_dim(); ++q) {
	  itf = hess_e.begin()
	    + ((j * c.N() + k) * target_dim() + q) * org_fem->nb_base(cv); 
	  for (size_type i = 0; i < org_fem->nb_base(cv); ++i, ++itf)
	    if (selected_dofs.is_in(i)) *it++ = *itf;
	}
      }
    }
    assert(it == t.end());
  }

}  /* end of namespace getfem.                                             */

