// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
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


#include <getfem_mesh_fem_sum.h>

namespace getfem {
    
  void fem_sum::init() {
    cvr = pfems[0]->ref_convex(cv);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = 1;

    std::stringstream nm;
    nm << "FEM_SUM(" << pfems[0]->debug_name() << ",";
    for (size_type i = 1; i < pfems.size(); ++i)
      nm << pfems[1]->debug_name() << ",";
    nm << cv << ")";
    debug_name_ = nm.str();
    
    init_cvs_node();
    for (size_type i = 0; i < pfems.size(); ++i) {
      if (pfems[i]->target_dim() != 1)
	DAL_THROW(to_be_done_error, "Vectorial fems not supported");

      for (size_type k = 0; k < pfems[i]->nb_dof(cv); ++k) {
	add_node(pfems[i]->dof_types()[k], pfems[i]->node_of_dof(cv,k));
      }
    }
  }

  size_type fem_sum::index_of_global_dof(size_type , size_type j) const {
    for (size_type i = 0; i < pfems.size(); ++i) {
      size_type nb = pfems[i]->nb_dof(cv);
      if (j >= nb) j -= pfems[i]->nb_dof(cv);
      else return pfems[i]->index_of_global_dof(cv, j);
    }
    DAL_THROW(failure_error, "incoherent global dof.");
  }

  void fem_sum::base_value(const base_node &, 
				 base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_sum::grad_base_value(const base_node &, 
				      base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_sum::hess_base_value(const base_node &, 
			     base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }

  void fem_sum::real_base_value(const fem_interpolation_context &c,
				      base_tensor &t, bool) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;

    fem_interpolation_context c0 = c;
    std::vector<base_tensor> val_e(pfems.size());
    for (size_type k = 0; k < pfems.size(); ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab()));
      } else { c0.set_pf(pfems[k]); }
      c0.base_value(val_e[k]);
    }
    
    for (dim_type q = 0; q < target_dim(); ++q) {
      for (size_type k = 0; k < pfems.size(); ++k) {
	itf = val_e[k].begin() + q * pfems[k]->nb_base(cv);
	for (size_type i = 0; i <  pfems[k]->nb_base(cv); ++i)
	  *it++ = *itf++;
      }
    }
    assert(it == t.end());
  }

  void fem_sum::real_grad_base_value(const fem_interpolation_context &c,
					   base_tensor &t, bool) const {
    bgeot::multi_index mi(3);
    mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin(), itf;
    
    fem_interpolation_context c0 = c;
    std::vector<base_tensor> grad_e(pfems.size());
    for (size_type k = 0; k < pfems.size(); ++k) {
      if (c0.have_pfp()) {
	c0.set_pfp(fem_precomp(pfems[k], &c0.pfp()->get_point_tab()));
      } else { c0.set_pf(pfems[k]); }
      c0.grad_base_value(grad_e[k]);
    }

    for (dim_type k = 0; k < c.N() ; ++k) {
      for (dim_type q = 0; q < target_dim(); ++q) {
	for (size_type f = 0; f < pfems.size(); ++f) {
	  itf = grad_e[f].begin()
	    + (k * target_dim() + q) * pfems[f]->nb_base(cv); 
	  for (size_type i = 0; i < pfems[f]->nb_base(cv); ++i)
	    *it++ = *itf++;
	}
      }
    }
    assert(it == t.end());
  }
  
  void fem_sum::real_hess_base_value(const fem_interpolation_context &,
				  base_tensor &, bool) const {
    DAL_THROW(to_be_done_error,
	      "Sorry order 2 derivatives for fem_sum to be done.");
  }


  void mesh_fem_sum::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }
  void mesh_fem_sum::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }
  void mesh_fem_sum::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void mesh_fem_sum::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    situations.clear();
    is_adapted = false;
  }
  
  DAL_SIMPLE_KEY(special_mflsum_key, pfem);
  
  void mesh_fem_sum::adapt(void) {
    context_check();
    clear();
    for (dal::bv_visitor i(linked_mesh().convex_index()); !i.finished(); ++i) {
      std::vector<pfem> pfems;
      bool is_cv_dep = false;
      for (size_type j = 0; j < mfs.size(); ++j) {
	if (mfs[j]->convex_index().is_in(i)) {
	  pfem pf = mfs[j]->fem_of_element(i);
	  if (pf->nb_dof(i)) {
	    pfems.push_back(pf);
	    if (pf->is_on_real_element()) is_cv_dep = true;
	  }
	}
      }
      if (pfems.size() == 1) {
	set_finite_element(i, pfems[0]);
      }
      else if (pfems.size() > 0) {
	if (situations.find(pfems) == situations.end() || is_cv_dep) {
	  pfem pf = new fem_sum(pfems, i);
	  dal::add_stored_object(new special_mflsum_key(pf), pf,
				 pf->ref_convex(0),
				 pf->node_tab(0));
	  build_methods.push_back(pf);
	  situations[pfems] = pf;
	}
	set_finite_element(i, situations[pfems]);
      }
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                             */

