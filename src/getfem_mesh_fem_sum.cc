// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_level_set.cc : definition of a finite element
//           method reprensenting a direct sum of two mesh_fem.
// Date    : March 18, 2005..
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


#include <getfem_fem_fem_sum.h>

namespace getfem {
    
  void fem_sum::init() {
    cvr = bfem->ref_convex(0);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = bfem->target_dim();
 
    ls_index.sup(0, mls.nb_sums());

    for (size_type i=0; i < mls.nb_sums(); ++i) {
      char c = '*';
      for (size_type k=0; k < bfem->nb_dof(0); ++k) {
	if (dofzones[k]) {
	  for (size_type j=0; j < dofzones[k]->size(); ++j) {
	    char d = (*(*dofzones[k])[j])[i];
	    if (c == '*') c = d;
	    else if (c != d) { ls_index.add(i); break; }
	  }
	}
      }
    }
    
    init_cvs_node();
    cout << "creating fem_level_sel:\n";
    for (size_type k = 0; k < bfem->nb_dof(0); ++k) {
      if (!dofzones[k]) {
	add_node(bfem->dof_types()[k], bfem->node_of_dof(0,k));
      } else {
	for (size_type j = 0; j < dofzones[k]->size(); ++j) {
	  cout << " -> +dof: '" << *(*dofzones[k])[j] << "'\n";
	  add_node(xfem_dof(bfem->dof_types()[k], j+1000), /* +1000 to avoid messing with the real xfem */
		   bfem->node_of_dof(0,k));
	}
      }
    }
    cout << "done (nb_dof=" << nb_dof(0) << ")\n";
  }

  void fem_sum::base_value(const base_node &x, 
				 base_tensor &t) const
  { bfem->base_value(x, t); }
  void fem_sum::grad_base_value(const base_node &x, 
				      base_tensor &t) const
  { bfem->grad_base_value(x, t); }
  void fem_sum::hess_base_value(const base_node &x, 
			     base_tensor &t) const
  { bfem->hess_base_value(x, t); }

  void fem_sum::real_base_value(const fem_interpolation_context &c,
				      base_tensor &t) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin();
    base_tensor tt; c.base_value(tt);
    base_tensor::const_iterator itf = tt.begin();
    std::vector<scalar_type> lsval(mls.nb_sums());
    for (dal::bv_visitor i(ls_index); !i.finished(); ++i) {
      mesher_sum eval = mls.get_sum(i)->
	mls_of_convex(c.convex_num());
      lsval[i] = eval(c.xref());
    }
    for (dim_type q = 0; q < target_dim(); ++q) {
      for (size_type d = 0; d < bfem->nb_base(0); ++d, ++itf) {
	if (dofzones[d]) { /* enriched dof ? */
	  const dof_ls_enrichment &de = *dofzones[d];
	  for (size_type k = 0; k < de.size(); ++k) {
	    scalar_type v = *itf;
	    for (dal::bv_visitor il(ls_index); !il.finished(); ++il) {
	      if ((*de[k])[il] == '0') continue;
	      if (((*de[k])[il] == '+' && lsval[il] < 0) ||
		  ((*de[k])[il] == '-' && lsval[il] > 0))
		v = 0;
	    }
	    *it++ = v;
	  }
	} else *it++ = *itf;
      }
    }
    assert(it == t.end());
  }

  void fem_sum::real_grad_base_value(const fem_interpolation_context &c,
					   base_tensor &t) const {
    bgeot::multi_index mi(3);
    mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);

    base_tensor tt; c.grad_base_value(tt);

    base_tensor::iterator it = t.begin();
    base_tensor::const_iterator itf = tt.begin();

    std::vector<scalar_type> lsval(mls.nb_sums());
    for (dal::bv_visitor i(ls_index); !i.finished(); ++i) {
      mesher_sum eval = mls.get_sum(i)->
	mls_of_convex(c.convex_num());
      lsval[i] = eval(c.xref());
      cout << "cv=" << c.convex_num() << ", xref=" << c.xref() << " (xreal=" << c.xreal() << " : lsval[" << i << "]=" << lsval[i] << "\n";
    }
    for (dim_type i = 0; i < c.N() ; ++i) {
      for (dim_type q = 0; q < target_dim(); ++q) {
	for (size_type d = 0; d < bfem->nb_base(0); ++d, ++itf) {
	  if (dofzones[d]) { /* enriched dof ? */
	    const dof_ls_enrichment &de = *dofzones[d];
	    for (size_type k = 0; k < de.size(); ++k) {
	      scalar_type v = *itf;
	      for (dal::bv_visitor il(ls_index); !il.finished(); ++il) {
		if ((*de[k])[il] == '0') continue;
		if (((*de[k])[il] == '+' && lsval[il] < 0) ||
		    ((*de[k])[il] == '-' && lsval[il] > 0))
		  v = 0;
	      }
	      *it++ = v;
	    }
	  } else *it++ = *itf;
	}
      }
    }
    assert(it == t.end());
  }
  
  void fem_sum::real_hess_base_value(const fem_interpolation_context &,
				  base_tensor &) const {
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
    is_adapted = false;
  }
  
  DAL_SIMPLE_KEY(special_mflsum_key, pfem);

  void mesh_fem_sum::build_method_of_convex(size_type cv) {
    ... a modifier;
      pfem pf = new fem_sum(index_ref_iterator
				  (dof_enrichments.begin(),
				   mf.ind_dof_of_element(cv).begin()) ,
				  mf.fem_of_element(cv), mls);
      dal::add_stored_object(new special_mflsum_key(pf), pf,
			     pf->ref_convex(0),
			     pf->node_tab(0));
      build_methods.push_back(pf);
      set_finite_element(cv, pf);
  }
  
  void mesh_fem_sum::adapt(void) {
    clear();

    for (dal::bv_visitor i(linked_mesh().convex_index()); !i.finished(); ++i) {
      if (enriched_elements[i]) build_method_of_convex(i); else
	set_finite_element(i, mf.fem_of_element(i));
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                             */

