// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_level_set.cc : Dealing with level set representation.
//           
// Date    : January 31, 2005.
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


#include <getfem_fem_level_set.h>

namespace getfem {
    
  void fem_level_set::init() {
    cvr = bfem->ref_convex(0);
    dim_ = cvr->structure()->dim();
    is_equiv = false; real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = bfem->target_dim();
 
    ls_index.sup(0, mls.nb_level_sets());

    for (size_type i=0; i < mls.nb_level_sets(); ++i) {
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
    cout << "creating fem_level_set:\n";
    for (size_type k = 0; k < bfem->nb_dof(0); ++k) {
      if (!dofzones[k]) {
	add_node(bfem->dof_types()[k], bfem->node_of_dof(0,k));
      } else {
	for (size_type j = 0; j < dofzones[k]->size(); ++j) {
	  cout << " -> +dof: '" << *(*dofzones[k])[j] << "'\n";
	  add_node(xfem_dof(bfem->dof_types()[k], j+XFEM_INDEX_START), /* +1000 to avoid messing with the real xfem */
		   bfem->node_of_dof(0,k));
	}
      }
    }
    cout << "done (nb_dof=" << nb_dof(0) << ")\n";
  }

  void fem_level_set::base_value(const base_node &, 
				 base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_level_set::grad_base_value(const base_node &, 
				      base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_level_set::hess_base_value(const base_node &, 
			     base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element.");  }

  void fem_level_set::real_base_value(const fem_interpolation_context &c,
				      base_tensor &t) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    base_tensor::iterator it = t.begin();
    fem_interpolation_context c0 = c;
    if (c0.have_pfp())
      c0.set_pfp(fem_precomp(bfem, &c0.pfp()->get_point_tab()));
    else  c0.set_pf(bfem); 
    base_tensor tt; c0.base_value(tt);
    base_tensor::const_iterator itf = tt.begin();
    std::vector<scalar_type> lsval(mls.nb_level_sets());
    for (dal::bv_visitor i(ls_index); !i.finished(); ++i) {
      mesher_level_set eval = mls.get_level_set(i)->
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

  void fem_level_set::real_grad_base_value(const fem_interpolation_context &c,
					   base_tensor &t) const {
    bgeot::multi_index mi(3);
    mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    fem_interpolation_context c0 = c;
    if (c0.have_pfp())
      c0.set_pfp(fem_precomp(bfem, &c0.pfp()->get_point_tab()));
    else  c0.set_pf(bfem); 
    base_tensor tt; c0.grad_base_value(tt);

    base_tensor::iterator it = t.begin();
    base_tensor::const_iterator itf = tt.begin();

    std::vector<scalar_type> lsval(mls.nb_level_sets());
    for (dal::bv_visitor i(ls_index); !i.finished(); ++i) {
      mesher_level_set eval = mls.get_level_set(i)->
	mls_of_convex(c.convex_num());
      lsval[i] = eval(c.xref());
      //cout << "cv=" << c.convex_num() << ", xref=" << c.xref() << " (xreal=" << c.xreal() << " : lsval[" << i << "]=" << lsval[i] << "\n";
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
  
  void fem_level_set::real_hess_base_value(const fem_interpolation_context &,
				  base_tensor &) const {
    DAL_THROW(to_be_done_error,
	      "Sorry order 2 derivatives for fem_level_set to be done.");
  }


}  /* end of namespace getfem.                                             */

