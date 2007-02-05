// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2003-2007 Yves Renard
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

#include "getfem/getfem_Xfem.h"

namespace getfem
{
  void Xfem::valid(void) {
    init_cvs_node();
    /* setup nodes of the base fem */
    if (pfb)
      for (size_type k = 0; k < pfb->nb_base(0); ++k)
	add_node(pfb->dof_types()[k], pfb->node_of_dof(0,k));
    
    /* setup nodes of the enriched fems */
    for (size_type k = 0; k < nb_func; ++k) {
      for (size_type j = 0; j < pfe(k)->nb_base(0); ++j) {
	add_node(xfem_dof(pfe(k)->dof_types()[j], func_indices[k]),
		 pfe(k)->node_of_dof(0,j));
      }
    }
    is_valid = true;
  }
  
  size_type Xfem::nb_dof(size_type) const {
    GMM_ASSERT1(is_valid, "Valid the Xfem element before using it");
    return dof_types_.size();
  }

  void Xfem::add_func(pfem pf, pXfem_func pXf, size_type ind) {
    if (!pfb) init(pf);
    nb_func ++;
    if (ind == size_type(-1)) ind = nb_func;
    funcs.resize(nb_func);
    func_indices.resize(nb_func);
    funcs[nb_func-1] = pXf;
    if (cvr != pf->ref_convex(0) || (pfb && pfb->target_dim() != pf->target_dim()))
      GMM_ASSERT1(false, "Incompatible Xfem fems");

    /* insert the new fem in the list */
    std::vector<pfem>::const_iterator it;
    if ((it=std::find(uniq_pfe.begin(), uniq_pfe.end(), pf)) == uniq_pfe.end()) {
      uniq_pfe.push_back(pf); func_pf.push_back(uniq_pfe.size()-1);
    } else {
      func_pf.push_back(it - uniq_pfe.begin());
    }

    func_indices[nb_func-1] = ind;
    is_valid = false;
  }
  
  /* create an interpolation_context array based on
     c0, for each fem of the Xfem. */
  void Xfem::get_fem_interpolation_context_tab(const fem_interpolation_context& c0,
					       std::vector<fem_interpolation_context>& vc) const {
    vc.resize(uniq_pfe.size());
    for (size_type k=0; k < uniq_pfe.size(); ++k) {
      vc[k] = c0; 
      if (c0.have_pfp()) {
	vc[k].set_pfp(fem_precomp(uniq_pfe[k], &c0.pfp()->get_point_tab()));
      } else { vc[k].set_pf(uniq_pfe[k]); }
    }
  }
  
  void Xfem::base_value(const base_node &, base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element.");  }
  void Xfem::grad_base_value(const base_node &, base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element.");  }
  void Xfem::hess_base_value(const base_node &, base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element.");  }

  void Xfem::real_base_value(const fem_interpolation_context &c,
			     base_tensor &t, bool) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    scalar_type a;
    Xfem_func_context ctx(c);
    base_tensor::iterator it = t.begin();
    fem_interpolation_context c0 = c;
    base_tensor tt; 
    if (pfb) {
      if (c0.have_pfp())
	c0.set_pfp(fem_precomp(pfb, &c0.pfp()->get_point_tab()));
      else  c0.set_pf(pfb); 
      c0.base_value(tt);
    }
    base_tensor::const_iterator itf = tt.begin();
    std::vector<fem_interpolation_context> vc; get_fem_interpolation_context_tab(c, vc);
    for (dim_type q = 0; q < target_dim(); ++q) {
      for (size_type i = 0; i < (pfb ? pfb->nb_base(0) : 0); ++i, ++itf, ++it)
          *it = *itf;
      for (size_type k = 0; k < nb_func; ++k) {
	base_tensor val_e; vc[func_pf[k]].base_value(val_e);
	ctx.pf = pfe(k);
	for (size_type i = 0; i < pfe(k)->nb_base(0); ++i, ++it) {
	  ctx.base_num = i; a = funcs[k]->val(ctx);
	  *it = val_e[i + q*pfe(k)->nb_base(0)] * a;
	}
      }
    }
  }

  void Xfem::real_grad_base_value(const fem_interpolation_context &c,
				  base_tensor &t, bool) const {
    bgeot::multi_index mi(3);
    mi[2] = c.N(); mi[1] = target_dim(); mi[0] = nb_base(0);
    t.adjust_sizes(mi);
    
    Xfem_func_context ctx(c);
    fem_interpolation_context c0 = c;
    base_tensor tt; 
    if (pfb) {
      if (c0.have_pfp())
	c0.set_pfp(fem_precomp(pfb, &c0.pfp()->get_point_tab()));
      else  c0.set_pf(pfb); 
      c0.grad_base_value(tt);
    }

    base_tensor::iterator it = t.begin();
    base_tensor::const_iterator itvf = tt.begin();
    std::vector<fem_interpolation_context> vc; get_fem_interpolation_context_tab(c, vc);
    std::vector<base_tensor> val_e(nb_func);
    std::vector<base_tensor> grad_e(nb_func);
    for (size_type i=0; i < uniq_pfe.size(); ++i) {
      vc[i].base_value(val_e[i]); vc[i].grad_base_value(grad_e[i]);
    }
    std::vector<std::vector<scalar_type> > vf(nb_func);
    std::vector<std::vector<base_small_vector> > gvf(nb_func);
    for (size_type f = 0; f < nb_func; ++f) {
      vf[f].resize(pfe(f)->nb_base(0));
      gvf[f].resize(pfe(f)->nb_base(0));
      ctx.pf = pfe(f);
      for (ctx.base_num=0; ctx.base_num < pfe(f)->nb_base(0); ++ctx.base_num) {
	vf[f][ctx.base_num] = funcs[f]->val(ctx); 
	gvf[f][ctx.base_num] = funcs[f]->grad(ctx); 
      }
    }

    //    cerr << "pfp->val(ii)={"; 
    //    for (size_type i=0; i < pfp->val(ii).size(); ++i) cerr << pfp->val(ii)[i] << " "; cerr << "}\n";
    
    for (dim_type k = 0; k < c.N() ; ++k) {
      for (dim_type q = 0; q < target_dim(); ++q) {
	for (size_type i = 0; i < (pfb ? pfb->nb_base(0) : 0); ++i, ++it)
	    *it = *itvf++;
	for (size_type f = 0; f < nb_func; ++f) {
          size_type posg = pfe(f)->nb_base(0)*(q + k*target_dim());
          size_type posv = pfe(f)->nb_base(0)*q;
	  for (size_type i = 0; i < pfe(f)->nb_base(0); ++i, ++it) {
	    *it = grad_e[func_pf[f]][i + posg] * vf[f][i];
	    *it += gvf[f][i][k] * (val_e[func_pf[f]])[i + posv];
	  }
	}
      }
    }
  }
  
  void Xfem::real_hess_base_value(const fem_interpolation_context &,
				  base_tensor &, bool) const
  { GMM_ASSERT1(false, "Sorry order 2 derivatives for Xfem to be done."); }

  void Xfem::init(pfem pf) {
    cvr = pf->ref_convex(0);
    dim_ = cvr->structure()->dim();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = pf->target_dim();
  }
  
  Xfem::Xfem(pfem pf) : pfb(pf), is_valid(false), nb_func(0) {
    if (pf) init(pfb);
  }

}  /* end of namespace getfem.                                            */
