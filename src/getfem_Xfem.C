/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_Xfem.C : definition of eXtended fems.                  */
/*           see for instance  "A finite element method for crack growth   */
/*           without remeshing", N. Moës, J. Dolbow and T. Belytschko,     */
/*           Int. J. Num. Meth. Engng. 46, 131-150 (1999).                 */
/*                                                                         */
/* Date : April 8, 2003.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#include <getfem_Xfem.h>

namespace getfem
{
  void Xfem::valid(void) {
    init_cvs_node();
    /* setup nodes of the base fem */
    for (size_type k = 0; k < pfb->nb_base(); ++k)
      add_node(pfb->dof_types()[k], pfb->node_of_dof(k));
    
    /* setup nodes of the enriched fems */
    for (size_type k = 0; k < nb_func; ++k) {
      for (size_type j = 0; j < pfe(k)->nb_base(); ++j) {
	add_node(xfem_dof(pfe(k)->dof_types()[j], func_indices[k]),
		 pfe(k)->node_of_dof(j));
      }
    }
    is_valid = true;
  }
  
  size_type Xfem::nb_dof(void) const {
    if (!is_valid)
      DAL_THROW(failure_error, "Valid the Xfem element before using it");
    return _dof_types.size();
  }

  void Xfem::add_func(pfem pf, pXfem_func pXf, size_type ind) {
    nb_func ++;
    if (ind == size_type(-1)) ind = nb_func;
    funcs.resize(nb_func);
    func_indices.resize(nb_func);
    funcs[nb_func-1] = pXf;
    if (cvr != pf->ref_convex() || pfb->target_dim() != pf->target_dim())
      DAL_THROW(failure_error, "Incompatible Xfem fems");

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
  
  // Interpolation method : call the interpolation of pfi for ordinary
  // base function and for each additional function and sum the
  // contributions.
  void Xfem::interpolation(const base_node &x, const base_matrix &G,
			   bgeot::pgeometric_trans pgt,
			   const base_vector &coeff, base_node &val) const {
    Xfem_func_context ctx(pgt,x,G);
    base_node val2(val.size());
    pfb->interpolation(x, G, pgt, coeff, val);
    base_vector coeff2;
    for (size_type k = 0, dofcnt = pfb->nb_base(); k < nb_func; ++k) {
      coeff2.resize(pfe(k)->nb_base());
      ctx.pf = pfe(k); 
      for (ctx.base_num = 0; ctx.base_num != pfe(k)->nb_base(); 
	   ++ctx.base_num, ++dofcnt) {
	coeff2[ctx.base_num] = coeff[dofcnt] * funcs[k]->val(ctx);
      }
      pfe(k)->interpolation(x, G, pgt, coeff2, val2);
      val += val2;
    }
  }

  void Xfem::get_fem_precomp_tab(pfem_precomp pfp, std::vector<pfem_precomp>& vpfp) const {
    vpfp.resize(uniq_pfe.size());
    for (size_type k=0; k < uniq_pfe.size(); ++k) 
      vpfp[k] = (uniq_pfe[k] == pfb) ? pfp : fem_precomp(uniq_pfe[k], pfp->get_point_tab());
  }

  // take into account the fact that the pfp is the same for pfi and
  // the Xfem.
  void Xfem::interpolation(pfem_precomp pfp, size_type ii,
			   const base_matrix &G,
			   bgeot::pgeometric_trans pgt, 
			   const base_vector &coeff, 
			   base_node &val, dim_type Qdim) const {
    size_type Qmult = size_type(Qdim) / target_dim();
    base_node val2(val.size());
    Xfem_func_context ctx(pgt,(*(pfp->get_point_tab()))[ii],G);
    pfb->interpolation(pfp, ii, G, pgt, coeff, val, Qdim);
    base_vector coeff2;
    std::vector<pfem_precomp> vpfp; get_fem_precomp_tab(pfp,vpfp);
    for (size_type k = 0, dofcnt = pfb->nb_base(); k < nb_func; ++k) {
      coeff2.resize(pfe(k)->nb_base() * Qmult);
      ctx.pf = pfe(k);
      for (size_type j = 0; j != pfe(k)->nb_base(); ++j, ++dofcnt) {
	ctx.base_num = j; 
	scalar_type a = funcs[k]->val(ctx);
	for (dim_type q = 0; q < Qmult; ++q)
	  coeff2[j*Qmult + q] = coeff[dofcnt*Qmult + q] * a;
      }
      pfe(k)->interpolation(vpfp[func_pf[k]], ii, G, pgt, coeff2, val2, Qdim);
      val += val2;
    }
  }
  
  void Xfem::interpolation_grad(const base_node &x,
				const base_matrix &G,
				bgeot::pgeometric_trans pgt,
				const base_vector &coeff,
				base_matrix &val) const {
    base_matrix val2(val.nrows(), val.ncols());
    base_node val3(ntarget_dim);
    Xfem_func_context ctx(pgt,x,G);
    pfb->interpolation_grad(x, G, pgt, coeff, val);
    base_vector coeff2;
    // func * grad(base)
    for (size_type k = 0, dofcnt = pfb->nb_base(); k < nb_func; ++k) {
      coeff2.resize(pfe(k)->nb_base());
      ctx.pf = pfe(k);
      for (size_type j = 0; j != pfe(k)->nb_base(); ++j, ++dofcnt) {
	ctx.base_num = j;
	scalar_type a = funcs[k]->val(ctx);
	coeff2[j] = coeff[dofcnt] * a;
      }
      pfe(k)->interpolation_grad(x, G, pgt, coeff2, val2);
      gmm::add(val2, val);
    }
    // grad(func) * base
    for (size_type k = 0, dofcnt = pfb->nb_base(); k < nb_func; ++k) {
      coeff2.resize(pfe(k)->nb_base());
      ctx.pf = pfe(k);
      for (size_type q = 0; q < G.nrows(); ++q) {
	for (size_type j = 0; j != pfe(k)->nb_base(); ++j, ++dofcnt) {
	  ctx.base_num = j;
	  base_vector v = funcs[k]->grad(ctx);
	  coeff2[j] = coeff[dofcnt] * v[q];
	}
	pfe(k)->interpolation(x, G, pgt, coeff2, val3);
	for (dim_type r = 0; r < ntarget_dim; ++r) val2(r, q) = val3[r];
      }
      gmm::add(val2, val);
    }
  }
  
  void Xfem::base_value(const base_node &x, base_tensor &t) const
  { pfb->base_value(x, t); }
  void Xfem::grad_base_value(const base_node &x, base_tensor &t) const
  { pfb->grad_base_value(x, t); }
  void Xfem::hess_base_value(const base_node &x, base_tensor &t) const
  { pfb->hess_base_value(x, t); }

  void Xfem::real_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			     size_type ii, const base_matrix &G,
			     base_tensor &t, size_type) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base();
    t.adjust_sizes(mi);
    scalar_type a;
    Xfem_func_context ctx(pgp->get_trans(), (*(pfp->get_point_tab()))[ii], G);
    //base_node xreal = pgp->transform(ii, G);
    base_tensor::iterator it = t.begin();
    base_tensor::const_iterator itf = pfp->val(ii).begin();
    std::vector<pfem_precomp> vpfp; get_fem_precomp_tab(pfp,vpfp);

    for (dim_type q = 0; q < target_dim(); ++q) {
      for (size_type i = 0; i < pfb->nb_base(); ++i, ++itf, ++it)
	  *it = *itf;

      for (size_type k = 0; k < nb_func; ++k) {
	const base_tensor& val_e = vpfp[func_pf[k]]->val(ii);
	ctx.pf = pfe(k);
	for (size_type i = 0; i < pfe(k)->nb_base(); ++i, ++it) {
	  ctx.base_num = i; a = funcs[k]->val(ctx);
	  *it = val_e[i + q*pfe(k)->nb_base()] * a;
	}
      }
    }
  }

  void Xfem::real_grad_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
				  size_type ii, const base_matrix &G,
				  const base_matrix &B, base_tensor &t,
				  size_type) const {
    bgeot::multi_index mi(3);
    dim_type n = G.nrows();
    mi[2] = n; mi[1] = target_dim(); mi[0] = nb_base();
    t.adjust_sizes(mi);
    
    Xfem_func_context ctx(pgp->get_trans(), (*(pfp->get_point_tab()))[ii], G);
    base_tensor tt; tt.mat_transp_reduction(pfp->grad(ii), B, 2);

    base_tensor::iterator it = t.begin();
    base_tensor::const_iterator itvf = tt.begin();

    std::vector<pfem_precomp> vpfp; get_fem_precomp_tab(pfp,vpfp);
    std::vector<const base_tensor*> val_e(nb_func);
    std::vector<base_tensor> grad_e(nb_func);
    for (size_type i=0; i < uniq_pfe.size(); ++i) {
      val_e[i] = &vpfp[i]->val(ii);
      grad_e[i].mat_transp_reduction(vpfp[i]->grad(ii), B, 2);
    }
    std::vector<std::vector<scalar_type> > vf(nb_func);
    std::vector<std::vector<base_vector> > gvf(nb_func);
    for (size_type f = 0; f < nb_func; ++f) {
      vf[f].resize(pfe(f)->nb_base());
      gvf[f].resize(pfe(f)->nb_base());
      ctx.pf = pfe(f);
      for (ctx.base_num=0; ctx.base_num < pfe(f)->nb_base(); ++ctx.base_num) {
	vf[f][ctx.base_num] = funcs[f]->val(ctx); 
	gvf[f][ctx.base_num] = funcs[f]->grad(ctx); 
      }
    }

    //    cerr << "pfp->val(ii)={"; 
    //    for (size_type i=0; i < pfp->val(ii).size(); ++i) cerr << pfp->val(ii)[i] << " "; cerr << "}\n";
    
    for (dim_type k = 0; k < n ; ++k) {      
      for (dim_type q = 0; q < target_dim(); ++q) {
	for (size_type i = 0; i < pfb->nb_base(); ++i, ++it)
	    *it = *itvf++;
	for (size_type f = 0; f < nb_func; ++f) {
          size_type posg = pfe(f)->nb_base()*(q + k*target_dim());
          size_type posv = pfe(f)->nb_base()*q;
	  for (size_type i = 0; i < pfe(f)->nb_base(); ++i, ++it) {
	    *it = grad_e[func_pf[f]][i + posg] * vf[f][i];
	    *it += gvf[f][i][k] * (*val_e[func_pf[f]])[i + posv];
	  }
	}
      }
    }
  }
  
  void Xfem::real_hess_base_value(pgeotrans_precomp, pfem_precomp,
				  size_type, const base_matrix &,
				  const base_matrix &, const base_matrix &,
				  base_tensor &, size_type) const {
    DAL_THROW(to_be_done_error,
	      "Sorry order 2 derivatives for Xfem to be done.");
  }
  
  Xfem::Xfem(pfem pf) : pfb(pf), is_valid(false), nb_func(0) {
    if (!(pfb->is_equivalent()))
      DAL_THROW(to_be_done_error,
		"Sorry, Xfem for non tau-equivalent elements to be done.");
    cvr = pfb->ref_convex();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5; /* humm ... */
    ntarget_dim = pfb->target_dim();
  }


}  /* end of namespace getfem.                                            */
