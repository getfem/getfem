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

  static virtual_Xfem_hess no_Xfem_hess_defined;
  pXfem_hess pno_Xfem_hess_defined = &no_Xfem_hess_defined;

  base_matrix virtual_Xfem_hess::operator()(const base_node &) {
    DAL_THROW(failure_error, "No hessian defined for this function");
  }

  void Xfem::valid(void) {
    init_cvs_node();
    for (size_type k = 0; k < pfi->nb_base(); ++k)
      add_node(pfi->dof_types()[k], pfi->node_of_dof(k));
    
    for (size_type k = 0; k < nb_func; ++k) {
      for (size_type j = 0; j < pfi->nb_base(); ++j) {
	add_node(xfem_dof(pfi->dof_types()[j], func_indices[k]),
		 pfi->node_of_dof(j));
      }
    }
    cout << "valid ok\n";
    is_valid = true;
  }
  
  size_type Xfem::nb_dof(void) const {
    if (!is_valid) DAL_THROW(failure_error,
			     "Valid the Xfem element before using it");
    return _dof_types.size();
  }

  void Xfem::add_func(pXfem_func pXf, pXfem_grad pXg, pXfem_hess pXh,
		      size_type ind) {
    nb_func ++;
    if (ind == size_type(-1)) ind = nb_func;
    funcs.resize(nb_func); grads.resize(nb_func); hess.resize(nb_func);
    func_indices.resize(nb_func);
    funcs[nb_func-1] = pXf;
    grads[nb_func-1] = pXg;
    hess[nb_func-1] = pXh;
    func_indices[nb_func-1] = ind;
    is_valid = false;
  }
  
  // Interpolation method : call the interpolation of pfi for ordinary
  // base function and for each additional function and sum the
  // contributions.
  void Xfem::interpolation(const base_node &x, const base_matrix &G,
			   bgeot::pgeometric_trans pgt,
			   const base_vector &coeff, base_node &val) const {
    base_node val2(val.size());
    // cerr << "coucou Xfem::interpolation(x=" << x << ", G=" << G << endl;
    base_node xreal = pgt->transform(x, G);
    size_type nbb = pfi->nb_base();
    pfi->interpolation(x, G, pgt, coeff, val);
    base_vector coeff2(nbb);
    for (size_type i = 0; i < nb_func; ++i) {
      scalar_type a = (*(funcs[i]))(xreal);
      for (size_type j = 0; j != nbb; ++j) {
	coeff2[j] = coeff[(i+1) * nbb + j] * a;
      }
      pfi->interpolation(x, G, pgt, coeff2, val2);
      val += val2;
    }
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
    base_node xreal = pgt->transform((*(pfp->get_point_tab()))[ii], G);
    size_type nbb = pfi->nb_base();
    pfi->interpolation(pfp, ii, G, pgt, coeff, val, Qdim);
    base_vector coeff2(pfi->nb_base());
    for (size_type i = 0; i < nb_func; ++i) {
      scalar_type a = (*(funcs[i]))(xreal);
      for (size_type j = 0; j != nbb; ++j) {
	for (dim_type q = 0; q < Qmult; ++q)
	  coeff2[j*Qmult + q] = coeff[((i+1)*nbb + j)* Qmult + q] * a;
      }
      pfi->interpolation(pfp, ii, G, pgt, coeff2, val2, Qdim);
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
    base_node xreal = pgt->transform(x, G);
    size_type nbb = pfi->nb_base();
    pfi->interpolation_grad(x, G, pgt, coeff, val);
    base_vector coeff2(nbb);
    for (size_type i = 0; i < nb_func; ++i) {
      scalar_type a = (*(funcs[i]))(xreal);
      for (size_type j = 0; j != nbb; ++j) {
	coeff2[j] = coeff[(i+1) * nbb + j] * a;
      }
      pfi->interpolation_grad(x, G, pgt, coeff2, val2);
      val += val2;
    }
    for (size_type i = 0; i < nb_func; ++i) {
      base_vector v = (*(grads[i]))(xreal);
      for (size_type q = 0; q < G.nrows(); ++q) {
	for (size_type j = 0; j != nbb; ++j) {
	  coeff2[j] = coeff[(i+1) * nbb + j] * v[q];
	}
	pfi->interpolation(x, G, pgt, coeff2, val3);
	for (dim_type r = 0; r < ntarget_dim; ++r) val2(r, q) = val3[r];
      }
      val += val2;
    }
  }
  
  void Xfem::base_value(const base_node &x, base_tensor &t) const
  { pfi->base_value(x, t); }
  void Xfem::grad_base_value(const base_node &x, base_tensor &t) const
  { pfi->grad_base_value(x, t); }
  void Xfem::hess_base_value(const base_node &x, base_tensor &t) const
  { pfi->hess_base_value(x, t); }

  void Xfem::real_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			     size_type ii, const base_matrix &G,
			     base_tensor &t) const {
    bgeot::multi_index mi(2);
    mi[1] = target_dim(); mi[0] = nb_base();
    t.adjust_sizes(mi);
    scalar_type a;
    base_node xreal = pgp->transform(ii, G);
    base_tensor::iterator it = t.begin();
    base_tensor::const_iterator itf = pfp->val(ii).begin(), itf2;
    
    for (dim_type q = 0; q < target_dim(); ++q) {
      for (size_type f = 0; f <= nb_func; ++f) {
	if (f == 0) a = scalar_type(1); else a = (*(funcs[f-1]))(xreal);
	itf2 = itf;
	for (size_type i = 0; i < pfi->nb_base(); ++i, ++itf2, ++it)
	  *it = *itf2 * a;
      }
      itf = itf2;
    }
  }

  void Xfem::real_grad_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
				  size_type ii, const base_matrix &G,
				  const base_matrix &B, base_tensor &t) const {
    
    bgeot::multi_index mi(3);
    dim_type n = G.nrows();
    mi[2] = n; mi[1] = target_dim(); mi[0] = nb_base();
    t.adjust_sizes(mi);
    base_vector v(n);
    scalar_type a;
    
    base_node xreal = pgp->transform(ii, G);
    base_tensor tt; tt.mat_transp_reduction(pfp->grad(ii), B, 2);

    base_tensor::iterator it = t.begin();
    base_tensor::const_iterator itvf = tt.begin(), itvf2;
    base_tensor::const_iterator itf = pfp->val(ii).begin(), itf2;
    
    //    cerr << "pfp->val(ii)={"; 
    //    for (size_type i=0; i < pfp->val(ii).size(); ++i) cerr << pfp->val(ii)[i] << " "; cerr << "}\n";

    std::vector<scalar_type> vf(nb_func);
    std::vector<base_vector> gvf(nb_func);
    for (size_type f = 0; f < nb_func; ++f)
      { vf[f] = (*(funcs[f]))(xreal); gvf[f] = (*(grads[f]))(xreal); }

    for (dim_type k = 0; k < n ; ++k) {
      itf = pfp->val(ii).begin();
      for (dim_type q = 0; q < target_dim(); ++q) {
	for (size_type f = 0; f <= nb_func; ++f) {
	  if (f == 0) { a = scalar_type(1); v.fill(1); }
	  else { a = vf[f-1]; v = gvf[f-1]; }
	  itvf2 = itvf; itf2 = itf;
	  for (size_type i = 0; i < pfi->nb_base(); ++i, ++it) {
	    /*	    cerr << "Xfem::real_grad_base_value(ii=" << ii << "): it=" << &(*it) << "=" << *it;
	    cerr << ", f=" << f << ", k=" << k << ", v[k]=" << v[k] << ", t.begin=" << &t[0] << ", tt.begin=" << &tt[0] 
		 << ", pfp->val(ii).size()=" << pfp->val(ii).size() << "tt.size()=" << tt.size() << ", t.size()=" << t.size() << endl;
	    cerr << "itf=" << itf - pfp->val(ii).begin() << ", itf2=" << itf2 - pfp->val(ii).begin() << ", itvf=" << itvf - tt.begin() << ", itvf2=" << itvf2-tt.begin() << endl;
	    cerr << "*itf=" << *itf << ", *itf2=" << *itf2 << ", *itvf=" << *itvf << ", *itvf2=" << *itvf2 << endl;*/
	    *it = *itvf2++ * a;
	    if (f > 0) *it += v[k] * (*itf2++);
	  }
	}
	itvf = itvf2; itf = itf2; 
      }
    }
  }
  
  void Xfem::real_hess_base_value(pgeotrans_precomp, pfem_precomp,
				  size_type, const base_matrix &,
				  const base_matrix &, const base_matrix &,
				  base_tensor &) const {
    DAL_THROW(to_be_done_error,
	      "Sorry order 2 derivatives for Xfem to be done.");
  }
  
  Xfem::Xfem(pfem pf) : pfi(pf), is_valid(false), nb_func(0) {
    if (!(pfi->is_equivalent())) 
      DAL_THROW(to_be_done_error,
		"Sorry, Xfem for non tau-equivalent elements to be done.");
    cvr = pfi->ref_convex();
    is_equiv = real_element_defined = true;
    is_polycomp = is_pol = is_lag = false;
    es_degree = 5;
    ntarget_dim = pfi->target_dim();
  }


}  /* end of namespace getfem.                                            */
