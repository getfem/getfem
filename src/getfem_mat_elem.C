/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_mat_elem.C : computation of elementary matrices.      */
/*     									   */
/*                                                                         */
/* Date : December 21, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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


#include <deque>
#include <getfem_mat_elem.h>
#include <getfem_precomp.h>

namespace getfem
{
  /* ********************************************************************* */
  /*       Elementary matrices computation.                                */
  /* ********************************************************************* */

  struct _emelem_comp_light {
    pmat_elem_type pmt;
    pintegration_method ppi;
    bgeot::pgeometric_trans pgt;
    bool operator < (const _emelem_comp_light &ls) const {
      if (pmt < ls.pmt) return true; if (pmt > ls.pmt) return false; 
      if (ppi < ls.ppi) return true; if (ppi > ls.ppi) return false; 
      if (pgt < ls.pgt) return true; return false;
    }
    _emelem_comp_light(pmat_elem_type pm, pintegration_method pi,
		       bgeot::pgeometric_trans pg)
    { pmt = pm; ppi = pi; pgt = pg; }
    _emelem_comp_light(void) { }
  };


  struct _emelem_comp_structure : public mat_elem_computation
  {
    pgeotrans_precomp pgp;
    ppoly_integration ppi;
    papprox_integration pai;
    bool is_ppi;
    std::vector<base_tensor> mref;
    std::vector<pfem_precomp> pfp;
    std::vector<base_tensor> elmt_stored;
    short_type nbf, dim; 
    base_matrix K, CS, TMP1, B, Htau, M, B2, B3, B32;
    std::deque<short_type> grad_reduction, hess_reduction, trans_reduction;
    std::deque<pfem> trans_reduction_pfi;
    base_vector un, up;
    bool faces_computed;
    bool volume_computed;
    bool is_linear;

    size_type memsize() const {
      size_type sz = sizeof(_emelem_comp_structure) +
	mref.capacity()*sizeof(base_tensor) +
	grad_reduction.size()*sizeof(short_type) +
	hess_reduction.size()*sizeof(short_type) +
	trans_reduction.size()*sizeof(short_type) +
	trans_reduction_pfi.size()*sizeof(pfem) + 
	K.memsize() + CS.memsize()+ TMP1.memsize()+ 
	B.memsize()+ Htau.memsize()+ M.memsize()+ 
	B2.memsize()+ B3.memsize()+ B32.memsize() +
	un.memsize() + up.memsize();

      for (size_type i=0; i < mref.size(); ++i) sz += mref[i].memsize();
      return sz;
    }


    void pre_tensors_for_linear_trans(bool volumic) {

      if ((volumic && volume_computed) || (!volumic && faces_computed)) return;
      scalar_type exectime = ftool::uclock_sec();

      bgeot::multi_index mi(pme->mi.size()), sizes = pme->mi;
      bgeot::multi_index::iterator mit = sizes.begin(), mite = sizes.end();
      size_type f = 1;
      for ( ; mit != mite; ++mit, ++f) f *= *mit;
      if (f > 1000000)
	DAL_WARNING(2, "Warning, very large elementary computations.\n" 
		    << "Be sure you need to compute this elementary matrix.\n"
		    << "(sizes = " << sizes << " )\n");

      base_tensor aux(sizes);
      std::fill(aux.begin(), aux.end(), 0.0);
      if (volumic) {
	volume_computed = true;
	mref[0] = aux;
      }
      else {
	faces_computed = true;
	std::fill(mref.begin()+1, mref.end(), aux);
      }

      if (volumic) cout << "precompute volumic\n"; else cout << "precompute faces\n";

      if (is_ppi)
      {
	base_poly P(dim, 0), Q;
	for ( ; !mi.finished(sizes); mi.incrementation(sizes)) {
	  P.one();
	  mat_elem_type::const_iterator it = pme->begin(), ite = pme->end(); 
	  mit = mi.begin();
	  
	  for (; it != ite; ++it) {
	    size_type ind = *mit; ++mit;

	    if ((*it).pfi->target_dim() > 1)
	      { ind += (*it).pfi->nb_base() * (*mit); ++mit; }
	    Q = ((ppolyfem)((*it).pfi))->base()[ind];

	    switch ((*it).t) {
	    case GETFEM__GRAD    : Q.derivative(*mit); ++mit; break;
	    case GETFEM__HESSIAN :
	      Q.derivative(*mit % dim); Q.derivative(*mit / dim);
	      ++mit; break;
	    case GETFEM__BASE : break;
	    }
	    P *= Q;

	  }
	  if (volumic) mref[0](mi) = ppi->int_poly(P);
	  for (f = 0; f < nbf && !volumic; ++f)
	    mref[f+1](mi) = ppi->int_poly_on_face(P, f);
	}
      }
      else
      { // very very very inefficient ... !!

	// il faudrait faire le calcul par point d'intégration et sommer ...
	// le calcul par point se fait en calc des produit tensoriel d'ordre 1
	// à copier sur le futur calcul en transformation non linéaire
	scalar_type V;

	for ( ; !mi.finished(sizes); mi.incrementation(sizes)) {
	  size_type ind_l = 0, nb_ptc = pai->nb_points_on_convex(), 
	    nb_pt_l = nb_ptc, nb_pt_tot =(volumic ? nb_ptc : pai->nb_points());
	  for (size_type ip = (volumic ? 0:nb_ptc); ip < nb_pt_tot; ++ip) {
	      while (ip == nb_pt_l
		     && ind_l < nbf) {
		nb_pt_l += pai->nb_points_on_face(ind_l);
		ind_l++;
	      }
	      V = 1.0;
	      
	      mat_elem_type::const_iterator it = pme->begin(), ite= pme->end();
	      mit = mi.begin();

	      for (size_type k = 0; it != ite; ++it, ++k) { 
		size_type ind = *mit; ++mit;
		if ((*it).pfi->target_dim() > 1)
		  { ind += (*it).pfi->nb_base() * (*mit); ++mit; }
		switch ((*it).t) {
		case GETFEM__BASE    :
		  V *= (pfp[k]->val(ip))[ind]; break;
		case GETFEM__GRAD    :
		  V *= (pfp[k]->grad(ip))[ind + (*it).pfi->nb_base() *
			 (*it).pfi->target_dim() * (*mit)];
		  ++mit; break;
		case GETFEM__HESSIAN :
		  V *= (pfp[k]->hess(ip))[ind + (*it).pfi->nb_base()
					  * (*it).pfi->target_dim() * (*mit)];
		  ++mit; break;
		}
	      }
	      mref[ind_l](mi) += V * pai->coeff(ip);
	    }
	}
      }

      cout << "precompute Mat elem computation time : "
	 << ftool::uclock_sec() - exectime << endl;

    }

    _emelem_comp_structure(const _emelem_comp_light &ls) {
      
      pgt = ls.pgt;
      pgp = geotrans_precomp(ls.pgt, &(ls.ppi->integration_points()));
      pme = ls.pmt;
      ppi = ls.ppi->method.ppi;
      pai = ls.ppi->method.pai;
      is_ppi = ls.ppi->is_ppi;
      faces_computed = volume_computed = false;
      is_linear = pgt->is_linear();
      nbf = pgt->structure()->nb_faces();
      dim = pgt->structure()->dim();
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
      
      for (size_type k = 0; it != ite; ++it, ++k) {

	if (is_ppi && (!((*it).pfi->is_polynomial()) || !is_linear))
	  DAL_THROW(std::invalid_argument, 
		    "Exact integration not allowed in this context");
	if((*it).pfi->basic_structure() != pgt->basic_structure())
	  DAL_THROW(std::invalid_argument, "incorrect computation");
	
	if (!((*it).pfi->is_equivalent())) {
	  trans_reduction.push_back(k);
	  trans_reduction_pfi.push_back((*it).pfi);
	}
	switch ((*it).t) {
	case GETFEM__BASE    : break;
	case GETFEM__GRAD    : ++k;
	  if ((*it).pfi->do_grad_reduction()) grad_reduction.push_back(k);
	  break;
	case GETFEM__HESSIAN : ++k; hess_reduction.push_back(k); break;
	}
      }

      if (!is_ppi) {
	pfp.resize(pme->size());
	it = pme->begin(), ite = pme->end();
	for (size_type k = 0; it != ite; ++it, ++k) {
	  pfp[k] = fem_precomp((*it).pfi, &(pai->integration_points()));
	}
	elmt_stored.resize(pme->size());
      }
      if (is_linear) mref.resize(nbf + 1);
    }

    // compute pour le cas lineaire uniquement pour le moment
    void compute(base_tensor &t, const base_matrix &G, size_type ir)
    {
      
      dim_type P = dim, N = G.nrows();
      short_type NP = pgt->nb_points();
      scalar_type J;
      if (G.ncols() != NP) DAL_THROW(dimension_error, "dimensions mismatch");
      
      K.resize(N, P); CS.resize(P, P); TMP1.resize(P, P); B.resize(N, P);
      if (hess_reduction.size() > 0) {
	B2.resize(P*P, P); B3.resize(N*N, P*P); Htau.resize(N, P*P);
	B32.resize(N*N, P); B2.fill(0.0);
      }
      if (ir > 0) {
	un.resize(P); up.resize(N);
	un = pgt->normals()[ir-1];
      }
      
      if (is_linear) { // linear transformation
      
      // if (is_ppi) {
	
	pre_tensors_for_linear_trans(ir == 0);
	base_tensor taux;
	
	// on peut simplifier les calculs pour N = P
	// cout << "mat G : " << G << endl;
	// cout << "mat grad : " << pgp->grad(ir) << endl;
	bgeot::mat_product(G, pgp->grad(0), K);
	bgeot::mat_product_tn(K, K, CS);
	J = ::sqrt(bgeot::mat_inv_cholesky(CS, TMP1));
	bgeot::mat_product(K, CS, B);
	
	if (ir > 0) {
	  bgeot::mat_vect_product(B, un, up);
	  J *= bgeot::vect_norm2(up);
	}
	
	bool flag = false;
	t = mref[ir]; t *= J;
	
	if (hess_reduction.size() > 0) {
	  for (short_type i = 0; i < P; ++i)
	    for (short_type j = 0; j < P; ++j)
	      for (short_type k = 0; k < N; ++k)
		for (short_type l = 0; l < N; ++l)
		  B3(k + N*l, i + P*j) = B(k, i) * B(l, j);
	}
	
	if (grad_reduction.size() > 0) {
	  std::deque<short_type>::const_iterator it = grad_reduction.begin(),
	    ite = grad_reduction.end();
	  for ( ; it != ite; ++it) {
	    (flag ? t:taux).mat_transp_reduction(flag ? taux:t, B, *it);
	    flag = !flag;
	  }
	}
	
	if (hess_reduction.size() > 0) {
	  std::deque<short_type>::const_iterator it = hess_reduction.begin(),
	    ite = hess_reduction.end();
	  for (short_type l = 1; it != ite; ++it, l *= 2) {
	    (flag ? t:taux).mat_transp_reduction(flag ? taux:t, B3, *it);
	    flag = !flag;
	  }
	}
	
	/* Applying linear transformation for non tau-equivalent elements.   */
	
	if (trans_reduction.size() > 0) {
	  std::deque<short_type>::const_iterator it = trans_reduction.begin(),
	    ite = trans_reduction.end();
	  std::deque<pfem>::const_iterator iti = trans_reduction_pfi.begin();
	  for ( ; it != ite; ++it, ++iti) { 
	    if ((*iti)->nb_dof() != M.nrows() || (*iti)->nb_base()!=M.ncols())
	      M.resize((*iti)->nb_base(), (*iti)->nb_dof());
	    (*iti)->mat_trans(M, G, pgt);
	    (flag ? t:taux).mat_reduction(flag ? taux:t, M, *it);
	    flag = !flag;
	  }
	}
	
	if (flag) t = taux;
      }
      else { // non linear transformation

	// cout << "on commnce ...\n";

	bool first = true;

	for (size_type ip=(ir == 0) ? 0 : pai->repart()[ir-1];
	     ip < pai->repart()[ir]; ++ip, first = false) {

	  bgeot::mat_product(G, pgp->grad(ip), K);
	  bgeot::mat_product_tn(K, K, CS);
	  J = ::sqrt(bgeot::mat_inv_cholesky(CS, TMP1));
	  bgeot::mat_product(K, CS, B);
	  
	  if (ir > 0) {
	    bgeot::mat_vect_product(B, un, up);
	    J *= bgeot::vect_norm2(up);
	  }

	  if (hess_reduction.size() > 0) {
	    for (short_type i = 0; i < P; ++i)
	      for (short_type j = 0; j < P; ++j)
		for (short_type k = 0; k < N; ++k)
		  for (short_type l = 0; l < N; ++l)
		    B3(k + N*l, i + P*j) = B(k, i) * B(l, j);
	    bgeot::mat_product(G, pgp->hessian(ip), Htau);
	    for (short_type i = 0; i < P; ++i)
	      for (short_type j = 0; j < P; ++j)
		for (short_type k = 0; k < P; ++k)
		  for (short_type l = 0; l < N; ++l)
		    B2(i + P*j, k) += Htau(l, i + P*j) * B(l,k);
	    bgeot::mat_product(B3, B2, B32);
	  }

	  mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
	  bgeot::multi_index mi(pme->mi.size()), sizes = pme->mi;
	  // cout << "sizes : " << sizes << endl;
	  bgeot::multi_index::iterator mit = sizes.begin();
	  
	  for (size_type k = 0; it != ite; ++it, ++k) {
	    ++mit; if ((*it).pfi->target_dim() > 1) ++mit;

	    switch ((*it).t) {
	    case GETFEM__BASE    :
	      elmt_stored[k] = pfp[k]->val(ip);  break;
	    case GETFEM__GRAD    :
	      if ((*it).pfi->do_grad_reduction()) {
		elmt_stored[k].mat_transp_reduction(pfp[k]->grad(ip), B, 2);
		*mit++ = N;
	      }
	      else
		elmt_stored[k] = pfp[k]->grad(ip);
	      // a verif
	      break;
	    case GETFEM__HESSIAN :
	      base_tensor tt = pfp[k]->hess(ip);
	      bgeot::multi_index mim(3);
	      mim[2] = dal::sqr(tt.sizes()[2]); mim[1] = tt.sizes()[1];
	      mim[0] = tt.sizes()[0];
	      tt.adjust_sizes(mim);
	      elmt_stored[k].mat_transp_reduction(tt, B3, 2);
	      tt.mat_transp_reduction(pfp[k]->grad(ip), B32, 2);
	      elmt_stored[k] -= tt;
	      *mit++ = N*N; 
	      break;
	    }
	  }
	  
	  if (first) {
	    t.adjust_sizes(sizes);
	    std::fill(t.begin(), t.end(), 0.0);
	  }
	  base_tensor::iterator pt = t.begin();
	  std::vector<base_tensor::iterator> pts(pme->size());


	  J *= pai->coeff(ip);
	  scalar_type V = J;
	  size_type k;
	  for (k = 0; k < pme->size(); ++k) {
	    pts[k] = elmt_stored[k].begin();
	    if (k != 0) V *= *(pts[k]); 
	    // cout << "tensor " << k << " : " << elmt_stored[k] << endl;
	  }

	  for (;;) {
	    *pt += V * (*(pts[0])); pt++; (pts[0])++;
	    if (pts[0] == elmt_stored[0].end()) {
	      pts[0] = elmt_stored[0].begin();
	      for  (k = 1; k < pme->size(); ++k) {
		(pts[k])++; 
		if (pts[k] == elmt_stored[k].end())
		  pts[k] = elmt_stored[k].begin();
		else break;
	      }
	      if (k == pme->size()) break;
	      V = J;
	      for  (k = 1; k < pme->size(); ++k) V *= *(pts[k]);
	    }
	  }
	  if (pt != t.end()) DAL_THROW(internal_error, "Internal error");
	}
      }
    }
    

    void compute(base_tensor &t, const base_matrix &G)
    { compute(t, G, 0); }

    void compute_on_face(base_tensor &t, const base_matrix &G, short_type f)
    { compute(t, G, f+1); }

  };

  static dal::FONC_TABLE<_emelem_comp_light, _emelem_comp_structure>
    *_tab__mat_elet = 0;
  
  pmat_elem_computation mat_elem(pmat_elem_type pm, pintegration_method pi,
				 bgeot::pgeometric_trans pg)
  { 
    if (_tab__mat_elet == 0)
      _tab__mat_elet = 
	new dal::FONC_TABLE<_emelem_comp_light, _emelem_comp_structure>();
    return _tab__mat_elet->add(_emelem_comp_light(pm, pi, pg));
  }


}  /* end of namespace getfem.                                            */

