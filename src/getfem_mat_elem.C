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
/* Copyright (C) 2000-2003  Yves Renard.                                   */
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

  struct emelem_comp_light_ {
    pmat_elem_type pmt;
    pintegration_method ppi;
    bgeot::pgeometric_trans pgt;
    bool operator < (const emelem_comp_light_ &ls) const {
      if (pmt < ls.pmt) return true; if (pmt > ls.pmt) return false; 
      if (ppi < ls.ppi) return true; if (ppi > ls.ppi) return false; 
      if (pgt < ls.pgt) return true; return false;
    }
    emelem_comp_light_(pmat_elem_type pm, pintegration_method pi,
		       bgeot::pgeometric_trans pg)
    { pmt = pm; ppi = pi; pgt = pg; }
    emelem_comp_light_(void) { }
  };


  struct emelem_comp_structure_ : public mat_elem_computation
  {
    pgeotrans_precomp pgp;
    ppoly_integration ppi;
    papprox_integration pai;
    bool is_ppi;
    std::vector<base_tensor> mref;
    std::vector<pfem_precomp> pfp;
    std::vector<base_tensor> elmt_stored;
    short_type nbf, dim; 
    base_matrix K, CS, B, Htau, M, B2, B3, B32;
    std::deque<short_type> grad_reduction, hess_reduction, trans_reduction;
    std::deque<pfem> trans_reduction_pfi;
    base_small_vector un, up;
    bool faces_computed;
    bool volume_computed;
    bool is_linear;
    bool computed_on_real_element;

    size_type memsize() const {
      size_type sz = sizeof(emelem_comp_structure_) +
	mref.capacity()*sizeof(base_tensor) +
	grad_reduction.size()*sizeof(short_type) +
	hess_reduction.size()*sizeof(short_type) +
	trans_reduction.size()*sizeof(short_type) +
	trans_reduction_pfi.size()*sizeof(pfem);

      for (size_type i=0; i < mref.size(); ++i) sz += mref[i].memsize();
      return sz;
    }

    emelem_comp_structure_(const emelem_comp_light_ &ls) {
      
      pgt = ls.pgt;
      pgp = geotrans_precomp(ls.pgt, &(ls.ppi->integration_points()));
      pme = ls.pmt;
      ppi = ls.ppi->method.ppi;
      pai = ls.ppi->method.pai;
      is_ppi = ls.ppi->is_ppi;
      faces_computed = volume_computed = false;
      is_linear = pgt->is_linear();
      computed_on_real_element = !is_linear;
      nbf = pgt->structure()->nb_faces();
      dim = pgt->structure()->dim();
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
      
      for (size_type k = 0; it != ite; ++it, ++k) {

	if ((*it).pfi->is_on_real_element()) computed_on_real_element = true;
	if (is_ppi && (!((*it).pfi->is_polynomial()) || !is_linear 
		       || computed_on_real_element))
	  DAL_THROW(std::invalid_argument, 
		    "Exact integration not allowed in this context");
	if((*it).pfi->basic_structure() != pgt->basic_structure())
	  DAL_THROW(std::invalid_argument, "incorrect computation");
	
	if (!((*it).pfi->is_equivalent())) {
	  trans_reduction.push_back(k);
	  trans_reduction_pfi.push_back((*it).pfi);
	}
	switch ((*it).t) {
	case GETFEM_BASE_    : break;
	case GETFEM_GRAD_    : ++k;
	  if (!((*it).pfi->is_on_real_element())) grad_reduction.push_back(k);
	  break;
	case GETFEM_HESSIAN_ : ++k; hess_reduction.push_back(k); break;
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
      if (!computed_on_real_element) mref.resize(nbf + 1);
    }

    void add_elem(base_tensor &t, const base_matrix &G, size_type ip,
		  scalar_type J, dim_type N, size_type elt, 
		  bool first, bool trans = true) {
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
      bgeot::multi_index mi(pme->mi.size()), sizes = pme->mi;
      bgeot::multi_index::iterator mit = sizes.begin();
      
      for (size_type k = 0; it != ite; ++it, ++k) {
	++mit; if ((*it).pfi->target_dim() > 1) ++mit;
	
	switch ((*it).t) {
	case GETFEM_BASE_    :
	  (*it).pfi->real_base_value(pgp, pfp[k], ip, G, elmt_stored[k], elt);
	  break;
	case GETFEM_GRAD_    :
	  if (trans) {
	    (*it).pfi->real_grad_base_value(pgp, pfp[k], ip, G, B, 
					    elmt_stored[k], elt);
	    *mit++ = N;
	  }
	  else
	    elmt_stored[k] = pfp[k]->grad(ip);
	  break;
	case GETFEM_HESSIAN_ :
	  if (trans) {
	    (*it).pfi->real_hess_base_value(pgp, pfp[k], ip, G, B3, 
					    B32, elmt_stored[k], elt);
	    *mit++ = N*N;
	  }
	  else {
	    base_tensor tt = pfp[k]->hess(ip);
	    bgeot::multi_index mim(3);
	    mim[2] = dal::sqr(tt.sizes()[2]); mim[1] = tt.sizes()[1];
	    mim[0] = tt.sizes()[0];
	    tt.adjust_sizes(mim);
	    elmt_stored[k] = tt;
	  }
	  break;
	}
      }
      
      if (first) {
	t.adjust_sizes(sizes);
	std::fill(t.begin(), t.end(), 0.0);
      }

      scalar_type V;
      size_type k;
      base_tensor::iterator pt = t.begin();
      std::vector<base_tensor::const_iterator> pts(pme->size());
      std::vector<scalar_type> Vtab(pme->size());
            
      J *= pai->coeff(ip);
      for (k = 0; k < pme->size(); ++k)
	pts[k] = elmt_stored[k].begin();
      base_tensor::const_iterator pts0 = pts[0];
      
      size_type n0 = elmt_stored[0].size();
      k = pme->size()-1; Vtab[k] = J;
      /* very heavy reduction .. takes much time */
      do {
        for (V = Vtab[k]; k; --k)
          Vtab[k-1] = V = *pts[k] * V;
        for (; k < n0/*elmt_stored[0].size()*/; ++k)
          *pt++ += V * pts0[k];
        for (k=1; k != pme->size() && ++pts[k] == elmt_stored[k].end(); ++k)
          pts[k] = elmt_stored[k].begin();
      } while (k != pme->size());
      if (pt != t.end()) DAL_THROW(internal_error, "Internal error");
    }

    void pre_tensors_for_linear_trans(bool volumic) {

      if ((volumic && volume_computed) || (!volumic && faces_computed)) return;
      // scalar_type exectime = ftool::uclock_sec();

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

      if (is_ppi) // pour accelerer, il faudrait précalculer les dérivées
      {
	base_poly P(dim, 0), Q(dim, 0), R(dim, 0);
	size_type old_ind = size_type(-1), ind; 
	for ( ; !mi.finished(sizes); mi.incrementation(sizes)) {
	  
	  mat_elem_type::const_iterator it = pme->begin(), ite = pme->end(); 
	  mit = mi.begin();

	  ind = *mit; ++mit;

	  if ((*it).pfi->target_dim() > 1)
	    { ind += (*it).pfi->nb_base() * (*mit); ++mit; }
	  Q = ((ppolyfem)((*it).pfi))->base()[ind];

	  switch ((*it).t) {
	  case GETFEM_GRAD_    : Q.derivative(*mit); ++mit; break;
	  case GETFEM_HESSIAN_ :
	    Q.derivative(*mit % dim); Q.derivative(*mit / dim);
	    ++mit; break;
	  case GETFEM_BASE_ : break;
	  }
	  ++it;

	  if (it != ite && *mit != old_ind) {
	    old_ind = *mit; 
	    P.one();
	    for (; it != ite; ++it) {
	      ind = *mit; ++mit;
	      
	      if ((*it).pfi->target_dim() > 1)
		{ ind += (*it).pfi->nb_base() * (*mit); ++mit; }
	      R = ((ppolyfem)((*it).pfi))->base()[ind];
	      
	      switch ((*it).t) {
	      case GETFEM_GRAD_    : R.derivative(*mit); ++mit; break;
	      case GETFEM_HESSIAN_ :
		R.derivative(*mit % dim); R.derivative(*mit / dim);
		++mit; break;
	      case GETFEM_BASE_ : break;
	      }
	      P *= R;   
	    }
	  }
	  R = P * Q;
	  if (volumic) mref[0](mi) = ppi->int_poly(R);
	  for (f = 0; f < nbf && !volumic; ++f)
	    mref[f+1](mi) = ppi->int_poly_on_face(R, f);
	}
      }
      else { 
	bool first = true;
	
	size_type ind_l = 0, nb_ptc = pai->nb_points_on_convex(), 
	  nb_pt_l = nb_ptc, nb_pt_tot =(volumic ? nb_ptc : pai->nb_points());
	for (size_type ip = (volumic ? 0:nb_ptc); ip < nb_pt_tot; ++ip) {
	  while (ip == nb_pt_l && ind_l < nbf)
	    { nb_pt_l += pai->nb_points_on_face(ind_l); ind_l++; }
	  add_elem(mref[ind_l], B, ip, 1.0, 0, 0, first, false);
	  first = false;
	}
      }
      // cout << "precompute Mat elem computation time : "
      //   << ftool::uclock_sec() - exectime << endl;
    }


    void compute(base_tensor &t, const base_matrix &G, size_type ir,
		 size_type elt) {
      dim_type P = dim, N = G.nrows();
      short_type NP = pgt->nb_points();
      scalar_type J;
      if (G.ncols() != NP) DAL_THROW(dimension_error, "dimensions mismatch");
      
      K.resize(N, P); CS.resize(P, P); B.resize(P, N);
      if (hess_reduction.size() > 0) {
	B2.resize(P*P, P); B3.resize(N*N, P*P); Htau.resize(N, P*P);
	B32.resize(N*N, P); B2.fill(0.0);
      }
      if (ir > 0) {
	un.resize(P); up.resize(N);
	un = pgt->normals()[ir-1];
      }
      base_tensor taux;
      bool flag = false;

      if (!computed_on_real_element) {
	
	pre_tensors_for_linear_trans(ir == 0);
	
	// computation of the pseudo inverse
	gmm::mult(gmm::transposed(pgp->grad(0)), gmm::transposed(G), K);
	if (P != N) {
	  gmm::mult(K, gmm::transposed(K), CS);
	  J = ::sqrt(gmm::lu_inverse(CS));
	  gmm::mult(gmm::transposed(K), CS, B);
	}
	else {
	  J = dal::abs(gmm::lu_inverse(K)); B = K;
	}
	
	if (ir > 0) {
	  gmm::mult(B, un, up);
	  J *= bgeot::vect_norm2(up);
	}
     
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
	
      }
      else { // non linear transformation and methods defined on real elements

	bool first = true;

	for (size_type ip=(ir == 0) ? 0 : pai->repart()[ir-1];
	     ip < pai->repart()[ir]; ++ip, first = false) {

	  // computation of the pseudo inverse
	  gmm::mult(gmm::transposed(pgp->grad(ip)), gmm::transposed(G), K);
	  if (P != N) {
	    gmm::mult(K, gmm::transposed(K), CS);
	    J = ::sqrt(gmm::lu_inverse(CS));
	    gmm::mult(gmm::transposed(K), CS, B);
	  }
	  else {
  	    J = dal::abs(gmm::lu_inverse(K)); B = K;
  	  }
	  
	  if (ir > 0) {
	    gmm::mult(B, un, up);
	    J *= bgeot::vect_norm2(up);
	  }

	  if (hess_reduction.size() > 0) {
	    for (short_type i = 0; i < P; ++i)
	      for (short_type j = 0; j < P; ++j)
		for (short_type k = 0; k < N; ++k)
		  for (short_type l = 0; l < N; ++l)
		    B3(k + N*l, i + P*j) = B(k, i) * B(l, j);
	    gmm::mult(G, pgp->hessian(ip), Htau);
	    for (short_type i = 0; i < P; ++i)
	      for (short_type j = 0; j < P; ++j)
		for (short_type k = 0; k < P; ++k)
		  for (short_type l = 0; l < N; ++l)
		    B2(i + P*j, k) += Htau(l, i + P*j) * B(l,k);
	    gmm::mult(B3, B2, B32);
	  }

	  add_elem(t, G,  ip, J, N, elt, first);
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
    
    void compute(base_tensor &t, const base_matrix &G, size_type elt)
    { compute(t, G, 0, elt); }

    void compute_on_face(base_tensor &t, const base_matrix &G,
			 short_type f, size_type elt)
    { compute(t, G, f+1, elt); }

  };

  static dal::FONC_TABLE<emelem_comp_light_, emelem_comp_structure_>
    *tab__mat_elet_ = 0;
  
  pmat_elem_computation mat_elem(pmat_elem_type pm, pintegration_method pi,
				 bgeot::pgeometric_trans pg)
  { 
    if (tab__mat_elet_ == 0)
      tab__mat_elet_ = 
	new dal::FONC_TABLE<emelem_comp_light_, emelem_comp_structure_>();
    return tab__mat_elet_->add(emelem_comp_light_(pm, pi, pg));
  }


}  /* end of namespace getfem.                                            */

