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

  struct _emelem_comp_light
  {
    pmat_elem_type pmt;
    pintegration_method ppi;
    bgeot::pgeometric_trans pgt;
    bool operator < (const _emelem_comp_light &ls) const
    {
      if (pmt < ls.pmt) return true; if (pmt > ls.pmt) return false; 
      if (ppi < ls.ppi) return true; if (ppi > ls.ppi) return false; 
      if (pgt < ls.pgt) return true; return false;
    }
    _emelem_comp_light(pmat_elem_type pm,
		   pintegration_method pi, bgeot::pgeometric_trans pg)
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
    std::vector<short_type>  mref_count;
    std::vector<scalar_type>  mref_coeff;
    std::vector<pfem_precomp> pfp;
    short_type nbf, nhess; 
    size_type nint;
    std::vector< std::vector<short_type> > indv;
    base_matrix K, CS, TMP1, B, Htau, M, B2, B3, B32;
    std::deque<short_type> grad_reduction;
    std::deque<short_type> hess_reduction;
    std::deque<short_type> trans_reduction;
    std::deque<pfem> trans_reduction_pfi;
    base_vector un, up;
    
    inline size_type indcomp(size_type ip, short_type h, short_type c)
    { return ip * nhess * 3 + size_type(h) * 3 + size_type(c); }

    size_type memsize() const {
      size_type sz = sizeof(_emelem_comp_structure) +
	mref.capacity()*sizeof(base_tensor) +
	mref_count.capacity()*sizeof(short_type) +
	mref_coeff.capacity()*sizeof(scalar_type) +
	grad_reduction.size()*sizeof(short_type) +
	hess_reduction.size()*sizeof(short_type) +
	trans_reduction.size()*sizeof(short_type) +
	trans_reduction_pfi.size()*sizeof(pfem) + 
	K.memsize() + CS.memsize()+ TMP1.memsize()+ 
	B.memsize()+ Htau.memsize()+ M.memsize()+ 
	B2.memsize()+ B3.memsize()+ B32.memsize() +
	un.memsize() + up.memsize() +
	indv.capacity() * sizeof(std::vector<short_type>);

      for (size_type i=0; i < mref.size(); ++i) sz += mref[i].memsize();
      for (size_type i=0; i < indv.size(); ++i) sz += indv[i].capacity() * sizeof(short_type);
      return sz;
    }

    _emelem_comp_structure(const _emelem_comp_light &ls)
    { // optimisable ... !!
      pgt = ls.pgt;
      pgp = geotrans_precomp(ls.pgt, &(ls.ppi->integration_points()));
      pme = ls.pmt;
      ppi = ls.ppi->method.ppi;
      pai = ls.ppi->method.pai;
      is_ppi = ls.ppi->is_ppi;
      size_type k;

      // cout << "debut de construction\n";
      
      nbf = pgt->structure()->nb_faces();
      nhess = 1;
      nint = is_ppi ? nbf + 1 : pai->nb_points();

      mat_elem_type::const_iterator it = ls.pmt->begin(), ite = ls.pmt->end();
      
      for (k = 0; it != ite; ++it, ++k) {

	if (is_ppi && (!((*it).pfi->is_polynomial()) || 
		       !(pgt->is_linear())))
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
	case GETFEM__HESSIAN : ++k; hess_reduction.push_back(k);
	  if (!(pgt->is_linear())) nhess *= 2; break;
	}
      }
      
      mref.resize(nhess * nint * 3);
      mref_count.resize(nhess * nint * 3);
      mref_coeff.resize(nhess * nint * 3);
      bgeot::multi_index mi(ls.pmt->mi.size()), sizes = ls.pmt->mi;
      bgeot::multi_index::iterator mit = sizes.begin(), mite = sizes.end();
      for (k = 1; mit != mite; ++mit, ++k) k *= *mit;
      if (k * nhess * nint * 3 > 1000000)
	cerr << "Warning, very large elementary computations.\n"
	     << "Be sure you need to compute this elementary matrix.\n"
	     << "(sizes = " << sizes << " times " << nhess*nint*3 << " )\n";
      std::fill(mref.begin(), mref.end(), base_tensor(sizes));
      
      for (k = 0; k < nhess; ++k)
	for (size_type j = 0; j < nint; ++j)
	  mref[indcomp(j, k, 0)] = base_tensor(sizes);
      
      if (is_ppi)
      {
	base_poly P(ls.pgt->structure()->dim(), 0), Q;
	for ( ; !mi.finished(sizes); mi.incrementation(sizes)) {
	  P.one();
	  it = ls.pmt->begin(), ite = ls.pmt->end();
	  mit = mi.begin();
	  short_type hc = 1;
	  
	  for (; it != ite; ++it) {
	    // cout << "P = " << P << endl;
	    size_type ind = *mit; ++mit;

	    if ((*it).pfi->target_dim() > 1)
	      { ind += (*it).pfi->nb_base() * (*mit); ++mit; }
	    Q = ((ppolyfem)((*it).pfi))->base()[ind];

	    switch ((*it).t) {
	    case GETFEM__BASE    : P *= Q; break;
	    case GETFEM__GRAD    : Q.derivative(*mit); ++mit; P *= Q; break;
	    case GETFEM__HESSIAN : 
	      Q.derivative(*mit % ls.pgt->structure()->dim());
	      Q.derivative(*mit / ls.pgt->structure()->dim());
	      ++mit; P *= Q; hc *= 2; break;
	    }

	  }
	  mref[indcomp(0, 0, 0)](mi) = ppi->int_poly(P);
	  for (short_type f = 0; f < ls.pgt->structure()->nb_faces(); ++f)
	    mref[indcomp(f+1, 0, 0)](mi) = ppi->int_poly_on_face(P, f);
	}
	// cout << "mref[0] = " << mref[0] << endl;
	// cout << "mref[indcomp(1, 0, 0)] = " << mref[indcomp(1, 0, 0)]<< endl;
      }
      else
      { // very inefficient ...
	scalar_type V;
	pfp.resize(ls.pmt->size());
	it = ls.pmt->begin(), ite = ls.pmt->end();
	for (k = 0; it != ite; ++it, ++k) {
	  pfp[k] = fem_precomp((*it).pfi, &(pai->integration_points()));
	}

	for (;!mi.finished(sizes);mi.incrementation(sizes)) {
	  for (short_type hi = 0; hi < nhess; ++hi) {
	    for (size_type ip = 0; ip < nint; ++ip) {
	      V = 1.0; mref_coeff[indcomp(ip, hi, 0)] = pai->coeff(ip);
	      it = ls.pmt->begin(), ite = ls.pmt->end();
	      mit = mi.begin();
	      short_type hc = 1;

	      for (k = 0; it != ite; ++it, ++k) { 
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
		  if (hi & hc) {
		    V *= -(pfp[k]->grad(ip))[ind + (*it).pfi->nb_base() *
			   (*it).pfi->target_dim() *
			   (*mit % ls.pgt->structure()->dim())];
		    ++mit;
		  }
		  else {
		    V *= (pfp[k]->hess(ip))[ind
			   + (*it).pfi->nb_base() * (*it).pfi->target_dim()
			   * (*mit)];
		    ++mit;
		  }
		  hc *= 2; break;
		}
	      }
	      mref[indcomp(ip, hi, 0)](mi) = V;
	    }
	  }
	  
	}
	// If the geometric transformation is linear, it is possible to
	// precompute the integrals
	if (pgt->is_linear())
	{
	  mref[indcomp(0, 0, 0)] *= pai->coeff(0);
	  size_type j = 0;
	  for (size_type i = 1; i < pai->nb_points_on_convex(); ++i)
	    mref[indcomp(0, 0, 0)].addmul(pai->coeff(i),
					  mref[indcomp(++j, 0, 0)]);
	  for (short_type f = 0; f < pgt->structure()->nb_faces(); ++f) {
	    mref[indcomp(f+1, 0, 0)] = mref[indcomp(++j, 0, 0)];
	    mref[indcomp(f+1, 0, 0)] *= pai->coeff_on_face(f, 0);
	    for (size_type i = 1; i < pai->nb_points_on_face(f); ++i)
	      mref[indcomp(f+1, 0, 0)].addmul(pai->coeff_on_face(f, i),
					      mref[indcomp(++j, 0, 0)]);
	  }
	  nint = nbf + 1;
	  mref.resize(nint * 3);
	  mref_count.resize(nint * 3);
	  mref_coeff.resize(nint * 3);

	  // cout << "mref[0] = " << mref[0] << endl;
	}
      }

      indv.resize(nbf+1);
      if (pgt->is_linear()) {
	indv[0].resize(1); indv[0][0] = 0; mref_coeff[indcomp(0, 0, 0)] = 1.0;
	for (short_type f = 0; f < nbf; ++f) {
	  indv[f+1].resize(1); indv[f+1][0] = f+1;
	  mref_coeff[indcomp(f+1, 0, 0)] = 1.0;
	}
      }
      else {
	indv[0].resize(pai->nb_points_on_convex());
	size_type j = 0;
	for (size_type i = 0; i < indv[0].size(); ++i) indv[0][i] = j++;
	for (short_type f = 0; f < nbf; ++f) {
	  indv[f+1].resize(pai->nb_points_on_face(f));
	  for (size_type i = 0; i < indv[f+1].size(); ++i) indv[f+1][i] = j++;
	}

      }

    }

    inline void compute(base_tensor &t, const base_matrix &G, size_type ir)
    {
      dim_type P = pgt->structure()->dim(), N = G.nrows();
      short_type NP = pgt->nb_points();
      scalar_type J;
      if (G.ncols() != NP) DAL_THROW(dimension_error, "dimensions mismatch");
      
      K.resize(N, P); CS.resize(P, P); TMP1.resize(P, P); B.resize(N, P);
      if (hess_reduction.size() > 0) {
	B2.resize(P*P, P*P); B3.resize(N*N, P*P); Htau.resize(N, P*P);
	B32.resize(N*N, P*P); B2.fill(0.0);
      }
      std::fill(mref_count.begin(), mref_count.end(), 0);
      if (ir > 0) {
	un.resize(P); up.resize(N);
	un = pgt->normals()[ir-1];
      }

      for (size_type ivt = 0; ivt < indv[ir].size(); ++ivt) {

	size_type ip = indv[ir][ivt];
	// on peut simplifier les calculs pour N = P
	// cout << "mat G : " << G << endl;
	// cout << "mat grad : " << pgp->grad(ip) << endl;
	if (pgt->is_linear())
	  bgeot::mat_product(G, pgp->grad(0), K);
	else
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
	  if (!(pgt->is_linear())) {
	    bgeot::mat_product(G, pgp->hessian(ip), Htau);
	    for (short_type i = 0; i < P; ++i)
	      for (short_type j = 0; j < P; ++j)
		for (short_type k = 0; k < P; ++k)
		  for (short_type l = 0; l < N; ++l)
		    B2(i + P*j, k) += Htau(l, i + P*j) * B(l,k);
	    bgeot::mat_product(B3, B2, B32);
	  }
	}

	for (short_type hi = 0; hi < nhess; ++hi) {
	  short_type *p = &(mref_count[indcomp(ip, hi, 0)]);
	  short_type j, k;

	  if (grad_reduction.size() > 0) {
	    std::deque<short_type>::const_iterator it = grad_reduction.begin(),
	      ite = grad_reduction.end();
	    for ( ; it != ite; ++it) {
	      j = *p; k = (j < 2) ? j + 1 : 1; *p = k;
	      mref[indcomp(ip, hi, k)].mat_transp_reduction
		(mref[indcomp(ip, hi, j)], B, *it);
	    }
	  }

	  if (hess_reduction.size() > 0) { // une mise en facteur peut être
	    // faite ... mais penible à réaliser.
	    std::deque<short_type>::const_iterator it = hess_reduction.begin(),
	      ite = hess_reduction.end();
	    for (short_type l = 1; it != ite; ++it, l *= 2) {
	      j = *p;  k = (j < 2) ? j + 1 : 1; *p = k;
	      if (hi & l)
		mref[indcomp(ip, hi, k)].mat_transp_reduction
		  (mref[indcomp(ip,hi,j)], B32, *it);
	      else
	        mref[indcomp(ip, hi, k)].mat_transp_reduction
		  (mref[indcomp(ip,hi,j)], B3, *it);
	    }
	  }
	  
	  if (*p == 0)
	    { mref[indcomp(ip, hi, 1)] = mref[indcomp(ip, hi, 0)]; *p = 1; }
	  mref[indcomp(ip, hi, *p)] *= J;

	}
      }

      bool ini = false;
      short_type j = mref_count[indcomp(indv[ir][0], 0, 0)];
      for (size_type ivt = 0; ivt < indv[ir].size(); ++ivt) {
	for (short_type hi = 0; hi < nhess; ++hi)
	  if (ini)
	    t.addmul(mref_coeff[indcomp(indv[ir][ivt], hi, 0)],
		     mref[indcomp(indv[ir][ivt], hi, j)]);
	  else {
	    t = mref[indcomp(indv[ir][ivt], hi, j)];
	    t *= mref_coeff[indcomp(indv[ir][ivt], hi, 0)];
	    ini = true;
	  }
      }

      /* Applying linear transformation for non tau-equivalent elements.   */

      if (trans_reduction.size() > 0) {
	std::deque<short_type>::const_iterator it = trans_reduction.begin(),
	  ite = trans_reduction.end();
	std::deque<pfem>::const_iterator iti = trans_reduction_pfi.begin();
	for ( ; it != ite; ++it, ++iti) { 
	  if ((*iti)->nb_dof() != M.nrows() || (*iti)->nb_base() != M.ncols())
	    M.resize((*iti)->nb_base(), (*iti)->nb_dof());
	  (*iti)->mat_trans(M, G, pgt);
	  base_tensor aux = t; // Optimisable (avoid copy).
	  t.mat_reduction(aux, M, *it);
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

