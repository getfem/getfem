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
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */


#include <deque>
#include <getfem_mat_elem.h>

namespace getfem
{

  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  struct _pre_geot_light
  {
    bgeot::pgeometric_trans pgt;
    pintegration_method ppi;
    bool operator < (const _pre_geot_light &ls) const
    {
      if (pgt < ls.pgt) return true; if (pgt > ls.pgt) return false; 
      if (ppi.method.ppi < ls.ppi.method.ppi) return true; return false;
    }
    _pre_geot_light(bgeot::pgeometric_trans pg, pintegration_method pi)
    { pgt = pg; ppi = pi; }
    _pre_geot_light(void) { }
   
  };

  class _geotrans_precomp
  {
    protected :

      bgeot::pgeometric_trans pgt;
      pintegration_method ppi;
      std::vector<base_matrix> pc;
      std::vector<base_matrix> hpc;

    public :

      inline const base_matrix &grad(size_type i) { return pc[i]; }
      inline const base_matrix &hessian(size_type i) { return hpc[i]; }

      _geotrans_precomp(const _pre_geot_light &ls)
      {
	if (ls.ppi.is_ppi) {
	  base_poly P;
	  assert(ls.pgt->basic_structure() == ls.ppi.method.ppi->structure());
	  assert(ls.pgt->is_linear());
	  pgt = ls.pgt; ppi = ls.ppi;
	  pc.resize(1);
	  hpc.resize(1);
	  pc[0] = base_matrix(pgt->nb_points() , pgt->structure()->dim());
	  hpc[0] = base_matrix(pgt->nb_points(),
			       dal::sqr(pgt->structure()->dim()));
	  hpc[0].fill(0.0);
	  for (size_type i = 0; i < pgt->nb_points(); ++i)
	    for (dim_type n = 0; n < pgt->structure()->dim(); ++n)
	      { P = pgt->poly_vector()[i]; P.derivative(n); pc[0](i,n)=P[0]; }
	}
	else {
	  base_poly P, Q;
	  dim_type N = ls.pgt->structure()->dim();
	  assert(ls.pgt->basic_structure() == ls.ppi.method.pai->structure());
	  pgt = ls.pgt; ppi = ls.ppi;
	  pc.resize(ppi.method.pai->nb_points());
	  hpc.resize(ppi.method.pai->nb_points());
	  std::fill(pc.begin(), pc.end(),
		    base_matrix(pgt->nb_points() , N));
	  std::fill(hpc.begin(), hpc.end(),
		    base_matrix(dal::sqr(N), pgt->nb_points()));
	  for (size_type i = 0; i < pgt->nb_points(); ++i)
	    for (dim_type n = 0; n < N; ++n) {
	      P = pgt->poly_vector()[i];
	      P.derivative(n);
	      for (size_type j = 0; j < ppi.method.pai->nb_points(); ++j)
		pc[j](i,n) = P.eval(ppi.method.pai->point(j).begin());
	      for (dim_type m = 0; m <= n; ++m) {
		Q = P; Q.derivative(m);
		for (size_type j = 0; j < ppi.method.pai->nb_points(); ++j)
		  hpc[j](m * N + n, i) = hpc[j](n * N + m, i)
		    = P.eval(ppi.method.pai->point(j).begin());
	      }
	    }
	}
      }
  };
  
  typedef _geotrans_precomp * pgeotrans_precomp;

  pgeotrans_precomp geotrans_precomp(bgeot::pgeometric_trans pg,
				     pintegration_method pi)
  { 
    static dal::FONC_TABLE<_pre_geot_light, _geotrans_precomp> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_pre_geot_light, _geotrans_precomp>();
      isinit = true;
    }
    return tab->add(_pre_geot_light(pg, pi));
  }

  /* ********************************************************************* */
  /*       Precomputation on fem.                                          */
  /* ********************************************************************* */

  struct _pre_fem_light
  {
    pfem pf;
    pintegration_method ppi;
    bool operator < (const _pre_fem_light &ls) const
    {
      if (pf < ls.pf) return true; if (pf > ls.pf) return false; 
      if (ppi.method.ppi < ls.ppi.method.ppi) return true; return false;
    }
    _pre_fem_light(pfem pff, pintegration_method pi)
    { pf = pff; ppi = pi; }
    _pre_fem_light(void) { }
   
  };

  class _fem_precomp
  {
    protected :
      
      std::vector<base_tensor> c;
      std::vector<base_tensor> pc;
      std::vector<base_tensor> hpc;

    public :

      inline const base_tensor &val(size_type i) { return c[i]; }
      inline const base_tensor &grad(size_type i) { return pc[i]; }
      inline const base_tensor &hess(size_type i) { return hpc[i]; }

      _fem_precomp(const _pre_fem_light &ls) {
	assert(!(ls.ppi.is_ppi));
	assert(ls.pf->basic_structure() == ls.ppi.method.pai->structure());
	pc.resize(ls.ppi.method.pai->nb_points());
	hpc.resize(ls.ppi.method.pai->nb_points());
	c.resize(ls.ppi.method.pai->nb_points());
	for (size_type i = 0; i < ls.ppi.method.pai->nb_points(); ++i) {
	  ls.pf->base_value(ls.ppi.method.pai->point(i), c[i]);
	  ls.pf->grad_base_value(ls.ppi.method.pai->point(i), pc[i]);
	  ls.pf->hess_base_value(ls.ppi.method.pai->point(i), hpc[i]);
	}
      }
  };
  
  typedef _fem_precomp * pfem_precomp;

  pfem_precomp fem_precomp(pfem pf, pintegration_method pi)
  { 
    static dal::FONC_TABLE<_pre_fem_light, _fem_precomp> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_pre_fem_light, _fem_precomp>();
      isinit = true;
    }
    return tab->add(_pre_fem_light(pf, pi));
  }


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
      if (ppi.method.ppi < ls.ppi.method.ppi) return true;
      if (ppi.method.ppi > ls.ppi.method.ppi) return false; 
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
    bgeot::ppoly_integration ppi;
    bgeot::papprox_integration pai;
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


    _emelem_comp_structure(const _emelem_comp_light &ls)
    { // optimisable ... !!
      pgt = ls.pgt;
      pgp = geotrans_precomp(ls.pgt, ls.ppi);
      pme = ls.pmt;
      ppi = ls.ppi.method.ppi;
      pai = ls.ppi.method.pai;
      is_ppi = ls.ppi.is_ppi;
      size_type k;

      // cout << "debut de construction\n";
      
      nbf = pgt->structure()->nb_faces();
      nhess = 1;
      nint = is_ppi ? nbf + 1 : pai->nb_points();

      mat_elem_type::const_iterator it = ls.pmt->begin(), ite = ls.pmt->end();
      
      for (k = 0; it != ite; ++it, ++k) {
	if (is_ppi) assert((*it).pfi->is_polynomial());
	assert((*it).pfi->basic_structure() == pgt->basic_structure());
	if (!((*it).pfi->is_equivalent())) {
	  trans_reduction.push_back(k);
	  trans_reduction_pfi.push_back((*ite).pfi);
	}
	switch ((*it).t) {
	case GETFEM__BASE    : break;
	case GETFEM__GRAD    : ++k; grad_reduction.push_back(k); break;
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
	STD_NEEDED cout << "Warning, very large elementary computations,\n"
	     << "Are you sure to want to compute this elementary matrix ?\n"
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
	      { ind += (*it).pfi->nb_dof() * (*mit); ++mit; }
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
	// cout << "mref[indcomp(1, 0, 0)] = " << mref[indcomp(1, 0, 0)] << endl;
      }
      else
      { // very inefficient ...
	scalar_type V;
	pfp.resize(ls.pmt->size());
	it = ls.pmt->begin(), ite = ls.pmt->end();
	for (k = 0; it != ite; ++it, ++k)
	  pfp[k] = fem_precomp((*it).pfi, pai);

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
		  { ind += (*it).pfi->nb_dof() * (*mit); ++mit; }
		switch ((*it).t) {
		case GETFEM__BASE    :
		  V *= (pfp[k]->val(ip))[ind]; break;
		case GETFEM__GRAD    :
		  V *= (pfp[k]->grad(ip))[ind + (*it).pfi->nb_dof() *
			 (*it).pfi->target_dim() * (*mit)];
		  ++mit; break;
		case GETFEM__HESSIAN :
		  if (hi & hc) {
		    V *= -(pfp[k]->grad(ip))[ind + (*it).pfi->nb_dof() *
			   (*it).pfi->target_dim() *
			   (*mit % ls.pgt->structure()->dim())];
		    ++mit;
		  }
		  else {
		    V *= (pfp[k]->hess(ip))[ind
			   + (*it).pfi->nb_dof() * (*it).pfi->target_dim()
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
      assert(G.ncols() == NP);
      K.resize(N, P); CS.resize(P, P); TMP1.resize(P, P); B.resize(N, P);
      if (hess_reduction.size() > 0) {
	B2.resize(P*P, P*P); B3.resize(N*N, P*P); Htau.resize(N, P*P);
	B32.resize(N*N, P*P); B2.fill(0.0);
      }
      std::fill(mref_count.begin(), mref_count.end(), 0);
      if (ir > 0)
      {
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

	if (ir > 0)
	{
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

	  if (grad_reduction.size() > 0)
	  {
	    std::deque<short_type>::const_iterator it = grad_reduction.begin(),
	      ite = grad_reduction.end();
	    for ( ; it != ite; ++it) {
	      j = *p; k = (j < 2) ? j + 1 : 1; *p = k;
	      mref[indcomp(ip, hi, k)].mat_reduction(mref[indcomp(ip, hi, j)],
						     B, *it);
	    }
	  }

	  if (hess_reduction.size() > 0) { // une mise en facteur peut être
	    // faite ... mais penible à réaliser.
	    std::deque<short_type>::const_iterator it = hess_reduction.begin(),
	      ite = hess_reduction.end();
	    for (short_type l = 1; it != ite; ++it, l *= 2) {
	      j = *p;  k = (j < 2) ? j + 1 : 1; *p = k;
	      if (hi & l)
		mref[indcomp(ip, hi, k)].mat_reduction(mref[indcomp(ip,hi,j)],
						       B32, *it);
	      else
	        mref[indcomp(ip, hi, k)].mat_reduction(mref[indcomp(ip,hi,j)],
						       B3, *it);
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
	    t.addmul(mref_coeff[indcomp(indv[ir][ivt], hi, j)],
		     mref[indcomp(indv[ir][ivt], hi, j)]);
	  else {
	    // cout << "here we put " << indcomp(indv[ir][ivt], hi, j) << " : " << j << "tensor : " << mref[indcomp(indv[ir][ivt], hi, j)] << endl;
	    t = mref[indcomp(indv[ir][ivt], hi, j)];
	    t *= mref_coeff[indcomp(indv[ir][ivt], hi, 0)];
	    ini = true;
	  }
      }

      /* Applying linear transformation for non tau-equivalent elements.   */

      if (trans_reduction.size() > 0)
      {
	std::deque<short_type>::const_iterator it = trans_reduction.begin(),
	  ite = trans_reduction.end();
	std::deque<pfem>::const_iterator
	  iti = trans_reduction_pfi.begin();
	for ( ; it != ite; ++it, ++iti)
	{ 
	  if (t.size(*it) != M.nrows() || t.size(*it) != M.ncols())
	    M.resize(t.size(*it), t.size(*it));
	  (*iti)->mat_trans(M, G);
	  t.mat_reduction(t, M, *it);
	}
      }
      
    }


    void compute(base_tensor &t, const base_matrix &G)
    { compute(t, G, 0); }

    void compute_on_face(base_tensor &t, const base_matrix &G, short_type f)
    { compute(t, G, f+1); }

  };

  pmat_elem_computation mat_elem(pmat_elem_type pm, pintegration_method pi,
				 bgeot::pgeometric_trans pg)
  { 
    static dal::FONC_TABLE<_emelem_comp_light, _emelem_comp_structure> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_emelem_comp_light, _emelem_comp_structure>();
      isinit = true;
    }
    return tab->add(_emelem_comp_light(pm, pi, pg));
  }


  /* norm of returned vector is the ratio between the face surface on
     the reel element and the face surface on the reference element 
     IT IS NOT UNITARY

     pt is the position of the evaluation point on the reference element
  */
  base_vector compute_normal(const base_matrix &G, size_type ir, bgeot::pgeometric_trans pgt, const base_node &pt)
    {
      dim_type P = pgt->structure()->dim(), N = G.nrows();
      short_type NP = pgt->nb_points();
      base_matrix K(N,P), CS(P,P), B(N,P), Grad(pgt->nb_points(),P), TMP1(P,P);
      base_vector un, up;
      base_poly Poly;

      assert(G.ncols() == NP);

      un.resize(P); up.resize(N);
      un = pgt->normals()[ir];
      //cout << "un=" << un << endl;

      for (size_type i = 0; i < pgt->nb_points(); ++i)
	{
	  for (dim_type n = 0; n < N; ++n) {
	    Poly = pgt->poly_vector()[i];
	    Poly.derivative(n);
	    Grad(i,n) = Poly.eval(pt.begin());
	  }
	}


      // on peut simplifier les calculs pour N = P
      // cout << "mat G : " << G << endl;
      // cout << "mat grad : " << Grad << endl;
      bgeot::mat_product(G, Grad, K);
      bgeot::mat_product_tn(K, K, CS);
      bgeot::mat_inv_cholesky(CS, TMP1);
      bgeot::mat_product(K, CS, B);
      bgeot::mat_vect_product(B, un, up);

      return up;
    }

}  /* end of namespace getfem.                                            */

