// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mat_elem.cc : computation of elementary matrices.
//           
// Date    : December 21, 2000.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Yves Renard
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


#include <deque>
#include <dal_singleton.h>
#include <getfem_fem.h>
#include <getfem_mat_elem.h>

extern "C" void daxpy_(const int *n, const double *alpha, const double *x,
		       const int *incx, double *y, const int *incy);
extern "C" void dger_(const int *m, const int *n, const double *alpha,
		      const double *x, const int *incx, const double *y,
		      const int *incy, double *A, const int *lda);

namespace getfem {
  /* ********************************************************************* */
  /*       Elementary matrices computation.                                */
  /* ********************************************************************* */

  struct emelem_comp_key_  : virtual public dal::static_stored_object_key {
    pmat_elem_type pmt;
    pintegration_method ppi;
    bgeot::pgeometric_trans pgt;
    /* prefer_comp_on_real_element: compute elementary matrices on the real
       element if possible (i.e. if no exact integration is used); this allow
       using inline reduction during the integration */
    bool prefer_comp_on_real_element; 
    virtual bool compare(const static_stored_object_key &oo) const {
      const emelem_comp_key_ &o = dynamic_cast<const emelem_comp_key_ &>(oo);  
      if (pmt < o.pmt) return true; if (o.pmt < pmt) return false; 
      if (ppi < o.ppi) return true; if (o.ppi < ppi) return false; 
      if (pgt < o.pgt) return true; if (o.pgt < pgt) return false;
      if (prefer_comp_on_real_element < o.prefer_comp_on_real_element)
	return true;
      return false;
    }
    emelem_comp_key_(pmat_elem_type pm, pintegration_method pi,
		       bgeot::pgeometric_trans pg, bool on_relt)
    { pmt = pm; ppi = pi; pgt = pg; prefer_comp_on_real_element = on_relt; }
    emelem_comp_key_(void) { }
  };
  
  struct emelem_comp_structure_ : public mat_elem_computation {
    bgeot::pgeotrans_precomp pgp;
    ppoly_integration ppi;
    papprox_integration pai;
    bool is_ppi;
    mutable std::vector<base_tensor> mref;
    mutable std::vector<pfem_precomp> pfp;
    mutable std::vector<base_tensor> elmt_stored;
    short_type nbf, dim; 
    std::deque<short_type> grad_reduction, hess_reduction, trans_reduction;
    std::deque<pfem> trans_reduction_pfi;
    mutable base_small_vector un, up;
    mutable bool faces_computed;
    mutable bool volume_computed;
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

    emelem_comp_structure_(pmat_elem_type pm, pintegration_method pi,
			   bgeot::pgeometric_trans pg, 
			   bool prefer_comp_on_real_element) {
      
      pgt = pg;
      pgp = bgeot::geotrans_precomp(pg, &(pi->integration_points()));
      pme = pm;
      switch (pi->type()) {
      case IM_EXACT: 
	ppi = pi->exact_method(); pai = 0;  is_ppi = true; break;
      case IM_APPROX: 
	ppi = 0; pai = pi->approx_method(); is_ppi = false; break;
      case IM_NONE: 
	DAL_THROW(dal::failure_error, 
		  "Attempt to use IM_NONE integration method in assembly!\n");
      }
      faces_computed = volume_computed = false;
      is_linear = pgt->is_linear();
      computed_on_real_element = !is_linear || (prefer_comp_on_real_element && !is_ppi);
      nbf = pgt->structure()->nb_faces();
      dim = pgt->structure()->dim();
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
      
      for (size_type k = 0; it != ite; ++it, ++k) {
	if ((*it).pfi) {
	  if ((*it).pfi->is_on_real_element()) computed_on_real_element = true;
	  if (is_ppi && (!((*it).pfi->is_polynomial()) || !is_linear 
			 || computed_on_real_element))
	    DAL_THROW(std::invalid_argument, 
		      "Exact integration not allowed in this context");
// 	  if((*it).pfi->basic_structure() != pgt->basic_structure())
// 	    DAL_THROW(std::invalid_argument, "incorrect computation");
	  
	  if (!((*it).pfi->is_equivalent()) && (*it).t != GETFEM_NONLINEAR_) {
	    // TODO : le numero d'indice à reduire peut changer ...
	    trans_reduction.push_back(k);
	    trans_reduction_pfi.push_back((*it).pfi);
	  }
	}
	switch ((*it).t) {
	case GETFEM_BASE_    : break;
	case GETFEM_UNIT_NORMAL_    : computed_on_real_element = true; break;
	case GETFEM_GRAD_    : ++k;
	  if (!((*it).pfi->is_on_real_element())) grad_reduction.push_back(k);
	  break;
	case GETFEM_HESSIAN_ : ++k; hess_reduction.push_back(k); break;
	case GETFEM_NONLINEAR_ :
	  if ((*it).nl_part == 0) {
	    for (dim_type ii = 1; ii < (*it).nlt->sizes().size(); ++ii) ++k;
	    if (is_ppi) DAL_THROW(failure_error, "For nonlinear terms you have "
				  "to use approximated integration");
	    computed_on_real_element = true;
	  }
	  break;
	}
      }

      if (!is_ppi) {
	pfp.resize(pme->size());
	it = pme->begin(), ite = pme->end();
	for (size_type k = 0; it != ite; ++it, ++k)
	  if ((*it).pfi)
	    pfp[k] = fem_precomp((*it).pfi, &(pai->integration_points()));
	  else pfp[k] = 0;
	elmt_stored.resize(pme->size());
      }
      if (!computed_on_real_element) mref.resize(nbf + 1);
    }

    void add_elem(base_tensor &t, fem_interpolation_context& ctx,
                  scalar_type J, bool first, bool trans, 
                  mat_elem_integration_callback *icb,
		  bgeot::multi_index sizes) const {
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();

      bgeot::multi_index::iterator mit = sizes.begin();
      for (size_type k = 0; it != ite; ++it, ++k) {
	if (pfp[k]) ctx.set_pfp(pfp[k]);
	++mit; if ((*it).pfi && (*it).pfi->target_dim() > 1) ++mit;
	
	switch ((*it).t) {
	case GETFEM_BASE_    :
	  (*it).pfi->real_base_value(ctx, elmt_stored[k]);
	  break;
	case GETFEM_GRAD_    :
	  if (trans) {
	    (*it).pfi->real_grad_base_value(ctx, elmt_stored[k]);
	    *mit++ = ctx.N();
	  }
	  else
	    elmt_stored[k] = pfp[k]->grad(ctx.ii());
	  break;
	case GETFEM_HESSIAN_ :
	  if (trans) {
	    (*it).pfi->real_hess_base_value(ctx, elmt_stored[k]);
	    *mit++ = gmm::sqr(ctx.N());
	  }
	  else {
	    base_tensor tt = pfp[k]->hess(ctx.ii());
	    bgeot::multi_index mim(3);
	    mim[2] = gmm::sqr(tt.sizes()[2]); mim[1] = tt.sizes()[1];
	    mim[0] = tt.sizes()[0];
	    tt.adjust_sizes(mim);
	    elmt_stored[k] = tt;
	  }
	  break;
	case GETFEM_UNIT_NORMAL_ :
	  *(mit-1) = ctx.N();
	  { 
	    bgeot::multi_index sz(1); sz[0] = ctx.N();
	    elmt_stored[k].adjust_sizes(sz);
	  }
	  std::copy(up.begin(), up.end(), elmt_stored[k].begin());
	  break;
	case GETFEM_NONLINEAR_ :
	  if ((*it).nl_part != 0) { /* for auxiliary fem of the nonlinear_term, */
	                            /* the "prepare" method is called           */
	     (*it).nlt->prepare(ctx, (*it).nl_part);
	    /* the dummy assistant multiplies everybody by 1
	       -> not efficient ! */
	    bgeot::multi_index sz(1); sz[0] = 1;
	    elmt_stored[k].adjust_sizes(sz); elmt_stored[k][0] = 1.;
	  } else {
	    elmt_stored[k].adjust_sizes((*it).nlt->sizes());
	    (*it).nlt->compute(ctx, elmt_stored[k]);
	    for (dim_type ii = 1; ii < (*it).nlt->sizes().size(); ++ii) ++mit;
	  }
	  break;
	}
      }
      
      //expand_product_old(t,J*pai->coeff(ctx.ii()), first);
      scalar_type c = J*pai->coeff(ctx.ii());
      if (!icb) {
	if (first) { t.adjust_sizes(sizes); }
	expand_product_daxpy(t, c, first);
      } else {
        icb->eltm.resize(0);
	for (unsigned k=0; k != pme->size(); ++k) {
	  if (icb && !((*pme)[k].t == GETFEM_NONLINEAR_ && (*pme)[k].nl_part != 0))
	    icb->eltm.push_back(&elmt_stored[k]);
	}
	icb->exec(t, first, c);
      }
    }


    void expand_product_old(base_tensor &t, scalar_type J, bool first) const {
      scalar_type V;
      size_type k;
      if (first) std::fill(t.begin(), t.end(), 0.0);
      base_tensor::iterator pt = t.begin();
      std::vector<base_tensor::const_iterator> pts(pme->size());
      std::vector<scalar_type> Vtab(pme->size());
      for (k = 0; k < pme->size(); ++k)
	pts[k] = elmt_stored[k].begin();
      
      size_type k0 = 0;
      unsigned n0 = elmt_stored[0].size();
      /*while (elmt_stored[k0].size() == 1 && k0+1 < pme->size()) {
        J *= elmt_stored[k0][0];
        ++k0; n0 = elmt_stored[k0].size();
        }*/
      base_tensor::const_iterator pts0 = pts[k0];


      k = pme->size()-1; Vtab[k] = J;
      /* very heavy expansion .. takes much time */
      do {
        for (V = Vtab[k]; k!=k0; --k)
          Vtab[k-1] = V = *pts[k] * V;
        for (k=0; k < n0; ++k)
          *pt++ += V * pts0[k];
        for (k=k0+1; k != pme->size() && ++pts[k] == elmt_stored[k].end(); ++k)
          pts[k] = elmt_stored[k].begin();
      } while (k != pme->size());
      if (pt != t.end()) DAL_THROW(internal_error, "Internal error");
    }

    /* do the tensorial product using the blas function daxpy (much more
       efficient than a loop).

       efficiency is maximized when the first tensor has a large dimension
     */
    void expand_product_daxpy(base_tensor &t, scalar_type J, bool first)const {
      size_type k;
      base_tensor::iterator pt = t.begin();
      static std::vector<base_tensor::const_iterator> pts, es_beg, es_end;
      static std::vector<scalar_type> Vtab;
      pts.resize(pme->size()); es_beg.resize(pme->size());
      es_end.resize(pme->size()); Vtab.resize(pme->size());
      size_type nm = 0;
      if (first) memset(&(*t.begin()), 0, t.size()*sizeof(*t.begin())); //std::fill(t.begin(), t.end(), 0.0);
      for (k = 0, nm = 0; k < pme->size(); ++k) {
        if (elmt_stored[k].size() != 1) {
          es_beg[nm] = elmt_stored[k].begin();
          es_end[nm] = elmt_stored[k].end();
          pts[nm] = elmt_stored[k].begin(); 
          ++nm;
        } else J *= elmt_stored[k][0];
      }
      if (nm == 0) {
        t[0] += J;
      } else {
        int n0 = es_end[0] - es_beg[0];
        base_tensor::const_iterator pts0 = pts[0];

        /* very heavy reduction .. takes much time */
        k = nm-1; Vtab[k] = J;
        int one = 1;
        scalar_type V;
        do {
          for (V = Vtab[k]; k; --k)
            Vtab[k-1] = V = *pts[k] * V;
          daxpy_(&n0, &V, const_cast<double*>(&(pts0[0])), &one,
		 (double*)&(*pt), &one); 
          pt+=n0;
          for (k=1; k != nm && ++pts[k] == es_end[k]; ++k)
            pts[k] = es_beg[k];
        } while (k != nm);
        if (pt != t.end()) DAL_THROW(internal_error, "Internal error");
      }
    }


    void pre_tensors_for_linear_trans(bool volumic) const {

      if ((volumic && volume_computed) || (!volumic && faces_computed)) return;
      // scalar_type exectime = ftool::uclock_sec();

      bgeot::multi_index sizes = pme->sizes(0), mi(sizes.size());
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

	  if ((*it).pfi) {
	    if ((*it).pfi->target_dim() > 1)
	      { ind += (*it).pfi->nb_base(0) * (*mit); ++mit; }
	    
	    Q = ((ppolyfem)((*it).pfi).get())->base()[ind];
	  }

	  switch ((*it).t) {
	  case GETFEM_GRAD_    : Q.derivative(*mit); ++mit; break;
	  case GETFEM_HESSIAN_ :
	    Q.derivative(*mit % dim); Q.derivative(*mit / dim);
	    ++mit; break;
	  case GETFEM_BASE_ : break;
	  case GETFEM_UNIT_NORMAL_ :
	  case GETFEM_NONLINEAR_ :
	    DAL_THROW(failure_error, "No nonlinear term allowed here");
	  }
	  ++it;

	  if (it != ite && *mit != old_ind) {
	    old_ind = *mit; 
	    P.one();
	    for (; it != ite; ++it) {
	      ind = *mit; ++mit;
	      
	      if ((*it).pfi->target_dim() > 1)
		{ ind += (*it).pfi->nb_base(0) * (*mit); ++mit; }
	      R = ((ppolyfem)((*it).pfi).get())->base()[ind];
	      
	      switch ((*it).t) {
	      case GETFEM_GRAD_    : R.derivative(*mit); ++mit; break;
	      case GETFEM_HESSIAN_ :
		R.derivative(*mit % dim); R.derivative(*mit / dim);
		++mit; break;
	      case GETFEM_BASE_ : break;
	      case GETFEM_UNIT_NORMAL_ :
	      case GETFEM_NONLINEAR_ :
		DAL_THROW(failure_error, "No nonlinear term allowed here");
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
	fem_interpolation_context ctx;
	size_type ind_l = 0, nb_ptc = pai->nb_points_on_convex(), 
	  nb_pt_l = nb_ptc, nb_pt_tot =(volumic ? nb_ptc : pai->nb_points());
	for (size_type ip = (volumic ? 0:nb_ptc); ip < nb_pt_tot; ++ip) {
	  while (ip == nb_pt_l && ind_l < nbf)
	    { nb_pt_l += pai->nb_points_on_face(ind_l); ind_l++; }
	  ctx.set_ii(ip); 
	  add_elem(mref[ind_l], ctx, 1.0, first, false, NULL, sizes);
	  first = false;
	}
      }
      // cout << "precompute Mat elem computation time : "
      //   << ftool::uclock_sec() - exectime << endl;
    }


    void compute(base_tensor &t, const base_matrix &G, size_type ir,
		 size_type elt, mat_elem_integration_callback *icb = 0) const {
      dim_type P = dim, N = G.nrows();
      short_type NP = pgt->nb_points();
      fem_interpolation_context ctx(pgp,0,0,G,elt);
      bgeot::multi_index sizes = pme->sizes(elt);

      if (G.ncols() != NP) DAL_THROW(dimension_error, "dimensions mismatch");
      if (ir > 0) {
	up.resize(N); un.resize(P);
	un = pgt->normals()[ir-1];
      }
      base_tensor taux;
      bool flag = false;

      if (!computed_on_real_element) {
	pre_tensors_for_linear_trans(ir == 0);
	const base_matrix& B = ctx.B(); // compute B and J
	scalar_type J=ctx.J();
	if (ir > 0) {
	  gmm::mult(B, un, up);
	  scalar_type nup = bgeot::vect_norm2(up);
	  J *= nup; up /= nup;
	}
     
	t = mref[ir]; t *= J;
	
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
	    (flag ? t:taux).mat_transp_reduction(flag ? taux:t, ctx.B3(), *it);
	    flag = !flag;
	  }
	}
	
      } else { // non linear transformation and methods defined on real elements
	bool first = true;

	for (size_type ip=(ir == 0) ? 0 : pai->repart()[ir-1];
	     ip < pai->repart()[ir]; ++ip, first = false) {
	  ctx.set_ii(ip);
	  const base_matrix& B = ctx.B(); // J computed as side-effect
	  scalar_type J = ctx.J();
	  if (ir > 0) {
	    gmm::mult(B, un, up);
	    scalar_type nup = bgeot::vect_norm2(up);
	    J *= nup; up /= nup;
	  }	  
	  add_elem(t, ctx, J, first, true, icb, sizes);
	}
      }

      /* Applying linear transformation for non tau-equivalent elements.   */
      
      if (trans_reduction.size() > 0) {
	if (icb) // Dans ce cas, il faudrait annuler la reduction finale (si
	  // l'indice des numerod de fonctions de base est réduit) et faire
	  // la reduction sur chaque point de Gauss.
	  DAL_INTERNAL_ERROR("Non tau-equivalent elements are not"
			     "working with this kind of assembly!");
	std::deque<short_type>::const_iterator it = trans_reduction.begin(),
	  ite = trans_reduction.end();
	std::deque<pfem>::const_iterator iti = trans_reduction_pfi.begin();
	for ( ; it != ite; ++it, ++iti) { 
	  ctx.set_pf(*iti);
	  (flag ? t:taux).mat_reduction(flag ? taux:t, ctx.M(), *it);
	  flag = !flag;
	}
      }
      
      if (flag) t = taux;
    }
    
    void compute(base_tensor &t, const base_matrix &G, size_type elt, 
		 mat_elem_integration_callback *icb) const   
    { compute(t, G, 0, elt, icb); }

    void compute_on_face(base_tensor &t, const base_matrix &G,
			 short_type f, size_type elt, 
			 mat_elem_integration_callback *icb) const
    { compute(t, G, f+1, elt, icb); }
  };

  pmat_elem_computation mat_elem(pmat_elem_type pm, pintegration_method pi,
				 bgeot::pgeometric_trans pg, 
                                 bool prefer_comp_on_real_element) { 
    dal::pstatic_stored_object o
      = dal::search_stored_object(emelem_comp_key_(pm, pi, pg,
						prefer_comp_on_real_element));
    if (o) return dal::stored_cast<mat_elem_computation>(o);
    pmat_elem_computation p = new emelem_comp_structure_(pm, pi, pg,
						prefer_comp_on_real_element);
    dal::add_stored_object(new emelem_comp_key_(pm, pi, pg,
					       prefer_comp_on_real_element),
			   p, pm, pi, pg);
    return p;
  }


}  /* end of namespace getfem.                                            */

