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
#include <dal_singleton.h>
#include <getfem_mat_elem.h>
#include <getfem_precomp.h>
#include <bgeot_precomp.h>

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
    bgeot::pgeotrans_precomp pgp;
    ppoly_integration ppi;
    papprox_integration pai;
    bool is_ppi;
    std::vector<base_tensor> mref;
    std::vector<pfem_precomp> pfp;
    std::vector<base_tensor> elmt_stored;
    short_type nbf, dim; 
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
      pgp = bgeot::geotrans_precomp(ls.pgt, &(ls.ppi->integration_points()));
      pme = ls.pmt;
      switch (ls.ppi->type()) {
      case IM_EXACT: 
	ppi = ls.ppi->exact_method(); pai = 0;  is_ppi = true; break;
      case IM_APPROX: 
	ppi = 0; pai = ls.ppi->approx_method(); is_ppi = false; break;
      case IM_NONE: 
	DAL_THROW(dal::failure_error, 
		  "Attempt to use IM_NONE integration method in assembly!\n");
      }
      faces_computed = volume_computed = false;
      is_linear = pgt->is_linear();
      computed_on_real_element = !is_linear;
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
	  if((*it).pfi->basic_structure() != pgt->basic_structure())
	    DAL_THROW(std::invalid_argument, "incorrect computation");
	  
	  if (!((*it).pfi->is_equivalent())) {
	    trans_reduction.push_back(k);
	    trans_reduction_pfi.push_back((*it).pfi);
	  }
	}
	switch ((*it).t) {
	case GETFEM_BASE_    : break;
	case GETFEM_GRAD_    : ++k;
	  if (!((*it).pfi->is_on_real_element())) grad_reduction.push_back(k);
	  break;
	case GETFEM_HESSIAN_ : ++k; hess_reduction.push_back(k); break;
	case GETFEM_NONLINEAR_ :
	  if ((*it).nl_part == 0) {
	    for (dim_type ii = 1; ii < (*it).nlt->sizes().size(); ++ii) ++k;
	    if (is_ppi)
	      DAL_THROW(failure_error,
			"For nonlinear terms you have to use approximated integration");
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
	elmt_stored.resize(pme->size());
      }
      if (!computed_on_real_element) mref.resize(nbf + 1);
    }

    void add_elem(base_tensor &t, fem_interpolation_context& ctx,
		  scalar_type J, bool first, bool trans = true) {
      mat_elem_type::const_iterator it = pme->begin(), ite = pme->end();
      bgeot::multi_index mi(pme->mi.size()), sizes = pme->mi;
      bgeot::multi_index::iterator mit = sizes.begin();
      for (size_type k = 0; it != ite; ++it, ++k) {
	ctx.set_pfp(pfp[k]);
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
	    *mit++ = dal::sqr(ctx.N());
	  }
	  else {
	    base_tensor tt = pfp[k]->hess(ctx.ii());
	    bgeot::multi_index mim(3);
	    mim[2] = dal::sqr(tt.sizes()[2]); mim[1] = tt.sizes()[1];
	    mim[0] = tt.sizes()[0];
	    tt.adjust_sizes(mim);
	    elmt_stored[k] = tt;
	  }
	  break;
	case GETFEM_NONLINEAR_ :
	  if ((*it).nl_part != 0) { /* for auxiliary fem of the nonlinear_term, the "prepare" method is called */
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
      
      if (first) {
	t.adjust_sizes(sizes);
	std::fill(t.begin(), t.end(), 0.0);
      }

      scalar_type V;
      size_type k;
      base_tensor::iterator pt = t.begin();
      std::vector<base_tensor::const_iterator> pts(pme->size());
      std::vector<scalar_type> Vtab(pme->size());
      J *= pai->coeff(ctx.ii());
      for (k = 0; k < pme->size(); ++k)
	pts[k] = elmt_stored[k].begin();
      base_tensor::const_iterator pts0 = pts[0];
      
      size_type n0 = elmt_stored[0].size();
      k = pme->size()-1; Vtab[k] = J;
      /* very heavy reduction .. takes much time */
      do {
        for (V = Vtab[k]; k; --k)
          Vtab[k-1] = V = *pts[k] * V;
        for (; k < n0; ++k)
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

	  if ((*it).pfi) {
	    if ((*it).pfi->target_dim() > 1)
	      { ind += (*it).pfi->nb_base() * (*mit); ++mit; }
	    
	    Q = ((ppolyfem)((*it).pfi))->base()[ind];
	  }

	  switch ((*it).t) {
	  case GETFEM_GRAD_    : Q.derivative(*mit); ++mit; break;
	  case GETFEM_HESSIAN_ :
	    Q.derivative(*mit % dim); Q.derivative(*mit / dim);
	    ++mit; break;
	  case GETFEM_BASE_ : break;
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
		{ ind += (*it).pfi->nb_base() * (*mit); ++mit; }
	      R = ((ppolyfem)((*it).pfi))->base()[ind];
	      
	      switch ((*it).t) {
	      case GETFEM_GRAD_    : R.derivative(*mit); ++mit; break;
	      case GETFEM_HESSIAN_ :
		R.derivative(*mit % dim); R.derivative(*mit / dim);
		++mit; break;
	      case GETFEM_BASE_ : break;
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
	  add_elem(mref[ind_l], ctx, 1.0, first, false);
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
      fem_interpolation_context ctx(pgp,0,0,G,elt);

      if (G.ncols() != NP) DAL_THROW(dimension_error, "dimensions mismatch");
      
      if (ir > 0) {
	un.resize(P); up.resize(N);
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
	  J *= bgeot::vect_norm2(up);
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
	    J *= bgeot::vect_norm2(up);
	  }	  
	  add_elem(t, ctx, J, first);
	}
      }

      /* Applying linear transformation for non tau-equivalent elements.   */
      
      if (trans_reduction.size() > 0) {
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
    
    void compute(base_tensor &t, const base_matrix &G, size_type elt)
    { compute(t, G, 0, elt); }

    void compute_on_face(base_tensor &t, const base_matrix &G,
			 short_type f, size_type elt)
    { compute(t, G, f+1, elt); }

  };
  
  struct emelem_comp_light_FUNC_TABLE : 
    public dal::FONC_TABLE<emelem_comp_light_, emelem_comp_structure_> { };

  size_type stored_mat_elem_memsize() {
    const emelem_comp_light_FUNC_TABLE & f = 
      dal::singleton<emelem_comp_light_FUNC_TABLE>::const_instance();
    size_type sz = 0;
    for (dal::bv_visitor i(f.index()); !i.finished(); ++i) {
      sz += f.table()[i]->memsize();
    }
    return sz;
  }

  pmat_elem_computation mat_elem(pmat_elem_type pm, pintegration_method pi,
				 bgeot::pgeometric_trans pg) { 
    return dal::singleton<emelem_comp_light_FUNC_TABLE>
      ::instance().add(emelem_comp_light_(pm, pi, pg));
  }

  /* remove all occurences of pm from the emelem_comp_light_FUNC_TABLE */
  void mat_elem_forget_mat_elem_type(pmat_elem_type pm) {
    emelem_comp_light_FUNC_TABLE& f = 
      dal::singleton<emelem_comp_light_FUNC_TABLE>::instance();
    for (dal::bv_visitor_c i(f.index()); !i.finished(); ++i) { 
      if (f.light_table()[i].pmt == pm) f.sup(f.light_table()[i]);
    }
  }
}  /* end of namespace getfem.                                            */

