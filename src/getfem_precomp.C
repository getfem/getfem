/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_precomp.C : pre-computations.                         */
/*     									   */
/*                                                                         */
/* Date : June 17, 2002.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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

#include <getfem_mat_elem.h>
#include <getfem_precomp.h>

namespace getfem
{



  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  struct pre_geot_light_
  {
    bgeot::pgeometric_trans pgt;
    bgeot::pstored_point_tab pspt;
    bool operator < (const pre_geot_light_ &ls) const
    {
      if (pgt < ls.pgt) return true; if (pgt > ls.pgt) return false; 
      if (pspt < ls.pspt) return true; return false;
    }
    pre_geot_light_(bgeot::pgeometric_trans pg, bgeot::pstored_point_tab ps)
    { pgt = pg; pspt = ps; }
    pre_geot_light_(void) { }
   
  };  

  geotrans_precomp_::geotrans_precomp_(const pre_geot_light_ &ls) {
    assign(ls);
  }

  geotrans_precomp_::geotrans_precomp_() {
    pgt = 0; pspt = 0;
  }

  void geotrans_precomp_::assign(const pre_geot_light_ &ls) {
    c.clear();
    pc.clear();
    hpc.clear();
    pgt = ls.pgt; pspt = ls.pspt;
  }

  void geotrans_precomp_::init_val() const {
    c.clear();  
    c.resize(pspt->size(), base_vector(pgt->nb_points()));
    for (size_type i = 0; i < pgt->nb_points(); ++i) {
      for (size_type j = 0; j < pspt->size(); ++j) {
	c[j][i] = pgt->poly_vector()[i].eval((*pspt)[j].begin());
      }
    }
  }

  void geotrans_precomp_::init_grad() const {
    dim_type N = pgt->structure()->dim();
    pc.clear(); 
    pc.resize(pspt->size(), base_matrix(pgt->nb_points() , N)); 
    for (size_type i = 0; i < pgt->nb_points(); ++i) {
      for (dim_type n = 0; n < N; ++n) {
	base_poly P = pgt->poly_vector()[i];
	P.derivative(n);
	for (size_type j = 0; j < pspt->size(); ++j) {
	  if ((*pspt)[j].size() != N)
	    DAL_THROW(dimension_error, "dimensions mismatch");
	  if(pgt->convex_ref()->is_in((*pspt)[j]) > 1.0E-7)
	    DAL_THROW(internal_error, "point " << j
		      << " mismatch the element");
	  pc[j](i,n) = P.eval((*pspt)[j].begin());
	}
      }
    }
  }

  void geotrans_precomp_::init_hess() const {
    base_poly P, Q;
    dim_type N = pgt->structure()->dim();
    hpc.clear();
    hpc.resize(pspt->size(), base_matrix(dal::sqr(N), pgt->nb_points()));
    for (size_type i = 0; i < pgt->nb_points(); ++i) {
      for (dim_type n = 0; n < N; ++n) {
	P = pgt->poly_vector()[i];
	P.derivative(n);
	for (dim_type m = 0; m <= n; ++m) {
	  Q = P; Q.derivative(m);
	  for (size_type j = 0; j < pspt->size(); ++j)
	    hpc[j](m * N + n, i) = hpc[j](n * N + m, i)
	      = P.eval((*pspt)[j].begin());
	}
      }
    }
  }

  base_node geotrans_precomp_::transform(size_type i,
					 const base_matrix &G) const {
    if (c.empty()) init_val();
    size_type N = G.nrows(), k = pgt->nb_points();
    base_node P(N); P.fill(0.0);
    base_matrix::const_iterator git = G.begin();
    for (size_type l = 0; l < k; ++l) {
      scalar_type a = c[i][l]; 
      base_node::iterator pit = P.begin(), pite = P.end();
      for (; pit != pite; ++git, ++pit) *pit += a * (*git);
    }
    return P;
  }

  pgeotrans_precomp geotrans_precomp(bgeot::pgeometric_trans pg,
				     bgeot::pstored_point_tab pspt)
  { 
    static dal::FONC_TABLE<pre_geot_light_, geotrans_precomp_> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<pre_geot_light_, geotrans_precomp_>();
      isinit = true;
    }
    return tab->add(pre_geot_light_(pg, pspt));
  }
  
  void geotrans_precomp_not_stored(bgeot::pgeometric_trans pg,
				   bgeot::pstored_point_tab pspt,
				   geotrans_precomp_& gp) {
    gp.assign(pre_geot_light_(pg, pspt));
  }

  /* ********************************************************************* */
  /*       Precomputation on fem.                                          */
  /* ********************************************************************* */

  struct pre_fem_light_
  {
    pfem pf;
    bgeot::pstored_point_tab pspt;
    bool operator < (const pre_fem_light_ &ls) const
    {
      if (pf < ls.pf) return true; if (pf > ls.pf) return false; 
      if (pspt < ls.pspt) return true; return false;
    }
    pre_fem_light_(pfem pff, bgeot::pstored_point_tab ps)
    { pf = pff; pspt = ps; }
    pre_fem_light_(void) { }
   
  };

  void fem_precomp_::assign(const pre_fem_light_ &ls) {
    c.clear();
    pc.clear();
    hpc.clear();
    pf = ls.pf; pspt = ls.pspt;
    for (size_type i = 0; i < pspt->size(); ++i)
      if ((*pspt)[i].size() != pf->structure()->dim())
	DAL_THROW(dimension_error, "dimensions mismatch");
  }

  fem_precomp_::fem_precomp_(const pre_fem_light_ &ls) { assign(ls); }

  fem_precomp_::fem_precomp_() : pf(0), pspt(0) {}
  
  void fem_precomp_::init_val() const {
    c.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i) 
      pf->base_value((*pspt)[i], c[i]);
  }

  void fem_precomp_::init_grad() const {
    pc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->grad_base_value((*pspt)[i], pc[i]);
  }

  void fem_precomp_::init_hess() const {
    hpc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->hess_base_value((*pspt)[i], hpc[i]);
  }

  typedef const fem_precomp_ * pfem_precomp;

  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt)
  { 
    static dal::FONC_TABLE<pre_fem_light_, fem_precomp_> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<pre_fem_light_, fem_precomp_>();
      isinit = true;
    }
    return tab->add(pre_fem_light_(pf, pspt));
  }

  
  void fem_precomp_not_stored(pfem pf, bgeot::pstored_point_tab pspt,
			      fem_precomp_& fp) {
    fp.assign(pre_fem_light_(pf,pspt));
  }
  

  /* fem_precomp_pool */
  class fem_precomp_pool_private : 
    public dal::FONC_TABLE<pre_fem_light_, fem_precomp_> {};
  fem_precomp_pool::fem_precomp_pool() : p(new fem_precomp_pool_private()) {}
  fem_precomp_pool::~fem_precomp_pool() { delete p; }
  pfem_precomp fem_precomp_pool::operator()(pfem pf, bgeot::pstored_point_tab pspt) {
    return p->add(pre_fem_light_(pf, pspt));
  }  
}  /* end of namespace getfem.                                            */

