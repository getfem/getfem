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

  struct _pre_geot_light
  {
    bgeot::pgeometric_trans pgt;
    bgeot::pstored_point_tab pspt;
    bool operator < (const _pre_geot_light &ls) const
    {
      if (pgt < ls.pgt) return true; if (pgt > ls.pgt) return false; 
      if (pspt < ls.pspt) return true; return false;
    }
    _pre_geot_light(bgeot::pgeometric_trans pg, bgeot::pstored_point_tab ps)
    { pgt = pg; pspt = ps; }
    _pre_geot_light(void) { }
   
  };  

  _geotrans_precomp::_geotrans_precomp(const _pre_geot_light &ls) {
    assign(ls);
  }

  _geotrans_precomp::_geotrans_precomp() {
    pgt = 0; pspt = 0;
  }

  void _geotrans_precomp::assign(const _pre_geot_light &ls) {
    c.clear();
    pc.clear();
    hpc.clear();
    pgt = ls.pgt; pspt = ls.pspt;
  }

  void _geotrans_precomp::init_val() const {
    c.clear();  
    c.resize(pspt->size(), base_vector(pgt->nb_points()));
    for (size_type i = 0; i < pgt->nb_points(); ++i) {
      for (size_type j = 0; j < pspt->size(); ++j) {
	c[j][i] = pgt->poly_vector()[i].eval((*pspt)[j].begin());
      }
    }
  }

  void _geotrans_precomp::init_grad() const {
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

  void _geotrans_precomp::init_hess() const {
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

  base_node _geotrans_precomp::transform(size_type i,
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
    static dal::FONC_TABLE<_pre_geot_light, _geotrans_precomp> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_pre_geot_light, _geotrans_precomp>();
      isinit = true;
    }
    return tab->add(_pre_geot_light(pg, pspt));
  }
  
  void geotrans_precomp_not_stored(bgeot::pgeometric_trans pg,
				   bgeot::pstored_point_tab pspt,
				   _geotrans_precomp& gp) {
    gp.assign(_pre_geot_light(pg, pspt));
  }

  /* ********************************************************************* */
  /*       Precomputation on fem.                                          */
  /* ********************************************************************* */

  struct _pre_fem_light
  {
    pfem pf;
    bgeot::pstored_point_tab pspt;
    bool operator < (const _pre_fem_light &ls) const
    {
      if (pf < ls.pf) return true; if (pf > ls.pf) return false; 
      if (pspt < ls.pspt) return true; return false;
    }
    _pre_fem_light(pfem pff, bgeot::pstored_point_tab ps)
    { pf = pff; pspt = ps; }
    _pre_fem_light(void) { }
   
  };

  void _fem_precomp::assign(const _pre_fem_light &ls) {
    c.clear();
    pc.clear();
    hpc.clear();
    pf = ls.pf; pspt = ls.pspt;
    for (size_type i = 0; i < pspt->size(); ++i)
      if ((*pspt)[i].size() != pf->structure()->dim())
	DAL_THROW(dimension_error, "dimensions mismatch");
  }

  _fem_precomp::_fem_precomp(const _pre_fem_light &ls) { assign(ls); }

  _fem_precomp::_fem_precomp() : pf(0), pspt(0) {}
  
  void _fem_precomp::init_val() const {
    c.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i) 
      pf->base_value((*pspt)[i], c[i]);
  }

  void _fem_precomp::init_grad() const {
    pc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->grad_base_value((*pspt)[i], pc[i]);
  }

  void _fem_precomp::init_hess() const {
    hpc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->hess_base_value((*pspt)[i], hpc[i]);
  }

  typedef const _fem_precomp * pfem_precomp;

  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt)
  { 
    static dal::FONC_TABLE<_pre_fem_light, _fem_precomp> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_pre_fem_light, _fem_precomp>();
      isinit = true;
    }
    return tab->add(_pre_fem_light(pf, pspt));
  }

  
  void fem_precomp_not_stored(pfem pf, bgeot::pstored_point_tab pspt,
			      _fem_precomp& fp) {
    fp.assign(_pre_fem_light(pf,pspt));
  }
  

  /* fem_precomp_pool */
  class fem_precomp_pool_private : 
    public dal::FONC_TABLE<_pre_fem_light, _fem_precomp> {};
  fem_precomp_pool::fem_precomp_pool() : p(new fem_precomp_pool_private()) {}
  fem_precomp_pool::~fem_precomp_pool() { delete p; }
  pfem_precomp fem_precomp_pool::operator()(pfem pf, bgeot::pstored_point_tab pspt) {
    return p->add(_pre_fem_light(pf, pspt));
  }  
}  /* end of namespace getfem.                                            */

