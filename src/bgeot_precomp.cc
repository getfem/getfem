//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : bgeot_precomp.cc : pre-computations.
//           
// Date    : 2004/01/11.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2002-2005 Yves Renard
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
#include <bgeot_precomp.h>

namespace bgeot {
  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  struct pre_geot_light_
  {
    bgeot::pgeometric_trans pgt;
    bgeot::pstored_point_tab pspt;
    geotrans_precomp_pool *pool;
    bool operator < (const pre_geot_light_ &ls) const
    {
      if (pgt < ls.pgt) return true; if (pgt > ls.pgt) return false; 
      if (pspt < ls.pspt) return true; return false;
    }
    pre_geot_light_(bgeot::pgeometric_trans pg, 
		    bgeot::pstored_point_tab ps, 
		    geotrans_precomp_pool *p) 
      : pgt(pg), pspt(ps), pool(p) {}
    pre_geot_light_(void) { }   
  };  

  geotrans_precomp_::geotrans_precomp_(const pre_geot_light_ &ls) 
    : pgt(ls.pgt), pspt(ls.pspt), pool_(ls.pool) {}

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
    hpc.resize(pspt->size(), base_matrix(gmm::sqr(N), pgt->nb_points()));
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

  /* geotrans_precomp_pool */
  class geotrans_precomp_pool_private : 
    public dal::FONC_TABLE<pre_geot_light_, geotrans_precomp_> {};

  geotrans_precomp_pool::geotrans_precomp_pool() : 
    p(new geotrans_precomp_pool_private()) {}

  void geotrans_precomp_pool::clear() { 
    delete p; p = new geotrans_precomp_pool_private(); 
  }
  geotrans_precomp_pool::~geotrans_precomp_pool() { delete p; }

  pgeotrans_precomp 
  geotrans_precomp_pool::operator()(bgeot::pgeometric_trans pg,
				    bgeot::pstored_point_tab pspt) { 
    return p->add(pre_geot_light_(pg, pspt, this));
  }

  pgeotrans_precomp geotrans_precomp(bgeot::pgeometric_trans pg,
				     bgeot::pstored_point_tab pspt) { 
    static geotrans_precomp_pool *pool = 0;
    if (!pool) pool = new geotrans_precomp_pool();
    return (*pool)(pg, pspt);
  }
}
