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
    base_poly P, Q;
    dim_type N = ls.pgt->structure()->dim();
    pgt = ls.pgt; pspt = ls.pspt;
    pc.resize(pspt->size());
    hpc.resize(pspt->size());
    std::fill(pc.begin(), pc.end(),
	      base_matrix(pgt->nb_points() , N));
    std::fill(hpc.begin(), hpc.end(),
	      base_matrix(dal::sqr(N), pgt->nb_points()));
    for (size_type i = 0; i < pgt->nb_points(); ++i)
      for (dim_type n = 0; n < N; ++n) {
	P = pgt->poly_vector()[i];
	P.derivative(n);
	for (size_type j = 0; j < pspt->size(); ++j) {
	  if ((*pspt)[j].size() != N)
	    DAL_THROW(dimension_error, "dimensions mismatch");
	  if(pgt->convex_ref()->is_in((*pspt)[j]) > 1.0E-7)
	    DAL_THROW(internal_error, "point " << j
		      << " mismatch the element");
	  pc[j](i,n) = P.eval((*pspt)[j].begin());
	}
	for (dim_type m = 0; m <= n; ++m) {
	  Q = P; Q.derivative(m);
	  for (size_type j = 0; j < pspt->size(); ++j)
	    hpc[j](m * N + n, i) = hpc[j](n * N + m, i)
	      = P.eval((*pspt)[j].begin());
	}
      }
    
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

  _fem_precomp::_fem_precomp(const _pre_fem_light &ls)
    : pf(ls.pf), pspt(ls.pspt) {
    dim_type N = ls.pf->structure()->dim();
    size_type npt = ls.pspt->size();
    // cout << "Dimension = " << N << endl;
    pc.resize(npt);
    hpc.resize(npt);
    c.resize(npt);
    // cout << "taille : " << npt << endl;
    for (size_type i = 0; i < npt; ++i) {
      if ((*ls.pspt)[i].size() != N)
	DAL_THROW(dimension_error, "dimensions mismatch");
      ls.pf->base_value((*(ls.pspt))[i], c[i]);
      ls.pf->grad_base_value((*(ls.pspt))[i], pc[i]);
      ls.pf->hess_base_value((*(ls.pspt))[i], hpc[i]);
      // cout << "---------------------------------------------------------\n";
      // cout << "point " << (*ls.pspt)[i] << " value : " << c[i] << " grad : " << pc[i] << endl;
    }
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

}  /* end of namespace getfem.                                            */

