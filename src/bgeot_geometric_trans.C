/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_geometric_trans.C : geometric transformations on convex*/
/*     									   */
/*                                                                         */
/* Date : December 20, 2000.                                               */
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


#include <bgeot_geometric_trans.h>
#include <dal_tree_sorted.h>

namespace bgeot
{

  /* ******************************************************************** */
  /* transformation on simplex.                                           */
  /* ******************************************************************** */

  struct _simplex_trans_light
  {
    dim_type N; short_type K;
    bool operator < (const _simplex_trans_light &l) const
    {
      if (N < l.N) return true; if (N > l.N) return false; 
      if (K < l.K) return true; return false;
    }
    _simplex_trans_light(dim_type n, short_type k) { N = n; K = k; }
    _simplex_trans_light(void) { }
  };

  struct _simplex_trans : public geometric_trans
  {
    void calc_base_func(base_poly &p, size_type i, short_type K) const
    {
      dim_type N = dim();
      base_poly l0(N, 0), l1(N, 0);
      power_index w(N+1);
      l0.one(); l1.one(); p = l0;
      for (int nn = 0; nn < N; ++nn) l0 -= base_poly(N, 1, nn);
      
      w[0] = K;
      for (int nn = 1; nn <= N; ++nn)
      { w[nn]=int(floor(0.5+((cvr->points()[i])[nn-1]*double(K)))); w[0]-=w[nn]; }
      
      for (int nn = 0; nn <= N; ++nn)
	for (int j = 0; j < w[nn]; ++j)
	  if (nn == 0)
	    p *= (l0 * (scalar_type(K) / scalar_type(j+1))) 
	       - (l1 * (scalar_type(j) / scalar_type(j+1)));
	  else
	    p *= (base_poly(N, 1, nn-1) * (scalar_type(K) / scalar_type(j+1))) 
	       - (l1 * (scalar_type(j) / scalar_type(j+1)));
    }

    _simplex_trans(const _simplex_trans_light &ls)
    {
      cvr = simplex_of_reference(ls.N, ls.K);
      size_type R = cvr->structure()->nb_points();
      is_lin = (ls.K == 1);
      trans.resize(R);
      
      for (size_type r = 0; r < R; ++r) calc_base_func(trans[r], r, ls.K);
    }
  };

  pgeometric_trans simplex_trans(dim_type n, short_type r)
  {
    static dal::FONC_TABLE<_simplex_trans_light, _simplex_trans> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_simplex_trans_light, _simplex_trans>();
      isinit = true;
    }
    return tab->add(_simplex_trans_light(n, r));
  }

  /* ******************************************************************** */
  /* direct product transformation                                        */
  /* ******************************************************************** */

  struct _cv_pr_t_light
  {
    pgeometric_trans cv1, cv2;
    
    bool operator < (const _cv_pr_t_light &ls) const
    {
      if (cv1 < ls.cv1) return true; if (cv1 > ls.cv1) return false; 
      if (cv2 < ls.cv2) return true; return false;
    }
    _cv_pr_t_light(pgeometric_trans a, pgeometric_trans b) { cv1=a; cv2=b; }
    _cv_pr_t_light(void) { }
  };

  struct _cv_pr_t : public geometric_trans
  {
    _cv_pr_t(const _cv_pr_t_light &ls)
    {
      cvr = convex_ref_product(ls.cv1->convex_ref(), ls.cv2->convex_ref());
      is_lin = false;

      size_type n1 = ls.cv1->nb_points(), n2 = ls.cv2->nb_points();
      trans.resize(n1 * n2);
      for (size_type i1 = 0; i1 < n1; ++i1)
	for (size_type i2 = 0; i2 < n2; ++i2)
	{
	  trans[i1 + i2 * n1] = ls.cv1->poly_vector()[i1];
	  trans[i1 + i2 * n1].direct_product(ls.cv2->poly_vector()[i2]);
	}
    }
  };

  pgeometric_trans convex_product_trans(pgeometric_trans a, pgeometric_trans b)
  {
    static dal::FONC_TABLE<_cv_pr_t_light, _cv_pr_t> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_cv_pr_t_light, _cv_pr_t>();
      isinit = true;
    }
    return tab->add(_cv_pr_t_light(a, b));
  }

  /* ******************************************************************** */
  /* linear direct product transformation.                                */
  /* ******************************************************************** */

  struct _cv_pr_tl_light
  {
    pgeometric_trans cv1, cv2;
    
    bool operator < (const _cv_pr_tl_light &ls) const
    {
      if (cv1 < ls.cv1) return true; if (cv1 > ls.cv1) return false; 
      if (cv2 < ls.cv2) return true; return false;
    }
    _cv_pr_tl_light(pgeometric_trans a, pgeometric_trans b) { cv1=a; cv2=b; }
    _cv_pr_tl_light(void) { }
  };

  struct _cv_pr_tl : public geometric_trans
  {
    _cv_pr_tl(const _cv_pr_tl_light &ls)
    {
      assert(ls.cv1->is_linear() && ls.cv2->is_linear());
      cvr = convex_ref_product(ls.cv1->convex_ref(), ls.cv2->convex_ref());
      is_lin = true;

      trans.resize(ls.cv1->nb_points() * ls.cv2->nb_points());

      for (size_type i = 0; i <= dim(); ++i)
	trans[cvr->structure()->ind_dir_points()[i]] 
	  = simplex_trans(dim(), 1)->poly_vector()[i];
    }
  };

  pgeometric_trans linear_product_trans(pgeometric_trans a,
					pgeometric_trans b)
  {
    static dal::FONC_TABLE<_cv_pr_tl_light, _cv_pr_tl> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_cv_pr_tl_light, _cv_pr_tl>();
      isinit = true;
    }
    return tab->add(_cv_pr_tl_light(a, b));
  }

  /* ******************************************************************** */
  /* parallelepiped transformation.                                       */
  /* ******************************************************************** */

  struct _cv_para_t_light
  {
    dim_type nc; short_type r;
    pgeometric_trans pgt;
    
    bool operator < (const _cv_para_t_light &ls) const
    {
      if (nc < ls.nc) return true; if (nc > ls.nc) return false; 
      if (r < ls.r) return true; return false;
    }
    _cv_para_t_light(dim_type nnc, short_type rr)
    { nc = nnc; r = rr; }
    _cv_para_t_light(void) { }
  };

  struct _cv_para_t
  {
    pgeometric_trans pgt;
    _cv_para_t(void);
    _cv_para_t(const _cv_para_t_light &l) {
      if (l.nc <= 1) 
	pgt = simplex_trans(l.nc, l.r);
      else
	pgt = convex_product_trans(parallelepiped_trans(l.nc-1, l.r),
				   simplex_trans(1, l.r));
    }
  };


  pgeometric_trans parallelepiped_trans(dim_type nc, short_type r)
  {
    static dal::FONC_TABLE<_cv_para_t_light, _cv_para_t> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_cv_para_t_light, _cv_para_t>();
      isinit = true;
    }
    return (tab->add(_cv_para_t_light(nc, r)))->pgt;
  }


  pgeometric_trans associated_trans(pconvex_structure cvs)
  {
    size_type n = cvs->dim(), nbp = cvs->nb_points();
    if (nbp = n+1)
      if (cvs == bgeot::simplex_structure(n))
	return simplex_trans(n, 1);

    if (nbp == (1 << n))
      if (cvs == bgeot::parallelepiped_structure(n))
	return parallelepiped_trans(n, 1);

    if (nbp == 2 * n)
      if (cvs == bgeot::prism_structure(n))
	return prism_trans(n, 1);
    
    // To be completed
    

    STD_NEEDED cerr << "This element is not taken into account. Contact us\n";
    assert(false);
    return NULL;
  }


}  /* end of namespace bgeot.                                            */

