/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_approx_integration.C : approximated  integration       */
/*               on convexes of reference.                                 */
/*                                                                         */
/* Date : December 17, 2000.                                               */
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


#include <bgeot_approx_integration.h>

namespace bgeot
{

  /* ********************************************************************* */
  /* method de Gauss.                                                      */
  /* ********************************************************************* */

  static dal::dynamic_array<base_poly> *Legendre_polynomials;
  static dal::dynamic_array< std::vector<scalar_type> > *Legendres_roots;
  static void init_legendre(short_type d)
  {
    static int nb_lp = -1;

    if (nb_lp < 0)
    {
      Legendre_polynomials = new dal::dynamic_array<base_poly>();
      Legendres_roots = new dal::dynamic_array< std::vector<scalar_type> >();
      (*Legendre_polynomials)[0] = base_poly(1,0);
      (*Legendre_polynomials)[0].one();
      (*Legendre_polynomials)[1] = base_poly(1,1,0);
      (*Legendres_roots)[1].resize(1);
      (*Legendres_roots)[1][0] = 0.0;
      nb_lp = 1;
    }
    while (nb_lp < d)
    {
      ++nb_lp;
      (*Legendre_polynomials)[nb_lp] =
	(base_poly(1,1,0) * (*Legendre_polynomials)[nb_lp-1]
	 * ((2.0 * scalar_type(nb_lp) - 1.0) / scalar_type(nb_lp)))
	+ ((*Legendre_polynomials)[nb_lp-2]
	* ((1.0 - scalar_type(nb_lp)) / scalar_type(nb_lp)));
      (*Legendres_roots)[nb_lp].resize(nb_lp);
      (*Legendres_roots)[nb_lp][nb_lp/2] = 0.0;
      scalar_type a = -1.0, b, c, d, e, cv, dv, ev, ecart, ecart2;
      for (int k = 0; k < nb_lp / 2; ++k) // + symetrie ...
      {
	b = (*Legendres_roots)[nb_lp-1][k];

	c = a, d = b;
	cv = (*Legendre_polynomials)[nb_lp].eval(&c);
	dv = (*Legendre_polynomials)[nb_lp].eval(&d);
	int imax = 10000;
	ecart = 2.0 * (d - c);
	while(c != d)
	{
	  --imax; if (imax == 0) break;
	  e = (c + d) / 2.0;
	  ecart2 = d - c; if (ecart2 >= ecart) break;
	  ecart = ecart2;
	  ev = (*Legendre_polynomials)[nb_lp].eval(&e);
	  if (ev * cv < 0) { d = e; dv = ev; } else { c = e; cv = ev; }
	}

	(*Legendres_roots)[nb_lp][k] = c;
	(*Legendres_roots)[nb_lp][nb_lp-k-1] = -c;
	a = b;
      }
    } 
  }

  struct _gauss_approx_integration : public approx_integration
  {

    _gauss_approx_integration(short_type nbpt)
    {
      assert(nbpt > 0);
      cvs = simplex_structure(1);
      stored_point_tab int_points(nbpt+2);
      int_coeffs.resize(nbpt+2);
      repartition.resize(3);
      repartition[0] = nbpt; 
      repartition[1] = nbpt + 1;
      repartition[2] = nbpt + 2; 

      init_legendre(nbpt);

      for (short_type i = 0; i < nbpt; ++i)
      {
	int_points[i].resize(1);
	scalar_type lr = (*Legendres_roots)[nbpt][i];
	int_points[i][0] = 0.5 + 0.5 * lr;
	int_coeffs[i] = (1.0 - dal::sqr(lr))
	  / dal::sqr( (nbpt) * ((*Legendre_polynomials)[nbpt-1].eval(&lr)));
      }
      
      int_points[nbpt].resize(1); /* a verifier. */
      int_points[nbpt][0] = 1.0; int_coeffs[nbpt] = 1.0;

      int_points[nbpt+1].resize(1);
      int_points[nbpt+1][0] = 0.0; int_coeffs[nbpt+1] = 1.0;
      pint_points = store_point_tab(int_points);

    }
  };

  papprox_integration Gauss_approx_integration(short_type nbpt)
  {
    static dal::dynamic_array<papprox_integration> *_gauss;
    static dal::bit_vector *_gauss_exists;
    static bool isinit = false;
    if (!isinit) {
      _gauss = new dal::dynamic_array<papprox_integration>();
      _gauss_exists = new dal::bit_vector();
      isinit = true;
    }

    if (!((*_gauss_exists)[nbpt]))
    {
      (*_gauss)[nbpt] = new _gauss_approx_integration(nbpt);
      (*_gauss_exists)[nbpt] = true;
    }
    return (*_gauss)[nbpt];
  }

  /* ********************************************************************* */
  /* integration on simplexes                                              */
  /* ********************************************************************* */

  struct _NC_apx_light
  {
    int n, k;
    bool operator < (const _NC_apx_light &ls) const
    {
      if (n < ls.n) return true; if (n > ls.n) return false; 
      if (k < ls.k) return true; return false;
    }
    _NC_apx_light(int nn, int kk) { n = nn; k = kk; }
    _NC_apx_light(void) { }
   
  };


  struct _Newton_Cotes_approx_integration : public approx_integration
  {

    void calc_base_func(base_poly &p, size_type i, short_type K, base_node &c)
      const
    {
      dim_type N = dim();
      base_poly l0(N, 0), l1(N, 0);
      power_index w(N+1);
      l0.one(); l1.one(); p = l0;
      for (int nn = 0; nn < N; ++nn) l0 -= base_poly(N, 1, nn);
      
      w[0] = K;
      for (int nn = 1; nn <= N; ++nn)
      { w[nn]=int(floor(0.5+(c[nn-1]*double(K)))); w[0]-=w[nn]; }
      
      for (int nn = 0; nn <= N; ++nn)
	for (int j = 0; j < w[nn]; ++j)
	  if (nn == 0)
	    p *= (l0 * (scalar_type(K) / scalar_type(j+1))) 
	       - (l1 * (scalar_type(j) / scalar_type(j+1)));
	  else
	    p *= (base_poly(N, 1, nn-1) * (scalar_type(K) / scalar_type(j+1))) 
	       - (l1 * (scalar_type(j) / scalar_type(j+1)));
    }

    _Newton_Cotes_approx_integration(const _NC_apx_light &ls)
    {
      cvs = simplex_structure(ls.n);
      size_type R = alpha(ls.n,ls.k);
      size_type R2 = (ls.n > 0) ? alpha(ls.n-1,ls.k) : 0;
      base_poly P;
      ppoly_integration ppi = simplex_poly_integration(ls.n);
      std::vector<size_type> fa(ls.n+1);
      stored_point_tab int_points;
      int_points.resize(R + (ls.n+1) * R2);
      int_coeffs.resize(R + (ls.n+1) * R2);
      repartition.resize(ls.n+2);
      repartition[0] = R;
      for (short_type f = 0; f <= ls.n; ++f)
      {
	fa[f] = repartition[f];
	repartition[f+1] = repartition[f] + R2;
      }

      size_type sum = 0, l;
      base_node c(ls.n); 
      c *= scalar_type(0.0);
      if (ls.k == 0) c.fill(1.0 / scalar_type(ls.n+1));

      for (size_type r = 0; r < R; ++r)
      {
	int_points[r] = c;
	calc_base_func(P, r, ls.k, c);
	int_coeffs[r] = ppi->int_poly(P);

	for (short_type f = 1; f <= ls.n; ++f)
	  if (c[f-1] == 0.0)
	  {
	    int_points[fa[f]] = c;
	    int_coeffs[fa[f]] = ppi->int_poly_on_face(P, f);
	    (fa[f])++;
	  }
	if (ls.k == 0)
	{
	  for (short_type f = 0; f <= ls.n; ++f)
	  {
	    int_points[fa[f]].resize(dim());
	    int_points[fa[f]].fill(1.0 / scalar_type(ls.n));
	    if (f > 0) int_points[fa[f]][f-1] = 0.0;
	    int_coeffs[fa[f]] = ppi->int_poly_on_face(P, f);
	  }
	}
	else
	if (sum == ls.k)
	{
	  int_points[fa[0]] = c;
	  int_coeffs[fa[0]] = ppi->int_poly_on_face(P, 0);
	  (fa[0])++;
	}  

	l = 0; c[l] += 1.0 / scalar_type(ls.k); sum++;
	while (sum > ls.k)
	{
	  sum -= int(floor(0.5+(c[l] * ls.k)));
	  c[l] = 0.0; l++; if (l == ls.n) break;
	  c[l] += 1.0 / scalar_type(ls.k); sum++;
	}
      }
      pint_points = store_point_tab(int_points);
    }
  };


  


  papprox_integration Newton_Cotes_approx_integration(dim_type n,
						      short_type k)
  { 
    static dal::FONC_TABLE<_NC_apx_light, _Newton_Cotes_approx_integration> *t;
    static bool isinit = false;
    if (!isinit) {
      t= new dal::FONC_TABLE<_NC_apx_light,_Newton_Cotes_approx_integration>();
      isinit = true;
    }
    return t->add(_NC_apx_light(n, k));
  }

  /* ********************************************************************* */
  /* integration on direct product of convex structures                    */
  /* ********************************************************************* */
  
  struct a_int_pro_light
  {
    papprox_integration cv1, cv2;
    bool operator < (const a_int_pro_light &ls) const
    {
      if (cv1 < ls.cv1) return true; if (cv1 > ls.cv1) return false; 
      if (cv2 < ls.cv2) return true; return false;
    }
    a_int_pro_light(papprox_integration a, papprox_integration b)
    { cv1 = a; cv2 = b; }
    a_int_pro_light(void) { }
   
  };


  struct a_int_pro_integration : public approx_integration
  {

    a_int_pro_integration(const a_int_pro_light &ls)
    {
      cvs = convex_product_structure(ls.cv1->structure(), ls.cv2->structure());
      size_type n1 = ls.cv1->nb_points_on_convex();
      size_type n2 = ls.cv2->nb_points_on_convex();
      stored_point_tab int_points;
      int_points.resize(n1 * n2);
      int_coeffs.resize(n1 * n2);
      repartition.resize(cvs->nb_faces()+1);
      repartition[0] = n1 * n2;
      for (size_type i1 = 0; i1 < n1; ++i1)
	for (size_type i2 = 0; i2 < n2; ++i2)
	{
	  int_coeffs[i1 + i2 * n1] = ls.cv1->coeff(i1) * ls.cv2->coeff(i2);
	  int_points[i1 + i2 * n1].resize(dim());
	  std::copy(ls.cv1->point(i1).begin(), ls.cv1->point(i1).end(),
		    int_points[i1 + i2 * n1].begin());
	  std::copy(ls.cv2->point(i2).begin(), ls.cv2->point(i2).end(),
		    int_points[i1 + i2 * n1].begin() + ls.cv1->dim());
	}
      short_type f = 0;
      for (short_type f1 = 0; f1 < ls.cv1->structure()->nb_faces(); ++f1, ++f)
      {
	n1 = ls.cv1->nb_points_on_face(f1);
	n2 = ls.cv2->nb_points_on_convex();
	size_type w = repartition[f];
	repartition[f+1] = w + n1 * n2;
	int_points.resize(repartition[f+1]);
	int_coeffs.resize(repartition[f+1]);
	for (size_type i1 = 0; i1 < n1; ++i1)
	  for (size_type i2 = 0; i2 < n2; ++i2)
	  {
	    int_coeffs[w + i1 + i2 * n1] = ls.cv1->coeff_on_face(f1, i1)
	                                 * ls.cv2->coeff(i2);
	    int_points[w + i1 + i2 * n1].resize(dim());
	    std::copy(ls.cv1->point_on_face(f1, i1).begin(),
		      ls.cv1->point_on_face(f1, i1).end(),
		      int_points[w + i1 + i2 * n1].begin());
	    std::copy(ls.cv2->point(i2).begin(), ls.cv2->point(i2).end(),
		      int_points[w + i1 + i2 * n1].begin() + ls.cv1->dim());
	  }
      }
      for (short_type f2 = 0; f2 < ls.cv2->structure()->nb_faces(); ++f2, ++f)
      {
	n1 = ls.cv1->nb_points_on_convex();
	n2 = ls.cv2->nb_points_on_face(f2);
	size_type w = repartition[f];
	repartition[f+1] = w + n1 * n2;
	int_points.resize(repartition[f+1]);
	int_coeffs.resize(repartition[f+1]);
	for (size_type i1 = 0; i1 < n1; ++i1)
	  for (size_type i2 = 0; i2 < n2; ++i2)
	  {
	    int_coeffs[w + i1 + i2 * n1] = ls.cv1->coeff(i1)
	                                 * ls.cv2->coeff_on_face(f2, i2);
	    int_points[w + i1 + i2 * n1].resize(dim());
	    std::copy(ls.cv1->point(i1).begin(), ls.cv1->point(i1).end(),
		      int_points[w + i1 + i2 * n1].begin());
	    std::copy(ls.cv2->point_on_face(f2, i2).begin(),
		      ls.cv2->point_on_face(f2, i2).end(),
		      int_points[w + i1 + i2 * n1].begin() + ls.cv1->dim());
	  }
      }
      pint_points = store_point_tab(int_points);
    }
  };

  papprox_integration convex_product_approx_integration(papprox_integration a,
							papprox_integration b)
  {
    static dal::FONC_TABLE<a_int_pro_light, a_int_pro_integration> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<a_int_pro_light, a_int_pro_integration>();
      isinit = true;
    }
    return tab->add(a_int_pro_light(a, b));
  }

}  /* end of namespace bgeot.                                            */

