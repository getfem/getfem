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
      scalar_type a = -1.0, b, c, d, e, cv, ev, ecart, ecart2;
      for (int k = 0; k < nb_lp / 2; ++k) // + symetrie ...
      {
	b = (*Legendres_roots)[nb_lp-1][k];

	c = a, d = b;
	cv = (*Legendre_polynomials)[nb_lp].eval(&c);
	int imax = 10000;
	ecart = 2.0 * (d - c);
	while(c != d)
	{
	  --imax; if (imax == 0) break;
	  e = (c + d) / 2.0;
	  ecart2 = d - c; if (ecart2 >= ecart) break;
	  ecart = ecart2;
	  ev = (*Legendre_polynomials)[nb_lp].eval(&e);
	  if (ev * cv < 0) { d = e; } else { c = e; cv = ev; }
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
      if (nbpt > 32000) DAL_THROW(std::out_of_range, "too much points");
      
      cvr = simplex_of_reference(1);
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
    dim_type n; short_type k;
    bool operator < (const _NC_apx_light &ls) const
    {
      if (n < ls.n) return true; if (n > ls.n) return false; 
      if (k < ls.k) return true; return false;
    }
    _NC_apx_light(dim_type nn, short_type kk) { n = nn; k = kk; }
    _NC_apx_light(void) { }
   
  };


  struct _Newton_Cotes_approx_integration : public approx_integration
  {

    void calc_base_func(base_poly &p, short_type K, base_node &c)
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
      cvr = simplex_of_reference(ls.n);
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
	calc_base_func(P, ls.k, c);
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

	if (ls.k != 0) {
	  l = 0; c[l] += 1.0 / scalar_type(ls.k); sum++;
	  while (sum > ls.k) {
	    sum -= int(floor(0.5+(c[l] * ls.k)));
	    c[l] = 0.0; l++; if (l == ls.n) break;
	    c[l] += 1.0 / scalar_type(ls.k); sum++;
	  }
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
      cvr = convex_ref_product(ls.cv1->ref_convex(), ls.cv2->ref_convex());
      size_type n1 = ls.cv1->nb_points_on_convex();
      size_type n2 = ls.cv2->nb_points_on_convex();
      stored_point_tab int_points;
      int_points.resize(n1 * n2);
      int_coeffs.resize(n1 * n2);
      repartition.resize(cvr->structure()->nb_faces()+1);
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

  /* ********************************************************************* */
  /* integration on parallelepiped with Newton Cotes formulae              */
  /* ********************************************************************* */
  
  struct a_par_Ne_light
  {
    dim_type N; short_type K;
    bool operator < (const a_par_Ne_light &ls) const
    {
      if (N < ls.N) return true; if (N > ls.N) return false; 
      if (K < ls.K) return true; return false;
    }
    a_par_Ne_light(dim_type NN, short_type KK) { N = NN; K = KK; }
    a_par_Ne_light(void) { }
  };

  struct a_par_Ne
  {
    papprox_integration pai;
    a_par_Ne(const a_par_Ne_light &ls)
    {
      papprox_integration aux;
      pai = aux = Newton_Cotes_approx_integration(1, ls.K);
      for (int i = 1; i < ls.N; ++i)
	pai = convex_product_approx_integration(pai, aux);
    }
  };

  papprox_integration parallelepiped_Newton_Cotes_approx_integration
  (dim_type N, short_type K)
  {
    static dal::FONC_TABLE<a_par_Ne_light, a_par_Ne> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<a_par_Ne_light, a_par_Ne>();
      isinit = true;
    }
    return (tab->add(a_par_Ne_light(N, K)))->pai;
  }

  /* ********************************************************************* */
  /* integration on parallelepiped with Gauss formulae                     */
  /* ********************************************************************* */
  
  struct a_par_Gauss_light
  {
    dim_type N; short_type K;
    bool operator < (const a_par_Gauss_light &ls) const
    {
      if (N < ls.N) return true; if (N > ls.N) return false; 
      if (K < ls.K) return true; return false;
    }
    a_par_Gauss_light(dim_type NN, short_type KK) { N = NN; K = KK;}
    a_par_Gauss_light(void) { }
  };

  struct a_par_Gauss
  {
    papprox_integration pai;
    a_par_Gauss(const a_par_Gauss_light &ls)
    {
      papprox_integration aux;
      pai = aux = Gauss_approx_integration(ls.K);
      for (int i = 1; i < ls.N; ++i)
	pai = convex_product_approx_integration(pai, aux);
    }
  };

  papprox_integration parallelepiped_Gauss_approx_integration
  (dim_type N, short_type K)
  {
    static dal::FONC_TABLE<a_par_Gauss_light, a_par_Gauss> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<a_par_Gauss_light, a_par_Gauss>();
      isinit = true;
    }
    return (tab->add(a_par_Gauss_light(N, K)))->pai;
  }

  /* ********************************************************************* */
  /*                                                                       */
  /*   Particular integration methods on dimension 2.                      */
  /*                                                                       */
  /* ********************************************************************* */

  /* ********************************************************************* */
  /*   triangle1 :    Integration on a triangle of order 1 with 1 point    */
  /* ********************************************************************* */

  struct _particular_approx : public approx_integration
  {
    friend papprox_integration triangle1_approx_integration(void);
    friend papprox_integration triangle2_approx_integration(void);
    friend papprox_integration triangle2bis_approx_integration(void);
    friend papprox_integration triangle3_approx_integration(void);
    friend papprox_integration triangle4_approx_integration(void);
    friend papprox_integration triangle5_approx_integration(void);
    friend papprox_integration triangle6_approx_integration(void);
    friend papprox_integration triangle7_approx_integration(void);
    friend papprox_integration quad2_approx_integration(void);
    friend papprox_integration quad3_approx_integration(void);    
    friend papprox_integration quad5_approx_integration(void);
    friend papprox_integration tetrahedron1_approx_integration(void);
    friend papprox_integration tetrahedron2_approx_integration(void);
    friend papprox_integration tetrahedron3_approx_integration(void);
    friend papprox_integration tetrahedron5_approx_integration(void);
  };

  papprox_integration triangle1_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(4);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      ptab[i][0] = 1.0 / 3.0; ptab[i][1] = 1.0 / 3.0;
      p->int_coeffs[i] = 0.5; 
      p->repartition[0] = 1;
      // face 0
      ptab[++i][0] = 0.5; ptab[i][1] = 0.5;
      p->int_coeffs[i] = ::sqrt(2.0); 
      p->repartition[1] = p->repartition[0] + 1;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = 0.5;
      p->int_coeffs[i] = 1.0; 
      p->repartition[2] = p->repartition[1] + 1;
      // face 2
      ptab[++i][0] = 0.5; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 1.0; 
      p->repartition[3] = p->repartition[2] + 1;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }

  /* ********************************************************************* */
  /*   triangle2 :    Integration on a triangle of order 2 with 3 points   */
  /* ********************************************************************* */

    papprox_integration triangle2_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(9);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      ptab[  i][0] = 1.0 / 6.0; ptab[i][1] = 1.0 / 6.0;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = 2.0 / 3.0; ptab[i][1] = 1.0 / 6.0;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = 1.0 / 6.0; ptab[i][1] = 2.0 / 3.0;
      p->int_coeffs[i] = 1.0 / 6.0; 
      p->repartition[0] = 3;
      // face 0
      double a = 0.5 - 0.5/::sqrt(3.0);
      double b = 0.5 + 0.5/::sqrt(3.0);
      ptab[++i][0] = a; ptab[i][1] = b;
      p->int_coeffs[i] = ::sqrt(2.0) * 0.5; 
      ptab[++i][0] = b; ptab[i][1] = a;
      p->int_coeffs[i] = ::sqrt(2.0) * 0.5; 
      p->repartition[1] = p->repartition[0] + 2;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = a;
      p->int_coeffs[i] = 0.5; 
      ptab[++i][0] = 0.0; ptab[i][1] = b;
      p->int_coeffs[i] = 0.5; 
      p->repartition[2] = p->repartition[1] + 2;
      // face 2
      ptab[++i][0] = a; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 0.5; 
      ptab[++i][0] = b; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 0.5; 
      p->repartition[3] = p->repartition[2] + 2;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }

  /* ********************************************************************* */
  /*   triangle2bis :   Integration on a triangle of order 2 with 3 points */
  /* ********************************************************************* */

    papprox_integration triangle2bis_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(9);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      ptab[  i][0] = 0.5; ptab[i][1] = 0.5;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = 0; ptab[i][1] = 0.5;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = 0.5; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 1.0 / 6.0; 
      p->repartition[0] = 3;
      // face 0
      double a = 0.5 - 0.5/::sqrt(3.0);
      double b = 0.5 + 0.5/::sqrt(3.0);
      ptab[++i][0] = a; ptab[i][1] = b;
      p->int_coeffs[i] = ::sqrt(2.0) * 0.5; 
      ptab[++i][0] = b; ptab[i][1] = a;
      p->int_coeffs[i] = ::sqrt(2.0) * 0.5; 
      p->repartition[1] = p->repartition[0] + 2;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = a;
      p->int_coeffs[i] = 0.5; 
      ptab[++i][0] = 0.0; ptab[i][1] = b;
      p->int_coeffs[i] = 0.5; 
      p->repartition[2] = p->repartition[1] + 2;
      // face 2
      ptab[++i][0] = a; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 0.5; 
      ptab[++i][0] = b; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 0.5; 
      p->repartition[3] = p->repartition[2] + 2;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }

  /* ********************************************************************* */
  /*   triangle3 :    Integration on a triangle of order 3 with 4 points   */
  /* ********************************************************************* */

  papprox_integration triangle3_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(10);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      ptab[  i][0] = 1.0 / 3.0; ptab[i][1] = 1.0 / 3.0;
      p->int_coeffs[i] = -27.0 / 96.0; 
      ptab[++i][0] = 1.0 / 5.0; ptab[i][1] = 1.0 / 5.0;
      p->int_coeffs[i] = 25.0 / 96.0; 
      ptab[++i][0] = 3.0 / 5.0; ptab[i][1] = 1.0 / 5.0;
      p->int_coeffs[i] = 25.0 / 96.0; 
      ptab[++i][0] = 1.0 / 5.0; ptab[i][1] = 3.0 / 5.0;
      p->int_coeffs[i] = 25.0 / 96.0; 
      p->repartition[0] = 4;
      // face 0
      double a = 0.5 - 0.5/::sqrt(3.0);
      double b = 0.5 + 0.5/::sqrt(3.0);
      ptab[++i][0] = a; ptab[i][1] = b;
      p->int_coeffs[i] = ::sqrt(2.0) * 0.5; 
      ptab[++i][0] = b; ptab[i][1] = a;
      p->int_coeffs[i] = ::sqrt(2.0) * 0.5; 
      p->repartition[1] = p->repartition[0] + 2;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = a;
      p->int_coeffs[i] = 0.5; 
      ptab[++i][0] = 0.0; ptab[i][1] = b;
      p->int_coeffs[i] = 0.5; 
      p->repartition[2] = p->repartition[1] + 2;
      // face 2
      ptab[++i][0] = a; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 0.5; 
      ptab[++i][0] = b; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 0.5; 
      p->repartition[3] = p->repartition[2] + 2;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }


  /* ********************************************************************* */
  /*   triangle4 :    Integration on a triangle of order 4 with 6 points   */
  /* ********************************************************************* */

  papprox_integration triangle4_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(15);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      double a = 0.445948490915965;
      double b = 0.091576213509771;
      double c = 0.111690794839005;
      double d = 0.054975871827661;
      ptab[  i][0] = a; ptab[i][1] = a;
      p->int_coeffs[i] = c; 
      ptab[++i][0] = 1.0 - 2.0 * a; ptab[i][1] = a;
      p->int_coeffs[i] = c; 
      ptab[++i][0] = a; ptab[i][1] = 1.0 - 2.0 * a;
      p->int_coeffs[i] = c; 
      ptab[++i][0] = b; ptab[i][1] = b;
      p->int_coeffs[i] = d; 
      ptab[++i][0] = 1.0 - 2.0 * b; ptab[i][1] = b;
      p->int_coeffs[i] = d; 
      ptab[++i][0] = b; ptab[i][1] = 1.0 - 2.0 * b;
      p->int_coeffs[i] = d; 
      p->repartition[0] = 6;
      // face 0
      double e = 0.5 - 0.5*::sqrt(3.0 / 5.0);
      double f = 0.5 + 0.5*::sqrt(3.0 / 5.0);
      ptab[++i][0] = e; ptab[i][1] = f;
      p->int_coeffs[i] = ::sqrt(2.0) * 5.0 / 18.0; 
      ptab[++i][0] = 0.5; ptab[i][1] = 0.5;
      p->int_coeffs[i] = ::sqrt(2.0) * 8.0 / 18.0; 
      ptab[++i][0] = f; ptab[i][1] = e;
      p->int_coeffs[i] = ::sqrt(2.0) * 5.0 / 18.0; 
      p->repartition[1] = p->repartition[0] + 3;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = e;
      p->int_coeffs[i] = 5.0 / 18.0; 
      ptab[++i][0] = 0.0; ptab[i][1] = 0.5;
      p->int_coeffs[i] = 8.0 / 18.0; 
      ptab[++i][0] = 0.0; ptab[i][1] = f;
      p->int_coeffs[i] = 5.0 / 18.0; 
      p->repartition[2] = p->repartition[1] + 3;
      // face 2
      ptab[++i][0] = e; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 5.0 / 18.0; 
      ptab[++i][0] = 0.5; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 8.0 / 18.0; 
      ptab[++i][0] = f; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 5.0 / 18.0; 
      p->repartition[3] = p->repartition[2] + 3;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }

  /* ********************************************************************* */
  /*   triangle5 :    Integration on a triangle of order 5 with 7 points   */
  /* ********************************************************************* */

  papprox_integration triangle5_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(16);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      double a = 0.470142064105115;
      double b = 0.101286507323456;
      double c = 0.0661970763942530;
      double d = 0.0629695902724135;
      ptab[  i][0] = a; ptab[i][1] = a;
      p->int_coeffs[i] = c; 
      ptab[++i][0] = 1.0 - 2.0 * a; ptab[i][1] = a;
      p->int_coeffs[i] = c; 
      ptab[++i][0] = a; ptab[i][1] = 1.0 - 2.0 * a;
      p->int_coeffs[i] = c; 
      ptab[++i][0] = b; ptab[i][1] = b;
      p->int_coeffs[i] = d; 
      ptab[++i][0] = 1.0 - 2.0 * b; ptab[i][1] = b;
      p->int_coeffs[i] = d; 
      ptab[++i][0] = b; ptab[i][1] = 1.0 - 2.0 * b;
      p->int_coeffs[i] = d; 
      ptab[++i][0] = 1.0 / 3.0; ptab[i][1] = 1.0 / 3.0;
      p->int_coeffs[i] = 9.0 / 80.0; 
      p->repartition[0] = 7;
      // face 0
      double e = 0.5 - 0.5*::sqrt(3.0 / 5.0);
      double f = 0.5 + 0.5*::sqrt(3.0 / 5.0);
      ptab[++i][0] = e; ptab[i][1] = f;
      p->int_coeffs[i] = ::sqrt(2.0) * 5.0 / 18.0; 
      ptab[++i][0] = 0.5; ptab[i][1] = 0.5;
      p->int_coeffs[i] = ::sqrt(2.0) * 8.0 / 18.0; 
      ptab[++i][0] = f; ptab[i][1] = e;
      p->int_coeffs[i] = ::sqrt(2.0) * 5.0 / 18.0; 
      p->repartition[1] = p->repartition[0] + 3;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = e;
      p->int_coeffs[i] = 5.0 / 18.0; 
      ptab[++i][0] = 0.0; ptab[i][1] = 0.5;
      p->int_coeffs[i] = 8.0 / 18.0; 
      ptab[++i][0] = 0.0; ptab[i][1] = f;
      p->int_coeffs[i] = 5.0 / 18.0; 
      p->repartition[2] = p->repartition[1] + 3;
      // face 2
      ptab[++i][0] = e; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 5.0 / 18.0; 
      ptab[++i][0] = 0.5; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 8.0 / 18.0; 
      ptab[++i][0] = f; ptab[i][1] = 0.0;
      p->int_coeffs[i] = 5.0 / 18.0; 
      p->repartition[3] = p->repartition[2] + 3;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }


  /* ********************************************************************* */
  /*   triangle6 :    Integration on a triangle of order 6 with 12 points  */
  /* ********************************************************************* */

  papprox_integration triangle6_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(24);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      size_type i = 0;

      double a1 = 0.063089104491502;
      double a2 = 0.249286745170910;
      double aa = 0.310352451033785;
      double bb = 0.053145049844816;

      ptab[  i][0] = a1; ptab[i][1] = a1;
      p->int_coeffs[i] = 0.050844906370206 * 0.5; 
      ptab[++i][0] = 1.0 - 2.0 * a1; ptab[i][1] = a1;
      p->int_coeffs[i] = 0.050844906370206 * 0.5; 
      ptab[++i][0] = a1; ptab[i][1] = 1.0 - 2.0 * a1;
      p->int_coeffs[i] = 0.050844906370206 * 0.5; 
      ptab[++i][0] = a2; ptab[i][1] = a2;
      p->int_coeffs[i] = 0.116786275726378 * 0.5; 
      ptab[++i][0] = 1.0 - 2.0 * a2; ptab[i][1] = a2;
      p->int_coeffs[i] = 0.116786275726378 * 0.5; 
      ptab[++i][0] = a2; ptab[i][1] = 1.0 - 2.0 * a2;
      p->int_coeffs[i] = 0.116786275726378 * 0.5; 
      ptab[++i][0] = aa; ptab[i][1] = bb;
      p->int_coeffs[i] = 0.082851075618374 * 0.5; 
      ptab[++i][0] = aa; ptab[i][1] = 1.0 - aa - bb;
      p->int_coeffs[i] = 0.082851075618374 * 0.5; 
      ptab[++i][0] = bb; ptab[i][1] = aa;
      p->int_coeffs[i] = 0.082851075618374 * 0.5; 
      ptab[++i][0] = bb; ptab[i][1] = 1.0 - aa - bb;
      p->int_coeffs[i] = 0.082851075618374 * 0.5; 
      ptab[++i][0] = 1.0 - aa - bb; ptab[i][1] = aa;
      p->int_coeffs[i] = 0.082851075618374 * 0.5; 
      ptab[++i][0] = 1.0 - aa - bb; ptab[i][1] = bb;
      p->int_coeffs[i] = 0.082851075618374 * 0.5;
      p->repartition[0] = 12;

      // face 0
      double a = 0.5 - 0.4305681557970265;
      double b = 0.5 - 0.1699905217924280;
      double c = 0.5 + 0.1699905217924280;
      double d = 0.5 + 0.4305681557970265;
      double e = 0.326072577431273;
      double f = 0.173927422568727;
      ptab[++i][0] = a; ptab[i][1] = d;
      p->int_coeffs[i] = ::sqrt(2.0) * f; 
      ptab[++i][0] = b; ptab[i][1] = c;
      p->int_coeffs[i] = ::sqrt(2.0) * e; 
      ptab[++i][0] = c; ptab[i][1] = b;
      p->int_coeffs[i] = ::sqrt(2.0) * e; 
      ptab[++i][0] = d; ptab[i][1] = a;
      p->int_coeffs[i] = ::sqrt(2.0) * f; 
      p->repartition[1] = p->repartition[0] + 4;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = a;
      p->int_coeffs[i] = f; 
      ptab[++i][0] = 0.0; ptab[i][1] = b;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = 0.0; ptab[i][1] = c;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = 0.0; ptab[i][1] = d;
      p->int_coeffs[i] = f; 
      p->repartition[2] = p->repartition[1] + 4;
      // face 2
      ptab[++i][0] = a; ptab[i][1] = 0.0;
      p->int_coeffs[i] = f; 
      ptab[++i][0] = b; ptab[i][1] = 0.0;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = c; ptab[i][1] = 0.0;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = d; ptab[i][1] = 0.0;
      p->int_coeffs[i] = f; 
      p->repartition[3] = p->repartition[2] + 4;

      p->pint_points = store_point_tab(ptab);
      if (++i != ptab.size()) DAL_THROW(internal_error, "internal error");
    }
    return p;
  }

  /* ********************************************************************* */
  /*   triangle7 :    Integration on a triangle of order 7 with 13 points  */
  /* ********************************************************************* */

  papprox_integration triangle7_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(25);
      std::fill(ptab.begin(), ptab.end(), base_node(2));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(2);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      size_type i = 0;
      double r1  = 0.0651301029022;
      double r2  = 0.8697397941956;
      double r4  = 0.3128654960049;
      double r5  = 0.6384441885698;
      double r6  = 0.0486903154253;
      double r10 = 0.2603459660790;
      double r11 = 0.4793080678419;
      double r13 = 0.3333333333333;
      double w1  = 0.0533472356088;
      double w4  = 0.0771137608903;
      double w10 = 0.1756152574332;
      double w13 = -0.1495700444677;
      ptab[  i][0] = r1; ptab[i][1] = r1;
      p->int_coeffs[i] = w1 * 0.5; 
      ptab[++i][0] = r2; ptab[i][1] = r1;
      p->int_coeffs[i] = w1 * 0.5; 
      ptab[++i][0] = r1; ptab[i][1] = r2;
      p->int_coeffs[i] = w1 * 0.5; 
      ptab[++i][0] = r4; ptab[i][1] = r6;
      p->int_coeffs[i] = w4 * 0.5; 
      ptab[++i][0] = r5; ptab[i][1] = r4;
      p->int_coeffs[i] = w4 * 0.5; 
      ptab[++i][0] = r6; ptab[i][1] = r5;
      p->int_coeffs[i] = w4 * 0.5; 
      ptab[++i][0] = r5; ptab[i][1] = r6;
      p->int_coeffs[i] = w4 * 0.5; 
      ptab[++i][0] = r4; ptab[i][1] = r5;
      p->int_coeffs[i] = w4 * 0.5; 
      ptab[++i][0] = r6; ptab[i][1] = r4;
      p->int_coeffs[i] = w4 * 0.5; 
      ptab[++i][0] = r10; ptab[i][1] = r10;
      p->int_coeffs[i] = w10 * 0.5; 
      ptab[++i][0] = r11; ptab[i][1] = r10;
      p->int_coeffs[i] = w10 * 0.5; 
      ptab[++i][0] = r10; ptab[i][1] = r11;
      p->int_coeffs[i] = w10 * 0.5; 
      ptab[++i][0] = r13; ptab[i][1] = r13;
      p->int_coeffs[i] = w13 * 0.5; 
      p->repartition[0] = 13;
      // face 0
      double a = 0.5 - 0.4305681557970265;
      double b = 0.5 - 0.1699905217924280;
      double c = 0.5 + 0.1699905217924280;
      double d = 0.5 + 0.4305681557970265;
      double e = 0.326072577431273;
      double f = 0.173927422568727;
      ptab[++i][0] = a; ptab[i][1] = d;
      p->int_coeffs[i] = ::sqrt(2.0) * f; 
      ptab[++i][0] = b; ptab[i][1] = c;
      p->int_coeffs[i] = ::sqrt(2.0) * e; 
      ptab[++i][0] = c; ptab[i][1] = b;
      p->int_coeffs[i] = ::sqrt(2.0) * e; 
      ptab[++i][0] = d; ptab[i][1] = a;
      p->int_coeffs[i] = ::sqrt(2.0) * f; 
      p->repartition[1] = p->repartition[0] + 4;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = a;
      p->int_coeffs[i] = f; 
      ptab[++i][0] = 0.0; ptab[i][1] = b;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = 0.0; ptab[i][1] = c;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = 0.0; ptab[i][1] = d;
      p->int_coeffs[i] = f; 
      p->repartition[2] = p->repartition[1] + 4;
      // face 2
      ptab[++i][0] = a; ptab[i][1] = 0.0;
      p->int_coeffs[i] = f; 
      ptab[++i][0] = b; ptab[i][1] = 0.0;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = c; ptab[i][1] = 0.0;
      p->int_coeffs[i] = e; 
      ptab[++i][0] = d; ptab[i][1] = 0.0;
      p->int_coeffs[i] = f; 
      p->repartition[3] = p->repartition[2] + 4;

      p->pint_points = store_point_tab(ptab);
      if (++i != ptab.size()) DAL_THROW(internal_error, "internal error");

    }
    return p;
  }

  /* ********************************************************************* */
  /* quad2 : Integration on quadrilaterals of order 2 with 3 points        */
  /* ********************************************************************* */

  papprox_integration quad2_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      const int NB_PER_VOL = 3;
      const int NB_PER_FA  = 2;
      const int NB_FA = 4;
      const int dim = 2;
      std::vector<base_node> ptab(NB_PER_VOL + NB_PER_FA * NB_FA);
      base_vector nullpt(dim); nullpt.fill(0);
      std::fill(ptab.begin(), ptab.end(), nullpt);
      p = new _particular_approx;
      p->cvr = parallelepiped_of_reference(dim);
      p->repartition.resize(NB_FA+1);
      p->int_coeffs.resize(ptab.size());
      std::vector<base_node>::iterator itp = ptab.begin();
      std::vector<scalar_type>::iterator itc = p->int_coeffs.begin(); 

      // volume
      *itp++ = base_vector(0.5 + ::sqrt(1.0/6.0), 0.5);
      *itc++ = 1.0 / 3.0;
      *itp++ = base_vector(0.5 - ::sqrt(1.0/24.0), 0.5 + ::sqrt(1.0/8.0));
      *itc++ = 1.0 / 3.0;
      *itp++ = base_vector(0.5 - ::sqrt(1.0/24.0), 0.5 - ::sqrt(1.0/8.0));
      *itc++ = 1.0 / 3.0;
      p->repartition[0] = NB_PER_VOL;

      
      double a = 0.5 + ::sqrt(1.0/3.0) * 0.5;
      double b = 0.5 - ::sqrt(1.0/3.0) * 0.5;
      for (int i = 0; i < NB_FA; ++i) {
	int i1 = (i < 2) ? 1 : 0;
	(*itp)[i1] = a; 
	(*itp)[1 - i1] = scalar_type(1 - (i & 1)); ++itp;
	*itc++ = 0.5;
	(*itp)[i1] = b; 
	(*itp)[1 - i1] = scalar_type(1 - (i & 1)); ++itp;
	*itc++ = 0.5;

	p->repartition[i+1] = p->repartition[i] + NB_PER_FA;
      }

      p->pint_points = store_point_tab(ptab);
      if (itp != ptab.end()) DAL_THROW(internal_error, "internal error");
    }
    return p;
  }

  /* ********************************************************************* */
  /* quad3 : Integration on quadrilaterals of order 3 with 4 points        */
  /* ********************************************************************* */

  papprox_integration quad3_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      const int NB_PER_VOL = 4;
      const int NB_PER_FA  = 2;
      const int NB_FA = 4;
      const int dim = 2;
      std::vector<base_node> ptab(NB_PER_VOL + NB_PER_FA * NB_FA);
      base_vector nullpt(dim); nullpt.fill(0);
      std::fill(ptab.begin(), ptab.end(), nullpt);
      p = new _particular_approx;
      p->cvr = parallelepiped_of_reference(dim);
      p->repartition.resize(NB_FA+1);
      p->int_coeffs.resize(ptab.size());
      std::vector<base_node>::iterator itp = ptab.begin();
      std::vector<scalar_type>::iterator itc = p->int_coeffs.begin(); 

      // volume
      *itp++ = base_vector(0.5 + ::sqrt(2.0/12.0), 0.5);
      *itc++ = 0.25;
      *itp++ = base_vector(0.5 - ::sqrt(2.0/12.0), 0.5);
      *itc++ = 0.25;
      *itp++ = base_vector(0.5, 0.5 + ::sqrt(2.0/12.0));
      *itc++ = 0.25;
      *itp++ = base_vector(0.5, 0.5 - ::sqrt(2.0/12.0));
      *itc++ = 0.25;
      p->repartition[0] = NB_PER_VOL;

      
      double a = 0.5 + ::sqrt(1.0/3.0) * 0.5;
      double b = 0.5 - ::sqrt(1.0/3.0) * 0.5;
      for (int i = 0; i < NB_FA; ++i) {
	int i1 = (i < 2) ? 1 : 0;
	(*itp)[i1] = a; 
	(*itp)[1 - i1] = scalar_type(1 - (i & 1)); ++itp;
	*itc++ = 0.5;
	(*itp)[i1] = b; 
	(*itp)[1 - i1] = scalar_type(1 - (i & 1)); ++itp;
	*itc++ = 0.5;

	p->repartition[i+1] = p->repartition[i] + NB_PER_FA;
      }

      p->pint_points = store_point_tab(ptab);
      if (itp != ptab.end()) DAL_THROW(internal_error, "internal error");
    }
    return p;
  }


  /* ********************************************************************* */
  /* quad5 : Integration on quadrilaterals of order 5 with 7 points        */
  /* ********************************************************************* */

  papprox_integration quad5_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      const int NB_PER_VOL = 7;
      const int NB_PER_FA  = 3;
      const int NB_FA = 4;
      const int dim = 2;
      std::vector<base_node> ptab(NB_PER_VOL + NB_PER_FA * NB_FA);
      base_vector nullpt(dim); nullpt.fill(0);
      std::fill(ptab.begin(), ptab.end(), nullpt);
      p = new _particular_approx;
      p->cvr = parallelepiped_of_reference(dim);
      p->repartition.resize(NB_FA+1);
      p->int_coeffs.resize(ptab.size());
      std::vector<base_node>::iterator itp = ptab.begin();
      std::vector<scalar_type>::iterator itc = p->int_coeffs.begin(); 

      // volume
      *itp++ = base_vector(0.5, 0.5);
      *itc++ = 2.0 / 7.0;
      *itp++ = base_vector(0.5, 0.5 + ::sqrt(14.0/60.0));
      *itc++ = 5.0 / 63.0;
      *itp++ = base_vector(0.5, 0.5 - ::sqrt(14.0/60.0));
      *itc++ = 5.0 / 63.0;
      *itp++ = base_vector(0.5 + ::sqrt(3.0/20.0), 0.5 + ::sqrt(3.0/20.0));
      *itc++ = 5.0 / 36.0;
      *itp++ = base_vector(0.5 - ::sqrt(3.0/20.0), 0.5 + ::sqrt(3.0/20.0));
      *itc++ = 5.0 / 36.0;
      *itp++ = base_vector(0.5 + ::sqrt(3.0/20.0), 0.5 - ::sqrt(3.0/20.0));
      *itc++ = 5.0 / 36.0;
      *itp++ = base_vector(0.5 - ::sqrt(3.0/20.0), 0.5 - ::sqrt(3.0/20.0));
      *itc++ = 5.0 / 36.0;
      p->repartition[0] = NB_PER_VOL;
      
      double a = 0.5 + ::sqrt(3.0/5.0) * 0.5;
      for (int i = 0; i < NB_FA; ++i) {
	int i1 = (i < 2) ? 1 : 0;
	(*itp)[i1] = 0.5; 
	(*itp)[1 - i1] = scalar_type(1 - (i & 1)); ++itp;
	
	*itc++ = 8.0 / 18.0;
	(*itp)[i1] = a; 
	(*itp)[1 - i1] = scalar_type(1 - (i & 1)); ++itp;
	*itc++ = 5.0 / 18.0;
	(*itp)[i1] = 1.0 - a; 
	(*itp)[1 - i1] = scalar_type(1 - (i & 1)); ++itp;
	*itc++ = 5.0 / 18.0;

	p->repartition[i+1] = p->repartition[i] + NB_PER_FA;
      }

      p->pint_points = store_point_tab(ptab);
      if (itp != ptab.end()) DAL_THROW(internal_error, "internal error");
    }
    return p;
  }


  /* ********************************************************************* */
  /*                                                                       */
  /*   Particular integration methods on dimension 3.                      */
  /*                                                                       */
  /* ********************************************************************* */


  /* ********************************************************************* */
  /*  tetrahedron1 : Integration on a tetrahedron of order 1 with 1 point  */
  /* ********************************************************************* */

  papprox_integration tetrahedron1_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(5);
      std::fill(ptab.begin(), ptab.end(), base_node(3));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(3);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      ptab[  i][0] = 0.25; ptab[i][1] = 0.25; ptab[i][2] = 0.25;
      p->int_coeffs[i] = 1.0 / 6.0; 
      p->repartition[0] = 1;
      // face 0
      ptab[++i][0] = 1.0 / 3.0; ptab[i][1] = 1.0 / 3.0; ptab[i][2] = 1.0 / 3.0;
      p->int_coeffs[i] = 0.5 * ::sqrt(3.0); 
      p->repartition[1] = p->repartition[0] + 1;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = 0.5; ptab[i][2] = 0.5;
      p->int_coeffs[i] = 0.5; 
      p->repartition[2] = p->repartition[1] + 1;
      // face 2
      ptab[++i][0] = 0.5; ptab[i][1] = 0.0; ptab[i][2] = 0.5;
      p->int_coeffs[i] = 0.5; 
      p->repartition[3] = p->repartition[2] + 1;
      // face 3
      ptab[++i][0] = 0.5; ptab[i][1] = 0.5; ptab[i][2] = 0.0;
      p->int_coeffs[i] = 0.5; 
      p->repartition[4] = p->repartition[3] + 1;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }

  /* ********************************************************************* */
  /*  tetrahedron2 : Integration on a tetrahedron of order 2 with 4 points */
  /* ********************************************************************* */

  papprox_integration tetrahedron2_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      std::vector<base_node> ptab(16);
      std::fill(ptab.begin(), ptab.end(), base_node(3));
      p = new _particular_approx;
      p->cvr = simplex_of_reference(3);
      p->repartition.resize(p->cvr->structure()->nb_faces()+1);
      p->int_coeffs.resize(ptab.size());
      // volume
      int i = 0;
      double c = 0.13819660112501052;
      ptab[  i][0] = c; ptab[i][1] = c; ptab[i][2] = c;
      p->int_coeffs[i] = 1.0 / 24.0; 
      ptab[++i][0] = 1 - 3.0 * c; ptab[i][1] = c; ptab[i][2] = c;
      p->int_coeffs[i] = 1.0 / 24.0; 
      ptab[++i][0] = c; ptab[i][1] = 1 - 3.0 * c; ptab[i][2] = c;
      p->int_coeffs[i] = 1.0 / 24.0; 
      ptab[++i][0] = c; ptab[i][1] = c; ptab[i][2] = 1 - 3.0 * c;
      p->int_coeffs[i] = 1.0 / 24.0; 
      p->repartition[0] = 4;

      // face 0
      double a = 1.0 / 6.0;
      double b = 2.0 / 3.0;
      ptab[++i][0] = 1.0 - a - a; ptab[i][1] = a; ptab[i][2] = a;
      p->int_coeffs[i] = ::sqrt(3.0) / 6.0; 
      ptab[++i][0] = 1.0 - b - a; ptab[i][1] = b; ptab[i][2] = a;
      p->int_coeffs[i] = ::sqrt(3.0) / 6.0; 
      ptab[++i][0] = 1.0 - a - b; ptab[i][1] = a; ptab[i][2] = b;
      p->int_coeffs[i] = ::sqrt(3.0) / 6.0; 
      p->repartition[1] = p->repartition[0] + 3;
      // face 1
      ptab[++i][0] = 0.0; ptab[i][1] = a; ptab[i][2] = a;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = 0.0; ptab[i][1] = b; ptab[i][2] = a;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = 0.0; ptab[i][1] = a; ptab[i][2] = b;
      p->int_coeffs[i] = 1.0 / 6.0; 
      p->repartition[2] = p->repartition[1] + 3;
      // face 2
      ptab[++i][0] = a; ptab[i][1] = 0.0; ptab[i][2] = a;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = b; ptab[i][1] = 0.0; ptab[i][2] = a;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = a; ptab[i][1] = 0.0; ptab[i][2] = b;
      p->int_coeffs[i] = 1.0 / 6.0; 
      p->repartition[3] = p->repartition[2] + 3;
      // face 3
      ptab[++i][0] = a; ptab[i][1] = a; ptab[i][2] = 0.0;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = b; ptab[i][1] = a; ptab[i][2] = 0.0;
      p->int_coeffs[i] = 1.0 / 6.0; 
      ptab[++i][0] = a; ptab[i][1] = b; ptab[i][2] = 0.0;
      p->int_coeffs[i] = 1.0 / 6.0; 
      p->repartition[4] = p->repartition[3] + 3;

      p->pint_points = store_point_tab(ptab);
      
    }
    return p;
  }


  /* ********************************************************************* */
  /*  tetrahedron3 : Integration on a tetrahedron of order 3 with 5 points */
  /* ********************************************************************* */

  papprox_integration tetrahedron3_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      const int NB_PER_VOL = 5;
      const int NB_PER_FA  = 4;
      const int NB_FA = 4;
      const int dim = 3;
      std::vector<base_node> ptab(NB_PER_VOL + NB_PER_FA * NB_FA);
      base_vector nullpt(dim); nullpt.fill(0);
      std::fill(ptab.begin(), ptab.end(), nullpt);
      p = new _particular_approx;
      p->cvr = simplex_of_reference(dim);
      p->repartition.resize(NB_FA+1);
      p->int_coeffs.resize(ptab.size());
      std::vector<base_node>::iterator itp = ptab.begin();
      std::vector<scalar_type>::iterator itc = p->int_coeffs.begin(); 
      // volume
      *itp++ = base_vector(0.25, 0.25, 0.25);
      *itc++ = - 4.0 / 30.0;
      *itp++ = base_vector(1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0);
      *itc++ = 9.0 / 120.0; 
      *itp++ = base_vector(1.0 / 2.0, 1.0 / 6.0, 1.0 / 6.0);
      *itc++ = 9.0 / 120.0; 
      *itp++ = base_vector(1.0 / 6.0, 1.0 / 2.0, 1.0 / 6.0);
      *itc++ = 9.0 / 120.0; 
      *itp++ = base_vector(1.0 / 6.0, 1.0 / 6.0, 1.0 / 2.0);
      *itc++ = 9.0 / 120.0; 
      p->repartition[0] = NB_PER_VOL;

      double a = 1.0 / 3.0;
      double b = 1.0 / 5.0;
      double c = 3.0 / 5.0;
      double surf = ::sqrt(3.0) * 0.5;
      for (int i = 0; i < NB_FA; ++i) {
	int i1 = (i < 2) ? 1 : 0;
	int i2 = (i < 3) ? 2 : 1;
	(*itp)[i1] = a; (*itp)[i2] = a;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = -surf * 9.0 / 16.0; ++itp;
	(*itp)[i1] = b; (*itp)[i2] = b;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 25.0 / 48.0; ++itp;
	(*itp)[i1] = c; (*itp)[i2] = b;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 25.0 / 48.0; ++itp;
	(*itp)[i1] = b; (*itp)[i2] = c;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 25.0 / 48.0; ++itp;

	p->repartition[i+1] = p->repartition[i] + NB_PER_FA;
	surf = 0.5;
      }

      p->pint_points = store_point_tab(ptab);
      if (itp != ptab.end()) DAL_THROW(internal_error, "internal error");
    }
    return p;
  }

  /* ********************************************************************* */
  /* tetrahedron5 : Integration on a tetrahedron of order 5 with 15 points */
  /* ********************************************************************* */

  papprox_integration tetrahedron5_approx_integration(void) {
    static _particular_approx *p = NULL;
    if (p == NULL)
    {
      const int NB_PER_VOL = 15;
      const int NB_PER_FA  = 7;
      const int NB_FA = 4;
      const int dim = 3;
      std::vector<base_node> ptab(NB_PER_VOL + NB_PER_FA * NB_FA);
      base_vector nullpt(dim); nullpt.fill(0);
      std::fill(ptab.begin(), ptab.end(), nullpt);
      p = new _particular_approx;
      p->cvr = simplex_of_reference(dim);
      p->repartition.resize(NB_FA+1);
      p->int_coeffs.resize(ptab.size());
      std::vector<base_node>::iterator itp = ptab.begin();
      std::vector<scalar_type>::iterator itc = p->int_coeffs.begin(); 

      // volume
      double b1 = 0.31979362782962991;
      double b2 = 0.091971078052723033;
      double c1 = 0.040619116511110275;
      double c2 = 0.7240867658418309;
      double d  = 0.056350832689629156;
      double e  = 0.44364916731037084;
      double w0 = 0.019753086419753086;
      double w1 = 0.011511367871045398;
      double w2 = 0.01198951396316977;
      double w3 = 0.008818342151675485;
      *itp++ = base_vector(0.25, 0.25, 0.25);
      *itc++ = w0;
      *itp++ = base_vector(b1, b1, b1);
      *itc++ = w1;
      *itp++ = base_vector(c1, b1, b1);
      *itc++ = w1;
      *itp++ = base_vector(b1, c1, b1);
      *itc++ = w1;
      *itp++ = base_vector(b1, b1, c1);
      *itc++ = w1;
      *itp++ = base_vector(b2, b2, b2);
      *itc++ = w2;
      *itp++ = base_vector(c2, b2, b2);
      *itc++ = w2;
      *itp++ = base_vector(b2, c2, b2);
      *itc++ = w2;
      *itp++ = base_vector(b2, b2, c2);
      *itc++ = w2;
      *itp++ = base_vector(d, d, e);
      *itc++ = w3;
      *itp++ = base_vector(d, e, d);
      *itc++ = w3;
      *itp++ = base_vector(e, d, d);
      *itc++ = w3;
      *itp++ = base_vector(d, e, e);
      *itc++ = w3; 
      *itp++ = base_vector(e, e, d);
      *itc++ = w3;
      *itp++ = base_vector(e, d, e);
      *itc++ = w3;
      p->repartition[0] = NB_PER_VOL;

      
      double aa1 = 0.10128650732345634;
      double aa2 = 0.47014206410511509;
      double surf = ::sqrt(3.0) * 0.5;
      for (int i = 0; i < NB_FA; ++i) {
	int i1 = (i < 2) ? 1 : 0;
	int i2 = (i < 3) ? 2 : 1;
	(*itp)[i1] = 1.0 / 3.0; (*itp)[i2] = 1.0 / 3.0;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 9.0 / 40.0; ++itp;
	(*itp)[i1] = aa1; (*itp)[i2] = aa1;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 0.12593918054482715; ++itp;
	(*itp)[i1] = aa1; (*itp)[i2] = 1.0 - 2.0 * aa1;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 0.12593918054482715; ++itp;
	(*itp)[i1] = 1.0 - 2.0 * aa1; (*itp)[i2] = aa1;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 0.12593918054482715; ++itp;
	(*itp)[i1] = aa2; (*itp)[i2] = aa2;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 0.13239415278850618; ++itp;
	(*itp)[i1] = aa2; (*itp)[i2] = 1.0 - 2.0 * aa2;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 0.13239415278850618; ++itp;
	(*itp)[i1] = 1.0 - 2.0 * aa2; (*itp)[i2] = aa2;
	if (i == 0) (*itp)[0] = 1.0 - (*itp)[1] - (*itp)[2];
	*itc++ = surf * 0.13239415278850618; ++itp;

	p->repartition[i+1] = p->repartition[i] + NB_PER_FA;
	surf = 0.5;
      }

      p->pint_points = store_point_tab(ptab);
      if (itp != ptab.end()) DAL_THROW(internal_error, "internal error");
    }
    return p;
  }


}  /* end of namespace bgeot.                                            */

