/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_poly_integration.C : exact integration of polynomial   */
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


#include <bgeot_integration.h>
#include <ftool_naming.h>

namespace bgeot
{

  typedef ftool::naming_system<integration_method>::param_list im_param_list;

  typedef polynomial<scalar_type> base_poly;

  scalar_type poly_integration::int_poly(const base_poly &P) const
  {
    scalar_type res = 0.0;
    if (P.size() > int_monomials.size())
    {
      std::vector<scalar_type>
	*hum = (std::vector<scalar_type> *)(&int_monomials);
      size_type i = P.size(), j = int_monomials.size();
      hum->resize(i);
      power_index mi(P.dim()); mi[P.dim()-1] = P.degree();
      for (size_type k = i; k > j; --k, --mi)
	(*hum)[k-1] = int_monomial(mi);
    }
    polynomial<scalar_type>::const_iterator it = P.begin(), ite = P.end();
    std::vector<scalar_type>::const_iterator itb = int_monomials.begin();
    for ( ; it != ite; ++it, ++itb) res += (*it) * (*itb);
    return res;
  }

  scalar_type
    poly_integration::int_poly_on_face(const base_poly &P, short_type f) const
  {
    scalar_type res = 0.0;
    std::vector<scalar_type>
	*hum = (std::vector<scalar_type> *)(&(int_face_monomials[f]));
    if (P.size() > hum->size())
    {
      size_type i = P.size(), j = hum->size();
      hum->resize(i);
      power_index mi(P.dim()); mi[P.dim()-1] = P.degree();
      for (size_type k = i; k > j; --k, --mi)
	(*hum)[k-1] = int_monomial_on_face(mi, f);
    }
    polynomial<scalar_type>::const_iterator it = P.begin(), ite = P.end();
    std::vector<scalar_type>::const_iterator itb = hum->begin();
    for ( ; it != ite; ++it, ++itb) res += (*it) * (*itb);
    return res;
  }

  /* ******************************************************************** */
  /* integration on simplex                                               */
  /* ******************************************************************** */

  struct _simplex_poly_integration : public poly_integration
  {
    scalar_type int_monomial(const power_index &power) const
    {
      scalar_type res = 1.0;
      short_type fa = 1;
      power_index::const_iterator itm = power.begin(), itme = power.end();
      for ( ; itm != itme; ++itm)
	for (int k = 1; k <= *itm; ++k, ++fa)
	  res *= scalar_type(k) / scalar_type(fa);
	
      for (int k = 0; k < cvs->dim(); k++) { res /= scalar_type(fa); fa++; }
      return res;
    }

    scalar_type int_monomial_on_face(const power_index &power, 
					       short_type f) const
    {
      scalar_type res = 0.0;
    
      if (f == 0 || power[f-1] == 0.0)
      {
	res = (f == 0) ? sqrt(scalar_type(cvs->dim())) : 1.0;
	short_type fa = 1;
	power_index::const_iterator itm = power.begin(), itme = power.end();
	for ( ; itm != itme; ++itm)
	  for (int k = 1; k <= *itm; ++k, ++fa)
	    res *= scalar_type(k) / scalar_type(fa);
	
	for (int k = 1; k < cvs->dim(); k++) { res /= scalar_type(fa); fa++; }
      }
      return res;
    }

    _simplex_poly_integration(pconvex_structure c)
    {
      cvs = c;
      int_face_monomials.resize(c->nb_faces());
    }
  };

  static pintegration_method exact_simplex(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    return new integration_method
      (new _simplex_poly_integration(simplex_structure(n)));
  }

  /* ******************************************************************** */
  /* integration on direct product of convex structures                   */
  /* ******************************************************************** */

  struct _plyint_mul_structure : public poly_integration
  {
    ppoly_integration cv1, cv2;

    scalar_type int_monomial(const power_index &power) const
    {
      power_index mi1(cv1->dim()), mi2(cv2->dim());
      std::copy(power.begin(), power.begin() + cv1->dim(), mi1.begin());
      std::copy(power.begin() + cv1->dim(), power.end(), mi2.begin());
      return cv1->int_monomial(mi1) * cv2->int_monomial(mi2);
    }

    scalar_type int_monomial_on_face(const power_index &power, 
				     short_type f) const
    {
      power_index mi1(cv1->dim()), mi2(cv2->dim());
      std::copy(power.begin(), power.begin() + cv1->dim(), mi1.begin());
      std::copy(power.begin() + cv1->dim(), power.end(), mi2.begin());
      short_type nfx = cv1->structure()->nb_faces();
      if (f < nfx)
	return cv1->int_monomial_on_face(mi1,f) * cv2->int_monomial(mi2);
      else
	return cv1->int_monomial(mi1) * cv2->int_monomial_on_face(mi2, f-nfx);
    }

    _plyint_mul_structure(ppoly_integration a, ppoly_integration b)
    {
      cv1 = a; cv2 = b;
      cvs = convex_product_structure(cv1->structure(), cv2->structure());
      int_face_monomials.resize(cvs->nb_faces());
    }
  };

  static pintegration_method product_exact(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pintegration_method a = params[0].method();
    pintegration_method b = params[1].method();
    if (!(a->is_ppi && b->is_ppi))
      DAL_THROW(failure_error, "Bad parameters");
    return new integration_method(new _plyint_mul_structure(a->method.ppi,
							     b->method.ppi));
  }

  /* ******************************************************************** */
  /* integration on parallelepiped.                                       */
  /* ******************************************************************** */

  static pintegration_method exact_parallelepiped(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");

    _STRINGSTREAM name;
    if (n == 1)
      name << "IM_EXACT_SIMPLEX(1)" << ends;
    else 
      name << "IM_PRODUCT(IM_EXACT_PARALLELEPIPED(" << n-1
	   << "),IM_EXACT_SIMPLEX(1)))" << ends;
    return int_method_descriptor(name.str());
  }

  static pintegration_method exact_prism(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 1 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");

    _STRINGSTREAM name;
    name << "IM_PRODUCT(IM_EXACT_SIMPLEX(" << n-1
	 << "),IM_EXACT_SIMPLEX(1)" << ends;
    return int_method_descriptor(name.str());
  }

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

  static pintegration_method gauss1d(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n < 0 || n >= 32000 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    if (n & 1) {
      _STRINGSTREAM name;
      name << "IM_GAUSS1D(" << n-1 << ")" << ends;
      return int_method_descriptor(name.str());
    }
    else
      return new integration_method(new _gauss_approx_integration(n/2 + 1));
  }

  /* ********************************************************************* */
  /* integration on simplexes                                              */
  /* ********************************************************************* */

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

    _Newton_Cotes_approx_integration(dim_type nc, short_type k)
    {
      cvr = simplex_of_reference(nc);
      size_type R = alpha(nc,k);
      size_type R2 = (nc > 0) ? alpha(nc-1,k) : 0;
      base_poly P;
      _STRINGSTREAM name;
      name << "IM_EXACT_SIMPLEX(" << int(nc) << ")" << ends;
      ppoly_integration ppi 
	= int_method_descriptor(name.str())->method.ppi;
      std::vector<size_type> fa(nc+1);
      stored_point_tab int_points;
      int_points.resize(R + (nc+1) * R2);
      int_coeffs.resize(R + (nc+1) * R2);
      repartition.resize(nc+2);
      repartition[0] = R;
      for (short_type f = 0; f <= nc; ++f)
      {
	fa[f] = repartition[f];
	repartition[f+1] = repartition[f] + R2;
      }

      size_type sum = 0, l;
      base_node c(nc); 
      c *= scalar_type(0.0);
      if (k == 0) c.fill(1.0 / scalar_type(nc+1));

      for (size_type r = 0; r < R; ++r)
      {
	int_points[r] = c;
	calc_base_func(P, k, c);
	int_coeffs[r] = ppi->int_poly(P);

	for (short_type f = 1; f <= nc; ++f)
	  if (c[f-1] == 0.0)
	  {
	    int_points[fa[f]] = c;
	    int_coeffs[fa[f]] = ppi->int_poly_on_face(P, f);
	    (fa[f])++;
	  }
	if (k == 0)
	{
	  for (short_type f = 0; f <= nc; ++f)
	  {
	    int_points[fa[f]].resize(dim());
	    int_points[fa[f]].fill(1.0 / scalar_type(nc));
	    if (f > 0) int_points[fa[f]][f-1] = 0.0;
	    int_coeffs[fa[f]] = ppi->int_poly_on_face(P, f);
	  }
	}
	else
	if (sum == k)
	{
	  int_points[fa[0]] = c;
	  int_coeffs[fa[0]] = ppi->int_poly_on_face(P, 0);
	  (fa[0])++;
	}  

	if (k != 0) {
	  l = 0; c[l] += 1.0 / scalar_type(k); sum++;
	  while (sum > k) {
	    sum -= int(floor(0.5+(c[l] * k)));
	    c[l] = 0.0; l++; if (l == nc) break;
	    c[l] += 1.0 / scalar_type(k); sum++;
	  }
	}
      }
      pint_points = store_point_tab(int_points);
    }
  };

  static pintegration_method Newton_Cotes(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
          "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");
    return new integration_method(new _Newton_Cotes_approx_integration(n, k));
  }

  /* ********************************************************************* */
  /* integration on direct product of convex structures                    */
  /* ********************************************************************* */

  struct a_int_pro_integration : public approx_integration
  {

    a_int_pro_integration(papprox_integration a, papprox_integration b)
    {
      cvr = convex_ref_product(a->ref_convex(), b->ref_convex());
      size_type n1 = a->nb_points_on_convex();
      size_type n2 = b->nb_points_on_convex();
      stored_point_tab int_points;
      int_points.resize(n1 * n2);
      int_coeffs.resize(n1 * n2);
      repartition.resize(cvr->structure()->nb_faces()+1);
      repartition[0] = n1 * n2;
      for (size_type i1 = 0; i1 < n1; ++i1)
	for (size_type i2 = 0; i2 < n2; ++i2)
	{
	  int_coeffs[i1 + i2 * n1] = a->coeff(i1) * b->coeff(i2);
	  int_points[i1 + i2 * n1].resize(dim());
	  std::copy(a->point(i1).begin(), a->point(i1).end(),
		    int_points[i1 + i2 * n1].begin());
	  std::copy(b->point(i2).begin(), b->point(i2).end(),
		    int_points[i1 + i2 * n1].begin() + a->dim());
	}
      short_type f = 0;
      for (short_type f1 = 0; f1 < a->structure()->nb_faces(); ++f1, ++f)
      {
	n1 = a->nb_points_on_face(f1);
	n2 = b->nb_points_on_convex();
	size_type w = repartition[f];
	repartition[f+1] = w + n1 * n2;
	int_points.resize(repartition[f+1]);
	int_coeffs.resize(repartition[f+1]);
	for (size_type i1 = 0; i1 < n1; ++i1)
	  for (size_type i2 = 0; i2 < n2; ++i2)
	  {
	    int_coeffs[w + i1 + i2 * n1] = a->coeff_on_face(f1, i1)
	                                 * b->coeff(i2);
	    int_points[w + i1 + i2 * n1].resize(dim());
	    std::copy(a->point_on_face(f1, i1).begin(),
		      a->point_on_face(f1, i1).end(),
		      int_points[w + i1 + i2 * n1].begin());
	    std::copy(b->point(i2).begin(), b->point(i2).end(),
		      int_points[w + i1 + i2 * n1].begin() + a->dim());
	  }
      }
      for (short_type f2 = 0; f2 < b->structure()->nb_faces(); ++f2, ++f)
      {
	n1 = a->nb_points_on_convex();
	n2 = b->nb_points_on_face(f2);
	size_type w = repartition[f];
	repartition[f+1] = w + n1 * n2;
	int_points.resize(repartition[f+1]);
	int_coeffs.resize(repartition[f+1]);
	for (size_type i1 = 0; i1 < n1; ++i1)
	  for (size_type i2 = 0; i2 < n2; ++i2)
	  {
	    int_coeffs[w + i1 + i2 * n1] = a->coeff(i1)
	                                 * b->coeff_on_face(f2, i2);
	    int_points[w + i1 + i2 * n1].resize(dim());
	    std::copy(a->point(i1).begin(), a->point(i1).end(),
		      int_points[w + i1 + i2 * n1].begin());
	    std::copy(b->point_on_face(f2, i2).begin(),
		      b->point_on_face(f2, i2).end(),
		      int_points[w + i1 + i2 * n1].begin() + a->dim());
	  }
      }
      pint_points = store_point_tab(int_points);
    }
  };

  static pintegration_method product_approx(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pintegration_method a = params[0].method();
    pintegration_method b = params[1].method();
    if (a->is_ppi || b->is_ppi)
      DAL_THROW(failure_error, "Bad parameters");
    return new integration_method(new a_int_pro_integration(a->method.pai,
							    b->method.pai));
  }

  static pintegration_method product_which(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
	  "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 1 || params[1].type() != 1)
      DAL_THROW(failure_error, "Bad type of parameters");
    pintegration_method a = params[0].method();
    pintegration_method b = params[1].method();
    if (a->is_ppi || b->is_ppi) return product_exact(params);
    else return product_approx(params);
  }


  /* ********************************************************************* */
  /* integration on parallelepiped with Newton Cotes formulae              */
  /* ********************************************************************* */

  static pintegration_method Newton_Cotes_para(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
          "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");

    _STRINGSTREAM name;
    if (n == 1)
      name << "IM_NC(1," << k << ")" << ends;
    else 
      name << "IM_PRODUCT(IM_NC_PARALLELEPIPED(" << n-1 << "," << k
	   << "),IM_NC(1," << k << "))" << ends;
    return int_method_descriptor(name.str());
  }

  static pintegration_method Newton_Cotes_prism(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
          "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 1 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");

    _STRINGSTREAM name;
     name << "IM_PRODUCT(IM_NC(" << n-1 << "," << k << "),IM_NC(1,"
	 << k << "))" << ends;
    return int_method_descriptor(name.str());
  }

  /* ********************************************************************* */
  /* integration on parallelepiped with Gauss formulae                     */
  /* ********************************************************************* */

  static pintegration_method Gauss_paramul(im_param_list &params) {
    if (params.size() != 2)
      DAL_THROW(failure_error, 
          "Bad number of parameters : " << params.size() << " should be 2.");
    if (params[0].type() != 0 || params[1].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    int k = int(::floor(params[1].num() + 0.01));
    if (n <= 0 || n >= 100 || k < 0 || k > 150 ||
	double(n) != params[0].num() || double(k) != params[1].num())
      DAL_THROW(failure_error, "Bad parameters");

    _STRINGSTREAM name;
    if (n == 1)
      name << "IM_GAUSS1D(" << k << ")" << ends;
    else 
      name << "IM_PRODUCT(IM_GAUSS_PARALLELEPIPED(" << n-1 << "," << k
	   << "),IM_GAUSS1D(" << k << "))" << ends;
    return int_method_descriptor(name.str());
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

  static pintegration_method approx_triangle(im_param_list &params) {
    if (params.size() <= 0 || params.size()>=2)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    int variante = 0;
    if (params.size() == 2) {
      variante = int(::floor(params[1].num() + 0.01));
      if (variante != 2 || n != 2 || 
	  double(variante) != params[1].num())
	DAL_THROW(failure_error, "Bad parameters");
    }
    switch (n) {
    case 1: return new integration_method(triangle1_approx_integration());
      break;
    case 2:
      if (variante == 0)
	return new integration_method(triangle2_approx_integration());
      else 
	return new integration_method(triangle2bis_approx_integration());
      break;
    case 3: return new integration_method(triangle3_approx_integration());
      break;
    case 4: return new integration_method(triangle4_approx_integration());
      break;
    case 5: return new integration_method(triangle5_approx_integration());
      break;
    case 6: return new integration_method(triangle6_approx_integration());
      break;
    case 7: return new integration_method(triangle7_approx_integration());
      break;
    default : DAL_THROW(failure_error, "Bad parameters");
    }
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

  static pintegration_method approx_quad(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    switch (n) {
    case 2: return new integration_method(quad2_approx_integration());
      break;
    case 3: return new integration_method(quad3_approx_integration());
      break;
    case 5: return new integration_method(quad5_approx_integration());
      break;
    default : DAL_THROW(failure_error, "Bad parameters");
    }
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

  static pintegration_method approx_tetra(im_param_list &params) {
    if (params.size() != 1)
      DAL_THROW(failure_error, 
	   "Bad number of parameters : " << params.size() << " should be 1.");
    if (params[0].type() != 0)
      DAL_THROW(failure_error, "Bad type of parameters");
    int n = int(::floor(params[0].num() + 0.01));
    if (n <= 0 || n >= 100 || double(n) != params[0].num())
      DAL_THROW(failure_error, "Bad parameters");
    switch (n) {
    case 1: return new integration_method(tetrahedron1_approx_integration());
      break;
    case 2: return new integration_method(tetrahedron2_approx_integration());
      break;
    case 3: return new integration_method(tetrahedron3_approx_integration());
      break;
    case 5: return new integration_method(tetrahedron5_approx_integration());
      break;
    default : DAL_THROW(failure_error, "Bad parameters");
    }
  }

  /* ******************************************************************** */
  /*    Naming system                                                     */
  /* ******************************************************************** */

  static ftool::naming_system<integration_method> *_im_naming_system = 0;
  
  static void init_im_naming_system(void) {
    _im_naming_system = new ftool::naming_system<integration_method>("IM");
    _im_naming_system->add_suffix("EXACT_SIMPLEX", exact_simplex);
    _im_naming_system->add_suffix("PRODUCT", product_which);
    _im_naming_system->add_suffix("EXACT_PARALLELEPIPED",exact_parallelepiped);
    _im_naming_system->add_suffix("EXACT_PRISM", exact_prism);
    _im_naming_system->add_suffix("GAUSS1D", gauss1d);
    _im_naming_system->add_suffix("NC", Newton_Cotes);
    _im_naming_system->add_suffix("NC_PARALLELEPIPED", Newton_Cotes_para);
    _im_naming_system->add_suffix("NC_PRISM", Newton_Cotes_prism);
    _im_naming_system->add_suffix("GAUSS_PARALLELEPIPED", Gauss_paramul);
    _im_naming_system->add_suffix("TRIANGLE", approx_triangle);
    _im_naming_system->add_suffix("QUAD", approx_quad);
    _im_naming_system->add_suffix("TETRAHEDRON", approx_tetra);
  }
  
  pintegration_method int_method_descriptor(std::string name) {
    if (_im_naming_system == 0) init_im_naming_system();
    size_type i = 0;
    return _im_naming_system->method(name, i);
  }

  std::string name_of_int_method(pintegration_method p) {
    if (_im_naming_system == 0) init_im_naming_system();
    return _im_naming_system->name_of_method(p);
  }

  /* Fonctions pour la ref. directe.                                     */
  
  pintegration_method exact_simplex_im(size_type n) {
    static pintegration_method pim = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      _STRINGSTREAM name;
      name << "IM_EXACT_SIMPLEX(" << n << ")" << ends;
      pim = int_method_descriptor(name.str());
      d = n;
    }
    return pim;
  }

  pintegration_method exact_parallelepiped_im(size_type n) {
    static pintegration_method pim = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      _STRINGSTREAM name;
      name << "IM_EXACT_PARALLELEPIPED(" << n << ")" << ends;
      pim = int_method_descriptor(name.str());
      d = n;
    }
    return pim;
  }

  pintegration_method exact_prism_im(size_type n) {
    static pintegration_method pim = 0;
    static size_type d = size_type(-2);
    if (d != n) {
      _STRINGSTREAM name;
      name << "IM_EXACT_PRISM(" << n << ")" << ends;
      pim = int_method_descriptor(name.str());
      d = n;
    }
    return pim;
  }


}  /* end of namespace bgeot.                                            */

