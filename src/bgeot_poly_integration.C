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


#include <bgeot_poly_integration.h>

namespace bgeot
{

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


  ppoly_integration simplex_poly_integration(dim_type nc)
  {
    static dal::dynamic_array<ppoly_integration> *_simplex;
    static dal::bit_vector *_simplex_exists;
    static bool isinit = false;
    if (!isinit) {
      _simplex = new dal::dynamic_array<ppoly_integration>();
      _simplex_exists = new dal::bit_vector();
      isinit = true;
    }

    if (!((*_simplex_exists)[nc]))
    {
      (*_simplex)[nc] = new _simplex_poly_integration(simplex_structure(nc));
      (*_simplex_exists)[nc] = true;
    }
    return (*_simplex)[nc];
  }

  /* ******************************************************************** */
  /* integration on direct product of convex structures                   */
  /* ******************************************************************** */


  struct _plyint_mul_light
  {
    ppoly_integration cv1, cv2;
    bool operator < (const _plyint_mul_light &ls) const
    {
      if (cv1 < ls.cv1) return true; if (cv1 > ls.cv1) return false; 
      if (cv2 < ls.cv2) return true; return false;
    }
    _plyint_mul_light(ppoly_integration a, ppoly_integration b)
    { cv1 = a; cv2 = b; }
    _plyint_mul_light(void) { }
  };


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

    _plyint_mul_structure(const _plyint_mul_light &ls)
    {
      cv1 = ls.cv1; cv2 = ls.cv2;
      cvs = convex_product_structure(cv1->structure(), cv2->structure());
      int_face_monomials.resize(cvs->nb_faces());
    }
  };

  
  ppoly_integration convex_product_poly_integration(ppoly_integration a,
						    ppoly_integration b)
  {
    static dal::FONC_TABLE<_plyint_mul_light, _plyint_mul_structure> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_plyint_mul_light, _plyint_mul_structure>();
      isinit = true;
    }
    return tab->add(_plyint_mul_light(a, b));
  }

  /* ******************************************************************** */
  /* integration on parallelepiped.                                       */
  /* ******************************************************************** */


  ppoly_integration parallelepiped_poly_integration(dim_type nc)
  {
    static dal::dynamic_array<ppoly_integration> *tab;
    static int _nb_parallelepiped = -1;

    static bool isinit = false;
    if (!isinit) {
      tab = new dal::dynamic_array<ppoly_integration>();
      isinit = true;
    }

    if (nc > _nb_parallelepiped)
    {
      if (_nb_parallelepiped < 0)
      {
	(*tab)[0] = simplex_poly_integration(0);
	(*tab)[1] = simplex_poly_integration(1);
      }
      for (int i = 1; i < nc; i++)
	(*tab)[i+1] = convex_product_poly_integration((*tab)[i],
					 simplex_poly_integration(1));
      _nb_parallelepiped = nc;
    }
    return (*tab)[nc];
  }

}  /* end of namespace bgeot.                                            */

