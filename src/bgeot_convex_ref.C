/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_convex_ref.C :  convexes of reference                  */
/*     									   */
/* Date : Septembre 28, 2001.                                              */
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

#include <bgeot_convex_ref.h>

namespace bgeot
{

  /* simplexes.                                                            */

  struct _K_simplex_ref_light
  {
    dim_type N; short_type K;
    bool operator < (const _K_simplex_ref_light &l) const
    {
      if (N < l.N) return true; if (N > l.N) return false; 
      if (K < l.K) return true; return false;
    }
    _K_simplex_ref_light(dim_type n, short_type k) { N = n; K = k; }
    _K_simplex_ref_light(void) { }
  };

  class _K_simplex_of_ref : public convex_of_reference 
  {
    public :
      scalar_type is_in(const base_node &pt) const
      { // return a negative or null number if pt is in the convex
        assert(pt.size() == cvs->dim());
	scalar_type e = -1.0, r = 0.0;
	base_node::const_iterator it = pt.begin(), ite = pt.end();
	for (; it != ite; e += *it, ++it) r = std::min(r, *it);
	return std::max(-r, e);
      }
      scalar_type is_in_face(short_type f, const base_node &pt) const
      { // return a null number if pt is in the face of the convex
	assert(pt.size() == cvs->dim());
	if (f > 0) return dal::abs(pt[f-1]);
	scalar_type e = -1.0;
	base_node::const_iterator it = pt.begin(), ite = pt.end();
	for (; it != ite; e += *it, ++it);
	return dal::abs(e);
      }
      _K_simplex_of_ref(const _K_simplex_ref_light &ls)
      {
	cvs = simplex_structure(ls.N, ls.K);
	size_type R = cvs->nb_points();
	points().resize(R);
	_normals.resize(ls.N+1);
	base_node null(ls.N); null.fill(0.0);
	std::fill(_normals.begin(), _normals.end(), null);
	std::fill(points().begin(), points().end(), null);
	for (size_type i = 1; i <= ls.N; ++i)
	  _normals[i][i-1] = -1.0;
	std::fill(_normals[0].begin(), _normals[0].end(),
		  scalar_type(1.0)/sqrt(scalar_type(ls.N)));
	base_node c(ls.N);  c.fill(0.0);
	
	if (ls.K == 0)
	{
	  c.fill(1.0/(ls.N+1));
	  points()[0] = c;
	}
	else
	{
	  size_type sum = 0, l;
	  for (size_type r = 0; r < R; ++r)
	  {
	    points()[r] = c;
	    l = 0; c[l] += 1.0 / scalar_type(ls.K); sum++;
	    while (sum > ls.K)
	    {
	      sum -= int(floor(0.5+(c[l] * ls.K)));
	      c[l] = 0.0; l++; if (l == ls.N) break;
	      c[l] += 1.0 / scalar_type(ls.K); sum++;
	    }
	  }
	}
      }
  };


  pconvex_ref simplex_of_reference(dim_type nc, short_type k)
  {
    static dal::FONC_TABLE<_K_simplex_ref_light, _K_simplex_of_ref> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_K_simplex_ref_light, _K_simplex_of_ref>();
      isinit = true;
    }
    return tab->add(_K_simplex_ref_light(nc, k));
  }


  /* products.                                                             */

  struct _product_ref_light
  {
    pconvex_ref cvr1, cvr2;
    bool operator < (const _product_ref_light &ls) const
    {
      if (cvr1 < ls.cvr1) return true; if (cvr1 > ls.cvr1) return false; 
      if (cvr2 < ls.cvr2) return true; return false;
    }
    _product_ref_light(pconvex_ref a, pconvex_ref b)
    { cvr1 = a; cvr2 = b; }
    _product_ref_light(void) { }
  };

  struct _product_ref : public convex_of_reference 
  {
    pconvex_ref cvr1, cvr2;
 
    scalar_type is_in(const base_node &pt) const
    {
      static base_node *pt1, *pt2;
      static bool isinit = false;
      if (!isinit) { pt1=new base_node(1); pt2=new base_node(1); isinit=true; }
      dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
      assert(pt.size() == cvs->dim());
      if (pt1->size() != n1) { delete pt1; pt1 = new base_node(n1); }
      if (pt2->size() != n2) { delete pt2; pt2 = new base_node(n2); }
      std::copy(pt.begin(), pt.begin()+n1, pt1->begin());
      std::copy(pt.begin()+n1,   pt.end(), pt2->begin());
      return std::max(cvr1->is_in(*pt1), cvr2->is_in(*pt2));
    }
    scalar_type is_in_face(short_type f, const base_node &pt) const
    { // ne controle pas si le point est dans le convexe mais si un point
      // supposé appartenir au convexe est dans une face donnée
      static base_node *pt1, *pt2;
      static bool isinit = false;
      if (!isinit) { pt1=new base_node(1); pt2=new base_node(1); isinit=true; }
      dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
      assert(pt.size() == cvs->dim());
      if (pt1->size() != n1) { delete pt1; pt1 = new base_node(n1); }
      if (pt2->size() != n2) { delete pt2; pt2 = new base_node(n2); }
      std::copy(pt.begin(), pt.begin()+n1, pt1->begin());
      std::copy(pt.begin()+n1,   pt.end(), pt2->begin());

      if (f < cvr1->structure()->nb_faces()) return cvr1->is_in_face(f, *pt1);
        else return cvr2->is_in_face(f - cvr1->structure()->nb_faces(), *pt2);
    }


    _product_ref(const _product_ref_light &ls)
    { 
      cvr1 = ls.cvr1; cvr2 = ls.cvr2;
      *((convex<base_node> *)(this)) = convex_product(*(ls.cvr1), *(ls.cvr2));
      _normals.resize(cvs->nb_faces());
      base_vector null(cvs->dim()); null.fill(0.0);
      std::fill(_normals.begin(), _normals.end(), null);
      for (size_type r = 0; r < cvr1->structure()->nb_faces(); r++)
	std::copy(cvr1->normals()[r].begin(), cvr1->normals()[r].end(),
		  _normals[r].begin());
      for (size_type r = 0; r < cvr2->structure()->nb_faces(); r++)
	std::copy(cvr2->normals()[r].begin(), cvr2->normals()[r].end(),
		  _normals[r+cvr1->structure()->nb_faces()].begin()
		  + cvr1->structure()->dim());
    }
  };

  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b)
  { 
    static dal::FONC_TABLE<_product_ref_light, _product_ref> *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::FONC_TABLE<_product_ref_light, _product_ref>();
      isinit = true;
    }
    return tab->add(_product_ref_light(a, b));
  }

  pconvex_ref parallelepiped_of_reference(dim_type nc)
  {
    static dal::dynamic_array<pconvex_ref> *ptab;
    static dim_type ncd = 1;
    static bool isinit = false;
    if (!isinit) {
      ptab = new dal::dynamic_array<pconvex_ref>();
      isinit = true;
    }

    if (nc <= 1) return simplex_of_reference(nc);
    if (nc > ncd)
    { 
      (*ptab)[nc] = convex_ref_product(parallelepiped_of_reference(nc-1),
				       simplex_of_reference(1));
      ncd = nc;
    }
    return (*ptab)[nc];
  }

  /* ******************************************************************** */
  /* multiple structures.                                                 */
  /* ******************************************************************** */

//   struct _cv_mul_ref_light
//   {
//     pconvex_ref cv;
//     dim_type mul;
//     bool operator < (const _cv_mul_ref_light &ls) const
//     {
//       if (cv < ls.cv) return true; if (cv > ls.cv) return false; 
//       if (mul < ls.mul) return true; return false;
//     }
//     _cv_mul_ref_light(pconvex_ref a, dim_type m) { cv = a; mul = m; }
//     _cv_mul_ref_light(void) { }
   
//   };


//   struct _cv_mul_ref : public convex_of_reference
//   {
//     pconvex_ref cv;
//     scalar_type is_in(const base_node &pt) const
//     { return cv->is_in(pt); }
//     scalar_type is_in_face(short_type f, const base_node &pt) const
//     { return cv->is_in_face(f, pt); }

//     _cv_mul_ref(const _cv_mul_ref_light &ls)
//     { 
//       cv = ls.cv;
//       *((convex<base_node> *)(this)) = convex_multiply(*(ls.cv), ls.mul);
//       _normals.resize(cvs->nb_faces());
//       std::copy(cv->normals().begin(), cv->normals().end(), _normals.begin());
//     }
//   };

//   pconvex_ref multiply_convex_of_reference(pconvex_ref a, dim_type n)
//   { 
//     static dal::FONC_TABLE<_cv_mul_ref_light, _cv_mul_ref> *tab;
//      static bool isinit = false;
//    if (!isinit) {
//      tab = new dal::FONC_TABLE<_cv_mul_ref_light, _cv_mul_ref>();
//      isinit = true;
//    }
//     assert(n != 0); if (n == 1) return a;
//     return tab->add(_cv_mul_ref_light(a, n));
//   }

  /* ******************************************************************** */
  /* nonconforming triangle structure.                                    */
  /* ******************************************************************** */

//   struct _ncf_triangle_ref : public _K_simplex_of_ref
//   {
//     _ncf_triangle_ref(void) : _K_simplex_of_ref(_K_simplex_ref_light(2, 1))
//     {
//       cvs = nonconforming_triangle_structure();
//       for (int i = 0; i < 3; i++) points()[i].fill(0.5);
//       points()[1][0] = 0.0; points()[2][1] = 0.0;
//     }
//   };

//   pconvex_ref nonconforming_triangle_ref(void)
//   {
//     static _ncf_triangle_ref *cs;
//     static bool initialized = false;
//     if (!initialized) { cs = new _ncf_triangle_ref(); initialized = true; }
//     return cs;
//   }

}  /* end of namespace bgeot.                                              */
