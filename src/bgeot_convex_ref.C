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
/* Copyright (C) 2001-2002  Yves Renard.                                   */
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


#include <bgeot_convex_ref.h>
#include <bgeot_simplexify.h>

namespace bgeot
{
  
  /* ********************************************************************* */
  /*       Point tab storage.                                              */
  /* ********************************************************************* */

  int comp_stored_point_tab::operator()(const stored_point_tab &x,
					const stored_point_tab &y) const {
    std::vector<base_node>::const_iterator it1 = x.begin(), it2 = y.begin();
    base_node::const_iterator itn1, itn2, itne;
    for ( ; it1 != x.end() && it2 != y.end() ; ++it1, ++it2) {
      if ((*it1).size() < (*it2).size()) return -1;
      if ((*it1).size() > (*it2).size()) return 1;
      itn1 = (*it1).begin(); itne = (*it1).end(); itn2 = (*it2).begin();
      for ( ; itn1 != itne ; ++itn1, ++itn2)
	if (*itn1 < *itn2) return -1;
	else if (*itn1 > *itn2) return 1;
    }
    if (it2 != y.end()) return -1;
    if (it1 != x.end()) return 1;
    return 0;
  }

  dal::dynamic_tree_sorted<stored_point_tab, comp_stored_point_tab>
    *_stored_point_tab_tab;
  bool isinit_stored_point_tab_tab = false;

  pstored_point_tab org_stored_point_tab(size_type n)
  {
    static std::vector<pstored_point_tab> *tab;
    static bool is_init = false;
    
    if (!is_init) { tab = new std::vector<pstored_point_tab>(); is_init = true; }
    while (n >= tab->size())
    {
      size_type i = tab->size();
      stored_point_tab spt;
      tab->resize(i+1); spt.resize(1); spt[0].resize(i); spt[0].fill(0.0); 
      (*tab)[i] = store_point_tab(spt);
    }
    return (*tab)[n];
  }

  /* should be called on the basic_convex_ref */
  const mesh_structure*
  convex_of_reference::simplexified_convex() const {    
    if (psimplexified_convex == NULL) {
      psimplexified_convex = new mesh_structure();
      if (this != basic_convex_ref()) 
	DAL_THROW(to_be_done_error, "always use simplexified_convex on the basic_convex_ref() [this=" << nb_points() << ", basic=" << basic_convex_ref()->nb_points());
      mesh_structure ms;
      std::vector<size_type> ipts(nb_points());
      for (size_type i=0; i < ipts.size(); ++i) ipts[i] = i;
      ms.add_convex(structure(), ipts.begin());
      ms.to_edges();
      bgeot::simplexify(ms,*psimplexified_convex, points(), 
			std::max(structure()->dim(),dim_type(1)), 1e-12);
    }
    return psimplexified_convex;
  }

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
	if (pt.size() != cvs->dim())
	  throw dimension_error(
		    "_K_simplex_of_ref::is_in : Dimension does not match");
	scalar_type e = -1.0, r = 0.0;
	base_node::const_iterator it = pt.begin(), ite = pt.end();
	for (; it != ite; e += *it, ++it) r = std::min(r, *it);
	return std::max(-r, e);
      }
      scalar_type is_in_face(short_type f, const base_node &pt) const
      { // return a null number if pt is in the face of the convex
	if (pt.size() != cvs->dim())
	  throw dimension_error(
		  "_K_simplex_of_ref::is_in_face : Dimension does not match");
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
	if (ls.N > 0)
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
	  for (size_type r = 0; r < R; ++r) {
	    points()[r] = c;
	    if (ls.K != 0 && ls.N > 0) {
	      l = 0; c[l] += 1.0 / scalar_type(ls.K); sum++;
	      while (sum > ls.K) {
		sum -= int(floor(0.5+(c[l] * ls.K)));
		c[l] = 0.0; l++; if (l == ls.N) break;
		c[l] += 1.0 / scalar_type(ls.K); sum++;
	      }
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
    bgeot::convex_of_reference * p1 = tab->add(_K_simplex_ref_light(nc, 1));
    bgeot::convex_of_reference * pk = tab->add(_K_simplex_ref_light(nc, k));
    p1->attach_basic_convex_ref(p1);
    pk->attach_basic_convex_ref(p1);
    return pk;
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
 
    scalar_type is_in(const base_node &pt) const;
    scalar_type is_in_face(short_type f, const base_node &pt) const;


    _product_ref(const _product_ref_light &ls);
  };

  scalar_type _product_ref::is_in(const base_node &pt) const {
    dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
    base_node pt1(n1), pt2(n2);
    if (pt.size() != cvs->dim())
      throw dimension_error(
			    "_product_ref::is_in : Dimension does not match");
    std::copy(pt.begin(), pt.begin()+n1, pt1.begin());
    std::copy(pt.begin()+n1,   pt.end(), pt2.begin());
    return std::max(cvr1->is_in(pt1), cvr2->is_in(pt2));
  }

  scalar_type _product_ref::is_in_face(short_type f, const base_node &pt) const
  { // ne controle pas si le point est dans le convexe mais si un point
    // supposé appartenir au convexe est dans une face donnée
    dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
    base_node pt1(n1), pt2(n2);
    if (pt.size() != cvs->dim())
      DAL_THROW(dimension_error, "Dimensions mismatch");
    std::copy(pt.begin(), pt.begin()+n1, pt1.begin());
    std::copy(pt.begin()+n1,   pt.end(), pt2.begin());
    if (f < cvr1->structure()->nb_faces()) return cvr1->is_in_face(f, pt1);
    else return cvr2->is_in_face(f - cvr1->structure()->nb_faces(), pt2);
  }
  
  
  _product_ref::_product_ref(const _product_ref_light &ls) { 
    if (ls.cvr1->structure()->dim() < ls.cvr2->structure()->dim())
      DAL_WARNING(1, "Illegal convex : swap your operands: dim(cv1)=" << int(ls.cvr1->structure()->dim()) << " < dim(cv2)=" << int(ls.cvr2->structure()->dim()));
    cvr1 = ls.cvr1; cvr2 = ls.cvr2;
    *((convex<base_node> *)(this)) 
      = convex_direct_product(*(ls.cvr1), *(ls.cvr2));
    _normals.resize(cvs->nb_faces());
    base_small_vector null(cvs->dim()); null.fill(0.0);
    std::fill(_normals.begin(), _normals.end(), null);
    for (size_type r = 0; r < cvr1->structure()->nb_faces(); r++)
      std::copy(cvr1->normals()[r].begin(), cvr1->normals()[r].end(),
		_normals[r].begin());
    for (size_type r = 0; r < cvr2->structure()->nb_faces(); r++)
      std::copy(cvr2->normals()[r].begin(), cvr2->normals()[r].end(),
		_normals[r+cvr1->structure()->nb_faces()].begin()
		+ cvr1->structure()->dim());
  }
  

  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b) { 
    static dal::FONC_TABLE<_product_ref_light, _product_ref> *tab = 0;
    if (!tab) tab = new dal::FONC_TABLE<_product_ref_light, _product_ref>();
    bgeot::convex_of_reference *bprod = tab->add(_product_ref_light(a->basic_convex_ref(), b->basic_convex_ref()));
    bgeot::convex_of_reference *prod = tab->add(_product_ref_light(a, b));
    bprod->attach_basic_convex_ref(bprod);
    prod->attach_basic_convex_ref(bprod);
    return prod;
  }

  pconvex_ref parallelepiped_of_reference(dim_type nc) {
    static dal::dynamic_array<pconvex_ref> *ptab = 0;
    static dim_type ncd = 1;
    if (!ptab) ptab = new dal::dynamic_array<pconvex_ref>();
     

    if (nc <= 1) return simplex_of_reference(nc);
    if (nc > ncd)
    { 
      (*ptab)[nc] = convex_ref_product(parallelepiped_of_reference(nc-1),
				       simplex_of_reference(1));
      ncd = nc;
    }
    return (*ptab)[nc];
  }

}  /* end of namespace bgeot.                                              */
