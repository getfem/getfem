// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_convex_ref.cc : convexes of reference
//           
// Date    : Septembre 28, 2001.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2001-2005 Yves Renard
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

#include <dal_singleton.h>
#include <bgeot_convex_ref.h>
#include <bgeot_simplexify.h>

namespace bgeot
{
  
  /* ********************************************************************* */
  /*       Point tab storage.                                              */
  /* ********************************************************************* */

  int compare_stored_point_tab::operator()(const stored_point_tab &x,
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

  stored_point_tab_tab &stored_point_tab_tab::instance() {
    return dal::singleton<stored_point_tab_tab>::instance();
  }

  struct org_stored_point_tab_data : public std::vector<pstored_point_tab> { };

  pstored_point_tab org_stored_point_tab(size_type n) {
    std::vector<pstored_point_tab> &tab = 
      dal::singleton<org_stored_point_tab_data>::instance();
    while (n >= tab.size()) {
      stored_point_tab spt(1);
      spt[0].resize(tab.size()); spt[0].fill(0.0); 
      tab.push_back(store_point_tab(spt));
    }
    return tab[n];
  }

  struct cleanup_simplexified_convexes : public dal::ptr_collection<mesh_structure> {};

  /* should be called on the basic_convex_ref */
  const mesh_structure*
  convex_of_reference::simplexified_convex() const {    
    if (psimplexified_convex == NULL) {
      psimplexified_convex = new mesh_structure();
      dal::singleton<cleanup_simplexified_convexes>::instance().push_back(psimplexified_convex);
      if (this != basic_convex_ref()) 
	DAL_THROW(to_be_done_error, 
                  "always use simplexified_convex on the basic_convex_ref() [this=" << 
                  nb_points() << ", basic=" << basic_convex_ref()->nb_points());
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

  struct K_simplex_ref_light_
  {
    dim_type N; short_type K;
    bool operator < (const K_simplex_ref_light_ &l) const
    {
      if (N < l.N) return true; if (N > l.N) return false; 
      if (K < l.K) return true; return false;
    }
    K_simplex_ref_light_(dim_type n, short_type k) { N = n; K = k; }
    K_simplex_ref_light_(void) { }
  };

  class K_simplex_of_ref_ : public convex_of_reference 
  {
    public :
      scalar_type is_in(const base_node &pt) const
      { // return a negative or null number if pt is in the convex
	if (pt.size() != cvs->dim())
	  throw dimension_error(
		    "K_simplex_of_ref_::is_in : Dimension does not match");
	scalar_type e = -1.0, r = 0.0;
	base_node::const_iterator it = pt.begin(), ite = pt.end();
	for (; it != ite; e += *it, ++it) r = std::min(r, *it);
	return std::max(-r, e);
      }
      scalar_type is_in_face(short_type f, const base_node &pt) const
      { // return a null number if pt is in the face of the convex
	if (pt.size() != cvs->dim())
	  throw dimension_error(
		  "K_simplex_of_ref_::is_in_face : Dimension does not match");
	if (f > 0) return gmm::abs(pt[f-1]);
	scalar_type e = -1.0;
	base_node::const_iterator it = pt.begin(), ite = pt.end();
	for (; it != ite; e += *it, ++it);
	return gmm::abs(e);
      }
      K_simplex_of_ref_(const K_simplex_ref_light_ &ls)
      {
	cvs = simplex_structure(ls.N, ls.K);
	size_type R = cvs->nb_points();
	points().resize(R);
	normals_.resize(ls.N+1);
	base_node null(ls.N); null.fill(0.0);
	std::fill(normals_.begin(), normals_.end(), null);
	std::fill(points().begin(), points().end(), null);
	for (size_type i = 1; i <= ls.N; ++i)
	  normals_[i][i-1] = -1.0;
	if (ls.N > 0)
	  std::fill(normals_[0].begin(), normals_[0].end(),
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

  struct simplex_of_reference_FONC_TABLE : 
    public dal::FONC_TABLE<K_simplex_ref_light_, K_simplex_of_ref_> {};

  pconvex_ref simplex_of_reference(dim_type nc, short_type k) {
    simplex_of_reference_FONC_TABLE &tab = 
      dal::singleton<simplex_of_reference_FONC_TABLE>::instance();
    bgeot::convex_of_reference * p1 = tab.add(K_simplex_ref_light_(nc, 1));
    bgeot::convex_of_reference * pk = tab.add(K_simplex_ref_light_(nc, k));
    p1->attach_basic_convex_ref(p1);
    pk->attach_basic_convex_ref(p1);
    return pk;
  }

  /* products.                                                             */

  struct product_ref_light_ {
    pconvex_ref cvr1, cvr2;
    bool operator < (const product_ref_light_ &ls) const
    {
      if (cvr1 < ls.cvr1) return true; if (cvr1 > ls.cvr1) return false; 
      if (cvr2 < ls.cvr2) return true; return false;
    }
    product_ref_light_(pconvex_ref a, pconvex_ref b)
    { cvr1 = a; cvr2 = b; }
    product_ref_light_(void) { }
  };

  struct product_ref_ : public convex_of_reference {
    pconvex_ref cvr1, cvr2;
 
    scalar_type is_in(const base_node &pt) const;
    scalar_type is_in_face(short_type f, const base_node &pt) const;


    product_ref_(const product_ref_light_ &ls);
  };

  scalar_type product_ref_::is_in(const base_node &pt) const {
    dim_type n1 = cvr1->structure()->dim(), n2 = cvr2->structure()->dim();
    base_node pt1(n1), pt2(n2);
    if (pt.size() != cvs->dim())
      throw dimension_error(
			    "product_ref_::is_in : Dimension does not match");
    std::copy(pt.begin(), pt.begin()+n1, pt1.begin());
    std::copy(pt.begin()+n1,   pt.end(), pt2.begin());
    return std::max(cvr1->is_in(pt1), cvr2->is_in(pt2));
  }

  scalar_type product_ref_::is_in_face(short_type f, const base_node &pt) const
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
  
  
  product_ref_::product_ref_(const product_ref_light_ &ls) { 
    if (ls.cvr1->structure()->dim() < ls.cvr2->structure()->dim())
      DAL_WARNING(1, "Illegal convex : swap your operands: dim(cv1)=" << 
		  int(ls.cvr1->structure()->dim()) << " < dim(cv2)=" << 
		  int(ls.cvr2->structure()->dim()));
    cvr1 = ls.cvr1; cvr2 = ls.cvr2;
    *((convex<base_node> *)(this)) 
      = convex_direct_product(*(ls.cvr1), *(ls.cvr2));
    normals_.resize(cvs->nb_faces());
    base_small_vector null(cvs->dim()); null.fill(0.0);
    std::fill(normals_.begin(), normals_.end(), null);
    for (size_type r = 0; r < cvr1->structure()->nb_faces(); r++)
      std::copy(cvr1->normals()[r].begin(), cvr1->normals()[r].end(),
		normals_[r].begin());
    for (size_type r = 0; r < cvr2->structure()->nb_faces(); r++)
      std::copy(cvr2->normals()[r].begin(), cvr2->normals()[r].end(),
		normals_[r+cvr1->structure()->nb_faces()].begin()
		+ cvr1->structure()->dim());
  }

  struct convex_direct_product_FONC_TABLE : 
    public dal::FONC_TABLE<product_ref_light_, product_ref_> {};

  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b) { 
    convex_direct_product_FONC_TABLE &tab = 
      dal::singleton<convex_direct_product_FONC_TABLE>::instance();
    bgeot::convex_of_reference *bprod = 
      tab.add(product_ref_light_(a->basic_convex_ref(), b->basic_convex_ref()));
    bgeot::convex_of_reference *prod = tab.add(product_ref_light_(a, b));
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

  /* equilateral ref convexes are used for estimation of convex quality */
  class equilateral_simplex_of_ref_ : public convex_of_reference {
  public:
    scalar_type is_in(const base_node &pt) const {
      scalar_type d;
      for (size_type f = 0; f < normals().size(); ++f) {
        const base_node &x0 = (f ? points()[f-1] : points().back());
        scalar_type v = gmm::vect_sp(pt-x0, normals()[f]);
        if (f == 0) d = v; else d = std::max(d,v);
      }
      return d;
    }
    scalar_type is_in_face(short_type f, const base_node &pt) const {
      const base_node &x0 = (f ? points()[f-1] : points().back());
      return gmm::abs(gmm::vect_sp(pt-x0, normals()[f])); 
    }
    equilateral_simplex_of_ref_(size_type N) {
      pconvex_ref prev = equilateral_simplex_of_reference(N-1);
      basic_convex_ref_ = this;
      cvs = simplex_structure(N, 1);
      points().resize(N+1);
      normals_.resize(N+1);
      base_node G(N); G.fill(0.);
      for (size_type i=0; i < N+1; ++i) {
        points()[i].resize(N);
        if (i != N) {
          std::copy(prev->points()[i].begin(), prev->points()[i].end(),
		    points()[i].begin());
          points()[i][N-1] = 0.;
        } else {
          points()[i] = 1./N * G;
          points()[i][N-1] = sqrt(1. - gmm::vect_norm2_sqr(points()[i]));
        }
        G += points()[i];
      }
      gmm::scale(G, 1./(N+1));
      for (size_type f=0; f < N+1; ++f) {
        normals_[f] = G - points()[f]; 
        gmm::scale(normals_[f], 1/gmm::vect_norm2(normals_[f]));
      }
    }
  };


  struct equilateral_simplex_list :
    public dal::ptr_collection<equilateral_simplex_of_ref_> { };

  pconvex_ref equilateral_simplex_of_reference(dim_type nc) {
    equilateral_simplex_list &stab = 
      dal::singleton<equilateral_simplex_list>::instance();
    if (nc <= 1) return simplex_of_reference(nc);
    if (nc >= stab.size()) stab.resize(nc+1);
    if (stab[nc] == 0) stab[nc] = new equilateral_simplex_of_ref_(nc);
    return stab[nc];
  }

  /* generic convex with n global nodes      */

  class generic_dummy_ : public convex_of_reference {
  public:
    scalar_type is_in(const base_node &) const
    { DAL_THROW(failure_error, "Information not available here"); }
    scalar_type is_in_face(short_type, const base_node &) const 
    { DAL_THROW(failure_error, "Information not available here"); }
  
    generic_dummy_(dim_type d, size_type n) {
      cvs = generic_dummy_structure(d, n);
      points().resize(n);
      normals_.resize(0);
      base_node P(d);
      std::fill(P.begin(), P.end(), scalar_type(1)/scalar_type(20));
      std::fill(points().begin(), points().end(), P);
    }
  };

  pconvex_ref generic_dummy_convex_ref(dim_type nc, size_type n) {
    static std::vector< std::vector<pconvex_ref> > tab;
    if (size_type(nc)+1 > tab.size()) tab.resize(nc+1);
    if (n+1 > tab[nc].size()) {
      size_type d = tab[nc].size();
      tab[nc].resize(n+1);
      for (size_type i = d; i <= n; ++i)
	tab[nc][i] = new generic_dummy_(nc, i);
    }
    return tab[nc][n];
  }

}  /* end of namespace bgeot.                                              */
