/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool (bgeot)                                 */
/* File    :  bgeot_geotrans_inv.h : Allows to inverse geometric transf-   */
/*     	      ormations and to localize a set of points.    	       	   */
/*                                                                         */
/* Date : December 20, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2003  Yves Renard.                                   */
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

#ifndef BGEOT_GEOTRANS_INV_H__
#define BGEOT_GEOTRANS_INV_H__

#include <bgeot_geometric_trans.h>
#include <bgeot_vector.h>
#include <bgeot_kdtree.h>

namespace bgeot {
  /** 
      does the inversion of the geometric transformation for a given convex
  */
  class geotrans_inv_convex {
    size_type N, P;
    base_poly PO;
    base_matrix a, pc, grad, B0, CS;
    pgeometric_trans pgt;
    std::vector<base_node> cvpts; /* used only for non-linear geotrans -- we should use the matrix a instead... */
    scalar_type EPS;
  public:
    geotrans_inv_convex(scalar_type e=10e-12) : N(0), P(0), pgt(0), EPS(e) {};
    template<class TAB> geotrans_inv_convex(const convex<base_node, TAB> &cv,
					    pgeometric_trans pgt_, 
                                            scalar_type e=10e-12)
      : N(0), P(0), pgt(0), EPS(e) { init(cv,pgt_); }
    template<class TAB> void init(const convex<base_node, TAB> &cv,
				  pgeometric_trans pgt_);
    
    /**
       given the node on the real element, returns the node
       on the reference element (even if it is outside of the reference convex)
       @return true if the n is inside the convex
       @param n node on the real element 
       @param n_ref computed node on the reference convex
    */
    bool invert(const base_node& n, base_node& n_ref, scalar_type IN_EPS=1e-12) {
      n_ref.resize(pgt->structure()->dim());
      if (pgt->is_linear()) {
        return invert_lin(n, n_ref,IN_EPS);
      } else return invert_nonlin(n, n_ref,IN_EPS);
    }
  private:
    bool invert_lin(const base_node& n, base_node& n_ref, scalar_type IN_EPS);
    bool invert_nonlin(const base_node& n, base_node& n_ref, scalar_type IN_EPS);
  };

  template<class TAB> void geotrans_inv_convex::init(const convex<base_node,
						     TAB> &cv, pgeometric_trans pgt_) {
    bool geotrans_changed = (pgt != pgt_); if (geotrans_changed) pgt = pgt_;
    if (!cv.points().size()) DAL_INTERNAL_ERROR("empty points!");
    if (P != cv.points()[0].size()) { P = cv.points()[0].size(); geotrans_changed = true; }
    if (geotrans_changed) {
      N = pgt->structure()->dim();
      pc.resize(pgt->nb_points() , N);
      grad.resize(P,N); B0.resize(N,P); CS.resize(N,N);
      a.resize(P, pgt->nb_points());
    }
    for (size_type j = 0; j < pgt->nb_points(); ++j) {// à optimiser !!
      base_node pt = cv.points()[j]; /* need a temporary storage since cv.points()[j] may not
					be a reference to a base_node, but a temporary base_node !!
				     */
      base_node::const_iterator it = pt.begin();
      for (size_type i = 0; i < P; ++i) { 
        this->a(i,j) = it[i];
      }
    }
    if (pgt->is_linear()) {
      if (geotrans_changed) {
	for (size_type i = 0; i < pgt->nb_points(); ++i) {
	  for (dim_type n = 0; n < N; ++n) { 
	    PO = pgt->poly_vector()[i]; PO.derivative(n); this->pc(i,n) = PO[0]; // optimisable
	  }
	}
      }
      // computation of the pseudo inverse
      gmm::mult(a, pc, grad);
      if (N != P) {
	gmm::mult(gmm::transposed(grad), grad, CS);
        gmm::lu_inverse(CS);
	gmm::mult(gmm::transposed(CS), gmm::transposed(grad), B0);
      }
      else {
        // L'inversion peut être optimisée par le non calcul global de B0
        // et la resolution d'un système linéaire.
        gmm::lu_inverse(grad); B0 = grad;
      }
    } else { /* not much to precompute for non-linear geometric transformations .. */
      cvpts.resize(cv.nb_points());
      for (size_type j = 0; j < pgt->nb_points(); ++j) 
        cvpts[j] = cv.points()[j];
    }
  }


  /**
     handles the geometric inversion for a given (supposedly quite large) set of points
  */
  class geotrans_inv
  {
  protected :
    mutable kdtree tree;
    scalar_type EPS;
    geotrans_inv_convex gic;
  public :
    void clear(void) { tree.clear(); }
    /// Add the points contained in c to the list of points.
    template<class CONT> void add_points(const CONT &c) {
      tree.reserve(std::distance(c.begin(),c.end()));
      typename CONT::const_iterator it = c.begin(), ite = c.end();
      for (; it != ite; ++it) tree.add_point(*it);
    }

    /// Number of points.
    size_type nb_points(void) const { return tree.nb_points(); }
    /// Add point p to the list of points.
    size_type add_point(base_node p) { return tree.add_point(p); }
    void add_point_with_id(base_node p,size_type id) { tree.add_point_with_id(p,id); }
      
    /// Find all the points present in the box between min and max.
    size_type points_in_box(kdtree_tab_type &ipts,
			    const base_node &min, 
			    const base_node &max) const {
      tree.points_in_box(ipts, min, max);
      return ipts.size();
    }

    /** Search all the points in the convex cv, which is the transformation
     *  of the convex cref via the geometric transformation pgt.
     *
     *  IMPORTANT : It is assumed that the whole convex is include in the
     *  minmax box of its nodes times a factor 1.2. If the transformation is
     *  linear, the factor is 1.0.
     *
     *  @param cv the convex points (as given by getfem_mesh::convex(ic))
     *
     *  @param pgt the geometric trans (as given by
     *  getfem_mesh::trans_of_convex(ic))
     *
     *  @param pftab container for the coordinates of points in the reference
     *  convex (should be of size nb_points())
     *
     *  @param itab the indices of points found in the convex
     *
     *  @return the number of points in the convex (i.e. the size of itab,
     *  and pftab)
     */
    template<class TAB, class CONT1, class CONT2>
    size_type points_in_convex(const convex<base_node, TAB> &cv,
			       pgeometric_trans pgt,
			       CONT1 &pftab, CONT2 &itab, bool bruteforce=false);
      
    geotrans_inv(scalar_type EPS_ = 10E-12) : EPS(EPS_) {}
  };



  template<class TAB, class CONT1, class CONT2>
  size_type geotrans_inv::points_in_convex(const convex<base_node, TAB> &cv,
					   pgeometric_trans pgt,
					   CONT1 &pftab, CONT2 &itab, bool bruteforce) {
    base_node min, max; /* bound of the box enclosing the convex */
    size_type nbpt = 0; /* nb of points in the convex */
    kdtree_tab_type boxpts;
    bounding_box(min, max, cv.points(), pgt);
    gic.init(cv,pgt);
    /* get the points in a box enclosing the convex */
    if (!bruteforce) points_in_box(boxpts, min, max);
    else boxpts = tree.points();
    /* and invert the geotrans, and check if the obtained point is 
       inside the reference convex */
    for (size_type l = 0; l < boxpts.size(); ++l) {
      base_node pt_ref;
      if (gic.invert(boxpts[l].n, pftab[nbpt], EPS)) {
	itab[nbpt++] = boxpts[l].i;
      }
    }
    return nbpt;
  }



    

}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_GEOMETRIC_TRANSFORMATION_H__                              */
