/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2000-2015 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/

/**@file bgeot_geotrans_inv.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date December 20, 2000.
   @brief Inversion of geometric transformations.

   Inversion means: given a set of convexes and a point, find:

     - a subset of candidate convexes, which are likely to contain the
     point (using bgeot::kdtree).

     - on these candidate convexes, invert the geometric
     transformation, i.e. find the corresponding coordinates on the
     reference element.

   Inversion of a geometric transformation is not a trivial task,
   especially with non-linear geometric transformations. This is the
   central part of interpolation routines from a getfem::mesh_fem onto
   another.
*/

#ifndef BGEOT_GEOTRANS_INV_H__
#define BGEOT_GEOTRANS_INV_H__

#include "bgeot_geometric_trans.h"
#include "bgeot_small_vector.h"
#include "bgeot_kdtree.h"

namespace bgeot {
  /** 
      does the inversion of the geometric transformation for a given convex
  */
  class geotrans_inv_convex {
    size_type N, P;
    base_matrix G, pc, K, B, CS;
    pgeometric_trans pgt;
    std::vector<base_node> cvpts; /* used only for non-linear geotrans
				  -- we should use the matrix a instead... */
    scalar_type EPS;
  public:
    const base_matrix &get_G() const { return G; }
    geotrans_inv_convex(scalar_type e=10e-12) : N(0), P(0), pgt(0), EPS(e) {};
    template<class TAB> geotrans_inv_convex(const convex<base_node, TAB> &cv,
					    pgeometric_trans pgt_, 
                                            scalar_type e=10e-12)
      : N(0), P(0), pgt(0), EPS(e) { init(cv.points(),pgt_); }
    geotrans_inv_convex(const std::vector<base_node> &nodes,
			pgeometric_trans pgt_, scalar_type e=10e-12)
      : N(0), P(0), pgt(0), EPS(e) { init(nodes,pgt_); }
    template<class TAB> void init(const TAB &nodes, pgeometric_trans pgt_);
    
    /**
       given the node on the real element, returns the node on the
       reference element (even if it is outside of the ref. convex).
       
       If the geometric transformation is not invertible at point n,
       an exception is thrown.

       @return true if the n is inside the convex

       @param n node on the real element 

       @param n_ref computed node on the reference convex

       @param IN_EPS a threshold.
    */
    bool invert(const base_node& n, base_node& n_ref,
		scalar_type IN_EPS=1e-12);

    /**
       given the node on the real element, returns the node
       on the reference element (even if it is outside of the ref. convex).
       
       This version will not throw an exception if the geometric
       transformation is not invertible at point n. 
       
       @return true if the n is inside the convex

       @param n node on the real element 

       @param n_ref computed node on the reference convex

       @param converged on output, will be set to true if the
       geometric transformation could be inverted.

       @param IN_EPS a threshold.
    */
    bool invert(const base_node& n, base_node& n_ref, bool &converged, 
		scalar_type IN_EPS=1e-12);
  private:
    bool invert_lin(const base_node& n, base_node& n_ref, scalar_type IN_EPS);
    bool invert_nonlin(const base_node& n, base_node& n_ref,
		       scalar_type IN_EPS, bool &converged, bool throw_except);
    void update_B();
    friend class geotrans_inv_convex_bfgs;
  };

  template<class TAB>
  void geotrans_inv_convex::init(const TAB &nodes,  pgeometric_trans pgt_) {
    bool geotrans_changed = (pgt != pgt_); if (geotrans_changed) pgt = pgt_;
    GMM_ASSERT3(!nodes.empty(), "empty points!");
    if (N != nodes[0].size())
      { N = nodes[0].size(); geotrans_changed = true; }
    if (geotrans_changed) {
      P = pgt->structure()->dim();
      pc.resize(pgt->nb_points() , P);
      K.resize(N,P); B.resize(N,P); CS.resize(P,P);
      G.resize(N, pgt->nb_points());
    }
    vectors_to_base_matrix(G, nodes);
    if (pgt->is_linear()) {
      if (geotrans_changed) {
	base_node Dummy(P);
	pgt->poly_vector_grad(Dummy, pc);
      }
      // computation of the pseudo inverse
      update_B();
    } else { /* not much to precompute for non-linear geometric
		transformations .. */
      cvpts.assign(nodes.begin(), nodes.end());
    }
  }


  /**
     handles the geometric inversion for a given (supposedly quite large)
     set of points
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
    void add_point_with_id(base_node p,size_type id)
    { tree.add_point_with_id(p,id); }
      
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
     *  @param cv the convex points (as given by getfem_mesh::convex(ic)).
     *
     *  @param pgt the geometric trans (as given by
     *  getfem_mesh::trans_of_convex(ic)).
     *
     *  @param pftab container for the coordinates of points in the reference
     *  convex (should be of size nb_points())
     *
     *  @param itab the indices of points found in the convex.
     *
     *  @param bruteforce use a brute force search (only for debugging purposes).
     *
     *  @return the number of points in the convex (i.e. the size of itab,
     *  and pftab)
     */
    template<class TAB, class CONT1, class CONT2>
    size_type points_in_convex(const convex<base_node, TAB> &cv,
			       pgeometric_trans pgt,
			       CONT1 &pftab, CONT2 &itab,
			       bool bruteforce=false);
      
    geotrans_inv(scalar_type EPS_ = 10E-12) : EPS(EPS_) {}
  };



  template<class TAB, class CONT1, class CONT2>
  size_type geotrans_inv::points_in_convex(const convex<base_node, TAB> &cv,
					   pgeometric_trans pgt,
					   CONT1 &pftab, CONT2 &itab,
					   bool bruteforce) {
    base_node min, max; /* bound of the box enclosing the convex */
    size_type nbpt = 0; /* nb of points in the convex */
    kdtree_tab_type boxpts;
    bounding_box(min, max, cv.points(), pgt);
    for (size_type k=0; k < min.size(); ++k) { min[k] -= EPS; max[k] += EPS; }
    gic.init(cv.points(),pgt);
    /* get the points in a box enclosing the convex */
    if (!bruteforce) points_in_box(boxpts, min, max);
    else boxpts = tree.points();
    /* and invert the geotrans, and check if the obtained point is 
       inside the reference convex */
    for (size_type l = 0; l < boxpts.size(); ++l) {
      // base_node pt_ref;
      if (gic.invert(boxpts[l].n, pftab[nbpt], EPS)) {
	itab[nbpt++] = boxpts[l].i;
      }
    }
    return nbpt;
  }



    

}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_GEOMETRIC_TRANSFORMATION_H__                              */
