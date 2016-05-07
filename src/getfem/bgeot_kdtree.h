/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2015 Julien Pommier

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

#ifndef BGEOT_KDTREE_H
#define BGEOT_KDTREE_H

/** @file bgeot_kdtree.h
    @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
    @date January 2004.
    @brief Simple implementation of a KD-tree.

    Basically, a KD-tree is a balanced N-dimensional tree.
*/
#include "bgeot_small_vector.h"

namespace bgeot {

  /* generic node for the kdtree */
  struct kdtree_elt_base {
    enum { PTS_PER_LEAF=8 };
    unsigned n; /* 0 => is a tree node, != 0 => tree leaf storing n points */
    bool isleaf() const { return (n != 0); }
    kdtree_elt_base(unsigned n_) : n(n_) {}
    virtual ~kdtree_elt_base() {}
  };

  /// store a point and the associated index for the kdtree. 
  /* std::pair<size_type,base_node> is not ok since it does not
     have a suitable overloaded swap function ...
  */
  struct index_node_pair {
    size_type i;
    base_node n;
    index_node_pair() {}
    index_node_pair(const index_node_pair& o) : i(o.i), n(o.n) {}
    index_node_pair(size_type i_, base_node n_) : i(i_), n(n_) {}
    void swap(index_node_pair& other) { std::swap(i,other.i); n.swap(other.n);}
  };

  /// store a set of points with associated indexes.
  typedef std::vector<index_node_pair> kdtree_tab_type;

  /** Balanced tree over a set of points.
   
  Once the tree have been built, it is possible to query very
  quickly for the list of points lying in a given box. Note that
  this is not a dynamic structure: once you start to call
  kdtree::points_in_box, you should not use anymore kdtree::add_point.
  
  Here is an example of use (which tries to find the mapping between
  the dof of the mesh_fem and the node numbers of its mesh):
  
  @code
  void dof_to_nodes(const getfem::mesh_fem &mf) {
    const getfem::getfem_mesh &m = mf.linked_mesh();
    bgeot::kdtree tree;
    for (dal::bv_visitor ip(m.points().index()); !ip.finished(); ++ip) {
      tree.add_point_with_id(m.points()[ip], ip);
    }
    bgeot::kdtree_tab_type t;
    for (size_type d = 0; d < mf.nb_dof(); ++d) {
      getfem::base_node P(mf.point_of_dof(d)), P2(P);
      for (unsigned i=0; i < P.size(); ++i) {
        P[i] -= 1e-12; P2[i]+= 1e-12;
      }
      tree.points_in_box(t, P, P2);
      if (t.size()) {
        assert(t.size() == 1);
	cout << " dof " << d << " maps to mesh node " << t[0].i << "\n";
      }
    } 
  }
  @endcode
  */
  class kdtree : public boost::noncopyable {
    dim_type N; /* dimension of points */
    std::unique_ptr<kdtree_elt_base> tree;
    kdtree_tab_type pts;
  public:
    kdtree() : N(0) {}
    /// reset the tree, remove all points
    void clear() { clear_tree(); pts = kdtree_tab_type(); N = 0; }
    void reserve(size_type n) { pts.reserve(n); }
    /// insert a new point
    size_type add_point(const base_node& n) { 
      size_type i = pts.size(); add_point_with_id(n,i); return i;
    }
    /// insert a new point, with an associated number.
    void add_point_with_id(const base_node& n, size_type i) {
      if (pts.size() == 0) N = n.size(); 
      else if (N != n.size()) 
	GMM_ASSERT2(false, "invalid dimension");
      if (tree) clear_tree();
      pts.push_back(index_node_pair(i, n));
    }
    size_type nb_points() const { return pts.size(); }
    const kdtree_tab_type &points() const { return pts; }
    /* fills ipts with the indexes of points in the box 
       [min,max] */
    void points_in_box(kdtree_tab_type &ipts,
		       const base_node &min, 
		       const base_node &max);
    /* assigns at ipt the index of the nearest neighbor at location
       pos and returns the square of the distance to this point*/
    scalar_type nearest_neighbor(index_node_pair &ipt,
                                 const base_node &pos);
  private:
    typedef std::vector<size_type>::const_iterator ITER;  
    void clear_tree();
  };
}

#endif
