// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool (bgeot)
// File    : bgeot_kdtree.h : implements a k-d-tree for window search on
//           a set of points.
// Date    : January 2004.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Julien Pommier
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

#ifndef BGEOT_KDTREE_H
#define BGEOT_KDTREE_H

/** \file bgeot_kdtree.h
    \brief Simple implementation of a KD-tree.

    Basically, a KD-tree is a balanced N-dimensional tree.
*/
#include <bgeot_vector.h>

namespace bgeot {
  struct kdtree_elt_base;
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
   *
   * Once the tree have been built, it is possible to query very
   * quickly for the list of points lying in a given box. Note that
   * this is not a dynamic structure: once you start to call
   * kdtree::points_in_box, you should not use anymore kdtree::add_point.
   */
  class kdtree : public boost::noncopyable {
    dim_type N; /* dimension of points */
    kdtree_elt_base *tree;
    kdtree_tab_type pts;
  public:
    kdtree() : N(0), tree(0) {}
    ~kdtree() { clear_tree(); }
    /// reset the tree, remove all points
    void clear() { clear_tree(); pts.clear(); N = 0; }
    void reserve(size_type n) { pts.reserve(n); }
    /// insert a new point
    size_type add_point(const base_node& n) { 
      size_type i = pts.size(); add_point_with_id(n,i); return i;
    }
    /// insert a new point, with an associated number.
    void add_point_with_id(const base_node& n, size_type i) {
      if (pts.size() == 0) N = n.size(); 
      else if (N != n.size()) 
	DAL_THROW(dal::failure_error, "invalid dimension");
      if (tree) clear_tree();
      pts.push_back(index_node_pair(i, n));
    }
    size_type nb_points() const { return pts.size(); }
    const kdtree_tab_type &points() const { return pts; }
    /* fills ipts with the indexes of points in the box 
       [min,max] */
    void points_in_box(kdtree_tab_type &inpts,
		       const base_node &min, 
		       const base_node &max);
  private:
    typedef std::vector<size_type>::const_iterator ITER;  
    void clear_tree();
  };
}
#endif
