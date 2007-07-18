// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2007-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#ifndef BGEOT_NODE_TAB_H
#define BGEOT_NODE_TAB_H

/** @file bgeot_node_tab.h
    @author  Yves Renard <Yves.Renard@insa-lyon.fr>
    @date January 2004.
    @brief Structure which dynamically collects points identifying points
    that are nearer than a certain very small distance.

*/
#include "bgeot_vector.h"
#include "dal_tree_sorted.h"
#include "set"

namespace bgeot {
  
  
  /** Store a set of points, identifying points
      that are nearer than a certain very small distance.
  */
  class node_tab : public dal::dynamic_tas<base_node> {
    
  protected :
    
    struct component_comp {
      const dal::dynamic_tas<base_node> *vbn;
      const scalar_type *c;
      std::vector<scalar_type> v;
      bool operator()(size_type i1, size_type i2) const;
      component_comp(const dal::dynamic_tas<base_node> &vbn_,
		     const scalar_type &c_, unsigned dim);
    };
    typedef std::set<size_type, component_comp> sorter;
    
    mutable std::vector<sorter> sorters;
    mutable scalar_type c;
    scalar_type eps, prec_factor, max_radius;
    unsigned dim_;

    void add_sorter(void) const;

  public :
    
    /// reset the array, remove all points
    void clear(void);
    
    /** Search a node in the array. return its index if it exists
	or size_type(-1) otherwise.
    */
    size_type search_node(const base_node &pt) const;
    /** Add a point to the array or identify it with a very close existing
	point.
    */
    size_type add_node(const base_node &pt);
    size_type add(const base_node &pt) { return add_node(pt); }
    void sup_node(size_type i);
    void sup(size_type i) { sup_node(i); }
    // void resort(void) { sorters = std::vector<sorter>(); }
    dim_type dim(void) const { return dim_; }
    void translation(const base_small_vector &V);
    void transformation(const base_matrix &M);

    void swap_points(size_type i, size_type j);
    
    node_tab(scalar_type prec_loose = scalar_type(10000));
    node_tab(const node_tab &t);
    node_tab &operator =(const node_tab &t);
  };


 
}
#endif
