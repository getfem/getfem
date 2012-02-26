/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
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
      const base_node *c;
      base_small_vector v;
      bool operator()(size_type i1, size_type i2) const;
      component_comp(const dal::dynamic_tas<base_node> &vbn_,
		     const base_node &c_, unsigned dim);
	  component_comp() : vbn(0), c(0) {}
    };
    typedef std::set<size_type, component_comp> sorter;

    mutable std::vector<sorter> sorters;
    mutable base_node c;
    scalar_type eps, prec_factor, max_radius;
    unsigned dim_;

    void add_sorter(void) const;

  public :

    /// reset the array, remove all points
    void clear(void);

    /** Search a node in the array. return its index if it exists
	or size_type(-1) otherwise.
    */
    size_type search_node(const base_node &pt, const scalar_type radius=0) const;
    /** Add a point to the array or identify it with a very close existing
	point.
    */
    size_type add_node(const base_node &pt);
    size_type add(const base_node &pt) { return add_node(pt); }
    void sup_node(size_type i);
    void sup(size_type i) { sup_node(i); }
    void resort(void) { sorters = std::vector<sorter>(); }
    dim_type dim(void) const { return dim_type(dim_); }
    void translation(const base_small_vector &V);
    void transformation(const base_matrix &M);

    void swap_points(size_type i, size_type j);
    void swap(size_type i, size_type j) { swap_points(i,j); }

    node_tab(scalar_type prec_loose = scalar_type(10000));
    node_tab(const node_tab &t);
    node_tab &operator =(const node_tab &t);
  };



}
#endif
