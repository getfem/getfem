/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2000-2020 Julien Pommier

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

#ifndef BGEOT_RTREE_H
#define BGEOT_RTREE_H

/** @file bgeot_rtree.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date January 2004.
    @brief region-tree for window/point search on a set of rectangles.
*/

#include <set>
#include "bgeot_small_vector.h"
#include "bgeot_node_tab.h"

namespace bgeot {

  struct box_index {
    size_type id;
    const base_node *min, *max;
  };

  struct box_index_id_compare {
    bool operator()(const box_index *plhs, const box_index *prhs) const {
      return plhs->id < prhs->id;
    }
  };

  struct box_index_topology_compare {
    const scalar_type EPS;
    bool is_less(const base_node &lhs, const base_node &rhs) const {
      GMM_ASSERT2(lhs.size() == rhs.size(), "size mismatch");
      for (size_type i = 0; i < lhs.size(); ++i)
        if (gmm::abs(lhs[i] - rhs[i]) > EPS) {
          return lhs[i] < rhs[i];
        }
      return false;
    }
    
    box_index_topology_compare(scalar_type EPS_) : EPS{EPS_} {}
 
    bool operator()(const box_index &lhs, const box_index &rhs) const {
      if (EPS == scalar_type(0))
        return lhs.id < rhs.id;
      else
        return is_less(*lhs.min, *rhs.min) ||
          (!is_less(*rhs.min, *lhs.min) && is_less(*lhs.max, *rhs.max));
    }
  };

  struct rtree_elt_base {
    enum { RECTS_PER_LEAF=8 };
    bool isleaf_;
    bool isleaf() const { return isleaf_; }
    base_node rmin, rmax;
    rtree_elt_base(bool leaf, const base_node& rmin_, const base_node& rmax_)
      : isleaf_(leaf), rmin(rmin_), rmax(rmax_) {}
    virtual ~rtree_elt_base() {}
  };

  /** Balanced tree of n-dimensional rectangles.
   *
   * This is not a dynamic structure. Once a query has been made on the
   * tree, new boxes should not be added.
   *
   * CAUTION : For EPS > 0, nearly identically boxes are eliminated
   *           For EPS = 0 all boxes are stored.
   */
  class rtree {
  public:
    using box_cont = std::set<box_index,box_index_topology_compare> ;
    using pbox_cont = std::vector<const box_index*>;
    using pbox_set = std::set<const box_index *, box_index_id_compare>;

    rtree(scalar_type EPS = 0);
    rtree(const rtree&) = delete;
    rtree& operator = (const rtree&) = delete;

    size_type add_box(const base_node &min, const base_node &max,
                      size_type id=size_type(-1));
    size_type nb_boxes() const { return boxes.size(); }
    void clear();

    void find_intersecting_boxes(const base_node& bmin, const base_node& bmax,
                                 pbox_set& boxlst) const;
    void find_containing_boxes(const base_node& bmin, const base_node& bmax,
                               pbox_set& boxlst) const;
    void find_contained_boxes(const base_node& bmin, const base_node& bmax,
                              pbox_set& boxlst) const;
    void find_boxes_at_point(const base_node& P, pbox_set& boxlst) const;
    void find_line_intersecting_boxes(const base_node& org,
                                      const base_small_vector& dirv,
                                      pbox_set& boxlst) const;
    void find_line_intersecting_boxes(const base_node& org,
                                      const base_small_vector& dirv,
                                      const base_node& bmin,
                                      const base_node& bmax,
                                      pbox_set& boxlst) const;

    void find_intersecting_boxes(const base_node& bmin, const base_node& bmax,
                                 std::vector<size_type>& idvec) {
      pbox_set bs;
      find_intersecting_boxes(bmin, bmax, bs);
      pbox_set_to_idvec(bs, idvec);
    }
    void find_containing_boxes(const base_node& bmin, const base_node& bmax,
                               std::vector<size_type>& idvec) {
      pbox_set bs;
      find_containing_boxes(bmin, bmax, bs);
      pbox_set_to_idvec(bs, idvec);
    }
    void find_contained_boxes(const base_node& bmin,
                              const base_node& bmax,
                              std::vector<size_type>& idvec) {
      pbox_set bs;
      find_contained_boxes(bmin, bmax, bs);
      pbox_set_to_idvec(bs, idvec);
    }
    void find_boxes_at_point(const base_node& P,
                             std::vector<size_type>& idvec) const
    { pbox_set bs; find_boxes_at_point(P, bs);  pbox_set_to_idvec(bs, idvec); }
    void find_line_intersecting_boxes(const base_node& org,
                                      const base_small_vector& dirv,
                                      std::vector<size_type>& idvec) {
      pbox_set bs;
      find_line_intersecting_boxes(org, dirv, bs);
      pbox_set_to_idvec(bs, idvec);
    }
    void find_line_intersecting_boxes(const base_node& org,
                                      const base_small_vector& dirv,
                                      const base_node& bmin,
                                      const base_node& bmax,
                                      std::vector<size_type>& idvec) {
      pbox_set bs;
      find_line_intersecting_boxes(org, dirv, bmin, bmax, bs);
      pbox_set_to_idvec(bs, idvec);
    }

    void dump();
    void build_tree();
  private:
    void pbox_set_to_idvec(pbox_set bs, std::vector<size_type>& idvec) const {
      idvec.reserve(bs.size()); idvec.resize(0);
      for (pbox_set::const_iterator it=bs.begin(); it != bs.end(); ++it)
        idvec.push_back((*it)->id);
    }

    const scalar_type EPS;
    node_tab nodes;
    box_cont boxes;
    std::unique_ptr<rtree_elt_base> root;
    bool tree_built;
    getfem::lock_factory locks_;
  };

}

#endif
