// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2000-2007 Julien Pommier
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

#include "getfem/bgeot_rtree.h"

namespace bgeot {
  struct rtree_elt_base {
    enum { RECTS_PER_LEAF=8 };
    bool isleaf_;
    bool isleaf() const { return isleaf_; }
    base_node rmin, rmax;
    rtree_elt_base(bool leaf, const base_node& rmin_, const base_node& rmax_) 
      : isleaf_(leaf), rmin(rmin_), rmax(rmax_) {}
  };
  
  struct rtree_node : public rtree_elt_base {
    rtree_elt_base *left, *right;
    rtree_node(const base_node& bmin, const base_node& bmax, 
	       rtree_elt_base *left_, rtree_elt_base *right_) 
      : rtree_elt_base(false, bmin, bmax), left(left_), right(right_) {}
  };

  struct rtree_leaf : public rtree_elt_base {
    rtree::pbox_cont lst;
    rtree_leaf(const base_node& bmin, const base_node& bmax, 
	       rtree::pbox_cont& lst_)
      : rtree_elt_base(true, bmin, bmax) { lst.swap(lst_); }
  };

  /* enlarge box to hold [a..b] */
  static void update_box(base_node& bmin, base_node& bmax, const base_node& a, const base_node& b) {
    base_node::iterator itmin=bmin.begin(), itmax=bmax.begin();
    for (size_type i=0; i < a.size(); ++i) {
      itmin[i] = std::min(itmin[i], a.at(i));
      itmax[i] = std::max(itmax[i], b.at(i));
    }
  }

  static bool r1_ge_r2(const base_node& min1, const base_node& max1,
		       const base_node& min2, const base_node& max2) {
    for (size_type i=0; i < min1.size(); ++i)
      if (!(min1[i] <= min2[i] && max1[i] >= max2[i])) return false; 
    return true;
  }

  static bool r1_inter_r2(const base_node& min1, const base_node& max1,
			  const base_node& min2, const base_node& max2) {
    for (size_type i=0; i < min1.size(); ++i) 
      if (max1[i] < min2[i] || min1[i] > max2[i]) return false; 
    return true;
  }

  /* some predicates for searches */
  struct intersection_p {
    const base_node min,max;
    intersection_p(const base_node& min_, const base_node& max_) : min(min_), max(max_) {}
    bool operator()(const base_node& min2, const base_node& max2)
    { return r1_inter_r2(min,max,min2,max2); }
    bool accept(const base_node& min2, const base_node& max2) 
    { return operator()(min2,max2); }
  };

  /* match boxes containing [min..max] */
  struct contains_p {
    const base_node min,max;
    contains_p(const base_node& min_, const base_node& max_) : min(min_), max(max_) {}
    bool operator()(const base_node& min2, const base_node& max2)
    { return r1_ge_r2(min2,max2,min,max); }
    bool accept(const base_node& min2, const base_node& max2) 
    { return r1_inter_r2(min,max,min2,max2); }
  };

  /* match boxes contained in [min..max] */
  struct contained_p {
    const base_node min,max;
    contained_p(const base_node& min_, const base_node& max_) : min(min_), max(max_) {}
    bool accept(const base_node& min2, const base_node& max2)
    { return r1_inter_r2(min,max,min2,max2); }
    bool operator()(const base_node& min2, const base_node& max2) 
    { return r1_ge_r2(min,max,min2,max2); }
  };

  struct has_point_p {
    const base_node P;
    has_point_p(const base_node& P_) : P(P_) {}
    bool operator()(const base_node& min2, const base_node& max2) {
      for (size_type i=0; i < P.size(); ++i) 
	if (P[i] < min2[i] || P[i] > max2[i]) return false;
      return true;
    }
    bool accept(const base_node& min2, const base_node& max2) 
    { return operator()(min2,max2); }
  };

  template <typename Predicate>
  static void find_matching_boxes_(rtree_elt_base *n, rtree::pbox_set& boxlst,
				   Predicate p) {
    //cout << "find_matching_boxes_: " << n->rmin << ".." << n->rmax << "\n";
    if (n->isleaf()) {
      //cout << "find_matching_boxes_ in leaf\n";
      const rtree_leaf *rl = static_cast<rtree_leaf*>(n);
      for (rtree::pbox_cont::const_iterator it = rl->lst.begin(); it != rl->lst.end(); ++it) {
	//cout << "  ->match(" << (*it)->id << "=" << (*it)->min << "," << (*it)->max << " -> " << p((*it)->min, (*it)->max) << "\n";
	if (p((*it)->min, (*it)->max)) { boxlst.insert(*it); }
      }
    } else {
      //cout << "find_matching_boxes_ in branch\n";
      const rtree_node *rn = static_cast<rtree_node*>(n);
      if (p.accept(rn->left->rmin,rn->left->rmax)) 
	bgeot::find_matching_boxes_(rn->left, boxlst, p);
      if (p.accept(rn->right->rmin,rn->right->rmax)) 
	bgeot::find_matching_boxes_(rn->right, boxlst, p);
    }
  }

  void rtree::find_intersecting_boxes(const base_node& bmin, const base_node& bmax, pbox_set& boxlst) {
    boxlst.clear(); if (!root) build_tree(); 
    if (root) find_matching_boxes_(root, boxlst, intersection_p(bmin,bmax));
    //cout << "find_intersecting_boxes : found " << boxlst.size() << " matches\n";
  }

  void rtree::find_containing_boxes(const base_node& bmin, const base_node& bmax, pbox_set& boxlst) {
    boxlst.clear(); if (!root) build_tree();
    if (root) find_matching_boxes_(root, boxlst, contains_p(bmin,bmax));
  }

  void rtree::find_contained_boxes(const base_node& bmin, const base_node& bmax, pbox_set& boxlst) {
    boxlst.clear(); if (!root) build_tree();
    if (root) find_matching_boxes_(root, boxlst, contained_p(bmin,bmax));
    //cout << "find_matching_boxes : found " << boxlst.size() << " matches\n";
  }

  void rtree::find_boxes_at_point(const base_node& P, pbox_set& boxlst) {
    boxlst.clear(); if (!root) build_tree();
    if (root) find_matching_boxes_(root, boxlst, has_point_p(P));
  }

  /* 
     try to split at the approximate center of the box. Could be much more
     sophisticated 
  */
  static bool split_test(const rtree::pbox_cont& b, 
			 const base_node& bmin, const base_node& bmax, 
			 unsigned dir, scalar_type& split_v) {
    scalar_type v = bmin[dir] + (bmax[dir] - bmin[dir])/2; split_v = v;
    size_type cnt = 0;
    //cout << "[enter]Split_test: dir=" << dir << ", split_v=" << v << ", bmin=" << bmin << ", bmax=" << bmax << "\n";
    for (rtree::pbox_cont::const_iterator it = b.begin(); it != b.end(); ++it) {
      if ((*it)->max[dir] < v) {
	if (cnt == 0) split_v = (*it)->max[dir]; 
	else split_v = std::max((*it)->max[dir],split_v);
	cnt++; 
      }
    }
    //cout << "[exit] Split_test cnt = " << cnt << ", b.size()=" << b.size() << ", split_v=" << split_v << "\n";
    return (cnt > 0 && cnt < b.size());
  }

  /* 
     there are many flavors of rtree ... this one is more or less a quadtree
     where splitting does not occurs at predefined positions (hence the split_test
     function above). Regions of the tree do not overlap (box are splitted).
  */
  static rtree_elt_base* build_tree_(rtree::pbox_cont b, 
				     const base_node& bmin, const base_node& bmax,
				     unsigned last_dir) {
    size_type N=bmin.size();
    scalar_type split_v; 
    unsigned split_dir = (last_dir+1)%N;
    //cout << " build_tree_ [b.size=" << b.size() << "], bmin=" << bmin << ", bmax=" << bmax << "\n";
    bool split_ok = false;
    if (b.size() > rtree_elt_base::RECTS_PER_LEAF) {
      for (size_type itry=0; itry < N; ++itry) {
	//cout << "split_test: dir=" << split_dir << "\n";
	if (split_test(b, bmin, bmax, split_dir, split_v)) { split_ok = true; break; }
	split_dir = (split_dir+1)%N;
      }
      //if (!split_ok && b.size() > rtree_elt_base::RECTS_PER_LEAF*2) cout << "FAILED TO SPLIT ...\n";
    }
    if (split_ok) {
      size_type cnt1=0,cnt2=0;
      //cout << "splitting with v=" << split_v << "\n";
      for (rtree::pbox_cont::const_iterator it = b.begin(); it != b.end(); ++it) {
	//cout << " . test box" << (*it)->min[split_dir] << ".." << (*it)->max[split_dir] << "\n";
	if ((*it)->min[split_dir] < split_v) cnt1++; 
	if ((*it)->max[split_dir] > split_v) cnt2++;
      }
      //cout << "  -> left : " << cnt1 << " boxes, right : " << cnt2 << " boxes\n";
      assert(cnt1); assert(cnt2); assert(cnt1+cnt2 >= b.size());
      rtree::pbox_cont v1(cnt1), v2(cnt2);
      base_node bmin1(bmax), bmax1(bmin); 
      base_node bmin2(bmax), bmax2(bmin);
      cnt1 = cnt2 = 0;
      for (rtree::pbox_cont::const_iterator it = b.begin(); it != b.end(); ++it) {
	if ((*it)->min[split_dir] < split_v) {
	  v1[cnt1++] = *it; 
	  update_box(bmin1,bmax1,(*it)->min,(*it)->max);
	  //cout << "update_box bmin1=" << bmin1 << ", bmax1=" << bmax1 << "\n";
	}
	if ((*it)->max[split_dir] > split_v) {
	  v2[cnt2++] = *it;
	  update_box(bmin2,bmax2,(*it)->min,(*it)->max);
	}
      }
      for (size_type k=0; k < N; ++k) {
	bmin1[k] = std::max(bmin1[k],bmin[k]);
	bmax1[k] = std::min(bmax1[k],bmax[k]);
	bmin2[k] = std::max(bmin2[k],bmin[k]);
	bmax2[k] = std::min(bmax2[k],bmax[k]);
      }
      bmax1[split_dir] = std::min(bmax1[split_dir], split_v);
      bmin2[split_dir] = std::max(bmin2[split_dir], split_v);
      assert(cnt1 == v1.size()); assert(cnt2 == v2.size());
      return new rtree_node(bmin,bmax, 
			    build_tree_(v1, bmin1, bmax1, split_dir),
			    build_tree_(v2, bmin2, bmax2, split_dir));
    } else {
      return new rtree_leaf(bmin,bmax,b);
    }
  }

 void rtree::build_tree() {
   //cout << "build tree\n";
   if (boxes.size() == 0) return;
   assert(root == 0);
   pbox_cont b(boxes.size());
   pbox_cont::iterator b_it = b.begin();
   base_node bmin(boxes.front().min), bmax(boxes.front().max);
   for (box_cont::const_iterator it=boxes.begin(); it != boxes.end(); ++it) {
     update_box(bmin,bmax,(*it).min,(*it).max);
     *b_it++ = &(*it);
   }
   root = build_tree_(b,bmin,bmax,0);
 }
  
  static void destroy_tree_(rtree_elt_base *n) {
    if (n->isleaf()) delete static_cast<rtree_leaf*>(n);
    else {
      const rtree_node *rn = static_cast<rtree_node*>(n);
      if (rn->left) { destroy_tree_(rn->left); }
      if (rn->right) { destroy_tree_(rn->right); }
      delete rn;
    }
  }

  void rtree::destroy_tree() {
    if (root) { destroy_tree_(root); }
    root = 0;
  }

  static void dump_tree_(rtree_elt_base *p, int level, size_type& count) {
    if (!p) return;
    cout << level << "|";
    for (int i=0; i < level; ++i) cout << "  ";
    cout << "span=" << p->rmin << ".." << p->rmax << " ";
    if (p->isleaf()) {
      rtree_leaf *rl = static_cast<rtree_leaf*>(p);
      cout << "Leaf [" << rl->lst.size() << " elts] = ";
      for (size_type i=0; i < rl->lst.size(); ++i) cout << " " << rl->lst[i]->id; cout << "\n";
      count += rl->lst.size();
    } else {
      cout << "Node\n";
      const rtree_node *rn = static_cast<rtree_node*>(p);
      if (rn->left) { dump_tree_(rn->left, level+1, count); }
      if (rn->right) { dump_tree_(rn->right, level+1, count); }
    }
  }

  void rtree::dump() {
    cout << "tree dump follows\n";
    if (!root) build_tree();
    size_type count = 0;
    dump_tree_(root, 0, count);
    cout << " --- end of tree dump, nb of rectangles: " << boxes.size() 
	 << ", rectangle ref in tree: " << count << "\n";
  }
}
