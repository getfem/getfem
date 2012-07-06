/*===========================================================================
 
 Copyright (C) 2004-2012 Julien Pommier
 
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
 
===========================================================================*/


#include "getfem/bgeot_kdtree.h"

namespace bgeot {
  /* generic node for the kdtree */
  struct kdtree_elt_base {
    enum { PTS_PER_LEAF=8 };
    unsigned n; /* 0 => is a tree node, != 0 => tree leaf storing n points */
    bool isleaf() const { return (n != 0); }
    kdtree_elt_base(unsigned n_) : n(n_) {}
  };
  
  /* leafs contains a list of points */
  struct kdtree_leaf : public kdtree_elt_base {
    kdtree_tab_type::const_iterator it;
    kdtree_leaf(kdtree_tab_type::const_iterator begin, 
		kdtree_tab_type::const_iterator end) : 
      kdtree_elt_base(unsigned(std::distance(begin,end))) { it = begin; }
  };
  
  struct kdtree_node : public kdtree_elt_base {
    scalar_type split_v;
    kdtree_elt_base *left, *right; /* left: <=v, right: >v */
    kdtree_node(scalar_type v, kdtree_elt_base *left_, kdtree_elt_base *right_) : 
      kdtree_elt_base(0), split_v(v), left(left_), right(right_) {}
  };

  typedef kdtree_tab_type::iterator ITER;
  
  /* sorting of kdtree_tab_type with respect to a component of the points */
  struct component_sort {
    unsigned dir; 
    component_sort(unsigned d) : dir(d) {}
    bool operator()(const index_node_pair& a, const index_node_pair& b) 
    { return a.n.at(dir) < b.n.at(dir); }
  };
  /* splitting of kdtree_tab_type with respect to a given value */
  /*struct component_split { 
    unsigned dir; scalar_type v;
    component_split(unsigned d, scalar_type v_) : dir(d), v(v_) {}
    bool operator()(const index_node_pair& a) 
    { return (a.n.at(dir) <= v); }
  };
  */

  static ITER partition(ITER begin, ITER end, unsigned dir, scalar_type median) {
    --end;
    do {
      while (begin <= end && (*begin).n.at(dir) <= median) ++begin; 
      while (begin <= end && (*end).n.at(dir) > median) --end;
      if (begin < end) { begin->swap(*end); ++begin; --end; } else break;
    } while (1);
    return begin;
  }
  /* 
     build (recursively) a complete tree for points in the interval
     [begin, end[ dir is the splitting direction.
     (may be ameliorated doing a single sort on each component at the begining)
  */
  static kdtree_elt_base *build_tree_(ITER begin, ITER end, unsigned dir) {
    if (begin == end) return 0;
    size_type npts = std::distance(begin,end);
    if (npts > kdtree_elt_base::PTS_PER_LEAF) {
      ITER itmedian;
      scalar_type median;
      size_type N = begin->n.size();
      unsigned ndir_tests = unsigned(dir/N); dir = unsigned(dir % N);
      if (npts > 50) {
	/* too much points for an exact median: estimation of the median .. */
        std::vector<index_node_pair> v(30);
        //random_sample(begin,end,v.begin(),v.end());
	for (size_type i=0; i < v.size(); ++i) 
	  v[i] = begin[rand() % npts];
        std::sort(v.begin(), v.end(), component_sort(dir));
        median = (v[v.size()/2-1].n.at(dir)+v[v.size()/2].n.at(dir))/2;
        //itmedian = std::partition(begin,end,component_split(dir, median));
	itmedian = partition(begin,end,dir,median);
	/*ITER it = begin; cout << "}"; while (it < end) { 
	  if (it == itmedian) cout << "<|>"; else cout << " ";
	  cout << (*it).n[dir];
	  ++it;
	  } cout << "}\n";*/

      } else {
	/* exact median computation */
        std::sort(begin, end, component_sort(dir));
        itmedian = begin + npts/2 - 1;
        median = (*itmedian).n[dir];
        while (itmedian < end && (*itmedian).n[dir] == median) itmedian++;
      }
      if (itmedian == end) /* could not split the set (all points have same value for component 'dir' !) */
	if (ndir_tests == N-1) /* tested all N direction ? so all points are strictly the same */
	  return new kdtree_leaf(begin,end); 
        else return new kdtree_node(median, 
				    build_tree_(begin, itmedian, unsigned((dir+1)%N + (ndir_tests+1)*N)), 0);
      else { /* the general case */
	assert((*itmedian).n[dir] >= median && median >= (*(itmedian-1)).n[dir]);
	return new kdtree_node(median, 
			       build_tree_(begin, itmedian, unsigned((dir+1)%N)), 
			       build_tree_(itmedian,end, unsigned((dir+1)%N)));
      }
    } else {
      return new kdtree_leaf(begin,end);
    }
  }
  
  static void destroy_tree_(kdtree_elt_base *t) {
    if (t == 0) return;
    if (!t->isleaf()) {
      kdtree_node *tn = static_cast<kdtree_node*>(t);
      destroy_tree_(tn->right);
      destroy_tree_(tn->left);
      delete tn;
    } else {
      kdtree_leaf *tl = static_cast<kdtree_leaf*>(t);
      delete tl;
    }
  }

  /* avoid pushing too much arguments on the stack for points_in_box_ */
  struct points_in_box_data_ {
    base_node::const_iterator bmin;
    base_node::const_iterator bmax;
    kdtree_tab_type *ipts;
    size_type N;
  };

  /* recursive lookup for points inside a given box */
  static void points_in_box_(const points_in_box_data_& p,
			     const kdtree_elt_base *t, unsigned dir) {
    if (!t->isleaf()) {
      const kdtree_node *tn = static_cast<const kdtree_node*>(t);
      if (p.bmin[dir] <= tn->split_v && tn->left)
        points_in_box_(p, tn->left, unsigned((dir+1)%p.N));
      if (p.bmax[dir] > tn->split_v && tn->right)
        points_in_box_(p, tn->right, unsigned((dir+1)%p.N));
    } else {
      const kdtree_leaf *tl = static_cast<const kdtree_leaf*>(t);
      kdtree_tab_type::const_iterator itpt = tl->it;
      for (size_type i=tl->n;i; --i, ++itpt) {
        bool is_in = true;
        base_node::const_iterator it=itpt->n.const_begin();
        for (size_type k=0; k < p.N; ++k) {
	  //cout << "test: k=" << k << ", " << it[k] << ", p.bmin[k]=" << p.bmin[k] << ", p.bmax[k]=" << p.bmax[k] << "\n";
          if (it[k] < p.bmin[k] || it[k] > p.bmax[k]) {
	    is_in = false; break; 
	  }
        }
        if (is_in) p.ipts->push_back(*itpt);
      }
    }
  }

  /* avoid pushing too much arguments on the stack for nearest_neighbor_ */
  struct nearest_neighbor_data_ {
    base_node::const_iterator pos;
    index_node_pair *ipt;
    size_type N;
    mutable scalar_type dist2;
    base_node::iterator vec_to_tree_elm;
  };

  /* recursive lookup for the nearest neighbor point at a given position */
  static void nearest_neighbor_assist(const nearest_neighbor_data_& p,
			              const kdtree_elt_base *t, unsigned dir) {
    scalar_type dist2(0);
    for (size_type k=0; k < p.N; ++k)
      dist2 += p.vec_to_tree_elm[k] * p.vec_to_tree_elm[k];
    if (dist2 > p.dist2 && p.dist2 > scalar_type(0)) return;

    if (!t->isleaf()) {
      const kdtree_node *tn = static_cast<const kdtree_node*>(t);
      scalar_type tmp = p.vec_to_tree_elm[dir];
      scalar_type dist = p.pos[dir] - tn->split_v;
      if (tn->left) {
        if (dist > tmp) p.vec_to_tree_elm[dir] = dist;
        nearest_neighbor_assist(p, tn->left, unsigned((dir+1)%p.N));
        p.vec_to_tree_elm[dir] = tmp;
      }
      if (tn->right) {
        if (-dist > tmp) p.vec_to_tree_elm[dir] = -dist;
        nearest_neighbor_assist(p, tn->right, unsigned((dir+1)%p.N));
        p.vec_to_tree_elm[dir] = tmp;
      }
    } else {
      // find the nearest neighbor inside the leaf
      const kdtree_leaf *tl = static_cast<const kdtree_leaf*>(t);
      kdtree_tab_type::const_iterator itpt = tl->it;
      for (size_type i=tl->n;i; --i, ++itpt) {
        dist2 = scalar_type(0);
        base_node::const_iterator it=itpt->n.const_begin();
        for (size_type k=0; k < p.N; ++k) {
          scalar_type dist = it[k] - p.pos[k];
          dist2 += dist * dist;
        }
        if (dist2 < p.dist2 || p.dist2 < scalar_type(0)){
          *(p.ipt) = *itpt;
          p.dist2 = dist2;
        }
      }
    }
  }

  static void nearest_neighbor_main(const nearest_neighbor_data_& p,
			            const kdtree_elt_base *t, unsigned dir) {
    if (!t->isleaf()) {
      const kdtree_node *tn = static_cast<const kdtree_node*>(t);
      scalar_type dist = p.pos[dir] - tn->split_v;
      if ((dist <= scalar_type(0) && tn->left) || !tn->right) {
        nearest_neighbor_main(p, tn->left, unsigned((dir+1)%p.N));
      } else if (tn->right) {
        nearest_neighbor_main(p, tn->right, unsigned((dir+1)%p.N));
      } else {
        assert(false);
      }
      // check the possibility of points at the opposite side of the current
      // tree node which are closer to pos as the current minimum distance
      if (dist * dist <= p.dist2) {
        for (size_type k=0; k < p.N; ++k) p.vec_to_tree_elm[k] = 0.;
        if ((dist <= scalar_type(0) && tn->right) || !tn->left) {
          p.vec_to_tree_elm[dir] = -dist;
          nearest_neighbor_assist(p, tn->right, unsigned((dir+1)%p.N));
        } else if (tn->left) {
          p.vec_to_tree_elm[dir] = dist;
          nearest_neighbor_assist(p, tn->left, unsigned((dir+1)%p.N));
        }
      }
    } else {
      // find the nearest neighbor inside the leaf which contains pos
      nearest_neighbor_assist(p, t, dir);
    }
  }

  void kdtree::clear_tree() {
    destroy_tree_(tree); tree = 0;
  }

  void kdtree::points_in_box(kdtree_tab_type &ipts,
			     const base_node &min, 
			     const base_node &max) {
    ipts.resize(0);
    if (tree == 0) { tree = build_tree_(pts.begin(), pts.end(), 0); if (!tree) return; }
    base_node bmin(min), bmax(max);
    for (size_type i=0; i < bmin.size(); ++i) if (bmin[i] > bmax[i]) return;
    points_in_box_data_ p; 
    p.bmin = bmin.const_begin(); p.bmax = bmax.const_begin();
    p.ipts = &ipts; p.N = N; 
    points_in_box_(p, tree, 0);
  }

   scalar_type kdtree::nearest_neighbor(index_node_pair &ipt,
                                        const base_node &pos) {

    ipt.i = -1;
    if (tree == 0) {
      tree = build_tree_(pts.begin(), pts.end(), 0);
      if (!tree) return scalar_type(-1);
    }
    nearest_neighbor_data_ p;
    p.pos = pos.const_begin();
    p.ipt = &ipt;
    p.N = N;
    p.dist2 = scalar_type(-1);
    base_node tmp(N);
    p.vec_to_tree_elm = tmp.begin();
    nearest_neighbor_main(p, tree, 0);
    return p.dist2;
  }
}
