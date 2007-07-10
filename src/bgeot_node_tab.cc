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


#include "getfem/bgeot_node_tab.h"

namespace bgeot {

  bool node_tab::component_comp::operator()(size_type i1, size_type i2) const {
    scalar_type a = (i1==size_type(-1)) ? *c : gmm::vect_sp((*vbn)[i1], v);
    scalar_type b = (i2==size_type(-1)) ? *c : gmm::vect_sp((*vbn)[i2], v);
    if (a == b) return (i1 < i2); else return a < b;
  }

  node_tab::component_comp::component_comp
  (const dal::dynamic_tas<base_node> &vbn_,
   const scalar_type &c_, unsigned dim)
    : vbn(&vbn_), c(&c_), v(dim) {
    do gmm::fill_random(v); while (gmm::vect_norm2(v) == 0);
    gmm::scale(v, scalar_type(1) / gmm::vect_norm2(v) );
  }

  void node_tab::add_sorter(void) const {
    if (sorters.size() > 0) GMM_WARNING3("Multiple sort need for node tab");
    sorters.push_back(sorter(component_comp(*this, c, dim_)));
    for (dal::bv_visitor i(index()); !i.finished(); ++i)
      sorters.back().insert(size_type(i));
  }
  
  size_type node_tab::search_node(const base_node &pt) const {
    for (size_type is = 0; ; ++is) {
      if (is >= sorters.size()) add_sorter();
      
      c = gmm::vect_sp(pt, sorters[is].key_comp().v) - eps;
      sorter::const_iterator it = sorters[is].lower_bound(size_type(-1));
      size_type count = 0;
      for (; it != sorters[is].end() && count < 20; ++it, ++count) {
	const base_node &pt2 = (*this)[*it];
	if (gmm::vect_dist2(pt, pt2) < eps) return *it;
	if (gmm::vect_sp(pt2, sorters[is].key_comp().v) > c+eps+eps)
	  return size_type(-1);
      }
      if (it == sorters[is].end()) return size_type(-1);
      GMM_ASSERT1(is < 10, "Problem in node structure");
    }
  }
  
  void node_tab::clear(void) {
    dal::dynamic_tas<base_node>::clear();
    sorters = std::vector<sorter>();
    max_radius = scalar_type(1e-60);
    eps = max_radius * prec_factor;
  }
  
  size_type node_tab::add_node(const base_node &pt) {
    scalar_type npt = gmm::vect_norm2(pt);
    max_radius = std::max(max_radius, npt);
    eps = max_radius * prec_factor;
    
    if (this->card() == 0) {
      dim_ = pt.size();
      return this->add(pt);
    }
    else {
      GMM_ASSERT1(dim_ == pt.size(), "Nodes should have the same dimension");
      size_type id = search_node(pt);
      if (id == size_type(-1)) {
	id = this->add(pt);
	for (size_type i = 0; i < sorters.size(); ++i) sorters[i].insert(id);
      }
      return id;
    }
  }

  void node_tab::swap_points(size_type i, size_type j) {
    if (i != j) {
      bool existi = index().is_in(i);
      bool existj = index().is_in(j);
      for (size_type is = 0; is < sorters.size(); ++is) {
	if (existi) sorters[is].erase(i);
	if (existj) sorters[is].erase(j);
      }
      dal::dynamic_tas<base_node>::swap(i, j);
      for (size_type is = 0; is < sorters.size(); ++is) {
	if (existi) sorters[is].insert(j);
	if (existj) sorters[is].insert(i);
      }
    }
  }

  void node_tab::translation(const base_small_vector &V) {
    for (dal::bv_visitor i(index()); !i.finished(); ++i)
      (*this)[i] += V;
    sorters = std::vector<sorter>();
  }

  void node_tab::transformation(const base_matrix &M) {
    base_small_vector w(M.nrows());
    GMM_ASSERT1(gmm::mat_nrows(M) != 0 && gmm::mat_ncols(M) == dim(),
		"invalid dimensions for the transformation matrix");
    dim_ = gmm::mat_nrows(M);
    for (dal::bv_visitor i(index()); !i.finished(); ++i) {
      w = (*this)[i];
      gmm::resize((*this)[i], dim_);
      gmm::mult(M,w,(*this)[i]);
    }
    sorters = std::vector<sorter>();
  }
  
  
  node_tab::node_tab(scalar_type prec_loose) {
    max_radius = scalar_type(1e-60);
    sorters.reserve(5);
    prec_factor = gmm::default_tol(scalar_type()) * prec_loose;
    eps = max_radius * prec_factor;
  }

 
}
