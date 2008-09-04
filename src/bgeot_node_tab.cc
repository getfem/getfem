// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================


#include "getfem/bgeot_node_tab.h"

namespace bgeot {
  
  bool node_tab::component_comp::operator()(size_type i1,size_type i2) const {
    if (i1 == i2) return false;
    const base_node &pt1((i1 == size_type(-1)) ? *c : (*vbn)[i1]);
    const base_node &pt2((i2 == size_type(-1)) ? *c : (*vbn)[i2]);
    unsigned d = pt1.size();
    base_small_vector::const_iterator it = v.begin();
    base_node::const_iterator it1 = pt1.begin(), it2 = pt2.begin();
    scalar_type a(0);
    for (size_type i = 0; i < d; ++i) a += (*it++) * (*it1++ - *it2++);
    if (a != scalar_type(0)) return a < 0;
    if (i1 == size_type(-1) || i2 == size_type(-1)) return false;
    return i1 < i2;
  }

  node_tab::component_comp::component_comp
  (const dal::dynamic_tas<base_node> &vbn_, const base_node &c_, unsigned dim)
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
    if (card() == 0) return size_type(-1);
    for (size_type is = 0; ; ++is) {
      if (is >= sorters.size()) add_sorter();
      c = pt - eps * sorters[is].key_comp().v;
      sorter::const_iterator it = sorters[is].lower_bound(size_type(-1));
      scalar_type up_bound
	= gmm::vect_sp(pt, sorters[is].key_comp().v) + eps + eps;
      size_type count = 0;
      for (; it != sorters[is].end() && count < 20; ++it, ++count) {
	const base_node &pt2 = (*this)[*it];
	if (gmm::vect_dist2(pt, pt2) < eps) return *it;
	if (gmm::vect_sp(pt2, sorters[is].key_comp().v) > up_bound)
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
    size_type id;
    if (this->card() == 0) {
      dim_ = pt.size();
      id = dal::dynamic_tas<base_node>::add(pt);
      for (size_type i = 0; i < sorters.size(); ++i) sorters[i].insert(id);
    }
    else {
      GMM_ASSERT1(dim_ == pt.size(), "Nodes should have the same dimension");
      id = search_node(pt);
      if (id == size_type(-1)) {
	id = dal::dynamic_tas<base_node>::add(pt);
	for (size_type i = 0; i < sorters.size(); ++i) {
	  sorters[i].insert(id);
	  GMM_ASSERT3(sorters[i].size() == card(), "internal error");
	}
      }
    }
    return id;
  }

  void node_tab::swap_points(size_type i, size_type j) {
    if (i != j) {
      bool existi = index().is_in(i), existj = index().is_in(j);
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

  void node_tab::sup_node(size_type i) {
    if (index().is_in(i)) {
      for (size_type is = 0; is < sorters.size(); ++is) sorters[is].erase(i);
      dal::dynamic_tas<base_node>::sup(i);
    }
  }

  void node_tab::translation(const base_small_vector &V) {
    for (dal::bv_visitor i(index()); !i.finished(); ++i) (*this)[i] += V;
    resort();
  }
  
  void node_tab::transformation(const base_matrix &M) {
    base_small_vector w(M.nrows());
    GMM_ASSERT1(gmm::mat_nrows(M) != 0 && gmm::mat_ncols(M) == dim(),
		"invalid dimensions for the transformation matrix");
    dim_ = unsigned(gmm::mat_nrows(M));
    for (dal::bv_visitor i(index()); !i.finished(); ++i) {
      w = (*this)[i];
      gmm::resize((*this)[i], dim_);
      gmm::mult(M,w,(*this)[i]);
    }
    resort();
  }
  
  node_tab::node_tab(scalar_type prec_loose) {
    max_radius = scalar_type(1e-60);
    sorters.reserve(5);
    prec_factor = gmm::default_tol(scalar_type()) * prec_loose;
    eps = max_radius * prec_factor;
  }

  node_tab::node_tab(const node_tab &t)
    : dal::dynamic_tas<base_node>(t), sorters(), eps(t.eps),
      prec_factor(t.prec_factor), max_radius(t.max_radius), dim_(t.dim_)  {}

  node_tab &node_tab::operator =(const node_tab &t) {
    dal::dynamic_tas<base_node>::operator =(t);
    sorters = std::vector<sorter>();
    eps = t.eps; prec_factor = t.prec_factor;
    max_radius = t.max_radius; dim_ = t.dim_;
    return *this;
  }

}
