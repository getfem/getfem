// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mat_elem_type.cc : precomputations on fem
//           interpolations.
// Date    : December 21, 2000.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Yves Renard
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

#include <dal_singleton.h>
#include <getfem_mat_elem_type.h>
#include <dal_tree_sorted.h>
#include <getfem_mat_elem.h> /* for mat_elem_forget_mat_elem_type */
namespace getfem {

  bool operator < (const constituant &m, const constituant &n) {
    if (m.t < n.t) return true; if (m.t > n.t) return false;
    if (m.pfi < n.pfi) return true; if (m.pfi > n.pfi) return false;
    if (m.nlt < n.nlt) return true; if (m.nlt > n.nlt) return false;
    return false;
  }

  typedef dal::dynamic_tree_sorted<mat_elem_type,
    dal::lexicographical_less<mat_elem_type> > mat_elem_type_tab;

  static pmat_elem_type add_to_met_tab(const mat_elem_type &f) {
    mat_elem_type_tab &tab = dal::singleton<mat_elem_type_tab>::instance();
    return &(tab[tab.add_norepeat(f)]);
  }

  /* on destruction, all occurences of the nonlinear term are removed
     from the mat_elem_type cache; */
  nonlinear_elem_term::~nonlinear_elem_term() {
    mat_elem_type_tab &tab = dal::singleton<mat_elem_type_tab>::instance();
    for (dal::bv_visitor_c i(tab.index()); !i.finished(); ++i) {
      for (size_type j=0; j < tab[i].size(); ++j) {
	if (tab[i][j].t == GETFEM_NONLINEAR_ && tab[i][j].nlt == this) {
	  /* remove all mat_elem structures pointing to the mat_elem_type */
	  mat_elem_forget_mat_elem_type(&tab[i]); 
	  tab.sup(i);
	  break;
	}
      }
    }
    /*
    for (mat_elem_type_tab::sorted_iterator it = tab.begin(); it != tab.end();) {
      mat_elem_type_tab::sorted_iterator itnext = it+1;
      if ((*it).t = GETFEM_NONLINEAR_ && ((*it).nlt == this)) {
	mat_elem_forget_mat_elem_type(&(*it)); 
	tab.sup((*it).index());
      }
      it = itnext;
    }*/
  }

  pmat_elem_type mat_elem_base(pfem pfi) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_BASE_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1)
      { f.get_mi().resize(1); f.get_mi()[0] = 1; }
    else {
      f.get_mi().resize(2); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->target_dim();
    }
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_unit_normal(void) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_UNIT_NORMAL_;
    f[0].pfi = 0; f[0].nlt = 0;
    f.get_mi().resize(1); f.get_mi()[0] = 1;
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_grad(pfem pfi) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_GRAD_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1) { 
      f.get_mi().resize(2); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->dim();
    }
    else {
      f.get_mi().resize(3); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->target_dim();
      f.get_mi()[2] = pfi->dim();
    }
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_hessian(pfem pfi) {
    mat_elem_type f; f.resize(1);  f[0].t = GETFEM_HESSIAN_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1) { 
      f.get_mi().resize(2); f.get_mi()[0] = 1;
      f.get_mi()[1] = gmm::sqr(pfi->dim());
    }
    else {
      f.get_mi().resize(3); f.get_mi()[0] = 1;
      f.get_mi()[1] = pfi->target_dim();
      f.get_mi()[2] = gmm::sqr(pfi->dim());
    }
    return add_to_met_tab(f);
  }

  static pmat_elem_type mat_elem_nonlinear_(pnonlinear_elem_term nlt,
					    pfem pfi, unsigned nl_part) {
    mat_elem_type f; f.resize(1); 
    f[0].t = GETFEM_NONLINEAR_; f[0].nl_part = nl_part;
    f[0].pfi = pfi;
    f[0].nlt = nlt;
    if (nl_part) {
      f.get_mi().resize(1); f.get_mi()[0] = 1;
    } else f.get_mi() = nlt->sizes();
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_nonlinear(pnonlinear_elem_term nlt,
				    std::vector<pfem> pfi) {
    if (pfi.size() == 0)
      DAL_THROW(dal::dimension_error, "mat_elem_nonlinear with no pfem!");
    pmat_elem_type me = mat_elem_nonlinear_(nlt, pfi[0], 0);
    for (size_type i=1; i < pfi.size(); ++i)
      me = mat_elem_product(mat_elem_nonlinear_(nlt, pfi[i], i),me);
    return me;
  }

  pmat_elem_type mat_elem_product(pmat_elem_type a, pmat_elem_type b) {
    mat_elem_type f; f.reserve(a->size() + b->size());
    f.get_mi().reserve(a->get_mi().size() + b->get_mi().size());
    f.insert(f.end(), (*a).begin(), (*a).end());
    f.insert(f.end(), (*b).begin(), (*b).end());
    f.get_mi().insert(f.get_mi().end(), (*a).get_mi().begin(), (*a).get_mi().end());
    f.get_mi().insert(f.get_mi().end(), (*b).get_mi().begin(), (*b).get_mi().end());

    /*    mat_elem_type f; f.resize(a->size() + b->size());
	  f.mi.resize(a->mi.size() + b->mi.size());
    mat_elem_type::const_iterator ita = a->begin(), itae = a->end();
    mat_elem_type::const_iterator itb = b->begin(), itbe = b->end(), it;
    mat_elem_type::iterator itf = f.begin();
    bgeot::multi_index::const_iterator itma = a->mi.begin();
    bgeot::multi_index::const_iterator itmb = b->mi.begin(), *itm;
    bgeot::multi_index::iterator itmf = f.mi.begin();
    for( ;  ita != itae || itb != itbe; ++itf ) {
      if (ita == itae)      { it = itb; ++itb; itm = &(itmb); }
      else                  { it = ita; ++ita; itm = &(itma); }
     
      *itf = *it;
      switch ((*it).t) { 
      case GETFEM_BASE_      : *itmf++ = *(*itm)++; break;
      case GETFEM_GRAD_      : *itmf++ = *(*itm)++; *itmf++ = *(*itm)++; break;
      case GETFEM_HESSIAN_   : *itmf++ = *(*itm)++; *itmf++ = *(*itm)++; break;
      case GETFEM_NONLINEAR_ :
	for (dim_type i = 0; i < (*it).nlt->sizes().size(); ++i) *itmf++ = *(*itm)++;
	break;
      }
    }
    */
    return add_to_met_tab(f);
  }

  bgeot::multi_index mat_elem_type::sizes(size_type cv) const {
    bgeot::multi_index mii = mi;
    for (size_type i = 0, j = 0; i < size(); ++i, ++j) {
      switch ((*this)[i].t) {
      case GETFEM_BASE_ :
	mii[j] = (*this)[i].pfi->nb_base(cv);
	if ((*this)[i].pfi->target_dim() != 1) ++j;
	break;
      case GETFEM_GRAD_ :
	mii[j] = (*this)[i].pfi->nb_base(cv); ++j;
	if ((*this)[i].pfi->target_dim() != 1) ++j;
	break;     
      case GETFEM_HESSIAN_   :
	mii[j] = (*this)[i].pfi->nb_base(cv); ++j;
	if ((*this)[i].pfi->target_dim() != 1) ++j;
	break;
      case GETFEM_UNIT_NORMAL_ :
	break;
      case GETFEM_NONLINEAR_ :
	if ((*this)[i].nl_part == 0)
	  { j+=(*this)[i].nlt->sizes().size(); --j; }
	break;
      }
    }
    return mii;
  }


}  /* end of namespace getfem.                                            */

