/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_mat_elem_type.C : precomputations on fem              */
/*              interpolations.                                            */
/*                                                                         */
/* Date : December 21, 2000.                                               */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2002  Yves Renard.                                   */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#include <getfem_mat_elem_type.h>
#include <dal_tree_sorted.h>

namespace getfem {

  bool operator < (const constituant &m, const constituant &n) {
    if (m.t < n.t) return true; if (m.t > n.t) return false;
    if (m.pfi < n.pfi) return true; if (m.pfi > n.pfi) return false;
    if (m.nlt < n.nlt) return true; if (m.nlt > n.nlt) return false;
    return false;
  }

  static pmat_elem_type add_to_met_tab(const mat_elem_type &f) {
    static dal::dynamic_tree_sorted<mat_elem_type,
      dal::lexicographical_less<mat_elem_type> > *tab;
    static bool isinit = false;
    if (!isinit) {
      tab = new dal::dynamic_tree_sorted<mat_elem_type,
	dal::lexicographical_less<mat_elem_type> >();
      isinit = true;
    }
    return &((*tab)[tab->add_norepeat(f)]);
  }

  pmat_elem_type mat_elem_base(pfem pfi) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_BASE_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1)
    { f.mi.resize(1); f.mi[0] = pfi->nb_base(); }
    else
    { f.mi.resize(2); f.mi[0] = pfi->nb_base(); f.mi[1] = pfi->target_dim(); }
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_grad(pfem pfi) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_GRAD_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1) { 
      f.mi.resize(2); f.mi[0] = pfi->nb_base();
      f.mi[1] = pfi->structure()->dim();
    }
    else {
      f.mi.resize(3); f.mi[0] = pfi->nb_base();
      f.mi[1] = pfi->target_dim(); f.mi[2] = pfi->structure()->dim();
    }
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_hessian(pfem pfi) {
    mat_elem_type f; f.resize(1);  f[0].t = GETFEM_HESSIAN_; f[0].pfi = pfi;
    f[0].nlt = 0;
    if (pfi->target_dim() == 1) { 
      f.mi.resize(2); f.mi[0] = pfi->nb_base();
      f.mi[1] = dal::sqr(pfi->structure()->dim());
    }
    else {
      f.mi.resize(3); f.mi[0] = pfi->nb_base();
      f.mi[1] = pfi->target_dim();
      f.mi[2] = dal::sqr(pfi->structure()->dim());
    }
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_nonlinear(const nonlinear_elem_term &nlt, pfem pfi) {
    mat_elem_type f; f.resize(1); f[0].t = GETFEM_NONLINEAR_;
    f[0].pfi = pfi;
    f[0].nlt = &nlt;
    f.mi = nlt.sizes();
    return add_to_met_tab(f);
  }

  pmat_elem_type mat_elem_product(pmat_elem_type a, pmat_elem_type b) {
    mat_elem_type f; f.resize(a->size() + b->size());
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
    return add_to_met_tab(f);
  }


}  /* end of namespace getfem.                                            */

