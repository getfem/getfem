/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_precomp.C : pre-computations.                         */
/*     									   */
/*                                                                         */
/* Date : June 17, 2002.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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

#include <dal_singleton.h>
#include <getfem_mat_elem.h>
#include <getfem_precomp.h>

namespace getfem
{
  /* ********************************************************************* */
  /*       Precomputation on fem.                                          */
  /* ********************************************************************* */

  struct pre_fem_light_
  {
    pfem pf;
    bgeot::pstored_point_tab pspt;
    fem_precomp_pool *pool;
    bool operator < (const pre_fem_light_ &ls) const
    {
      if (pf < ls.pf) return true; if (pf > ls.pf) return false; 
      if (pspt < ls.pspt) return true; return false;
    }
    pre_fem_light_(pfem pff, bgeot::pstored_point_tab ps, fem_precomp_pool *p) :
      pf(pff), pspt(ps), pool(p) {}
    pre_fem_light_(void) { }   
  };

  fem_precomp_::fem_precomp_(const pre_fem_light_ &ls) :
    pf(ls.pf), pspt(ls.pspt), pool_(ls.pool) {
      for (size_type i = 0; i < pspt->size(); ++i)
	if ((*pspt)[i].size() != pf->dim())
	  DAL_THROW(dimension_error, "dimensions mismatch");
    }

  //  fem_precomp_::fem_precomp_() : pf(0), pspt(0) {}
  
  void fem_precomp_::init_val() const {
    c.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i) 
      pf->base_value((*pspt)[i], c[i]);
  }

  void fem_precomp_::init_grad() const {
    pc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->grad_base_value((*pspt)[i], pc[i]);
  }

  void fem_precomp_::init_hess() const {
    hpc.resize(pspt->size());
    for (size_type i = 0; i < pspt->size(); ++i)
      pf->hess_base_value((*pspt)[i], hpc[i]);
  }

  typedef const fem_precomp_ * pfem_precomp;


  /*  void fem_precomp_not_stored(pfem pf, bgeot::pstored_point_tab pspt,
			      fem_precomp_& fp) {
    fp.assign(pre_fem_light_(pf,pspt));
  }
  */

  /* fem_precomp_pool */
  class fem_precomp_pool_private : 
    public dal::FONC_TABLE<pre_fem_light_, fem_precomp_> {};

  fem_precomp_pool::fem_precomp_pool() : 
    p(new fem_precomp_pool_private()) {}

  void fem_precomp_pool::clear() { 
    delete p; p = new fem_precomp_pool_private(); 
  }
  fem_precomp_pool::~fem_precomp_pool() { delete p; }

  pfem_precomp 
  fem_precomp_pool::operator()(pfem pf, bgeot::pstored_point_tab pspt) {
    return p->add(pre_fem_light_(pf, pspt, this));
  }

  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt)
  { 
    return dal::singleton<fem_precomp_pool>::instance()(pf, pspt);
  }
}  /* end of namespace getfem.                                            */

