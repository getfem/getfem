/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_precomp.h : pre-computations.                         */
/*     									   */
/*                                                                         */
/* Date : 2004/01/11.                                                      */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Julien Pommier, pommier@gmm.insa-tlse.fr                       */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004 Yves Renard, Julien Pommier.                    */
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

#ifndef GETFEM_PRECOMP_H__
#define GETFEM_PRECOMP_H__

#include <bgeot_geometric_trans.h>

namespace getfem
{
  /**
     Pre-computations on a fem.     
  */
  class virtual_fem;
  typedef const virtual_fem * pfem;

  struct pre_fem_light_; 
  class fem_precomp_pool;
  
  class fem_precomp_ {
  protected:      
    pfem pf;
    bgeot::pstored_point_tab pspt;
    mutable std::vector<base_tensor> c;
    mutable std::vector<base_tensor> pc;
    mutable std::vector<base_tensor> hpc;
    fem_precomp_pool *pool_;
  public:    
    inline const base_tensor &val(size_type i) const
    { if (c.empty()) init_val(); return c[i]; }
    inline const base_tensor &grad(size_type i) const
    { if (pc.empty()) init_grad(); return pc[i]; }
    inline const base_tensor &hess(size_type i) const
    { if (hpc.empty()) init_hess(); return hpc[i]; }
    inline pfem get_pfem() const { return pf; }
    inline const bgeot::stored_point_tab& get_point_tab() const { return *pspt; }

    fem_precomp_(const pre_fem_light_ &ls);
    fem_precomp_pool *pool() const { return pool_; }
  private:
    void init_val() const;
    void init_grad() const;
    void init_hess() const;    
  };
  
  typedef const fem_precomp_ * pfem_precomp;

  /**
     statically allocates a fem-precomputation object, and returns a
     pointer to it. The precomputations are "cached", i.e. they are
     stored in a global pool and if this function is called two times
     with the same arguments, a pointer to the same object will be
     returned.

     @param pf a pointer to the fem object.
     @param pspt a pointer to a list of points in the reference convex.CAUTION:
     this array must not be destroyed as long as the fem_precomp is used!!. 

     Moreover pspt is supposed to identify uniquely the set of
     points. This means that you should NOT alter its content at any
     time after using this function.

     If you need a set of "temporary" fem_precomp, create them
     via a geotrans_precomp_pool structure. All memory will be freed
     when this structure will be destroyed.  */
  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt);

  /**
     handles a pool (i.e. a set) of fem_precomp. The difference with
     the global fem_precomp function is that these fem_precomp objects
     are freed when the fem_precomp_pool is destroyed (they can eat
     much memory). An example of use can be found in the
     interpolation_solution functions of getfem_export.h

     @param pf a pointer to the fem object.
     @param pspt a pointer to a list of points in the reference convex.CAUTION:
     this array must not be destroyed as long as the fem_precomp is used!!

     Moreover pspt is supposed to identify uniquely the set of
     points. This means that you should NOT alter its content until
     the fem_precomp_pool is destroyed.
  */
  class fem_precomp_pool_private;
  class fem_precomp_pool {
    fem_precomp_pool_private * p;
  public:
    fem_precomp_pool();
    ~fem_precomp_pool();
    pfem_precomp operator()(pfem pf, bgeot::pstored_point_tab pspt);
    void clear();
  };
}  /* end of namespace getfem.                                            */

#endif
