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

#ifndef BGEOT_PRECOMP_H__
#define BGEOT_PRECOMP_H__

#include <bgeot_geometric_trans.h>

namespace bgeot
{

  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  struct pre_geot_light_;
  class geotrans_precomp_pool;
  /**
     precomputed geometric transformation operations use this for
     repetitive evaluation of a geometric transformations on a set of
     points "pspt" in the the reference convex which do not change.
  */
  class geotrans_precomp_ {
  protected:      
    pgeometric_trans pgt;
    pstored_point_tab pspt;  /* a set of points in the reference elt*/
    mutable std::vector<base_vector> c;  /* precomputed values for the     */
                                         /* transformation                 */
    mutable std::vector<base_matrix> pc; /* precomputed values for gradient*/
                                         /* of the transformation.         */
    mutable std::vector<base_matrix> hpc; /* precomputed values for hessian*/
                                          /*  of the transformation.       */
    geotrans_precomp_pool *pool_;
  public:
    inline const base_vector &val(size_type i) const 
    { if (c.empty()) init_val(); return c[i]; }
    inline const base_matrix &grad(size_type i) const
    { if (pc.empty()) init_grad(); return pc[i]; }
    inline const base_matrix &hessian(size_type i) const
    { if (hpc.empty()) init_hess(); return hpc[i]; }
    
    /**
       Apply the geometric transformation from the reference convex to
       the convex whose vertices are stored in G, to the set of points
       listed in pspt. @param G any container of vertices of the transformed
       convex
       @result pt_tab transformed points
    */
    template <typename CONT> 
    void transform(const CONT& G,
		   stored_point_tab& pt_tab) const;
    template <typename CONT, typename VEC> 
    void transform(const CONT& G, size_type ii, VEC& pt) const;

    base_node transform(size_type i, const base_matrix &G) const;
    pgeometric_trans get_trans() const { return pgt; }
    inline const stored_point_tab& get_point_tab() const { return *pspt; }
    geotrans_precomp_(const pre_geot_light_ &ls);
    geotrans_precomp_pool *pool() const { return pool_; }
  private:
    void init_val() const;
    void init_grad() const;
    void init_hess() const;
  };


  template <typename CONT, typename VEC> 
  void geotrans_precomp_::transform(const CONT& G, size_type j, VEC& pt) const {
    size_type k = 0;
    if (c.empty()) init_val();
    for (typename CONT::const_iterator itk = G.begin(); 
         itk != G.end(); ++itk, ++k) {
      typename CONT::value_type::const_iterator Gk = (*itk).begin();
      typename VEC::iterator ipt = pt.begin();
      for (size_type i=0; i < (*itk).size(); ++i) {
        ipt[i] += Gk[i] * c[j][k];
      }
    }
  }

  template <typename CONT> 
  void geotrans_precomp_::transform(const CONT& G,
                                    stored_point_tab& pt_tab) const {
    if (c.empty()) init_val();
    pt_tab.clear(); pt_tab.resize(c.size(), base_node(G[0].size()));
    for (size_type j = 0; j < c.size(); ++j) {
      transform(G, j, pt_tab[j]);
    }
  }
  
  typedef const geotrans_precomp_ * pgeotrans_precomp;

  /**
     precomputes a geometric transformation for a fixed set of 
     points in the reference convex. The result is stored in a static
     list of geotrans_precomp (global pool)
     @param pg the geometric transformation
     @param pspt a pointer to a list of points in the reference convex.CAUTION:
     this array must not be destroyed as long as the fem_precomp is used!!

     If you need a set of "temporary" geotrans_precomp, create them
     via a geotrans_precomp_pool structure. All memory will be freed
     when this structure will be destroyed.
  */
  pgeotrans_precomp geotrans_precomp(bgeot::pgeometric_trans pg,
				     bgeot::pstored_point_tab pspt);
  
  class geotrans_precomp_pool_private;
  /* handles a set of geotrans_precomp */
  class geotrans_precomp_pool {
    geotrans_precomp_pool_private * p;
  public:
    geotrans_precomp_pool();
    ~geotrans_precomp_pool();
    pgeotrans_precomp operator()(bgeot::pgeometric_trans pg,
				 bgeot::pstored_point_tab pspt);
    void clear();
  };
}
#endif
