/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_precomp.h : pre-computations.                         */
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

#ifndef __GETFEM_PRECOMP_H
#define __GETFEM_PRECOMP_H

#include <bgeot_geometric_trans.h>

namespace getfem
{

  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  struct _pre_geot_light;

  /**
     precomputed geometric transformation operations use this for
     repetitive evaluation of a geometric transformations on a set of
     points "pspt" in the the reference convex which do not change.
  */
  class _geotrans_precomp {
  protected:      
    bgeot::pgeometric_trans pgt;
    bgeot::pstored_point_tab pspt;  /* a set of points in the reference elt*/
    mutable std::vector<base_vector> c;  /* precomputed values for the     */
                                         /* transformation                 */
    mutable std::vector<base_matrix> pc; /* precomputed values for gradient*/
                                         /* of the transformation.         */
    mutable std::vector<base_matrix> hpc; /* precomputed values for hessian*/
                                          /*  of the transformation.       */
    
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
    template <typename CONT> void transform(const CONT& G,
					    bgeot::stored_point_tab& pt_tab);
    template <typename CONT, typename VEC> void transform(const CONT& G, size_type ii, VEC& pt);

    base_node transform(size_type i, const base_matrix &G) const;
    bgeot::pgeometric_trans get_trans() const { return pgt; }
    void assign(const _pre_geot_light &ls);
    _geotrans_precomp(const _pre_geot_light &ls);
    _geotrans_precomp();
  private:
    void init_val() const;
    void init_grad() const;
    void init_hess() const;
  };


  template <typename CONT, typename VEC> 
  void _geotrans_precomp::transform(const CONT& G, size_type j, VEC& pt) {
    size_type k = 0;
    if (c.empty()) init_val();
    for (typename CONT::const_iterator itk = G.begin(); 
         itk != G.end(); ++itk, ++k) {
      for (size_type i=0; i < (*itk).size(); ++i) {
        pt[i] += (*itk)[i] * c[j][k];
      }
    }
  }

  template <typename CONT> 
  void _geotrans_precomp::transform(const CONT& G,
                                    bgeot::stored_point_tab& pt_tab) {
    if (c.empty()) init_val();
    pt_tab.clear(); pt_tab.resize(c.size(), base_node(G[0].size()));
    for (size_type j = 0; j < c.size(); ++j) {
      transform(G, j, pt_tab[j]);
    }
  }
  
  typedef const _geotrans_precomp * pgeotrans_precomp;

  /**
     precomputes a geometric transformation for a fixed set of 
     points in the reference convex. The result is stored in a static
     "cache": do not delete it!
     @param pg the geometric transformation
     @param pspt a pointer to a list of points in the reference convex.CAUTION:
     this array must not be destroyed as long as the fem_precomp is used!!
  */
  pgeotrans_precomp geotrans_precomp(bgeot::pgeometric_trans pg,
				     bgeot::pstored_point_tab pspt);
  
  /**
     precomputes a geometric transformation for a fixed set of 
     points in the reference convex. 
     @param pg the geometric transformation
     @param pspt a pointer to a list of points in the reference convex.CAUTION:
     this array must not be destroyed as long as the fem_precomp is used!!
  */
  void geotrans_precomp_not_stored(bgeot::pgeometric_trans pg,
				   bgeot::pstored_point_tab pspt,
				   _geotrans_precomp& gp);

  /**
     Pre-computations on a fem.     
  */
  class virtual_fem;
  typedef const virtual_fem * pfem;

  struct _pre_fem_light;
  
  class _fem_precomp {
  protected:      
    pfem pf;
    bgeot::pstored_point_tab pspt;
    mutable std::vector<base_tensor> c;
    mutable std::vector<base_tensor> pc;
    mutable std::vector<base_tensor> hpc;

  public:    
    inline const base_tensor &val(size_type i) const
    { if (c.empty()) init_val(); return c[i]; }
    inline const base_tensor &grad(size_type i) const
    { if (pc.empty()) init_grad(); return pc[i]; }
    inline const base_tensor &hess(size_type i) const
    { if (hpc.empty()) init_hess(); return hpc[i]; }
    inline pfem get_pfem() const { return pf; }
    inline bgeot::pstored_point_tab get_point_tab() const { return pspt; }

    void assign(const _pre_fem_light &ls);
    _fem_precomp(const _pre_fem_light &ls);
    _fem_precomp();

  private:
    void init_val() const;
    void init_grad() const;
    void init_hess() const;    
  };
  
  typedef const _fem_precomp * pfem_precomp;

  /**
     statically allocates a fem-precomputation object, and returns a pointer
     to it. The precomputations are "cached", i.e. if this function is called
     two times with the same arguments, a pointer to the same object will be
     returned.
     @param pf a pointer to the fem object.
     @param pspt a pointer to a list of points in the reference convex.CAUTION:
     this array must not be destroyed as long as the fem_precomp is used!!
  */
  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt);

  /**
     fills a _fem_precomp object.
     @param pf a pointer to the fem object.
     @param pspt a pointer to a list of points in the reference convex.CAUTION:
     this array must not be destroyed as long as the fem_precomp is used!!
   */
  void fem_precomp_not_stored(pfem pf, bgeot::pstored_point_tab pspt,
			      _fem_precomp& fp);
}  /* end of namespace getfem.                                            */

#endif
