/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_mat_elem.h : computation of elementary matrices.      */
/*     									   */
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


#ifndef GETFEM_MAT_ELEM_H__
#define GETFEM_MAT_ELEM_H__

#include <bgeot_geometric_trans.h>
#include <getfem_mat_elem_type.h>
#include <getfem_fem.h>

namespace getfem
{
  /** (optional) callback to be called for each point of the
      integration (i.e. only with approximate integrations). It is
      used by getfem_assembling_tensors to perform reductions before
      integration.
  */
  struct mat_elem_integration_callback {
    /* this vector is filled by mat_elem, it contains a pointer to each tensor
       build from a term in the mat_elem_type (base functions, gradients,
       nonlinear terms etc)
    */
    std::vector<const bgeot::base_tensor*> eltm;
    /** executes the callback
        @param t the destination tensor
        @param first indicates if this is the first integration point (in that case,
               t should be set to the correct size and filled with zeros)
        @param c the current coefficient (contains the norm of Jacobian and the integration weight)
    */
    virtual void exec(bgeot::base_tensor &t, bool first, scalar_type c) = 0;
    virtual ~mat_elem_integration_callback() {}
  };

  /** 
      this class (whose intances are returned by the mat_elem function, see below)
      holds all computations of elementary integrals over convexes or faces of convexes
  */
  class mat_elem_computation {
  protected : 
    
    bgeot::pgeometric_trans pgt;
    pmat_elem_type pme;
    base_matrix pa;
    
  public :
    
    virtual void compute(base_tensor &t, const base_matrix &a,
                         size_type elt, 
			 mat_elem_integration_callback *icb = 0) = 0;
    virtual void compute_on_face(base_tensor &t, const base_matrix &a,
                                 short_type f, size_type elt, 
				 mat_elem_integration_callback *icb = 0) = 0;
    /**
       perform the integration on the volume of a convex.
       @param t     the destination tensor
       @param a     the list of vertices of the convex
       @param elt   the convex number
       @param icb   an optional callback which will be called for each integration point
    */
    template <class CONT> void 
    gen_compute(base_tensor &t, const CONT &a,  size_type elt, 
		mat_elem_integration_callback *icb = 0) { 
      bgeot::vectors_to_base_matrix(pa, a); 
      compute(t, pa, elt, icb); 
    }
    /** 
        perform the integration on a face of the convex
    */
    template <class CONT> void 
    gen_compute_on_face(base_tensor &t,
                        const CONT &a, short_type f, size_type elt, 
			mat_elem_integration_callback *icb = 0) {
      bgeot::vectors_to_base_matrix(pa, a); 
      compute_on_face(t, pa, f, elt, icb); 
    }
    
    virtual ~mat_elem_computation() {}
    virtual size_type memsize() const = 0;
  };

  typedef mat_elem_computation *pmat_elem_computation;

  /** 
      allocate a structure for computation (integration over elements or faces
      of elements) of elementary tensors. Internally this structure is linked to a 
      "cache" which stores some pre-computed data.
  */ 
  pmat_elem_computation mat_elem(pmat_elem_type pm, 
				 pintegration_method pi,
				 bgeot::pgeometric_trans pg,
                                 bool prefer_comp_on_real_element = false);

  /** 
      return the number of bytes used for storage of all cached data associated
      to pmat_elem_computation objects
  */
  size_type stored_mat_elem_memsize();

  /**
     free all caches associated to the given type of elementary matrix (when
     using high degree elements and large mat_elem_type, these may become quite
     memory consuming).
  */
  void mat_elem_forget_mat_elem_type(pmat_elem_type pm);


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MAT_ELEM_H__                                              */
