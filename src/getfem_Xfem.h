/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_Xfem.h : definition of eXtended fems.                  */
/*           see for instance "A finite element method for crack growth    */
/*           without remeshing", N. Moës, J. Dolbow and T. Belytschko,     */
/*           Int. J. Num. Meth. Engng. 46, 131-150 (1999).                 */
/*                                                                         */
/* Date : April 8, 2003.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
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

#include <getfem_fem.h>
#include <getfem_mesh_fem.h>

/* Works only for tau-equivalent elements                                  */

namespace getfem
{

  // Object representing global functions. To be derived.

  struct virtual_Xfem_func {
    virtual scalar_type operator()(const base_node &) = 0;
  };
  typedef virtual_Xfem_func *pXfem_func;

  struct virtual_Xfem_grad {
    virtual base_vector operator()(const base_node &) = 0;
  };
  typedef virtual_Xfem_grad *pXfem_grad;

  struct virtual_Xfem_hess {
    virtual base_matrix operator()(const base_node &);
  };
  typedef virtual_Xfem_hess *pXfem_hess;
  
  extern pXfem_hess pno_Xfem_hess_defined;
  

  // Xfem definition
  
  class Xfem : public virtual_fem {
  
  protected:
    pfem pfi;
    bool is_valid;
    size_type nb_func;
    std::vector<pXfem_func> funcs; // List of functions to be added
    std::vector<pXfem_grad> grads; // Gradients of theses functions
    std::vector<pXfem_hess> hess;  // Hessians of theses functions
    std::vector<size_type> func_indices;


  public:

    void valid(void);

    virtual size_type nb_dof(void) const;

    void add_func(pXfem_func pXf, pXfem_grad pXg,
		  pXfem_hess pXh = pno_Xfem_hess_defined,
		  size_type ind = size_type(-1));
    
    void interpolation(const base_node &x, const base_matrix &G,
		       bgeot::pgeometric_trans pgt,
		       const base_vector &coeff, base_node &val) const;

    void interpolation(pfem_precomp pfp, size_type ii,
		       const base_matrix &G,
		       bgeot::pgeometric_trans pgt, 
		       const base_vector &coeff, 
		       base_node &val, dim_type Qdim=1) const;

    virtual void interpolation_grad(const base_node &x,
				    const base_matrix &G,
				    bgeot::pgeometric_trans pgt,
				    const base_vector &coeff,
				    base_matrix &val) const;

    void base_value(const base_node &x, base_tensor &t) const;
    void grad_base_value(const base_node &x, base_tensor &t) const;
    void hess_base_value(const base_node &x, base_tensor &t) const;

    void real_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			 size_type ii, const base_matrix &G,
			 base_tensor &t, size_type elt) const;
    
    void real_grad_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			      size_type ii, const base_matrix &G,
			      const base_matrix &B, base_tensor &t,
			      size_type elt) const;
    
    void real_hess_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			      size_type ii, const base_matrix &G,
			      const base_matrix &B3, const base_matrix &B32,
			      base_tensor &t, size_type elt) const;
    
    Xfem(pfem pf);
  };


}  /* end of namespace getfem.                                            */
