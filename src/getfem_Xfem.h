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

/*
  what is an Xfem ?

  It is a "base" fem (for example PK(2,2)), with additional base functions.
  These additionnal base functions are the product of:
   - a global function (the virtual_Xfem_func)
   - base functions of another fem (for example PK(2,1))

  The Xfem is built using the add_func member function, which takes as
  parameters, a global function and a fem.

  example of use: enrichment of the finite elements space with
  particular functions, which may represent discontinuities of the
  field in some elements, or singularities in the field (crack tip
  functions..)
*/

#include <getfem_fem.h>
#include <getfem_mesh_fem.h>

/* Works only for tau-equivalent elements                                  */

namespace getfem
{

  // Object representing global functions. To be derived.

  struct virtual_Xfem_func {
    virtual scalar_type val(const base_node &) { DAL_THROW(dal::failure_error,"this Xfem_func has no value"); }
    virtual base_vector grad(const base_node&) { DAL_THROW(dal::failure_error,"this Xfem_func has no gradient"); }
    virtual base_matrix hess(const base_node&) { DAL_THROW(dal::failure_error,"this Xfem_func has no hessian"); }
    virtual ~virtual_Xfem_func() {}
  };
  typedef virtual_Xfem_func *pXfem_func;
  
  // Xfem definition
  
  class Xfem : public virtual_fem {
  
  protected:
    pfem pfb; // base fem
    std::vector<pfem> uniq_pfe; 
    std::vector<size_type> func_pf; // nb_func fems which are enriched (indexes in the array uniq_pfe)
    bool is_valid;
    size_type nb_func;
    std::vector<pXfem_func> funcs; // List of functions to be added
    std::vector<size_type> func_indices;

    void get_fem_precomp_tab(pfem_precomp pfp, std::vector<pfem_precomp>& vpfp) const;
    pfem pfe(size_type k) const { return uniq_pfe[func_pf[k]]; }
  public:

    void valid(void);

    virtual size_type nb_dof(void) const;

    /* ind should be > 0 */
    void add_func(pfem pf, pXfem_func pXf,
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
    
    Xfem(pfem pfb);
  };


}  /* end of namespace getfem.                                            */
