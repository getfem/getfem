/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_Xfem.C : definition of eXtended fems.                  */
/*           see for instance  "A finite element method for crack growth   */
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
#include <bgeot_geotrans_inv.h>

// Works only for tau-equivalent elements

namespace getfem
{

  // Object representing global functionsFunctional object. To be derived.

  class virtual_Xfem_func {
    virtual long_scalar_type operator(const base_node &) = 0;
  };
  typedef virtual_Xfem_func *pXfem_func;

  class virtual_Xfem_grad {
    virtual base_vector operator(const base_node &) = 0;
  };
  typedef virtual_Xfem_grad *pXfem_grad;

  class virtual_Xfem_hess {
    virtual base_matrix operator(const base_node &) = 0;
  };
  typedef virtual_Xfem_hess *pXfem_hess;
  
  class Xfem : public virtual_fem {
  
  protected:
    size_type nb_func;
    std::vector<pXfem_func> funcs; // List of functions to be added
    std::vector<pXfem_grad> grads; // Gradients of theses functions
    std::vector<pXfem_hess> hess;  // Hessians of theses functions
    std::vector<size_type> func_indices;
    papprox_integration pai;
    pfem pfi;
    bool is_valid;


  public:
    void valid(void) {
      init_cvs_node();
      for (size_type k = 0; k < pfi->nb_base(); ++k)
	add_node(pfi->dof_types()[k], pfi->node_of_dof(k));

      for (size_type k = 0; k < nb_func; ++k) {
	for (size_type j = 0; j < pfi->nb_base(); ++j) {
	  add_node(xfem_dof(pfi->dof_types()[j], func_indices[k]),
		   pfi->node_of_dof(j));
	}
      }
      is_valid = true;
    }

    void add_func(pXfem_func pXf, pXfem_grad pXg, pXfem_hess pXh,
		  size_type ind) {
      nb_func ++;
      funcs.resize(nb_func); grads.resize(nb_func); hess.resize(nb_func);
      dof_concerned.resize(nb_func); func_indices.resize(nb_func);
      funcs[nb_func-1] = pXf;
      grads[nb_func-1] = pXg;
      hess[nb_func-1] = pHh;
      func_indices[nb_func-1] = ind;
      is_valid = false;
    }
    
    // Interpolation method : call the interpolation of pfi for ordinary
    // base function and for each additional function and sum the
    // contributions.
    void interpolation(const base_node &x, const base_matrix &G,
		       bgeot::pgeometric_trans pgt,
		       const base_vector coeff, base_node &val) const {
      base_node val2(val.size());
      base_node xreal = pgt->transform(x, G);
      size_type nbb = pfi->nb_base();
      pfi->interpolation(x, G, pgt, coeff, val);
      std::vector<scalar_type> coeff2(nbb);
      for (size_type i = 0; i < nb_func; ++i) {
	for (size_type j = 0; j != nbb; ++j) {
	  coeff2[j] = coeff[(i+1) * nbb + j] * (*(funcs[i]))(xreal);
	}
	pfi->interpolation(x, G, pgt, coeff2, val2);
	val += val2;
      }
    }

    void interpolation(pfem_precomp pfp, size_type ii,
		       const base_matrix &G,
		       bgeot::pgeometric_trans pgt, 
		       const base_vector &coeff, 
		       base_node &val, dim_type Qdim=1) const {

      // + gérer Qdim ... voir virtual_fem::interpolation
      // gérer le pfp nouveau ...
      base_node val2(val.size());
      base_node xreal = pgt->transform((*(pfp->get_point_tab()))[ii], G);
      size_type nbb = pfi->nb_base();
      pfi->interpolation(pfp..., ii, G, pgt, coeff, val);
      std::vector<scalar_type> coeff2(pfi->nb_base());
      for (size_type i = 0; i < nb_func; ++i) {
	for (size_type j = 0; j != nbb; ++j) {
	  coeff2[j] = coeff[(i+1) * nbb + j] * (*(funcs[i]))(xreal);
	}
	pfi->interpolation(pfp..., ii, G, pgt, coeff2, val2);
	val += val2;
      }
      // + completer par les autres ddl ...
    }

    virtual void interpolation_grad(const base_node &x,
				    const base_matrix &G,
				    bgeot::pgeometric_trans pgt,
				    const base_vector &coeff,
				    base_matrix &val) const {
      // ...
    }

    void base_value(const base_node &x, base_tensor &t) const {
      const bgeot::stored_point_tab *p = &(pai->integration_points());
      for (size_type i = 0; i < p->size(); ++i)
	if (&((*p)[i]) == &x) {
	  bgeot::multi_index mi(2);
	  mi[1] = target_dim(); mi[0] = nb_base();
	  t.adjust_sizes(mi);
	  std::fill(t.begin(), t.end(), 0.0);
	  t[i] = 1.0;
	  return;
	}
      // à changer ... completer par zéro ...
      DAL_THROW(internal_error,
	  "You cannot interpolate this element, use the original element.");
    }
    void grad_base_value(const base_node &x, base_tensor &t) const {
      // il faudrait vérifier dans mat_elem que le fait de passer di comme
      // dimension donne le bon calcul quand des dimensions differentes 
      // interviennent.
      const bgeot::stored_point_tab *p = &(pai->integration_points());
      for (size_type i = 0; i < p->size(); ++i)
	if (&((*p)[i]) == &x) { 
	  bgeot::multi_index mi(3);
	  mi[2] = di; mi[1] = target_dim(); mi[0] = nb_base();
	  t.adjust_sizes(mi);
	  std::fill(t.begin(), t.end(), 0.0);
	  if (with_grad)
	    for (dim_type k = 0; k < di; ++k)
	      t[k * mi[0] + i + pai->nb_points() * (k+1)] = 1.0;
	  return;
	}
       // à changer ...
      DAL_THROW(internal_error,
	  "You cannot interpolate this element, use the original element.");
    }
    void hess_base_value(const base_node &x, base_tensor &t) const {
      const bgeot::stored_point_tab *p = &(pai->integration_points());
      for (size_type i = 0; i < p->size(); ++i)
	if (&((*p)[i]) == &x) {
	  bgeot::multi_index mi(4);
	  mi[2] = di; mi[2] = di; mi[1] = target_dim(); mi[0] = nb_base();
	  t.adjust_sizes(mi);
	  std::fill(t.begin(), t.end(), 0.0);
	  return;
	}
       // à changer ...
      DAL_THROW(internal_error,
	  "You cannot interpolate this element, use the original element.");
    }


    virtual void hess_base_value(const base_node &x, base_tensor &t) const = 0;
    /** Gives the value of all components of the base functions at the point
     *  ii of the pgp possibly using information on pfp and G.
     *  Used by elementary computations.
     *  Matrix M for non tau-equivalent elements not taken into account.
     */
    virtual void real_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
				 size_type ii, const base_matrix &G,
				 base_tensor &t) const {
      // ...
    }
    /** Gives the value of all gradients on the real element of the components
     *  of the base functions at the point ii of the pgp possibly using
     *  information on pfp, G and B. B is the matrix wich transforms the
     *  gradient on the reference element to the gradient on the real element.
     *  Used by elementary computations.
     *  Matrix M for non tau-equivalent elements not taken into account.
     */
    virtual void real_grad_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
				      size_type ii, const base_matrix &G,
				      const base_matrix &B,
				      base_tensor &t) const {
      // ...
    }
    /** Gives the value of all hessians on the real element of the components
     *  of the base functions at the point ii of the pgp possibly using
     *  information on pfp, G, B2 and B32. B2 and B32 are the matrices used
     *  to transform a Hessian on reference element to a Hessian on real
     *  element. Used by elementary computations.
     *  Matrix M for non tau-equivalent elements not taken into account.
     */
    virtual void real_hess_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
				      size_type ii, const base_matrix &G,
				      const base_matrix &B3,
				      const base_matrix &B32,
				      base_tensor &t) const {
      // ...
    }


    
    Xfem(...) {
      is_equiv = is_pol = is_lag = false; es_degree = 5; // a refaire ...
      real_element_defined = true;
      cvr = pai->ref_convex();
      di = ls.pmflf->mf_target().linked_mesh().dim();
      is_valid = false;
    }

  };




}  /* end of namespace getfem.                                            */
