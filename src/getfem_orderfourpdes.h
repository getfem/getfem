// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_orderfourspdes.h : assembly procedures and bricks 
//                                     for order four pdes.
// Date    : January 6, 2006.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard
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

/**@file getfem_orderfourpdes.h
   @brief assembly procedures and bricks for order four pdes.
*/
#ifndef GETFEM_ORDERFOURPDES_H_
#define GETFEM_ORDERFOURPDES_H__

#include <getfem_modeling.h>
#include <getfem_assembling_tensors.h>

namespace getfem {

  /**
     assembly of @f$\int_\Gamma{\partial_n u f}@f$.
     @ingroup asm
  */
  template<typename VECT1, typename VECT2>
  void asm_normal_derivative_source_term(VECT1 &B, const mesh_im &mim,
					 const mesh_fem &mf,
					 const mesh_fem &mf_data,
					 const VECT2 &F,
					 const mesh_region &rg) {
    if (mf_data.get_qdim() != 1)
      DAL_THROW(invalid_argument,
		"invalid data mesh fem (Qdim=1 required)");

    std::string s;
    if (mf.get_qdim() == 1)
      s = "F=data(#2);"
	"V(#1)+=comp(Grad(#1).Normal().Base(#2))(:,i,i,j).F(j);";
    else
      s = "F=data(qdim(#1),#2);"
	"V(#1)+=comp(vGrad(#1).Normal.Base(#2))(:,i,k,k,j).F(i,j);";
    asm_real_or_complex_1_param(B, mim, mf, mf_data, F, rg, s);
  }
  
  /**
     assembly of @f$\int_\Omega \Delta u.\Delta v@f$.
     
     @ingroup asm
  */
  template<typename MAT, typename VECT>
  void asm_stiffness_matrix_for_bilaplacian
  (const MAT &M, const mesh_im &mim, const mesh_fem &mf,
   const mesh_fem &mf_data, const VECT &A,
   const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem
      ("a=data$1(#2);"
       "M(#1,#1)+=sym(comp(Hess(#1).Hess(#1).Base(#2))(:,i,i,:,j,j,k).a(k))");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_data);
    assem.push_data(A);
    assem.push_mat(const_cast<MAT &>(M));
    assem.assembly(rg);
  }

  /* ******************************************************************** */
  /*		Bilaplacian brick.                                        */
  /* ******************************************************************** */

# define MDBRICK_BILAPLACIAN 783465
  
  /** Bilaplacian brick @f$ \Delta \Delta u @f$.

  @see asm_stiffness_matrix_for_bilaplacian
  @ingroup bricks 
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_bilaplacian
    : public mdbrick_abstract_linear_pde<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    /* coeff_ a scalar field */
    mdbrick_parameter<VECTOR> coeff_;

    void proper_update_K(void) {
      asm_stiffness_matrix_for_bilaplacian
	(this->K, this->mim, this->mf_u, coeff().mf(),  coeff().get(),
	 this->mf_u.linked_mesh().get_mpi_region());
    }
  public :

    /** accessor to the coefficient k */
    mdbrick_parameter<VECTOR> &coeff() { return coeff_; }
    const mdbrick_parameter<VECTOR> &coeff() const { return  coeff_; }

    /** Switch between a scalar coefficient, a NxN matrix field (with
	N = mf_u.linked_mesh().dim()), and a NxNxNxN tensor field. */
    void set_coeff_dimension(unsigned d) { coeff_.redim(d); }

    /** Constructor, the default coeff is a scalar equal to one
	(i.e. it gives the Laplace operator).

        The coeff can be later changed.

	@param mim the integration method that is used. 
	@param mf_u the mesh_fem for the unknown u.
	@param k the scalar value of the coefficient.
    */
    mdbrick_bilaplacian(const mesh_im &mim_,
			const mesh_fem &mf_u_, value_type k = 1.)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_,
						 MDBRICK_BILAPLACIAN),
	coeff_("coeff", mf_u_.linked_mesh(), this) {
      coeff_.set(k);
    }
  };



}  /* end of namespace getfem.                                             */


#endif /* GETFEM_ORDERFOURPDES_H__  */
