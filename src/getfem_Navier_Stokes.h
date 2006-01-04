// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_Navier_Stokes.h : 
//           
// Date    : April 15, 2005.
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2000-2005 Yves Renard
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

/**@file getfem_Navier_Stokes.h
   @brief Navier-Stokes dynamic brick.
*/
#ifndef GETFEM_NAVIER_STOKES_H__
#define GETFEM_NAVIER_STOKES_H__

#include <getfem_modeling.h>
#include <getfem_assembling_tensors.h>

namespace getfem {

  /**
     assembly of Tangent matrix for Navier-Stokes.
     @ingroup asm
  */
  template<typename MAT, typename VECT>
  void asm_navier_stokes_tgm(const MAT &M, 
			     const mesh_im &mim, 
			     const mesh_fem &mf,
			     const VECT &U,
			     const mesh_region &rg = mesh_region::all_convexes()) {
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    getfem::generic_assembly assem;
    assem.set("u=data(#1);"
	      "t=comp(vGrad(#1).vBase(#1).vBase(#1));"
	      "M(#1, #1) += u(i).t(i,k,j,:,k,:,j);"
	      "M(#1, #1) += u(i).t(:,j,k,:,k,i,j);"
	      "M(#1, #1) += u(i).t(i,j,j,:,k,:,k)/2;"
	      "M(#1, #1) += u(i).t(:,k,k,:,j,i,j)/2;");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mat(const_cast<MAT&>(M));
    assem.push_data(U);
    assem.assembly(rg);
  }

  /**
     assembly of right hand side for Navier-Stokes.
     @ingroup asm
  */
  template<typename VECT1, typename VECT2> 
  void asm_navier_stokes_rhs(const VECT1 &V, 
			     const mesh_im &mim, 
			     const mesh_fem &mf,
			     const VECT2 &U,
			     const mesh_region &rg = mesh_region::all_convexes()) {
    if (mf.get_qdim() != mf.linked_mesh().dim())
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");

    getfem::generic_assembly assem;
    assem.set("u=data(#1);"
	      "t=comp(vBase(#1).vGrad(#1).vBase(#1));"
	      "V(#1) += u(i).u(j).t(i,k,j,k,l,:,l);"
	      "V(#1) += u(i).u(j).t(i,k,j,l,l,:,k)/2;");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_vec(const_cast<VECT1&>(V));
    assem.push_data(U);
    assem.assembly(rg);
  }


# define MDBRICK_NAVIER_STOKES 394329

  /** Internal brick for mdbrick_pre_navier_stokes */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_pre_navier_stokes : public mdbrick_abstract_linear_pde<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    value_type nu;

    virtual void proper_update_K(void) {
      DAL_TRACE2("Assembling laplacian for mdbrick_pre_navier_stokes");
      asm_stiffness_matrix_for_homogeneous_laplacian_componentwise(this->K, this->mim, this->mf_u);
      gmm::scale(this->K, nu);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, this->nb_dof());
      gmm::copy(this->get_K(), gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      asm_navier_stokes_tgm(gmm::sub_matrix(MS.tangent_matrix(), SUBI),
			    this->mim, this->mf_u, gmm::sub_vector(MS.state(), SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, this->nb_dof());
      gmm::mult(this->get_K(), gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residual(), SUBI));
      asm_navier_stokes_rhs(gmm::sub_vector(MS.residual(), SUBI), this->mim,
			    this->mf_u, gmm::sub_vector(MS.state(), SUBI));
    }

    mdbrick_pre_navier_stokes(const mesh_im &mim_, const mesh_fem &mf_u_,
			      value_type nu_)
      : mdbrick_abstract_linear_pde<MODEL_STATE>(mim_, mf_u_, MDBRICK_NAVIER_STOKES),
	nu(nu_) {
      this->proper_is_linear_ = false;
      this->proper_is_coercive_ = this->proper_is_symmetric_ = false;
      this->force_update();
    }

  };



  /**
     Incompressible Navier-Stokes brick.
     @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_navier_stokes : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::sub_vector_type<VECTOR *,
				 gmm::sub_interval>::vector_type SUBVECTOR;

    mdbrick_pre_navier_stokes<MODEL_STATE> velocity_part;
    mdbrick_linear_incomp<MODEL_STATE> sub_problem;

    virtual void proper_update(void) {}

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					size_type) {}
    virtual void do_compute_residual(MODEL_STATE &, size_type, size_type) {}

    SUBVECTOR get_velocity(MODEL_STATE &MS) 
    { return velocity_part.get_solution(MS); }

    SUBVECTOR get_pressure(MODEL_STATE &MS) 
    { return sub_problem.get_pressure(MS); }

    mdbrick_navier_stokes(const mesh_im &mim, const mesh_fem &mf_u,
			  const mesh_fem &mf_p, value_type nu)
      : velocity_part(mim, mf_u, nu), sub_problem(velocity_part, mf_p) {
      this->add_sub_brick(sub_problem);
      this->force_update();
    }
    
  };






 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NAVIER_STOKES_H__ */
