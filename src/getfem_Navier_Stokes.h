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


#ifndef GETFEM_NAVIER_STOKES_H__
#define GETFEM_NAVIER_STOKES_H__

#include <getfem_modeling.h>
#include <getfem_assembling_tensors.h>

namespace getfem {

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


  /* ******************************************************************** */
  /*	Incompressible Navier-Stokes brick.                               */
  /* ******************************************************************** */

# define MDBRICK_NAVIER_STOKES 394329

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_pre_navier_stokes : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    mesh_im &mim;
    mesh_fem &mf_u;
    value_type nu;
    T_MATRIX K;

    virtual void proper_update(void) {
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      gmm::clear(K);
      asm_stiffness_matrix_for_homogeneous_laplacian_componentwise(K, mim, mf_u);
      gmm::scale(K, nu);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, this->nb_dof());
      gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
      asm_navier_stokes_tgm(gmm::sub_matrix(MS.tangent_matrix(), SUBI),
			    mim, mf_u, gmm::sub_vector(MS.state(), SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0, size_type) {
      gmm::sub_interval SUBI(i0, this->nb_dof());
      gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
      asm_navier_stokes_rhs(gmm::sub_vector(MS.residu(), SUBI), mim,
			    mf_u, gmm::sub_vector(MS.state(), SUBI));
    }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), this->nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      this->add_proper_mesh_fem(mf_u, MDBRICK_NAVIER_STOKES);
      this->add_proper_mesh_im(mim);
      this->proper_is_linear_ = false;
      this->proper_is_coercive_ = this->proper_is_symmetric_ = false;
      this->update_from_context();
    }

    mdbrick_pre_navier_stokes(mesh_im &mim_, mesh_fem &mf_u_, value_type nu_)
      : mim(mim_), mf_u(mf_u_), nu(nu_) { init_(); }

  };



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
    virtual void do_compute_residu(MODEL_STATE &, size_type, size_type) {}

    SUBVECTOR get_velocity(MODEL_STATE &MS) 
    { return velocity_part.get_solution(MS); }

    SUBVECTOR get_pressure(MODEL_STATE &MS) 
    { return sub_problem.get_pressure(MS); }

    mdbrick_navier_stokes(mesh_im &mim, mesh_fem &mf_u, mesh_fem &mf_p,
			  value_type nu)
      : velocity_part(mim, mf_u, nu), sub_problem(velocity_part, mf_p) {
      this->add_sub_brick(sub_problem);
      this->update_from_context();
    }
    
  };






 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_NAVIER_STOKES_H__ */
