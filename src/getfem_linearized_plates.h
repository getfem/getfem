/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_linearized_plates.h                                   */
/*     									   */
/* Date : November 1, 2004.                                                */
/* Authors : Yves Renard, Yves.Renard@insa-toulouse.fr                     */
/*           Michel Salaun, msalaun@ensica.fr                              */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2000-2004  Yves Renard.                                   */
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


#ifndef GETFEM_LINEARIZED_PLATES_H__
#define GETFEM_LINEARIZED_PLATES_H__

iiiiik#include <getfem_modeling.h>
#include <getfem_assembling_tensors.h>

namespace getfem {

  /* ******************************************************************** */
  /*		Linear plate specific assembly procedures.                */
  /* ******************************************************************** */

  template<class MAT, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear
  (const MAT &RM, const mesh_fem &mf_u3, const mesh_fem &mf_theta,
   const mesh_fem &mfdata, const VECT &LAMBDA,const VECT &MU) {
    gmm::sub_interval I1(0, mf_u3.nb_dof());
    gmm::sub_interval I2(mf_u3.nb_dof(), mf_theta.nb_dof());
    
    asm_stiffness_matrix_for_plate_transverse_shear
      (gmm::sub_matrix(RM, I1), gmm::sub_matrix(RM, I1, I2),
       gmm::sub_matrix(RM, I2, I1), gmm::sub_matrix(RM, I2), mf_u3, mf_theta,
       mfdata, LAMBDA, MU);
  }

  template<class MAT, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear
  (const MAT &RM1, const MAT &RM2, const MAT &RM3, const MAT &RM4,
   const mesh_fem &mf_u3, const mesh_fem &mf_theta,
   const mesh_fem &mfdata, const VECT &LAMBDA,const VECT &MU) {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    
    if (mf_u3.get_qdim() != 1 || mf_theta.get_qdim() != 2)
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    /* e = strain tensor,
       M = 2*mu*e(u):e(v) + lambda*tr(e(u))*tr(e(v))
    */
    generic_assembly assem("lambda=data$1(#3); mu=data$2(#3);"
			   "t1=comp(Grad(#1).Grad(#1).Base(#3));"
			   "t2=comp(vBase(#2).Grad(#1).Base(#3));"
			   "t3=comp(vBase(#2).vBase(#2).Base(#3));"
                           "M$1(#1,#1)+= sym(t(:,i,:,i,j).mu(j)"
			   "+t1(:,j,i,:,i,j,k).mu(k)" // incoherent
                             " + t(:,i,i,:,j,j,k).lambda(k))"
			   "M$2(#2,#1)+=t2(:,i,:,i,j).mu(j)" // incoherent
			   "M$3(#2,#2)+=t3(:,i,:,i,j).mu(j)" // incoherent
			   );
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mf(mfdata);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_mat(const_cast<MAT &>(RM1));
    assem.push_mat(const_cast<MAT &>(RM2));
    assem.push_mat(const_cast<MAT &>(RM4));
    assem.volumic_assembly();
    gmm::copy(gmm::transposed(RM2), RM3);
  }


  /* ******************************************************************** */
  /*		Linear plate model brick.                                 */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_isotropic_linearized_plate
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    mesh_fem &mf_ut;
    mesh_fem &mf_u3;
    mesh_fem &mf_theta;
    mesh_fem &mf_data;
    VECTOR lambda_, mu_;
    bool homogeneous;
    bool matrix_stored;
    T_MATRIX K;

    void compute_K(void) {
      gmm::clear(K);
      gmm::resize(K, nb_dof(), nb_dof());
      VECTOR lambda(mf_data.nb_dof()), mu(mf_data.nb_dof());
      if (homogeneous) {
	std::fill(lambda.begin(), lambda.end(), value_type(lambda_[0]));
	std::fill(mu.begin(), mu.end(), value_type(mu_[0]));
      }
      else { gmm::copy(lambda_, lambda); gmm::copy(mu_, mu); }
      gmm::sub_interval I1(0, mf_ut.nb_dof());
      gmm::sub_interval I2(mf_ut.nb_dof(), mf_u3.nb_dof()+mf_theta.nb_dof());
      gmm::sub_interval I3(mf_ut.nb_dof() + mf_u3.nb_dof(), mf_theta.nb_dof());
      // Il faut multiplier les coeff par epsilon ou epsilon^2 suivant les cas
      asm_stiffness_matrix_for_linear_elasticity
	(gmm::sub_matrix(K, I1), mf_ut, mf_data, lambda, mu);
      asm_stiffness_matrix_for_linear_elasticity
	(gmm::sub_matrix(K, I3), mf_theta, mf_data, lambda, mu);
      asm_stiffness_matrix_for_plate_transverse_shear
	(gmm::sub_matrix(K, I2), mf_u3, mf_theta, mf_data, lambda3, mu3);
      this->computed();
    }
    
  public :

    virtual bool is_linear(void) { return true; }
    virtual bool is_coercive(void) { return true; }
    virtual bool is_symmetric(void) { return true; }
    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_dof(void)
    { return mf_ut.nb_dof() + mf_u3.nb_dof() + mf_theta.nb_dof(); }
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false) {
      if (modified && !matrix_stored) 
	DAL_THROW(failure_error, "The residu will not be consistant. "
		  "Use this brick with the stiffness matrix stored option");
      react(MS, i0, modified);
      gmm::sub_interval SUBI(i0, nb_dof());
      if (this->to_be_computed()
	  || (!matrix_stored && this->to_be_transferred()))
	compute_K();
      if (this->to_be_transferred()) { 
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
	this->transferred();
      }
      if (!matrix_stored) gmm::clear(K);
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0) {
      react(MS, i0, false);
      gmm::sub_interval SUBI(i0, nb_dof());
      if (this->to_be_computed()) { 
	compute_K();
	if (!matrix_stored) {
	  gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	  gmm::clear(K);
	}
      }
      if (matrix_stored) {
	gmm::mult(K, gmm::sub_vector(MS.state(), SUBI),
		  gmm::sub_vector(MS.residu(), SUBI));
      } else {
	gmm::mult(gmm::sub_matrix(MS.tangent_matrix(), SUBI),
		  gmm::sub_vector(MS.state(), SUBI),
		  gmm::sub_vector(MS.residu(), SUBI));
      }
    }

    virtual mesh_fem &main_mesh_fem(void) { return mf_ut; }

    void set_Lame_coeff(value_type lambdai, value_type mui) {
      homogeneous = true;
      gmm::resize(lambda_, 1); lambda_[0] = lambdai;
      gmm::resize(mu_, 1); mu_[0] = mui;
      this->force_recompute();
    }

    void set_Lame_coeff(const VECTOR &lambdai, const VECTOR &mui) {
      homogeneous = false;
      gmm::resize(lambda_, mf_data.nb_dof()); gmm::copy(lambdai, lambda_);
      gmm::resize(mu_, mf_data.nb_dof()); gmm::copy(mui, mu_);
      this->force_recompute();
    }

    template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    // constructor for a homogeneous material (constant lambda and mu)
    mdbrick_isotropic_linearized_plate
    (mesh_fem &mf_ut_, mesh_fem &mf_u3_, mesh_fem &mf_theta_,
     mesh_fem &mf_data_, value_type lambdai, value_type mui,
     bool mat_stored = false)
      : mf_ut(mf_ut_), mf_u3(mf_u3_), mf_theta(mf_theta_), mf_data(mf_data_),
	matrix_stored(mat_stored) {
      set_Lame_coeff(lambdai, mui);
      this->add_dependency(mf_ut); this->add_dependency(mf_u3);
      this->add_dependency(mf_theta); this->add_dependency(mf_data);
    }

    // constructor for a non-homogeneous material
    mdbrick_isotropic_linearized_plate
    (mesh_fem &mf_ut_, mesh_fem &mf_u3_, mesh_fem &mf_theta_,
     mesh_fem &mf_data_, const VECTOR &lambdai, const VECTOR &mui,
     bool mat_stored = false)
      : mf_ut(mf_ut_), mf_u3(mf_u3_), mf_theta(mf_theta_), mf_data(mf_data_),
	matrix_stored(mat_stored) {
      set_Lame_coeff(lambdai, mui);
      this->add_dependency(mf_ut); this->add_dependency(mf_u3);
      this->add_dependency(mf_theta); this->add_dependency(mf_data);
    }
 
  };

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LINEARIZED_PLATES_H__ */
