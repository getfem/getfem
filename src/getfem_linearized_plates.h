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

#include <getfem_modeling.h>

namespace getfem {

  /* ******************************************************************** */
  /*		Linear plate specific assembly procedures.                */
  /* ******************************************************************** */

  template<class MAT, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear
  (const MAT &RM, const mesh_fem &mf_u3, const mesh_fem &mf_theta,
   const mesh_fem &mfdata, const VECT &MU) {
    gmm::sub_interval I1(0, mf_u3.nb_dof());
    gmm::sub_interval I2(mf_u3.nb_dof(), mf_theta.nb_dof());
    
    asm_stiffness_matrix_for_plate_transverse_shear
      (gmm::sub_matrix(RM, I1), gmm::sub_matrix(RM, I1, I2),
       gmm::transposed(gmm::sub_matrix(RM, I2, I1)),
       gmm::sub_matrix(RM, I2), mf_u3, mf_theta, mfdata, MU);
  }

  template<class MAT, class MAT3, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear
  (const MAT &RM1, const MAT &RM2, const MAT3 &RM3, const MAT &RM4,
   const mesh_fem &mf_u3, const mesh_fem &mf_theta,
   const mesh_fem &mfdata, const VECT &MU) {
    if (mfdata.get_qdim() != 1)
      DAL_THROW(invalid_argument, "invalid data mesh fem (Qdim=1 required)");
    
    if (mf_u3.get_qdim() != 1 || mf_theta.get_qdim() != 2)
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    generic_assembly assem("mu=data$1(#3);"
			   "t1=comp(Grad(#1).Grad(#1).Base(#3));"
			   "M$1(#1,#1)+=sym(t1(:,i,:,i,j).mu(j));"
			   // "t2=comp(vBase(#2).vBase(#2).Base(#3));"
			   "t2=comp(#1,vBase(#2).vBase(#2).Base(#3));"
			   "M$4(#2,#2)+=sym(t2(:,i,:,i,j).mu(j));"
			   // "M$4(#2,#2)+=sym(t2(k,:,i,:,i,j).mu(j));"
			   "t3=comp(Grad(#1).vBase(#2).Base(#3));"
			   "M$2(#1,#2)+=t3(:,i,:,i,j).mu(j);"
			   "M$3(#1,#2)+=t3(:,i,:,i,j).mu(j);"
			   );
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mf(mfdata);
    assem.push_data(MU);
    assem.push_mat(const_cast<MAT &>(RM1));
    assem.push_mat(const_cast<MAT &>(RM2));
    assem.push_mat(const_cast<MAT3 &>(RM3));
    assem.push_mat(const_cast<MAT &>(RM4));
    assem.volumic_assembly();
  }


  /* ******************************************************************** */
  /*		Linear plate model brick.                                 */
  /* ******************************************************************** */

# define MDBRICK_LINEAR_PLATE 897523

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
    value_type epsilon;
    VECTOR lambda_, mu_;
    bool homogeneous;
    bool matrix_stored;
    T_MATRIX K;

    void compute_K(void) {
      gmm::clear(K);
      gmm::resize(K, this->nb_dof(), this->nb_dof());
      VECTOR lambda(mf_data.nb_dof()), mu(mf_data.nb_dof());
      if (homogeneous) {
	std::fill(lambda.begin(), lambda.end(), value_type(lambda_[0]));
	std::fill(mu.begin(), mu.end(), value_type(mu_[0]));
      }
      else { gmm::copy(lambda_, lambda); gmm::copy(mu_, mu); }
      gmm::sub_interval I1(0, mf_ut.nb_dof());
      gmm::sub_interval I2(mf_ut.nb_dof(), mf_u3.nb_dof()+mf_theta.nb_dof());
      gmm::sub_interval I3(mf_ut.nb_dof() + mf_u3.nb_dof(), mf_theta.nb_dof());
      gmm::scale(lambda, value_type(2) * epsilon);
      gmm::scale(mu, value_type(2) * epsilon);
      asm_stiffness_matrix_for_linear_elasticity
	(gmm::sub_matrix(K, I1), mf_ut, mf_data, lambda, mu);
      gmm::scale(mu, value_type(1) / value_type(2));
      asm_stiffness_matrix_for_plate_transverse_shear
	(gmm::sub_matrix(K, I2), mf_u3, mf_theta, mf_data, mu);
      gmm::scale(lambda, epsilon * epsilon / value_type(3));
      gmm::scale(mu, value_type(2) * epsilon * epsilon / value_type(3));
      asm_stiffness_matrix_for_linear_elasticity
	(gmm::sub_matrix(K, I3), mf_theta, mf_data, lambda, mu);
      this->computed();
    }
    
  public :

    virtual void mixed_variables(dal::bit_vector &, size_type = 0) {}
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type = 0, bool modified = false) {
      if (modified && !matrix_stored) 
	DAL_THROW(failure_error, "The residu will not be consistant. "
		  "Use this brick with the stiffness matrix stored option");
      react(MS, i0, modified);
      gmm::sub_interval SUBI(i0, this->nb_dof());
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
      gmm::sub_interval SUBI(i0, this->nb_dof());
      if (this->to_be_computed()) { 
	compute_K();
	if (!matrix_stored) {
	  gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI)); 
	  gmm::clear(K);
	  this->transferred();
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

    void set_elastic_coeff(value_type E, value_type nu) {
      set_Lame_coeff(nu * E / (value_type(1) - nu*nu),
		     E/(value_type(2)*(value_type(1)+nu)));
    }

    template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), this->nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    template<typename VECT> void get_ut(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), mf_ut.nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }
    template<typename VECT> void get_u3(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index() + mf_ut.nb_dof(),
			     mf_u3.nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }
    template<typename VECT> void get_theta(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index() + mf_ut.nb_dof()
			     + mf_u3.nb_dof(), mf_theta.nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    void init_(void) {
      if (mf_ut.get_qdim() != 2)
	DAL_THROW(failure_error, "Qdim of mf_ut should be 2.");
      if (mf_u3.get_qdim() != 1)
	DAL_THROW(failure_error, "Qdim of mf_ut should be 1.");
      if (mf_theta.get_qdim() != 2)
	DAL_THROW(failure_error, "Qdim of mf_theta should be 2.");
      this->add_proper_mesh_fem(mf_ut, MDBRICK_LINEAR_PLATE, 1);
      this->add_proper_mesh_fem(mf_u3, MDBRICK_LINEAR_PLATE, 0);
      this->add_proper_mesh_fem(mf_theta, MDBRICK_LINEAR_PLATE, 0);
      this->add_dependency(mf_data);
      this->update_from_context();
    }

    // constructor for a homogeneous material (constant lambda and mu)
    mdbrick_isotropic_linearized_plate
    (mesh_fem &mf_ut_, mesh_fem &mf_u3_, mesh_fem &mf_theta_,
     mesh_fem &mf_data_, value_type lambdai, value_type mui,
     double epsilon_, bool mat_stored = false)
      : mf_ut(mf_ut_), mf_u3(mf_u3_), mf_theta(mf_theta_), mf_data(mf_data_),
	epsilon(epsilon_), matrix_stored(mat_stored)
    { set_Lame_coeff(lambdai, mui); init_(); }

    // constructor for a non-homogeneous material
    mdbrick_isotropic_linearized_plate
    (mesh_fem &mf_ut_, mesh_fem &mf_u3_, mesh_fem &mf_theta_,
     mesh_fem &mf_data_, const VECTOR &lambdai, const VECTOR &mui,
     double epsilon_, bool mat_stored = false)
      : mf_ut(mf_ut_), mf_u3(mf_u3_), mf_theta(mf_theta_), mf_data(mf_data_),
	epsilon(epsilon_), matrix_stored(mat_stored)
    { set_Lame_coeff(lambdai, mui); init_(); }
 
  };


  /* ******************************************************************** */
  /*		Mixed linear plate specific assembly procedures.          */
  /* ******************************************************************** */

  template<class MAT>
  void asm_coupling_u3theta(const MAT &RM, const mesh_fem &mf_u3,
			    const mesh_fem &mf_theta) {
    
    if (mf_u3.get_qdim() != 1 || mf_theta.get_qdim() != 2)
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    generic_assembly assem("t1=comp(Grad(#1).vBase(#2));"
			   "M$1(#1,#2)+=t1(:,i,:,i);");
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mat(const_cast<MAT &>(RM));
    assem.volumic_assembly();
  }

  template<class MAT>
  void asm_coupling_psitheta(const MAT &RM, const mesh_fem &mf_u3,
			     const mesh_fem &mf_theta) {
    
    if (mf_u3.get_qdim() != 1 || mf_theta.get_qdim() != 2)
      DAL_THROW(std::logic_error, "wrong qdim for the mesh_fem");
    generic_assembly assem("t1=comp(Base(#1).vGrad(#2));"
			   "M$1(#1,#2)+=t1(:,:,2,1)-t1(:,:,1,2);");
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mat(const_cast<MAT &>(RM));
    assem.volumic_assembly();
  }

  //   template<class MAT> void affiche_moi_valp(const MAT &M) {
  //     size_type nrows = gmm::mat_nrows(M);
  //     size_type ncols = gmm::mat_ncols(M);
  //     size_type nreddof = std::min(nrows, ncols);
  //     gmm::dense_matrix<double> MM(nreddof, nreddof);
  //     if (nrows < ncols)
  //       gmm::mult(M, gmm::transposed(M), MM);
  //     else if (nrows > ncols)
  //       gmm::mult(gmm::transposed(M), M, MM);
  //     else 
  //       gmm::copy(M, MM);
  //     std::vector<double> eigval(nreddof);
  //     gmm::symmetric_qr_algorithm(MM, eigval);
  //     //     std::sort(eigval.begin(), eigval.end(), Esort);
  //     cout << "eival = " << eigval << endl;
  //   }


  /* ******************************************************************** */
  /*		Mixed linear plate model brick.                           */
  /* ******************************************************************** */

# define MDBRICK_MIXED_LINEAR_PLATE 213456

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_mixed_isotropic_linearized_plate
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    mesh_fem &mf_ut;
    mesh_fem &mf_u3;
    mesh_fem &mf_theta;
    mesh_fem &mf_data;
    value_type epsilon;
    VECTOR lambda_, mu_;
    bool homogeneous;
    bool matrix_stored, symmetrized;
    T_MATRIX K;

    void compute_K(void) {
      gmm::clear(K);
      gmm::resize(K, this->nb_dof(), this->nb_dof());
      VECTOR lambda(mf_data.nb_dof()), mu(mf_data.nb_dof());
      if (homogeneous) {
	std::fill(lambda.begin(), lambda.end(), value_type(lambda_[0]));
	std::fill(mu.begin(), mu.end(), value_type(mu_[0]));
      }
      else { gmm::copy(lambda_, lambda); gmm::copy(mu_, mu); }
      size_type nd1=mf_ut.nb_dof(), nd2=mf_u3.nb_dof(), nd3=mf_theta.nb_dof();
      gmm::sub_interval I1(0, nd1), I2(nd1, nd2), I3(nd1 + nd2, nd3);
      gmm::sub_interval I4(nd1 + nd2 + nd3, nd2), I5(nd1 + 2*nd2 + nd3, nd2);
      
      asm_stiffness_matrix_for_linear_elasticity
	(gmm::sub_matrix(K, I1), mf_ut, mf_data, lambda, mu);
      gmm::scale(gmm::sub_matrix(K, I1), value_type(2) * epsilon);
      
      
      asm_stiffness_matrix_for_homogeneous_laplacian(gmm::sub_matrix(K, I2),
						     mf_u3);
      gmm::scale(gmm::sub_matrix(K, I2),
		 value_type(2) * epsilon * epsilon * epsilon / value_type(3));
      
      asm_stiffness_matrix_for_linear_elasticity
	(gmm::sub_matrix(K, I3), mf_theta, mf_data, lambda, mu);
      //       gmm::scale(gmm::sub_matrix(K, I3),
      //  		 value_type(2) * epsilon * epsilon * epsilon / value_type(3));
      
      
      asm_coupling_u3theta(gmm::sub_matrix(K, I2, I3), mf_u3, mf_theta);
      gmm::scale(gmm::sub_matrix(K, I2, I3),
		 value_type(2) * epsilon * epsilon * epsilon / value_type(3));

      //       cout << "\n\nval p de I2 I3 ";
      //       affiche_moi_valp(gmm::sub_matrix(K, I2, I3));


      asm_coupling_psitheta(gmm::sub_matrix(K, I4, I3), mf_u3, mf_theta);
      gmm::scale(gmm::sub_matrix(K, I4, I3), epsilon*epsilon/value_type(3));

      //       cout << "\n\nval p de I4 I3 ";
      //       affiche_moi_valp(gmm::sub_matrix(K, I4, I3));
     

      asm_coupling_psitheta(gmm::transposed(gmm::sub_matrix(K, I3, I4)),
			    mf_u3, mf_theta);
      gmm::scale(gmm::sub_matrix(K, I3, I4), epsilon*epsilon/value_type(3));


      asm_coupling_u3theta(gmm::transposed(gmm::sub_matrix(K, I3, I5)),
			   mf_u3, mf_theta);
      gmm::scale(gmm::sub_matrix(K, I3, I5), epsilon*epsilon/value_type(3));

      if (!symmetrized)
	asm_stiffness_matrix_for_homogeneous_laplacian(gmm::sub_matrix(K, I5),
						       mf_u3);
      if (symmetrized) {
	asm_mass_matrix(gmm::sub_matrix(K, I3), mf_theta);
	asm_coupling_u3theta(gmm::transposed(gmm::sub_matrix(K, I3, I2)),
			     mf_u3, mf_theta);
	gmm::scale(gmm::sub_matrix(K, I3, I2),
		   value_type(2) * epsilon * epsilon * epsilon / value_type(3));

	asm_stiffness_matrix_for_homogeneous_laplacian
	  (gmm::sub_matrix(K, I2, I5), mf_u3);
	gmm::scale(gmm::sub_matrix(K, I2, I5), epsilon*epsilon/value_type(3));
	asm_stiffness_matrix_for_homogeneous_laplacian
	  (gmm::sub_matrix(K, I5, I2), mf_u3);
	gmm::scale(gmm::sub_matrix(K, I5, I2), epsilon*epsilon/value_type(3));
	asm_coupling_u3theta(gmm::sub_matrix(K, I5, I3), mf_u3, mf_theta);
	gmm::scale(gmm::sub_matrix(K, I5, I3), epsilon*epsilon/value_type(3));
      }
      gmm::scale(gmm::sub_matrix(K, I3),
 		 value_type(2) * epsilon * epsilon * epsilon / value_type(3));

      this->computed();
    }
    
  public :
    virtual void mixed_variables(dal::bit_vector &bv, size_type i0 = 0) {
      bv.add(i0+ mf_ut.nb_dof() + mf_u3.nb_dof()+mf_theta.nb_dof(),
	     mf_u3.nb_dof()*2);
    }
    //    virtual size_type nb_constraints(void) { return 1; } 
    virtual size_type nb_constraints(void) { return 0; }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0= 0, bool modified = false) {
      if (modified && !matrix_stored) 
	DAL_THROW(failure_error, "The residu will not be consistant. "
		  "Use this brick with the stiffness matrix stored option");
      react(MS, i0, modified);
      gmm::sub_interval SUBI(i0, this->nb_dof());
      if (this->to_be_computed()
	  || (!matrix_stored && this->to_be_transferred()))
	compute_K();
      if (this->to_be_transferred()) { 
	gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
	this->transferred();
      }
      if (!matrix_stored) gmm::clear(K);
      gmm::clear(gmm::sub_matrix(MS.constraints_matrix(),
				 gmm::sub_interval(j0, 1),
				 gmm::sub_interval(i0,  mf_ut.nb_dof()
						   + 3 * mf_u3.nb_dof() + mf_theta.nb_dof())));
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0,
				size_type = 0) {
      react(MS, i0, false);
      gmm::sub_interval SUBI(i0, this->nb_dof());
      if (this->to_be_computed()) { 
	compute_K();
	if (!matrix_stored) {
	  gmm::copy(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
	  this->transferred();
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

    void set_elastic_coeff(value_type E, value_type nu) {
      set_Lame_coeff(nu * E / (value_type(1) - nu*nu),
		     E/(value_type(2)*(value_type(1)+nu)));
    }

    template<typename VECT> void get_solution(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), this->nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }
    template<typename VECT> void get_ut(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index(), mf_ut.nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }
    template<typename VECT> void get_u3(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index() + mf_ut.nb_dof(),
			     mf_u3.nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }
    template<typename VECT> void get_theta(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(this->first_index() + mf_ut.nb_dof()
			     + mf_u3.nb_dof(), mf_theta.nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    void init_(void) {
      size_type info = 1 + (symmetrized ? 2 : 0);
      if (mf_ut.get_qdim() != 2)
	DAL_THROW(failure_error, "Qdim of mf_ut should be 2.");
      if (mf_u3.get_qdim() != 1)
	DAL_THROW(failure_error, "Qdim of mf_ut should be 1.");
      if (mf_theta.get_qdim() != 2)
	DAL_THROW(failure_error, "Qdim of mf_theta should be 2.");
      this->add_proper_mesh_fem(mf_ut, MDBRICK_MIXED_LINEAR_PLATE, info );
      this->add_proper_mesh_fem(mf_u3, MDBRICK_MIXED_LINEAR_PLATE, 0);
      this->add_proper_mesh_fem(mf_theta, MDBRICK_MIXED_LINEAR_PLATE, 0); 
      this->add_proper_mesh_fem(mf_u3, MDBRICK_MIXED_LINEAR_PLATE, 0);
      this->add_proper_mesh_fem(mf_u3, MDBRICK_MIXED_LINEAR_PLATE, 0);
      this->add_dependency(mf_data);
      this->update_from_context();
      this->proper_is_symmetric_ = symmetrized;
      this->proper_is_coercive_ = false;
    }

    // constructor for a homogeneous material (constant lambda and mu)
    mdbrick_mixed_isotropic_linearized_plate
    (mesh_fem &mf_ut_, mesh_fem &mf_u3_, mesh_fem &mf_theta_,
     mesh_fem &mf_data_, value_type lambdai, value_type mui,
     double epsilon_, bool sym = false, bool mat_stored = false)
      : mf_ut(mf_ut_), mf_u3(mf_u3_), mf_theta(mf_theta_), mf_data(mf_data_),
	epsilon(epsilon_), matrix_stored(mat_stored), symmetrized(sym)
    { set_Lame_coeff(lambdai, mui); init_(); }

    // constructor for a non-homogeneous material
    mdbrick_mixed_isotropic_linearized_plate
    (mesh_fem &mf_ut_, mesh_fem &mf_u3_, mesh_fem &mf_theta_,
     mesh_fem &mf_data_, const VECTOR &lambdai, const VECTOR &mui,
     double epsilon_, bool sym = false, bool mat_stored = false)
      : mf_ut(mf_ut_), mf_u3(mf_u3_), mf_theta(mf_theta_), mf_data(mf_data_),
	epsilon(epsilon_), matrix_stored(mat_stored), symmetrized(sym)
    { set_Lame_coeff(lambdai, mui); init_(); }
 
  };


  /* ******************************************************************** */
  /*		plate source term model brick.                            */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_source_term : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;

    mdbrick_source_term<MODEL_STATE> *ut_part,*u3_part,*phi_part,*sub_problem;

  public :

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem->mixed_variables(b, i0); }
    virtual size_type nb_constraints(void) 
    { return sub_problem->nb_constraints(); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0 = 0, bool modified=false)
    { sub_problem->compute_tangent_matrix(MS, i0, j0, modified); }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0=0,
				size_type j0=0)
    { sub_problem->compute_residu(MS, i0, j0); }

    mdbrick_plate_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			      mesh_fem &mf_data, const VECTOR &B,
			      size_type bound = size_type(-1),
			      size_type num_fem = 0) {
      ut_part = phi_part = u3_part = 0;
      bool mixed = false, symmetrized = false;
      if (problem.get_mesh_fem_info(num_fem).brick_ident
	  == MDBRICK_LINEAR_PLATE)
	{ mixed = false; symmetrized = false; } 
      else if (problem.get_mesh_fem_info(num_fem).brick_ident 
	       == MDBRICK_MIXED_LINEAR_PLATE) {
	mixed=true;
	symmetrized = ((problem.get_mesh_fem_info(num_fem).info) & 2);
      }
      else DAL_THROW(failure_error,
		     "This brick should only be applied to a plate problem");
      if ((!(problem.get_mesh_fem_info(num_fem).info & 1))
	  || (num_fem + (mixed ? 4 : 2) >= problem.nb_mesh_fems()))
	DAL_THROW(failure_error, "The mesh_fem number is not correct");

      size_type n = gmm::vect_size(B) / 3;
      VECTOR Bt(2*n);
      gmm::copy(gmm::sub_vector(B, gmm::sub_slice(0, n, 3)),
		gmm::sub_vector(Bt, gmm::sub_slice(0, n, 2)));
      gmm::copy(gmm::sub_vector(B, gmm::sub_slice(1, n, 3)),
		gmm::sub_vector(Bt, gmm::sub_slice(1, n, 2)));
      ut_part = sub_problem = new mdbrick_source_term<MODEL_STATE>
	(problem, mf_data, Bt, bound, num_fem);
      VECTOR Bn(n);
      gmm::copy(gmm::sub_vector(B, gmm::sub_slice(2, n, 3)), Bn);
      if (!mixed || symmetrized)
	sub_problem = u3_part = new mdbrick_source_term<MODEL_STATE>
	  (*ut_part, mf_data, Bn, bound, num_fem+1);
      
      if (mixed && !symmetrized)
	sub_problem = phi_part = new mdbrick_source_term<MODEL_STATE>
	  (*sub_problem, mf_data, Bn, bound, num_fem+4);

      this->add_sub_brick(*sub_problem);
      if (bound != size_type(-1)) {
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
	this->add_proper_boundary_info(num_fem+1, bound, MDBRICK_NEUMANN);
      }
      this->update_from_context();
    }

    ~mdbrick_plate_source_term() {
      delete ut_part;
      if (u3_part) delete u3_part;
      if (phi_part) delete phi_part;
    }
    
  };

  /* ******************************************************************** */
  /*		Simple support condition for plate model brick.           */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_simple_support : public mdbrick_abstract<MODEL_STATE>  {
    
    mdbrick_Dirichlet<MODEL_STATE> *ut_part, *u3_part;
    mdbrick_Dirichlet<MODEL_STATE> *phi_part, *sub_problem;

  public :

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem->mixed_variables(b, i0); }
    virtual size_type nb_constraints(void) 
    { return sub_problem->nb_constraints(); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0 = 0, bool modified=false)
    { sub_problem->compute_tangent_matrix(MS, i0, j0, modified); }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0=0,size_type j0=0)
    { sub_problem->compute_residu(MS, i0, j0); }

    mdbrick_plate_simple_support(mdbrick_abstract<MODEL_STATE> &problem,
				 mesh_fem &mf_data, size_type bound,
				 size_type num_fem = 0,
				 bool with_mult = false) : phi_part(0) {
      ut_part = new  mdbrick_Dirichlet<MODEL_STATE>
	(problem, mf_data, bound, num_fem, with_mult);
      u3_part = new  mdbrick_Dirichlet<MODEL_STATE>
	(*ut_part, mf_data, bound, num_fem+1, with_mult);
      bool mixed = false, symmetrized = false;
      if (problem.get_mesh_fem_info(num_fem).brick_ident
	  == MDBRICK_LINEAR_PLATE)
	{ mixed = false; symmetrized = false; } 
      else if (problem.get_mesh_fem_info(num_fem).brick_ident 
	       == MDBRICK_MIXED_LINEAR_PLATE) {
	mixed=true;
	symmetrized = ((problem.get_mesh_fem_info(num_fem).info) & 2);
      }
      else DAL_THROW(failure_error,
		     "This brick should only be applied to a plate problem");
      if ((!(problem.get_mesh_fem_info(num_fem).info & 1))
	  || (num_fem + (mixed ? 4 : 2) >= problem.nb_mesh_fems()))
	DAL_THROW(failure_error, "The mesh_fem number is not correct");

      if (mixed)
	sub_problem = phi_part = new  mdbrick_Dirichlet<MODEL_STATE>
	  (*u3_part, mf_data, bound, num_fem+4, with_mult);
      else sub_problem = u3_part;
      this->add_sub_brick(*sub_problem);
      this->add_proper_boundary_info(num_fem, bound, MDBRICK_SIMPLE_SUPPORT);
      this->add_proper_boundary_info(num_fem+1, bound, MDBRICK_SIMPLE_SUPPORT);
      this->add_proper_boundary_info(num_fem+2, bound, MDBRICK_SIMPLE_SUPPORT);
      this->update_from_context();
    }

    virtual ~mdbrick_plate_simple_support() {
      delete ut_part; delete u3_part;
      if (phi_part) delete phi_part;
    }
    
  };

  /* ******************************************************************** */
  /*		Clamped condition for plate model brick.                  */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_clamped_support: public mdbrick_abstract<MODEL_STATE> {
    
    mdbrick_Dirichlet<MODEL_STATE> ut_part, u3_part, theta_part;
    mdbrick_Dirichlet<MODEL_STATE> *phi_part, *sub_problem;

  public :

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0)
    { sub_problem->mixed_variables(b, i0); }
    virtual size_type nb_constraints(void) 
    { return sub_problem->nb_constraints(); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0 = 0, bool modified=false)
    { sub_problem->compute_tangent_matrix(MS, i0, j0, modified); }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0=0,
				size_type j0=0)
    { sub_problem->compute_residu(MS, i0, j0); }

    mdbrick_plate_clamped_support(mdbrick_abstract<MODEL_STATE> &problem,
				  mesh_fem &mf_data, size_type bound,
				  size_type num_fem = 0,
				  bool with_mult = false)
      : ut_part(problem, mf_data, bound, num_fem, with_mult),
	u3_part(ut_part, mf_data, bound, num_fem+1, with_mult),
	theta_part(u3_part, mf_data, bound, num_fem+2, with_mult),
	phi_part(0) {

      bool mixed = false, symmetrized = false;
      if (problem.get_mesh_fem_info(num_fem).brick_ident
	  == MDBRICK_LINEAR_PLATE)
	{ mixed = false; symmetrized = false; } 
      else if (problem.get_mesh_fem_info(num_fem).brick_ident 
	       == MDBRICK_MIXED_LINEAR_PLATE) {
	mixed=true;
	symmetrized = ((problem.get_mesh_fem_info(num_fem).info) & 2);
      }
      else DAL_THROW(failure_error,
		     "This brick should only be applied to a plate problem");
      if ((!(problem.get_mesh_fem_info(num_fem).info & 1))
	  || (num_fem + (mixed ? 4 : 2) >= problem.nb_mesh_fems()))
	DAL_THROW(failure_error, "The mesh_fem number is not correct");

      if (mixed) {
	sub_problem = phi_part = new  mdbrick_Dirichlet<MODEL_STATE>
	  (theta_part, mf_data, bound, num_fem+4, with_mult);
	this->add_sub_brick(*phi_part);
      }
      else { 
	this->add_sub_brick(theta_part);
	sub_problem = &theta_part;
      }
      this->add_proper_boundary_info(num_fem, bound, MDBRICK_CLAMPED_SUPPORT);
      this->add_proper_boundary_info(num_fem+1, bound, MDBRICK_CLAMPED_SUPPORT);
      this->add_proper_boundary_info(num_fem+2, bound, MDBRICK_CLAMPED_SUPPORT);

      this->update_from_context();
    }

    ~mdbrick_plate_clamped_support() { if (phi_part) delete phi_part; }
    
  };


  /* ******************************************************************** */
  /*		Free edges condition for mixed plate model brick.         */
  /* ******************************************************************** */

  template<class VEC>
  void asm_constraint_on_theta(const VEC &V, const mesh_fem &mf_theta,
			       size_type boundary) {
    generic_assembly assem("t=comp(vBase(#1).Normal());"
			   "V(#1)+= t(:,2,1) - t(:,1,2);");
    assem.push_mf(mf_theta);
    assem.push_vec(const_cast<VEC &>(V));
    assem.boundary_assembly(boundary);
  }

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_closing: public mdbrick_abstract<MODEL_STATE> {
 
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;
    
    mdbrick_abstract<MODEL_STATE> *sub_problem;
    mesh_fem *mf_theta;
    gmm::row_matrix<gmm::rsvector<value_type> > CO;
    size_type num_fem;
    bool mixed, symmetrized, allclamped, with_multipliers;

    void compute_constraints(void) {
      if (!mixed) {  allclamped = false; gmm::resize(CO, 0, 0); return; }
      allclamped = true;
      std::vector<size_type> cv_nums;
      std::vector<short_type> face_nums;
      
      getfem_mesh *mesh = &(mf_theta->linked_mesh());
    
      getfem::convex_face_ct border_faces;
      getfem::outer_faces_of_mesh(*mesh, border_faces);
      dal::bit_vector vb = (this->mesh_fems[num_fem])->get_valid_boundaries();
      
      for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
	   it != border_faces.end(); ++it) {
	bool add = true;
	// cout << "face " << it->f << " of cv " << it->cv << "boundaries : ";
	for (dal::bv_visitor i(vb); !i.finished(); ++i) {
	  if ((this->mesh_fems[num_fem])->is_face_on_boundary(i,it->cv,it->f)) {
	    // cout << i << endl;
	    bound_cond_type bct = this->boundary_type(num_fem, i);
	    if (bct != MDBRICK_UNDEFINED && bct != MDBRICK_NEUMANN) add = false;
	    if (bct != MDBRICK_CLAMPED_SUPPORT) allclamped = false;
	  }
	}
	
	if (add) {
	  cv_nums.push_back(it->cv); face_nums.push_back(it->f); allclamped = false;
	  // cout << " adding";
	}
	// cout << endl;
      }
      cout << "allclamped = " << allclamped << endl;

      std::vector<size_type> comp_conns(cv_nums.size(), size_type(-1)); 
      size_type nbmax = mesh->points().ind_last() + 1, p1, p2;
      std::vector<size_type> E1(nbmax, size_type(-1)), E2(nbmax, size_type(-1));
      for (size_type j = 0; j < cv_nums.size(); ++j) {
	p1 = mesh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[0];
	p2 = mesh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[1];
	if (E1[p1] == size_type(-1)) E1[p1] = j; else E2[p1] = j;
	if (E1[p2] == size_type(-1)) E1[p2] = j; else E2[p2] = j;	
      }

      size_type comp_conn = 0;
      for (size_type i = 0; i < comp_conns.size(); ++i) {
	if (comp_conns[i] == size_type(-1)) {
	  
	  comp_conns[i] = comp_conn;
	  p1 = mesh->ind_points_of_face_of_convex(cv_nums[i],face_nums[i])[0];
	  p2 = mesh->ind_points_of_face_of_convex(cv_nums[i],face_nums[i])[1];
	  size_type j1 = (E1[p1] == i) ? E2[p1] :  E1[p1];
	  size_type j2 = (E1[p2] == i) ? E2[p2] :  E1[p2];
	  
	  for (unsigned k = 0; k < 2; ++k) {
	    size_type j = (k == 0) ? j1 : j2;
	    
	    while (j != size_type(-1) && comp_conns[j] == size_type(-1)) {
	      comp_conns[j] = comp_conn;
	      p1 = mesh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[0];
	      p2 = mesh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[1];
	      size_type i1 = (E1[p1] == j) ? E2[p1] :  E1[p1];
	      size_type i2 = (E1[p2] == j) ? E2[p2] :  E1[p2];
	      if (i1 == size_type(-1) || comp_conns[i1] != size_type(-1))
		j = i2; else j = i1;
	    }
	  }
	  
	  ++comp_conn;
	}
      }

      cout << "Number of comp conn : " << comp_conn << endl;
      if (comp_conn < 2) {
	gmm::resize(CO, 0, 0);
      }
      else {
	vb = mf_theta->get_valid_boundaries();
	size_type boundary = vb.first_false();
	gmm::resize(CO, comp_conn, mf_theta->nb_dof());
	for (size_type k = 0; k < comp_conn; ++k) {
	  for (size_type i = 0; i < comp_conns.size(); ++i)
	    if (comp_conns[i] == k)
	      mf_theta->add_boundary_elt(boundary, cv_nums[i], face_nums[i]);

	  std::vector<value_type> V(mf_theta->nb_dof());
	  asm_constraint_on_theta(V, *mf_theta, boundary);
	  gmm::copy(V, gmm::mat_row(CO, k));
	  mf_theta->sup_boundary(boundary);
	}
	cout << "CO = " << CO << endl;
      }
      this->computed();
    }

  public :

    virtual void mixed_variables(dal::bit_vector &b, size_type i0 = 0) {
      this->context_check(); if (this->to_be_computed()) compute_constraints();
      sub_problem->mixed_variables(b, i0);
      if (with_multipliers)
	b.add(i0+sub_problem->nb_dof(), gmm::mat_nrows(CO)+(allclamped ? 1:0));
    }
    virtual size_type nb_constraints(void)  {
      this->context_check(); if (this->to_be_computed()) compute_constraints();
      if (with_multipliers)
	return sub_problem->nb_constraints();
      else 
	return sub_problem->nb_constraints() 
	  + gmm::mat_nrows(CO) + (allclamped ? 1:0);
    }
    virtual size_type nb_dof(void) {
      this->context_check(); if (this->to_be_computed()) compute_constraints();
      if (with_multipliers)
	return sub_problem->nb_dof() + gmm::mat_nrows(CO) + (allclamped ? 1:0);
      return sub_problem->nb_dof();
    }

    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0,
					size_type j0=0, bool modified=false) {
      sub_problem->compute_tangent_matrix(MS, i0, j0, modified);
      react(MS, i0, modified);
      if (this->to_be_computed()) compute_constraints();
      if (this->to_be_transferred()) {
	gmm::sub_interval SUBJ(i0+this->mesh_fem_positions[num_fem+2],
			       mf_theta->nb_dof());
	size_type nbd = sub_problem->nb_dof();
	if (with_multipliers) {
	  if (gmm::mat_nrows(CO) > 0) {
	    gmm::sub_interval SUBI(i0+nbd, gmm::mat_nrows(CO));
	    gmm::copy(CO, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	    gmm::copy(gmm::transposed(CO),
		      gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	  }
	  if (allclamped) {
	    size_type i = i0 + nbd + gmm::mat_nrows(CO);
	    size_type j = i0 + this->mesh_fem_positions[num_fem+3];
	    MS.tangent_matrix()(i, j) = value_type(1);
	    MS.tangent_matrix()(j, i) = value_type(1);
	  }
	}
	else {
	  size_type ncs = sub_problem->nb_constraints();
	  if (gmm::mat_nrows(CO) > 0) {
	    gmm::sub_interval SUBI(j0 + ncs, gmm::mat_nrows(CO));
	    gmm::copy(CO,gmm::sub_matrix(MS.constraints_matrix(),SUBI,SUBJ));
	  }
	  if (allclamped) {
	    MS.constraints_matrix()(j0+ncs+gmm::mat_nrows(CO),
				    i0 + this->mesh_fem_positions[num_fem+3])
	      = value_type(1);
	  }
	}
	this->transferred();
      }
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0=0,
				size_type j0=0) {
      sub_problem->compute_residu(MS, i0, j0);
      react(MS, i0, false);
      if (this->to_be_computed()) compute_constraints();
      gmm::sub_interval SUBJ(i0+this->mesh_fem_positions[num_fem+2],
			       mf_theta->nb_dof());
      if (with_multipliers) {
	  size_type nbd = sub_problem->nb_dof();
	if (gmm::mat_nrows(CO) > 0) {
	  gmm::sub_interval SUBI(i0 + nbd, gmm::mat_nrows(CO));
	  gmm::mult(CO, gmm::sub_vector(MS.state(), SUBJ),
		    gmm::sub_vector(MS.residu(), SUBI));
	  gmm::mult_add(gmm::transposed(CO), gmm::sub_vector(MS.state(), SUBI),
			gmm::sub_vector(MS.residu(), SUBJ));
	}
	if (allclamped) {
	  size_type i = i0 + nbd + gmm::mat_nrows(CO);
	  size_type j = i0 + this->mesh_fem_positions[num_fem+3];
	  MS.residu()[i] = MS.state()[j];
	  MS.residu()[j] += MS.state()[i];
	}
      }
      else {
	size_type ncs = sub_problem->nb_constraints();
	if (gmm::mat_nrows(CO) > 0) {
	  gmm::sub_interval SUBI(j0+ncs,gmm::mat_nrows(CO));
	  gmm::mult(CO, gmm::scaled(gmm::sub_vector(MS.state(), SUBJ), 
				    value_type(-1)),
		    gmm::sub_vector(MS.constraints_rhs(), SUBI));
	}
	if (allclamped) {
	  (MS.constraints_rhs())[j0+ncs+gmm::mat_nrows(CO)] =
	    -(MS.state())[i0 + this->mesh_fem_positions[num_fem+3]];
	}
      }
    }

    mdbrick_plate_closing(mdbrick_abstract<MODEL_STATE> &problem,
			  size_type num_fem_ = 0, int with_mult = -1)
      : sub_problem(&problem), num_fem(num_fem_),
	with_multipliers(with_mult!=0) {

      if (with_mult == -1)
	with_multipliers = (sub_problem->nb_constraints() == 0);

      mixed = false; symmetrized = false;
      if (problem.get_mesh_fem_info(num_fem).brick_ident
	  == MDBRICK_LINEAR_PLATE)
	{ mixed = false; symmetrized = false; } 
      else if (problem.get_mesh_fem_info(num_fem).brick_ident 
	       == MDBRICK_MIXED_LINEAR_PLATE) {
	mixed=true;
	symmetrized = ((problem.get_mesh_fem_info(num_fem).info) & 2);
      }
      else DAL_THROW(failure_error,
		     "This brick should only be applied to a plate problem");
      if ((!(problem.get_mesh_fem_info(num_fem).info & 1))
	  || (num_fem + (mixed ? 4 : 2) >= problem.nb_mesh_fems()))
	DAL_THROW(failure_error, "The mesh_fem number is not correct");


      this->add_sub_brick(problem);
      this->update_from_context();
      mf_theta = this->mesh_fems[num_fem+2];
      compute_constraints();
    }
    
  };





}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LINEARIZED_PLATES_H__ */
