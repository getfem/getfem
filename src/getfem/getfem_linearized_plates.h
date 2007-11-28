// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================
/**@file getfem_linearized_plates.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, Michel Salaun, msalaun@ensica.fr
   @date November 1, 2004.
   @brief Define a linear plate model brick.
*/

#ifndef GETFEM_LINEARIZED_PLATES_H__
#define GETFEM_LINEARIZED_PLATES_H__

#include "getfem_modeling.h"

namespace getfem {

  /* ******************************************************************** */
  /*		Linear plate specific assembly procedures.                */
  /* ******************************************************************** */

  /**@ingroup asm*/
  template<class MAT, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear
  (const MAT &RM, const mesh_im &mim, const mesh_fem &mf_u3,
   const mesh_fem &mf_theta, const mesh_fem &mfdata, const VECT &MU,
   const mesh_region &rg = mesh_region::all_convexes()) {
    gmm::sub_interval I1(0, mf_u3.nb_dof());
    gmm::sub_interval I2(mf_u3.nb_dof(), mf_theta.nb_dof());
    
    asm_stiffness_matrix_for_plate_transverse_shear
      (gmm::sub_matrix(RM, I1), gmm::sub_matrix(RM, I1, I2),
       gmm::transposed(gmm::sub_matrix(RM, I2, I1)),
       gmm::sub_matrix(RM, I2), mim, mf_u3, mf_theta, mfdata, MU, rg);
  }

  /**@ingroup asm*/
  template<class MAT, class MAT3, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear
  (const MAT &RM1, const MAT &RM2, const MAT3 &RM3, const MAT &RM4,
   const mesh_im &mim, const mesh_fem &mf_u3, const mesh_fem &mf_theta,
   const mesh_fem &mfdata, const VECT &MU,
   const mesh_region &rg = mesh_region::all_convexes()) {
    GMM_ASSERT1(mfdata.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
    
    GMM_ASSERT1(mf_u3.get_qdim() == 1 && mf_theta.get_qdim() == 2,
		"wrong qdim for the mesh_fem");
    generic_assembly assem("mu=data$1(#3);"
			   "t1=comp(Grad(#1).Grad(#1).Base(#3));"
			   "M$1(#1,#1)+=sym(t1(:,i,:,i,j).mu(j));"
			   "t2=comp(vBase(#2).vBase(#2).Base(#3));"
			   "M$4(#2,#2)+=sym(t2(:,i,:,i,j).mu(j));"
			   "t3=comp(Grad(#1).vBase(#2).Base(#3));"
			   "M$2(#1,#2)+=t3(:,i,:,i,j).mu(j);"
			   "M$3(#1,#2)+=t3(:,i,:,i,j).mu(j);");
    assem.push_mi(mim);
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mf(mfdata);
    assem.push_data(MU);
    assem.push_mat(const_cast<MAT &>(RM1));
    assem.push_mat(const_cast<MAT &>(RM2));
    assem.push_mat(const_cast<MAT3 &>(RM3));
    assem.push_mat(const_cast<MAT &>(RM4));
    assem.assembly(rg);
  }
  
  /**@ingroup asm*/
  template<class MAT, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear_mitc
  (const MAT &RM, const mesh_im &mim, const mesh_fem &mf_u3,
   const mesh_fem &mf_theta, const mesh_fem &mfdata, const VECT &MU,
   const mesh_region &rg = mesh_region::all_convexes()) {
    gmm::sub_interval I1(0, mf_u3.nb_dof());
    gmm::sub_interval I2(mf_u3.nb_dof(), mf_theta.nb_dof());
    
    asm_stiffness_matrix_for_plate_transverse_shear_mitc_new
      (gmm::sub_matrix(RM, I1), gmm::sub_matrix(RM, I1, I2),
       gmm::transposed(gmm::sub_matrix(RM, I2, I1)),
       gmm::sub_matrix(RM, I2), mim, mf_u3, mf_theta, mfdata, MU, rg);
  }




  class mitc4_projection_term  : public getfem::nonlinear_elem_term {

    bgeot::multi_index sizes_;

  public:

     mitc4_projection_term(void) : sizes_(8,8)  { }
     const bgeot::multi_index &sizes() const {  return sizes_; }
     virtual void compute(getfem::fem_interpolation_context & ctx,
			  bgeot::base_tensor &t) {
       
       //     ctx.G()  --> coordonées des noeuds
       //     ctx.B()  --> (grad tau) ^{-T}

       for (size_type i = 0; i < 8; ++i)
	  for (size_type j = 0; j < 8; ++j) 
	      t(i, j) = 0 ;  // Initialisation
       
//        // Remplissage des termes non nuls
	// lignes 0 et 1
	t(0,0) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(0,0) + ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(0,1) ;
	t(0,1) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(0,0) + ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(0,1) ;
	t(0,2) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(0,0) ; 
	t(0,3) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(0,0) ;
	t(0,4) = ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(0,1) ;
	t(0,5) = ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(0,1) ;
	
	t(1,0) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(1,0) + ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(1,1) ;
	t(1,1) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(1,0) + ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(1,1) ;
	t(1,2) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(1,0) ; 
	t(1,3) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(1,0) ;
	t(1,4) = ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(1,1) ;
	t(1,5) = ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(1,1) ;
	// lignes 2 et 3
	t(2,0) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(0,0) ; 
	t(2,1) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(0,0) ;
	t(2,2) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(0,0) + ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(0,1) ;
	t(2,3) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(0,0) + ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(0,1) ;
	t(2,6) = ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(0,1) ;
	t(2,7) = ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(0,1) ;
	
	t(3,0) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(1,0) ; 
	t(3,1) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(1,0) ;
	t(3,2) = ( ctx.G()(0,1) - ctx.G()(0,0) ) * ctx.B()(1,0) + ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(1,1) ;
	t(3,3) = ( ctx.G()(1,1) - ctx.G()(1,0) ) * ctx.B()(1,0) + ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(1,1) ;
	t(3,6) = ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(1,1) ;
	t(3,7) = ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(1,1) ;
	// lignes 4 et 5
	t(4,0) = ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(0,1) ; 
	t(4,1) = ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(0,1) ;
	t(4,4) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(0,0) + ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(0,1) ;
	t(4,5) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(0,0) + ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(0,1) ;
	t(4,6) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(0,0) ;
	t(4,7) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(0,0) ;
	
	t(5,0) = ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(1,1) ; 
	t(5,1) = ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(1,1) ;
	t(5,4) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(1,0) + ( ctx.G()(0,2) - ctx.G()(0,0) ) * ctx.B()(1,1) ;
	t(5,5) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(1,0) + ( ctx.G()(1,2) - ctx.G()(1,0) ) * ctx.B()(1,1) ;
	t(5,6) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(1,0) ;
	t(5,7) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(1,0) ;
	// lignes 6 et 7
	t(6,2) = ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(0,1) ; 
	t(6,3) = ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(0,1) ;
	t(6,4) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(0,0) ;
	t(6,5) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(0,0) ;
	t(6,6) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(0,0) + ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(0,1) ; 
	t(6,7) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(0,0) + ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(0,1) ;
	
	t(7,2) = ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(1,1) ;
	t(7,3) = ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(1,1) ;
	t(7,4) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(1,0) ;
	t(7,5) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(1,0) ;
	t(7,6) = ( ctx.G()(0,3) - ctx.G()(0,2) ) * ctx.B()(1,0) + ( ctx.G()(0,3) - ctx.G()(0,1) ) * ctx.B()(1,1) ;
	t(7,7) = ( ctx.G()(1,3) - ctx.G()(1,2) ) * ctx.B()(1,0) + ( ctx.G()(1,3) - ctx.G()(1,1) ) * ctx.B()(1,1) ;
	
        t *= 0.5 ;

	    }
      
    
  };


  /**@ingroup asm*/
  template<class MAT, class MAT3, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear_mitc_new
  (const MAT &RM1, const MAT &RM2, const MAT3 &RM3, const MAT &RM4,
   const mesh_im &mim, const mesh_fem &mf_u3, const mesh_fem &mf_theta,
   const mesh_fem &mfdata, const VECT &MU,
   const mesh_region &rg = mesh_region::all_convexes()) {
    typedef typename gmm::linalg_traits<VECT>::value_type value_type;

    GMM_ASSERT1(mfdata.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
    
    GMM_ASSERT1(mf_u3.get_qdim() == 1 && mf_theta.get_qdim() == 2,
		"wrong qdim for the mesh_fem");

    mitc4_projection_term mitc4;

    generic_assembly assem("mu=data$1(#3);"
			   "t1=comp(Grad(#1).Grad(#1).Base(#3));"
			   "M$1(#1,#1)+=sym(t1(:,i,:,i,j).mu(j));"
			   "M$4(#2,#2)+=sym(comp(NonLin(#2)(k,:).vBase(#2)(k,i).vBase(#2)(l,i).Base(#3)(:).NonLin(#2)(l,:))(:,j,:).mu(j));"
			   "M$2(#1,#2)+=comp(Grad(#1)(:,i).vBase(#2)(l,i).Base(#3)(:).NonLin(#2)(l,:))(:,j,:).mu(j);"   
			   "M$3(#1,#2)+=comp(Grad(#1)(:,i).vBase(#2)(l,i).Base(#3)(:).NonLin(#2)(l,:))(:,j,:).mu(j);"
			   );

    assem.push_mi(mim);
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mf(mfdata);
    assem.push_data(MU);
    assem.push_nonlinear_term(&mitc4);
    assem.push_mat(const_cast<MAT &>(RM1));
    assem.push_mat(const_cast<MAT &>(RM2));
    assem.push_mat(const_cast<MAT3 &>(RM3));
    assem.push_mat(const_cast<MAT &>(RM4));
    assem.assembly(rg);
    //cout << "RM3 = " << RM3 << endl; getchar();
  }

  /**@ingroup asm*/
  template<class MAT, class MAT3, class VECT>
  void asm_stiffness_matrix_for_plate_transverse_shear_mitc
  (const MAT &RM1, const MAT &RM2, const MAT3 &RM3, const MAT &RM4,
   const mesh_im &mim, const mesh_fem &mf_u3, const mesh_fem &mf_theta,
   const mesh_fem &mfdata, const VECT &MU,
   const mesh_region &rg = mesh_region::all_convexes()) {
    typedef typename gmm::linalg_traits<VECT>::value_type value_type;

    GMM_ASSERT1(mfdata.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
    
    GMM_ASSERT1(mf_u3.get_qdim() == 1 && mf_theta.get_qdim() == 2,
		"wrong qdim for the mesh_fem");

    generic_assembly assem("mu=data$1(#3);"
			   "A=data$2(8,8);"
			   "t1=comp(Grad(#1).Grad(#1).Base(#3));"
			   "M$1(#1,#1)+=sym(t1(:,i,:,i,j).mu(j));"
			   "t2=comp(vBase(#2).vBase(#2).Base(#3));"
			   "M$4(#2,#2)+=sym(A(k,:).t2(k,i,l,i,j).mu(j).A(l,:));"
			   "t3=comp(Grad(#1).vBase(#2).Base(#3));"
			   "M$2(#1,#2)+=t3(:,i,l,i,j).mu(j).A(l,:);"
			   "M$3(#1,#2)+=t3(:,i,l,i,j).mu(j).A(l,:);"
			   );

    std::vector<value_type> A(64);
    // remplissage de A :
    std::fill(A.begin(), A.end(), 0.) ;
    A[ 0] = 0.5 ;   A[16] = 0.5 ;   A[36] = 0.5 ;   A[52] = 0.5 ;
    A[ 2] = 0.5 ;   A[18] = 0.5 ;   A[38] = 0.5 ;   A[54] = 0.5 ;
    A[ 9] = 0.5 ;   A[27] = 0.5 ;   A[41] = 0.5 ;   A[59] = 0.5 ;
    A[13] = 0.5 ;   A[31] = 0.5 ;   A[45] = 0.5 ;   A[63] = 0.5 ;
    
    assem.push_mi(mim);
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mf(mfdata);
    assem.push_data(MU);
    assem.push_data(A);
    assem.push_mat(const_cast<MAT &>(RM1));
    assem.push_mat(const_cast<MAT &>(RM2));
    assem.push_mat(const_cast<MAT3 &>(RM3));
    assem.push_mat(const_cast<MAT &>(RM4));
    assem.assembly(rg);
    //cout << "RM3 = " << RM3 << endl; getchar();
  }


  /* ******************************************************************** */
  /*		Linear plate model brick.                                 */
  /* ******************************************************************** */

# define MDBRICK_LINEAR_PLATE 897523

  /**
     Linear plate model brick (for moderately thick plates, using the
     Reissner-Mindlin model).

     For very thin plates, see the
     mdbrick_mixed_isotropic_linearized_plate or the
     mdbrick_bilaplacian.  @ingroup bricks
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_isotropic_linearized_plate
    : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    const mesh_im &mim, &mim_subint;
    const mesh_fem &mf_ut, &mf_u3, &mf_theta;
    mdbrick_parameter<VECTOR> lambda_, mu_;
    value_type epsilon;
    bool homogeneous, mitc, K_uptodate;
    T_MATRIX K;
    size_type nbdof;

    virtual void proper_update(void) {
      K_uptodate = false;
      nbdof = mf_ut.nb_dof() + mf_u3.nb_dof() + mf_theta.nb_dof();
    }

  public :

    const T_MATRIX &get_K(void) {
      this->context_check(); 
      if (!K_uptodate || this->parameters_is_any_modified()) {
	GMM_ASSERT1(&lambda_.mf() == &mu_.mf(),
		    "lambda and mu should share the same mesh_fem");
	gmm::resize(K, nbdof, nbdof);
	gmm::clear(K);
	gmm::sub_interval I1(0, mf_ut.nb_dof());
	gmm::sub_interval I2(mf_ut.nb_dof(), mf_u3.nb_dof()+mf_theta.nb_dof());
	gmm::sub_interval I3(mf_ut.nb_dof()+mf_u3.nb_dof(),mf_theta.nb_dof());
	VECTOR vlambda(lambda_.get()), vmu(mu_.get());
	gmm::scale(vlambda, value_type(2) * epsilon);
	gmm::scale(vmu, value_type(2) * epsilon);
	asm_stiffness_matrix_for_linear_elasticity
	  (gmm::sub_matrix(K, I1), mim, mf_ut, lambda_.mf(), vlambda, vmu,
	   mf_ut.linked_mesh().get_mpi_region());
	// gmm::scale(mu, value_type(1) / value_type(2));
	if (mitc) 
	  asm_stiffness_matrix_for_plate_transverse_shear_mitc
	    (gmm::sub_matrix(K, I2), mim_subint, mf_u3, mf_theta, lambda_.mf(),
	     vmu, mf_ut.linked_mesh().get_mpi_region());
	else
	  asm_stiffness_matrix_for_plate_transverse_shear
	    (gmm::sub_matrix(K, I2), mim_subint, mf_u3, mf_theta, lambda_.mf(),
	     vmu, mf_ut.linked_mesh().get_mpi_region());
	gmm::scale(vlambda, epsilon * epsilon / value_type(3));
	// gmm::scale(mu, value_type(2) * epsilon * epsilon / value_type(3));
	gmm::scale(vmu, epsilon * epsilon / value_type(3));
	asm_stiffness_matrix_for_linear_elasticity
	  (gmm::sub_matrix(K, I3), mim, mf_theta, lambda_.mf(), vlambda, vmu,
	   mf_ut.linked_mesh().get_mpi_region());
	K_uptodate = true;
	this->parameters_set_uptodate();
      }
      return K;
    }

    mdbrick_parameter<VECTOR> &lambda(void) { return lambda_; }
    const mdbrick_parameter<VECTOR> &lambda(void) const { return lambda_; }
    mdbrick_parameter<VECTOR> &mu(void) { return mu_; }
    const mdbrick_parameter<VECTOR> &mu(void) const { return mu_; }


    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0, nbdof);
      gmm::copy(get_K(), gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				size_type) {
      gmm::sub_interval SUBI(i0, nbdof);
      gmm::mult(get_K(), gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residual(), SUBI));
    }

    void set_elastic_coeff(value_type E, value_type nu) {
      lambda_.set(nu * E / (value_type(1) - nu*nu));
      mu_.set(E/(value_type(2)*(value_type(1)+nu)));
    }

    void set_mitc(void) {
      if (!mitc) {
	mitc = true;  K_uptodate = false;
	this->change_context();
      }
    }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), nbdof);
      return gmm::sub_vector(MS.state(), SUBU);
    }
    SUBVECTOR get_ut(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), mf_ut.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }
    SUBVECTOR get_u3(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index() + mf_ut.nb_dof(),
			     mf_u3.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }
    SUBVECTOR get_theta(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index() + mf_ut.nb_dof()
			     + mf_u3.nb_dof(), mf_theta.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      GMM_ASSERT1(mf_ut.get_qdim() == 2, "Qdim of mf_ut should be 2.");
      GMM_ASSERT1(mf_u3.get_qdim() == 1, "Qdim of mf_u3 should be 1.");
      GMM_ASSERT1(mf_theta.get_qdim() == 2, "Qdim of mf_theta should be 2.");
      mitc = false;
      this->add_proper_mesh_im(mim);
      this->add_proper_mesh_im(mim_subint);
      this->add_proper_mesh_fem(mf_ut, MDBRICK_LINEAR_PLATE, 1);
      this->add_proper_mesh_fem(mf_u3, MDBRICK_LINEAR_PLATE, 0);
      this->add_proper_mesh_fem(mf_theta, MDBRICK_LINEAR_PLATE, 0);
      this->force_update();
    }

    /** 
	Constructor for a homogeneous material (constant lambda and
	mu). In order to avoid locking phenomena it is preferable for
	mf_theta to be of lower order compared to mf_u3. Another
	possibility is to use a sub-integration of the transverse
	shear term (with mim_subint).

	@param mf_ut the finite element method for the membrane displacement
	@param mf_u3 the finite element method for the transverse displacement.
	@param mf_theta the finite element method for the rotation of the normal (section rotations).
	@param epsilon the thickness of the plate.
    */
    mdbrick_isotropic_linearized_plate
    (const mesh_im &mim_, const mesh_fem &mf_ut_, const mesh_fem &mf_u3_,
     const mesh_fem &mf_theta_,
     value_type lambdai, value_type mui, double epsilon_)
      : mim(mim_), mim_subint(mim_), mf_ut(mf_ut_), mf_u3(mf_u3_),
	mf_theta(mf_theta_), lambda_("lambda", mf_ut_.linked_mesh(), this),
	mu_("mu", mf_ut_.linked_mesh(), this), epsilon(epsilon_)
    { lambda_.set(lambdai); mu_.set(mui); init_(); }

    /** Constructor for a homogeneous material (constant lambda and
	mu) with sub integration for the transverse shear term.

	@param mf_ut the finite element method for the membrane displacement
	@param mf_u3 the finite element method for the transverse displacement.
	@param mf_theta the finite element method for the rotation of the normal (section rotations).
	@param epsilon the thickness of the plate.
     */
    mdbrick_isotropic_linearized_plate
    (const mesh_im &mim_, const mesh_im &mim_subint_, const mesh_fem &mf_ut_,
     const mesh_fem &mf_u3_, const mesh_fem &mf_theta_, value_type lambdai,
     value_type mui, double epsilon_)
      : mim(mim_), mim_subint(mim_subint_), mf_ut(mf_ut_), mf_u3(mf_u3_),
	mf_theta(mf_theta_), lambda_("lambda", mf_ut_.linked_mesh(), this),
	mu_("mu", mf_ut_.linked_mesh(), this), epsilon(epsilon_)
    {  lambda_.set(lambdai); mu_.set(mui); init_(); }
 
  };


  /* ******************************************************************** */
  /*		Mixed linear plate specific assembly procedures.          */
  /* ******************************************************************** */

  /**@ingroup asm*/
  template<class MAT>
  void asm_coupling_u3theta(const MAT &RM, const mesh_im &mim,
			    const mesh_fem &mf_u3,
			    const mesh_fem &mf_theta,
			    const mesh_region &rg
			    = mesh_region::all_convexes()) {
    
    GMM_ASSERT1(mf_u3.get_qdim() == 1 && mf_theta.get_qdim() == 2,
		"wrong qdim for the mesh_fem");
    generic_assembly assem("t1=comp(Grad(#1).vBase(#2));"
			   "M$1(#1,#2)+=t1(:,i,:,i);");
    assem.push_mi(mim);
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mat(const_cast<MAT &>(RM));
    assem.assembly(rg);
  }

  /**@ingroup asm*/
  template<class MAT>
  void asm_coupling_psitheta(const MAT &RM,  const mesh_im &mim,
			     const mesh_fem &mf_u3,
			     const mesh_fem &mf_theta,
			     const mesh_region &rg
			     = mesh_region::all_convexes()) {
    
    GMM_ASSERT1(mf_u3.get_qdim() == 1 && mf_theta.get_qdim() == 2,
		"wrong qdim for the mesh_fem");
    generic_assembly assem("t1=comp(Base(#1).vGrad(#2));"
			   "M$1(#1,#2)+=t1(:,:,2,1)-t1(:,:,1,2);");
    assem.push_mi(mim);
    assem.push_mf(mf_u3);
    assem.push_mf(mf_theta);
    assem.push_mat(const_cast<MAT &>(RM));
    assem.assembly(rg);
  }

  /* ******************************************************************** */
  /*		Mixed linear plate model brick.                           */
  /* ******************************************************************** */

# define MDBRICK_MIXED_LINEAR_PLATE 213456

  /** 
      Mixed linear plate model brick (for thin plates, using Kirchhoff-Love model).

      Do NOT forget to use the mdbrick_plate_closing with this one !

      @ingroup bricks
      @see mdbrick_bilaplacian, mdbrick_mixed_isotropic_linearized_plate
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_mixed_isotropic_linearized_plate
    : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;

    const mesh_im &mim;
    const mesh_fem &mf_ut, &mf_u3, &mf_theta;
    mdbrick_parameter<VECTOR> lambda_, mu_;
    value_type epsilon;
    bool homogeneous, symmetrized, K_uptodate;
    T_MATRIX K;
    size_type nbdof;

    virtual void proper_update(void) {
      K_uptodate = false;
      nbdof = mf_ut.nb_dof() + mf_u3.nb_dof()*3 + mf_theta.nb_dof();
    }

    const T_MATRIX &get_K(void) {
      this->context_check(); 
      if (!K_uptodate || this->parameters_is_any_modified()) {
	gmm::clear(K);
	gmm::resize(K, nbdof, nbdof);
	size_type nd1=mf_ut.nb_dof(), nd2=mf_u3.nb_dof(),nd3=mf_theta.nb_dof();
	gmm::sub_interval I1(0, nd1), I2(nd1, nd2), I3(nd1 + nd2, nd3);
	gmm::sub_interval I4(nd1 + nd2 + nd3, nd2), I5(nd1 + 2*nd2 + nd3, nd2);
	
	asm_stiffness_matrix_for_linear_elasticity
	  (gmm::sub_matrix(K, I1), mim, mf_ut, lambda_.mf(),
	   lambda_.get(), mu_.get());
	gmm::scale(gmm::sub_matrix(K, I1), value_type(2) * epsilon);
	
	
	asm_stiffness_matrix_for_homogeneous_laplacian(gmm::sub_matrix(K, I2),
						       mim, mf_u3);
	gmm::scale(gmm::sub_matrix(K, I2),
		   value_type(2)* epsilon * epsilon * epsilon / value_type(3));
	
	asm_stiffness_matrix_for_linear_elasticity
	  (gmm::sub_matrix(K, I3), mim, mf_theta, lambda_.mf(), lambda_.get(),
	   mu_.get());
	//   gmm::scale(gmm::sub_matrix(K, I3),
	//  	 value_type(2) * epsilon * epsilon * epsilon / value_type(3));
	
	
	asm_coupling_u3theta(gmm::sub_matrix(K, I2, I3), mim, mf_u3, mf_theta);
	gmm::scale(gmm::sub_matrix(K, I2, I3),
		   value_type(2)* epsilon * epsilon * epsilon / value_type(3));
	
	//       cout << "\n\nval p de I2 I3 ";
	//       affiche_moi_valp(gmm::sub_matrix(K, I2, I3));
	
	
	asm_coupling_psitheta(gmm::sub_matrix(K, I4, I3),mim, mf_u3, mf_theta);
	gmm::scale(gmm::sub_matrix(K, I4, I3), epsilon*epsilon/value_type(3));
	
	//       cout << "\n\nval p de I4 I3 ";
	//       affiche_moi_valp(gmm::sub_matrix(K, I4, I3));
	
	
	asm_coupling_psitheta(gmm::transposed(gmm::sub_matrix(K, I3, I4)), mim,
			      mf_u3, mf_theta);
	gmm::scale(gmm::sub_matrix(K, I3, I4), epsilon*epsilon/value_type(3));
	
	
	asm_coupling_u3theta(gmm::transposed(gmm::sub_matrix(K, I3, I5)), mim,
			     mf_u3, mf_theta);
	gmm::scale(gmm::sub_matrix(K, I3, I5), epsilon*epsilon/value_type(3));
	
	if (!symmetrized)
	  asm_stiffness_matrix_for_homogeneous_laplacian
	    (gmm::sub_matrix(K, I5), mim, mf_u3);
	if (symmetrized) {
	  asm_mass_matrix(gmm::sub_matrix(K, I3), mim, mf_theta);
	  asm_coupling_u3theta(gmm::transposed(gmm::sub_matrix(K, I3, I2)),
			       mim, mf_u3, mf_theta);
	  gmm::scale(gmm::sub_matrix(K, I3, I2),
		     value_type(2)*epsilon*epsilon*epsilon / value_type(3));
	  
	  asm_stiffness_matrix_for_homogeneous_laplacian
	    (gmm::sub_matrix(K, I2, I5), mim, mf_u3);
	  gmm::scale(gmm::sub_matrix(K, I2, I5),epsilon*epsilon/value_type(3));
	  asm_stiffness_matrix_for_homogeneous_laplacian
	    (gmm::sub_matrix(K, I5, I2), mim, mf_u3);
	  gmm::scale(gmm::sub_matrix(K, I5, I2),epsilon*epsilon/value_type(3));
	  asm_coupling_u3theta(gmm::sub_matrix(K, I5, I3),mim,mf_u3, mf_theta);
	  gmm::scale(gmm::sub_matrix(K, I5, I3),epsilon*epsilon/value_type(3));
	}
	gmm::scale(gmm::sub_matrix(K, I3),
		   value_type(2)* epsilon * epsilon * epsilon / value_type(3));
	
	this->proper_mixed_variables.clear();
	this->proper_mixed_variables.add(nbdof - mf_u3.nb_dof()*2,
					 mf_u3.nb_dof()*2);
	K_uptodate = true;
	this->parameters_set_uptodate();
      }
      return K;
    }
    
  public :
    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					size_type) {
      gmm::sub_interval SUBI(i0, nbdof);
      gmm::copy(get_K(), gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				size_type) {
      gmm::sub_interval SUBI(i0, nbdof);
      gmm::mult(get_K(), gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residual(), SUBI));
    }

    void set_elastic_coeff(value_type E, value_type nu) {
      lambda_.set(nu * E / (value_type(1) - nu*nu));
      mu_.set(E/(value_type(2)*(value_type(1)+nu)));
    }

    SUBVECTOR get_solution(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), nbdof);
      return gmm::sub_vector(MS.state(), SUBU);
    }
    SUBVECTOR get_ut(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index(), mf_ut.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }
    SUBVECTOR get_u3(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index() + mf_ut.nb_dof(),
			     mf_u3.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }
    SUBVECTOR get_theta(MODEL_STATE &MS) {
      gmm::sub_interval SUBU(this->first_index() + mf_ut.nb_dof()
			     + mf_u3.nb_dof(), mf_theta.nb_dof());
      return gmm::sub_vector(MS.state(), SUBU);
    }

    void init_(void) {
      size_type info = 1 + (symmetrized ? 2 : 0);
      GMM_ASSERT1(mf_ut.get_qdim() == 2, "Qdim of mf_ut should be 2.");
      GMM_ASSERT1(mf_u3.get_qdim() == 1, "Qdim of mf_u3 should be 1.");
      GMM_ASSERT1(mf_theta.get_qdim() == 2, "Qdim of mf_theta should be 2.");
      this->add_proper_mesh_im(mim);
      this->add_proper_mesh_fem(mf_ut, MDBRICK_MIXED_LINEAR_PLATE, info );
      this->add_proper_mesh_fem(mf_u3, MDBRICK_MIXED_LINEAR_PLATE, 0);
      this->add_proper_mesh_fem(mf_theta, MDBRICK_MIXED_LINEAR_PLATE, 0); 
      this->add_proper_mesh_fem(mf_u3, MDBRICK_MIXED_LINEAR_PLATE, 0);
      this->add_proper_mesh_fem(mf_u3, MDBRICK_MIXED_LINEAR_PLATE, 0);
      this->proper_is_symmetric_ = symmetrized;
      this->proper_is_coercive_ = false;
      this->force_update();
    }


    /** constructor for a homogeneous material (constant lambda and
	mu). An inf-sup condition between mf_u3 and mf_theta has to be
	satisfied (mf_u3 should be of lower order than mf_theta).
     
	@param mf_ut the finite element method for the membrane displacement
	@param mf_u3 the finite element method for the transverse displacement.
	@param mf_theta the finite element method for the rotation of the normal (section rotations).
	@param epsilon the thickness of the plate.
    */
    mdbrick_mixed_isotropic_linearized_plate
    (const mesh_im &mim_, const mesh_fem &mf_ut_, const mesh_fem &mf_u3_,
     const mesh_fem &mf_theta_,
     value_type lambdai, value_type mui, double epsilon_, bool sym = false)
      : mim(mim_), mf_ut(mf_ut_), mf_u3(mf_u3_), mf_theta(mf_theta_),
	lambda_("lambda", mf_ut_.linked_mesh(), this),
	mu_("mu", mf_ut_.linked_mesh(), this), epsilon(epsilon_), symmetrized(sym)
    { lambda_.set(lambdai); mu_.set(mui); init_(); }
 
  };

  /**
     Plate source term brick (apply a classical source term on the
     @f$ut@f$, @f$u3@f$, and @f$\theta@f$ fields)
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_source_term : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_source_term<MODEL_STATE> *ut_part, *theta_part, *u3_part,
      *phi_part, *sub_problem;
    
    mdbrick_parameter<VECTOR> B_;
    bool mixed, symmetrized;

    virtual void proper_update(void) {
      size_type n = B_.mf().nb_dof();
      VECTOR Bt(2*n);

      gmm::copy(gmm::sub_vector(B_.get(), gmm::sub_slice(0, n, 3)),
		gmm::sub_vector(Bt, gmm::sub_slice(0, n, 2)));
      gmm::copy(gmm::sub_vector(B_.get(), gmm::sub_slice(1, n, 3)),
		gmm::sub_vector(Bt, gmm::sub_slice(1, n, 2)));
     
      ut_part->source_term().set(B_.mf(), Bt);
      
      VECTOR Bn(n);
      gmm::copy(gmm::sub_vector(B_.get(), gmm::sub_slice(2, n, 3)), Bn);
      if (!mixed || symmetrized)
	u3_part->source_term().set(B_.mf(), Bn);
      
      if (mixed && !symmetrized)
	phi_part->source_term().set(B_.mf(), Bn);
    }

  public :

    /** displacement (ut and u3) source term */
    mdbrick_parameter<VECTOR> &B(void) { return B_; }
    const mdbrick_parameter<VECTOR> &B(void) const { return B_; }
    /** moment source term (i.e. the source term on the rotation of
	the normal). */
    mdbrick_parameter<VECTOR> &M(void) { return theta_part->source_term(); }
    const mdbrick_parameter<VECTOR> &M(void) const
    { return  theta_part->source_term(); }



    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					   size_type) { }
    virtual void do_compute_residual(MODEL_STATE &, size_type, size_type) { }

    mdbrick_plate_source_term(mdbrick_abstract<MODEL_STATE> &problem,
			      const mesh_fem &mf_data, const VECTOR &B__,
			      const VECTOR &M__,
			      size_type bound = size_type(-1),
			      size_type num_fem = 0)
      :  B_("B", mf_data, this, 3) {
      B_.set(B__);
      ut_part = phi_part = u3_part = theta_part = 0;
      mixed = false; symmetrized = false;
      if (problem.get_mesh_fem_info(num_fem).brick_ident
	  == MDBRICK_LINEAR_PLATE)
	{ mixed = false; symmetrized = false; } 
      else if (problem.get_mesh_fem_info(num_fem).brick_ident 
	       == MDBRICK_MIXED_LINEAR_PLATE) {
	mixed=true;
	symmetrized = ((problem.get_mesh_fem_info(num_fem).info) & 2);
      }
      else GMM_ASSERT1(false, "This brick should only be applied to a "
		       "plate problem");
      GMM_ASSERT1((problem.get_mesh_fem_info(num_fem).info & 1)
		  && (num_fem + (mixed ? 4 : 2) < problem.nb_mesh_fems()),
		  "The mesh_fem number is not correct");

      theta_part = new mdbrick_source_term<MODEL_STATE>
	(problem, mf_data, M__, bound, num_fem+2);

      // alias the source_term param of the theta_part brick into this brick
      this->parameters["M"] = &M();

      ut_part = sub_problem = new mdbrick_source_term<MODEL_STATE>
	(*theta_part, mf_data, VECTOR(), bound, num_fem);
    
      if (!mixed || symmetrized)
	sub_problem = u3_part = new mdbrick_source_term<MODEL_STATE>
	  (*ut_part, mf_data, VECTOR(), bound, num_fem+1);
      
      if (mixed && !symmetrized)
	sub_problem = phi_part = new mdbrick_source_term<MODEL_STATE>
	  (*sub_problem, mf_data, VECTOR(), bound, num_fem+4);

      this->add_sub_brick(*sub_problem);

      if (bound != size_type(-1)) {
	this->add_proper_boundary_info(num_fem, bound, MDBRICK_NEUMANN);
	this->add_proper_boundary_info(num_fem+1, bound, MDBRICK_NEUMANN);
      }

      this->force_update();
    }

    ~mdbrick_plate_source_term() {
      delete ut_part;
      if (u3_part) delete u3_part;
      if (phi_part) delete phi_part;
      if (theta_part) delete theta_part;
    }
    
  };


  /**
     Simple support condition for plate model brick (Dirichlet
     condition on @f$ut@f$ and @f$u3@f$, free rotation)

     @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_simple_support : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_Dirichlet<MODEL_STATE> *ut_part, *u3_part;
    mdbrick_Dirichlet<MODEL_STATE> *phi_part, *sub_problem;

    virtual void proper_update(void) {}

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					size_type) { }
    virtual void do_compute_residual(MODEL_STATE &, size_type ,size_type) {}

    mdbrick_plate_simple_support(mdbrick_abstract<MODEL_STATE> &problem,
				 size_type bound,
				 size_type num_fem = 0,
				 constraints_type cot = AUGMENTED_CONSTRAINTS)
      : phi_part(0) {
      ut_part = new mdbrick_Dirichlet<MODEL_STATE>(problem, bound,
						   dummy_mesh_fem(), num_fem);
      ut_part->set_constraints_type(cot);

      u3_part = new mdbrick_Dirichlet<MODEL_STATE>
	(*ut_part, bound, dummy_mesh_fem(), num_fem+1);
      u3_part->set_constraints_type(cot);

      bool mixed = false, symmetrized = false;
      if (problem.get_mesh_fem_info(num_fem).brick_ident
	  == MDBRICK_LINEAR_PLATE)
	{ mixed = false; symmetrized = false; } 
      else if (problem.get_mesh_fem_info(num_fem).brick_ident 
	       == MDBRICK_MIXED_LINEAR_PLATE) {
	mixed=true;
	symmetrized = ((problem.get_mesh_fem_info(num_fem).info) & 2);
      }
      else GMM_ASSERT1(false, "This brick should only be applied to "
		       "a plate problem");
      GMM_ASSERT1((problem.get_mesh_fem_info(num_fem).info & 1)
		  && (num_fem + (mixed ? 4 : 2) < problem.nb_mesh_fems()),
		  "The mesh_fem number is not correct");

      if (mixed) {
	sub_problem = phi_part = new mdbrick_Dirichlet<MODEL_STATE>(*u3_part, bound, dummy_mesh_fem(), num_fem+4);
	sub_problem->set_constraints_type(cot);
      }
      else sub_problem = u3_part;
      this->add_sub_brick(*sub_problem);
      this->add_proper_boundary_info(num_fem, bound, MDBRICK_SIMPLE_SUPPORT);
      this->add_proper_boundary_info(num_fem+1, bound, MDBRICK_SIMPLE_SUPPORT);
      this->add_proper_boundary_info(num_fem+2, bound, MDBRICK_SIMPLE_SUPPORT);
      this->force_update();
    }

    virtual ~mdbrick_plate_simple_support() {
      delete ut_part; delete u3_part;
      if (phi_part) delete phi_part;
    }
    
  };


  /**
    Clamped condition for plate model brick (Dirichlet condition on
    the displacement and the rotation).
    @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_clamped_support: public mdbrick_abstract<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_Dirichlet<MODEL_STATE> ut_part, u3_part, theta_part;
    mdbrick_Dirichlet<MODEL_STATE> *phi_part, *sub_problem;

    virtual void proper_update(void) {}

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &, size_type,
					   size_type) {}
    virtual void do_compute_residual(MODEL_STATE &, size_type, size_type) {}

    mdbrick_plate_clamped_support(mdbrick_abstract<MODEL_STATE> &problem,
				  size_type bound, size_type num_fem = 0,
				  constraints_type cot = AUGMENTED_CONSTRAINTS)
      : ut_part(problem, bound, dummy_mesh_fem(), num_fem),
	u3_part(ut_part, bound, dummy_mesh_fem(), num_fem+1),
	theta_part(u3_part, bound, dummy_mesh_fem(), num_fem+2), phi_part(0) {

      ut_part.set_constraints_type(cot);
      u3_part.set_constraints_type(cot);
      theta_part.set_constraints_type(cot);

      bool mixed = false, symmetrized = false;
      if (problem.get_mesh_fem_info(num_fem).brick_ident
	  == MDBRICK_LINEAR_PLATE)
	{ mixed = false; symmetrized = false; } 
      else if (problem.get_mesh_fem_info(num_fem).brick_ident 
	       == MDBRICK_MIXED_LINEAR_PLATE) {
	mixed=true;
	symmetrized = ((problem.get_mesh_fem_info(num_fem).info) & 2);
      }
      else GMM_ASSERT1(false, "This brick should only be applied to "
		       "a plate problem");
      GMM_ASSERT1((problem.get_mesh_fem_info(num_fem).info & 1)
		  && (num_fem + (mixed ? 4 : 2) < problem.nb_mesh_fems()),
		  "The mesh_fem number is not correct");

      if (mixed) {
	sub_problem = phi_part = new  mdbrick_Dirichlet<MODEL_STATE>
	  (theta_part, bound, dummy_mesh_fem(), num_fem+4);
	sub_problem->set_constraints_type(cot);
	this->add_sub_brick(*phi_part);
      }
      else { 
	this->add_sub_brick(theta_part);
	sub_problem = &theta_part;
      }
      this->add_proper_boundary_info(num_fem, bound, MDBRICK_CLAMPED_SUPPORT);
      this->add_proper_boundary_info(num_fem+1, bound, MDBRICK_CLAMPED_SUPPORT);
      this->add_proper_boundary_info(num_fem+2, bound, MDBRICK_CLAMPED_SUPPORT);

      this->force_update();
    }

    ~mdbrick_plate_clamped_support() { if (phi_part) delete phi_part; }
    
  };


  /* ******************************************************************** */
  /*		Free edges condition for mixed plate model brick.         */
  /* ******************************************************************** */

  /**@ingroup asm*/
  template<class VEC>
  void asm_constraint_on_theta(const VEC &V, const mesh_im &mim, 
			       const mesh_fem &mf_theta, const mesh_region &boundary) {
    generic_assembly assem("t=comp(vBase(#1).Normal());"
			   "V(#1)+= t(:,2,1) - t(:,1,2);");
    assem.push_mi(mim);
    assem.push_mf(mf_theta);
    assem.push_vec(const_cast<VEC &>(V));
    assem.assembly(boundary);
  }

  /**
     Free edges condition for mixed plate model brick. This brick has to be added
     for the mixed linearized plate brick after all other boundary conditions.

     (the reason is that the brick has to inspect all other boundary
     conditions to determine the number of disconnected boundary parts
     which are free edges).

     @ingroup bricks
   */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_plate_closing: public mdbrick_abstract<MODEL_STATE> {
 
    TYPEDEF_MODEL_STATE_TYPES;
    
    mdbrick_abstract<MODEL_STATE> *sub_problem;
    const mesh_fem *mf_theta;
    gmm::row_matrix<gmm::rsvector<value_type> > CO;
    size_type num_fem;
    bool mixed, symmetrized, allclamped, with_multipliers;

    virtual void proper_update(void) {
      mf_theta = this->mesh_fems[num_fem+2];
      if (!mixed) {  allclamped = false; gmm::resize(CO, 0, 0); return; }
      allclamped = true;
      std::vector<size_type> cv_nums;
      std::vector<short_type> face_nums;
      
      const mesh *msh = &(mf_theta->linked_mesh());
    
      getfem::mesh_region border_faces;
      getfem::outer_faces_of_mesh(*msh, border_faces);
      dal::bit_vector vb = msh->regions_index();
      
      for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
	bool add = true;
	// cout << "face " << it->f << " of cv " << it->cv << "boundaries : ";
	for (dal::bv_visitor i(vb); !i.finished(); ++i) {
	  if (msh->region(i).is_in(it.cv(),it.f())) {
	    // cout << i << endl;
	    bound_cond_type bct = this->boundary_type(num_fem, i);
	    if (bct != MDBRICK_UNDEFINED && bct != MDBRICK_NEUMANN) add = false;
	    if (bct != MDBRICK_CLAMPED_SUPPORT) allclamped = false;
	  }
	}
	
	if (add) {
	  cv_nums.push_back(it.cv()); 
	  face_nums.push_back(it.f()); 
	  allclamped = false;
	  // cout << " adding";
	}
	// cout << endl;
      }
      //cout << "allclamped = " << allclamped << endl;

      std::vector<size_type> comp_conns(cv_nums.size(), size_type(-1)); 
      size_type nbmax = msh->points().ind_last() + 1, p1, p2;
      std::vector<size_type> E1(nbmax, size_type(-1)), E2(nbmax, size_type(-1));
      for (size_type j = 0; j < cv_nums.size(); ++j) {
	p1 = msh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[0];
	p2 = msh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[1];
	if (E1[p1] == size_type(-1)) E1[p1] = j; else E2[p1] = j;
	if (E1[p2] == size_type(-1)) E1[p2] = j; else E2[p2] = j;	
      }

      size_type comp_conn = 0;
      for (size_type i = 0; i < comp_conns.size(); ++i) {
	if (comp_conns[i] == size_type(-1)) {
	  
	  comp_conns[i] = comp_conn;
	  p1 = msh->ind_points_of_face_of_convex(cv_nums[i],face_nums[i])[0];
	  p2 = msh->ind_points_of_face_of_convex(cv_nums[i],face_nums[i])[1];
	  size_type j1 = (E1[p1] == i) ? E2[p1] :  E1[p1];
	  size_type j2 = (E1[p2] == i) ? E2[p2] :  E1[p2];
	  
	  for (unsigned k = 0; k < 2; ++k) {
	    size_type j = (k == 0) ? j1 : j2;
	    
	    while (j != size_type(-1) && comp_conns[j] == size_type(-1)) {
	      comp_conns[j] = comp_conn;
	      p1 = msh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[0];
	      p2 = msh->ind_points_of_face_of_convex(cv_nums[j],face_nums[j])[1];
	      size_type i1 = (E1[p1] == j) ? E2[p1] :  E1[p1];
	      size_type i2 = (E1[p2] == j) ? E2[p2] :  E1[p2];
	      if (i1 == size_type(-1) || comp_conns[i1] != size_type(-1))
		j = i2; else j = i1;
	    }
	  }
	  
	  ++comp_conn;
	}
      }

      //cout << "Number of comp conn : " << comp_conn << endl;
      if (comp_conn < 2) {
	gmm::resize(CO, 0, 0);
      }
      else {
	mesh_region boundary;
	gmm::resize(CO, comp_conn, mf_theta->nb_dof());
	for (size_type k = 0; k < comp_conn; ++k) {
	  for (size_type i = 0; i < comp_conns.size(); ++i)
	    if (comp_conns[i] == k) {
	      boundary.add(cv_nums[i], face_nums[i]);
	    }
	  std::vector<value_type> V(mf_theta->nb_dof());
	  asm_constraint_on_theta(V, *(this->mesh_ims[0]), *mf_theta,
				  boundary);
	  gmm::copy(V, gmm::mat_row(CO, k));
	}
	// cout << "CO = " << CO << endl;
      }

      size_type nb_const = gmm::mat_nrows(CO) + (allclamped ? 1:0);
      this->proper_mixed_variables.clear();
      this->proper_additional_dof = with_multipliers ? nb_const : 0;
      this->proper_nb_constraints = with_multipliers ? 0 : nb_const;
      if (with_multipliers)
	this->proper_mixed_variables.add(sub_problem->nb_dof(), nb_const);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) {
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
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				   size_type j0) {
      gmm::sub_interval SUBJ(i0+this->mesh_fem_positions[num_fem+2],
			     mf_theta->nb_dof());
      if (with_multipliers) {
	size_type nbd = sub_problem->nb_dof();
	if (gmm::mat_nrows(CO) > 0) {
	  gmm::sub_interval SUBI(i0 + nbd, gmm::mat_nrows(CO));
	  gmm::mult(CO, gmm::sub_vector(MS.state(), SUBJ),
		    gmm::sub_vector(MS.residual(), SUBI));
	  gmm::mult_add(gmm::transposed(CO), gmm::sub_vector(MS.state(), SUBI),
			gmm::sub_vector(MS.residual(), SUBJ));
	}
	if (allclamped) {
	  size_type i = i0 + nbd + gmm::mat_nrows(CO);
	  size_type j = i0 + this->mesh_fem_positions[num_fem+3];
	  MS.residual()[i] = MS.state()[j];
	  MS.residual()[j] += MS.state()[i];
	}
      }
      else {
	size_type ncs = sub_problem->nb_constraints();
	if (gmm::mat_nrows(CO) > 0) {
	  gmm::sub_interval SUBI(j0+ncs,gmm::mat_nrows(CO));
	  gmm::mult(CO, gmm::sub_vector(MS.state(), SUBJ),
		    gmm::sub_vector(MS.constraints_rhs(), SUBI));
	}
	if (allclamped) {
	  (MS.constraints_rhs())[j0+ncs+gmm::mat_nrows(CO)] =
	    (MS.state())[i0 + this->mesh_fem_positions[num_fem+3]];
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
      else GMM_ASSERT1(false, "This brick should only be applied to "
		       "a plate problem");
      GMM_ASSERT1((problem.get_mesh_fem_info(num_fem).info & 1)
		  && (num_fem + (mixed ? 4 : 2) < problem.nb_mesh_fems()),
		  "The mesh_fem number is not correct");

      this->add_sub_brick(problem);
      this->force_update();
    }
    
  };





}  /* end of namespace getfem.                                             */


#endif /* GETFEM_LINEARIZED_PLATES_H__ */
