// -*- c++ -*- (enables emacs c++ mode)
#ifndef NAVIER_STOKES_H_
#define NAVIER_STOKES_H_

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <gmm.h>

namespace getfem {

  /*
   * div(uuT) term. 
   */
  
  template<typename MAT, typename VEC>
  void asm_NS_uuT(const MAT &M, const mesh_im &mim,
		  const mesh_fem &mf, const VEC &U0,
		  const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;    
    assem.set("u=data(#1); "
	      "t = comp(vBase(#1).vBase(#1).vGrad(#1));"
	      "M(#1, #1)+=u(i).t(i,k,:,j,:,j,k)"
	      // ";M(#1, #1)+=u(i).t(i,k,:,k,:,j,j)" // optional (or *0.5)
	      );
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U0);
    assem.push_mat(const_cast<MAT&>(M));
    assem.assembly(rg);
  }


  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_NS_uuT : public mdbrick_abstract<MODEL_STATE> {

    TYPEDEF_MODEL_STATE_TYPES;
    
    mdbrick_abstract<MODEL_STATE> &sub_problem;
    VECTOR U0;
    size_type num_fem, i1, nbd;
    T_MATRIX K;

    virtual void proper_update(void) {}

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					size_type) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      gmm::sub_interval SUBI(i0+this->mesh_fem_positions[num_fem],
			     mf_u.nb_dof());
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      gmm::clear(K);
      asm_NS_uuT(K, *(this->mesh_ims[0]), mf_u, U0);
      gmm::add(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				size_type) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      gmm::sub_interval SUBI(i0+this->mesh_fem_positions[num_fem],
			     mf_u.nb_dof());
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residu(), SUBI);
      gmm::mult_add(K, gmm::sub_vector(MS.state(), SUBI), SUBV);
    }

    template <class VEC> void set_U0(const VEC &q) {
      gmm::resize(U0, gmm::vect_size(q));
      gmm::copy(q, U0);
    }

    void init_(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->proper_is_symmetric_ = false;
      this->force_update();
    }

    // Constructor which homogeneous diagonal Q
    mdbrick_NS_uuT(mdbrick_abstract<MODEL_STATE> &problem,
		   size_type num_fem_=0) 
      : sub_problem(problem), num_fem(num_fem_) { init_(); }

  };


  /* ******************************************************************** */
  /*		A non-reflective condition.                               */
  /* ******************************************************************** */


# define MDBRICK_NS_NONREF1 3212435


  template<typename MAT, typename VECT>
  void asm_mass_matrix_param_un(MAT &M, const mesh_im &mim,
				const mesh_fem &mf_u,
				const mesh_fem &mfdata,
				const VECT &F, scalar_type nudt,		       
				const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;

    std::stringstream s;
    s << "F=data(#2);"
      "M(#1,#1)+=sym(comp(vBase(#1).vBase(#1).vBase(#2).Normal())(:,i,:,i,j,k,k).F(j))*("
      << nudt << ");";

    assem.set(s.str());
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mfdata);
    assem.push_data(F);
    assem.push_mat(M);
    assem.assembly(rg);
  }

  /*
    du/dt + u.n * du/dn = 0
    
    handled as:
      du/dt - u.n * lambda = 0
                    lambda = -du/dn
  */
  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_NS_nonref1 : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_lambda; /* lambda = -du/dn on the boundary */
    T_MATRIX B, M;
    VECTOR Un;
    size_type num_fem, i1, nbd, boundary;
    value_type dt, nu;
    dal::bit_vector doflambda;
    std::vector<size_type> indices;
    gmm::sub_index SUBK;

    void proper_update(void) {
      mesh_fem &mf_u = *(this->mesh_fems.at(num_fem));
      i1 = this->mesh_fem_positions.at(num_fem);
      nbd = mf_u.nb_dof();
      size_type nd = mf_u.nb_dof(), ndd = mf_lambda.nb_dof();
      T_MATRIX Baux(nd, ndd);
    
      asm_mass_matrix(Baux, *(this->mesh_ims.at(0)), mf_lambda, mf_u, boundary);
      doflambda = mf_lambda.dof_on_set(boundary);
      for (dal::bv_visitor i(doflambda); !i.finished(); ++i)
	indices.push_back(i);
      gmm::resize(M, doflambda.card(), doflambda.card());
      SUBK = gmm::sub_index(indices);
      gmm::resize(B, doflambda.card(), ndd);
      gmm::copy(gmm::sub_matrix(Baux, SUBK, gmm::sub_interval(0, ndd)), B);
      this->proper_mixed_variables.clear();
      this->proper_mixed_variables.add(sub_problem.nb_dof(), doflambda.card());
      this->proper_additional_dof = doflambda.card();
    }

  public :
    
    const T_MATRIX &get_B(void) const { this->context_check(); return B; }
    
    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type) {
      gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), doflambda.card());
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::copy(B, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
      gmm::copy(gmm::transposed(B),
		gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
      gmm::copy(M, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type) {
      gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), doflambda.card());
      gmm::sub_interval SUBJ(i0+i1, nbd);
      gmm::mult(B, gmm::sub_vector(MS.state(), SUBJ),
		gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(gmm::transposed(B), gmm::sub_vector(MS.state(), SUBI),
		    gmm::sub_vector(MS.residu(), SUBJ));
      gmm::mult_add(M, gmm::sub_vector(MS.state(), SUBI),
		    gmm::sub_vector(MS.residu(), SUBI));
      gmm::mult_add(B, gmm::scaled(Un, -value_type(1)),
		    gmm::sub_vector(MS.residu(), SUBI));
    }
    
    void set_Un(const VECTOR &B__) {
      mesh_fem &mf_u = *(this->mesh_fems.at(num_fem));
      this->context_check();
      gmm::resize(Un, gmm::vect_size(B__));
      gmm::copy(B__, Un);
      T_MATRIX Maux(mf_lambda.nb_dof(), mf_lambda.nb_dof());
      asm_mass_matrix_param_un(Maux, *(this->mesh_ims.at(0)), mf_lambda, mf_u,
 			       Un, -dt/nu, boundary);
      // asm_mass_matrix(Maux, *(this->mesh_ims.at(0)), mf_lambda, mf_u, boundary);
      gmm::copy(gmm::sub_matrix(Maux, SUBK, SUBK), M);
      gmm::scale(M, dt);
    }

    void set_dt(value_type dt_) {
      this->context_check();
      gmm::scale(M, dt / dt_); dt = dt_; 
    }

    void init_(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = false;
      this->force_update();
    }

    mdbrick_NS_nonref1(mdbrick_abstract<MODEL_STATE> &problem,
		       mesh_fem &mf_lambda_, size_type bound, value_type dt_,
		       value_type nu_, size_type num_fem_=0)
      : sub_problem(problem), mf_lambda(mf_lambda_), num_fem(num_fem_), boundary(bound), dt(dt_), nu(nu_)
    { init_(); }

  };






}

#endif
