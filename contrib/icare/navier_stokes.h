// -*- c++ -*- (enables emacs c++ mode)
#ifndef NAVIER_STOKES_H_
#define NAVIER_STOKES_H_

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <gmm.h>

namespace getfem {

  /**
   * div(uuT) term. 
   */
  
  template<typename MAT, typename VEC>
  void asm_NS_uuT(const MAT &M, const mesh_im &mim,
			 const mesh_fem &mf, const VEC &U0) {
    generic_assembly assem;    
    assem.set("u=data(#1); "
	      "t = comp(vGrad(#1).vBase(#1).vBase(#1));"
	      "M(#1, #1)+=u(i).t(:,j,k,i,j,:,k);"
	      "M(#1, #1)+=u(i).t(i,j,j,:,k,:,k)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U0);
    assem.push_mat(const_cast<MAT&>(M));
    assem.volumic_assembly();
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
      this->update_from_context();
    }

    // Constructor which homogeneous diagonal Q
    mdbrick_NS_uuT(mdbrick_abstract<MODEL_STATE> &problem,
		   size_type num_fem_=0) 
      : sub_problem(problem), num_fem(num_fem_) { init_(); }

  };
}

#endif
