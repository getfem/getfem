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
		  const mesh_fem &mf, const VEC &U0,
		  const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;    
    assem.set("u=data(#1); "
	      "t = comp(vGrad(#1).vBase(#1).vBase(#1));"
	      "M(#1, #1)+=u(i).t(:,j,k,i,j,:,k);"
	      "M(#1, #1)+=u(i).t(i,j,j,:,k,:,k)");
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





  /* ******************************************************************** */
  /*		A non-reflective condition.                               */
  /* ******************************************************************** */


  template<typename MAT, typename VEC>
  void asm_NS_uuT(const MAT &M, const mesh_im &mim,
		  const VEC &U0, const mesh_fem &mf,
		  const mesh_region &rg = mesh_region::all_convexes()) {
    generic_assembly assem;    
    assem.set("u=data(#1); "
	      "t = comp(vGrad(#1).vBase(#1).vBase(#1).Normal().Normal());"
	      "M(#1, #1)+=u(i).t(i,k,j,:,l,:,k,l,j);");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_data(U0);
    assem.push_mat(const_cast<MAT&>(M));
    assem.assembly(rg);
  }

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_NS_nonref1 : public mdbrick_abstract<MODEL_STATE>  {
    
    TYPEDEF_MODEL_STATE_TYPES;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    VECTOR H_, Vn;
    C_MATRIX G;
    T_MATRIX M;
    VECTOR CRHS;
    size_type boundary, nb_const, num_fem;
    bool with_H, with_multipliers;
    gmm::sub_index SUB_CT;
    size_type i1, nbd;
    value_type dt;

    void compute_constraints(void) {
      size_type nb_const_old = nb_const;
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      i1 = this->mesh_fem_positions[num_fem];
      nbd = mf_u.nb_dof();
      size_type nd = mf_u.nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > A(nd, nd);
      VECTOR V(nd);

      asm_NS_nonref1(A, Vn, *(this->mesh_ims[0]), mf_u, boundary);
      gmm::add(M, A);

      gmm::mult(M, Vn, V);
      
      R tol=gmm::mat_maxnorm(A)*gmm::default_tol(value_type())*R(100);
      gmm::clean(A, tol);
      std::vector<size_type> ind_ct;
      dal::bit_vector nn = mf_u.dof_on_set(boundary);
      // The following filter is not sufficient for an arbitrary matrix field
      // H for the multipliers version. To be ameliorated.
      for (size_type i = nn.take_first(); i != size_type(-1); i << nn)
	if (!with_multipliers || gmm::vect_norm2(gmm::mat_row(A, i)) > tol)
	  ind_ct.push_back(i);
      nb_const = ind_ct.size();
      SUB_CT = gmm::sub_index(ind_ct);
      gmm::resize(G, nb_const, nd);
      gmm::copy(gmm::sub_matrix(A, SUB_CT, gmm::sub_interval(0, nd)), G);

      gmm::resize(CRHS, nb_const);
      gmm::copy(gmm::sub_vector(V, SUB_CT), CRHS);
      if (nb_const != nb_const_old && nb_const != size_type(-1))
	this->force_update();
    }

    virtual void proper_update(void) {
      mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      size_type nd = mf_u.nb_dof();
      gmm::resize(M, nd, nd);
      gmm::resize(Vn, nd);
      asm_mass_matrix(M, *(this->mesh_ims[0]), mf_u, boundary);
      gmm::scale(M, 1./dt);
      nb_const = size_type(-1);
      compute_constraints();
      this->proper_mixed_variables.clear();
      this->proper_additional_dof = with_multipliers ? nb_const : 0;
      this->proper_nb_constraints = with_multipliers ? 0 : nb_const;
      if (with_multipliers)
	this->proper_mixed_variables.add(sub_problem.nb_dof(), nb_const);
    }

  public :

    virtual void do_compute_tangent_matrix(MODEL_STATE &MS, size_type i0,
					   size_type j0) {
      if (with_multipliers) {
	gmm::sub_interval SUBI(i0+sub_problem.nb_dof(), nb_const);
	gmm::sub_interval SUBJ(i0+i1, nbd);
	gmm::copy(G, gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBJ));
	gmm::copy(gmm::transposed(G),
		  gmm::sub_matrix(MS.tangent_matrix(), SUBJ, SUBI));
	gmm::clear(gmm::sub_matrix(MS.tangent_matrix(), SUBI, SUBI));
      }
      else {	  
	size_type ncs = sub_problem.nb_constraints();
	gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0+i1, nbd);
	gmm::copy(G, gmm::sub_matrix(MS.constraints_matrix(), SUBI, SUBJ));
      }
    }
    virtual void do_compute_residu(MODEL_STATE &MS, size_type i0,
				   size_type j0) {
      if (with_multipliers) {
	gmm::sub_interval SUBI(i0 + sub_problem.nb_dof(), nb_const);
	gmm::sub_interval SUBJ(i0+i1, nbd);
	
	gmm::mult(G, gmm::sub_vector(MS.state(), SUBJ),
		  gmm::scaled(CRHS, value_type(-1)),
		  gmm::sub_vector(MS.residu(), SUBI));
	
	gmm::mult_add(gmm::transposed(G), gmm::sub_vector(MS.state(), SUBI),
		      gmm::sub_vector(MS.residu(), SUBJ));
      }
      else {
	size_type ncs = sub_problem.nb_constraints();
	gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0+i1, nbd);
	gmm::mult(G, gmm::scaled(gmm::sub_vector(MS.state(), SUBJ), 
				 value_type(-1)),
		  CRHS, gmm::sub_vector(MS.constraints_rhs(), SUBI));
      }
    }

    void set_Vn(const VECTOR &B__) {
      this->context_check();
      gmm::resize(Vn, gmm::vect_size(B__));
      gmm::copy(B__, Vn);
      compute_constraints();
    }

    void init_(void) {
      this->add_sub_brick(sub_problem);
      this->proper_is_coercive_ = !with_multipliers;
      this->add_proper_boundary_info(num_fem, boundary,
				     MDBRICK_NAVIERSTOKESNONREF1);
      this->update_from_context();
    }

    // Constructor which does not define the rhs (0 rhs in fact)
    mdbrick_NS_nonref1(mdbrick_abstract<MODEL_STATE> &problem,
		       size_type bound, value_type dt_,
		       size_type num_fem_=0, bool with_mult = false)
      : sub_problem(problem), dt(dt_), boundary(bound),
	num_fem(num_fem_), with_multipliers(with_mult) {
      init_();
    }
  };















}

#endif
