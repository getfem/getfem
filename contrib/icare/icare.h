// -*- c++ -*- (enables emacs c++ mode)
#ifndef NAVIER_STOKES_H_
#define NAVIER_STOKES_H_

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_modeling.h"
#include "gmm/gmm.h"

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

 template<typename VEC, typename VECTOR>
 void traineePortance2D(VEC &Cd, VEC &Cl, const mesh_im &mim,
			const mesh_fem &mf, const VECTOR &lambda,
			const mesh_region &rg) {

   generic_assembly assem;
   assem.set("l=data$1(#1); V$1()+=comp(vBase(#1).Normal())(i,1,1).l(i);V$2()+=comp(vBase(#1).Normal())(i,1,2).l(i);");
   
   assem.push_mi(mim);
   
   assem.push_mf(mf);
   assem.push_data(lambda);

   assem.push_vec(Cd);
   assem.push_vec(Cl);

   assem.assembly(rg); //region = bord du cylindre
						 
 }




 template<typename VEC, typename VECTOR>
 void traineePortance3D(VEC &Cxn, VEC &Cxp,VEC &Cyn, VEC &Cyp, const mesh_im &mim,
			const mesh_fem &mf1, const mesh_fem &mf2,
			const VECTOR &nuDxU, const VECTOR &nuDyU, const VECTOR &nuDzU,
			const VECTOR &nuDxV, const VECTOR &nuDyV, const VECTOR &nuDzV, 
			const VECTOR &nuDxW, const VECTOR &nuDyW, const VECTOR &nuDzW, 
			const VECTOR &PP, const mesh_region &rg) {

   generic_assembly assem;

   assem.set("nuDxU=data$1(#1); nuDyU=data$2(#1); nuDzU=data$3(#1); nuDxV=data$4(#1); nuDyU=data$5(#1); nuDzU=data$6(#1); nuDxW=data$7(#1); nuDyW=data$8(#1);  nuDzW=data$9(#1);    p=data$10(#2);"
	     "t1=comp(vBase(#1).Normal());"
             "V$1()+=2*t1(i,1,1).nuDxU(i) + t1(i,2,2).nuDxU(i) + t1(i,1,2).nuDxV(i) + t1(i,1,3).nuDxW(i) + t1(i,3,3).nuDxU(i);"
	     "t2=comp(vBase(#2).Normal());"
             "V$2()+=t2(i,1,1).p(i);"
	     "t3=comp(vBase(#1).Normal());"
             "V$3()+=2*t3(i,2,2).nuDxV(i) + t3(i,2,1).nuDxU(i) + t3(i,1,1).nuDxV(i) + t3(i,2,3).nuDxW(i) + t3(i,3,3).nuDxV(i);"
	     "t4=comp(vBase(#2).Normal());"
             "V$4()+=t4(i,1,2).p(i);");




   assem.set("nuDxU=data$1(#1); nuDyU=data$2(#1); nuDxV=data$3(#1); nuDyV=data$4(#1); p=data$5(#2);"
	     "t1=comp(vBase(#1).Normal());"
             "V$1()+=2*t1(i,1,1).nuDxU(i) + t1(i,1,2).nuDyU(i) + t1(i,2,2).nuDxV(i);"
	     "t2=comp(vBase(#2).Normal());"
             "V$2()+=-1*t2(i,1,1).p(i);"
	     "t3=comp(vBase(#1).Normal());"
             "V$3()+=2*t3(i,2,2).nuDyV(i) + t3(i,1,1).nuDyU(i) + t3(i,2,1).nuDxV(i);"
	     "t4=comp(vBase(#2).Normal());"
             "V$4()+=-1*t4(i,1,2).p(i);");


   assem.push_mi(mim);

   assem.push_mf(mf1);
   assem.push_mf(mf2);

   assem.push_data(nuDxU);   
   assem.push_data(nuDyU);
   assem.push_data(nuDzU);

   assem.push_data(nuDxV);   
   assem.push_data(nuDyV);
   assem.push_data(nuDzV);

   assem.push_data(nuDxW);   
   assem.push_data(nuDyW);
   assem.push_data(nuDzW);
   
   assem.push_data(PP);

   assem.push_vec(Cxn); 
   assem.push_vec(Cxp);
   assem.push_vec(Cyn); 
   assem.push_vec(Cyp);
  
   assem.assembly(rg); //region = bord du cylindre
						 
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
      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      gmm::sub_interval SUBI(i0+this->mesh_fem_positions[num_fem],
			     mf_u.nb_dof());
      gmm::resize(K, mf_u.nb_dof(), mf_u.nb_dof());
      gmm::clear(K);
      asm_NS_uuT(K, *(this->mesh_ims[0]), mf_u, U0, 
		 mf_u.linked_mesh().get_mpi_region());
      gmm::add(K, gmm::sub_matrix(MS.tangent_matrix(), SUBI));
    }
    virtual void do_compute_residual(MODEL_STATE &MS, size_type i0,
				size_type) {
      const mesh_fem &mf_u = *(this->mesh_fems[num_fem]);
      gmm::sub_interval SUBI(i0+this->mesh_fem_positions[num_fem],
			     mf_u.nb_dof());
      typename gmm::sub_vector_type<VECTOR *, gmm::sub_interval>::vector_type
	SUBV = gmm::sub_vector(MS.residual(), SUBI);
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
  // idée: construire la condition non-réflective comme une condition de dirichlet en explicitant tous les termes: Unp1= Un+dt*Un.N*lamda




// initialisation des CL
// résolution du pb -> nouveau lambda<--| 
//mise à jour de CL de sortie        -->|  



# define MDBRICK_NS_NONREF1 3212435

 template<typename VECTOR, typename VECT>
  void asm_nonref_right_hand_side(VECTOR &R,
				  const mesh_im &mim,
				  const mesh_fem &mf_u,
				  const mesh_fem &mf_mult,
				  const VECT &Un,
				  const mesh_region &rg) {
    generic_assembly assem;
    // construction du terme de droite dans [M]*Unp1=F
    
    // mise en place de Un.N + Un.N*(dUn/dn).N

    std::stringstream ss;
    ss << "u=data$1(#1); "
     "V(#1)+=comp(vBase(#1).vBase(#2))(j,i,:,i).u(j)";
 
    assem.set(ss.str());
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);
    assem.push_data(Un);
    assem.push_vec(R);
    assem.assembly(rg);
  }

    

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_NS_nonref1 : public mdbrick_constraint<MODEL_STATE> {
    
    TYPEDEF_MODEL_STATE_TYPES;

    size_type boundary;
    bool mfdata_set, B_to_be_computed, Un_is_modified;
    gmm::sub_index SUB_CT;
    const mesh_fem *mf_mult;
    VECTOR Un;
    scalar_type dt;
    
    const mesh_fem &mf_u() { return *(this->mesh_fems[this->num_fem]); }
    const mesh_im  &mim() { return *(this->mesh_ims[0]); }

    void compute_constraints(unsigned version) {
      size_type ndu = mf_u().nb_dof(), ndm = mf_mult->nb_dof();
      gmm::row_matrix<gmm::rsvector<value_type> > M(ndm, ndu);
      VECTOR V(ndm);
      GMM_TRACE2("Assembling Dirichlet constraints, version " << version);
      if ((version | ASMDIR_BUILDH)) {
	asm_mass_matrix(M, mim(), mf_u(), *mf_mult,
			mf_u().linked_mesh().get_mpi_region());
      }
      if ((version | ASMDIR_BUILDR)) {
	asm_nonref_right_hand_side(V, mim(), mf_u(), *mf_mult, Un, mf_u().linked_mesh().get_mpi_sub_region(boundary));
      }
      if (version & ASMDIR_BUILDH)
	gmm::copy(gmm::sub_matrix(M, SUB_CT, gmm::sub_interval(0, ndu)), 
		  this->B);
      gmm::copy(gmm::sub_vector(V, SUB_CT), this->CRHS);
    }

    virtual void recompute_B_sizes(void) {
      size_type nd = mf_u().nb_dof();
      dal::bit_vector dof_on_bound;
      if (mf_mult->is_reduced())
	dof_on_bound.add(0, nd);
      else
	dof_on_bound = mf_mult->basic_dof_on_region(boundary);
      size_type nb_const = dof_on_bound.card();
      std::vector<size_type> ind_ct;
      for (dal::bv_visitor i(dof_on_bound); !i.finished(); ++i)
	ind_ct.push_back(i);
      SUB_CT = gmm::sub_index(ind_ct);
      gmm::resize(this->B, nb_const, nd);
      gmm::resize(this->CRHS, nb_const);
      B_to_be_computed = true;
    }

    virtual void recompute_B(void) {
      unsigned version = 0;
      if (Un_is_modified) { version = ASMDIR_BUILDR; }
      if (B_to_be_computed) { version = ASMDIR_BUILDR | ASMDIR_BUILDH; }
      if (version) { 
	compute_constraints(version);
	this->parameters_set_uptodate();
	B_to_be_computed = false; Un_is_modified = false;
      }
    }

  public :

    void set_Un(const VECTOR &U) {
      Un = U;
      Un_is_modified = true;
    }

    mdbrick_NS_nonref1(mdbrick_abstract<MODEL_STATE> &problem,
		       size_type bound, scalar_type dt_,
		       const mesh_fem &mf_mult_ = dummy_mesh_fem(), 
		       size_type num_fem_=0)
      : mdbrick_constraint<MODEL_STATE>(problem, num_fem_), 
	boundary(bound) {
      dt = dt_;
      mf_mult = (&mf_mult_ == &dummy_mesh_fem()) ? &(mf_u()) : &mf_mult_;
      GMM_ASSERT1(mf_mult->get_qdim() == mf_u().get_qdim(), "The lagrange multipliers mesh fem "
		  "for the mdbrick_NS_nonref1 brick should have the same Qdim as "
		  "the main mesh_fem");

      this->add_proper_boundary_info(this->num_fem, boundary, 
				     MDBRICK_DIRICHLET);
      this->add_dependency(*mf_mult);
      mfdata_set = false; B_to_be_computed = true; Un_is_modified = false;
      this->force_update();
    }
  };



  /* non reflective boundary conditions */


  // construction du terme de droite dans [M]*Unp1=F
  // mise en place de Un + Un.N*(dUn/dn).N
  template <typename VEC1, typename VEC2>
  void asm_basic_non_reflective_bc(VEC1 &VV, 
				   const getfem::mesh_im &mim, 
				   const getfem::mesh_fem &mf_u, 
				   const VEC2 &Un0, 
				   const getfem::mesh_fem &mf_mult,
				   scalar_type dt, scalar_type nu,
				   const getfem::mesh_region &rg) {
    getfem::generic_assembly assem;
    std::stringstream ss;
    ss << "u=data$1(#1); "
      "V(#2)+="
       << -dt << "*"
      " comp(vBase(#1).Normal().vGrad(#1).Normal().vBase(#2))(l,i,i,m,k,j,j,:,k).u(l).u(m)+"     //"(l,i,i,m,j,k,j,:,k).u(l).u(m)+"
      "comp(vBase(#1).vBase(#2))(j,i,:,i).u(j)"                                                  //"comp(vBase(#1).vBase(#2).u(j))(j,i,:,i)";
       <<-nu<<"*" << dt << "*"
      "comp(vGrad(#1).vGrad(#2))(l,i,j,:,i,j).u(l)"
       <<-nu<<"*" << dt << "*"
      "comp(vGrad(#1).Normal().vGrad(#2).Normal())(l,i,j,j,:,i,k,k).u(l)";
    assem.set(ss.str());
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);  
    assem.push_data(Un0);
    assem.push_vec(VV);
    assem.assembly(rg);
  }


  template<typename VECT1> class improved_non_reflective_bc_nonlinear_term 
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf;
    std::vector<scalar_type> U;
    scalar_type dt, nu;
    size_type N;
    base_vector coeff, valU;
    base_matrix gradU, hessU;
    bgeot::multi_index sizes_;
  public:
    improved_non_reflective_bc_nonlinear_term
    (const mesh_fem &mf_, const VECT1 &U_, scalar_type dt_, scalar_type nu_)
      : mf(mf_), U(mf_.nb_basic_dof()), dt(dt_), nu(nu_), N(mf_.get_qdim()), 
	valU(N), gradU(N, N), hessU(N,N*N) {
      sizes_.resize(1); sizes_[0] = short_type(N); /*assert(N == 2);*/
      mf.extend_vector(U_, U);
    }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_basic_dof_of_element(cv))),
		coeff);
      ctx.pf()->interpolation(ctx, coeff, valU, mf.get_qdim());
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf.get_qdim());
      ctx.pf()->interpolation_hess(ctx, coeff, hessU, mf.get_qdim());

      if (N==2){
	t[0] = (valU[0] +  nu *dt* hessU(0,3) ) / (1 + gradU(0,0)*dt);
	t[1] = valU[1] +   nu *dt* hessU(1,3) - dt * t[0] * gradU(1,0);
      }
      if (N==3) {
	t[0] = (valU[0] + nu *dt* hessU(0,4) ) / (1 + gradU(0,0)*dt);
	t[1] =  valU[1] + nu *dt* hessU(1,4) - dt * t[0] * gradU(1,0);
	t[2] =  valU[2] + nu *dt* hessU(2,4)  - dt * t[0] * gradU(2,0);
      }
    }
  };


  template <typename VEC1, typename VEC2>
  void asm_improved_non_reflective_bc(VEC1 &VV, 
				      const getfem::mesh_im &mim, 
				      const getfem::mesh_fem &mf_u, 
				      const VEC2 &Un0, 
				      const getfem::mesh_fem &mf_mult,
				      scalar_type dt, scalar_type nu,
				      const getfem::mesh_region &rg) {
    getfem::generic_assembly assem;
    improved_non_reflective_bc_nonlinear_term<VEC2> nlterm(mf_u, Un0, dt, nu);
    assem.set("V(#2)+=comp(NonLin(#1).vBase(#2))(i,:,i)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_mult);  
    assem.push_vec(VV);
    assem.push_nonlinear_term(&nlterm);
    assem.assembly(rg);
  }



}
#endif
