/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_modeling.h :  Defines model bricks to build           */
/*                                 complete models.                        */
/*     									   */
/* Date : June 15, 2004.                                                   */
/* Authors : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                     */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Yves Renard.                                        */
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


#ifndef GETFEM_MODELING_H__
#define GETFEM_MODELING_H__

#include <getfem_assembling.h>
#include <gmm_solver_cg.h>
#include <gmm_solver_gmres.h>
#include <gmm_precond_ildlt.h>
#include <gmm_precond_ilu.h>
namespace getfem {

  // TODO : gestion des changements : changement dans un mesh fem,
  //   changement d'un second membre, ...

  /* ******************************************************************** */
  /*		Generic definitions.                                      */
  /* ******************************************************************** */
 
  template<typename MODEL_STATE> class mdbrick_abstract;

  template<typename T_MATRIX, typename C_MATRIX, typename VECTOR>
  class model_state {

  protected :
    T_MATRIX tangent_matrix_;
    C_MATRIX constraints_matrix_;
    VECTOR state_, residu_, constraints_rhs_;
    long ident_;

  public :

    typedef T_MATRIX tangent_matrix_type;
    typedef C_MATRIX constraints_matrix_type;
    typedef VECTOR vector_type;
    typedef typename gmm::linalg_traits<VECTOR>::value_type value_type;

    const T_MATRIX &tangent_matrix(void) const 
    { return tangent_matrix_; }
    T_MATRIX &tangent_matrix(void) { return tangent_matrix_; }
    const C_MATRIX &constraints_matrix(void) const 
    { return constraints_matrix_; }
    C_MATRIX &constraints_matrix(void) { return constraints_matrix_; }
    const VECTOR &constraints_rhs(void) const  { return constraints_rhs_; }
    VECTOR &constraints_rhs(void)  { return constraints_rhs_; }
    const VECTOR &state(void) const  { return state_; }
    VECTOR &state(void)  { return state_; }
    const VECTOR &residu(void) const  { return residu_; }
    VECTOR &residu(void) { return residu_; }
    long ident(void) { return ident_; }
    void touch(void) { ident_ = context_dependencies::new_ident(); }
    void adapt_sizes(mdbrick_abstract<model_state> &problem) {
      size_type ndof = problem.nb_dof(), nc = problem.nb_constraints();
      gmm::resize(tangent_matrix_, ndof, ndof);
      gmm::resize(constraints_matrix_, nc, ndof);
      gmm::resize(constraints_rhs_, nc);
      gmm::resize(state_, ndof);
      gmm::resize(residu_, ndof);
    } 

    model_state(void) { ident_ = context_dependencies::new_ident(); }
  };

  template<typename MODEL_STATE>
  class mdbrick_abstract : public context_dependencies {
  public :
    virtual size_type nb_dof(void) = 0;
    virtual size_type nb_constraints(void) = 0;
    virtual void constraints_system(MODEL_STATE &MS, size_type i0 = 0,
				    size_type j0 = 0)= 0;
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0) = 0;
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0) = 0;
    virtual mesh_fem &main_mesh_fem(void) = 0;
    virtual bool is_linear(void) = 0;
    virtual bool is_coercive(void) = 0;
    virtual ~mdbrick_abstract() {}
  };

  typedef gmm::rsvector<scalar_type> modeling_standard_sparse_vector;
  typedef gmm::col_matrix<modeling_standard_sparse_vector>
                                    modeling_standard_sparse_matrix;
  typedef std::vector<scalar_type> modeling_standard_plain_vector;
  typedef model_state<modeling_standard_sparse_matrix,
		      modeling_standard_sparse_matrix,
		      modeling_standard_plain_vector > standard_model_state;


  /* ******************************************************************** */
  /*		Linearized elasticity bricks.                             */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Hooke_linearized_elasticity
    : public mdbrick_abstract<MODEL_STATE> {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;
    mesh_fem &mf_u;
    mesh_fem &mf_data;
    VECTOR lambda_, mu_;
    long ident_ms;
    size_type i0_stored;
  public :

    virtual bool is_linear(void) { return true; }
    virtual bool is_coercive(void) { return true; }
    virtual size_type nb_dof(void) { return mf_u.nb_dof(); }
    virtual size_type nb_constraints(void) { return 0; }
    virtual void constraints_system(MODEL_STATE &, size_type = 0,
				    size_type = 0) {}
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0) {
      if (ident_ms != MS.ident() || this->context_changed()) {
	gmm::sub_interval SUBI(i0, nb_dof());
	asm_stiffness_matrix_for_linear_elasticity
	  (gmm::sub_matrix(MS.tangent_matrix(), SUBI), mf_u, mf_data, 
	   lambda_, mu_);	
	ident_ms = MS.ident();
	i0_stored = i0;
      }
    }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0) {
      compute_tangent_matrix(MS, i0);
      gmm::sub_interval SUBI(i0, nb_dof());
      gmm::mult(gmm::sub_matrix(MS.tangent_matrix(), SUBI),
		gmm::sub_vector(MS.state(), SUBI),
		gmm::sub_vector(MS.residu(), SUBI));
      i0_stored = i0;
    }
    virtual mesh_fem &main_mesh_fem(void) { return mf_u; }

    VECTOR &lambda(void) { return lambda_; }
    VECTOR &mu(void) { return mu_; }
    const VECTOR &lambda(void) const { return lambda_; }
    const VECTOR &mu(void) const { return mu_; }

    template<typename VECT> void get_displacement(MODEL_STATE &MS, VECT &V) {
      gmm::sub_interval SUBI(i0_stored, nb_dof());
      gmm::copy(gmm::sub_vector(MS.state(), SUBI), V);
    }

    mdbrick_Hooke_linearized_elasticity(mesh_fem &mf_u_, mesh_fem &mf_data_,
       value_type lambdai = value_type(1), value_type mui = value_type(1))
      : mf_u(mf_u_), mf_data(mf_data_), ident_ms(-1), i0_stored(0) {
      gmm::resize(lambda_, mf_data_.nb_dof());
      gmm::resize(mu_, mf_data_.nb_dof());
      std::fill(gmm::vect_begin(lambda_), gmm::vect_end(lambda_), lambdai);
      std::fill(gmm::vect_begin(mu_), gmm::vect_end(mu_), mui);
      this->add_dependency(mf_u); this->add_dependency(mf_data); 
    }

  };

  /* A faire : tenseur d'élasticité quelconque */


  /* ******************************************************************** */
  /*		Source term brick.                                        */
  /* ******************************************************************** */

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_source_term : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data;
    VECTOR B_;
    VECTOR F_;
    size_type boundary, qmult;
    int ident_ms;
    
    void fixing_dimensions(void) {
      size_type q = mf_data.get_qdim();
      size_type qdim = sub_problem.main_mesh_fem().get_qdim();
      if (qdim != q && q != 1)
	DAL_THROW(dimension_error,"incompatible dimension of mesh_fem"
		  " structure for source term ");
      qmult = qdim / q;
      gmm::resize(B_, mf_data.nb_dof() * qmult);
    }

  public :
    
    virtual bool is_linear(void) { return sub_problem.is_linear(); }
    virtual bool is_coercive(void) { return sub_problem.is_coercive(); }
    virtual size_type nb_dof(void) { return sub_problem.nb_dof(); }
    virtual size_type nb_constraints(void)
    { return sub_problem.nb_constraints(); }
    virtual void constraints_system(MODEL_STATE &MS, size_type i0 = 0,
				    size_type j0 = 0) 
    { sub_problem.constraints_system(MS, i0, j0); }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0)
    { sub_problem.compute_tangent_matrix(MS, i0); }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0) {
      sub_problem.compute_residu(MS, i0);
      if (ident_ms != MS.ident() || this->context_changed()) {
	qmult = sub_problem.main_mesh_fem().get_qdim() / mf_data.get_qdim();
	if (gmm::vect_size(B_) != mf_data.nb_dof() * qmult) 
	  DAL_THROW(failure_error, "The data mesh fem structure has changed, "
		    " You have to change the rhs in that case.");
	gmm::resize(F_, sub_problem.main_mesh_fem().nb_dof());
	asm_source_term(F_, sub_problem.main_mesh_fem(),mf_data, B_,boundary);
	ident_ms = MS.ident();
      }
      gmm::add(gmm::scaled(F_, value_type(-1)), gmm::sub_vector(MS.residu(),
	       gmm::sub_interval(i0, sub_problem.main_mesh_fem().nb_dof())));
    }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    void changing_rhs(const VECTOR &B__)
    { fixing_dimensions(); gmm::copy(B__, B_); ident_ms = -1; }

    // Constructor which does not define the rhs
    mdbrick_source_term(mdbrick_abstract<MODEL_STATE> &problem,
		       mesh_fem &mf_data_, size_type bound = size_type(-1))
      : sub_problem(problem), mf_data(mf_data_), boundary(bound), ident_ms(-1){
      fixing_dimensions();
      gmm::clear(B_);
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
    }

    // Constructor defining the rhs
    mdbrick_source_term(mdbrick_abstract<MODEL_STATE> &problem,
		       mesh_fem &mf_data_, const VECTOR &B__,
		       size_type bound = size_type(-1))
      : sub_problem(problem), mf_data(mf_data_), boundary(bound), ident_ms(-1){
      fixing_dimensions();
      gmm::copy(B__, B_);
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
    }

  };


  /* ******************************************************************** */
  /*		Dirichlet condition bricks.                               */
  /* ******************************************************************** */

  // faire la version plus complète avec les matrices locales.

  template<typename MODEL_STATE = standard_model_state>
  class mdbrick_Dirichlet : public mdbrick_abstract<MODEL_STATE>  {
    
    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type R;

    mdbrick_abstract<MODEL_STATE> &sub_problem;
    mesh_fem &mf_data;
    VECTOR B_;
    C_MATRIX G;
    VECTOR CRHS;
    size_type boundary, nb_const;
    int ident_ms;

    void fixing_dimensions(void) {
      size_type q = mf_data.get_qdim();
      size_type qdim = sub_problem.main_mesh_fem().get_qdim();
      if (qdim != q && q != 1)
	DAL_THROW(dimension_error,"incompatible dimension of mesh_fem"
		  " structure for source term ");
      size_type qmult = qdim / q;
      gmm::resize(B_, mf_data.nb_dof() * qmult);
    }

  public :
    
    virtual bool is_linear(void) { return sub_problem.is_linear(); }
    virtual bool is_coercive(void) { return sub_problem.is_coercive(); }
    virtual size_type nb_dof(void) { return sub_problem.nb_dof(); }
    void compute_constraints(void) {
      mesh_fem &mf_u = sub_problem.main_mesh_fem();
      size_type nd = mf_u.nb_dof();
      C_MATRIX H(nd, nd);
      VECTOR V(nd);
      asm_dirichlet_constraints(H, V, sub_problem.main_mesh_fem(),
				mf_data, B_, boundary);
      gmm::clean(H, gmm::mat_maxnorm(H) * gmm::default_tol(value_type())
		 * value_type(100));
      gmm::rsvector<value_type> ei(nd), hi(nd);
      std::vector<size_type> ind(0);
      
      dal::bit_vector bdof = mf_u.dof_on_boundary(boundary);
      for (size_type i = bdof.take_first(); i != size_type(-1); i << bdof)
        ind.push_back(i);
      nb_const = ind.size();
      gmm::resize(G, nb_const, nd);
      gmm::copy(gmm::sub_matrix(H, gmm::sub_index(ind),
				gmm::sub_interval(0, nd)), G);
      gmm::resize(CRHS, nb_const);
      gmm::copy(gmm::sub_vector(V, gmm::sub_index(ind)), CRHS);
    }

    virtual size_type nb_constraints(void) {
      return sub_problem.nb_constraints() + nb_const;
    }
    virtual void constraints_system(MODEL_STATE &MS, size_type i0 = 0,
				    size_type j0 = 0) {
      sub_problem.constraints_system(MS, i0, j0);
      if (this->context_changed()) compute_constraints();

      size_type nd = sub_problem.main_mesh_fem().nb_dof();
      size_type ncs = sub_problem.nb_constraints();
      gmm::sub_interval SUBI(j0+ncs,nb_const), SUBJ(i0, nd);
      gmm::copy(G, gmm::sub_matrix(MS.constraints_matrix(), SUBI, SUBJ));
      gmm::mult(G, gmm::scaled(gmm::sub_vector(MS.state(), SUBJ),
			       value_type(-1)),
		CRHS, gmm::sub_vector(MS.constraints_rhs(), SUBI));
    }
    virtual void compute_tangent_matrix(MODEL_STATE &MS, size_type i0 = 0)
    { sub_problem.compute_tangent_matrix(MS, i0); }
    virtual void compute_residu(MODEL_STATE &MS, size_type i0 = 0)
    { sub_problem.compute_residu(MS, i0); }
    virtual mesh_fem &main_mesh_fem(void)
    { return sub_problem.main_mesh_fem(); }

    void changing_rhs(const VECTOR &B__) {
      if (this->context_changed()) fixing_dimensions();
      gmm::copy(B__, B_); compute_constraints();
    }

    // Constructor which does not define the rhs
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		     mesh_fem &mf_data_, size_type bound)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	ident_ms(-1) {
      fixing_dimensions();
      gmm::clear(B_); 
      compute_constraints();
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
    }

    // Constructor defining the rhs
    mdbrick_Dirichlet(mdbrick_abstract<MODEL_STATE> &problem,
		     mesh_fem &mf_data_, const VECTOR &B__,
		     size_type bound)
      : sub_problem(problem), mf_data(mf_data_), boundary(bound),
	ident_ms(-1) {
      fixing_dimensions();
      gmm::copy(B__, B_); 
      compute_constraints();
      this->add_dependency(mf_data);
      this->add_dependency(sub_problem.main_mesh_fem());
    }
    
  };

  
  /* ******************************************************************** */
  /*		Generic solvers.                                          */
  /* ******************************************************************** */

  // faire une version avec using_cg, using_gmres ... (appelée par celle-ci)

  template <typename MODEL_STATE> void
  standard_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
	gmm::iteration &iter) {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;

    size_type ndof = problem.nb_dof(), nc = problem.nb_constraints();
    
    MS.adapt_sizes(problem);
    gmm::fill_random(MS.state());

    problem.compute_tangent_matrix(MS);
    problem.compute_residu(MS);
    problem.constraints_system(MS);
    assert(problem.is_linear());

    VECTOR d(ndof);

    if (nc > 0) { // Take the constraints into account if any.
 
      gmm::col_matrix< gmm::rsvector<value_type> > NS(ndof, ndof);
      VECTOR Ud(ndof);
      
      size_type nbcols=getfem::Dirichlet_nullspace(MS.constraints_matrix(),
						NS, MS.constraints_rhs(), Ud);
      cout << "Nombre d'inconnues sur le systeme reduit : "<< nbcols<< endl;
      VECTOR dr(nbcols), f(nbcols);
      gmm::resize(NS, ndof, nbcols);
      T_MATRIX SM(nbcols, nbcols);
      if (nbcols != ndof) {
	VECTOR RHaux(ndof);
	gmm::mult(MS.tangent_matrix(), Ud, MS.residu(), RHaux);
	gmm::mult(gmm::transposed(NS), gmm::scaled(RHaux, value_type(-1)), f);
	T_MATRIX SMaux(nbcols, ndof);
	gmm::col_matrix< gmm::rsvector<value_type> >
	  NST(gmm::mat_ncols(NS), gmm::mat_nrows(NS));
	gmm::copy(gmm::transposed(NS), NST);
	gmm::mult(NST, MS.tangent_matrix(), SMaux);
	gmm::mult(SMaux, NS, SM);
      }
      else gmm::copy(MS.tangent_matrix(), SM);
      
      gmm::clear(dr);
      if (problem.is_coercive()) {
	gmm::ildlt_precond<T_MATRIX> P(SM);
	gmm::cg(SM, dr, f, P, iter);
      } else {
	gmm::ilu_precond<T_MATRIX> P(SM);
	gmm::gmres(SM, dr, f, P, 100, iter);
      }
      gmm::mult(NS, dr, Ud, d);
    }
    else {
      gmm::clear(d);
      if (problem.is_coercive()) {
	gmm::ildlt_precond<T_MATRIX> P(MS.tangent_matrix());
	gmm::cg(MS.tangent_matrix(), d,
		gmm::scaled(MS.residu(), value_type(-1)), P, iter);
      } else {
	gmm::ilu_precond<T_MATRIX> P(MS.tangent_matrix());
	gmm::gmres(MS.tangent_matrix(), d,
		   gmm::scaled(MS.residu(), value_type(-1)), P, 100, iter);
      }
    }
    
    gmm::add(d, MS.state());
  }


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MODELING_H__  */
