/*===========================================================================
 
 Copyright (C) 2000-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/


#include "getfem/getfem_models.h"
#include "getfem/getfem_plasticity.h"
#include "getfem/getfem_interpolation.h"


namespace getfem {


  /** Compute the projection of D*e + sigma_bar_ 
      on the dof of sigma. */
  class elastoplasticity_nonlinear_term : public nonlinear_elem_term {
    
  protected:
    base_vector params;
    base_vector coeff_precalc;
    size_type N, previous_cv;
    const mesh_im &mim;
    const mesh_fem &mf_u;
    const mesh_fem &mf_sigma;
    const mesh_fem *mf_data;
    std::vector<scalar_type> Sigma_n; 
    std::vector<scalar_type> *Sigma_np1;
    std::vector<scalar_type> U_n;
    std::vector<scalar_type> U_np1;
    std::vector<scalar_type> threshold, lambda, mu;  
    bgeot::multi_index sizes_;
    const abstract_constraints_projection &t_proj;
    std::vector<scalar_type> *saved_plast;
    fem_precomp_pool fppool;
    std::vector<scalar_type> stored_proj;
    const size_type flag_proj;
    bool write_sigma_np1;
    bool write_plast;
    

  
  public:  


    // constructor
    elastoplasticity_nonlinear_term(const mesh_im &mim_,
       	    const mesh_fem &mf_u_,
            const mesh_fem &mf_sigma_,
            const mesh_fem *mf_data_,
            const std::vector<scalar_type> &U_n_, 
	    const std::vector<scalar_type> &U_np1_,
            const std::vector<scalar_type> &Sigma_n_, 
            std::vector<scalar_type> *Sigma_np1_, 
            const std::vector<scalar_type> &threshold_, 
            const std::vector<scalar_type> &lambda_,
            const std::vector<scalar_type> &mu_, 
	    const abstract_constraints_projection  &t_proj_,
	    std::vector<scalar_type> *saved_plast_,
	    const size_type flag_proj_, bool write_sigma_np1_, 
	    bool write_plast_) : 
      mim(mim_), mf_u(mf_u_), mf_sigma(mf_sigma_), 
      Sigma_n(Sigma_n_), Sigma_np1(Sigma_np1_),
      t_proj(t_proj_), saved_plast(saved_plast_), flag_proj(flag_proj_),
      write_sigma_np1(write_sigma_np1_), write_plast(write_plast_) {
      
      params = base_vector(3); 
      N = mf_u_.linked_mesh().dim();
      coeff_precalc = base_vector(N*N*N*N);
      gmm::resize(U_n, mf_u_.nb_basic_dof());
      gmm::resize(U_np1, mf_u_.nb_basic_dof());
      gmm::resize(Sigma_n, mf_sigma_.nb_basic_dof());

      sizes_ = bgeot::multi_index(N, N, N, N);
      mf_u.extend_vector(gmm::sub_vector
			 (U_n_, gmm::sub_interval
			  (0,mf_u_.nb_dof())), U_n);
      mf_u.extend_vector(gmm::sub_vector
			 (U_np1_, gmm::sub_interval
			  (0,mf_u_.nb_dof())), U_np1);
      mf_sigma.extend_vector(gmm::sub_vector
			     (Sigma_n_, gmm::sub_interval
			      (0,mf_sigma_.nb_dof())), Sigma_n);
      

      if (mf_data_ != NULL) {
	gmm::resize(mu, mf_data_->nb_basic_dof());
	gmm::resize(lambda, mf_data_->nb_basic_dof());
	gmm::resize(threshold, mf_data_->nb_basic_dof());
	mf_data = mf_data_;
	mf_data->extend_vector(threshold_, threshold);
	mf_data->extend_vector(lambda_, lambda);
	mf_data->extend_vector(mu_, mu);

	
      } else {
	gmm::resize(mu, 1); mu[0]  =  mu_[0];
	gmm::resize(lambda, 1); lambda[0]  =  lambda_[0];
	gmm::resize(threshold, 1); 
	threshold[0] =  threshold_[0];
	mf_data = mf_data_;
	
      }
      GMM_ASSERT1(mf_u.get_qdim() == N, 
		  "wrong qdim for the mesh_fem"); 
      
      if (flag_proj==0) sizes_.resize(2);

      // used to know if the current element is different 
      // than the previous one and so if a new computation 
      // is necessary or not.
      previous_cv = size_type(-1);

    }



    const bgeot::multi_index &sizes() const { return sizes_; }



    // method from nonlinear_elem_term, 
    // gives on output the tensor
    virtual void compute(fem_interpolation_context& ctx,
			 bgeot::base_tensor &t){ 
      size_type cv = ctx.convex_num();//index of current element
      size_type qdim = mf_u.get_qdim();
      size_type qdim_sigma = mf_sigma.get_qdim();
      pfem pf_sigma = ctx.pf();
      size_type nbd_sigma = pf_sigma->nb_dof(cv);
      
      GMM_ASSERT1(pf_sigma->is_lagrange(),
		  "Sorry, works only for Lagrange fems");
      
      // size_proj is different if we compute the projection 
      // or the gradient of the projection
      size_type size_proj = qdim * qdim * 
	(flag_proj == 1 ? qdim * qdim : 1);

      // if the current element is different than the previous one
      if(previous_cv != cv) {

	stored_proj.resize(nbd_sigma * size_proj);
	base_matrix G;
	bgeot::vectors_to_base_matrix
	  (G, mf_u.linked_mesh().points_of_convex(cv));
	bgeot::pgeometric_trans pgt=
	  mf_u.linked_mesh().trans_of_convex(cv);

	fem_interpolation_context ctx_data;
	pfem pf_data;
	base_vector coeff_data;

	// if the Lame coefficient are vector fields
	if (mf_data) {
	 
	  pf_data = mf_data->fem_of_element(cv); 
	  size_type nbd_data = pf_data->nb_dof(cv);
	  coeff_data.resize(nbd_data*3);

	  // Definition of the Lame coeff
	  mesh_fem::ind_dof_ct::const_iterator itdof
	    = mf_data->ind_basic_dof_of_element(cv).begin();

	  for (size_type k = 0; k < nbd_data; ++k, ++itdof) {
	    coeff_data[k*3] = lambda[*itdof];
	    coeff_data[k*3+1] = mu[*itdof];
	    coeff_data[k*3+2] = threshold[*itdof];
	  } 
	  GMM_ASSERT1(pf_data->target_dim() == 1,
		      "won't interpolate on a vector FEM... ");

	  pfem_precomp pfp_data = 
	    fppool(pf_data, pf_sigma->node_tab(cv));
	  ctx_data = fem_interpolation_context
	    (pgt,pfp_data,size_type(-1), G, cv,size_type(-1));

	} else {
	  
	  params[0] = lambda[0];
	  params[1] = mu[0];
	  params[2] = threshold[0];

	}

	// for each dof of the current element
	for (size_type ii = 0; ii < nbd_sigma; ++ii) {

	  // interpolation of the data on sigma dof
	  if (mf_data) {
	    ctx_data.set_ii(ii);
	    pf_data->interpolation(ctx_data, coeff_data, params, 3);
	  } 
	  
	  std::vector<scalar_type> coeff_u_n, coeff_u_np1;
	  
	  pfem pf_u = mf_u.fem_of_element(cv);
	  size_type cvnbdof_u = 
	      mf_u.nb_basic_dof_of_element(cv); 
	  
	  
	  // Definition of the coeff for u_n and u_np1	  
	  coeff_u_n.resize(cvnbdof_u);
	  coeff_u_np1.resize(cvnbdof_u);
	  mesh_fem::ind_dof_ct::const_iterator itdof
	    = mf_u.ind_basic_dof_of_element(cv).begin();
	  for (size_type k = 0; k < cvnbdof_u; ++k, ++itdof) {
	    coeff_u_n[k] = U_n[*itdof];
	    coeff_u_np1[k] = U_np1[*itdof];
	  }
	  
	  base_matrix G_u_n(qdim, qdim), G_u_np1(qdim, qdim);
	  pfem_precomp pfp_u = fppool(pf_u, pf_sigma->node_tab(cv));
	  fem_interpolation_context ctx_u
	    (pgt,pfp_u,size_type(-1), G, cv,size_type(-1));
	    
	  // interpolation of the gradient of u_n and u_np1 
	  // on sigma dof
	  ctx_u.set_ii(ii);
	  pf_u->interpolation_grad
	    (ctx_u, coeff_u_n, G_u_n, dim_type(qdim));
	  pf_u->interpolation_grad
	    (ctx_u, coeff_u_np1, G_u_np1, dim_type(qdim));
	  
	  // Compute sigma_hat = D*esp_np1 - D*eps_n + sigma_n
	  // where D represent the elastic stiffness tensor
	  base_matrix sigma_hat(qdim, qdim);
	    
	  
	  // Compute lambda*tr(esp_n) and lambda*tr(esp_np1)
	  scalar_type ltrace_eps_n 
	    =  params[0]*gmm::mat_trace(G_u_n);
	  scalar_type ltrace_eps_np1 
	    = params[0]*gmm::mat_trace(G_u_np1);
	  
	  // Compute sigma_hat
	  size_type idof_sigma
	    = mf_sigma.ind_basic_dof_of_element(cv)[ii*qdim_sigma];
	  for(dim_type i = 0; i < qdim; ++i) {
	    for(dim_type j = 0; j < qdim; ++j) {
	      sigma_hat(i,j) = Sigma_n[idof_sigma + j*qdim +i]
		+ params[1]*(G_u_np1(i,j) + G_u_np1(j,i))
		- params[1]*(G_u_n(i,j) + G_u_n(j,i));
	      if (i==j)
		sigma_hat(i,j) += ltrace_eps_np1 - ltrace_eps_n;
	    }
	  }
	    
	  
	  base_matrix proj;

	  // Compute the projection or its grad
	  t_proj.do_projection(sigma_hat,params[2],proj,flag_proj);
	  
	  // Retrieve the projection or its grad
	  std::copy(proj.begin(), proj.end(),
		    stored_proj.begin() + proj.size() * ii);
	  
	  // Retrieve the new stress constraints values
	  // used only by the function 'write_sigma(...)'
	  if (flag_proj == 0 && write_sigma_np1) {
	    
	    for(dim_type i = 0; i < qdim; ++i){
	      for(dim_type j = 0; j < qdim; ++j){
		(*Sigma_np1)[idof_sigma + j*qdim + i] = proj(i,j);
	      }
	    }
	  }

	  // Compute the plastic part
	  if (flag_proj == 0 && write_plast) {
	    
	    for(dim_type i = 0; i < qdim; ++i){
	      for(dim_type j = 0; j < qdim; ++j){
		(*saved_plast)[idof_sigma + j*qdim +i] =
		  proj(i,j) - params[1]*(G_u_np1(i,j) + G_u_np1(j,i));
		if (i==j)
		  (*saved_plast)[idof_sigma + j*qdim +i] -= ltrace_eps_np1;
	      }
	    }
	  }
	}
	// upload previous_cv
	previous_cv = cv;
      }

      // interpolation of the projection on sigma dof
      coeff_precalc.resize(size_proj);
      pf_sigma->interpolation(ctx, stored_proj, coeff_precalc,
			      dim_type(size_proj));
      
      t.adjust_sizes(sizes_);
      
      // copy the result into the tensor returned t
      std::copy(coeff_precalc.begin(), coeff_precalc.end(), 
		t.begin());

    }

};





  /** 
     Right hand side vector for elastoplasticity 
      @ingroup asm
  */
  template<typename VECT> 
  void asm_elastoplasticity_rhs (VECT &V, 
	const mesh_im &mim, 
	const mesh_fem &mf_u, 
	const mesh_fem &mf_sigma, 
	const mesh_fem &mf_data, 
	const VECT &u_n,
	const VECT &u_np1, 
	const VECT &sigma_n, 
	VECT *sigma_np1, 
	const VECT &lambda, 
	const VECT &mu, 
	const VECT &threshold, 
        const abstract_constraints_projection  &t_proj,
	VECT *saved_plast,
	bool write_sigma_np1,
	bool write_plast,
	const mesh_region &rg = mesh_region::all_convexes()) {

    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");


    elastoplasticity_nonlinear_term plast(mim, mf_u, mf_sigma,
					  &mf_data, u_n, u_np1,
					  sigma_n, sigma_np1, 
					  threshold, lambda, mu, 
					  t_proj, saved_plast, 0, 
					  write_sigma_np1, write_plast);


    generic_assembly assem("V(#1) + =comp(NonLin(#2).vGrad(#1))(i,j,:,i,j);");

 
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_sigma);
    assem.push_nonlinear_term(&plast);
    assem.push_vec(V);
    assem.assembly(rg);


  }
  




  /** 
      Tangent matrix for elastoplasticity
      @ingroup asm
  */
  template<typename MAT,typename VECT> 
  void asm_elastoplasticity_tangent_matrix(MAT &H, 
	const mesh_im &mim, 
	const mesh_fem &mf_u, 
	const mesh_fem &mf_sigma,
	const mesh_fem &mf_data, 
	const VECT &u_n,
	const VECT &u_np1, 
	const VECT &sigma_n, 
	const VECT &lambda, 
	const VECT &mu, 
	const VECT &threshold, 
        const abstract_constraints_projection &t_proj,
        const mesh_region &rg = mesh_region::all_convexes()) {


    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    elastoplasticity_nonlinear_term gradplast(mim, mf_u, mf_sigma,
					&mf_data, u_n, u_np1,
					sigma_n, 0, 
					threshold, lambda, mu,
					      t_proj, 0, 1, false, false);

    generic_assembly assem;

    if (&(mf_data)!=NULL) {

      assem.set("lambda=data$1(#3); mu=data$2(#3);"
		"t=comp(NonLin(#2).vGrad(#1).vGrad(#1).Base(#3))(i,j,:,:,:,:,:,:,i,j,:);"
		"M(#1,#1)+=  sym(t(k,l,:,l,k,:,m).mu(m)+t(k,l,:,k,l,:,m).mu(m)+t(k,k,:,l,l,:,m).lambda(m))");
      
    } else {
      
      assem.set("lambda=data$1(1); mu=data$2(1);"
		"t=comp(NonLin(#2).vGrad(#1).vGrad(#1))(i,j,:,:,:,:,:,:,i,j);"
		"M(#1,#1)+= sym(t(k,l,:,l,k,:).mu(1)+t(k,l,:,k,l,:).mu(1)+t(k,k,:,l,l,:).lambda(1))");
    }
    
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_sigma);
    if (&(mf_data)!=NULL)
      assem.push_mf(mf_data);
    assem.push_data(lambda);
    assem.push_data(mu);
    assem.push_nonlinear_term(&gradplast);
    assem.push_mat(H);
    assem.assembly(rg);
    
  }



  //=================================================================
  //
  //  Elastoplasticity Brick
  //
  //=================================================================

  struct elastoplasticity_brick : public virtual_brick {
    
    const abstract_constraints_projection  &t_proj;
	

    virtual void asm_real_tangent_terms(const model &md, 
					size_type /* ib */,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,   
					model::real_veclist &vecl,
					model::real_veclist &,
					size_type region,
					build_version version)const {


      GMM_ASSERT1(mims.size() == 1,
		  "Elastoplasticity brick need a single mesh_im");
      GMM_ASSERT1(vl.size() == 1,
		  "Elastoplasticity brick need one variable"); 
      /** vl[0] = u */
      
      GMM_ASSERT1(dl.size() == 4,
		  "Wrong number of data for elastoplasticity brick, "
                  << dl.size() << " should be 4.");
      GMM_ASSERT1(matl.size() == 1,  "Wrong number of terms for "
		  "elastoplasticity brick");


      const model_real_plain_vector &u_np1 = 
	md.real_variable(vl[0], 0);
      const model_real_plain_vector &u_n = 
	md.real_variable(vl[0], 1);
      const mesh_fem &mf_u = 
	*(md.pmesh_fem_of_variable(vl[0]));

      const model_real_plain_vector &lambda = 
	md.real_variable(dl[0]);
      const model_real_plain_vector &mu = 
	md.real_variable(dl[1]);
      const model_real_plain_vector &threshold = 
	md.real_variable(dl[2]);
      const mesh_fem *mf_data = 
	md.pmesh_fem_of_variable(dl[0]);

      const model_real_plain_vector &sigma_n = 
	md.real_variable(dl[3]);
      const mesh_fem &mf_sigma = 
	*(md.pmesh_fem_of_variable(dl[3]));
      GMM_ASSERT1(!(mf_sigma.is_reduced()),
		  "Works only for pure Lagrange fems");
	
      const mesh_im &mim = *mims[0];
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	asm_elastoplasticity_tangent_matrix
	  (matl[0], mim, mf_u, mf_sigma, *mf_data, u_n,
	   u_np1, sigma_n, lambda, mu, threshold, t_proj, rg);
      }
      
      if (version & model::BUILD_RHS) {
	asm_elastoplasticity_rhs
	  (vecl[0], mim, mf_u, mf_sigma, *mf_data, u_n,
	   u_np1, sigma_n, (model_real_plain_vector *)(0), 
	   lambda, mu, threshold, t_proj, (model_real_plain_vector *)(0),
	   false, false, rg);
	gmm::scale(vecl[0], scalar_type(-1));
      }

    }




    elastoplasticity_brick(const abstract_constraints_projection &t_proj_)
      : t_proj(t_proj_){
      set_flags("Elastoplasticity brick", false /* is linear*/,
		true /* is symmetric */, false /* is coercive */,
		true /* is real */, false /* is complex */);
    }

  };
  





  //=================================================================
  //  Add a elastoplasticity brick
  //=================================================================

  size_type add_elastoplasticity_brick(model &md, 
	const mesh_im &mim, 
	const abstract_constraints_projection &ACP,
	const std::string &varname,
	const std::string &datalambda,
	const std::string &datamu,
	const std::string &datathreshold, 
	const std::string &datasigma,
	size_type region) {
    pbrick pbr = new elastoplasticity_brick(ACP);
    
    model::termlist tl;
    tl.push_back(model::term_description
		 (varname, varname, true));
    model::varnamelist dl(1, datalambda);
    dl.push_back(datamu); 
    dl.push_back(datathreshold); 
    dl.push_back(datasigma);
    model::varnamelist vl(1, varname);
    
    return md.add_brick(pbr, vl, dl, tl, 
			model::mimlist(1,&mim), region);
  }





  //=================================================================
  //  New stress constraints values computation and saved 
  //  Update of u and sigma on time iterates :
  //  u_np1 -> u_n         sigma_np1 -> sigma_n
  //=================================================================
  
  
  void elastoplasticity_next_iter(model &md, 
		   const mesh_im &mim,
		   const std::string &varname,
		   const abstract_constraints_projection &ACP,
		   const std::string &datalambda,
		   const std::string &datamu, 
		   const std::string &datathreshold, 
		   const std::string &datasigma) {
   

    const model_real_plain_vector &u_np1 = 
      md.real_variable(varname, 0);
    model_real_plain_vector &u_n = 
      md.set_real_variable(varname, 1);
    const mesh_fem &mf_u = 
      *(md.pmesh_fem_of_variable(varname));
    
    const model_real_plain_vector &lambda = 
      md.real_variable(datalambda);
    const model_real_plain_vector &mu = 
      md.real_variable(datamu);
    const model_real_plain_vector &threshold = 
      md.real_variable(datathreshold);
    const mesh_fem *mf_data = 
    md.pmesh_fem_of_variable(datalambda);

    const model_real_plain_vector &sigma_n = 
      md.real_variable(datasigma);
    const mesh_fem &mf_sigma = 
      *(md.pmesh_fem_of_variable(datasigma));
    
    unsigned N = unsigned(mf_sigma.linked_mesh().dim());

    std::vector<scalar_type> sigma_np1
      (mf_sigma.nb_dof()*N*N/mf_sigma.get_qdim());

    std::vector<scalar_type> V(mf_u.nb_dof());

    mesh_region rg = mim.linked_mesh().get_mpi_region();

    asm_elastoplasticity_rhs
      (V, mim, mf_u, mf_sigma, *mf_data, u_n,
       u_np1, sigma_n, &sigma_np1, 
       lambda, mu, threshold, ACP, (model_real_plain_vector *)(0), true, false, rg);

    // upload sigma and u : u_np1 -> u_n, sigma_np1 -> sigma_n 
    // be careful to use this function 
    // only if the computation is over
    MPI_SUM_VECTOR(sigma_np1, md.set_real_variable(datasigma));
    gmm::copy(u_np1, u_n);

 
  }



  //=================================================================
  //  Von Mises or Tresca stress computation for elastoplasticity 
  //=================================================================

  void compute_elastoplasticity_Von_Mises_or_Tresca(model &md, 
	const std::string &datasigma,
	const mesh_fem &mf_vm,
	model_real_plain_vector &VM,
	bool tresca) {

    GMM_ASSERT1(gmm::vect_size(VM) == mf_vm.nb_dof(),
		"The vector has not the good size");

    const model_real_plain_vector &sigma_np1 = 
      md.real_variable(datasigma, 0);
    const mesh_fem &mf_sigma = 
      *(md.pmesh_fem_of_variable(datasigma));
    
    // dimension of the finite element used
    unsigned N = unsigned(mf_sigma.linked_mesh().dim());
    
    GMM_ASSERT1(mf_vm.get_qdim() == 1,
		"Target dimension of mf_vm should be 1");

    base_matrix sigma(N, N), Id(N, N);
    base_vector eig(N);
    base_vector sigma_vm(mf_vm.nb_dof()*N*N);

    gmm::copy(gmm::identity_matrix(), Id);

    interpolation(mf_sigma, mf_vm, sigma_np1, sigma_vm);
    
    // for each dof we compute the Von Mises or Tresca stress
    for (size_type ii = 0; ii < mf_vm.nb_dof(); ++ii) {

      /* we retrieve the matrix sigma_vm on this dof */
      std::copy(sigma_vm.begin()+ii*N*N, sigma_vm.begin()+(ii+1)*N*N,
		sigma.begin());
      
      if (!tresca) {
	/* von mises: norm(deviator(sigma)) */
	gmm::add(gmm::scaled(Id, -gmm::mat_trace(sigma) / N), sigma);
	
	/* von mises stress=sqrt(3/2)* norm(sigma) */
	VM[ii] = sqrt(3.0/2.)*gmm::mat_euclidean_norm(sigma);
      } else {
	/* else compute the tresca criterion */
	gmm::symmetric_qr_algorithm(sigma, eig);
	std::sort(eig.begin(), eig.end());
	VM[ii] = eig.back() - eig.front();
      }
    }
  }
  


  //=================================================================
  //  Compute the plastic part  
  //=================================================================
  

  void compute_plastic_part(model &md, 
		   const mesh_im &mim,
		   const mesh_fem &mf_pl,
		   const std::string &varname,
		   const abstract_constraints_projection &ACP,
		   const std::string &datalambda,
		   const std::string &datamu, 
		   const std::string &datathreshold, 
		   const std::string &datasigma,
		   model_real_plain_vector &plast) {
   

    const model_real_plain_vector &u_np1 = 
      md.real_variable(varname, 0);
    model_real_plain_vector &u_n = 
      md.set_real_variable(varname, 1);
    const mesh_fem &mf_u = 
      *(md.pmesh_fem_of_variable(varname));
    
    const model_real_plain_vector &lambda = 
      md.real_variable(datalambda);
    const model_real_plain_vector &mu = 
      md.real_variable(datamu);
    const model_real_plain_vector &threshold = 
      md.real_variable(datathreshold);
    const mesh_fem *mf_data = 
      md.pmesh_fem_of_variable(datalambda);

    const model_real_plain_vector &sigma_n = 
      md.real_variable(datasigma);
    const mesh_fem &mf_sigma = 
      *(md.pmesh_fem_of_variable(datasigma));
    
    unsigned N = unsigned(mf_sigma.linked_mesh().dim());

    mesh_region rg = mesh_region::all_convexes(); // mim.linked_mesh().get_mpi_region();

    std::vector<scalar_type> V(mf_u.nb_dof());
    std::vector<scalar_type> saved_plast(mf_sigma.nb_dof());

    asm_elastoplasticity_rhs
      (V, mim, mf_u, mf_sigma, *mf_data, u_n,
       u_np1, sigma_n, (model_real_plain_vector *)(0), 
       lambda, mu, threshold, ACP, &saved_plast, false, true, rg);

    //MPI_SUM_VECTOR(saved_plast);

    /* Retrieve and save the plastic part */
    GMM_ASSERT1(gmm::vect_size(plast) == mf_pl.nb_dof(),
		"The vector has not the good size");
    
    GMM_ASSERT1(mf_pl.get_qdim() == 1,
		"Target dimension of mf_vm should be 1");

    base_matrix plast_tmp(N, N), Id(N, N);
    base_vector eig(N);
    base_vector saved_pl(mf_pl.nb_dof()*N*N);

    gmm::copy(gmm::identity_matrix(), Id);

    interpolation(mf_sigma, mf_pl, saved_plast, saved_pl);
       
    // for each dof we compute the norm of the plastic part
    for (size_type ii = 0; ii < mf_pl.nb_dof(); ++ii) {

      /* we retrieve the matrix sigma_pl on this dof */
      std::copy(saved_pl.begin()+ii*N*N, saved_pl.begin()+(ii+1)*N*N,
		plast_tmp.begin());
      
      plast[ii] = gmm::mat_euclidean_norm(plast_tmp);
      
    }
 
  }


  
}  /* end of namespace getfem.  */

