// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2010 Yves Renard
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


#include "getfem/getfem_models.h"
#include "getfem/getfem_plasticity.h"


namespace getfem {


  /** Compute the projection of D*e + sigma_bar_ 
      on the dof of sigma. */
  class plasticity_nonlinear_term : public nonlinear_elem_term {
  protected:
    base_vector params;
    base_matrix coeff_precalc;
    size_type N, previous_cv;
    const mesh_im &mim;
    const mesh_fem &mf_u;
    const mesh_fem &mf_sigma;
    const mesh_fem *mf_data;
    std::vector<scalar_type> U_n;
    std::vector<scalar_type> U_np1;
    std::vector<scalar_type> Sigma_n;
    std::vector<scalar_type> Sigma_np1;
    std::vector<scalar_type> threshold, lambda, mu;  
    bgeot::multi_index sizes_;
    const abstract_constraints_projection  &t_proj; 

 
    const size_type flag_proj;

//    bool fill_sigma_bar;

  
  public:  




    // constructor
    plasticity_nonlinear_term(const mesh_im &mim_,
	          const mesh_fem &mf_u_,
	          const mesh_fem &mf_sigma_,
	          const mesh_fem *mf_data_,
	          const std::vector<scalar_type> &U_n_, 
	          const std::vector<scalar_type> &U_np1_,
	          const std::vector<scalar_type> &Sigma_n_, 
		  const std::vector<scalar_type> &Sigma_np1_, 
	     	  const std::vector<scalar_type> &threshold_, 
	      	  const std::vector<scalar_type> &lambda_,
	       	  const std::vector<scalar_type> &mu_, 
	       	  const abstract_constraints_projection  &t_proj_,
		  const size_type flag_proj_):mim(mim_), mf_u(mf_u_), mf_sigma(mf_sigma_), t_proj(t_proj_), flag_proj(flag_proj_){
      GMM_TRACE2("Building the plasticity non linear term");
      

      params = base_vector(3); 
      N = mf_u_.linked_mesh().dim();
      gmm::resize(U_n, mf_u_.nb_basic_dof());
      gmm::resize(U_np1, mf_u_.nb_basic_dof());
      gmm::resize(Sigma_n, mf_sigma_.nb_basic_dof());
      gmm::resize(Sigma_np1, mf_sigma_.nb_basic_dof());
      sizes_ = bgeot::multi_index(N, N, N, N);


      if (mf_data_ != NULL) {
	gmm::resize(mu, mf_data_->nb_basic_dof());
	gmm::resize(lambda, mf_data_->nb_basic_dof());
	gmm::resize(threshold, mf_data_->nb_basic_dof());
	
	mf_u.extend_vector(gmm::sub_vector
	   (U_n_, gmm::sub_interval(0,mf_u_.nb_dof())), U_n);
	mf_u.extend_vector(gmm::sub_vector
	   (U_np1_, gmm::sub_interval(0,mf_u_.nb_dof())), U_np1);
	mf_sigma.extend_vector(gmm::sub_vector
	   (Sigma_n_, gmm::sub_interval(0, mf_sigma_.nb_dof())),Sigma_n);
	mf_sigma.extend_vector(gmm::sub_vector
	   (Sigma_np1_, gmm::sub_interval(0,mf_sigma_.nb_dof())),Sigma_np1);
	mf_data->extend_vector(threshold_, threshold);
	mf_data->extend_vector(lambda_, lambda);
	mf_data->extend_vector(mu_, mu);
	mf_data = mf_data_;
	
      } else {
	gmm::resize(mu, 1); mu[0]  =  mu_[0];
	gmm::resize(lambda, 1); lambda[0]  =  lambda_[0];
	gmm::resize(threshold, 1); 
	threshold[0] =  threshold_[0];
	mf_data = mf_data_;

      }
      

   //   fill_sigma_bar = fill_sigma;   /* always false during resolution, */
      /*                      true when called from compute_constraints */

      GMM_ASSERT1(mf_u.get_qdim() == N, 
		  "wrong qdim for the mesh_fem");      
      if (flag_proj==0) sizes_.resize(2);


      previous_cv = size_type(-1);
	GMM_TRACE2(
	   "End of building the plasticity non linear term");

    }





    const bgeot::multi_index &sizes() const { return sizes_; }





    // compute() method from nonlinear_elem, gives on output the tensor
    virtual void compute(fem_interpolation_context& ctx,
			 bgeot::base_tensor &t){ 

      GMM_TRACE2("Computing the plasticity");

      size_type cv = ctx.convex_num();//index of current element


      size_type qdim = mf_u.get_qdim();
      size_type qdim_sigma = mf_sigma.get_qdim();
      size_type ii = ctx.ii();
      pfem pf = ctx.pf();
      size_type nbd_sigma = pf->nb_dof(cv);
       cout<<"nbd sigma : "<<nbd_sigma<<endl;
      // = 3 : noeuds locaux

      size_type qqdim_u = size_type(gmm::vect_size(U_n)/mf_u.nb_dof()); //same for U_np1   =1

      size_type qmult_u = 1;
      size_type R = nbd_sigma*qdim;

      base_matrix G_u_n(R, R); 
      base_matrix G_u_np1(R, R);

      // cout<<"N : "<<N<<endl;
      // = 2

      if(previous_cv != cv){
	GMM_TRACE2("previous");
	
	if(mf_data){

	  GMM_TRACE2("data");
	  base_vector coeff;
	  fem_precomp_pool fppool_data;
 
	  bgeot::pgeometric_trans pgt=
	    mf_data->linked_mesh().trans_of_convex(cv);

	  pfem pf_data = mf_data->fem_of_element(cv); 

	  size_type nbd_data = pf_data->nb_dof(cv);
	  
 
	  mesh_fem::ind_dof_ct::const_iterator itdof;

	  size_type cvnbdof_data = 
	    mf_data->nb_basic_dof_of_element(cv);

	  // Definition of the coeff of Lame

	  coeff.resize(nbd_data*3);

	  itdof = mf_data->ind_basic_dof_of_element(cv).begin();
	  for (size_type k = 0; k < cvnbdof_data; ++k, ++itdof) {
	    coeff[k*3] = lambda[*itdof];
	    coeff[k*3+1] = mu[*itdof];
	    coeff[k*3+2] = threshold[*itdof];

	  } 

	  base_matrix G;
	  if (pf_data->need_G()) 
	    bgeot::vectors_to_base_matrix
	      (G, mf_data->linked_mesh().points_of_convex(cv));

	  GMM_ASSERT1(pf->target_dim() == 1,
                   "won't interpolate on a vector FEM... ");


	  pfem_precomp pfp_data = 
	    fppool_data(pf_data, pf->node_tab(cv));
	  fem_interpolation_context ctx_data(
	    pgt,pfp_data,size_type(-1), G, cv,size_type(-1));



	  for (size_type i = 0; i < nbd_sigma; ++i) {

	    ctx.set_ii(i);
	    pf_data->interpolation(ctx_data, 
				   coeff, params, qdim);
	  }
	}else{
	  GMM_TRACE2("non data");
	  for (size_type i = 0; i < nbd_sigma; ++i) {

	    ctx.set_ii(i);
	    params[i*3] = lambda[0];
	    params[i*3+1] = mu[0];
	    params[i*3+2] = threshold[0];
	  }
	  cout<<"nbd sigma : "<<nbd_sigma<<endl;
	  // cout<<"taille params : "<<sizeof(params)<<endl;
	  // = 24 = 3*2*4 ???



	  



	  std::vector<std::vector<scalar_type> > coeff_u_n, coeff_u_np1;
	  fem_precomp_pool fppool_u; // same for U_np1
	  
	   
	  
	  bgeot::pgeometric_trans pgt=
	    mf_u.linked_mesh().trans_of_convex(cv);
	  
	  pfem pf_u = mf_u.fem_of_element(cv);
	  pfem pf_sigma = mf_sigma.fem_of_element(cv);
	  
	  size_type nbd_u = pf_u->nb_dof(cv); // = 6
	  
	  mesh_fem::ind_dof_ct::const_iterator itdof;
	  
	 

	  size_type cvnbdof_u = 
	    mf_u.nb_basic_dof_of_element(cv); // = 12

	  
	  // Definition of the coeff u_n and u_np1
	  
	  coeff_u_n.resize(qqdim_u);
	  coeff_u_np1.resize(qqdim_u);

	  
	  for(size_type qq=0; qq<qqdim_u; ++qq){
	    coeff_u_n[qq].resize(nbd_u*qdim);
	    coeff_u_np1[qq].resize(nbd_u*qdim);
	    itdof = mf_u.ind_basic_dof_of_element(cv).begin();
	    for (size_type k = 0; k < cvnbdof_u; ++k, ++itdof) {
	      coeff_u_n[qq][k] = U_n[(*itdof)*qqdim_u+qq];
	      coeff_u_np1[qq][k] = U_np1[(*itdof)*qqdim_u+qq];
	    }
	  }

	
	  base_matrix G;
	  base_matrix G_u_n_temp(qdim,qdim), G_u_np1_temp(qdim, qdim);
 
	  bgeot::vectors_to_base_matrix
	    (G, mf_u.linked_mesh().points_of_convex(cv));
	  //dim G = 2*3

	  GMM_ASSERT1(pf->target_dim() == 1,
		      "won't interpolate on a vector FEM... ");

	  pfem_precomp pfp_u = 
	    fppool_u(pf_u, pf->node_tab(cv));
	  

	  fem_interpolation_context ctx_u(
	     pgt,pfp_u,size_type(-1), G, cv,size_type(-1));


	  // dal::bit_vector dof_sigma_done;
	  //  dof_sigma_done.sup(0, mf_sigma.nb_basic_dof());
	  // itdof = mf_sigma.ind_basic_dof_of_element(cv).begin();//taille = 12 = 3noeuds*4(qdim_sigma)
	 
	  size_type itdof2 = 0;

	  for (size_type i = 0; i<nbd_sigma; ++i, itdof2 += qdim){
	    size_type dof_sigma = itdof2*qmult_u;
	    //    dof_sigma_done.add(itdof);
	    ctx_u.set_ii(i);
	    for(size_type qq=0; qq<qqdim_u; ++qq){
	      pf_u->interpolation_grad(ctx_u, 
		coeff_u_n[qq], G_u_n_temp, qdim);
	      pf_u->interpolation_grad(ctx_u, 
		coeff_u_np1[qq], G_u_np1_temp, qdim);

	      cout<<"col = "<<gmm::mat_ncols(G_u_n_temp)<<", lig = "<<gmm::mat_nrows(G_u_np1_temp)<<endl;

	      for(size_type k = 0; k<qdim; ++k){
		for(size_type j = 0; j<qdim; ++j){
		  G_u_n((dof_sigma+k)*qqdim_u+qq,
			(dof_sigma+j)*qqdim_u+qq) 
		    = G_u_n_temp(k,j);
		  G_u_np1((dof_sigma+k)*qqdim_u+qq,
			  (dof_sigma+j)*qqdim_u+qq) 
		    = G_u_np1_temp(k,j);
		  cout<<"k = "<<k<<", j = "<<j<<", dof_sigma = "<<dof_sigma<<endl;
		}
	      }
	    }
	  } 
	  cout<<"nbd sigma : "<<nbd_sigma<<endl;
	}
	cout<<"nbd sigma : "<<nbd_sigma<<endl;
    GMM_TRACE2("non previous");

	// Compute sigma_hat = D*esp_np1 - D*eps_n + sigma_n
    base_matrix sigma_hat(R,R); //(N,N);

    cout<<"nbd sigma : "<<nbd_sigma<<endl;


	// Compute lambda*tr(esp_n) and lambda*tr(esp_np1)
	scalar_type ltrace_eps_n = 
	  params[0]*gmm::mat_trace(G_u_n);
	scalar_type ltrace_eps_np1 = 
	  params[0]*gmm::mat_trace(G_u_np1);
	cout<<"nbd sigma : "<<nbd_sigma<<endl;

	cout<<"col = "<<gmm::mat_ncols(G_u_n)<<", lig = "<<gmm::mat_nrows(G_u_n)<<", R = "<<R<<", nbd_sigma = "<<nbd_sigma<<endl;
	

	for(dim_type i = 0; i<R; ++i){
	  for(dim_type j = 0; j<R; ++j){
	    sigma_hat(i,j) = 2*params[1]*
	      (G_u_np1(i,j) + G_u_np1(j,i))/2 
	      - 2*params[1]*
	      (G_u_n(i,j) + G_u_n(j,i))/2 + 
	      Sigma_n[ii*qdim_sigma*qdim_sigma + j*qdim_sigma +i];
	    if(i==j) 
	      sigma_hat(i,i) += ltrace_eps_np1 
		- ltrace_eps_n;
	  }
	}
	//cout<<"sigma_hat(2,0) = "<<sigma_hat(2,0)<<endl;
	GMM_TRACE2("End of computing sigma_hat")
      
	base_matrix proj;

	//#if 0
	size_type itdof2 = 0;
	
	base_matrix sigma_temp(qdim, qdim);
	for(size_type i = 0; i<R; ++i,itdof2 += qdim){
	  size_type incj = 0;
	  for(size_type j = itdof2; j<itdof2+qdim; ++j, ++incj){
	    size_type inck = 0;
	    for(size_type k = itdof2; k<itdof2+qdim; ++k, ++inck){
	      sigma_temp(incj,inck) = sigma_hat(j,k);
	    }
	  }

	  t_proj.do_projection(sigma_temp, params[2], 
			       proj, flag_proj);
      
	  GMM_TRACE2("End of computing projection")
	    for(dim_type j = 0; i<qdim_sigma; ++i){
	      for(dim_type k = 0; j<qdim_sigma; ++j){
		Sigma_np1[ii*qdim_sigma*qdim_sigma + k*qdim_sigma +j] = proj(j,k);
	      }
	    }
	  GMM_TRACE2("End of computing Sigma_np1")
	    }
	//#endif
	previous_cv = cv;

      }
  
      t.adjust_sizes(sizes_);
      
      /* interpolation of proj(sigma_hat) on gauss points */
      //   pf->interpolation(ctx, Sigma_np1, 
      //	   coeff_precalc, gmm::sqr(mf_sigma.get_qdim())); 


      /* copy the result into the tensor returned t */
      // std::copy(coeff_precalc.begin(), 
      //	coeff_precalc.end(), t.begin());
      std::copy(Sigma_np1.begin(), 
		Sigma_np1.end(), t.begin());


      //==================================================
      //==================================================

      GMM_TRACE2("End of computing the plasticity");

    }

};

  /** 
     Right hand side vector for plasticity 
      @ingroup asm
  */
  template<typename VECT> 
  void asm_plasticity_rhs (VECT &V, 
			   const mesh_im &mim, 
			   const mesh_fem &mf_u, 
			   const mesh_fem &mf_sigma, 
			   const mesh_fem &mf_data, 
			   const VECT &u_n,
			   const VECT &u_np1, 
			   const VECT &sigma_n, 
			   const VECT &sigma_np1, 
			   const VECT &lambda, 
			   const VECT &mu, 
			   const VECT &threshold, 
        const abstract_constraints_projection  &t_proj, 
			   const size_type flag_proj,
	const mesh_region &rg = mesh_region::all_convexes()) {



    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    GMM_TRACE2("Assembling the plasticity rhs");

    plasticity_nonlinear_term plast(mim, mf_u, mf_sigma,
				    &mf_data, u_n, u_np1,
				    sigma_n, sigma_np1, 
				    threshold, lambda, mu, 
				    t_proj, flag_proj);

    GMM_TRACE2("Assembling the plasticity rhs : end of building plast");



    generic_assembly assem("t=comp(NonLin(#2).vGrad(#1));"
			   "e=(t{:,:,:,4,5}+t{:,:,:,5,4})/2;"
			   "V(#1) += e(i,j,:,i,j)");


    GMM_TRACE2("Assembling the plasticity rhs : end of assembling");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_sigma);
    if (&(mf_data)!=NULL)
      assem.push_mf(mf_data);
    assem.push_nonlinear_term(&plast);
    assem.push_vec(V);
    assem.assembly(rg);



    GMM_TRACE2("End of assembling the plasticity rhs");
  }
  

  /** 
      Tangent matrix for plasticity
      @ingroup asm
  */
  template<typename MAT,typename VECT> 
  void asm_plasticity_tangent_matrix(MAT &H, 
				     const mesh_im &mim, 
				     const mesh_fem &mf_u, 
				     const mesh_fem &mf_sigma,
				     const mesh_fem &mf_data, 
				     const VECT &u_n,
				     const VECT &u_np1, 
				     const VECT &sigma_n, 
				     const VECT &sigma_np1, 
				     const VECT &lambda, 
				     const VECT &mu, 
				     const VECT &threshold, 
        const abstract_constraints_projection &t_proj, 
				    const size_type flag_proj,
        const mesh_region &rg = mesh_region::all_convexes()) {



    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    GMM_TRACE2("Assembling the plasticity tangent matrix");


    plasticity_nonlinear_term gradplast(mim, mf_u, mf_sigma,
					&mf_data, u_n, u_np1,
					sigma_n, sigma_np1, 
					threshold, lambda, mu,
					t_proj, flag_proj);

    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(NonLin(#2).vGrad(#1).vGrad(#1).Base(#2))(i,j,:,:,:,:,:,:,i,j,:);"
			   "M(#1,#1)+=  sym(t(k,l,:,l,k,:,m).mu(m)+t(k,l,:,k,l,:,m).mu(m)+t(k,k,:,l,l,:,m).lambda(m))");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_sigma);
    if (&(mf_data)!=NULL)
      assem.push_mf(mf_data);;
    assem.push_data(lambda);
    assem.push_data(mu);
    assem.push_nonlinear_term(&gradplast);
    assem.push_mat(H);
    assem.assembly(rg);

    GMM_TRACE2("End of assembling the plasticity tangent matrix");

  }



  //=========================================================================
  //
  //  Plasticity Brick
  //
  //=========================================================================

  struct plasticity_brick : public virtual_brick {

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
                                build_version version) const {


      GMM_ASSERT1(mims.size() == 1,
		  "Plasticity brick need a single mesh_im");
      GMM_ASSERT1(vl.size() == 1,
		  "Plasticity brick need one variable"); /** vl[0] = u */
      GMM_ASSERT1(dl.size() == 4,
		  "Wrong number of data for plasticity brick, "
                  << dl.size() << " should be 4.");
      GMM_ASSERT1(matl.size() == 1,  "Wrong number of terms for "
		  "plasticity brick");


      GMM_TRACE2("Assembling the tangent terms");



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


      // ...

      const model_real_plain_vector &sigma_np1 = 
	md.real_variable(dl[3], 0);
      const model_real_plain_vector &sigma_n = 
	md.real_variable(dl[3], 1);
      const mesh_fem &mf_sigma = 
	*(md.pmesh_fem_of_variable(dl[3]));
	
      const mesh_im &mim = *mims[0];

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	asm_plasticity_tangent_matrix
  		(matl[0], mim, mf_u, mf_sigma, *mf_data, u_n,
  		 u_np1, sigma_n, sigma_np1, 
		 lambda, mu, threshold, t_proj, 1);
      }

      if (version & model::BUILD_RHS) {
	asm_plasticity_rhs
  		(vecl[0], mim, mf_u, mf_sigma, *mf_data, u_n,
  		 u_np1, sigma_n, sigma_np1, 
  		 lambda, mu, threshold, t_proj, 0);
	gmm::scale(vecl[0], scalar_type(-1));
      }

      GMM_TRACE2("End of assembling the tangent terms");
    }



    plasticity_brick(const abstract_constraints_projection &t_proj_):t_proj(t_proj_){

      set_flags("Plasticity brick", false /* is linear*/,
            true /* is symmetric */, false /* is coercive */,
	    true /* is real */, false /* is complex */);
      GMM_TRACE2("End of creating the plasticity brick");
    }

  };
  
  //=========================================================================
  //  Add a plasticity brick
  //=========================================================================

  size_type add_plasticity_brick (model &md, 
		            const mesh_im &mim, 
			    const std::string &varname,
			    const std::string &datalambda,
			    const std::string &datamu,
			    const std::string &datathreshold, 
			    const std::string &datasigma,
			    size_type region) {

    VM_projection proj(0);
    pbrick pbr = new plasticity_brick(proj);

    model::termlist tl;
    tl.push_back(model::term_description
		 (varname, varname, true));
    model::varnamelist dl(1, datalambda);
    dl.push_back(datamu); 
    dl.push_back(datathreshold); 
    dl.push_back(datasigma);
    model::varnamelist vl(1, varname);
	

    GMM_TRACE2("End of adding the plasticity brick");

    return md.add_brick(pbr, vl, dl, tl, 
			model::mimlist(1,&mim), region);
  }

}  /* end of namespace getfem.  */

