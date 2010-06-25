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
    base_vector coeff_precalc;
    size_type N, previous_cv;
    const mesh_im &mim;
    const mesh_fem &mf_u;
    const mesh_fem &mf_sigma;
    const mesh_fem *mf_data;
    const std::vector<scalar_type> &Sigma_n; 
    std::vector<scalar_type> *Sigma_np1;
    std::vector<scalar_type> U_n;
    std::vector<scalar_type> U_np1;
    std::vector<scalar_type> threshold, lambda, mu;  
    bgeot::multi_index sizes_;
    const abstract_constraints_projection &t_proj;
    fem_precomp_pool fppool;
    std::vector<scalar_type> stored_proj;
 
    const size_type flag_proj;
    bool write_sigma_np1;

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
		  std::vector<scalar_type> *Sigma_np1_, 
	     	  const std::vector<scalar_type> &threshold_, 
	      	  const std::vector<scalar_type> &lambda_,
	       	  const std::vector<scalar_type> &mu_, 
	       	  const abstract_constraints_projection  &t_proj_,
		  const size_type flag_proj_, bool write_sigma_np1_)
      : mim(mim_), mf_u(mf_u_), mf_sigma(mf_sigma_), Sigma_n(Sigma_n_),
	Sigma_np1(Sigma_np1_), t_proj(t_proj_), flag_proj(flag_proj_),
	write_sigma_np1(write_sigma_np1_) {
      GMM_TRACE2("Building the plasticity non linear term");
      //      cout << "t_proj flag_hyp debut nonlinear term = " << t_proj.flag_hyp << endl;

      params = base_vector(3); 
      N = mf_u_.linked_mesh().dim();
      coeff_precalc = base_vector(N*N*N*N);
      gmm::resize(U_n, mf_u_.nb_basic_dof());
      gmm::resize(U_np1, mf_u_.nb_basic_dof());
      sizes_ = bgeot::multi_index(N, N, N, N);
      mf_u.extend_vector(gmm::sub_vector
			 (U_n_, gmm::sub_interval(0,mf_u_.nb_dof())), U_n);
      mf_u.extend_vector(gmm::sub_vector
			 (U_np1_, gmm::sub_interval(0,mf_u_.nb_dof())), U_np1);


      if (mf_data_ != NULL) {
	gmm::resize(mu, mf_data_->nb_basic_dof());
	gmm::resize(lambda, mf_data_->nb_basic_dof());
	gmm::resize(threshold, mf_data_->nb_basic_dof());
	
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

      GMM_ASSERT1(mf_u.get_qdim() == N, "wrong qdim for the mesh_fem");      
      if (flag_proj==0) sizes_.resize(2);

      previous_cv = size_type(-1);
      GMM_TRACE2("End of building the plasticity non linear term");
    }


    const bgeot::multi_index &sizes() const { return sizes_; }



    // compute() method from nonlinear_elem, gives on output the tensor
    virtual void compute(fem_interpolation_context& ctx,
			 bgeot::base_tensor &t){ 

      GMM_TRACE2("Computing the plasticity");
      // cout << "t_proj flag_hyp debut compute = " << t_proj.flag_hyp << endl;


      size_type cv = ctx.convex_num();//index of current element
      size_type qdim = mf_u.get_qdim();
      size_type qdim_sigma = mf_sigma.get_qdim();
      // size_type qqdim_sigma = gmm::vect_size(Sigma_n) / mf_sigma.nb_dof();
      pfem pf_sigma = ctx.pf();
      size_type nbd_sigma = pf_sigma->nb_dof(cv);
      GMM_ASSERT1(pf_sigma->is_lagrange(),
		  "Sorry, works only for Lagrange fems");
      size_type size_proj = qdim * qdim * (flag_proj == 1 ? qdim * qdim : 1);

      if(previous_cv != cv) {
	GMM_TRACE2("previous");

	stored_proj.resize(nbd_sigma * size_proj);
	base_matrix G;
	bgeot::vectors_to_base_matrix
	  (G, mf_u.linked_mesh().points_of_convex(cv));
	bgeot::pgeometric_trans pgt=
	    mf_u.linked_mesh().trans_of_convex(cv);

	fem_interpolation_context ctx_data;
	pfem pf_data;
	base_vector coeff_data;

	if (mf_data) {

	  GMM_TRACE2("data");
	 
	  pf_data = mf_data->fem_of_element(cv); 
	  size_type nbd_data = pf_data->nb_dof(cv);
	  coeff_data.resize(nbd_data*3);

	  // Definition of the coeff of Lame

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
	  GMM_TRACE2("non data");
	  
	  params[0] = lambda[0];
	  params[1] = mu[0];
	  params[2] = threshold[0];
	}


	for (size_type ii = 0; ii < nbd_sigma; ++ii) {

	  if (mf_data) {
	    ctx_data.set_ii(ii);
	    pf_data->interpolation(ctx_data, coeff_data, params, 3);
	  } 

	  std::vector<scalar_type> coeff_u_n, coeff_u_np1;
	  
	  pfem pf_u = mf_u.fem_of_element(cv);
	  // size_type nbd_u = pf_u->nb_dof(cv); // = 6
	  
	  size_type cvnbdof_u = 
	    mf_u.nb_basic_dof_of_element(cv); // = 12
	  
	  
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
	  
	  
	  // dal::bit_vector dof_sigma_done;
	  //  dof_sigma_done.sup(0, mf_sigma.nb_basic_dof());
	  // itdof = mf_sigma.ind_basic_dof_of_element(cv).begin();//taille = 12 = 3noeuds*4(qdim_sigma)
	  
	   
	  ctx_u.set_ii(ii);
	  
	  pf_u->interpolation_grad(ctx_u, coeff_u_n, G_u_n, dim_type(qdim));
	  pf_u->interpolation_grad(ctx_u, coeff_u_np1, G_u_np1,
				   dim_type(qdim));
	  
	  GMM_TRACE2("non previous");
	  
	  // Compute sigma_hat = D*esp_np1 - D*eps_n + sigma_n
	  base_matrix sigma_hat(qdim, qdim);

	
	  // Compute lambda*tr(esp_n) and lambda*tr(esp_np1)
	  scalar_type ltrace_eps_n =  params[0]*gmm::mat_trace(G_u_n);
	  scalar_type ltrace_eps_np1 = params[0]*gmm::mat_trace(G_u_np1);
	  
	  // cout<<"col = "<<gmm::mat_ncols(G_u_n)<<", lig = "<<gmm::mat_nrows(G_u_n)<<", R = "<<R<<", nbd_sigma = "<<nbd_sigma<<endl;
	  
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
	  
	  GMM_TRACE2("End of computing sigma_hat");
	
	  base_matrix proj;
	  // cout<<*t_proj<<endl;

	  VM_projection pproj(0);

	  // cout << "Adress of t_proj " << &t_proj << endl;
	  // cout << "t_proj flag_hyp = " << t_proj.flag_hyp << endl;
	  // cout << "t_proj.do_projection " << &(t_proj.do_projection) << endl;

	  t_proj.do_projection(sigma_hat,params[2],proj,flag_proj);

	  cout<<"here ok"<<endl;
	   std::copy(proj.begin(), proj.end(),
	      stored_proj.begin() + proj.size() * ii);
	
	  if (flag_proj == 0 && write_sigma_np1) {
	  
	    GMM_TRACE2("End of computing projection");
	    for(dim_type i = 0; i < qdim; ++i){
	      for(dim_type j = 0; j < qdim; ++j){
		(*Sigma_np1)[idof_sigma + j*qdim_sigma + i] = proj(i,j);
	      }
	    }
	    GMM_TRACE2("End of computing Sigma_np1");
	  }
	  
	  
	}
	previous_cv = cv;
      }
      coeff_precalc.resize(size_proj);
      pf_sigma->interpolation(ctx, stored_proj, coeff_precalc,
			      dim_type(size_proj));
      
      t.adjust_sizes(sizes_);
      
      /*  copy the result into the tensor returned t */
       std::copy(coeff_precalc.begin(), coeff_precalc.end(), t.begin());


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
			   VECT *sigma_np1, 
			   const VECT &lambda, 
			   const VECT &mu, 
			   const VECT &threshold, 
        const abstract_constraints_projection  &t_proj,
			   bool write_sigma_np1,
	const mesh_region &rg = mesh_region::all_convexes()) {


    // cout << "t_proj flag_hyp debut asm_rhs = " << t_proj.flag_hyp << endl;

    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    GMM_TRACE2("Assembling the plasticity rhs");
    
    if(&t_proj == NULL)
	cout<<"pb ds plast asm rhs $$$$$$$$$$$$$$$$$"<<endl;

    plasticity_nonlinear_term plast(mim, mf_u, mf_sigma,
				    &mf_data, u_n, u_np1,
				    sigma_n, sigma_np1, 
				    threshold, lambda, mu, 
				    t_proj, 0, write_sigma_np1);

    GMM_TRACE2("Assembling the plasticity rhs : end of building plast");



//     generic_assembly assem("t=comp(NonLin(#2).vGrad(#1));"
// 			   "e=(t{:,:,:,4,5}+t{:,:,:,5,4})/2;"
// 			   "V(#1) += e(i,j,:,i,j)");

    generic_assembly assem("V(#1) + =comp(NonLin(#2).vGrad(#1))(i,j,:,i,j);");

    GMM_TRACE2("Assembling the plasticity rhs : end of assembling");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_sigma);
//  if (&(mf_data)!=NULL)
//       assem.push_mf(mf_data);
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
				     const VECT &lambda, 
				     const VECT &mu, 
				     const VECT &threshold, 
        const abstract_constraints_projection &t_proj,
        const mesh_region &rg = mesh_region::all_convexes()) {



    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    GMM_TRACE2("Assembling the plasticity tangent matrix");

    plasticity_nonlinear_term gradplast(mim, mf_u, mf_sigma,
					&mf_data, u_n, u_np1,
					sigma_n, 0, 
					threshold, lambda, mu,
					t_proj, 1, false);

    generic_assembly assem;

    if (&(mf_data)!=NULL) {
      assem.set("lambda=data$1(#2); mu=data$2(#2);"
		"t=comp(NonLin(#2).vGrad(#1).vGrad(#1).Base(#2))(i,j,:,:,:,:,:,:,i,j,:);"
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

      // cout << "t_proj flag_hyp debut asm_real = " << t_proj.flag_hyp << endl;


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

      // const model_real_plain_vector &sigma_np1 = 
      // md.real_variable(dl[3], 0);
      const model_real_plain_vector &sigma_n = 
	md.real_variable(dl[3], 1);
      const mesh_fem &mf_sigma = 
	*(md.pmesh_fem_of_variable(dl[3]));
      GMM_ASSERT1(!(mf_sigma.is_reduced()),
		  "Works only for pure Lagrange fems");
	
      const mesh_im &mim = *mims[0];
      
      if(&(t_proj) == NULL)
	cout<<"pb ds asm tangent term $$$$$$$$$$$$$$$$$"<<endl;

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	asm_plasticity_tangent_matrix
  		(matl[0], mim, mf_u, mf_sigma, *mf_data, u_n,
  		 u_np1, sigma_n, lambda, mu, threshold, t_proj, region);
      }

      if (version & model::BUILD_RHS) {
	asm_plasticity_rhs
  		(vecl[0], mim, mf_u, mf_sigma, *mf_data, u_n,
  		 u_np1, sigma_n, (model_real_plain_vector *)(0), 
  		 lambda, mu, threshold, t_proj, false, region);
	gmm::scale(vecl[0], scalar_type(-1));
      }

      GMM_TRACE2("End of assembling the tangent terms");
    }



    plasticity_brick(const abstract_constraints_projection &t_proj_)
      : t_proj(t_proj_){

      // cout << "t_proj_ flag_hyp brick = " << t_proj_.flag_hyp << endl;
      // cout << "t_proj flag_hyp brick = " << t_proj.flag_hyp << endl;

      set_flags("Plasticity brick", false /* is linear*/,
            true /* is symmetric */, false /* is coercive */,
	    true /* is real */, false /* is complex */);
      GMM_TRACE2("End of creating the plasticity brick");
    }

  };
  
  //=========================================================================
  //  Add a plasticity brick
  //=========================================================================

  size_type add_plasticity_brick(model &md, 
		            const mesh_im &mim, 
			    const std::string &varname,
			    const std::string &datalambda,
			    const std::string &datamu,
			    const std::string &datathreshold, 
			    const std::string &datasigma,
			    size_type region) {

    static VM_projection proj(0);

    // cout << "Adress of t_proj " << &proj << endl;
    // cout << "t_proj flag_hyp = " << proj.flag_hyp << endl;

    pbrick pbr = new plasticity_brick(proj);

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

}  /* end of namespace getfem.  */

