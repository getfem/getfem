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


#if 0



  /** Compute the projection of D*e + sigma_bar_ on a Gauss point. */
  class plasticity_projection : public nonlinear_elem_term {
  protected:
    base_vector params, coeff;
    size_type N, previous_cv;
    const mesh_im &mim;
    const mesh_fem &mf;
    const mesh_fem &mf_sigma;
    const mesh_fem &mf_data;
    std::vector<scalar_type> U;
    std::vector<scalar_type> stress_threshold;
    std::vector<scalar_type> lambda, mu;  
    bgeot::multi_index sizes_;
    const abstract_constraints_projection  *t_proj; 
    std::vector<std::vector<scalar_type> > &sigma_bar_;

    // to save the projection
    std::vector<std::vector<scalar_type> > &saved_proj_;
 
    const size_type flag_proj;
    bool fill_sigma_bar;

  
  public:  

    std::vector<std::vector<scalar_type> > &sigma_bar() { return sigma_bar_; }
    scalar_type &sigma_bar(size_type cv, size_type ii, int i, int j)
    { return sigma_bar_[cv][ii*N*N + j*N + i]; }
    std::vector<std::vector<scalar_type> > &saved_proj()
    { return saved_proj_; }
    scalar_type &saved_proj(size_type cv, size_type ii, int i, int j)
    { return saved_proj_[cv][ii*N*N + j*N + i]; }

    // constructor
    plasticity_projection(const mesh_im &mim_,
			  const mesh_fem &mf_,
			  const mesh_fem &mf_sigma_,
			  const mesh_fem &mf_data_,
			  const std::vector<scalar_type> &U_, 
			  const std::vector<scalar_type> &stress_threshold_, 
			  const std::vector<scalar_type> &lambda_,
			  const std::vector<scalar_type> &mu_, 
			  const abstract_constraints_projection  *t_proj_,
			  std::vector<std::vector<scalar_type> > &sigma_bar__, 
			  std::vector<std::vector<scalar_type> > &saved_proj__,
			  const size_type flag_proj_,
			  const bool fill_sigma) :
      params(3), N(mf_.linked_mesh().dim()), mim(mim_),
      mf(mf_), mf_sigma(mf_sigma_), mf_data(mf_data_),
      U(mf_.nb_basic_dof()),  stress_threshold(mf_data_.nb_basic_dof()),
      lambda(mf_data_.nb_basic_dof()), mu(mf_data_.nb_basic_dof()),
      sizes_(N, N, N, N), t_proj(t_proj_),
      sigma_bar_( sigma_bar__),saved_proj_(saved_proj__),
      flag_proj(flag_proj_)  {
    
      mf.extend_vector
	(gmm::sub_vector(U_, gmm::sub_interval(0, mf_.nb_dof())), U);
      mf_data.extend_vector(stress_threshold_, stress_threshold);
      mf_data.extend_vector(lambda_, lambda);
      mf_data.extend_vector(mu_, mu);
  
      fill_sigma_bar = fill_sigma;   /* always false during resolution, */
      /*                      true when called from compute_constraints */

      GMM_ASSERT1(mf.get_qdim() == N, "wrong qdim for the mesh_fem");      
      if (flag_proj==0) sizes_.resize(2);
    
      sigma_bar_.resize(mf.linked_mesh().convex_index().last_true()+1);    
      saved_proj_.resize(mf.linked_mesh().convex_index().last_true()+1);




      previous_cv = size_type(-1);

    }

    const bgeot::multi_index &sizes() const { return sizes_; }

    // compute() method from nonlinear_elem, gives on output the tensor
    virtual void compute(fem_interpolation_context& ctx,
			 bgeot::base_tensor &t){ 

      size_type cv = ctx.convex_num();

      if (cv != previous_cv) {









	previous_cv = cv;
      }

      t.adjust_sizes(sizes_);
      pf->interpolation(ctx, coeff_precalculés, t, gmm::sqr(mf.get_qdim()));


      /* 

      size_type ii = ctx.ii(); 

      pfem pf = ctx.pf();
      coeff.resize(mf.nb_basic_dof_of_element(cv));
      base_matrix gradU(N, N), sigma(N,N);
      gmm::copy(gmm::sub_vector
		(U, gmm::sub_index(mf.ind_basic_dof_of_element(cv))), coeff);
      pf->interpolation(ctx, coeff, gradU, mf.get_qdim());
      
      scalar_type ltrace_eps;
      ltrace_eps = params[0]*gmm::mat_trace(gradU);

      // if needed, we give sigma_bar[cv] and saved_proj[cv] a size equal
      // to the number of integration points on the convexe. Seems that
      // this is rarely needed. 
      if (sigma_bar_[cv].size() == 0){
	size_type nb_interp_pt = mim.int_method_of_element(cv)
	  ->approx_method()->nb_points_on_convex();
	sigma_bar_[cv].resize(N*N*nb_interp_pt);
	gmm::clear(sigma_bar_[cv]);
	saved_proj_[cv].resize(N*N*nb_interp_pt);
	gmm::clear(saved_proj_[cv]);
      }

      t.adjust_sizes(sizes_);
   
      for (dim_type i=0; i < N; ++i) {
	for (dim_type j=0; j < N; ++j) {
	  sigma(i,j) = 2*params[1]*(gradU(i,j)+gradU(j,i))/2.
	    + sigma_bar(cv,ii,i,j);
	  if(i==j) sigma(i,i) += ltrace_eps;
	}
      }
    
      base_matrix tau_star(N,N), gradproj(N,N), proj;
      t_proj->do_projection(sigma, params[2], proj, flag_proj);

      // we fill sigma_bar only when called from compute_constraints
      // (ie, when fill_sigma_bar is set)
      if (fill_sigma_bar && flag_proj==0) {

	for (dim_type i=0; i < N; ++i)
	  for (dim_type j=0; j < N; ++j)
	    saved_proj(cv,ii,i,j) = proj(i,j);

	
	gmm::add(gmm::scaled(gradU, -params[1]), proj);
	gmm::add(gmm::scaled(gmm::transposed(gradU), -params[1]), proj);

	for (dim_type i=0; i < N; ++i) {
	  proj(i,i) += ltrace_eps;
	  for (dim_type j=0; j < N; ++j)
	    sigma_bar(cv,ii,i,j) = proj(i,j);
	}
      }
      std::copy(proj.begin(),proj.end(), t.begin());


      */

    }

  };

  /** 
     Right hand side vector for plasticity 
      @ingroup asm
  */
  template<typename VECT> 
  void asm_rhs_for_plasticity
  (VECT &V, const mesh_im &mim, const mesh_fem &mf, const VECT &u_n,
   const VECT &u_np1, const mesh_fem &mfdata, const VECT &lambda,
   const VECT &mu, const mesh_fem &mf_sigma,
   const VECT &threshold, const VECT &sigma_n, const VECT &sigma_np1,
   const mesh_region &rg = mesh_region::all_convexes()) {

    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    plasticity_projection plast(u_n, u_np1, mfdata, lambda, mu,
				threshold, mf_sigma, sigma_n, sigma_np1);

    generic_assembly assem("sigma=data(#1);"
			   "t=comp(NonLin(#2).vGrad(#1));"
			   "e=(t{:,:,:,4,5}+t{:,:,:,5,4})/2;"
			   "V(#1) += e(i,j,:,i,j)");
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_sigma);
    assem.push_mf(mfdata);
    assem.push_data(SIGMA);
    assem.push_nonlinear_term(plast);
    assem.push_vec(V);
    assem.assembly(rg);
  }
  
  /** 
      Tangent matrix for plasticity
      @ingroup asm
  */
  template<typename MAT,typename VECT> 
  void asm_lhs_for_plasticity
  (MAT &H, const mesh_im &mim, const mesh_fem &mf, const mesh_fem &mf_sigma, const mesh_fem &mfdata,
   const VECT &LAMBDA, const VECT &MU, nonlinear_elem_term *gradplast,
   const mesh_region &rg = mesh_region::all_convexes()) {
    GMM_ASSERT1(mf.get_qdim() == mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");
    /*generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
				   "t=comp(NonLin(#1,#2).vGrad(#1).vGrad(#1).Base(#2));"
				   "e=(t{:,:,:,:,:,6,7,:,9,10,:}+t{:,:,:,:,:,7,6,:,9,10,:}+t{:,:,:,:,:,6,7,:,10,9,:}+t{:,:,:,:,:,7,6,:,10,9,:})/4;"
				   "M(#1,#1)+=  sym(2*e(i,j,k,l,:,k,l,:,i,j,m).mu(m)+e(i,j,k,k,:,l,l,:,i,j,m).lambda(m))");
    */
    generic_assembly assem("lambda=data$1(#2); mu=data$2(#2);"
			   "t=comp(NonLin(#1,#2).vGrad(#1).vGrad(#1).Base(#2))(i,j,:,:,:,:,:,:,i,j,:);"
			   //"t=comp(NonLin(#1,#2)(i,j,:,:).vGrad(#1).vGrad(#1)(:,i,j).Base(#2));"
			   "M(#1,#1)+=  sym(t(k,l,:,l,k,:,m).mu(m)+t(k,l,:,k,l,:,m).mu(m)+t(k,k,:,l,l,:,m).lambda(m))");
			   /*"a=comp(NonLin(#1,#2)(i,j,k,l).vGrad(#1)(:,k,l).vGrad(#1)(:,i,j).Base(#2)(:))"
			   " +comp(NonLin(#1,#2)(i,j,k,l).vGrad(#1)(:,l,k).vGrad(#1)(:,i,j).Base(#2)(:));"
			   "b=comp(NonLin(#1,#2)(i,j,k,k).vGrad(#1)(:,l,l).vGrad(#1)(:,i,j).Base(#2)(:));"
			   "M(#1,#1)+=  sym(a(:,:,m).mu(m) + b(:,:,m).lambda(m));");*/
    // comp()  to be optimized !!
    assem.push_mi(mim);
    assem.push_mf(mf);
    assem.push_mf(mf_sigma);
    assem.push_mf(mfdata);
    assem.push_data(LAMBDA);
    assem.push_data(MU);
    assem.push_nonlinear_term(gradplast);
    assem.push_mat(H);
    assem.assembly(rg);
  }








#endif



  //=========================================================================
  //
  //  Plasticity Brick
  //
  //=========================================================================

  struct plasticity_brick : public virtual_brick {



    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
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

      const model_real_plain_vector &u_np1 = md.real_variable(vl[0], 0);
      const model_real_plain_vector &u_n = md.real_variable(vl[0], 1);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));

      const model_real_plain_vector &lambda = md.real_variable(dl[0]);
      const mesh_fem *mf_lambda = md.pmesh_fem_of_variable(dl[0]);

      // ...

      const model_real_plain_vector &sigma_np1 = md.real_variable(dl[3], 0);
      const model_real_plain_vector &sigma_n = md.real_variable(dl[3], 1);
      const mesh_fem &mf_sigma = *(md.pmesh_fem_of_variable(dl[1]));
      

  /*    const mesh_fem &mf_params =*(md.pmesh_fem_of_variable(dl[0]));
      const model_real_plain_vector &lambda = md.real_variable(dl[0]);
      const model_real_plain_vector &mu = md.real_variable(dl[1]);
      const model_real_plain_vector &stress_threshold = md.real_variable(dl[2]);
      const model_real_plain_vector &sigma = md.real_variable(dl[3]);	
      const mesh_im &mim = *mims[0];
*/

#if 0

	size-type N = mf_sigma.linked_mesh().dim());
      std::vector<std::vector<scalar_type> > sigma_bar ;
	for(dim_type i = 0; i<N; ++i)
	sigma_bar[mf_sigma.convex_index()][] = 0;
 
  //  std::vector<std::vector<scalar_type> > &saved_proj ;
	size_type flag_hyp;
     //   const abstract_constraints_projection *t_proj;
    //  mesh_region rg(region);
    //  mf_u.linked_mesh().intersect_with_mpi_region(rg);

    
	
      

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	GMM_TRACE2("Plasticity stiffness matrix assembly");

	
        flag_hyp = 1;
	//VM_projection gradproj(flag_hyp);
	//plasticity_projection gradproj(mim, mf_u, mf_params, u,stress_threshold, lambda, mu, t_proj, sigma_bar,  saved_proj, flag_hyp , false);
	//asm_lhs_for_plasticity(matl[0], mim, mf_u, mf_params, lambda, mu, &gradproj );
      }


      if (version & model::BUILD_RHS) {
        flag_hyp = 0;
      //  VM_projection proj(flag_hyp);
	//plasticity_projection proj(mim, mf_u, mf_params, u,stress_threshold, lambda, mu, t_proj, sigma_bar,  saved_proj, flag_hyp, false);
	//asm_rhs_for_plasticity(vecl[0], mim, mf_u, mf_params, &proj );
	gmm::scale(vecl[0], scalar_type(-1));
      }

#endif

    }



    plasticity_brick(void){
      set_flags("Plasticity brick", false /* is linear*/,
                true /* is symmetric */, false /* is coercive */,
		true /* is real */, false /* is complex */);
    }

  };
  
  //=========================================================================
  //  Add a plasticity brick
  //=========================================================================

  size_type add_plasticity_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &datalambda, const std::string &datamu,
   const std::string &datathreshold, const std::string &datasigma,
   size_type region) {
    pbrick pbr = new plasticity_brick();

    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, datalambda);
    dl.push_back(datamu); dl.push_back(datathreshold); dl.push_back(datasigma);
    model::varnamelist vl(1, varname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1,&mim), region);
  }

}  /* end of namespace getfem.  */

