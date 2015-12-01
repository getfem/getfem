/*===========================================================================

 Copyright (C) 2000-2015 Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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
#include "getfem/getfem_generic_assembly.h"
#include "gmm/gmm_dense_matrix_functions.h"


namespace getfem {

  enum elastoplasticity_nonlinear_term_version { PROJ,
                                                 GRADPROJ,
                                                 PLAST
  };

  /** Compute the projection of D*e + sigma_bar_
      on the dof of sigma. */
  class elastoplasticity_nonlinear_term : public nonlinear_elem_term {

  protected:
    const mesh_im &mim;
    const mesh_fem &mf_u;
    const mesh_fem &mf_sigma;
    const mesh_fem *pmf_data;
    model_real_plain_vector U_n,U_np1;
    model_real_plain_vector Sigma_n;
    model_real_plain_vector threshold, lambda, mu;
    const abstract_constraints_projection &t_proj;
    const size_type option;
    const size_type flag_proj;
    const bool store_sigma;

    bgeot::multi_index sizes_;

    size_type N, size_proj;

    // temporary variables
    base_vector params;
    size_type current_cv;
    model_real_plain_vector convex_coeffs, interpolated_val;

    // storage variables
    model_real_plain_vector cumulated_sigma; // either the projected stress (option==PROJ)
                                             // or the plastic stress (option==PLAST)
    model_real_plain_vector cumulated_count;

    fem_precomp_pool fppool;


    // computes stresses or stress projections on all sigma dofs of a convex
    void compute_convex_coeffs(size_type cv) {

      current_cv = cv;

      pfem pf_sigma = mf_sigma.fem_of_element(cv);
      size_type nbd_sigma = pf_sigma->nb_dof(cv);
      size_type qdim_sigma = mf_sigma.get_qdim();

      gmm::resize(convex_coeffs, size_proj*nbd_sigma);

      base_matrix G;
      bgeot::vectors_to_base_matrix
        (G, mf_u.linked_mesh().points_of_convex(cv));
      bgeot::pgeometric_trans pgt =
        mf_u.linked_mesh().trans_of_convex(cv);

      // if the Lame coefficient are vector fields
      base_vector coeff_data;
      pfem pf_data;
      fem_interpolation_context ctx_data;
      if (pmf_data) {
        pf_data = pmf_data->fem_of_element(cv);
        size_type nbd_data = pf_data->nb_dof(cv);
        coeff_data.resize(nbd_data*3);

        // Definition of the Lame coeff
        mesh_fem::ind_dof_ct::const_iterator itdof
          = pmf_data->ind_basic_dof_of_element(cv).begin();
        for (size_type k = 0; k < nbd_data; ++k, ++itdof) {
          coeff_data[k*3] = lambda[*itdof];
          coeff_data[k*3+1] = mu[*itdof];
          coeff_data[k*3+2] = threshold[*itdof];
        }
        GMM_ASSERT1(pf_data->target_dim() == 1,
                    "won't interpolate on a vector FEM... ");

        pfem_precomp pfp_data = fppool(pf_data, pf_sigma->node_tab(cv));
        ctx_data = fem_interpolation_context
          (pgt, pfp_data, size_type(-1), G, cv, short_type(-1));
      }

      // Definition of the coeff for du = u_n-u_np1 and optionally for u_np1
      size_type cvnbdof_u = mf_u.nb_basic_dof_of_element(cv);
      model_real_plain_vector coeff_du(cvnbdof_u);
      model_real_plain_vector coeff_u_np1(cvnbdof_u);
      mesh_fem::ind_dof_ct::const_iterator itdof
        = mf_u.ind_basic_dof_of_element(cv).begin();
      for (size_type k = 0; k < cvnbdof_u; ++k, ++itdof) {
        coeff_du[k] = U_np1[*itdof] - U_n[*itdof];
        coeff_u_np1[k] = U_np1[*itdof];
      }

      pfem pf_u = mf_u.fem_of_element(cv);
      pfem_precomp pfp_u = fppool(pf_u, pf_sigma->node_tab(cv));
      fem_interpolation_context
        ctx_u(pgt, pfp_u, size_type(-1), G, cv, short_type(-1));

      size_type qdim = mf_u.get_qdim();
      base_matrix G_du(qdim, qdim), G_u_np1(qdim, qdim); // G_du = G_u_np1 - G_u_n

      for (size_type ii = 0; ii < nbd_sigma; ++ii) {

        if (pmf_data) {
          // interpolation of the data on sigma dof
          ctx_data.set_ii(ii);
          pf_data->interpolation(ctx_data, coeff_data, params, 3);
        }

        // interpolation of the gradient of du and u_np1 on sigma dof
        ctx_u.set_ii(ii);
        pf_u->interpolation_grad(ctx_u, coeff_du, G_du, dim_type(qdim));
        if (option == PLAST)
          pf_u->interpolation_grad(ctx_u, coeff_u_np1, G_u_np1, dim_type(qdim));

        // Compute lambda*(tr(eps_np1)-tr(eps_n)) and lambda*tr(eps_np1)
        scalar_type ltrace_deps = params[0]*gmm::mat_trace(G_du);
        scalar_type ltrace_eps_np1 = (option == PLAST) ?
                                     params[0]*gmm::mat_trace(G_u_np1) : 0.;

        // Compute sigma_hat = D*(eps_np1 - eps_n) + sigma_n
        // where D represents the elastic stiffness tensor
        base_matrix sigma_hat(qdim, qdim);
        size_type sigma_dof = mf_sigma.ind_basic_dof_of_element(cv)[ii*qdim_sigma];
        for (dim_type j = 0; j < qdim; ++j) {
          for (dim_type i = 0; i < qdim; ++i)
            sigma_hat(i,j) = Sigma_n[sigma_dof++]
                             + params[1]*(G_du(i,j) + G_du(j,i));
          sigma_hat(j,j) += ltrace_deps;
        }

        // Compute the projection or its grad
        base_matrix proj;
        t_proj.do_projection(sigma_hat, params[2], proj, flag_proj);

        // Compute the plastic part if required
        if (option == PLAST)
          for (dim_type i = 0; i < qdim; ++i) {
            for (dim_type j = 0; j < qdim; ++j)
              proj(i,j) -= params[1]*(G_u_np1(i,j) + G_u_np1(j,i));
            proj(i,i) -= ltrace_eps_np1;
          }

        // Fill in convex_coeffs with sigma or its grad
        std::copy(proj.begin(), proj.end(),
                  convex_coeffs.begin() + proj.size() * ii);

        // Store the projected or plastic sigma
        if (store_sigma) {
          sigma_dof = mf_sigma.ind_basic_dof_of_element(cv)[ii*qdim_sigma];
          for (dim_type j = 0; j < qdim; ++j) {
            for (dim_type i = 0; i < qdim; ++i) {
              cumulated_count[sigma_dof] += 1;
              cumulated_sigma[sigma_dof++] += proj(i,j);
            }
          }
        }

      } // ii = 0:nbd_sigma-1

    }

  public:

    // constructor
    elastoplasticity_nonlinear_term
      (const mesh_im &mim_,
       const mesh_fem &mf_u_,
       const mesh_fem &mf_sigma_,
       const mesh_fem *pmf_data_,
       const model_real_plain_vector &U_n_,
       const model_real_plain_vector &U_np1_,
       const model_real_plain_vector &Sigma_n_,
       const model_real_plain_vector &threshold_,
       const model_real_plain_vector &lambda_,
       const model_real_plain_vector &mu_,
       const abstract_constraints_projection  &t_proj_,
       size_type option_,
       bool store_sigma_) :
      mim(mim_), mf_u(mf_u_), mf_sigma(mf_sigma_), pmf_data(pmf_data_),
      Sigma_n(Sigma_n_), t_proj(t_proj_), option(option_),
      flag_proj(option == GRADPROJ ? 1 : 0),
      store_sigma(option == GRADPROJ ? false : store_sigma_) {

      params.resize(3);
      N = mf_u.linked_mesh().dim();

      sizes_ = (flag_proj == 0 ? bgeot::multi_index(N,N)
                               : bgeot::multi_index(N,N,N,N));

      // size_proj is different if we compute the projection
      // or the gradient of the projection
      size_proj = (flag_proj == 0 ? N*N : N*N*N*N);

      gmm::resize(U_n, mf_u.nb_basic_dof());
      gmm::resize(U_np1, mf_u.nb_basic_dof());
      gmm::resize(Sigma_n, mf_sigma.nb_basic_dof());
      mf_u.extend_vector(gmm::sub_vector(U_n_,
                                         gmm::sub_interval(0,mf_u.nb_dof())),
                         U_n);
      mf_u.extend_vector(gmm::sub_vector(U_np1_,
                                         gmm::sub_interval(0,mf_u.nb_dof())),
                         U_np1);
      mf_sigma.extend_vector(gmm::sub_vector(Sigma_n_,
                                             gmm::sub_interval(0,mf_sigma.nb_dof())),
                             Sigma_n);

      if (pmf_data != NULL) {
        gmm::resize(mu, pmf_data->nb_basic_dof());
        gmm::resize(lambda, pmf_data->nb_basic_dof());
        gmm::resize(threshold, pmf_data->nb_basic_dof());
        pmf_data->extend_vector(threshold_, threshold);
        pmf_data->extend_vector(lambda_, lambda);
        pmf_data->extend_vector(mu_, mu);
      } else {
        gmm::resize(mu, 1); mu[0]  =  mu_[0];
        gmm::resize(lambda, 1); lambda[0]  =  lambda_[0];
        gmm::resize(threshold, 1); threshold[0] =  threshold_[0];
        params[0] = lambda[0];
        params[1] = mu[0];
        params[2] = threshold[0];
      }
      GMM_ASSERT1(mf_u.get_qdim() == N,
                  "wrong qdim for the mesh_fem");

      gmm::resize(interpolated_val, size_proj);

      if (store_sigma) {
        cumulated_sigma.resize(mf_sigma.nb_dof());
        cumulated_count.resize(mf_sigma.nb_dof());
      }

      // used to know if the current element is different
      // than the previous one and so if a new computation
      // is necessary or not.
      current_cv = size_type(-1);

    }


    const bgeot::multi_index &sizes(size_type) const { return sizes_; }


    // method from nonlinear_elem_term, gives on output the tensor
    virtual void compute(fem_interpolation_context& ctx,
                         bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num(); //index of current element
      pfem pf_sigma = ctx.pf();
      GMM_ASSERT1(pf_sigma->is_lagrange(),
                  "Sorry, works only for Lagrange fems");

      // if the current element is different than the previous one
      if (cv != current_cv)
        compute_convex_coeffs(cv);

      // interpolation of the sigma or its grad on sigma dof
      pf_sigma->interpolation(ctx, convex_coeffs, interpolated_val, dim_type(size_proj));

      // copy the result into the returned tensor t
      t.adjust_sizes(sizes_);
      std::copy(interpolated_val.begin(), interpolated_val.end(), t.begin());
    }


    // method to get the averaged sigma stored during the assembly
    void get_averaged_sigmas(model_real_plain_vector &sigma) {
       model_real_plain_vector glob_cumulated_count(mf_sigma.nb_dof());
       MPI_SUM_VECTOR(cumulated_sigma, sigma);
       MPI_SUM_VECTOR(cumulated_count, glob_cumulated_count);
       size_type imax = mf_sigma.nb_dof();
       for (size_type i = 0; i < imax; ++i)
         sigma[i] /= glob_cumulated_count[i];
    }

};



  /**
     Right hand side vector for elastoplasticity
      @ingroup asm
  */
  void asm_elastoplasticity_rhs
    (model_real_plain_vector &V,
     model_real_plain_vector *saved_sigma,
     const mesh_im &mim,
     const mesh_fem &mf_u,
     const mesh_fem &mf_sigma,
     const mesh_fem &mf_data,
     const model_real_plain_vector &u_n,
     const model_real_plain_vector &u_np1,
     const model_real_plain_vector &sigma_n,
     const model_real_plain_vector &lambda,
     const model_real_plain_vector &mu,
     const model_real_plain_vector &threshold,
     const abstract_constraints_projection  &t_proj,
     size_type option_sigma,
     const mesh_region &rg = mesh_region::all_convexes()) {

    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
                "wrong qdim for the mesh_fem");
    GMM_ASSERT1(option_sigma == PROJ || option_sigma == PLAST,
                "wrong option parameter");

    elastoplasticity_nonlinear_term plast(mim, mf_u, mf_sigma, &mf_data,
                                          u_n, u_np1, sigma_n,
                                          threshold, lambda, mu,
                                          t_proj, option_sigma, (saved_sigma != NULL));

    generic_assembly assem("V(#1) + =comp(NonLin(#2).vGrad(#1))(i,j,:,i,j);");

    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_mf(mf_sigma);
    assem.push_nonlinear_term(&plast);
    assem.push_vec(V);
    assem.assembly(rg);

    if (saved_sigma)
      plast.get_averaged_sigmas(*saved_sigma);
  }


  /**
      Tangent matrix for elastoplasticity
      @ingroup asm
  */
  void asm_elastoplasticity_tangent_matrix
    (model_real_sparse_matrix &H,
     const mesh_im &mim,
     const mesh_fem &mf_u,
     const mesh_fem &mf_sigma,
     const mesh_fem &mf_data,
     const model_real_plain_vector &u_n,
     const model_real_plain_vector &u_np1,
     const model_real_plain_vector &sigma_n,
     const model_real_plain_vector &lambda,
     const model_real_plain_vector &mu,
     const model_real_plain_vector &threshold,
     const abstract_constraints_projection &t_proj,
     const mesh_region &rg = mesh_region::all_convexes()) {

    GMM_ASSERT1(mf_u.get_qdim() == mf_u.linked_mesh().dim(),
                "wrong qdim for the mesh_fem");

    elastoplasticity_nonlinear_term gradplast(mim, mf_u, mf_sigma, &mf_data,
                                              u_n, u_np1, sigma_n,
                                              threshold, lambda, mu,
                                              t_proj, GRADPROJ, false);

    generic_assembly assem;

    if (&(mf_data)!=NULL)
      assem.set("lambda=data$1(#3); mu=data$2(#3);"
                "t=comp(NonLin(#2).vGrad(#1).vGrad(#1).Base(#3))(i,j,:,:,:,:,:,:,i,j,:);"
                "M(#1,#1)+=  sym(t(k,l,:,l,k,:,m).mu(m)+t(k,l,:,k,l,:,m).mu(m)+t(k,k,:,l,l,:,m).lambda(m))");
    else
      assem.set("lambda=data$1(1); mu=data$2(1);"
                "t=comp(NonLin(#2).vGrad(#1).vGrad(#1))(i,j,:,:,:,:,:,:,i,j);"
                "M(#1,#1)+= sym(t(k,l,:,l,k,:).mu(1)+t(k,l,:,k,l,:).mu(1)+t(k,k,:,l,l,:).lambda(1))");

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

    pconstraints_projection  t_proj;

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

      const model_real_plain_vector &u_np1 = md.real_variable(vl[0], 0);
      const model_real_plain_vector &u_n = md.real_variable(vl[0], 1);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));

      const model_real_plain_vector &lambda = md.real_variable(dl[0]);
      const model_real_plain_vector &mu = md.real_variable(dl[1]);
      const model_real_plain_vector &threshold = md.real_variable(dl[2]);
      const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);

      const model_real_plain_vector &sigma_n = md.real_variable(dl[3]);
      const mesh_fem &mf_sigma = *(md.pmesh_fem_of_variable(dl[3]));
      GMM_ASSERT1(!(mf_sigma.is_reduced()),
                  "Works only for pure Lagrange fems");

      const mesh_im &mim = *mims[0];
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
        gmm::clear(matl[0]);
        asm_elastoplasticity_tangent_matrix
          (matl[0], mim, mf_u, mf_sigma, *mf_data, u_n,
           u_np1, sigma_n, lambda, mu, threshold, *t_proj, rg);
      }

      if (version & model::BUILD_RHS) {
        model_real_plain_vector *dummy = 0;
        asm_elastoplasticity_rhs(vecl[0], dummy,
                                 mim, mf_u, mf_sigma, *mf_data,
                                 u_n, u_np1, sigma_n,
                                 lambda, mu, threshold, *t_proj, PROJ, rg);
        gmm::scale(vecl[0], scalar_type(-1));
      }

    }

    // constructor
    elastoplasticity_brick(const pconstraints_projection &t_proj_)
      : t_proj(t_proj_) {
      set_flags("Elastoplasticity brick", false /* is linear*/,
                true /* is symmetric */, false /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };


  //=================================================================
  //  Add a elastoplasticity brick
  //=================================================================

  size_type add_elastoplasticity_brick
    (model &md,
     const mesh_im &mim,
     const pconstraints_projection &ACP,
     const std::string &varname,
     const std::string &datalambda,
     const std::string &datamu,
     const std::string &datathreshold,
     const std::string &datasigma,
     size_type region) {

    pbrick pbr = std::make_shared<elastoplasticity_brick>(ACP);

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
                                  const pconstraints_projection &ACP,
                                  const std::string &datalambda,
                                  const std::string &datamu,
                                  const std::string &datathreshold,
                                  const std::string &datasigma) {

    const model_real_plain_vector &u_np1 = md.real_variable(varname, 0);
    model_real_plain_vector &u_n = md.set_real_variable(varname, 1);
    const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(varname));

    const model_real_plain_vector &lambda = md.real_variable(datalambda);
    const model_real_plain_vector &mu = md.real_variable(datamu);
    const model_real_plain_vector &threshold = md.real_variable(datathreshold);
    const mesh_fem *mf_data = md.pmesh_fem_of_variable(datalambda);

    const model_real_plain_vector &sigma_n = md.real_variable(datasigma);
    const mesh_fem &mf_sigma = *(md.pmesh_fem_of_variable(datasigma));

    // dim_type N = mf_sigma.linked_mesh().dim();

    mesh_region rg = mim.linked_mesh().get_mpi_region();

    model_real_plain_vector sigma_np1(mf_sigma.nb_dof());
    model_real_plain_vector dummyV(mf_u.nb_dof());
    asm_elastoplasticity_rhs(dummyV, &sigma_np1,
                             mim, mf_u, mf_sigma, *mf_data,
                             u_n, u_np1, sigma_n,
                             lambda, mu, threshold, *ACP, PROJ, rg);

    // upload sigma and u : u_np1 -> u_n, sigma_np1 -> sigma_n
    // be careful to use this function
    // only if the computation is over
    gmm::copy(sigma_np1, md.set_real_variable(datasigma));
    gmm::copy(u_np1, u_n);

  }



  //=================================================================
  //  Von Mises or Tresca stress computation for elastoplasticity
  //=================================================================

  void compute_elastoplasticity_Von_Mises_or_Tresca
    (model &md,
     const std::string &datasigma,
     const mesh_fem &mf_vm,
     model_real_plain_vector &VM,
     bool tresca) {

    GMM_ASSERT1(gmm::vect_size(VM) == mf_vm.nb_dof(),
                "The vector has not the right size");

    const model_real_plain_vector &sigma_np1 = md.real_variable(datasigma, 0);
    const mesh_fem &mf_sigma = *(md.pmesh_fem_of_variable(datasigma));

    // dimension of the finite element used
    dim_type N = mf_sigma.linked_mesh().dim();

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
                            const pconstraints_projection &ACP,
                            const std::string &datalambda,
                            const std::string &datamu,
                            const std::string &datathreshold,
                            const std::string &datasigma,
                            model_real_plain_vector &plast) {

    const model_real_plain_vector &u_np1 = md.real_variable(varname, 0);
    model_real_plain_vector &u_n = md.set_real_variable(varname, 1);
    const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(varname));

    const model_real_plain_vector &lambda = md.real_variable(datalambda);
    const model_real_plain_vector &mu = md.real_variable(datamu);
    const model_real_plain_vector &threshold = md.real_variable(datathreshold);
    const mesh_fem *pmf_data = md.pmesh_fem_of_variable(datalambda);

    const model_real_plain_vector &sigma_n = md.real_variable(datasigma);
    const mesh_fem &mf_sigma = *(md.pmesh_fem_of_variable(datasigma));

    dim_type N = mf_sigma.linked_mesh().dim();

    mesh_region rg = mim.linked_mesh().get_mpi_region();

    model_real_plain_vector dummyV(mf_u.nb_dof());
    model_real_plain_vector saved_plast(mf_sigma.nb_dof());
    asm_elastoplasticity_rhs(dummyV, &saved_plast,
                             mim, mf_u, mf_sigma, *pmf_data,
                             u_n, u_np1, sigma_n,
                             lambda, mu, threshold, *ACP, PLAST, rg);

    /* Retrieve and save the plastic part */
    GMM_ASSERT1(gmm::vect_size(plast) == mf_pl.nb_dof(),
                "The vector has not the right size");
    GMM_ASSERT1(mf_pl.get_qdim() == 1,
                "Target dimension of mf_pl should be 1");

    base_vector saved_pl(mf_pl.nb_dof()*N*N);
    interpolation(mf_sigma, mf_pl, saved_plast, saved_pl);

    // for each dof we compute the norm of the plastic part
    base_matrix plast_tmp(N, N);
    for (size_type ii = 0; ii < mf_pl.nb_dof(); ++ii) {

      /* we retrieve the matrix sigma_pl on this dof */
      std::copy(saved_pl.begin()+ii*N*N, saved_pl.begin()+(ii+1)*N*N,
                plast_tmp.begin());

      plast[ii] = gmm::mat_euclidean_norm(plast_tmp);

    }

  }





  //=========================================================================
  //
  //  Specific nonlinear operators of the high-level generic assembly
  //  language, useful for plasticity modeling
  //
  //=========================================================================

  static void ga_init_scalar(bgeot::multi_index &mi) { mi.resize(0); }
  static void ga_init_vector(bgeot::multi_index &mi, size_type N)
    { mi.resize(1); mi[0] = N; }
  static void ga_init_matrix(bgeot::multi_index &mi, size_type M, size_type N)
    { mi.resize(2); mi[0] = M; mi[1] = N; }
  static void ga_init_square_matrix(bgeot::multi_index &mi, size_type N)
    { mi.resize(2); mi[0] = mi[1] = N; }


  bool expm(const base_matrix &a_, base_matrix &aexp, scalar_type tol=1e-15) {

    const size_type itmax = 40;
    base_matrix a(a_);
    // scale input matrix a
    int e;
    frexp(gmm::mat_norminf(a), &e);
    e = std::max(0, std::min(1023, e));
    gmm::scale(a, pow(scalar_type(2),-scalar_type(e)));

    base_matrix atmp(a), an(a);
    gmm::copy(a, aexp);
    gmm::add(gmm::identity_matrix(), aexp);
    scalar_type factn(1);
    bool success(false);
    for (size_type n=2; n < itmax; ++n) {
      factn /= scalar_type(n);
      gmm::mult(an, a, atmp);
      gmm::copy(atmp, an);
      gmm::scale(atmp, factn);
      gmm::add(atmp, aexp);
      if (gmm::mat_euclidean_norm(atmp) < tol) {
        success = true;
        break;
      }
    }
    // unscale result
    for (int i=0; i < e; ++i) {
      gmm::mult(aexp, aexp, atmp);
      gmm::copy(atmp, aexp);
    }
    return success;
  }

  bool expm_deriv(const base_matrix &a_, base_tensor &daexp,
                  base_matrix *paexp=NULL, scalar_type tol=1e-15) {

    const size_type itmax = 40;
    size_type N = gmm::mat_nrows(a_);
    size_type N2 = N*N;

    base_matrix a(a_);
    // scale input matrix a
    int e;
    frexp(gmm::mat_norminf(a), &e);
    e = std::max(0, std::min(1023, e));
    scalar_type scale = pow(scalar_type(2),-scalar_type(e));
    gmm::scale(a, scale);

    base_vector factnn(itmax);
    base_matrix atmp(a), an(a), aexp(a);
    base_tensor ann(bgeot::multi_index(N,N,itmax));
    gmm::add(gmm::identity_matrix(), aexp);
    gmm::copy(gmm::identity_matrix(), atmp);
    std::copy(atmp.begin(), atmp.end(), ann.begin());
    factnn[1] = 1;
    std::copy(a.begin(), a.end(), ann.begin()+N2);
    size_type n;
    bool success(false);
    for (n=2; n < itmax; ++n) {
      factnn[n] = factnn[n-1]/scalar_type(n);
      gmm::mult(an, a, atmp);
      gmm::copy(atmp, an);
      std::copy(an.begin(), an.end(), ann.begin()+n*N2);
      gmm::scale(atmp, factnn[n]);
      gmm::add(atmp, aexp);
      if (gmm::mat_euclidean_norm(atmp) < tol) {
        success = true;
        break;
      }
    }

    if (!success)
      return false;

    gmm::clear(daexp.as_vector());
    gmm::scale(factnn, scale);
    for (--n; n >= 1; --n) {
	   scalar_type factn = factnn[n];
      for (size_type m=1; m <= n; ++m)
        for (size_type l=0; l < N; ++l)
          for (size_type k=0; k < N; ++k)
            for (size_type j=0; j < N; ++j)
              for (size_type i=0; i < N; ++i)
                daexp(i,j,k,l) += factn*ann(i,k,m-1)*ann(l,j,n-m);
	 }

    // unscale result
    base_matrix atmp1(a), atmp2(a);
    for (int i=0; i < e; ++i) {
      for (size_type l=0; l < N; ++l)
        for (size_type k=0; k < N; ++k) {
          std::copy(&daexp(0,0,k,l), &daexp(0,0,k,l)+N*N, atmp.begin());
          gmm::mult(atmp, aexp, atmp1);
          gmm::mult(aexp, atmp, atmp2);
          gmm::add(atmp1, atmp2, atmp);
          std::copy(atmp.begin(), atmp.end(), &daexp(0,0,k,l));
        }
      gmm::mult(aexp, aexp, atmp);
      gmm::copy(atmp, aexp);
    }

    if (paexp) gmm::copy(aexp, *paexp);
    return true;
  }

  // numerical differantiation of logm
  // not used becaused it caused some issues and was slower than
  // simply inverting the derivative of expm
  bool logm_deriv(const base_matrix &a, base_tensor &dalog,
                  base_matrix *palog=NULL) {

    size_type N = gmm::mat_nrows(a);
    base_matrix a1(a), alog1(a), alog(a);
    logm(a, alog);
    scalar_type eps(1e-8);
    for (size_type k=0; k < N; ++k)
      for (size_type l=0; l < N; ++l) {
        gmm::copy(a, a1);
        a1(k,l) += eps;
        gmm::logm(a1, alog1);
        for (size_type i=0; i < N; ++i)
          for (size_type j=0; j < N; ++j)
            dalog(i,j,k,l) = (alog1(i,j) - alog(i,j))/eps;
      }
    if (palog) gmm::copy(alog, *palog);
    return true;
  }


  // Matrix exponential
  struct matrix_exponential_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() != 2
          || args[0]->sizes()[0] != args[0]->sizes()[1]) return false;
      ga_init_matrix(sizes, args[0]->sizes()[0], args[0]->sizes()[1]);
      return true;
    }

    // Value:
    void value(const arg_list &args, base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix inpmat(N,N), outmat(N,N);
      gmm::copy(args[0]->as_vector(), inpmat.as_vector());
      bool info = expm(inpmat, outmat);
      GMM_ASSERT1(info, "Matrix exponential calculation "
                        "failed to converge");
      gmm::copy(outmat.as_vector(), result.as_vector());
    }

    // Derivative:
    void derivative(const arg_list &args, size_type /*nder*/,
                    base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix inpmat(N,N);
      gmm::copy(args[0]->as_vector(), inpmat.as_vector());
      bool info = expm_deriv(inpmat, result);
      GMM_ASSERT1(info, "Matrix exponential derivative calculation "
                        "failed to converge");
    }

    // Second derivative : not implemented
    void second_derivative(const arg_list &, size_type, size_type,
                           base_tensor &) const {
      GMM_ASSERT1(false, "Sorry, second derivative not implemented");
    }
  };


  // Matrix logarithm
  struct matrix_logarithm_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() != 2
          || args[0]->sizes()[0] != args[0]->sizes()[1]) return false;
      ga_init_matrix(sizes, args[0]->sizes()[0], args[0]->sizes()[1]);
      return true;
    }

    // Value:
    void value(const arg_list &args, base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix inpmat(N,N), outmat(N,N);
      gmm::copy(args[0]->as_vector(), inpmat.as_vector());
      gmm::logm(inpmat, outmat);
      gmm::copy(outmat.as_vector(), result.as_vector());
    }

    // Derivative:
    void derivative(const arg_list &args, size_type /*nder*/,
                    base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix inpmat(N,N), outmat(N,N), tmpmat(N*N,N*N);
      gmm::copy(args[0]->as_vector(), inpmat.as_vector());
      gmm::logm(inpmat, outmat);
      bool info = expm_deriv(outmat, result);
      if (info) {
        gmm::copy(result.as_vector(), tmpmat.as_vector());
        scalar_type det = gmm::lu_inverse(tmpmat);
        if (det <= 0) gmm::copy(gmm::identity_matrix(), tmpmat);
        gmm::copy(tmpmat.as_vector(), result.as_vector());
      }
      GMM_ASSERT1(info, "Matrix logarithm derivative calculation "
                        "failed to converge");
    }

    // Second derivative : not implemented
    void second_derivative(const arg_list &, size_type, size_type,
                           base_tensor &) const {
      GMM_ASSERT1(false, "Sorry, second derivative not implemented");
    }
  };


  // Normalized vector/matrix operator : Vector/matrix divided by its Frobenius norm
  struct normalized_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() > 2
          || args[0]->sizes().size() < 1) return false;
      if (args[0]->sizes().size() == 1)
        ga_init_vector(sizes, args[0]->sizes()[0]);
      else
        ga_init_matrix(sizes, args[0]->sizes()[0], args[0]->sizes()[1]);
      return true;
    }
    
    // Value : u/|u|
    void value(const arg_list &args, base_tensor &result) const
    {
      scalar_type no = gmm::vect_norm2(args[0]->as_vector());
      if (no < 1E-15)
        gmm::clear(result.as_vector());
      else
        gmm::copy(gmm::scaled(args[0]->as_vector(), scalar_type(1)/no),
                  result.as_vector());
    }

    // Derivative : (|u|^2 Id - u x u)/|u|^3
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const {
      const base_tensor &t = *args[0];
      size_type N = t.size();
      scalar_type no = gmm::vect_norm2(t.as_vector());

      gmm::clear(result.as_vector());
      if (no >= 1E-15) {
        scalar_type no3 = no*no*no;
        for (size_type i = 0; i < N; ++i) {
          result[i*N+i] += scalar_type(1)/no;
          for (size_type j = 0; j < N; ++j)
            result[j*N+i] -= t[i]*t[j] / no3;
        }
      }
    }

    // Second derivative : not implemented
    void second_derivative(const arg_list &/*args*/, size_type, size_type,
                           base_tensor &/*result*/) const {
      GMM_ASSERT1(false, "Sorry, second derivative not implemented");
    }
  };


  //=================================================================
  // Von Mises projection
  //=================================================================


  struct Von_Mises_projection_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 2 || args[0]->sizes().size() > 2
          || args[1]->size() != 1) return false;
      size_type N = (args[0]->sizes().size() == 2) ?  args[0]->sizes()[0] : 1;
      if (args[0]->sizes().size() == 2 && args[0]->sizes()[1] != N) return false;
      if (args[0]->sizes().size() != 2 && args[0]->size() != 1)  return false;
      if (N > 1) ga_init_square_matrix(sizes, N); else ga_init_scalar(sizes);
      return true;
    }

    // Value:
    void value(const arg_list &args, base_tensor &result) const {
      size_type N = (args[0]->sizes().size() == 2) ? args[0]->sizes()[0] : 1;
      base_matrix tau(N, N), tau_D(N, N);
      gmm::copy(args[0]->as_vector(), tau.as_vector());
      scalar_type s = (*(args[1]))[0];

      scalar_type tau_m = gmm::mat_trace(tau) / scalar_type(N);
      gmm::copy(tau, tau_D);
      for (size_type i = 0; i < N; ++i) tau_D(i,i) -= tau_m;

      scalar_type norm_tau_D = gmm::mat_euclidean_norm(tau_D);
      
      if (norm_tau_D > s) gmm::scale(tau_D, s / norm_tau_D);

      for (size_type i = 0; i < N; ++i) tau_D(i,i) += tau_m;

      gmm::copy(tau_D.as_vector(), result.as_vector());
    }

    // Derivative:
    void derivative(const arg_list &args, size_type nder,
                    base_tensor &result) const {
      size_type N = (args[0]->sizes().size() == 2) ? args[0]->sizes()[0] : 1;
      base_matrix tau(N, N), tau_D(N, N);
      gmm::copy(args[0]->as_vector(), tau.as_vector());
      scalar_type s = (*(args[1]))[0];
      scalar_type tau_m = gmm::mat_trace(tau) / scalar_type(N);
      gmm::copy(tau, tau_D);
      for (size_type i = 0; i < N; ++i) tau_D(i,i) -= tau_m;
      scalar_type norm_tau_D = gmm::mat_euclidean_norm(tau_D);

      if (norm_tau_D != scalar_type(0))
        gmm::scale(tau_D, scalar_type(1)/norm_tau_D);

      switch(nder) {
      case 1: 
        if (norm_tau_D <= s) {
          gmm::clear(result.as_vector());
          for (size_type i = 0; i < N; ++i)
            for (size_type j = 0; j < N; ++j)
              result(i,j,i,j) = scalar_type(1);
        } else {
          for (size_type i = 0; i < N; ++i)
            for (size_type j = 0; j < N; ++j)
              for (size_type m = 0; m < N; ++m)
                for (size_type n = 0; n < N; ++n)
                  result(i,j,m,n)
                    = s * (-tau_D(i,j) * tau_D(m,n)
                           + ((i == m && j == n) ? scalar_type(1) : scalar_type(0))
                           - ((i == j && m == n) ? scalar_type(1)/scalar_type(N)
                              : scalar_type(0))) / norm_tau_D;
          for (size_type i = 0; i < N; ++i)
            for (size_type j = 0; j < N; ++j)
              result(i,i,j,j) += scalar_type(1)/scalar_type(N);
        }
        break;
      case 2:
        if (norm_tau_D < s)
          gmm::clear(result.as_vector());
        else
          gmm::copy(tau_D.as_vector(), result.as_vector());
        break;
      }
    }

    // Second derivative : not implemented
    void second_derivative(const arg_list &, size_type, size_type,
                           base_tensor &) const {
      GMM_ASSERT1(false, "Sorry, second derivative not implemented");
    }
  };

  static bool init_predef_operators(void) {

    ga_predef_operator_tab &PREDEF_OPERATORS
      = dal::singleton<ga_predef_operator_tab>::instance();

    PREDEF_OPERATORS.add_method("Expm",
                                std::make_shared<matrix_exponential_operator>());
    PREDEF_OPERATORS.add_method("Logm",
                                std::make_shared<matrix_logarithm_operator>());
    PREDEF_OPERATORS.add_method("Normalized",
                                std::make_shared<normalized_operator>());
    PREDEF_OPERATORS.add_method("Von_Mises_projection",
                                std::make_shared<Von_Mises_projection_operator>());
    return true;
   }

  // declared in getfem_generic_assembly.cc
  extern bool predef_operators_plasticity_initialized
    = init_predef_operators();

}  /* end of namespace getfem.  */

