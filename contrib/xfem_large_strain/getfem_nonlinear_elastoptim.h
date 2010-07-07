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


#include "getfem/getfem_modeling.h"
#include "getfem/getfem_models.h"
#include "getfem/getfem_assembling_tensors.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_assembling.h" 
#include "getfem/getfem_export.h"   
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_nonlinear_elasticity.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_spider_fem.h"
#include "getfem/getfem_mesh_fem_sum.h"
#include "getfem/getfem_superlu.h"
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"



namespace getfem {



  inline scalar_type mat_euclidean_sp(base_matrix &Sigma, base_matrix &Ealpha) {
    return gmm::vect_sp((std::vector<scalar_type> &)Sigma, (std::vector<scalar_type> &)Ealpha);
  }


  //==========================================================================
  //=               /*************************/                              =
  //=               /* Term 1 dL/dbeta       */                              =
  //=               /*************************/                              =
  //==========================================================================



  template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term1
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem &mf_p;
    std::vector<scalar_type> P, P_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old,  NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, valP;
    base_matrix E, Sigma, gradU;
    base_tensor tt;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term1
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem &mf_p_, const VECT1 &P_, const VECT1 &P_ls_,
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_p(mf_p_),
	P(mf_p_.nb_basic_dof()), P_ls(mf_p_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), E(N, N), Sigma(N, N), gradU(NFem, N),
	tt(N, N, N, N), sizes_(NFem, N, NFem, N) {
      sizes_.resize(1); sizes_[0] = 1;
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_p.extend_vector(P_, P);
      mf_u.extend_vector(U_ls_, U_ls);
      mf_p.extend_vector(P_ls_, P_ls);
      if (gmm::vect_size(PARAMS) == AHL_.nb_params())
	gmm::copy(PARAMS, params);
    }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf_u.get_qdim());

      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)+= scalar_type(1);
      
      scalar_type J = gmm::lu_det(gradU);

      t[0] = (J - scalar_type(1)) * valP[0];
    }
    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
      if (nb == 1) {

	valP.resize(1);
	size_type cv = ctx.convex_num();
	coeff.resize(mf_p.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector
	     (P_ls, gmm::sub_index(mf_p.ind_basic_dof_of_element(cv))), coeff);
	ctx.pf()->interpolation(ctx, coeff, valP, 1);
	

	if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
	}
	scalar_type x = mls_x(ctx.xref());
	scalar_type y = mls_y(ctx.xref());

	valP[0] *= log(x*x+y*y) / scalar_type(2);
      } else if (mf_data) {
	size_type cv = ctx.convex_num();
	size_type nbp = AHL.nb_params();
	coeff.resize(mf_data->nb_basic_dof_of_element(cv)*nbp);
	for (size_type i = 0; i < mf_data->nb_basic_dof_of_element(cv); ++i)
	  for (size_type k = 0; k < nbp; ++k)
	    coeff[i * nbp + k]
	      = PARAMS[mf_data->ind_basic_dof_of_element(cv)[i]*nbp+k];
	ctx.pf()->interpolation(ctx, coeff, params, dim_type(nbp));
      }
    } 
    
  };


//============================================================================
//=                /******************************/                          =
//=                /*      Term 2 dL/dalpha      */                          =
//=                /******************************/                          =
//============================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term2
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem &mf_p;
    std::vector<scalar_type> P, P_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, valP, u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, U_ls_theta;// U_ls_theta le terme U_ls multiplier par le vecteur 
    base_tensor tt;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term2
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem &mf_p_, const VECT1 &P_, const VECT1 &P_ls_,
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_p(mf_p_),
	P(mf_p_.nb_basic_dof()), P_ls(mf_p_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N), gradU(NFem, N),
	tt(N, N, N, N), sizes_(NFem, N, NFem, N) {
      sizes_.resize(1); sizes_[0] = 1;
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_p.extend_vector(P_, P);
      mf_u.extend_vector(U_ls_, U_ls);
      mf_p.extend_vector(P_ls_, P_ls);
      if (gmm::vect_size(PARAMS) == AHL_.nb_params())
	gmm::copy(PARAMS, params);
    }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf_u.get_qdim());


      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));


      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)+= scalar_type(1);
      
      scalar_type J = gmm::lu_det(gradU);


      // Computation of d/dalpha (grad U)

      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U_ls, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU_ls, mf_u.get_qdim());
      scalar_type x = mls_x(ctx.xref());
      scalar_type y = mls_y(ctx.xref());
      scalar_type r2 = x*x+y*y;
      gmm::scale(gradU_ls, log(r2) / scalar_type(2));

      ctx.pf()->interpolation(ctx, coeff, u_ls, mf_u.get_qdim());
      base_vector Runit(2);
      Runit[0] = x / r2; Runit[1] = y / r2; 
      gmm::rank_one_update(gradU_ls, Runit, u_ls);

      // Computation of d(E)/d(alpha)
      
      base_matrix Ealpha(N, N);
      gmm::mult(gmm::transposed(gradU_ls), gradU, Sigma);
      gmm::add(Sigma, gmm::transposed(Sigma), Ealpha);
      gmm::scale(Ealpha, 0.5);

      // The term sigma:dE/dalpha
      AHL.sigma(E, Sigma, params);
      t[0] = mat_euclidean_sp(Sigma, Ealpha);

      // The term pJF^{-T}:dgradu/dalpha
      gmm::lu_inverse(gradU);
      base_matrix FmT(N,N);
      gmm::copy(gmm::transposed(gradU), FmT);
      t[0] += valP[0] * J * mat_euclidean_sp(gradU_ls, FmT);
    }
    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
      if (nb == 1) {
	size_type cv = ctx.convex_num();

	valP.resize(1);
	coeff.resize(mf_p.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector
	     (P_ls, gmm::sub_index(mf_p.ind_basic_dof_of_element(cv))), coeff);
	ctx.pf()->interpolation(ctx, coeff, valP, 1);

	if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
	}

      } else if (mf_data) {
	size_type cv = ctx.convex_num();
	size_type nbp = AHL.nb_params();
	coeff.resize(mf_data->nb_basic_dof_of_element(cv)*nbp);
	for (size_type i = 0; i < mf_data->nb_basic_dof_of_element(cv); ++i)
	  for (size_type k = 0; k < nbp; ++k)
	    coeff[i * nbp + k]
	      = PARAMS[mf_data->ind_basic_dof_of_element(cv)[i]*nbp+k];
	ctx.pf()->interpolation(ctx, coeff, params, dim_type(nbp));
      }
    } 
    
  };

/******************************************************************************************************************************/

  //===================================================================
  //=                                                                 =
  //=             Terms of the tangent matrix for optimization        =
  //=                                                                 =
  //===================================================================

//==================================================================================
//=     /***********************************************************/              =
//=     /* Term 3 d^2(L)/d(alpha)d(u^h)    The part [*]:Grad(v^h)  */              =
//=     /***********************************************************/              =
//==================================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term3
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem &mf_p;
    std::vector<scalar_type> P, P_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, valP, u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, U_ls_theta, GsEalpha;// U_ls_theta le terme U_ls multiplier par le vecteur 
    base_tensor GSigma;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term3
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem &mf_p_, const VECT1 &P_, const VECT1 &P_ls_,
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_p(mf_p_),
	P(mf_p_.nb_basic_dof()), P_ls(mf_p_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N), gradU(NFem, N),
	GSigma(N, N, N, N) ,sizes_(NFem, N, NFem, N) {
      sizes_.resize(2); sizes_[0] = N; sizes_[1] = N;
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_p.extend_vector(P_, P);
      mf_u.extend_vector(U_ls_, U_ls);
      mf_p.extend_vector(P_ls_, P_ls);
      if (gmm::vect_size(PARAMS) == AHL_.nb_params())
	gmm::copy(PARAMS, params);
    }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf_u.get_qdim());


      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));


      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)=gradU(alpha, alpha)+ scalar_type(1);
      scalar_type J = gmm::lu_det(gradU);
      

      // Computation of d(grad U)/d(alpha) 

      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U_ls, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU_ls, mf_u.get_qdim());
      scalar_type x = mls_x(ctx.xref());
      scalar_type y = mls_y(ctx.xref());
      scalar_type r2 = x*x+y*y;
      gmm::scale(gradU_ls, log(r2) / scalar_type(2));

      ctx.pf()->interpolation(ctx, coeff, u_ls, mf_u.get_qdim());
      base_vector Runit(2);
      Runit[0] = x / r2; Runit[1] = y / r2; 
      gmm::rank_one_update(gradU_ls, Runit, u_ls);

      // Computation of d(E)/d(alpha) 
      
      base_matrix Ealpha(N, N);
      gmm::mult(gmm::transposed(gradU_ls), gradU, Sigma);
      gmm::add(Sigma, gmm::transposed(Sigma), Ealpha);
      gmm::scale(Ealpha, 0.5);

      // The term sigma.d(E)/d(alpha)
      base_matrix eas(N,N);
      AHL.sigma(E, Sigma, params);
      gmm::mult(Ealpha, Sigma, eas);

      // The term grad_sigma:d(E)/d(alpha)
      scalar_type temp ;
      AHL.grad_sigma(E, GSigma, params);
      for (size_type i = 0; i < N; ++i)
	for (size_type j = 0; j < N; ++j){
	  temp = 0.0;
	  for (size_type k = 0; k < N; ++k)
	    for (size_type l = 0; l < N; ++l) {
	      temp += GSigma(i,j,k,l)*Ealpha(l,k);	     	
	    }
	  GsEalpha(i,j)=temp;
	}

      // The term F.grad_sigma:d(E)/d(alpha)
      gmm::mult(gradU, GsEalpha, Sigma);
      gmm::add(Sigma, eas);// on a plus besoin du contenue de GsEalpha
      
      // The term pJ[F^{-T}:d(grad U)/d(alpha)].F^{-T} - JF^{-T}.d(grad U)/d(alpha).F^{-T}
      gmm::lu_inverse(gradU);
      base_matrix FmT(N,N);
      gmm::copy(gmm::transposed(gradU), FmT);
      temp = valP[0] * J * mat_euclidean_sp(gradU_ls, gradU);
      gmm::add(gmm::scaled(FmT,temp), eas);
      gmm::mult(FmT, gmm::transposed(gradU_ls), GsEalpha);
      gmm::mult(GsEalpha, FmT, Sigma);
      gmm::add(gmm::scaled(Sigma,-1.0 * J * valP[0]), eas);
      
      for (size_type i = 0; i < N; ++i)
	for (size_type j = 0; j < N; ++j)
	  t(i,j) = eas(i, j);
    }

    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
      if (nb == 1) {
	size_type cv = ctx.convex_num();

	valP.resize(1);
	coeff.resize(mf_p.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector
	     (P_ls, gmm::sub_index(mf_p.ind_basic_dof_of_element(cv))), coeff);
	ctx.pf()->interpolation(ctx, coeff, valP, 1);

	if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
	}

      } else if (mf_data) {
	size_type cv = ctx.convex_num();
	size_type nbp = AHL.nb_params();
	coeff.resize(mf_data->nb_basic_dof_of_element(cv)*nbp);
	for (size_type i = 0; i < mf_data->nb_basic_dof_of_element(cv); ++i)
	  for (size_type k = 0; k < nbp; ++k)
	    coeff[i * nbp + k]
	      = PARAMS[mf_data->ind_basic_dof_of_element(cv)[i]*nbp+k];
	ctx.pf()->interpolation(ctx, coeff, params, dim_type(nbp));
      }
    } 
    
  };


//==============================================================================================
//=     /***********************************************************************/              =
//=     /* Term 4 d^2(L)/d(alpha)d(u^h)    The part [*]:d(Grad(v^h))/d(alpha)  */              =
//=     /***********************************************************************/              =
//==============================================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term4
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem &mf_p;
    std::vector<scalar_type> P, P_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, valP, u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, U_ls_theta, GsEalpha;// U_ls_theta le terme U_ls multiplier par le vecteur 
    base_tensor GSigma;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term4
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem &mf_p_, const VECT1 &P_, const VECT1 &P_ls_,
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_p(mf_p_),
	P(mf_p_.nb_basic_dof()), P_ls(mf_p_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N), gradU(NFem, N),
	GSigma(N, N, N, N) ,sizes_(NFem, N, NFem, N) {
      sizes_.resize(2); sizes_[0] = N; sizes_[1] = N;
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_p.extend_vector(P_, P);
      mf_u.extend_vector(U_ls_, U_ls);
      mf_p.extend_vector(P_ls_, P_ls);
      if (gmm::vect_size(PARAMS) == AHL_.nb_params())
	gmm::copy(PARAMS, params);
    }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf_u.get_qdim());

      scalar_type x = mls_x(ctx.xref());
      scalar_type y = mls_y(ctx.xref());
      scalar_type r2 = x*x+y*y;

      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));


      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)=gradU(alpha, alpha)+ scalar_type(1);
      scalar_type J = gmm::lu_det(gradU);
      

      // Computation of F.Sigma 
      AHL.sigma(E, Sigma, params);
      base_matrix eas(N,N);
      gmm::mult(gradU, Sigma, eas);
      // Computation of p.J.F^{-T}
      gmm::lu_inverse(gradU);
      gmm::add(gmm::scaled(gmm::transposed(gradU), valP[0] * J), eas);     
      gmm::scale(eas, log(r2)*0.5);      
    }

    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
      if (nb == 1) {
	size_type cv = ctx.convex_num();

	valP.resize(1);
	coeff.resize(mf_p.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector
	     (P_ls, gmm::sub_index(mf_p.ind_basic_dof_of_element(cv))), coeff);
	ctx.pf()->interpolation(ctx, coeff, valP, 1);

	if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
	}

      } else if (mf_data) {
	size_type cv = ctx.convex_num();
	size_type nbp = AHL.nb_params();
	coeff.resize(mf_data->nb_basic_dof_of_element(cv)*nbp);
	for (size_type i = 0; i < mf_data->nb_basic_dof_of_element(cv); ++i)
	  for (size_type k = 0; k < nbp; ++k)
	    coeff[i * nbp + k]
	      = PARAMS[mf_data->ind_basic_dof_of_element(cv)[i]*nbp+k];
	ctx.pf()->interpolation(ctx, coeff, params, dim_type(nbp));
      }
    } 
    
  };


//==============================================================================================
//=     /***********************************************************************/              =
//=     /* Term 5 d^2(L)/d(alpha)d(u^h)    The part [*]:d(Grad(v^h))/d(alpha)  */              =
//=     /***********************************************************************/              =
//==============================================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term5
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem &mf_p;
    std::vector<scalar_type> P, P_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, valP, u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, U_ls_theta, GsEalpha;// U_ls_theta le terme U_ls multiplier par le vecteur 
    base_tensor GSigma;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term5
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem &mf_p_, const VECT1 &P_, const VECT1 &P_ls_,
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_p(mf_p_),
	P(mf_p_.nb_basic_dof()), P_ls(mf_p_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N), gradU(NFem, N),
	GSigma(N, N, N, N) ,sizes_(NFem, N, NFem, N) {
      sizes_.resize(1); sizes_[0] = N;
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_p.extend_vector(P_, P);
      mf_u.extend_vector(U_ls_, U_ls);
      mf_p.extend_vector(P_ls_, P_ls);
      if (gmm::vect_size(PARAMS) == AHL_.nb_params())
	gmm::copy(PARAMS, params);
    }
    const bgeot::multi_index &sizes() const {  return sizes_; }
    virtual void compute(getfem::fem_interpolation_context& ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU, mf_u.get_qdim());

      scalar_type x = mls_x(ctx.xref());
      scalar_type y = mls_y(ctx.xref());
      scalar_type r2 = x*x+y*y;

      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));


      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)=gradU(alpha, alpha)+ scalar_type(1);
      scalar_type J = gmm::lu_det(gradU);
      

      // Computation of F.Sigma 
      AHL.sigma(E, Sigma, params);
      base_matrix eas(N,N);
      gmm::mult(gradU, Sigma, eas);
      // Computation of p.J.F^{-T}
      gmm::lu_inverse(gradU);
      gmm::add(gmm::scaled(gmm::transposed(gradU), valP[0] * J), eas); 
      base_small_vector V(2), W(2); V[0] = x; V[1] = y;
      gmm::mult(gmm::transposed(eas), gmm::scaled(V, 1./r2), W);
    }

    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
      if (nb == 1) {
	size_type cv = ctx.convex_num();

	valP.resize(1);
	coeff.resize(mf_p.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector
	     (P_ls, gmm::sub_index(mf_p.ind_basic_dof_of_element(cv))), coeff);
	ctx.pf()->interpolation(ctx, coeff, valP, 1);

	if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
	}

      } else if (mf_data) {
	size_type cv = ctx.convex_num();
	size_type nbp = AHL.nb_params();
	coeff.resize(mf_data->nb_basic_dof_of_element(cv)*nbp);
	for (size_type i = 0; i < mf_data->nb_basic_dof_of_element(cv); ++i)
	  for (size_type k = 0; k < nbp; ++k)
	    coeff[i * nbp + k]
	      = PARAMS[mf_data->ind_basic_dof_of_element(cv)[i]*nbp+k];
	ctx.pf()->interpolation(ctx, coeff, params, dim_type(nbp));
      }
    } 
    
  };



  /***********************************************************/


  //=========================================================================
  //
  //  Nonlinear elasticity Brick
  //
  //=========================================================================


  //=========================================================================
  //  Assembling the tangent of the additional term
  //=========================================================================


template<typename VECT0, typename VECT1, typename VECT2> 
  
  void asm_nonlinear_elasticity_optim_tangent_matrix
  (const VECT0 &V_, const mesh_im &mim, 
   const getfem::mesh_fem &mf_u, const VECT1 &U, const VECT2 &U_ls,
   dal::bit_vector &u_enriched_dof,
   const getfem::mesh_fem &mf_p, const VECT1 &P, const VECT2 &P_ls,
   dal::bit_vector &p_enriched_dof,
   const VECT2 &alpha, const VECT2 &beta,
   const getfem::mesh_fem *mf_data, const VECT2 &PARAMS,
   const abstract_hyperelastic_law &AHL, const getfem::level_set &ls,
   const mesh_region &rg = mesh_region::all_convexes()) {
    VECT0 &V = const_cast<VECT0 &>(V_);
    GMM_ASSERT1(mf_u.get_qdim() >= mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    elasticity_nonlinear_optim_term3<VECT1, VECT2>
      nterm1(mf_u, U, U_ls, mf_p, P, P_ls, mf_data, PARAMS, AHL, ls);
    elasticity_nonlinear_optim_term4<VECT1, VECT2>
      nterm2(mf_u, U, U_ls, mf_p, P, P_ls, mf_data, PARAMS, AHL, ls);
    elasticity_nonlinear_optim_term5<VECT1, VECT2>
      nterm3(mf_u, U, U_ls, mf_p, P, P_ls, mf_data, PARAMS, AHL, ls);



    getfem::generic_assembly assem2;
    if (mf_data)
      assem2.set("V(#1)+=comp(NonLin$1(#1,#2)(i,j).vGrad(#1)(:,i,j));"
		 "V(#1)+=comp(NonLin$2(#1,#2)(i).vBase(#1)(:,i))");
    else
      assem2.set("V(#1)+=comp(NonLin$1(#1)(i,j).vGrad(#1)(:,i,j));"
		 "V(#1)+=comp(NonLin$2(#1)(i).vBase(#1)(:,i))");
    assem2.push_mi(mim);
    assem2.push_mf(mf_u);
    if (mf_data) assem2.push_mf(*mf_data);
    assem2.push_nonlinear_term(&nterm2);
    assem2.push_nonlinear_term(&nterm3);
    assem2.push_vec(V);
    assem2.assembly(rg);

    for (size_type i = 0; i < mf_u.nb_dof(); ++i)
      if (!(u_enriched_dof.is_in(i))) V[i] = 0.0;


    getfem::generic_assembly assem1;
    
    if (mf_data)
      assem1.set("V(#1)+=comp(NonLin(#1,#2)(i,j).vGrad(#1)(:,i,j))");
    else
      assem1.set("V(#1)+=comp(NonLin(#1)(i,j).vGrad(#1)(:,i,j))");
    assem1.push_mi(mim);
    assem1.push_mf(mf_u);
    if (mf_data) assem1.push_mf(*mf_data);
    assem1.push_nonlinear_term(&nterm1);
    assem1.push_vec(V);
    assem1.assembly(rg);

  }



  //========================================================================
  //=       Assembling the right hand side of the additional term          =
  //========================================================================

  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_elasticity_optim_rhs
  (const VECT1 &R1_, const VECT1 &R2_, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U, const VECT3 &U_ls,
   dal::bit_vector &u_enriched_dof,
   const getfem::mesh_fem &mf_p, const VECT2 &P, const VECT3 &P_ls,
   dal::bit_vector &p_enriched_dof,
   const VECT2 &alpha, const VECT2 &beta,
   const getfem::mesh_fem *mf_data, const VECT3 &PARAMS,
   const abstract_hyperelastic_law &AHL, const getfem::level_set &ls,
   const mesh_region &rg = mesh_region::all_convexes()) {
    VECT1 &R1 = const_cast<VECT1 &>(R1_);
    VECT1 &R2 = const_cast<VECT1 &>(R2_);
    GMM_ASSERT1(mf_u.get_qdim() >= mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");


    GMM_ASSERT1(!mf_u.is_reduced(), "Do not work for reduced fems");
    GMM_ASSERT1(!mf_p.is_reduced(), "Do not work for reduced fems");

    {
      
      elasticity_nonlinear_optim_term1<VECT2, VECT3>
	nterm(mf_u, U, U_ls, mf_p, P, P_ls, mf_data, PARAMS, AHL, ls);
      
      getfem::generic_assembly assem;
      if (mf_data)
	assem.set("t=comp(NonLin(#1,#2,#3))); V() += t(i)");
      else
	assem.set("t=comp(NonLin(#1,#2)); V() += t(i)");
      // comp() to be optimized ?
      assem.push_mi(mim);
      assem.push_mf(mf_u);
      assem.push_mf(mf_p);
      if (mf_data) assem.push_mf(*mf_data);
      assem.push_nonlinear_term(&nterm);
      assem.push_vec(R2);
      assem.assembly(rg);
      
    }
    
    {
      elasticity_nonlinear_optim_term2<VECT2, VECT3>
	nterm(mf_u, U, U_ls, mf_p, P, P_ls, mf_data, PARAMS, AHL, ls);
      
      getfem::generic_assembly assem;
      if (mf_data)
	assem.set("t=comp(NonLin(#1,#2,#3))); V() += t(i)");
      else
	assem.set("t=comp(NonLin(#1,#2)); V() += t(i)");
      // comp() to be optimized ?
      assem.push_mi(mim);
      assem.push_mf(mf_u);
      assem.push_mf(mf_p);
      if (mf_data) assem.push_mf(*mf_data);
      assem.push_nonlinear_term(&nterm);
      assem.push_vec(R1);
      assem.assembly(rg);

    }
    

  }









}  /* end of namespace getfem.                                             */

