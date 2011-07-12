//===========================================================================
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



  //  inline scalar_type mat_euclidean_sp(base_matrix &Sigma, base_matrix &Ealpha) {
  //  return gmm::vect_sp((std::vector<scalar_type> &)Sigma, (std::vector<scalar_type> &)Ealpha);
  //
  //  }


  

//=======================================================================
//=                /******************************/                     =
//=  compressible  /*   Term 1 d(L)/d(alpha)     */                     =
//=                /******************************/                     =
//=======================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term_compressible1
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
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
   elasticity_nonlinear_optim_term_compressible1 
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_,
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N),
	gradU(NFem, N), gradU_ls(NFem, N),
	tt(N, N, N, N), sizes_(NFem, N, NFem, N) {
      sizes_.resize(1); sizes_[0] = 1;
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_u.extend_vector(U_ls_, U_ls);
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
      
      //scalar_type J = gmm::lu_det(gradU);


      // Computation of d(E)/d(alpha)
      
      base_matrix Ealpha(N, N);
      gmm::mult(gmm::transposed(gradU_ls), gradU, Sigma);
      gmm::add(Sigma, gmm::transposed(Sigma), Ealpha);
      gmm::scale(Ealpha, 0.5);

      // The term sigma:dE/dalpha
      AHL.sigma(E, Sigma, params);
      t[0] = mat_euclidean_sp(Sigma, Ealpha);

     
    }
    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");

      if (mf_data) {
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

/***********************************************************************************************************************/
/***********************************************************************************************************************/
/***********************************************************************************************************************/
/***********************************************************************************************************************/
/***********************************************************************************************************************/

  //===================================================================
  //=                                                                 =
  //= Compressible Terms of the tangent matrix for optimization       =
  //=                                                                 =
  //===================================================================

//==================================================================================
//=     /***********************************************************/              =
//=     /* Term 2 d^2(L)/d(alpha)d(u^h)   The part [*]:Grad(v^h)   */              = Compressible
//=     /***********************************************************/              =
//==================================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term_compressible2
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff,  u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, U_ls_theta;// U_ls_theta le terme U_ls multiplier par le vecteur 
    base_tensor GSigma;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term_compressible2
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N),
	gradU(NFem, N), gradU_ls(NFem, N),
	GSigma(N, N, N, N) ,sizes_(NFem, N, NFem, N) {
      sizes_.resize(2); sizes_[0] = short_type(N); sizes_[1] = short_type(N);
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_u.extend_vector(U_ls_, U_ls);
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

      // Computation of tensor E 

      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));

      // Computation of J determinant of tensor F 

      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)=gradU(alpha, alpha)+ scalar_type(1); // gradU contient le tenseur F
      // scalar_type J = gmm::lu_det(gradU);
      

      // Computation of d(gradU)/d(alpha) 

      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U_ls, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU_ls, mf_u.get_qdim());

      if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
      }

      scalar_type x = mls_x(ctx.xref());
      
      scalar_type y = mls_y(ctx.xref());
      

      scalar_type r2 = x*x+y*y;
      gmm::scale(gradU_ls, 0.5*log(r2));

      ctx.pf()->interpolation(ctx, coeff, u_ls, mf_u.get_qdim());
      base_vector Runit(2);
      Runit[0] = x / r2; Runit[1] = y / r2; 
      gmm::rank_one_update(gradU_ls, Runit, u_ls); // gradU_ls contain d(gradU)/d(alpha) 

      // Computation of d(E)/d(alpha) 
      
      base_matrix Ealpha(N, N), GsEalpha(N, N);
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
	      temp += GSigma(i,j,k,l) * Ealpha(l,k);	     	
	    }
	  GsEalpha(i,j)=temp;
	}

      // The term F.grad_sigma:d(E)/d(alpha)
      gmm::mult(gradU, GsEalpha, Sigma);
      gmm::add(Sigma, eas);// on a plus besoin du contenue de GsEalpha
            
      
      for (size_type i = 0; i < N; ++i)
	for (size_type j = 0; j < N; ++j)
	  t(i,j) = eas(i, j);
    }

    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
     if (mf_data) {
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
//=     /* Term 3 d^2(L)/d(alpha)d(u^h)    The part [*]:Ln(r).Grad(v^h)        */              = Compressible
//=     /***********************************************************************/              =
//==============================================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term_compressible3
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, U_ls_theta;// U_ls_theta le terme U_ls multiplier par le vecteur 
    base_tensor GSigma;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term_compressible3
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N),
	gradU(NFem, N), gradU_ls(NFem, N),
	GSigma(N, N, N, N) ,sizes_(NFem, N, NFem, N) {
      sizes_.resize(2); sizes_[0] = short_type(N); sizes_[1] = short_type(N);
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      
      mf_u.extend_vector(U_ls_, U_ls);
      
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
      
      if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
      }

      scalar_type x = mls_x(ctx.xref());
      scalar_type y = mls_y(ctx.xref());
      scalar_type r2 = x*x+y*y;
      // Computation of tensor E
      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));

      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)=gradU(alpha, alpha)+ scalar_type(1); // GradU contient le tenseur F maintenant
      // scalar_type J = gmm::lu_det(gradU);

      /*****************************/      
      /*****************************/

      // Computation of F.Sigma 
      AHL.sigma(E, Sigma, params);
      base_matrix eas(N,N);
      gmm::mult(gradU, Sigma, eas);
      // cout << "Avant multiplication log :" << eas << endl ;
      // cout << "Valeur du rayon :" << r2 << endl ;
      // cout << "*****************************************************************\n" ;
      // cout << "Avant multiplication log" << eas << endl;
      gmm::scale(eas, 0.5*log(r2));
      // cout << "Aprés multiplication log" << eas << endl;
      for (size_type i = 0; i < N; ++i)
	for (size_type j = 0; j < N; ++j)
	  t(i,j) = eas(i, j);    
      }

    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
      if (mf_data) {
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


//============================================================================================================
//=     /*************************************************************************************/              =
//=     /* Term 4 d^2(L)/d(alpha)d(u^h)   The part [*]:1/r(cos\theta, sin\theta)(v^h)_ls)    */              = Compressible
//=     /*************************************************************************************/              =
//============================================================================================================


 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term_compressible4
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, U_ls_theta;// U_ls_theta le terme U_ls multiplier par le vecteur 
    base_tensor GSigma;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
        elasticity_nonlinear_optim_term_compressible4
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N),
	gradU(NFem, N), gradU_ls(NFem, N), 
	GSigma(N, N, N, N) ,sizes_(NFem, N, NFem, N) {
      sizes_.resize(1); sizes_[0] = short_type(N);
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_u.extend_vector(U_ls_, U_ls);
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

      if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
      }
      scalar_type x = mls_x(ctx.xref());
      scalar_type y = mls_y(ctx.xref());
      scalar_type r2 = x*x+y*y;
      // Computation of tensor E
      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));


      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)=gradU(alpha, alpha)+ scalar_type(1);
       // scalar_type J = gmm::lu_det(gradU);
      

      // Computation of F.Sigma 
      AHL.sigma(E, Sigma, params);
      base_matrix eas(N,N);
      gmm::mult(gradU, Sigma, eas); // eas contient le tenseur F fois sigma
       base_small_vector V(2), W(2); V[0] = x; V[1] = y;
      gmm::mult(gmm::transposed(eas), gmm::scaled(V, 1./r2), W);
      for (size_type i = 0; i < N; ++i) t[i] = W[i];
    }

    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
      if (mf_data) {
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


//========================================================================================
//=             /*********************************************************/              =
//=             /*              Term 5 d^2(L)/d^2(alpha)                 */              = Compressible
//=             /*********************************************************/              =
//========================================================================================

 template<typename VECT1, typename VECT2>
  class elasticity_nonlinear_optim_term_compressible5
    : public getfem::nonlinear_elem_term {
    const mesh_fem &mf_u;
    std::vector<scalar_type> U, U_ls;
    const mesh_fem *mf_data;
    const VECT2 &PARAMS;
    const getfem::level_set &ls;
    size_type N, cv_old, NFem;
    const abstract_hyperelastic_law &AHL;
    base_vector params, coeff, u_ls;
    base_matrix E, Sigma, gradU, gradU_ls, gradU2_ls;
    //  gradU_ls pour contenir du/dalpha et gradU2_ls pour contenir d^2u/dalpha
    base_tensor GSigma;
    bgeot::multi_index sizes_;
    mesher_level_set mls_x, mls_y;
  public:
    elasticity_nonlinear_optim_term_compressible5
    (const mesh_fem &mf_u_, const VECT1 &U_, const VECT1 &U_ls_, 
     const mesh_fem *mf_data_, const VECT2 &PARAMS_,
     const abstract_hyperelastic_law &AHL_, const getfem::level_set &ls_) 
      : mf_u(mf_u_), U(mf_u_.nb_basic_dof()), U_ls(mf_u_.nb_basic_dof()),
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_),
	N(mf_u_.linked_mesh().dim()), NFem(mf_u_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), u_ls(N), E(N, N), Sigma(N, N),
	gradU(NFem, N), gradU_ls(NFem, N), gradU2_ls(NFem, N),
	GSigma(N, N, N, N) ,sizes_(NFem, N, NFem, N) {
      sizes_.resize(1); sizes_[0] = 1;
      cv_old = size_type(-1);
      mf_u.extend_vector(U_, U);
      mf_u.extend_vector(U_ls_, U_ls);
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

      if (cv != cv_old) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
	  cv_old = cv;
      }

      scalar_type x = mls_x(ctx.xref());
      scalar_type y = mls_y(ctx.xref());
      scalar_type r2 = x*x+y*y;
      // Computation of E
      gmm::mult(gmm::transposed(gradU), gradU, E);
      gmm::add(gradU, E);
      gmm::add(gmm::transposed(gradU), E);
      gmm::scale(E, scalar_type(0.5));

      
      for (unsigned int alpha = 0; alpha < N; ++alpha)
	gradU(alpha, alpha)=gradU(alpha, alpha)+ scalar_type(1);
      // scalar_type J = gmm::lu_det(gradU);

      // Computation of d(grad U)/d(alpha) 
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
		(U_ls, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation_grad(ctx, coeff, gradU_ls, mf_u.get_qdim()); 
      ctx.pf()->interpolation_grad(ctx, coeff, gradU2_ls, mf_u.get_qdim());
      gmm::scale(gradU_ls, log(r2) / scalar_type(2));
      gmm::scale(gradU2_ls, log(r2) * log(r2) * 0.25);
      ctx.pf()->interpolation(ctx, coeff, u_ls, mf_u.get_qdim());
      base_vector Runit(2),Runitlog(2);
      Runit[0] = x / r2; Runit[1] = y / r2;
      Runitlog[0] = (x * log(r2))/ (0.5 * r2); Runitlog[1] = (y * log(r2)) / (0.5 * r2);
      gmm::rank_one_update(gradU_ls, Runit, u_ls);        // gradU_ls contient d(grad U)/d(alpha)
      gmm::rank_one_update(gradU2_ls, Runitlog, u_ls);    // gradU2_ls contient d^2(grad U)/(d(alpha))^2 
      
      // Computation of d(E)/d(alpha) 
      
      base_matrix Ealpha(N, N);
      gmm::mult(gmm::transposed(gradU_ls), gradU, Sigma);
      gmm::add(Sigma, gmm::transposed(Sigma), Ealpha);
      gmm::scale(Ealpha, 0.5);
      
      // Computation of d2(E)/d2(alpha)
      base_matrix Ealpha2(N, N);
      gmm::mult(gmm::transposed(gradU2_ls), gradU, Sigma);
      gmm::add(Sigma, gmm::transposed(Sigma), Ealpha2);
      gmm::scale(Ealpha2, 0.5);
      
      gmm::mult(gmm::transposed(gradU_ls), gradU_ls, Sigma);
      gmm::add(Sigma, Ealpha2);

      // Computation of Sigma:d^2(E)/d^2(alpha)  
      AHL.sigma(E, Sigma, params);
      t[0]=mat_euclidean_sp(Sigma, Ealpha2);
       
      // Computation of grad_Sigma:d(E)/d(alpha)
      // scalar_type temp ;
      AHL.grad_sigma(E, GSigma, params);
      for (size_type i = 0; i < N; ++i)
	for (size_type j = 0; j < N; ++j)
	  for (size_type k = 0; k < N; ++k)
	    for (size_type l = 0; l < N; ++l)
	      t[0]  += GSigma(i,j,k,l) * Ealpha(l,k) * Ealpha(i,j);

      
    }

    virtual void prepare(fem_interpolation_context& ctx, size_type nb) {
      GMM_ASSERT1(nb != 0, "Oops");
            
      if (mf_data) {
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




  /********************************************************************************************************************/
  /********************************************************************************************************************/
  /********************************************************************************************************************/
  /********************************************************************************************************************/
  /********************************************************************************************************************/
  /********************************************************************************************************************/


 
  //=========================================================================
  //
  //  Nonlinear elasticity Brick                                       Compressible
  //
  //=========================================================================


  //=========================================================================
  //
  //  Assembling the tangent of the additional term                    Compressible
  //
  //=========================================================================


template<typename VECT0, typename VECT1, typename VECT2> 
  
  void asm_nonlinear_elasticity_optim_compressible_tangent_matrix_alpha_u
  (const VECT0 &V_, const mesh_im &mim, 
   const getfem::mesh_fem &mf_u, const VECT1 &U, const VECT2 &U_ls,
   dal::bit_vector &u_enriched_dof, const VECT2 &/*alpha*/, 
   const getfem::mesh_fem *mf_data, const VECT2 &PARAMS,
   const abstract_hyperelastic_law &AHL, const getfem::level_set &ls,
   const mesh_region &rg = mesh_region::all_convexes()) {
    VECT0 &V = const_cast<VECT0 &>(V_);
    GMM_ASSERT1(mf_u.get_qdim() >= mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    elasticity_nonlinear_optim_term_compressible2 <VECT1, VECT2>
      nterm1(mf_u, U, U_ls, mf_data, PARAMS, AHL, ls);
    elasticity_nonlinear_optim_term_compressible3 <VECT1, VECT2>
      nterm2(mf_u, U, U_ls, mf_data, PARAMS, AHL, ls);
    elasticity_nonlinear_optim_term_compressible4 <VECT1, VECT2>
      nterm3(mf_u, U, U_ls, mf_data, PARAMS, AHL, ls);

    gmm::clear(V);

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





  template<typename VECT1, typename VECT2> 
  scalar_type asm_nonlinear_elasticity_optim_compressible_tangent_matrix_alpha_alpha
  (const mesh_im &mim, 
   const getfem::mesh_fem &mf_u, const VECT1 &U, const VECT2 &U_ls,
   dal::bit_vector &/*u_enriched_dof*/, const VECT2 &/*alpha*/, 
   const getfem::mesh_fem *mf_data, const VECT2 &PARAMS,
   const abstract_hyperelastic_law &AHL, const getfem::level_set &ls,
   const mesh_region &rg = mesh_region::all_convexes()) {
    std::vector<scalar_type> V(1);
    GMM_ASSERT1(mf_u.get_qdim() >= mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");

    elasticity_nonlinear_optim_term_compressible5<VECT1, VECT2>
      nterm1(mf_u, U, U_ls, mf_data, PARAMS, AHL, ls);

    getfem::generic_assembly assem1;
    
    if (mf_data)
      assem1.set("V()+=comp(NonLin(#1,#2))(1)");
    else
      assem1.set("V()+=comp(NonLin(#1))(1)");
    assem1.push_mi(mim);
    assem1.push_mf(mf_u);
    if (mf_data) assem1.push_mf(*mf_data);
    assem1.push_nonlinear_term(&nterm1);
    assem1.push_vec(V);
    assem1.assembly(rg);
    return V[0];
  }


  //=========================================================================
  //=========================================================================
  //==       Assembling the right hand side of the additional term         == Compressible
  //=========================================================================
  //=========================================================================


  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_elasticity_optim_compressible_rhs
  (const VECT1 &R1_, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT1 &U, const VECT3 &U_ls,
   dal::bit_vector &/*u_enriched_dof*/,
   const VECT2 &/*alpha*/, const getfem::mesh_fem *mf_data, const VECT3 &PARAMS,
   const abstract_hyperelastic_law &AHL, const getfem::level_set &ls,
   const mesh_region &rg = mesh_region::all_convexes()) {
    VECT1 &R1 = const_cast<VECT1 &>(R1_);
    GMM_ASSERT1(mf_u.get_qdim() >= mf_u.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");


    GMM_ASSERT1(!mf_u.is_reduced(), "Do not work for reduced fems");

    
    
    {
      elasticity_nonlinear_optim_term_compressible1<VECT2, VECT3>
	nterm(mf_u, U, U_ls, mf_data, PARAMS, AHL, ls);
      
      getfem::generic_assembly assem;
      if (mf_data)
	assem.set("t=comp(NonLin(#1,#2))); V() += t(i)");
      else
	assem.set("t=comp(NonLin(#1)); V() += t(i)");
      // comp() to be optimized ?
      assem.push_mi(mim);
      assem.push_mf(mf_u);
      if (mf_data) assem.push_mf(*mf_data);
      assem.push_nonlinear_term(&nterm);
      assem.push_vec(R1);
      assem.assembly(rg);

    }
    

  }



}  /* end of namespace getfem.                                             */

