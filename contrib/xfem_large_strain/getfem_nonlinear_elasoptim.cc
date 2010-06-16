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
#include "getfem/getfem_nonlinear_elasticity.h"

namespace getfem {




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
    size_type N, cv_old;
    size_type NFem;
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
	mf_data(mf_data_), PARAMS(PARAMS_), ls(ls_);
	N(mf_.linked_mesh().dim()), NFem(mf_.get_qdim()), AHL(AHL_),
	params(AHL_.nb_params()), E(N, N), Sigma(N, N), gradU(NFem, N),
	tt(N, N, N, N), sizes_(NFem, N, NFem, N) {
      switch (version) {
	sizes_.resize(1); sizes_[0] = 1; break;
      }
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
      GMM_ASSERT(nb != 0, "Oops");
      if (nb == 1) {

	valP.resize(1);
	size_type cv = ctx.convex_num();
	coeff.resize(mf_p.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector
	     (P_ls, gmm::sub_index(mf_p.ind_basic_dof_of_element(cv))), coeff);
	ctx.pf()->interpolation(ctx, coeff, valP);
	

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
	size_type nb = AHL.nb_params();
	coeff.resize(mf_data->nb_basic_dof_of_element(cv)*nb);
	for (size_type i = 0; i < mf_data->nb_basic_dof_of_element(cv); ++i)
	  for (size_type k = 0; k < nb; ++k)
	    coeff[i * nb + k]
	      = PARAMS[mf_data->ind_basic_dof_of_element(cv)[i]*nb+k];
	ctx.pf()->interpolation(ctx, coeff, params, dim_type(nb));
      }
    } 
    
  };




  //=========================================================================
  //
  //  Nonlinear elasticity Brick
  //
  //=========================================================================


  template<typename VECT1, typename VECT2, typename VECT3> 
  void asm_nonlinear_elasticity_optim_rhs
  (const VECT1 &R1_, const VECT1 &R2_, const mesh_im &mim,
   const getfem::mesh_fem &mf_u, const VECT2 &U,
   const getfem::mesh_fem &mf_p, const VECT2 &P,
   const VECT2 &alpha, const VECT2 &beta,
   const getfem::mesh_fem *mf_data, const VECT3 &PARAMS,
   const abstract_hyperelastic_law &AHL, const getfem::level_set &ls,
   const mesh_region &rg = mesh_region::all_convexes()) {
    VECT1 &R1 = const_cast<VECT1 &>(R1_);
    VECT1 &R2 = const_cast<VECT1 &>(R2_);
    GMM_ASSERT1(mf.get_qdim() >= mf.linked_mesh().dim(),
		"wrong qdim for the mesh_fem");


    GMM_ASSERT1(!mf_u.is_reduced(), "Do not work for reduced fems");
    GMM_ASSERT1(!mf_p.is_reduced(), "Do not work for reduced fems");


    std::vector<scalar_type> P_ls(gmm::vect_size(P));

    int d = mf_u.dim();
    dal::bit_vector p_enriched_dof;
    for (dal::bv_visitor cv(mf_p.linked_mesh().convex_index());
	 !cv.finished(); ++cv) {
      pfem pf = mf_p.fem_of_element(cv);
      for (size_type j = 0; j< pf->nb_dof(cv); ++j)
	if (pf->dof_type()[j] == global_dof(d)) {
	  size_type dof = mf_p.ind_basic_dof_of_element(cv)[j];
	  p_enriched_dof.add(dof);
	  P_ls[dof] = P[dof];
	  
	}
    }

    
    std::vector<scalar_type> U_ls(gmm::vect_size(U));

    dal::bit_vector u_enriched_dof;
    for (dal::bv_visitor cv(mf_u.linked_mesh().convex_index());
	 !cv.finished(); ++cv) {
      pfem pf = mf_u.fem_of_element(cv);
      for (size_type j = 0; j< pf->nb_dof(cv); ++j)
	if (pf->dof_type()[j] == global_dof(d)) {
	  for (k = 0; k < d; ++k) {
	    size_type dof = mf_u.ind_basic_dof_of_element(cv)[j*d+k];
	    u_enriched_dof.add(dof);
	    U_ls[dof] = U[dof];
	  }
	}
    }

  
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


  struct nonlinear_elasticity_optim_brick : public virtual_brick {

    const abstract_hyperelastic_law &AHL;
    const getfem::level_set &ls;
    
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
		  "Nonlinear elasticity brick need a single mesh_im");
      GMM_ASSERT1(vl.size() == 4,
		  "Nonlinear elasticity brick need a single variable");
      GMM_ASSERT1(dl.size() == 1,
		  "Wrong number of data for nonlinear elasticity brick, "
                  << dl.size() << " should be 1 (vector).");
      GMM_ASSERT1(matl.size() == 7,  "Wrong number of terms for nonlinear "
		  "elasticity brick");

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));

      const model_real_plain_vector &p = md.real_variable(vl[1]);
      const mesh_fem &mf_p = *(md.pmesh_fem_of_variable(vl[1]));


      const model_real_plain_vector &alpha = md.real_variable(vl[2]);
      const model_real_plain_vector &beta = md.real_variable(vl[3]);

      const mesh_fem *mf_params = md.pmesh_fem_of_variable(dl[0]);
      const model_real_plain_vector &params = md.real_variable(dl[0]);
      const mesh_im &mim = *mims[0];

      size_type sl = gmm::vect_size(params);
      if (mf_params) sl = sl * mf_params->get_qdim() / mf_params->nb_dof();
      GMM_ASSERT1(sl == AHL.nb_params(), "Wrong number of coefficients for the "
		  "nonlinear constitutive elastic law");

      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	GMM_TRACE2("Nonlinear elasticity stiffness matrix assembly");
	asm_nonlinear_elasticity_tangent_matrix
	  (matl[0], mim, mf_u, u, mf_params, params, AHL, rg);
      }


      if (version & model::BUILD_RHS) {
	asm_nonlinear_elasticity_optim_rhs(vecl[0], vecl[3], mim, 
					   mf_u, u, mf_p, p, alpha, beta,
					   mf_params, params, AHL, ls, rg);
	gmm::scale(vecl[0], scalar_type(-1));
      }

    }


    nonlinear_elasticity_optim_brick(const abstract_hyperelastic_law &AHL_,
				     const getfem::level_set &ls_)
      : AHL(AHL_), ls(ls_) {
      set_flags("Nonlinear elasticity brick", false /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
		true /* is real */, false /* is complex */);
    }

  };
  
  //=========================================================================
  //  Add a nonlinear elasticity brick.  
  //=========================================================================

  size_type add_nonlinear_elasticity_optim_brick
  (model &md, const mesh_im &mim, const std::string &varname_u,
   const std::string &varname_p,  const std::string &varname_alpha,
   const std::string &varname_beta,
   const abstract_hyperelastic_law &AHL, const std::string &dataname_law,
   const std::string &dataname_ls, const getfem::level_set &ls,
   size_type region) {
    pbrick pbr = new nonlinear_elasticity_brick(AHL, ls);

    model::termlist tl;
    tl.push_back(model::term_description(varname_alpha, varname_u, true));
    tl.push_back(model::term_description(varname_alpha, varname_p, true));
    tl.push_back(model::term_description(varname_alpha, varname_alpha, true));
    tl.push_back(model::term_description(varname_beta, varname_u, true));
    tl.push_back(model::term_description(varname_beta, varname_p, true));
    tl.push_back(model::term_description(varname_beta, varname_alpha, true));
    tl.push_back(model::term_description(varname_beta, varname_beta, true));
    model::varnamelist dl(1, dataname);
    model::varnamelist vl(1, varname_u);
    vl.push_back(varname_p);
    vl.push_back(varname_alpha);
    vl.push_back(varname_beta);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1,&mim), region);
  }






}  /* end of namespace getfem.                                             */

