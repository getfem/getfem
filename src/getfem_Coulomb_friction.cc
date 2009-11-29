// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2004-2010 Yves Renard
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


#include "getfem/getfem_Coulomb_friction.h"

namespace getfem {

  //=========================================================================
  //
  //  Projection on a ball and gradient of the projection.
  //
  //=========================================================================

  template<typename VEC> static void ball_projection(const VEC &x,
						     scalar_type radius) {
    scalar_type a = gmm::vect_norm2(x);
    if (radius <= 0) gmm::clear(const_cast<VEC&>(x));
    else if (a > radius) gmm::scale(const_cast<VEC&>(x), radius/a); 
  }
  
  template<class VEC, class VECR>
  static void ball_projection_grad_r(const VEC &x, scalar_type radius,
				     VECR &g) {
    scalar_type a = gmm::vect_norm2(x);
    if (radius > 0 && a >= radius)
      gmm::copy(gmm::scaled(x, scalar_type(1)/a), g);
    else gmm::clear(g);
  }
  
  template <class VEC, class MAT>
  static void ball_projection_grad(const VEC &x, double radius, MAT &g) {
    if (radius <= scalar_type(0)) { gmm::clear(g); return; }
    gmm::copy(gmm::identity_matrix(), g);
    scalar_type a = gmm::vect_norm2(x);
    if (a >= radius) { 
      gmm::scale(g, radius/a);
      // gmm::rank_one_update(g, gmm::scaled(x, -radius/(a*a*a)), x);
      for (size_type i = 0; i < x.size(); ++i)
	for (size_type j = 0; j < x.size(); ++j)
	  g(i,j) -= radius*x[i]*x[j] / (a*a*a);
    }
  }

  //=========================================================================
  //
  //  Basic Brick (with given BN, BT, gap)
  //
  //=========================================================================

  // To be done:
  // - Rédiger utilisation basique + interface matlab/python
  //      (+ acces à BN, BT pour pouvoir les changer ...)
  // - Fonctions virtuelles pour update en fonctions de bords donnés
  // - Ajout d'une pénalisation optionnelle ?
  // - Un moyen de récuperer BN, BT, ...


  struct Coulomb_friction_brick : public virtual_brick {

    typedef gmm::row_matrix<gmm::rsvector<scalar_type> > RT_MATRIX;
    RT_MATRIX BN, BT;
    mutable RT_MATRIX BBN, BBT;
    mutable model_real_plain_vector gap, threshold, friction_coeff, alpha;
    mutable model_real_plain_vector RLN, RLT; 
    mutable scalar_type r, gamma;
    mutable bool is_init;
    bool Tresca_version, symmetrized, contact_only;
    bool really_stationary, friction_dynamic_term;

    void init_BBN_BBT(void) const {
      gmm::resize(BBN, gmm::mat_nrows(BN), gmm::mat_ncols(BN));
      gmm::resize(BBT, gmm::mat_nrows(BT), gmm::mat_ncols(BT));
      gmm::copy(BN, BBN);
      gmm::copy(BT, BBT);
      size_type nbc = gmm::mat_nrows(BN);
      size_type d = gmm::mat_nrows(BT)/nbc;
      for (size_type i = 0; i < nbc; ++i) {
	gmm::scale(gmm::mat_row(BBN, i), alpha[i]);
	if (!contact_only)
	  for (size_type k = 0; k < d; ++k)
	    gmm::scale(gmm::mat_row(BBT, d*i+k), alpha[i]);
      }
      is_init = true;
    }

    void precomp(const model_real_plain_vector &u,
		 const model_real_plain_vector &lambda_n,
		 const model_real_plain_vector &lambda_t,
		 const model_real_plain_vector &wt) const {
      gmm::resize(RLN, gmm::mat_nrows(BN));
      gmm::resize(RLT, gmm::mat_nrows(BT));

      gmm::copy(gmm::scaled(gap, r), RLN);
      for (size_type i = 0; i < gmm::mat_nrows(BN); ++i) RLN[i] *= alpha[i];
      gmm::add(lambda_n, RLN);
      gmm::mult_add(BBN, gmm::scaled(u, -r), RLN);
      if (!contact_only) {
	gmm::copy(lambda_t, RLT);
	if (friction_dynamic_term)
	  gmm::mult_add(BBT, gmm::scaled(wt, -r*gamma), RLT);
	if (!really_stationary)
	  gmm::mult_add(BBT, gmm::scaled(u, -r), RLT);
      }
    }

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &vecl,
					model::real_veclist &,
					size_type /* region */,
					build_version version) const {
      GMM_ASSERT1(matl.size() == 4 || matl.size() == 7,
		  "Coulomb friction brick has four or seven terms");
      GMM_ASSERT1(mims.size() == 0, "Coulomb friction brick need no mesh_im");
      GMM_ASSERT1(vl.size() >= 2 && vl.size() <= 2
		  /* a adapter */ && dl.size() >= 2 && dl.size() <= 3,
		  "Wrong number of variables for Coulomb friction brick");
      const scalar_type vt1(1);
      bool friction = (vl.size() > 2);
      size_type nbc = gmm::mat_nrows(BN);
      size_type d = gmm::mat_nrows(BT)/nbc;

      // Matrices and rhs to be filled
      model_real_sparse_matrix &T_u_u = matl[0], &T_u_n = matl[1];
      model_real_sparse_matrix &T_n_u = matl[2], &T_n_n = matl[3];
      model_real_sparse_matrix &T_u_t = matl[friction ? 4 : 0];
      model_real_sparse_matrix &T_t_u = matl[friction ? 5 : 0];
      model_real_sparse_matrix &T_t_t = matl[friction ? 6 : 0];
      model_real_sparse_matrix &T_t_n = matl[friction ? 7 : 0];
      model_real_plain_vector &ru = vecl[0];
      model_real_plain_vector &rlambda_n = vecl[3];
      model_real_plain_vector &rlambda_t = vecl[friction ? 6 : 0];
      
      // Variables
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const model_real_plain_vector &lambda_n = md.real_variable(vl[1]);
      const model_real_plain_vector &lambda_t
	= md.real_variable(vl[friction ? 2 : 0]);

      // Parameters
      // (order : r, gap, alpha, friction_coeff, gamma, wt, threshold)
      size_type np = 0, np_wt = 0, np_alpha = 0;
      const model_real_plain_vector &vr = md.real_variable(dl[np++]);
      GMM_ASSERT1(gmm::vect_size(vr) == 1, "Parameter r should be a scalar");
      r = vr[0];
      const model_real_plain_vector &vgap = md.real_variable(dl[np++]);
      GMM_ASSERT1(gmm::vect_size(vgap) == 1 || gmm::vect_size(vgap) == nbc,
		  "Parameter gap has a wrong size");
      gmm::resize(gap, nbc);
      if (gmm::vect_size(vgap) == 1)
	gmm::fill(gap, vgap[0]);
      else
	gmm::copy(vgap, gap);
      np_alpha = np++;
      const model_real_plain_vector &valpha = md.real_variable(dl[np_alpha]);
      GMM_ASSERT1(gmm::vect_size(valpha)== 1 || gmm::vect_size(valpha) == nbc,
		  "Parameter alpha has a wrong size");
      gmm::resize(alpha, nbc);
      if (gmm::vect_size(valpha) == 1)
	gmm::fill(alpha, valpha[0]);
      else
	gmm::copy(valpha, alpha);
      if (!contact_only) {
	const model_real_plain_vector &vfr = md.real_variable(dl[np++]);
	GMM_ASSERT1(gmm::vect_size(vfr)==1 || gmm::vect_size(vfr) == nbc,
		    "Parameter friction_coeff has a wrong size");
	gmm::resize(friction_coeff, nbc);
	if (gmm::vect_size(vfr) == 1)
	  gmm::fill(friction_coeff, vfr[0]);
	else
	  gmm::copy(vfr, friction_coeff);
	if (friction_dynamic_term) {
	  const model_real_plain_vector &vg = md.real_variable(dl[np++]);
	  GMM_ASSERT1(gmm::vect_size(vg) == 1,
		      "Parameter gamma should be a scalar");
	  gamma = vg[0];
	  np_wt = np++;
	}
	if (Tresca_version) {
	  const model_real_plain_vector &vth = md.real_variable(dl[np++]);
	  GMM_ASSERT1(gmm::vect_size(vth) == 1 || gmm::vect_size(vth) == nbc,
		      "Parameter threshold has a wrong size");
	  gmm::resize(threshold, nbc);
	  if (gmm::vect_size(vth) == 1)
	    gmm::fill(threshold, vth[0]);
	  else
	    gmm::copy(vth, threshold);
	}
      }

      // pre-computations
      if (md.is_var_newer_than_brick(dl[np_alpha], ib)) is_init = false;
      if (!is_init) init_BBN_BBT();
      precomp(u, lambda_n, lambda_t, md.real_variable(dl[np_wt]));
      
      if (version & model::BUILD_MATRIX) {
	// Unilateral contact
	gmm::clear(T_n_n); gmm::clear(T_u_u);
	gmm::copy(gmm::scaled(gmm::transposed(BBN), -vt1), T_u_n);
	for (size_type i=0; i < nbc; ++i) {
	  if (RLN[i] >= scalar_type(0)) {
	    gmm::clear(gmm::mat_col(T_u_n, i));
	    T_n_n(i, i) = -vt1/r;
	  }
	}
	gmm::copy(gmm::transposed(T_u_n), T_n_u);
      
	// Friction
	if (!contact_only) {
	  base_matrix pg(d, d);
	  base_vector vg(d);
	  gmm::clear(T_u_t); gmm::clear(T_t_n); gmm::clear(T_t_t);

	  for (size_type i=0; i < nbc; ++i) {
	    gmm::sub_interval SUBI(i*d, d);
	    scalar_type th = Tresca_version ? threshold[i]
	      : - lambda_n[i] * friction_coeff[i];
	    ball_projection_grad(gmm::sub_vector(RLT, SUBI), th, pg);
	    if (!really_stationary) {
	      for (size_type k = 0; k < d; ++k)
		gmm::copy(gmm::scaled(gmm::mat_row(BBT, i*d+k), -pg[k]),
			  gmm::mat_col(T_u_t, i*d+k));
	    }
	    
	    if (!Tresca_version) {
	      ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
	      for (size_type k = 0; k < d; ++k)
		T_t_n(i*d+k, i) = - friction_coeff[i] * vg[k] / r;
	    }
	    for (size_type k = 0; k < d; ++k) pg(k,k) -= vt1;
	    gmm::copy(gmm::scaled(pg, vt1/r), gmm::sub_matrix(T_t_t, SUBI));
	  }
	  gmm::copy(gmm::transposed(T_u_t), T_t_u);
	}

	if (symmetrized) {
	  // gmm::copy(gmm::transposed(T_n_u), T_u_n);  // already done
	  model_real_sparse_matrix tmp(gmm::mat_ncols(BN), gmm::mat_ncols(BN));
	  gmm::mult(gmm::transposed(gmm::scaled(BBN,-r)), T_n_u, tmp);
	  gmm::add(tmp, T_u_u);
	  
	  if (!contact_only) {
	    // gmm::copy(gmm::transposed(T_t_u), T_u_t);  // already done
	    gmm::mult(gmm::transposed(gmm::scaled(BBT,-r)), T_t_u, tmp);
	    gmm::add(tmp, T_u_u);
	  }
	}
	else {
	  gmm::copy(gmm::scaled(gmm::transposed(BN), -vt1), T_u_n);
	  if (!contact_only)
	    gmm::copy(gmm::scaled(gmm::transposed(BT), -vt1), T_u_t);
	}
      }
      
    if (version & model::BUILD_RHS) {
	for (size_type i=0; i < nbc; ++i) {
	  RLN[i] = std::min(scalar_type(0), RLN[i]);
	  if (!contact_only) {
	    scalar_type radius = Tresca_version ? threshold[i]
	      : -friction_coeff[i]*lambda_n[i];
	    ball_projection
	      (gmm::sub_vector(RLT, gmm::sub_interval(i*d,d)), radius);
	  }
	}
	
	if (symmetrized) {
	  gmm::mult_add(gmm::transposed(BN), RLN, ru);
	  if (!contact_only)
	    gmm::mult_add(gmm::transposed(BT), RLT, ru);
	} else {
	  gmm::mult_add(gmm::transposed(BN), lambda_n, ru);
	  if (!contact_only)
	    gmm::mult_add(gmm::transposed(BT), lambda_t, ru);
	}
	
	gmm::add(gmm::scaled(lambda_n, vt1/r), gmm::scaled(RLN, -vt1/r),
		 rlambda_n);

	if (!contact_only)
	  gmm::add(gmm::scaled(lambda_t, vt1/r), gmm::scaled(RLT, -vt1/r),
		   rlambda_t);
      }
    }

    Coulomb_friction_brick(bool symmetrized_, bool contact_only_) {
      symmetrized = symmetrized_;
      contact_only = contact_only_;
      is_init = false;
      Tresca_version = false;   // future version ...
      really_stationary = false;   // future version ...
      friction_dynamic_term = false;  // future version ...
      set_flags("isotropic linearized elasticity", false /* is linear*/,
		/* is symmetric */
		symmetrized && (contact_only || Tresca_version),
		false /* is coercive */, true /* is real */,
		false /* is complex */);
    }

    
    void set_BN(model_real_sparse_matrix &BN_) {
      gmm::resize(BN, gmm::mat_nrows(BN_), gmm::mat_ncols(BN_));
      gmm::copy(BN_, BN);
      is_init = false;
    }

    void set_BT(model_real_sparse_matrix &BT_) {
      gmm::resize(BT, gmm::mat_nrows(BT_), gmm::mat_ncols(BT_));
      gmm::copy(BT_, BT);
      is_init = false;
    }


  };
  
  size_type add_basic_contact_brick
  (model &md, const std::string &varname_u, const std::string &multname_n,
   const std::string &dataname_r, model_real_sparse_matrix &BN,
   std::string dataname_gap, std::string dataname_alpha,
   bool symmetrized) {
    Coulomb_friction_brick *pbr_=new Coulomb_friction_brick(symmetrized,true);
    pbr_->set_BN(BN);
    pbrick pbr = pbr_;

    model::termlist tl;
    tl.push_back(model::term_description(varname_u, varname_u, false));
    tl.push_back(model::term_description(varname_u, multname_n, false));
    tl.push_back(model::term_description(multname_n, varname_u, false));
    tl.push_back(model::term_description(multname_n, multname_n, false));
    model::varnamelist dl(1, dataname_r);

    if (dataname_gap.size() == 0) {
      dataname_gap = md.new_name("contact_gap_on_" + varname_u);
      md.add_initialized_fixed_size_data
	(dataname_gap, model_real_plain_vector(1, scalar_type(0)));
    }
    dl.push_back(dataname_gap);
    
    if (dataname_alpha.size() == 0) {
      dataname_alpha = md.new_name("contact_parameter_alpha_on_"+ multname_n);
      md.add_initialized_fixed_size_data
	(dataname_alpha, model_real_plain_vector(1, scalar_type(1)));
    }
    dl.push_back(dataname_alpha);

    model::varnamelist vl(1, varname_u);
    vl.push_back(multname_n);
    
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }
  
//   size_type add_Coulomb_friction_brick
//   (model &md, const mesh_im &mim, const std::string &varname,
//    const std::string &dataname1, const std::string &dataname2,
//    size_type region, const std::string &dataname3) {
//     pbrick pbr = new Coulomb_friction_brick(true, true);
//     model::termlist tl;
//     tl.push_back(model::term_description(varname, varname, true));
//     model::varnamelist dl(1, dataname1);
//     dl.push_back(dataname2);
//     if (dataname3.size()) dl.push_back(dataname3);
//     return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
// 			model::mimlist(1, &mim), region);
//   }

  






}  /* end of namespace getfem.                                             */
