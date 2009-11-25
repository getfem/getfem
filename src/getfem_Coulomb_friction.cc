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


  // To be done:
  // - Fabriquer brique "au plus près"
  // - Rédiger utilisation basique
  // - Fonctions virtuelles pour update en fonctions de bords donnés
  // - Ajout d'une pénalisation optionnelle ?
  // - r doit être transformé en paramêtre du modèle.
  // - Un moyen de récuperer BN, BT, gap, r ...


  struct Coulomb_friction_brick : public virtual_brick {

    T_MATRIX BN, BT;
    typedef gmm::row_matrix<gmm::rsvector<value_type> > RT_MATRIX;
    // RT_MATRIX AUG_M; // For Hansbo augmentation.
    VECTOR gap, threshold, WT, WN, friction_coef, RLN, RLT; // to be transformed into parameters ...
    value_type r, alpha, beta; // to be transformed into parameters ...
    size_type d, nbc;

    // const mesh_fem *mf_u;
    // gmm::sub_interval SUBU, SUBN, SUBT;

    bool Tresca_version, symmetrized, contact_only, really_stationary;


    template<typename VEC> static void ball_projection(const VEC &x,
						       value_type radius) {
      value_type a = gmm::vect_norm2(x);
      if (radius <= 0) gmm::clear(const_cast<VEC&>(x));
      else if (a > radius) gmm::scale(const_cast<VEC&>(x), radius/a); 
    }

    template<class VEC, class VECR>
    static void ball_projection_grad_r(const VEC &x, value_type radius,
				       VECR &g) {
      value_type a = gmm::vect_norm2(x);
      if (radius > 0 && a >= radius)
	gmm::copy(gmm::scaled(x, value_type(1)/a), g);
      else gmm::clear(g);
    }

    template <class VEC, class MAT>
    static void ball_projection_grad(const VEC &x, double radius, MAT &g) {
      if (radius <= value_type(0)) { gmm::clear(g); return; }
      gmm::copy(gmm::identity_matrix(), g);
      value_type a = gmm::vect_norm2(x);
      if (a >= radius) { 
	gmm::scale(g, radius/a);
	// gmm::rank_one_update(g, gmm::scaled(x, -radius/(a*a*a)), x);
	for (size_type i = 0; i < x.size(); ++i)
	  for (size_type j = 0; j < x.size(); ++j)
	    g(i,j) -= radius*x[i]*x[j] / (a*a*a);
      }
    }

    void precomp(const model_real_plain_vector *u,
		 const model_real_plain_vector *lambda_n,
		 const model_real_plain_vector *lambda_t) {
      gmm::resize(RLN, gmm::mat_nrows(BN));
      gmm::resize(RLT, gmm::mat_nrows(BT));
   
      gmm::add(*lambda_n, gmm::scaled(gap, r), RLN);
      if (gmm::vect_size(WN) > 0)
	gmm::mult_add(BN, gmm::scaled(WN, -r*alpha), RLN);
      gmm::mult_add(BN, gmm::scaled(*u, -r*alpha), RLN);
      // if (gmm::mat_nrows(AUG_M) > 0)
      //   gmm::mult_add(AUG_M, gmm::scaled(*lambda_n,-r), RLN);
      if (!contact_only) {
	gmm::copy(*lambda_t, RLT);
	if (gmm::vect_size(WT) > 0)
	  gmm::mult_add(BT, gmm::scaled(WT, -r*beta), RLT);
	if (!really_stationary)
	  gmm::mult_add(BT, gmm::scaled(*u, -r*beta), RLT);
      }
    }

    // Il faut trois (ou deux) variables : u , lambda_n, lambda_t
    // Il y a pleins de termes : matl[0], vecl[0] pour u,
    //                           matl[1] pour u/lambda_n
    //                           matl[2] pour u/lambda_t
    //                           matl[3] pour lambda_n/lambda_n
    //                           matl[4] pour lambda_t/lambda_t
    //                           matl[5] pour lambda_t/lambda_n
    
    

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &vecl,
					model::real_veclist &,
					size_type region,
					build_version version) const {
      GMM_ASSERT1(/* a adapter */ matl.size() == 1,
		  "Coulomb friction brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
		  "Coulomb friction brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() >= 2 && vl.size() <= 2
		  /* a adapter */ && dl.size() >= 2 && dl.size() <= 3,
		  "Wrong number of variables for Coulomb friction brick");

      const model_real_plain_vector *u = &(md.real_variable(vl[0]));
      model_real_plain_vector *ru = &(vecl[0]);
      const model_real_plain_vector *lambda_n = &(md.real_variable(vl[1]));
      model_real_plain_vector *rlambda_n = &(vecl[1]); // bon ??????????????
      const model_real_plain_vector *lambda_t = 0;
      model_real_plain_vector *rlambda_t = &(vecl[2]); // bon ??????????????
      if (vl.size() > 2)
	lambda_t = &(md.real_variable(vl[2]));
      precomp(u, lambda_n, lambda_t);
      
      if (version --> matrix) {
       
	// BBN et MM les matrices passées en paramêtre (matl) ? : adapter.
	// Pb, accès aux lignes.
	RT_MATRIX BBN(gmm::mat_nrows(BN), gmm::mat_ncols(BN));
	// à la place de BBN, on peut utiliser une matl correspondant à une partie transposé puis faire une copie ... si c'est symétrique ...
	// RT_MATRIX MM(gmm::mat_nrows(BN), gmm::mat_nrows(BN));
	gmm::clear(matl[3]);

	gmm::copy(gmm::scaled(BN, -alpha), BBN);
	// if (gmm::mat_nrows(AUG_M) > 0)
	//  gmm::copy(gmm::scaled(AUG_M, -value_type(1)), MM);
	for (size_type i=0; i < gmm::mat_nrows(BN); ++i) {
	  if (RLN[i] >= value_type(0)) {
	    gmm::clear(BBN[i]);
	    // if (gmm::mat_nrows(AUG_M) > 0) gmm::clear(MM[i]);
	    matl[3](i, i) = -value_type(1)/r;
	  }
	}
	
	// gmm::copy(BBN, gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU));
	gmm::copy(BBN, matl[1]); // bon ?

	if (!contact_only) {
	  base_matrix pg(d-1, d-1);
	  base_vector vg(d-1);

	  RT_MATRIX BBT(gmm::mat_nrows(BT), gmm::mat_ncols(BT));
	  gmm::dense_matrix<value_type> BTi(d-1,  gmm::mat_ncols(BT)); // matrice pleine ??? peuton éviter ?
	  
	  for (size_type i=0; i < nb_contact_nodes(); ++i) {
	    gmm::sub_interval SUBI(i*(d-1), d-1);
	    gmm::sub_interval SUBJ(SUBT.first()+i*(d-1),(d-1));
	    gmm::sub_interval SUBJJ(i*(d-1),(d-1));
	    value_type th = Tresca_version ? threshold[i]
	      : - (MS.state())[SUBN.first()+i] * friction_coef[i];
	    
	    ball_projection_grad(gmm::sub_vector(RLT, SUBI), th, pg);
	    if (!really_stationary) {
	      gmm::mult(gmm::scaled(pg, -beta), 
			gmm::sub_matrix(BT, SUBI,
				      gmm::sub_interval(0,gmm::mat_ncols(BT))),
			BTi);
	      gmm::copy(BTi, gmm::sub_matrix(BBT, SUBJJ, SUBU));
	    }
	    
	    if (!Tresca_version) {
	      ball_projection_grad_r(gmm::sub_vector(RLT, SUBI), th, vg);
	      for (size_type k = 0; k < d-1; ++k)
		MS.tangent_matrix()(SUBT.first()+i*(d-1)+k, SUBN.first()+i)
		  = - friction_coef[i] * vg[k] / r;
	    }
	    for (size_type j = 0; j < d-1; ++j) pg(j,j) -= value_type(1);
	    gmm::copy(gmm::scaled(pg,value_type(1)/r), 
		      gmm::sub_matrix(MS.tangent_matrix(), SUBJ));
	  }
	  T_MATRIX BBBT(gmm::mat_nrows(BT), gmm::mat_ncols(BT));
	  gmm::copy(BBT, BBBT);
	  gmm::copy(BBBT, gmm::sub_matrix(MS.tangent_matrix(), SUBT, SUBU));
	}
	



	if (symmetrized) {
	  T_MATRIX tmp(mf_u->nb_dof(), mf_u->nb_dof());
	  
	  gmm::resize(tmp, mf_u->nb_dof(), gmm::mat_nrows(BN));
	  gmm::copy(gmm::transposed(gmm::sub_matrix(MS.tangent_matrix(),
						    SUBN, SUBU)), tmp);
	  gmm::copy(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBN));
	  gmm::resize(tmp, mf_u->nb_dof(), mf_u->nb_dof());
	  gmm::mult(gmm::transposed(gmm::scaled(BN,-r*alpha)),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBN, SUBU), tmp);
	  gmm::add(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU));
	  
	  if (!contact_only) {
	    gmm::mult(gmm::transposed(gmm::scaled(BT,-r*beta)),
		      gmm::sub_matrix(MS.tangent_matrix(), SUBT, SUBU), tmp);
	    gmm::add(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU));
	    gmm::resize(tmp, mf_u->nb_dof(), gmm::mat_nrows(BT));
	    gmm::copy(gmm::transposed(gmm::sub_matrix(MS.tangent_matrix(),
						      SUBT, SUBU)), tmp);
	    gmm::copy(tmp, gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBT));
	  }
	}
	else {
	  gmm::copy(gmm::scaled(gmm::transposed(BN), value_type(-1)),
		    gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBN));
	  if (!contact_only)
	    gmm::copy(gmm::scaled(gmm::transposed(BT), value_type(-1)), 
		      gmm::sub_matrix(MS.tangent_matrix(), SUBU, SUBT));
	}
	
	
      }
      
      if (version == RHS) {

	value_type c1(1);
	
	for (size_type i=0; i < nb_contact_nodes(); ++i) {
	  RLN[i] = std::min(value_type(0), RLN[i]);
	  if (!contact_only)
	    scalar_type radius = Tresca_version ? threshold[i]
	      : -friction_coef[i]*((*lambda_n)[i]);
	    ball_projection
	      (gmm::sub_vector(RLT, gmm::sub_interval(i*(d-1),d-1)), radius);
	}
	
	if (symmetrized) {
	  gmm::mult_add(gmm::transposed(BN), gmm::scaled(RLN, -c1), *ru);
	  if (!contact_only)
	    gmm::mult_add(gmm::transposed(BT), gmm::scaled(RLT, -c1), *ru);
	} else {
	  gmm::mult_add(gmm::transposed(BN), gmm::scaled(*lambda_n,-c1), *ru);
			gmm::sub_vector(MS.residual(), SUBU));
	  if (!contact_only)
	    gmm::mult_add(gmm::transposed(BT), gmm::scaled(*lambda_t,-c1),*ru);
	}
	
	/* residual on LN */
	gmm::add(gmm::scaled(lambda_n, -c1/r), gmm::scaled(RLN, c1/r),
		 *rlambda_n);
	
	/* residual on LT */
	if (!contact_only)
	  gmm::add(gmm::scaled(lambda_t, -c1/r), gmm::scaled(RLT, c1/r),
		   *rlambda_t);
      }
    }

    Coulomb_friction_brick(bool symmetrized_, bool contact_only_) {
      symmetrized = symmetrized_;
      contact_only = contact_only_;
      set_flags("isotropic linearized elasticity", false /* is linear*/,
		/* is symmetric */
		symmetrized && (contact_only || Tresca_version),
		false /* is coercive */, true /* is real */,
		false /* is complex */);
    }

  };

  // to be done
  size_type add_Coulomb_friction_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname1, const std::string &dataname2,
   size_type region, const std::string &dataname3) {
    pbrick pbr = new iso_lin_elasticity_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, dataname1);
    dl.push_back(dataname2);
    if (dataname3.size()) dl.push_back(dataname3);
    return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
			model::mimlist(1, &mim), region);
  }

  






}  /* end of namespace getfem.                                             */


#endif /* GETFEM_COULOMB_FRICTION_H__ */
