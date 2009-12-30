// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2000-2008 Yves Renard
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


  int check_symmetry(const base_tensor &t) {
    int flags = 7; size_type N = 3;
    for (size_type n = 0; n < N; ++n)
      for (size_type m = 0; m < N; ++m)
	for (size_type l = 0; l < N; ++l)
	  for (size_type k = 0; k < N; ++k) {
	    if (gmm::abs(t(n,m,l,k) - t(l,k,n,m))>1e-10) flags &= (~1); 
	    if (gmm::abs(t(n,m,l,k) - t(m,n,l,k))>1e-10) flags &= (~2); 
	    if (gmm::abs(t(n,m,l,k) - t(n,m,k,l))>1e-10) flags &= (~4);
	  }
    return flags;
  }

  void abstract_hyperelastic_law::random_E(base_matrix &E) {
    size_type N = gmm::mat_nrows(E);
    base_matrix Phi(N,N); gmm::fill_random(Phi);
    gmm::mult(gmm::transposed(Phi),Phi,E);
    gmm::scale(E,-1.); gmm::add(gmm::identity_matrix(),E); 
    gmm::scale(E,-0.5);
  }

  void abstract_hyperelastic_law::test_derivatives
  (size_type N, scalar_type h, const base_vector& param) const {
    base_matrix E(N,N), E2(N,N), DE(N,N); 
    random_E(E); random_E(DE);
    gmm::scale(DE,h);
    gmm::add(E,DE,E2);
    
    base_matrix sigma1(N,N), sigma2(N,N);
    getfem::base_tensor tdsigma(N,N,N,N);
    base_matrix dsigma(N,N);
    gmm::copy(E,E2); gmm::add(DE,E2);
    sigma(E, sigma1, param); sigma(E2, sigma2, param);
    
    scalar_type d = strain_energy(E2, param) - strain_energy(E, param);
    scalar_type d2 = 0;
    for (size_type i=0; i < N; ++i) 
      for (size_type j=0; j < N; ++j) d2 += sigma1(i,j)*DE(i,j);
    if (gmm::abs(d-d2) > h*1e-5) 
      cout << "wrong derivative of strain_energy, d=" << d
	   << ", d2=" << d2 << "\n";
    
    grad_sigma(E,tdsigma,param);
    for (size_type i=0; i < N; ++i) {
      for (size_type j=0; j < N; ++j) {
	dsigma(i,j) = 0;
	for (size_type k=0; k < N; ++k) {
	  for (size_type m=0; m < N; ++m) {
	    dsigma(i,j) += tdsigma(i,j,k,m)*DE(k,m);
	  }
	}
	sigma2(i,j) -= sigma1(i,j);
	if (gmm::abs(dsigma(i,j) - sigma2(i,j)) > h*1e-5) {
	  cout << "wrong derivative of sigma, i=" << i << ", j=" 
	       << j << ", dsigma=" << dsigma(i,j) << ", var sigma = " 
	       << sigma2(i,j) << "\n";
	}
      }
    }
  }
    
  scalar_type SaintVenant_Kirchhoff_hyperelastic_law::strain_energy
  (const base_matrix &E, const base_vector &params) const {
    return gmm::sqr(gmm::mat_trace(E)) * params[0] / scalar_type(2)
      + gmm::mat_euclidean_norm_sqr(E) * params[1];
  }
  
  void SaintVenant_Kirchhoff_hyperelastic_law::sigma
  (const base_matrix &E, base_matrix &result,const base_vector &params) const {
    gmm::copy(gmm::identity_matrix(), result);
    gmm::scale(result, params[0] * gmm::mat_trace(E));
    gmm::add(gmm::scaled(E, 2 * params[1]), result);
  }
  void SaintVenant_Kirchhoff_hyperelastic_law::grad_sigma
  (const base_matrix &E, base_tensor &result,const base_vector &params) const {
    std::fill(result.begin(), result.end(), scalar_type(0));
    size_type N = gmm::mat_nrows(E);
    for (size_type i = 0; i < N; ++i)
      for (size_type l = 0; l < N; ++l) {
	result(i, i, l, l) = params[0];
	result(i, l, i, l) += params[1];
	  result(i, l, l, i) += params[1];
      }
  }

  scalar_type membrane_elastic_law::strain_energy
  (const base_matrix & /* E */, const base_vector & /* params */) const {
    // to be done if needed
    GMM_ASSERT1(false, "To be done");
    return 0;
  }
  
  void membrane_elastic_law::sigma
  (const base_matrix &E, base_matrix &result,const base_vector &params) const {
    // should be optimized, maybe deriving sigma from strain energy
    base_tensor tt(2,2,2,2);
    size_type N = gmm::mat_nrows(E);
    grad_sigma(E,tt,params);
    for (size_type i = 0; i < N; ++i)
      for (size_type j = 0; j < N; ++j) {
	result(i,j)=0.0;
	for (size_type k = 0; k < N; ++k)
	  for (size_type l = 0; l < N; ++l) 
	    result(i,j)+=tt(i,j,k,l)*E(k,l);
      }
    // add pretension in X' direction
    if(params[4]!=0) result(0,0)+=params[4];	
    // add pretension in Y' direction
    if(params[5]!=0) result(1,1)+=params[5];
    //	cout<<"sigma="<<result<<endl;
  }
  
  void membrane_elastic_law::grad_sigma
  (const base_matrix & /* E */, base_tensor &result,
   const base_vector &params) const {
    // to be optimized!!
    std::fill(result.begin(), result.end(), scalar_type(0));
    scalar_type poisonXY=params[0]*params[1]/params[2];	//Ex*vYX=Ey*vXY
    scalar_type Ghalf=( params[3] == 0) ? params[0]/(4*(1+params[1])) : params[3]/2;	//if no G entered, compute G=E/(2*(1+v))	to be cfmd!!
    std::fill(result.begin(), result.end(), scalar_type(0));
    result(0,0,0,0) = params[0]/(1-params[1]*poisonXY);
    // result(0,0,0,1) = 0;
    // result(0,0,1,0) = 0;
    result(0,0,1,1) = params[1]*params[0]/(1-params[1]*poisonXY);
    result(1,1,0,0) = params[1]*params[0]/(1-params[1]*poisonXY);
    // result(1,1,0,1) = 0;
    // result(1,1,1,0) = 0;
    result(1,1,1,1) = params[2]/(1-params[1]*poisonXY);
    // result(0,1,0,0) = 0;
    result(0,1,0,1) = Ghalf;
    result(0,1,1,0) = Ghalf;
    // result(0,1,1,1) = 0;
    // result(1,0,0,0) = 0;
    result(1,0,0,1) = Ghalf;
    result(1,0,1,0) = Ghalf;
    // result(1,0,1,1) = 0;
  }

  scalar_type Mooney_Rivlin_hyperelastic_law::strain_energy
  (const base_matrix &E, const base_vector &params) const {
    scalar_type C1 = params[0], C2 = params[1];
    return scalar_type(2) *
      (gmm::mat_trace(E) * (C1 + scalar_type(2)*C2)
       + C2*(gmm::sqr(gmm::mat_trace(E)) - gmm::mat_euclidean_norm_sqr(E)));
  }

  void Mooney_Rivlin_hyperelastic_law::sigma
  (const base_matrix &E, base_matrix &result,const base_vector &params) const {
    scalar_type C12 = scalar_type(2) * params[0];
    scalar_type C24 = scalar_type(4) * params[1];
    gmm::copy(gmm::identity_matrix(), result);
    gmm::scale(result, C24*(gmm::mat_trace(E)+scalar_type(1)) + C12);
    gmm::add(gmm::scaled(E, -C24), result);
  }
  void Mooney_Rivlin_hyperelastic_law::grad_sigma
  (const base_matrix &E, base_tensor &result,const base_vector &params) const {
    scalar_type C22 = scalar_type(2) * params[1];
    std::fill(result.begin(), result.end(), scalar_type(0));
    size_type N = gmm::mat_nrows(E);
    for (size_type i = 0; i < N; ++i)
      for (size_type l = 0; l < N; ++l) {
	result(i, i, l, l) = scalar_type(2) * C22;
	result(i, l, i, l) -= C22;
	result(i, l, l, i) -= C22;
      }
  }
 
  scalar_type Ciarlet_Geymonat_hyperelastic_law::strain_energy
  (const base_matrix &E, const base_vector &params) const {
    size_type N = gmm::mat_nrows(E);
    scalar_type a = params[1] + params[2] / scalar_type(2);
    scalar_type b = -(params[1] + params[2]) / scalar_type(2);
    scalar_type c = params[0]/scalar_type(4)  - b;
    scalar_type d = params[0]/scalar_type(2) + params[1];
    //scalar_type d = params[0] - scalar_type(2)*params[2] - scalar_type(4)*b;
    scalar_type e = -(scalar_type(3)*(a+b) + c);
    base_matrix C(N, N);
    gmm::copy(gmm::scaled(E, scalar_type(2)), C);
    gmm::add(gmm::identity_matrix(), C);
    scalar_type det = gmm::lu_det(C);
    return a * gmm::mat_trace(C)
      + b * (gmm::sqr(gmm::mat_trace(C)) - 
	     gmm::mat_euclidean_norm_sqr(C))/scalar_type(2)
      + c * det - d * log(det) / scalar_type(2) + e;
  }
  void Ciarlet_Geymonat_hyperelastic_law::sigma
  (const base_matrix &E, base_matrix &result,const base_vector &params) const {
    size_type N = gmm::mat_nrows(E);
    scalar_type a = params[1] + params[2] / scalar_type(2);
    scalar_type b = -(params[1] + params[2]) / scalar_type(2);
    scalar_type c = params[0]/scalar_type(4)  - b;
    scalar_type d = params[0]/scalar_type(2) + params[1]; 
    //d=params[0] - scalar_type(2)*params[2] - scalar_type(4)*b;
    base_matrix C(N, N);
    assert(gmm::abs(2*a+4*b+2*c-d)<1e-5);
    gmm::copy(gmm::scaled(E, scalar_type(2)), C);
    gmm::add(gmm::identity_matrix(), C);
    gmm::copy(gmm::identity_matrix(), result);
    gmm::scale(result, scalar_type(2) * (a + b * gmm::mat_trace(C)));
    gmm::add(gmm::scaled(C, -scalar_type(2) * b), result);
    scalar_type det = gmm::lu_inverse(C);
    gmm::add(gmm::scaled(C, scalar_type(2) * c * det - d), result);
  }
  void Ciarlet_Geymonat_hyperelastic_law::grad_sigma
  (const base_matrix &E, base_tensor &result,const base_vector &params) const {
    size_type N = gmm::mat_nrows(E);
    scalar_type b2 = -(params[1] + params[2]); // b * 2
    scalar_type c = (params[0]  - 2*b2) / scalar_type(4);
    //scalar_type d = params[0] - scalar_type(2)*params[2] - 2*b2;
    scalar_type d = params[0]/scalar_type(2) + params[1]; 
    base_matrix C(N, N);
    gmm::copy(gmm::scaled(E, scalar_type(2)), C);
    gmm::add(gmm::identity_matrix(), C);
    scalar_type det = gmm::lu_inverse(C);
    std::fill(result.begin(), result.end(), scalar_type(0));
    for (size_type i = 0; i < N; ++i)
      for (size_type j = 0; j < N; ++j) {
	result(i, i, j, j) += 2*b2;
	result(i, j, i, j) -= b2;
	result(i, j, j, i) -= b2;
	for (size_type  k = 0; k < N; ++k)
	  for (size_type  l = 0; l < N; ++l)
	    result(i, j, k, l) += 
	      (C(i, k)*C(l, j) + C(i, l)*C(k, j)) * (d-scalar_type(2)*det*c)
	      + (C(i, j) * C(k, l)) * det*c*scalar_type(4);
      }
  }


  int levi_civita(int i, int j, int k) {
    int ii=i+1;
    int jj=j+1;
    int kk=k+1;	//i,j,k from 0 to 2 !
    return static_cast<int>
      (int(- 1)*(static_cast<int>(pow(double(ii-jj),2.))%3)
       * (static_cast<int> (pow(double(ii-kk),2))%3 )
       * (static_cast<int> (pow(double(jj-kk),2))%3)
       * (pow(double(jj-(ii%3))-double(0.5),2)-double(1.25)));
  }







  //=========================================================================
  //
  //  Nonlinear elasticity Brick
  //
  //=========================================================================

  struct nonlinear_elasticity_brick : public virtual_brick {

    const abstract_hyperelastic_law &AHL;
    
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
      GMM_ASSERT1(vl.size() == 1,
		  "Nonlinear elasticity brick need a single variable");
      GMM_ASSERT1(dl.size() == 1,
		  "Wrong number of data for nonlinear elasticity brick, "
                  << dl.size() << " should be 1 (vector).");
      GMM_ASSERT1(matl.size() == 1,  "Wrong number of terms for nonlinear "
		  "elasticity brick");

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));

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
	asm_nonlinear_elasticity_rhs(vecl[0], mim,
				     mf_u, u, mf_params, params, AHL, rg);
	gmm::scale(vecl[0], scalar_type(-1));
      }

    }

    nonlinear_elasticity_brick(const abstract_hyperelastic_law &AHL_)
      : AHL(AHL_) {
      set_flags("Nonlinear elasticity brick", false /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
		true /* is real */, false /* is complex */);
    }

  };
  
  //=========================================================================
  //  Add a nonlinear elasticity brick.  
  //=========================================================================

  size_type add_nonlinear_elasticity_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const abstract_hyperelastic_law &AHL, const std::string &dataname,
   size_type region) {
    pbrick pbr = new nonlinear_elasticity_brick(AHL);

    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, dataname);
    model::varnamelist vl(1, varname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1,&mim), region);
  }

  //=========================================================================
  //  Von Mises or Tresca stress computation.  
  //=========================================================================

  void compute_Von_Mises_or_Tresca(model &md,
				   const std::string &varname, 
				   const abstract_hyperelastic_law &AHL,
				   const std::string &dataname,
				   const mesh_fem &mf_vm,
				   model_real_plain_vector &VM,
				   bool tresca) {
    GMM_ASSERT1(gmm::vect_size(VM) == mf_vm.nb_dof(),
		"The vector has not the good size");
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    const model_real_plain_vector &u = md.real_variable(varname);
    const mesh_fem *mf_params = md.pmesh_fem_of_variable(dataname);
    const model_real_plain_vector &params = md.real_variable(dataname);
    
    size_type sl = gmm::vect_size(params);
    if (mf_params) sl = sl * mf_params->get_qdim() / mf_params->nb_dof();
    GMM_ASSERT1(sl == AHL.nb_params(), "Wrong number of coefficients for "
		"the nonlinear constitutive elastic law");
    
    unsigned N = unsigned(mf_u.linked_mesh().dim());
    unsigned NP = unsigned(AHL.nb_params()), NFem = mf_u.get_qdim();
    model_real_plain_vector GRAD(mf_vm.nb_dof()*NFem*N);
    model_real_plain_vector PARAMS(mf_vm.nb_dof()*NP);
    if (mf_params) interpolation(*mf_params, mf_vm, params, PARAMS);
    compute_gradient(mf_u, mf_vm, u, GRAD);
    base_matrix E(N, N), gradphi(NFem,N),gradphit(N,NFem), Id(N, N),
      sigmahathat(N,N),aux(NFem,N), sigma(NFem,NFem),
      IdNFem(NFem, NFem);
    base_vector p(NP);
    if (!mf_params) gmm::copy(params, p);
    base_vector eig(NFem);
    base_vector ez(NFem);	// vector normal at deformed surface, (ex X ey)
    double normEz(0);	//norm of ez
    gmm::copy(gmm::identity_matrix(), Id);
    gmm::copy(gmm::identity_matrix(), IdNFem);
    for (size_type i = 0; i < mf_vm.nb_dof(); ++i) {
      gmm::resize(gradphi,NFem,N);
      std::copy(GRAD.begin()+i*NFem*N, GRAD.begin()+(i+1)*NFem*N,
		gradphit.begin());
      gmm::copy(gmm::transposed(gradphit),gradphi);
      for (unsigned int alpha = 0; alpha <N; ++alpha)
	gradphi(alpha, alpha)+=1;
      gmm::mult(gmm::transposed(gradphi), gradphi, E);
      gmm::add(gmm::scaled(Id, -scalar_type(1)), E);
      gmm::scale(E, scalar_type(1)/scalar_type(2));
      if (mf_params)
	gmm::copy(gmm::sub_vector(PARAMS, gmm::sub_interval(i*NP,NP)), p);
      AHL.sigma(E, sigmahathat, p);
      if (NFem == 3 && N == 2) {
	//jyh : compute ez, normal on deformed surface
	for (unsigned int l = 0; l <NFem; ++l)  {
	  ez[l]=0;
	  for (unsigned int m = 0; m <NFem; ++m) 
	    for (unsigned int n = 0; n <NFem; ++n){
	      ez[l]+=levi_civita(l,m,n)*gradphi(m,0)*gradphi(n,1);
	    }
	  normEz= gmm::vect_norm2(ez);
	}
	//jyh : end compute ez
      }
      gmm::mult(gradphi, sigmahathat, aux);
      gmm::mult(aux, gmm::transposed(gradphi), sigma);
      
      /* jyh : complete gradphi for virtual 3rd dim (perpendicular to
	 deformed surface, same thickness) */
      if (NFem == 3 && N == 2) {
	gmm::resize(gradphi,NFem,NFem);
	for (unsigned int ll = 0; ll <NFem; ++ll) 
	  for (unsigned int ii = 0; ii <NFem; ++ii) 
	    for (unsigned int jj = 0; jj <NFem; ++jj) 
	      gradphi(ll,2)+=(levi_civita(ll,ii,jj)*gradphi(ii,0)
			      *gradphi(jj,1))/normEz;
	//jyh : end complete graphi
      }
      
      gmm::scale(sigma, scalar_type(1) / gmm::lu_det(gradphi));
      
      if (!tresca) {
	/* von mises: norm(deviator(sigma)) */
	gmm::add(gmm::scaled(IdNFem, -gmm::mat_trace(sigma) / NFem), sigma);
	
	//jyh : von mises stress=sqrt(3/2)* norm(sigma) ?
	VM[i] = sqrt(3.0/2)*gmm::mat_euclidean_norm(sigma);
      } else {
	/* else compute the tresca criterion */
	//jyh : to be adapted for membrane if necessary
	gmm::symmetric_qr_algorithm(sigma, eig);
	std::sort(eig.begin(), eig.end());
	VM[i] = eig.back() - eig.front();
      }
    }
  }
  

  // ----------------------------------------------------------------------
  //
  // Nonlinear incompressibility brick
  //
  // ----------------------------------------------------------------------

  struct nonlinear_incompressibility_brick : public virtual_brick {
    
    virtual void asm_real_tangent_terms(const model &md, size_type,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &vecl,
					model::real_veclist &,
					size_type region,
					build_version version) const {
      
      GMM_ASSERT1(matl.size() == 2,  "Wrong number of terms for nonlinear "
		  "incompressibility brick");
      GMM_ASSERT1(dl.size() == 0, "Nonlinear incompressibility brick need no "
		  "data");
      GMM_ASSERT1(mims.size() == 1, "Nonlinear incompressibility brick need a "
		  "single mesh_im");
      GMM_ASSERT1(vl.size() == 2, "Wrong number of variables for nonlinear "
		  "incompressibility brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_p = md.mesh_fem_of_variable(vl[1]);
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const model_real_plain_vector &p = md.real_variable(vl[1]);
      const mesh_im &mim = *mims[0];
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      GMM_TRACE2("Stokes term assembly");
      gmm::clear(matl[0]);
      asm_stokes_B(matl[0], mim, mf_u, mf_p, rg);

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	gmm::clear(matl[1]);
	asm_nonlinear_incomp_tangent_matrix(matl[0], matl[1],
					    mim, mf_u, mf_p, u, p, rg);
      }

      if (version & model::BUILD_RHS) {
	asm_nonlinear_incomp_rhs(vecl[0], vecl[1], mim, mf_u, mf_p, u, p, rg);
	gmm::scale(vecl[0], scalar_type(-1));
	gmm::scale(vecl[1], scalar_type(-1));
      }

    }

    nonlinear_incompressibility_brick(void) {
      set_flags("Nonlinear incompressibility brick",
		false /* is linear*/,
		true /* is symmetric */, false /* is coercive */,
		true /* is real */, false /* is complex */);
    }


  };

  size_type add_nonlinear_incompressibility
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region) {
    pbrick pbr = new nonlinear_incompressibility_brick();
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }





}  /* end of namespace getfem.                                             */

