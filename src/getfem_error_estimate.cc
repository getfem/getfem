/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2013-2014 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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


#include <getfem/getfem_error_estimate.h>
#include <getfem/getfem_contact_and_friction_common.h>
#include <getfem/getfem_mesher.h>

namespace getfem {

#ifdef EXPERIMENTAL_PURPOSE_ONLY

 
  void error_estimate_nitsche(const mesh_im & mim,
                              const mesh_fem &mf_u,
                              const base_vector &U,
                              int GAMMAC,
                              int GAMMAN,
                              scalar_type lambda,
                              scalar_type mu,
                              scalar_type gamma0,
                              scalar_type f_coeff,
			       scalar_type vertical_force,
                              base_vector &ERR) {
   
    
   
    
    // static double lambda, mu;
    const mesh &m = mf_u.linked_mesh();
    size_type N = m.dim();
    size_type qdim = mf_u.get_qdim();
    gmm::clear(ERR);
    std::vector<scalar_type> coeff1, coeff2;
    base_matrix grad1(qdim, N), grad2(qdim, N), E(N, N), S1(N, N), S2(N, N);
    base_matrix hess1(qdim, N*N);
    base_matrix G1, G2;
    bgeot::geotrans_inv_convex gic;
    base_node xref2(N);
    base_small_vector up(N), jump(N), U1(N), sig(N), sigt(N);
    // scalar_type young_modulus = 4*mu*(lambda + mu)/(lambda+2*mu);
    scalar_type Pr, scal, sign, Un ;
    scalar_type eta1 = 0, eta2 = 0, eta3 = 0, eta4 = 0;
    
    
    scalar_type force_coeff= f_coeff; // inutile pour l'instant
    force_coeff =0;
    
    //vertical force
    
    base_small_vector F(N);
    for (unsigned ii=0; ii < N-1; ++ii)
      F[ii]=0;
      F[N-1]=-vertical_force; 

    GMM_ASSERT1(!mf_u.is_reduced(), "To be adapted");
     
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv) {
      
      bgeot::pgeometric_trans pgt1 = m.trans_of_convex(cv);
      getfem::papprox_integration pai1 = 
      get_approx_im_or_fail(mim.int_method_of_element(cv));
      getfem::pfem pf1 = mf_u.fem_of_element(cv);
      scalar_type radius = m.convex_radius_estimate(cv);
      bgeot::vectors_to_base_matrix(G1, m.points_of_convex(cv));
      coeff1.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(cv))), coeff1);
      getfem::fem_interpolation_context ctx1(pgt1, pf1, base_node(N), G1, cv);
       
      // Residual on the element
      
      for (unsigned ii=0; ii < pai1->nb_points_on_convex(); ++ii) {
        base_small_vector res(N); 
        ctx1.set_xref(pai1->point(ii));
        pf1->interpolation_hess(ctx1, coeff1, hess1, dim_type(qdim));

        for (size_type i = 0; i < N; ++i)
          for (size_type j = 0; j < N; ++j)
            res[i] += (lambda + mu) * hess1(j, i*N+j) + mu * hess1(i, j*N+j)+F[i];
	
	  ERR[cv] += radius*radius*ctx1.J()*pai1->coeff(ii)*gmm::vect_norm2_sqr(res); //norme carrée
	  eta1 += (radius*radius*ctx1.J()*pai1->coeff(ii)*gmm::vect_norm2_sqr(res));
      }    
        //if (ERR[cv] > 100)
        //cout << "Erreur en résidu sur element " << cv << " : " << ERR[cv] << endl;    
        
        
       // jump of the stress between the element ant its neighbours.
	
	for (short_type f1=0; f1 < m.structure_of_convex(cv)->nb_faces(); ++f1) {
            
        size_type cvn = m.neighbour_of_convex(cv, f1);
        if (cvn == size_type(-1)) continue;
	
        bgeot::pgeometric_trans pgt2 = m.trans_of_convex(cvn);
        getfem::pfem pf2 = mf_u.fem_of_element(cvn);
        bgeot::vectors_to_base_matrix(G2, m.points_of_convex(cvn));
        coeff2.resize(mf_u.nb_basic_dof_of_element(cvn));
        gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(cvn))), coeff2);
        getfem::fem_interpolation_context ctx2(pgt2, pf2, base_node(N), G2, cvn);
        gic.init(m.points_of_convex(cvn), pgt2);
	for (unsigned ii=0; ii < pai1->nb_points_on_face(f1); ++ii) {
          
          ctx1.set_xref(pai1->point_on_face(f1, ii));
          gmm::mult(ctx1.B(), pgt1->normals()[f1], up);
          scalar_type norm = gmm::vect_norm2(up);
          up /= norm;
          scalar_type coefficient = pai1->coeff_on_face(f1, ii) * ctx1.J() * norm; 
          pf1->interpolation_grad(ctx1, coeff1, grad1, dim_type(qdim));
          gmm::copy(grad1, E); gmm::add(gmm::transposed(grad1), E);
          gmm::scale(E, 0.5);
          scalar_type trace = gmm::mat_trace(E);
          gmm::copy(gmm::identity_matrix(), S1);
          gmm::scale(S1, lambda * trace);
          gmm::add(gmm::scaled(E, 2*mu), S1);
          bool converged;
          gic.invert(ctx1.xreal(), xref2, converged);
          GMM_ASSERT1(converged, "geometric transformation not well inverted ... !");
          ctx2.set_xref(xref2);
          pf2->interpolation_grad(ctx2, coeff2, grad2, dim_type(qdim));
          gmm::copy(grad2, E); gmm::add(gmm::transposed(grad2), E);
          gmm::scale(E, 0.5);
          trace = gmm::mat_trace(E);
          gmm::copy(gmm::identity_matrix(), S2);
          gmm::scale(S2, lambda * trace);
          gmm::add(gmm::scaled(E, 2*mu), S2);
          gmm::mult(S1, up, jump);
          gmm::mult_add(S2, gmm::scaled(up, -1.0), jump);
	  
          ERR[cv] +=radius * coefficient * gmm::vect_norm2_sqr(jump);
	  eta2 += (radius * coefficient * gmm::vect_norm2_sqr(jump));
	  

	  }
	  
      	              
      }
      
   
    };

    
       {

      int bnum = GAMMAN; 
   
      getfem::mesh_region region = m.region(bnum);
      for (getfem::mr_visitor v(region,m);  !v.finished(); ++v) {
  
	bgeot::pgeometric_trans pgt1 = m.trans_of_convex(v.cv());
	getfem::papprox_integration pai1 = 
        get_approx_im_or_fail(mim.int_method_of_element(v.cv()));
	getfem::pfem pf1 = mf_u.fem_of_element(v.cv());
	scalar_type radius = m.convex_radius_estimate(v.cv());
      
	bgeot::vectors_to_base_matrix(G1, m.points_of_convex(v.cv()));
      
	coeff1.resize(mf_u.nb_basic_dof_of_element(v.cv()));
	gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(v.cv()))), coeff1);
      
	getfem::fem_interpolation_context ctx1(pgt1, pf1, base_node(N), G1, v.cv());
	
        short_type f = v.f(); 
	   
       	for (unsigned ii=0; ii < pai1->nb_points_on_face(f); ++ii) {
      
	  ctx1.set_xref(pai1->point_on_face(f, ii));
	  gmm::mult(ctx1.B(), pgt1->normals()[f], up);
	  scalar_type norm = gmm::vect_norm2(up);
	  up /= norm;
	  scalar_type coefficient = pai1->coeff_on_face(f, ii) * ctx1.J() * norm;
	  pf1->interpolation_grad(ctx1, coeff1, grad1, dim_type(qdim));
	  gmm::copy(grad1, E); gmm::add(gmm::transposed(grad1), E);
	  gmm::scale(E, 0.5);
	  scalar_type trace = gmm::mat_trace(E);
	  gmm::copy(gmm::identity_matrix(), S1);
	  gmm::scale(S1, lambda * trace);
	  gmm::add(gmm::scaled(E, 2*mu), S1);    
	  gmm::mult(S1, up, jump);
	  
	  ERR[v.cv()] +=radius * coefficient * gmm::vect_norm2_sqr(jump);
	  eta2 += (radius * coefficient * gmm::vect_norm2_sqr(jump));
        }
        //cout << "Erreur en neumann sur " << v.cv() << " ERR[v.cv()] = " << ERR[v.cv()]-eee << endl; 

        //if (ERR[v.cv()] > 100)
         // cout << "Erreur en neumann sur element " << v.cv() << " : " << ERR[v.cv()] << endl;
        
      }
      
    };
    
        
    {
      int bnum = GAMMAC; 
      
      getfem::mesh_region region = m.region(bnum);
      for (getfem::mr_visitor v(region, m);  !v.finished(); ++v) {
        
	bgeot::pgeometric_trans pgt1 = m.trans_of_convex(v.cv());
	getfem::papprox_integration pai1 = 
        get_approx_im_or_fail(mim.int_method_of_element(v.cv()));
	getfem::pfem pf1 = mf_u.fem_of_element(v.cv());
	//scalar_type radius = m.convex_radius_estimate(v.cv());//inutile
      
	bgeot::vectors_to_base_matrix(G1, m.points_of_convex(v.cv()));
      
	coeff1.resize(mf_u.nb_basic_dof_of_element(v.cv()));
	gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_u.ind_basic_dof_of_element(v.cv()))), coeff1);
      
	getfem::fem_interpolation_context ctx1(pgt1, pf1, base_node(N), G1, v.cv());
  
        // computation of h for gamma = gamma0*h
        scalar_type emax, emin, gamma;
        gmm::condition_number(ctx1.K(),emax,emin);
        gamma = gamma0 * emax * sqrt(scalar_type(N));
	//test:gamma=radius*gamma0;
	short_type f = v.f();
        for (unsigned ii=0; ii < pai1->nb_points_on_face(f); ++ii) {
         
	  ctx1.set_xref(pai1->point_on_face(f, ii));
	  gmm::mult(ctx1.B(), pgt1->normals()[f], up);
	  scalar_type norm = gmm::vect_norm2(up);
	  up /= norm;
	  scalar_type coefficient = pai1->coeff_on_face(f, ii) * ctx1.J() * norm;
	  pf1->interpolation_grad(ctx1, coeff1, grad1, dim_type(qdim));
	  pf1->interpolation(ctx1, coeff1, U1, dim_type(qdim));
	  gmm::copy(grad1, E); gmm::add(gmm::transposed(grad1), E);
	  gmm::scale(E, 0.5);
	  scalar_type trace = gmm::mat_trace(E);
	  gmm::copy(gmm::identity_matrix(), S1);
	  gmm::scale(S1, lambda * trace);
	  gmm::add(gmm::scaled(E, 2*mu), S1);  
	  gmm::mult(S1, up, sig); // sig = sigma(u)n
	  sign = gmm::vect_sp(sig,up);// sign = sigma_n(u)
	 Un = gmm::vect_sp(U1,up);// un = u_n
 	  scal = Un-gamma*sign;
	  if (scal<0)
	    Pr = sign;
	  else
	    Pr = (scal/gamma + sign);
	  
	  ERR[v.cv()] += coefficient*sqrt(gamma) *Pr*Pr; 
	  eta4 +=  coefficient*sqrt(gamma)* Pr*Pr;
	  
	  gmm::copy(up,sigt);
	  gmm::scale(sigt, - sign);
	  gmm::add(sig,sigt);

	  ERR[v.cv()] += coefficient *gmm::vect_norm2_sqr(sigt); 
	  eta3 +=  coefficient *gmm::vect_norm2_sqr(sigt);
	}	  	  
      	//  cout << "Erreur en contact sur " << v.cv() << " ERR[v.cv()] = " << ERR[v.cv()]-ee << endl;   
        
       // if (ERR[v.cv()] > 100)
         // cout << "Erreur en contact sur element " << v.cv() << " : " << ERR[v.cv()] << endl;
        
      }
      
    }

    cout << "eta1, eta2, eta3, eta4 = " << sqrt(eta1) << endl;  
    cout <<  sqrt(eta1) << endl;  
    cout <<  sqrt(eta2) << endl;  
    cout <<  sqrt(eta3) << endl;  
    cout <<  sqrt(eta4) << endl;
    cout <<  sqrt(eta1+eta2+eta3+eta4) << endl;

  }  
  
#endif
 
}

