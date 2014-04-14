/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 1999-2012 Yves Renard, Julien Pommier
 
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
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/


/**
   @file getfem_error_estimate.h 
   @author Yves Renard <Yves.Renard@insa-lyon.fr>
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date February 10, 2006.
   @brief Definition of a posteriori error estimates.
*/

#ifndef GETFEM_ERROR_ESTIMATE
#define GETFEM_ERROR_ESTIMATE

#include "getfem_mesh_im.h"
#include "getfem_mesh_fem.h"
#include "getfem_inter_element.h"

namespace getfem {

  
  template <typename VECT1, typename VECT2>
  class inter_element_normal_derivative_jump
    : public getfem::compute_on_inter_element {
    
  protected :

    typedef typename gmm::linalg_traits<VECT1>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;
    
    const VECT1 &U;
    VECT2 &err;

    std::vector<T> coeff1, coeff2, gradn, up;
    gmm::dense_matrix<T> grad1, grad2;

    virtual void compute_on_gauss_point
    (getfem::fem_interpolation_context ctx1, getfem::pfem pf1,
     getfem::fem_interpolation_context ctx2, getfem::pfem pf2,
     getfem::papprox_integration pai1) {
      
      size_type cv1 = ctx1.convex_num();
      size_type cv2 = ctx2.convex_num();
      
      if (cv1 > cv2) {

	unsigned qdim = mf.get_qdim(), N = mf.linked_mesh().dim();
        slice_vector_on_basic_dof_of_element(mf, U, cv1, coeff1);
	// coeff1.resize(mf.nb_basic_dof_of_element(cv1));
	// gmm::copy(gmm::sub_vector
	//	  (U, gmm::sub_index(mf.ind_basic_dof_of_element(cv1))),
	//	  coeff1);
        slice_vector_on_basic_dof_of_element(mf, U, cv2, coeff2);
	// coeff2.resize(mf.nb_basic_dof_of_element(cv2));
	// gmm::copy(gmm::sub_vector
	//	  (U, gmm::sub_index(mf.ind_basic_dof_of_element(cv2))),
	//	  coeff2);
	
	gmm::resize(grad1, qdim, N); gmm::resize(grad2, qdim, N);
	pf1->interpolation_grad(ctx1, coeff1, grad1, dim_type(qdim));
	pf2->interpolation_grad(ctx2, coeff2, grad2, dim_type(qdim));
	
	gradn.resize(qdim); up.resize(N);
	const base_matrix& B = ctx1.B();
	gmm::mult(B, pgt1->normals()[f1], up);
	scalar_type norm = gmm::vect_norm2(up);
	scalar_type J = ctx1.J() * norm;
	gmm::scale(up, R(1) / norm);
	gmm::mult(grad1, up, gradn);
	gmm::mult_add(grad2, gmm::scaled(up, R(-1)), gradn);
	R w = pai1->integration_coefficients()[ctx1.ii()];
	R a = gmm::vect_norm2_sqr(gradn) * w * J;
	err[cv1] += a; err[cv2] += a;
      }
    }


  public :
  
    inter_element_normal_derivative_jump
    (const VECT1 &UU, VECT2 &errr, const getfem::mesh_im &mmim,
     const getfem::mesh_fem &mmf)
      : compute_on_inter_element(mmim, mmf), U(UU), err(errr) {
//       GMM_ASSERT1(mf.get_qdim() <= 1,
// 		  "Vectorial elements not taken into account ... to be done");
    }

  };

  template <typename VECT1, typename VECT2>
  void error_estimate(const mesh_im &mim, const mesh_fem &mf,
		      const VECT1 &UU, VECT2 &err,
		      mesh_region rg = mesh_region::all_convexes()) {
    
    typedef typename gmm::linalg_traits<VECT1>::value_type T;
    std::vector<T> U(mf.nb_basic_dof());
    mf.extend_vector(UU, U);

    const mesh &m = mim.linked_mesh();
    GMM_ASSERT3(&m == &mf.linked_mesh() &&
		gmm::vect_size(err) >= m.convex_index().last_true()+1, "");
    rg.from_mesh(m);
    GMM_ASSERT1(rg.is_only_convexes(), "Invalid mesh region");
    mesh_region sub_rg = rg;
    m.intersect_with_mpi_region(sub_rg);
    gmm::clear(err);
    inter_element_normal_derivative_jump<std::vector<T>, VECT2> iendj(U, err, mim, mf);

    for (mr_visitor cv1(sub_rg); !cv1.finished(); ++cv1)
      for (short_type f=0; f<m.structure_of_convex(cv1.cv())->nb_faces(); ++f)
	iendj.compute_on_face(cv1.cv(), f);
    
    for (mr_visitor cv1(sub_rg); !cv1.finished(); ++cv1)
      err[cv1.cv()] *= m.convex_radius_estimate(cv1.cv());
    
    MPI_SUM_VECTOR(err);
  }




#ifdef EXPERIMENTAL_PURPOSE_ONLY



/*subprogramme*/

  template <typename VECT1, typename VECT2>
  class inter_element_stress_tensor_jump
    : public getfem::compute_on_inter_element_nitsche {
/*
 *
 *
 *
 *
     
  protected :

    typedef typename gmm::linalg_traits<VECT1>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;
    
    const VECT1 &U;
    VECT2 &err;

    std::vector<T> coeff1, coeff2, tGn, up;
    gmm::dense_matrix<T> tG1, tG2;

    virtual void compute_on_gauss_point
    (getfem::fem_interpolation_context ctx1, getfem::pfem pf1,
     getfem::fem_interpolation_context ctx2, getfem::pfem pf2,
     getfem::papprox_integration pai1) {
      
      size_type cv1 = ctx1.convex_num();
      size_type cv2 = ctx2.convex_num();
      
      if (cv1 > cv2) {

	unsigned qdim = mf.get_qdim(), N = mf.linked_mesh().dim();
        slice_vector_on_basic_dof_of_element(mf, U, cv1, coeff1);
        slice_vector_on_basic_dof_of_element(mf, U, cv2, coeff2);
	gmm::resize(tG1, qdim, N); gmm::resize(tG2, qdim, N);
	typedef typename gmm::linalg_traits<VECT1>::value_type T;
	std::vector<T> U(mf.nb_basic_dof());
	mf.extend_vector(UU, U);
	
	tGn.resize(qdim); up.resize(N);
	const mesh &m = mim.linked_mesh();
	const base_matrix& B = ctx1.B();	
	gmm::mult(B, pgt1->normals()[f1], up);
	scalar_type norm = gmm::vect_norm2(up);	
	
	m.compute_Neumann_terms(1, U, mf, UU, coeff1, up, tG1); // base element fini U ?
	m.compute_Neumann_terms(1, U, mf, UU, coeff2, up, tG2);	
	pf1->interpolation(ctx1, coeff1, tG1, dim_type(qdim));
	pf2->interpolation(ctx2, coeff2, tG2, dim_type(qdim));	
	scalar_type J = ctx1.J() * norm;
	gmm::scale(up, R(1) / norm);
	gmm::mult(tG1, up, tGn);
	gmm::mult_add(tG2, gmm::scaled(up, R(-1)), tGn);
	R w = pai1->integration_coefficients()[ctx1.ii()];
	R a = gmm::vect_norm2_sqr(tGn) * w * J;
	err[cv1] += a; err[cv2] += a;
	
      }
    }

    
    
    
    
    
    
    
*/
  public :
  
    inter_element_stress_tensor_jump
    (const VECT1 &UU, VECT2 &errr, const getfem::mesh_im &mmim,
     const getfem::mesh_fem &mmf)
      : compute_on_inter_element_nitsche(mmim, mmf), U(UU), err(errr) {
//       GMM_ASSERT1(mf.get_qdim() <= 1,
// 		  "Vectorial elements not taken into account ... to be done");
    }

  };


/*programme*/

 template <typename VECT1, typename VECT2>
  void error_estimate_nitsche(const mesh_im &mim, const mesh_fem &mf,
		      const VECT1 &UU, VECT2 &err,
		      mesh_region rg = mesh_region::all_convexes(),scalar_type GAMMAC,scalar_type GAMMAN) {
    
  /* 
    
    
    
    
    
    
    
    
    typedef typename gmm::linalg_traits<VECT1>::value_type T;
    std::vector<T> U(mf.nb_basic_dof());
    mf.extend_vector(UU, U);
    const mesh &m = mim.linked_mesh();
    GMM_ASSERT3(&m == &mf.linked_mesh() &&
		gmm::vect_size(err) >= m.convex_index().last_true()+1, "");   
    rg.from_mesh(m);
    MM_ASSERT1(rg.is_only_convexes(), "Invalid mesh region"); 
    mesh_region sub_rg = rg;
    m.intersect_with_mpi_region(sub_rg);
    const VECT2 &error = err;
    
    
    //error on inter element
    gmm::clear(err); gmm::clear(error) ;  
    inter_element_stress_tensor_jump<std::vector<T>, VECT2> iendj(U, error, mim, mf);
    
    
    for (mr_visitor cv1(sub_rg); !cv1.finished(); ++cv1)
    for (short_type f=0; f<m.structure_of_convex(cv1.cv())->nb_faces(); ++f)
    iendj.compute_on_face(cv1.cv(), f);
    
    for (mr_visitor cv1(sub_rg); !cv1.finished(); ++cv1)
    error[cv1.cv()] *= m.convex_radius_estimate(cv1.cv());  

      // v.cv(); // numéro de l'élément
      // v.f();  // numéro de la face   
     err = error; 
     
     
     // A faire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     //error on contact  
//     gmm::clear(error) ;     
     
    int bnum = GAMMAC;
    getfem::mesh_region region = m.region(bnum);
for (getfem::mr_visitor v(region,mesh);  !v.finished(); ++v) 
  for (short_type f=0; f<m.structure_of_convex(v.cv())->nb_faces(); ++f)
 	iendj.compute_on_face(v.cv(), f); 
  
for (getfem::mr_visitor v(region,mesh);  !v.finished(); ++v) 
        error[v.cv()] *= m.convex_radius_estimate(v.cv());

     err += error; 

     /*error on neumann boundary*/   
/*    int bnum = GAMMAN;
    getfem::mesh_region region = m.region(bnum);
for (getfem::mr_visitor v(region,mesh);  !v.finished(); ++v) 
  for (short_type f=0; f<m.structure_of_convex(v.cv())->nb_faces(); ++f)
 	iendj.compute_on_face(v.cv(), f); 
  
for (getfem::mr_visitor v(region,mesh);  !v.finished(); ++v) 
        error[v.cv()] *= m.convex_radius_estimate(v.cv());
         
      err += error; 

    MPI_SUM_VECTOR(err);
    
}








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#endif
}




/*


void crack_problem::error_estimate_nitsche(const plain_vector &U, plain_vector &ERR) {
  
  
  /* utile ??? */
/*  size_type N = mesh.dim();
  size_type qdim = mf_u().get_qdim();
  gmm::clear(ERR);
  std::vector<scalar_type> coeff1, coeff2;
  base_matrix grad1(qdim, N), grad2(qdim, N), E(N, N), S1(N, N), S2(N, N);
  base_matrix hess1(qdim, N*N);
  base_matrix G1, G2;
  bgeot::geotrans_inv_convex gic;
  base_node xref2(N);
  base_small_vector up(N), jump(N);
  
  GMM_ASSERT1(!mf_u().is_reduced(), "To be adapted");




    
    
    
    
    
    
    
    






  for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv) {
    
    getfem::mesher_level_set mmls = ls.mls_of_convex(cv, 0); // level wet
    bgeot::pgeometric_trans pgt1 = mesh.trans_of_convex(cv);
    getfem::papprox_integration pai1 = 
      get_approx_im_or_fail(mim.int_method_of_element(cv));
    getfem::pfem pf1 = mf_u().fem_of_element(cv);
    scalar_type radius = mesh.convex_radius_estimate(cv);

    bgeot::vectors_to_base_matrix(G1, mesh.points_of_convex(cv));

    coeff1.resize(mf_u().nb_basic_dof_of_element(cv));
    gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_u().ind_basic_dof_of_element(cv))), coeff1);

    getfem::fem_interpolation_context ctx1(pgt1, pf1, base_node(N), G1, cv);
     
    
    // Residual on the element  => à adapter

    
    for (unsigned ii=0; ii < pai1->nb_points_on_convex(); ++ii) {
      
      base_small_vector res = sol_f(pai1->point(ii));
      ctx1.set_xref(pai1->point(ii));
      pf1->interpolation_hess(ctx1, coeff1, hess1, dim_type(qdim));
      for (size_type i = 0; i < N; ++i)
	for (size_type j = 0; j < N; ++j)
	  res[i] += (lambda + mu) * hess1(j, i*N+j) + mu * hess1(i, j*N+j);
          ERR[cv] += radius*radius*ctx1.J()*pai1->coeff(ii)*gmm::vect_norm2(res);
    }

    if (ERR[cv] > 100)
      cout << "Erreur en residu sur element " << cv << " : " << ERR[cv] << endl;

    
    // Stress on the level set => DF ici à suppr
   
      
    getfem::pintegration_method pim = mimbound.int_method_of_element(cv);

    if (pim->type() == getfem::IM_APPROX) {
      getfem::papprox_integration pai_crack = pim->approx_method();
      
      base_small_vector gradls;
      for (unsigned ii=0; ii < pai_crack->nb_points(); ++ii) {
	
	ctx1.set_xref(pai_crack->point(ii));
	mmls.grad(pai_crack->point(ii), gradls);
	gradls /= gmm::vect_norm2(gradls);
	gmm::mult(ctx1.B(), gradls, up);
	scalar_type norm = gmm::vect_norm2(up);
	up /= norm;
	scalar_type coefficient = pai_crack->coeff(ii)*ctx1.J(); 
	
	for (scalar_type e = -1.0; e < 2.0; e += 2.0) {
	  
	  base_node ptref = pai_crack->point(ii) + e * 1.0E-7 * gradls;
	  if (pgt1->convex_ref()->is_in(ptref) > 0.) continue;
	  ctx1.set_xref(ptref);
	  pf1->interpolation_grad(ctx1, coeff1, grad1, dim_type(qdim));
	  gmm::copy(grad1, E); gmm::add(gmm::transposed(grad1), E);
	  gmm::scale(E, 0.5);
	  scalar_type trace = gmm::mat_trace(E);
	  gmm::copy(gmm::identity_matrix(), S1);
	  gmm::scale(S1, lambda * trace);
	  gmm::add(gmm::scaled(E, 2*mu), S1);
	  gmm::mult(S1, up, jump);
	
	  ERR[cv] += radius * coefficient * gmm::vect_norm2_sqr(jump);

	}
      }
    }

 
 // jump of the stress between the element ant its neighbours.
 
 
 for (short_type f1=0; f1 < mesh.structure_of_convex(cv)->nb_faces(); ++f1) {

      if (gmm::abs(mmls(mesh.trans_of_convex(cv)->convex_ref()->points_of_face(f1)[0])) < 1E-7 * radius) continue;

      size_type cvn = mesh.neighbour_of_convex(cv, f1);
      if (cvn == size_type(-1)) continue;
	
      bgeot::pgeometric_trans pgt2 = mesh.trans_of_convex(cvn);
      getfem::pfem pf2 = mf_u().fem_of_element(cvn);
      bgeot::vectors_to_base_matrix(G2, mesh.points_of_convex(cvn));
      coeff2.resize(mf_u().nb_basic_dof_of_element(cvn));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf_u().ind_basic_dof_of_element(cvn))), coeff2);
      getfem::fem_interpolation_context ctx2(pgt2, pf2, base_node(N), G2, cvn);
      gic.init(mesh.points_of_convex(cvn), pgt2);

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

      }
      
    }
    
  }
  




