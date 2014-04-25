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


  void error_estimate_nitsche(const mesh_im & mim,
                              const mesh_fem &mf_u,
                              const base_vector &U,
                              int GAMMAC,
                              int GAMMAN,
                              scalar_type lambda,
                              scalar_type mu,
                              scalar_type gamma0,
                              scalar_type f_coeff,
                              base_vector &ERR);

  

#endif


}

#endif

