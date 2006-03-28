// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================


#include <getfem_error_estimate.h>



namespace getfem {

  papprox_integration get_approx_im_or_fail(pintegration_method pim) {
    if (pim->type() != IM_APPROX) 
      DAL_THROW(failure_error, "error estimate work only with "
		"approximate integration methods");
    return pim->approx_method();
  }

  interelt_boundary_integration_::interelt_boundary_integration_
    (pintegration_method pa1, pintegration_method pa2)
    : pai1(get_approx_im_or_fail(pa1)), pai2(get_approx_im_or_fail(pa2)),
      warn_msg(false) {
    if (pai1->structure()->dim() !=  pai2->structure()->dim())
      DAL_THROW(failure_error, "dimensions mismatch");
    indices.resize(pai1->structure()->nb_faces()
		   * pai2->structure()->nb_faces());
  }


  std::vector<size_type> &
    interelt_boundary_integration_::face_indices(size_type f1,
						 size_type f2) const {
    if (f1 > pai1->structure()->nb_faces()
	|| f2 > pai2->structure()->nb_faces())
      DAL_THROW(invalid_argument, "face number invalid");
    std::vector<size_type> &ind=indices[f1*pai2->structure()->nb_faces()+f2];
    
    if (ind.size() == 0) {
      ind.resize(pai1->nb_points_on_face(f1));
      bgeot::pconvex_ref cvr1 = pai1->ref_convex();
      bgeot::pconvex_ref cvr2 = pai2->ref_convex();
      dim_type N = pai1->structure()->dim() - 1;
      base_matrix A1(N+1, N), A2(N+1, N), A(N, N), B(N, N+1);
      base_node o1 = (cvr1->dir_points_of_face(f1))[0];
      base_node o2 = (cvr2->dir_points_of_face(f2))[0];
      for (size_type i = 0; i <= N; ++i)
	for (size_type j = 0; j < N; ++j) {
	  A1(i, j) = (cvr1->dir_points_of_face(f1))[j+1][i] - o1[i];
	  A2(i, j) = (cvr2->dir_points_of_face(f2))[j+1][i] - o2[i];	      
	}
      
      if (N) {
	gmm::mult(gmm::transposed(A1), A1, A);
	gmm::lu_inverse(A);
	gmm::mult(A, gmm::transposed(A1), B);
      }
      
      for (size_type i = 0; i < pai1->nb_points_on_face(f1); ++i) {
	base_node pt = pai1->point_on_face(f1, i) - o1, ptaux(N);
	gmm::mult(B, pt, ptaux);
	gmm::mult(A2, ptaux, o2, pt);
	ind[i] = size_type(-1);
	for (size_type j = 0; j < pai2->nb_points_on_face(f2); ++j)
	  if (gmm::vect_dist2(pai2->point_on_face(f2, j), pt) < 1e-10)
	    { ind[i] = j; break; }
	if (ind[i] == size_type(-1)) {
	  if (!warn_msg)
	    DAL_WARNING2("Integration on a face between two elements with "
			 "non-conforming integration methods");
	  warn_msg = true;
	  ind[i] = add_points.size() + pai2->nb_points_on_face(f2);
	  add_points.push_back(pt);
	}
      }
    }
    return ind;
  }
  
  DAL_DOUBLE_KEY(intboundint_key_, pintegration_method, pintegration_method);

  pinterelt_boundary_integration interelt_boundary_integration
    (pintegration_method pa1, pintegration_method pa2) {
    
    dal::pstatic_stored_object o
      = dal::search_stored_object(intboundint_key_(pa1, pa2));
    if (o) return dal::stored_cast<interelt_boundary_integration_>(o);
    pinterelt_boundary_integration p =
      new interelt_boundary_integration_(pa1, pa2);
    
    dal::add_stored_object(new intboundint_key_(pa1, pa2), p,
			   pa1, pa2, dal::AUTODELETE_STATIC_OBJECT);
    return p;
  }
  




}

