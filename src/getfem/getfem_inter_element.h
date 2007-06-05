// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 1999-2007 Yves Renard
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


/**
   @file getfem_inter_element.h 
   @author Yves Renard <Yves.Renard@insa-lyon.fr>
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date February 10, 2006.
   @brief Definition of a posteriori error estimates.
*/

#ifndef GETFEM_ERROR_ESTIMATE
#define GETFEM_ERROR_ESTIMATE

#include "getfem_mesh_im.h"
#include "getfem_mesh_fem.h"

namespace getfem {

  papprox_integration get_approx_im_or_fail(pintegration_method pim);

  class interelt_boundary_integration_
    : virtual public dal::static_stored_object {
    
    papprox_integration pai1, pai2;
    mutable std::vector<base_node> add_points;
    mutable std::vector< std::vector<size_type> > indices;
    mutable bool warn_msg;

  public :
    
    interelt_boundary_integration_(pintegration_method pa1,
				   pintegration_method pa2);
    
    std::vector<size_type> &face_indices(size_type f1, size_type f2) const;
    const base_node &additional_point(size_type i) const
    { return add_points[i]; }

  };

  typedef boost::intrusive_ptr<const getfem::interelt_boundary_integration_>
    pinterelt_boundary_integration;

  pinterelt_boundary_integration interelt_boundary_integration
    (pintegration_method pa1, pintegration_method pa2);


  template <class BOUNDARY_COMP>
  void compute_on_inter_element(const mesh_im &mim, const mesh_fem &mf,
				size_type cv, size_type f1,
			        BOUNDARY_COMP &BC) {
    
    static pfem pf1_old = 0, pf2_old = 0;
    static papprox_integration pai_old = 0, pai1_old = 0, pai2_old = 0;
    static pfem_precomp pfp1 = 0, pfp2 = 0;
    static pinterelt_boundary_integration pibi = 0;

    base_matrix G1, G2;

    const getfem::mesh &m(mf.linked_mesh());

    GMM_ASSERT1(mf.convex_index().is_in(cv) &&
		mim.convex_index().is_in(cv), "");
    bgeot::mesh_structure::ind_set neighbours;

    bgeot::pgeometric_trans pgt1 = m.trans_of_convex(cv);
    papprox_integration pai1 =
      get_approx_im_or_fail(mim.int_method_of_element(cv));
    pfem pf1 = mf.fem_of_element(cv);
      
    if (pf1 != pf1_old || pai1 != pai_old) {
      pfp1 = fem_precomp(pf1, &pai1->integration_points());
      pf1_old = pf1; pai_old = pai1;
    }

    bgeot::vectors_to_base_matrix(G1, m.points_of_convex(cv));
      
    fem_interpolation_context ctx1(pgt1, pfp1, 0, G1, cv, size_type(-1));

    m.neighbours_of_convex(cv, f1, neighbours);
    for (bgeot::mesh_structure::ind_set::const_iterator
	   it = neighbours.begin(); it != neighbours.end(); ++it) {
      size_type cv2 = *it;

      papprox_integration pai2 = 
	get_approx_im_or_fail(mim.int_method_of_element(cv2));
      pfem pf2 = mf.fem_of_element(cv2);

      if (pai1 != pai1_old || pai2 != pai2_old || pf2 != pf2_old) {
	pfp2 = fem_precomp(pf2, &pai2->integration_points());
	pibi = interelt_boundary_integration
	  (mim.int_method_of_element(cv),
	   mim.int_method_of_element(cv2));
	pai1_old = pai1; pai2_old = pai2; pf2_old = pf2;
      }

      /* look for the face number of the second convex */
      bgeot::pgeometric_trans pgt2 = m.trans_of_convex(cv2);
      size_type f2 = 0;
      for (; f2 < pgt2->structure()->nb_faces(); ++f2) {
	if (m.is_convex_face_having_points
	    (cv2, f2, pgt1->structure()->nb_points_of_face(f1),
	     m.ind_points_of_face_of_convex(cv,f1).begin()))
	  break;
      }
      assert(f2 < pgt2->structure()->nb_faces());
      
      std::vector<size_type> &ind = pibi->face_indices(f1, f2);
      
//       coeff2.resize(mf.nb_dof_of_element(cv2));
//       gmm::copy(gmm::sub_vector(U, gmm::sub_index
// 				(mf.ind_dof_of_element(cv2))), coeff2);

      bgeot::vectors_to_base_matrix(G2, m.points_of_convex(cv2));

      fem_interpolation_context ctx2(pgt2, pfp2, 0, G2, cv2,
				     size_type(-1));
      fem_interpolation_context ctx3(pgt2, pf2,
				     base_node(pgt2->dim()), G2, cv2,
				     size_type(-1));

      for (unsigned ii=0; ii < pai1->nb_points_on_face(f1); ++ii) {
	ctx1.set_ii(pai1->ind_first_point_on_face(f1) + ii);
	
	if (ind[ii] < pai2->nb_points_on_face(f2)) {
	  ctx2.set_ii(pai2->ind_first_point_on_face(f2) + ind[ii]);
	  BC(ctx1, pf1, ctx2, pf2, pai1);
	}
	else {
	  ctx3.set_xref(pibi->additional_point
			(ind[ii] - pai2->nb_points_on_face(f2)));
	  BC(ctx1, pf1, ctx3, pf2, pai1);
	}
      }
    }
  }
  




 
}


#endif
