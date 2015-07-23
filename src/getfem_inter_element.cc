/*===========================================================================

 Copyright (C) 1999-2015 Yves Renard

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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


#include "getfem/getfem_inter_element.h"

namespace getfem {

  interelt_boundary_integration_::interelt_boundary_integration_
    (pintegration_method pa1, pintegration_method pa2)
    : pai1(get_approx_im_or_fail(pa1)), pai2(get_approx_im_or_fail(pa2)),
      warn_msg(false) {
      GMM_ASSERT1(pai1->structure()->dim() ==  pai2->structure()->dim(),
                  "dimensions mismatch");
    indices.resize(pai1->structure()->nb_faces()
                   * pai2->structure()->nb_faces());
  }


  std::vector<size_type> &
    interelt_boundary_integration_::face_indices(short_type f1,
                                                 short_type f2) const {
    GMM_ASSERT1(f1 <= pai1->structure()->nb_faces()
                && f2 <= pai2->structure()->nb_faces(), "face number invalid");
    std::vector<size_type> &ind=indices[f1*pai2->structure()->nb_faces()+f2];
    
    if (ind.size() == 0) {
      ind.resize(pai1->nb_points_on_face(f1));
      bgeot::pconvex_ref cvr1 = pai1->ref_convex();
      bgeot::pconvex_ref cvr2 = pai2->ref_convex();
      dim_type N = dim_type(pai1->structure()->dim() - 1);
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
            GMM_WARNING2("Integration on a face between two elements with "
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
  

  void compute_on_inter_element::compute_on_face(size_type cv,short_type ff1) {
    
    f1 = ff1;
    const getfem::mesh &m(mf.linked_mesh());
    
    GMM_ASSERT1(mf.convex_index().is_in(cv) &&
                mim.convex_index().is_in(cv), "");
    bgeot::mesh_structure::ind_set neighbours;
    
    pgt1 = m.trans_of_convex(cv);
    pintegration_method pi1 = mim.int_method_of_element(cv);
    papprox_integration pai1 = get_approx_im_or_fail(pi1);
    pfem pf1 = mf.fem_of_element(cv);
    
    if (pf1 != pf1_old || pai1 != pai_old) {
      pfp1 = fem_precomp(pf1, &pai1->integration_points(), pi1);
      pf1_old = pf1; pai_old = pai1;
    }
    
    bgeot::vectors_to_base_matrix(G1, m.points_of_convex(cv));
    
    fem_interpolation_context ctx1(pgt1, pfp1, 0, G1, cv, short_type(-1));
    
    m.neighbours_of_convex(cv, f1, neighbours);
    for (bgeot::mesh_structure::ind_set::const_iterator
           it = neighbours.begin(); it != neighbours.end(); ++it) {
      size_type cv2 = *it;
      
      pintegration_method pi2 = mim.int_method_of_element(cv2);
      papprox_integration pai2 = get_approx_im_or_fail(pi2);
      pfem pf2 = mf.fem_of_element(cv2);
      
      if (pai1 != pai1_old || pai2 != pai2_old || pf2 != pf2_old) {
        pfp2 = fem_precomp(pf2, &pai2->integration_points(), pi2);
        pibi = interelt_boundary_integration
          (mim.int_method_of_element(cv),
           mim.int_method_of_element(cv2));
        pai1_old = pai1; pai2_old = pai2; pf2_old = pf2;
      }
      
      /* look for the face number of the second convex */
      pgt2 = m.trans_of_convex(cv2);
      f2 = 0;
      for (; f2 < pgt2->structure()->nb_faces(); ++f2) {
        if (m.is_convex_face_having_points
            (cv2, f2, pgt1->structure()->nb_points_of_face(f1),
             m.ind_points_of_face_of_convex(cv,f1).begin()))
          break;
      }
      assert(f2 < pgt2->structure()->nb_faces());
      
      std::vector<size_type> &ind = pibi->face_indices(f1, f2);
      
      bgeot::vectors_to_base_matrix(G2, m.points_of_convex(cv2));
      
      fem_interpolation_context ctx2(pgt2, pfp2, 0, G2, cv2,
                                     short_type(-1));
      fem_interpolation_context ctx3(pgt2, pf2,
                                     base_node(pgt2->dim()), G2, cv2,
                                     short_type(-1));
      
      for (unsigned ii=0; ii < pai1->nb_points_on_face(f1); ++ii) {
        ctx1.set_ii(pai1->ind_first_point_on_face(f1) + ii);
        
        if (ind[ii] < pai2->nb_points_on_face(f2)) {
          ctx2.set_ii(pai2->ind_first_point_on_face(f2) + ind[ii]);
          compute_on_gauss_point(ctx1, pf1, ctx2, pf2, pai1);
        }
        else {
          ctx3.set_xref(pibi->additional_point
                        (ind[ii] - pai2->nb_points_on_face(f2)));
          compute_on_gauss_point(ctx1, pf1, ctx3, pf2, pai1);
        }
      }
    }
    
    
  }
  

}

