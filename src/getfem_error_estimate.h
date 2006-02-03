#ifndef GETFEM_ERROR_ESTIMATE
#define GETFEM_ERROR_ESTIMATE

#include <getfem_mesh_im.h>
#include <getfem_mesh_fem.h>



namespace getfem {

  namespace {
    papprox_integration get_approx_im_or_fail(pintegration_method pim) {
      if (pim->type() != IM_APPROX) 
	DAL_THROW(failure_error, "error estimate work only with "
		  "approximate integration methods");
      return pim->approx_method();
    }
  }

  template <typename VECT1, typename VECT2>
    void error_estimate(const mesh_im &mim,
			const mesh_fem &mf,
			const VECT1 &U, 
			VECT2 &err,
			const dal::bit_vector &cvlst) {
    typedef typename gmm::linalg_traits<VECT1>::value_type T;
    typedef typename gmm::number_traits<T>::magnitude_type R;

    const mesh &m = mim.linked_mesh();
    if (&m != &mf.linked_mesh() || 
	gmm::vect_size(U) < mf.nb_dof() ||
	gmm::vect_size(err) < m.convex_index().last_true()+1)
      DAL_INTERNAL_ERROR("");
    gmm::clear(err);
    std::vector<unsigned> nb_cedges(m.convex_index().last_true()+1);
    pfem pf_old = 0;
    papprox_integration pai_old = 0;
    pfem_precomp pfp = 0;
    unsigned qdim = mf.get_qdim(), N = m.dim();
    
    std::vector<T> coeff1, coeff2, gradn(qdim);
    base_vector up(N);
    gmm::dense_matrix<T> grad1(qdim, N), grad2(qdim, N);
    base_matrix G1, G2;
    bgeot::geotrans_inv_convex gic;

    for (dal::bv_visitor cv1(cvlst); !cv1.finished(); ++cv1) {
      if (!mf.convex_index().is_in(cv1) || !mim.convex_index().is_in(cv1))
	DAL_INTERNAL_ERROR("");
      bgeot::mesh_structure::ind_set neighbours;

      bgeot::pgeometric_trans pgt1 = m.trans_of_convex(cv1);
      papprox_integration pai = 
	get_approx_im_or_fail(mim.int_method_of_element(cv1));
      pfem pf1 = mf.fem_of_element(cv1);
      
      if (pf1 != pf_old || pai != pai_old) {
	pfp = fem_precomp(pf1, &pai->integration_points());
      }
      pf_old = pf1; pai_old = pai;

      bgeot::vectors_to_base_matrix(G1, m.points_of_convex(cv1));
      
      fem_interpolation_context ctx1(pgt1, pfp, 0, G1, cv1);

      coeff1.resize(mf.nb_dof_of_element(cv1));
      gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv1))), coeff1);

      for (unsigned f1=0; f1 < m.structure_of_convex(cv1)->nb_faces(); ++f1) {
	m.neighbours_of_convex(cv1, f1, neighbours);
	for (bgeot::mesh_structure::ind_set::const_iterator it = neighbours.begin();
	     it != neighbours.end(); ++it) {
	  size_type cv2 = *it;
	  if (cv2 > cv1) continue; // avoid dealing twice with the same face

	  /* look for the face number of the second convex */
	  bgeot::pgeometric_trans pgt2 = m.trans_of_convex(cv2);
	  size_type f2 = 0;
	  for (; f2 < pgt2->structure()->nb_faces(); ++f2) {
	    if (m.is_convex_face_having_points
		(cv2, f2, pgt1->structure()->nb_points_of_face(f1),
		 m.ind_points_of_face_of_convex(cv1,f1).begin()))
	      break;
	  }
	  assert(f2 < pgt2->structure()->nb_faces());

          coeff2.resize(mf.nb_dof_of_element(cv2));
          gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_dof_of_element(cv2))), coeff2);

	  bgeot::vectors_to_base_matrix(G2, m.points_of_convex(cv2));
	  pfem pf2 = mf.fem_of_element(cv2);

	  gic.init(m.points_of_convex(cv2), pgt2);
	  fem_interpolation_context ctx2(m.trans_of_convex(cv2), pf2,
					 base_node(pgt2->dim()), G2, cv2);

	  for (unsigned ii=0; ii < pai->nb_points_on_face(f1); ++ii) {
	    ctx1.set_ii(pai->ind_first_point_on_face(f1) + ii);
	    base_node xref2;
	    gic.invert(ctx1.xreal(), xref2);
	    ctx2.set_xref(xref2);
	    pf1->interpolation_grad(ctx1, coeff1, grad1, qdim);
	    pf2->interpolation_grad(ctx2, coeff2, grad2, qdim);

	    const base_matrix& B = ctx1.B();
	    scalar_type J=ctx1.J();
	    gmm::mult(B, pgt1->normals()[f1], up);
	    J /= gmm::vect_norm2(up);
	    gmm::mult(grad1, up, gradn);
	    gmm::mult_add(grad2, gmm::scaled(up, R(-1)), gradn);
	    R a = gmm::vect_norm2_sqr(gradn)
	      * pai->integration_coefficients()[ctx1.ii()] * J;
	    err[cv1] += gmm::sqrt(a); err[cv2] += gmm::sqrt(a);
	    nb_cedges[cv1]++; nb_cedges[cv2]++; 
	  }
	}
      }
    }
    for (dal::bv_visitor cv1(cvlst); !cv1.finished(); ++cv1)
      if (nb_cedges[cv1]) err[cv1] /= R(nb_cedges[cv1]);
  }
}


#endif
