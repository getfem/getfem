// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2005-2006 Yves Renard
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

#include <getfem_mesh_im_level_set.h>
#include <getfem_mesher.h>
#include <bgeot_kdtree.h>

namespace getfem {

  void mesh_im_level_set::update_from_context(void) const { 
    is_adapted = false; 
  }

  void mesh_im_level_set::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; 
  }

  void mesh_im_level_set::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; 

  }
  void mesh_im_level_set::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
    cut_im.clear();
  }
  void mesh_im_level_set::clear(void) {
    mesh_im::clear();
    clear_build_methods();
    is_adapted = false;
  }
  

  mesh_im_level_set::mesh_im_level_set(mesh_level_set &me,
				       int integrate_where_, 
				       pintegration_method reg,
				       pintegration_method sing) 
    : mesh_im(me.linked_mesh()), mls(me), cut_im(me.linked_mesh()),
      integrate_where(integrate_where_) {
    set_simplex_im(reg, sing);
    this->add_dependency(mls);
    is_adapted = false;
  }

  pintegration_method 
  mesh_im_level_set::int_method_of_element(size_type cv) const {
    if (!is_adapted) const_cast<mesh_im_level_set *>(this)->adapt();
    if (cut_im.convex_index().is_in(cv)) 
      return cut_im.int_method_of_element(cv); 
    else {
      if (ignored_im.is_in(cv)) //integrate_where == INTEGRATE_BOUNDARY)
	return getfem::im_none();

      return mesh_im::int_method_of_element(cv);
    }
  }

  DAL_SIMPLE_KEY(special_imls_key, papprox_integration);

  /* only for INTEGRATE_INSIDE or INTEGRATE_OUTSIDE */
  bool mesh_im_level_set::is_point_in_selected_area
  (const std::vector<mesher_level_set> &mesherls0,
   const std::vector<mesher_level_set> &mesherls1,
				  const base_node& P) {
    bool isin = true;
    for (unsigned i = 0; isin && (i < mls.nb_level_sets()); ++i) {
      isin = isin && ((mesherls0[i])(P) < 0);
      if (mls.get_level_set(i)->has_secondary())
	isin = isin && ((mesherls1[i])(P) < 0);
    }
    return ((integrate_where & INTEGRATE_INSIDE)) ? isin : !isin;
  }

  void mesh_im_level_set::build_method_of_convex(size_type cv) {
    const mesh &msh(mls.mesh_of_convex(cv));
    if (msh.convex_index().card() == 0) DAL_INTERNAL_ERROR("");
    base_matrix G;
    base_node B;

    std::vector<mesher_level_set> mesherls0(mls.nb_level_sets());
    std::vector<mesher_level_set> mesherls1(mls.nb_level_sets());
    dal::bit_vector convexes_arein;

    //std::fstream totof("totof", std::ios::out | std::ios::app);
    for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
      mesherls0[i] =  mls.get_level_set(i)->mls_of_convex(cv, 0, false);
      if (mls.get_level_set(i)->has_secondary())
	mesherls1[i] =  mls.get_level_set(i)->mls_of_convex(cv, 1, false);
    }

    if (integrate_where != (INTEGRATE_INSIDE | INTEGRATE_OUTSIDE)) {
      for (dal::bv_visitor scv(msh.convex_index()); !scv.finished(); ++scv) {
	B = dal::mean_value(msh.points_of_convex(scv));
	convexes_arein[scv] = is_point_in_selected_area(mesherls0, mesherls1, B);
      }
    }
    
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    bgeot::pgeometric_trans pgt2
      = msh.trans_of_convex(msh.convex_index().first_true());
    dim_type n = pgt->dim();

    if (base_singular_pim) {
      if ((n == 2 && base_singular_pim->structure()
	   != bgeot::parallelepiped_structure(2))
	  || (n == 3 && base_singular_pim->structure()
	      != bgeot::prism_structure(3)) || (n < 2) || (n > 3))
	DAL_THROW(failure_error,
	 "Base integration method for quasi polar integration not convenient");
    }


    approx_integration *new_approx = new approx_integration(pgt->convex_ref());
    base_matrix KK(n,n), CS(n,n);
    base_matrix pc(pgt2->nb_points(), n);
    std::vector<size_type> ptsing;

    for (dal::bv_visitor i(msh.convex_index()); !i.finished(); ++i) {
      papprox_integration pai = regular_simplex_pim->approx_method();
      
      if ((integrate_where != (INTEGRATE_INSIDE | INTEGRATE_OUTSIDE)) &&
	  !convexes_arein[i]) continue;
      
      if (base_singular_pim && mls.crack_tip_convexes().is_in(cv)) {
	ptsing.resize(0);
	unsigned sing_ls = unsigned(-1);
	
	// cout << "cv no " << cv << endl;

	for (unsigned ils = 0; ils < mls.nb_level_sets(); ++ils)
	  if (mls.get_level_set(ils)->has_secondary()) {
	    for (unsigned ipt = 0; ipt <= n; ++ipt) {
	      if (gmm::abs((mesherls0[ils])(msh.points_of_convex(i)[ipt])) < 1E-10
		  && gmm::abs((mesherls1[ils])(msh.points_of_convex(i)[ipt])) < 1E-10) {
		if (sing_ls == unsigned(-1)) sing_ls = ils;
		if (sing_ls != ils)
		  DAL_THROW(failure_error,
			    "Two singular point in one sub element. To be done");
		ptsing.push_back(ipt);
	      }
	    }
	  }
	assert(ptsing.size() < n);
	
	if (ptsing.size() > 0) {
	  std::stringstream sts;
	  sts << "IM_QUASI_POLAR(" << name_of_int_method(base_singular_pim) << ", " << ptsing[0];
	  if (ptsing.size() > 1) sts << ", " <<  ptsing[1];
	  sts << ")";
	  //cout << "Singular int method : " << sts.str() << endl;
	  pai = int_method_descriptor(sts.str())->approx_method();
	}
      }

      base_matrix G2;
      vectors_to_base_matrix(G2, linked_mesh().points_of_convex(cv));
      bgeot::geotrans_interpolation_context cc(linked_mesh().trans_of_convex(cv),
					       pai->point(0), G2);

      if (integrate_where & (INTEGRATE_INSIDE | INTEGRATE_OUTSIDE)) {
		
	vectors_to_base_matrix(G, msh.points_of_convex(i));
	bgeot::geotrans_interpolation_context c(msh.trans_of_convex(i),
						pai->point(0), G);

	for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	  c.set_xref(pai->point(j));
	  pgt2->poly_vector_grad(pai->point(j), pc);
	  gmm::mult(G,pc,KK);
	  scalar_type J = gmm::lu_det(KK);
	  new_approx->add_point(c.xreal(), pai->coeff(j) * gmm::abs(J));

	  /*if (integrate_where == INTEGRATE_INSIDE) {
	    cc.set_xref(c.xreal());
	    totof << cc.xreal()[0] << "\t" << cc.xreal()[1] << "\n";
	  }
	  */
	}
      }




      // pgt2 = msh.trans_of_convex(i);

      for (unsigned f = 0; f < pgt2->structure()->nb_faces(); ++f) {

	short_type ff = short_type(-1);
	unsigned isin = unsigned(-1);

	if (integrate_where == INTEGRATE_BOUNDARY) {
	  for (unsigned ils = 0; ils < mls.nb_level_sets(); ++ils) {
	    bool lisin = true;
	    for (unsigned ipt = 0;
		 ipt < pgt2->structure()->nb_points_of_face(f); ++ipt) {
	      lisin = lisin && (gmm::abs((mesherls0[ils])
		     (msh.points_of_face_of_convex(i, f)[ipt])) < 1E-7);
	    }
	    if (lisin) { isin = ils; break; }
	  }
	  if (isin ==  unsigned(-1)) continue;
	} else {
	  B = dal::mean_value(msh.points_of_face_of_convex(i, f));
	  if (pgt->convex_ref()->is_in(B) < -1E-7) continue;
	  for (short_type fi = 0; fi < pgt->structure()->nb_faces(); ++fi)
	    if (gmm::abs(pgt->convex_ref()->is_in_face(fi, B)) < 1E-6) ff = fi;
	  if (ff == short_type(-1)) DAL_INTERNAL_ERROR("");
	}
	  
	vectors_to_base_matrix(G, msh.points_of_convex(i));
	bgeot::geotrans_interpolation_context c(msh.trans_of_convex(i),
						pai->point(0), G);
	
	
	for (size_type j = 0; j < pai->nb_points_on_face(f); ++j) {
	  c.set_xref(pai->point_on_face(f, j));
	  base_small_vector un = pgt2->normals()[f], up(msh.dim());
	  gmm::mult(c.B(), un, up);
	  scalar_type nup = gmm::vect_norm2(up);

	  scalar_type nnup(1);
	  if (integrate_where == INTEGRATE_BOUNDARY) {
	    cc.set_xref(c.xreal());
	    mesherls0[isin].grad(c.xreal(), un);
	    un /= gmm::vect_norm2(un);
	    gmm::mult(cc.B(), un, up);
	    nnup = gmm::vect_norm2(up);
// 	    cout << "adding coeff " << pai->coeff_on_face(f, j)
// 	      * gmm::abs(c.J()) * nup * nnup << endl;
	  }
	  new_approx->add_point(c.xreal(), pai->coeff_on_face(f, j)
				* gmm::abs(c.J()) * nup * nnup, ff);


	  /*if (integrate_where == INTEGRATE_BOUNDARY) {
	    static double ssum = 0.0;
	    ssum += pai->coeff_on_face(f, j) * gmm::abs(c.J()) * nup * nnup;
	    cout << "add crack point " << c.xreal() << " : "
		 << pai->coeff_on_face(f, j) * gmm::abs(c.J()) * nup * nnup << " sum = " << ssum << endl;
	  }*/
	  /*if (integrate_where == INTEGRATE_BOUNDARY) {
	    cc.set_xref(c.xreal());
	    totof << cc.xreal()[0] << "\t" << cc.xreal()[1] << "\n"; // << cc.xreal()[2] << "\n";
	  }
	  */
	} 
      }
    }
    
    new_approx->valid_method();

    if (new_approx->nb_points()) {
      pintegration_method pim = new integration_method(new_approx);
      dal::add_stored_object(new special_imls_key(new_approx), pim,
			     new_approx->ref_convex(),
			     &(new_approx->integration_points()));
      build_methods.push_back(pim);
      cut_im.set_integration_method(cv, pim);
    }
    else delete new_approx;
  }


  void mesh_im_level_set::adapt(void) {
    context_check();
    clear_build_methods();
    ignored_im.clear();
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      if (mls.is_convex_cut(cv)) build_method_of_convex(cv);
      else {
	if (integrate_where == INTEGRATE_BOUNDARY) {
	  ignored_im.add(cv);
	} else if (integrate_where != (INTEGRATE_OUTSIDE|INTEGRATE_INSIDE)) {
	  /* remove convexes that are not in the integration area */
	  std::vector<mesher_level_set> mesherls0(mls.nb_level_sets());
	  std::vector<mesher_level_set> mesherls1(mls.nb_level_sets());
	  for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
	    mesherls0[i] = mls.get_level_set(i)->mls_of_convex(cv, 0, false);
	    if (mls.get_level_set(i)->has_secondary())
	      mesherls1[i] = mls.get_level_set(i)->mls_of_convex(cv, 1, false);
	  }

	  base_node B(dal::mean_value(linked_mesh().trans_of_convex(cv)->convex_ref()->points()));
	  if (!is_point_in_selected_area(mesherls0, mesherls1, B))
	    ignored_im.add(cv);
	}
      }
    }
    cerr << "mesh_im_level_set: integrate = " << integrate_where << ", ignored = " << ignored_im << "\n";
    is_adapted = true; touch();
  }

  
}  /* end of namespace getfem.                                             */



