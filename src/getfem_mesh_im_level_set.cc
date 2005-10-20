// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mesh_im_level_set.cc : adaptable integration methods
//           on convex meshes.
// Date    : February 02, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2005-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

#include <getfem_mesh_im_level_set.h>
#include <getfem_mesher.h>
#include <bgeot_kdtree.h>

namespace getfem {

  
  void mesh_im_level_set::receipt(const MESH_CLEAR &) { clear(); is_adapted = false; }
  void mesh_im_level_set::receipt(const MESH_DELETE &) { clear(); is_adapted = false; }
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
    regular_simplex_pim = reg;
    singular_simplex_pim = (sing == 0) ? reg : sing;
    this->add_dependency(mls);
    is_adapted = false;
  }

  pintegration_method 
  mesh_im_level_set::int_method_of_element(size_type cv) const {
    if (!is_adapted) adapt();
    if (cut_im.convex_index().is_in(cv)) 
      return cut_im.int_method_of_element(cv); 
    else {
      if (integrate_where == INTEGRATE_BOUNDARY)
	return getfem::im_none();

      return mesh_im::int_method_of_element(cv);
    }
  }

  DAL_SIMPLE_KEY(special_imls_key, papprox_integration);

  void mesh_im_level_set::build_method_of_convex(size_type cv) const {
    const getfem_mesh &mesh(mls.mesh_of_convex(cv));
    if (mesh.convex_index().card() == 0) DAL_INTERNAL_ERROR("");
    base_matrix G;
    base_node B;

    std::vector<mesher_level_set> mesherls0(mls.nb_level_sets());
    std::vector<mesher_level_set> mesherls1(mls.nb_level_sets());
    dal::bit_vector convexes_arein;

    std::fstream totof("totof", std::ios::out | std::ios::app);

    if (integrate_where != (INTEGRATE_INSIDE | INTEGRATE_OUTSIDE)) {
      for (unsigned i = 0; i < mls.nb_level_sets(); ++i) {
	mesherls0[i] =  mls.get_level_set(i)->mls_of_convex(cv, 0, false);
	if (mls.get_level_set(i)->has_secondary())
	  mesherls1[i] =  mls.get_level_set(i)->mls_of_convex(cv, 1, false);
      }
      for (dal::bv_visitor scv(mesh.convex_index()); !scv.finished(); ++scv) {
	B = dal::mean_value(mesh.points_of_convex(scv));
	bool isin = true;
	for (unsigned i = 0; isin && (i < mls.nb_level_sets()); ++i) {
	  isin = isin && ((integrate_where & INTEGRATE_OUTSIDE) ? ((mesherls0[i])(B) > 0)
			  : ((mesherls0[i])(B) < 0));
	  if (mls.get_level_set(i)->has_secondary())
	    isin = isin && ((integrate_where & INTEGRATE_OUTSIDE) ? ((mesherls1[i])(B) > 0)
			  : ((mesherls1[i])(B) < 0));
	}
	convexes_arein[scv] = isin;
      }
    }
    
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    bgeot::pgeometric_trans pgt2
      = mesh.trans_of_convex(mesh.convex_index().first_true());
    dim_type n = pgt->dim();
    papprox_integration pai = regular_simplex_pim->approx_method();
    approx_integration *new_approx = new approx_integration(pgt->convex_ref());
    base_matrix KK(n,n), CS(n,n);
    base_matrix pc(pgt2->nb_points(), n);

    for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
      
      if ((integrate_where != (INTEGRATE_INSIDE | INTEGRATE_OUTSIDE)) &&
	  !convexes_arein[i]) continue;

      if (integrate_where & (INTEGRATE_INSIDE | INTEGRATE_OUTSIDE)) {
		
	vectors_to_base_matrix(G, mesh.points_of_convex(i));
	bgeot::geotrans_interpolation_context c(mesh.trans_of_convex(i),
						pai->point(0), G);
	for (size_type j = 0; j < pai->nb_points_on_convex(); ++j) {
	  c.set_xref(pai->point(j));
	  pgt2->gradient(pai->point(j), pc);
	  gmm::mult(G,pc,KK);
	  scalar_type J = gmm::lu_det(KK);
	  new_approx->add_point(c.xreal(), pai->coeff(j) * gmm::abs(J));
	}
      }

      base_matrix G2;
      vectors_to_base_matrix(G2, linked_mesh().points_of_convex(cv));
      bgeot::geotrans_interpolation_context cc(linked_mesh().trans_of_convex(cv),
					       pai->point(0), G2);

      // pgt2 = mesh.trans_of_convex(i);

      for (unsigned f = 0; f < pgt2->structure()->nb_faces(); ++f) {

	short_type ff = short_type(-1);
	unsigned isin = unsigned(-1);

	if (integrate_where == INTEGRATE_BOUNDARY) {
	  for (unsigned ils = 0; ils < mls.nb_level_sets(); ++ils) {
	    bool lisin = true;
	    for (unsigned ipt = 0;
		 ipt < pgt2->structure()->nb_points_of_face(f); ++ipt) {
	      lisin = lisin && (gmm::abs((mesherls0[ils])
		     (mesh.points_of_face_of_convex(i, f)[ipt])) < 1E-7);
	    }
	    cout << "lisin = " << lisin << endl;
	    if (lisin) { isin = ils; break; }
	  }
	  if (isin ==  unsigned(-1)) continue;
	} else {
	  B = dal::mean_value(mesh.points_of_face_of_convex(i, f));
	  if (pgt->convex_ref()->is_in(B) < -1E-7) continue;
	  for (short_type fi = 0; fi < pgt->structure()->nb_faces(); ++fi)
	    if (gmm::abs(pgt->convex_ref()->is_in_face(fi, B)) < 1E-6) ff = fi;
	  if (ff == short_type(-1)) DAL_INTERNAL_ERROR("");
	}
	  
	vectors_to_base_matrix(G, mesh.points_of_convex(i));
	bgeot::geotrans_interpolation_context c(mesh.trans_of_convex(i),
						pai->point(0), G);
	
	base_small_vector un = pgt2->normals()[f], up(mesh.dim());
	gmm::mult(c.B(), un, up);
	scalar_type nup = bgeot::vect_norm2(up);
	
	for (size_type j = 0; j < pai->nb_points_on_face(f); ++j) {
	  c.set_xref(pai->point_on_face(f, j));
	  scalar_type nnup(1);
	  if (integrate_where == INTEGRATE_BOUNDARY) {
	    cc.set_xref(c.xreal());
	    mesherls0[isin].grad(c.xreal(), un);
	    un /= gmm::vect_norm2(un);
	    gmm::mult(cc.B(), un, up);
	    nnup = gmm::vect_norm2(up);
	    cout << "nnup = " << nnup << endl;
	  }
	  new_approx->add_point(c.xreal(), pai->coeff_on_face(f, j)
				* gmm::abs(c.J()) * nup * nnup, ff);

	  if (integrate_where == INTEGRATE_BOUNDARY) {
	    static double ssum = 0.0;
	    ssum += pai->coeff_on_face(f, j) * gmm::abs(c.J()) * nup * nnup;
	    cout << "add crack point " << c.xreal() << " : "
		 << pai->coeff_on_face(f, j) * gmm::abs(c.J()) * nup * nnup << " sum = " << ssum << endl;

	    cc.set_xref(c.xreal());
	    totof << cc.xreal()[0] << "\t" << cc.xreal()[1] << "\t" << cc.xreal()[2] << "\n";
	  }
	} 
      }
    }
    
    new_approx->valid_method();
    
    pintegration_method pim = new integration_method(new_approx);
    dal::add_stored_object(new special_imls_key(new_approx), pim,
			   new_approx->ref_convex(),
			   &(new_approx->integration_points()));
    build_methods.push_back(pim);
    cut_im.set_integration_method(cv, pim);
  }


  void mesh_im_level_set::adapt(void) const {
    for (dal::bv_visitor cv(linked_mesh().convex_index()); 
	 !cv.finished(); ++cv) {
      if (mls.is_convex_cut(cv)) build_method_of_convex(cv);
    }

	// + éléments non coupés intérieurs.

    is_adapted = true; touch();
  }

  
}  /* end of namespace getfem.                                             */



