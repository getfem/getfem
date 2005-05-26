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
				       pintegration_method reg,
				       pintegration_method sing) 
    : mesh_im(me.linked_mesh()), mls(me), cut_im(me.linked_mesh()) {
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
    else return mesh_im::int_method_of_element(cv);
  }

  DAL_SIMPLE_KEY(special_imls_key, papprox_integration);

  void mesh_im_level_set::build_method_of_convex(size_type cv) const {
    const getfem_mesh &mesh(mls.mesh_of_convex(cv));
    if (mesh.convex_index().card() == 0) DAL_INTERNAL_ERROR("");
    base_matrix G;
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    bgeot::pgeometric_trans pgt2
      = mesh.trans_of_convex(mesh.convex_index().first_true());
    dim_type n = pgt->dim();
    papprox_integration pai = regular_simplex_pim->approx_method();
    approx_integration *new_approx = new approx_integration(pgt->convex_ref());
    base_matrix KK(n,n), CS(n,n);
    base_matrix pc(pgt2->nb_points(), n); 
    for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
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
    
    getfem::mesh_region border_faces;
    getfem::outer_faces_of_mesh(mesh, border_faces);
    for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
      vectors_to_base_matrix(G, mesh.points_of_convex(it.cv()));
      bgeot::geotrans_interpolation_context c(mesh.trans_of_convex(it.cv()),
					      pai->point(0), G);
      for (size_type j = 0; j < pai->nb_points_on_face(it.f()); ++j) {
	c.set_xref(pai->point_on_face(it.f(), j));
	new_approx->add_point(c.xreal(), pai->coeff_on_face(it.f(), j)
			     * gmm::abs(c.J()), it.f());
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
    is_adapted = true; touch();
  }

  
}  /* end of namespace getfem.                                             */



