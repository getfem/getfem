// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2007 Yves Renard
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

#include "getfem/getfem_partial_fem.h"
#include "getfem/getfem_partial_mesh_fem.h"

namespace getfem {
  
  void partial_mesh_fem::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }
  void partial_mesh_fem::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }
  void partial_mesh_fem::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void partial_mesh_fem::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    is_adapted = false;
  }
  
  partial_mesh_fem::partial_mesh_fem(const mesh_fem &mef)
    : mesh_fem(mef.linked_mesh()), mf(mef) {
    GMM_ASSERT1(mf.get_qdim() == 1, "base mesh_fem for partial_mesh_fem has "
		"to be of qdim one for the moment ...");
    is_adapted = false;
  }
  
  dal::bit_vector select_dofs_from_im(const mesh_fem &mf, const mesh_im &mim,
				      unsigned P) {
    const mesh &m = mf.linked_mesh();
    unsigned N = m.dim();
    if (P == unsigned(-1)) P = N;
    base_matrix G;
    bgeot::pgeometric_trans pgt_old = 0;
    bgeot::pgeotrans_precomp pgp2 = 0;
    getfem::pfem pf_old = 0;
    getfem::pfem_precomp pfp = 0;
    papprox_integration pai1 = 0;
    
    std::vector<scalar_type> areas(mf.nb_dof()), area_supports(mf.nb_dof());
    dal::bit_vector kept_dofs;

    for (dal::bv_visitor cv(mim.convex_index()); !cv.finished(); ++cv) {
      bgeot::vectors_to_base_matrix(G, m.points_of_convex(cv));
      bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
      pintegration_method pim = mim.int_method_of_element(cv);
      if (pim == im_none()) continue;
      getfem::pfem pf = mf.fem_of_element(cv);
      GMM_ASSERT1(pim->type() == IM_APPROX,
		  "Works only with approximate integration");
      papprox_integration pai2= pim->approx_method();
      static papprox_integration pai2_old = 0;
      if (pgt_old != pgt || pai2 != pai2_old) {
	pai1 = getfem::classical_approx_im(pgt, 2)->approx_method();
      	pgp2 = bgeot::geotrans_precomp(pgt, &(pai2->integration_points()));
      }
      if (pai2 != pai2_old || pf != pf_old) {
	pf_old = pf;
	pfp = getfem::fem_precomp(pf, &(pai2->integration_points()));
      }
      pai2_old = pai2;
      pgt_old = pgt;

      bgeot::geotrans_interpolation_context c2(pgp2, 0, G);
      scalar_type area1 = convex_area_estimate(pgt, G, pai1);

      for (size_type i = 0; i < pai2->nb_points_on_convex(); ++i) {
	for (unsigned d = 0; d < pf->nb_dof(cv); ++d) {
	  // gerer Qdim eventuel ...
	  if (i == 0) areas[mf.ind_dof_of_element(cv)[d]] += area1;
	  c2.set_ii(i);
	  area_supports[mf.ind_dof_of_element(cv)[d]]
	    += pai2->coeff(i) * c2.J() * gmm::sqr(pfp->val(i)[d]);
	  //	    * ((gmm::abs(pfp->val(i)[d]) < 1e-10) ? 0.0 : 1.0);
	}
      }
    }

    for (size_type i = 0; i < mf.nb_dof(); ++i) {
      //cout << "area " << i << " : " << area_supports[i] << " : " << areas[i];
      if (area_supports[i] > pow(1e-14 * areas[i], scalar_type(P) / N))
	kept_dofs.add(i);
      // else cout << " eliminated"; cout << endl;
    }

    return kept_dofs;

  }

  DAL_SIMPLE_KEY(special_partialmf_key, pfem);
  void partial_mesh_fem::adapt(const dal::bit_vector &kept_dofs,
			       const dal::bit_vector &rejected_elt) {
    context_check();
    clear();
    
    dal::bit_vector kept_dof;
      
    for (dal::bv_visitor cv(linked_mesh().convex_index());
	 !cv.finished(); ++cv)
      if (!rejected_elt.is_in(cv)) {
	dal::bit_vector selected_dofs;
	for (size_type i = 0; i < mf.nb_dof_of_element(cv); ++i)
	  if (kept_dofs.is_in(mf.ind_dof_of_element(cv)[i])) 
	    selected_dofs.add(i);
	if (selected_dofs.card() == mf.nb_dof_of_element(cv))
	  set_finite_element(cv, mf.fem_of_element(cv));
	else if (selected_dofs.card()) {
	  assert(mf.fem_of_element(cv) != 0);
	  pfem pf = new partial_fem(mf.fem_of_element(cv), selected_dofs, cv);
	  dal::add_stored_object(new special_partialmf_key(pf), pf,
				 pf->ref_convex(0),
				 pf->node_tab(0));
	  build_methods.push_back(pf);
	  set_finite_element(cv, pf);
	}
      }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                            */

