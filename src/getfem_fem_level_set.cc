// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_fem_level_set.cc : definition of a finite element
//           method reprensenting a discontinous field across some level sets.
// Date    : March 09, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           Julien Pommier <Julien.Pommier@insa-toulouse.fr>           
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
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

#include <getfem_fem_level_set.h>

namespace getfem {
 
#ifdef lmdfjkdf
 
  void fem_level_set::update_from_context(void) const {
    adapt(); // mettre un is_adapted ...
  }

  void adapt(void) {

    dal::bit_vector enriched_dofs, enriched_elements;
    
    for (size_type i = 0; i < mf.nb_dof(); ++i) {
      bgeot::mesh_convex_ind_ct ct = mf.convex_to_dof(i);
      bool touch_cutted = false;
      for (bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin();
	   it != ct.end(); ++it)
	if (mls.convex_is_cutted(*it)) { touch_cutted = true; break; }
      

      if (touch_cutted) {
	std::vector<const std::string *> zoneset;
	
	bgeot::mesh_convex_ind_ct ct = mf.convex_to_dof(i);
	for (bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin();
	     it != ct.end(); ++it) {
	  if (mls.convex_is_cutted(*it))
	    mls.merge_zonesets(zoneset, mls.zoneset_of_element(*it));
	  else
	    mls.merge_zonesets(zoneset, primary_zone_of_element(*it));
	}
	
	
	if (zoneset.size() != 1) { // stockage dans un map + bitset
	  cout << "number of zones for dof " << i << " : "
	       << zoneset.size() << endl;
	  enriched_dofs.add(i);
	  for (bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin();
	       it != ct.end(); ++it) enriched_elements.add(*it);
	}
      }

      

    }
    // + fabrication des méthodes éléments finis sur les convexes enrichis





  }

  size_type fem_level_set::nb_dof(size_type cv) const
  { context_check(); return elements[cv].nb_dof; }
  
  size_type fem_level_set::index_of_global_dof
  (size_type cv, size_type i) const
  { return elements[cv].inddof[i]; }
  
  bgeot::pconvex_ref fem_level_set::ref_convex(size_type cv) const
  { return mim.int_method_of_element(cv)->approx_method()->ref_convex(); }
  
  const bgeot::convex<base_node> &
  fem_level_set::node_convex(size_type cv) const {
    return *(bgeot::generic_dummy_convex_ref
	     (dim(), nb_dof(cv),
	      mim.linked_mesh().structure_of_convex(cv)->nb_faces()));
  }
  bgeot::pstored_point_tab fem_level_set::node_tab(size_type)
    const { 
    if (!pspt_valid)
      { pspt = bgeot::store_point_tab(node_tab_); pspt_valid = true; }
    return pspt;
  }
  
  void fem_level_set::base_value(const base_node &, base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  void fem_level_set::grad_base_value(const base_node &,
					 base_tensor &) const
  { DAL_THROW(internal_error, "No grad values, real only element."); }
  void fem_level_set::hess_base_value(const base_node &,
					 base_tensor &) const
  { DAL_THROW(internal_error, "No hess values, real only element."); }
  
  inline void fem_level_set::actualize_fictx(pfem pf, size_type cv,
						const base_node &ptr) const {
    if (fictx_cv != cv) {
      if (pf->need_G()) 
	bgeot::vectors_to_base_matrix
	  (G, mf.linked_mesh().points_of_convex(cv));
      fictx = fem_interpolation_context
	(mf.linked_mesh().trans_of_convex(cv), pf, base_node(), G, cv);
      fictx_cv = cv;
    }
    fictx.set_xref(ptr);
  }
  
  void fem_level_set::real_base_value(const fem_interpolation_context& c, 
					 base_tensor &t) const {
    size_type nbdof = elements[c.convex_num()].nb_dof, cv;
    mi2[1] = target_dim(); mi2[0] = nbdof;
    t.adjust_sizes(mi2);
    std::fill(t.begin(), t.end(), scalar_type(0));
    if (nbdof == 0) return;
    
    if (c.have_pgp() && 
	(&c.pgp()->get_point_tab() == 
	 &mim.int_method_of_element(c.convex_num())->approx_method()->integration_points())) { 
      gausspt_interpolation_data &gpid
	= elements[c.convex_num()].gausspt[c.ii()];
      if (gpid.flags & 1) {
	cv = gpid.elt;
	pfem pf = mf.fem_of_element(cv);
	if (gpid.flags & 2) { t = gpid.base_val; return; }
	actualize_fictx(pf, cv, gpid.ptref);
	pf->real_base_value(fictx, taux);
	for (size_type i = 0; i < pf->nb_dof(cv); ++i)
	  if (gpid.local_dof[i] != size_type(-1))
	    for (size_type j = 0; j < target_dim(); ++j)
	      t(gpid.local_dof[i],j) = taux(i, j);
	if (store_values) { gpid.base_val = t; gpid.flags |= 2; }
      }
    }
    else {
      if (find_a_point(c.xreal(), ptref, cv)) {
	pfem pf = mf.fem_of_element(cv);
	actualize_fictx(pf, cv, ptref);
	pf->real_base_value(fictx, taux);
	for (size_type i = 0; i < elements[cv].nb_dof; ++i)
	  ind_dof[elements[cv].inddof[i]] = i;
	for (size_type i = 0; i < pf->nb_dof(cv); ++i)
	  for (size_type j = 0; j < target_dim(); ++j)
	    if (ind_dof[mf.ind_dof_of_element(cv)[i]] != size_type(-1))
	      t(ind_dof[mf.ind_dof_of_element(cv)[i]], j) = taux(i, j);
	for (size_type i = 0; i < elements[cv].nb_dof; ++i)
	  ind_dof[elements[cv].inddof[i]] = size_type(-1);
      }
    }
    
  }
  
  void fem_level_set::real_grad_base_value
  (const fem_interpolation_context& c, base_tensor &t) const {
    size_type nbdof = elements[c.convex_num()].nb_dof, cv;
    mi3[2] = dim(); mi3[1] = target_dim(); mi3[0] = nbdof;
    t.adjust_sizes(mi3);
    std::fill(t.begin(), t.end(), scalar_type(0));
    if (nbdof == 0) return;
    
    if (c.have_pgp()) { 
      gausspt_interpolation_data &gpid
	= elements[c.convex_num()].gausspt[c.ii()];
      if (gpid.flags & 1) {
	cv = gpid.elt;
	pfem pf = mf.fem_of_element(cv);
	if (gpid.flags & 4) { t = gpid.grad_val; return; }
	actualize_fictx(pf, cv, gpid.ptref);
	pf->real_grad_base_value(fictx, taux);

	if (pif) {
	  pif->grad(c.xreal(), trans);
	  for (size_type i = 0; i < pf->nb_dof(cv); ++i)
	    if (gpid.local_dof[i] != size_type(-1))
	      for (size_type j = 0; j < target_dim(); ++j)
		for (size_type k = 0; k < dim(); ++k) {
		  scalar_type e(0);
		  for (size_type l = 0; l < dim(); ++l)
		    e += trans(l, k) * taux(i, j, l);
		  t(gpid.local_dof[i], j, k) = e;
		}
	}
	else {
	  for (size_type i = 0; i < pf->nb_dof(cv); ++i)
	    if (gpid.local_dof[i] != size_type(-1))
	      for (size_type j = 0; j < target_dim(); ++j)
		for (size_type k = 0; k < dim(); ++k)
		  t(gpid.local_dof[i], j, k) = taux(i, j, k);
	  if (store_values) { gpid.grad_val = t; gpid.flags |= 4; }
	}
      }
    }
    else {
      if (find_a_point(c.xreal(), ptref, cv)) {
	pfem pf = mf.fem_of_element(cv);
	actualize_fictx(pf, cv, ptref);
	fictx.set_xref(ptref);
	pf->real_grad_base_value(fictx, taux);
	for (size_type i = 0; i < nbdof; ++i)
	  ind_dof[elements[cv].inddof[i]] = i;
	if (pif) {
	  pif->grad(c.xreal(), trans);
	  for (size_type i = 0; i < pf->nb_dof(cv); ++i)
	    for (size_type j = 0; j < target_dim(); ++j)
	      for (size_type k = 0; k < dim(); ++k)
		if (ind_dof[mf.ind_dof_of_element(cv)[i]] != size_type(-1)) {
		  scalar_type e(0);
		  for (size_type l = 0; l < dim(); ++l)
		    e += trans(l, k) * taux(i, j, l);
		  t(ind_dof[mf.ind_dof_of_element(cv)[i]],j,k) = e;
		}
	}
	else {
	  for (size_type i = 0; i < pf->nb_dof(cv); ++i)
	    for (size_type j = 0; j < target_dim(); ++j)
	      for (size_type k = 0; k < dim(); ++k)
		if (ind_dof[mf.ind_dof_of_element(cv)[i]] != size_type(-1))
		  t(ind_dof[mf.ind_dof_of_element(cv)[i]],j,k) = taux(i,j,k);
	}
	  for (size_type i = 0; i < nbdof; ++i)
	    ind_dof[elements[cv].inddof[i]] = size_type(-1);
      }
    }
  }
  
  void fem_level_set::real_hess_base_value
  (const fem_interpolation_context&, base_tensor &) const
  { DAL_THROW(internal_error, "Sorry, to be done."); }
  
  
  fem_level_set::fem_level_set(const mesh_fem &mef_,
			       const mesh_level_set &mls_)
    : mf(mef_), mls(mls_) {
    if (mef.get_qdim() != 1) 
      DAL_THROW(dal::to_be_done_error,"fem_level_set do not handle qdim != 1");
    this->add_dependency(mls);
    this->add_dependency(mf);
    is_pol = is_lag = false; es_degree = 5;
    is_equiv = real_element_defined = true;
    ntarget_dim = 1;
    update_from_context();
  }

  DAL_SIMPLE_KEY(special_femlevs_key, pfem);

  pfem new_fem_level_set(const mesh_fem &mef, const mesh_level_set &mls) {
    pfem pf = new fem_level_set(mef, mls);
    dal::add_stored_object(new special_femlevs_key(pf), pf);
    return pf;
  }

#endif

}  /* end of namespace getfem.                                            */

