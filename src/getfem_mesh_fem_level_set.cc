// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_fem_level_set.cc : definition of a finite element
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

#include <getfem_mesh_fem_level_set.h>

namespace getfem {
  
  void mesh_fem_level_set::receipt(const MESH_CLEAR &)
  { clear(); is_adapted = false; }
  void mesh_fem_level_set::receipt(const MESH_DELETE &)
  { clear(); is_adapted = false; }
  void mesh_fem_level_set::clear_build_methods() {
    for (size_type i = 0; i < build_methods.size(); ++i)
      del_stored_object(build_methods[i]);
    build_methods.clear();
  }
  void mesh_fem_level_set::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
    is_adapted = false;
  }
  
  mesh_fem_level_set::mesh_fem_level_set(const mesh_level_set &me,
					 const mesh_fem &mef)
    : mesh_fem(mef.linked_mesh()), mls(me), mf(mef) {
    if (mf.get_qdim() != 1)
      DAL_THROW(to_be_done_error, "base mesh_fem for mesh_fem_level_set has "
		"to be of qdim one for the moment ...");
    this->add_dependency(mls);
    is_adapted = false;
  }

  
  DAL_SIMPLE_KEY(special_mfls_key, pfem);

  void mesh_fem_level_set::build_method_of_convex(size_type cv) {
      pfem pf = new fem_level_set(index_ref_iterator
				  (dof_enrichments.begin(),
				   mf.ind_dof_of_element(cv).begin()) ,
				  mf.fem_of_element(cv), mls);
      dal::add_stored_object(new special_mfls_key(pf), pf,
			     pf->ref_convex(0),
			     pf->node_tab(0));
      build_methods.push_back(pf);
      set_finite_element(cv, pf);
  }
  
  void mesh_fem_level_set::adapt(void) {
    clear_build_methods();
    enriched_dofs.clear(); enriched_elements.clear();
    dof_enrichments.resize(0);
    dof_enrichments.resize(mf.nb_dof(), 0);

    for (size_type i = 0; i < mf.nb_dof(); ++i) {
      bgeot::mesh_convex_ind_ct ct = mf.convex_to_dof(i);
      bool touch_cut = false;
      for (bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin();
	   it != ct.end(); ++it)
	if (mls.convex_is_cut(*it)) { touch_cut = true; break; }
      

      if (touch_cut) {
	std::vector<const std::string *> zoneset;
	
	for (bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin();
	     it != ct.end(); ++it) {
	  if (mls.convex_is_cut(*it))
	    mls.merge_zonesets(zoneset, mls.zoneset_of_convex(*it));
	  else
	    mls.merge_zoneset(zoneset, mls.primary_zone_of_convex(*it));
	}
	
	cout << "number of zones for dof " << i << "  " << zoneset.size()
	     << endl;
	cout << "zoneset[0] = " << *(zoneset[0]) << endl; 

	if (zoneset.size() != 1) { // stockage dans un set et map
	  std::sort(zoneset.begin(), zoneset.end());
	  dof_enrichments[i] = &(*(enrichments.insert(zoneset).first));
	  cout << "number of zones for dof " << i << " : "
	       << zoneset.size() << endl;
	  enriched_dofs.add(i);
	  for (bgeot::mesh_convex_ind_ct::const_iterator it = ct.begin();
	       it != ct.end(); ++it) enriched_elements.add(*it);
	}
      }
    }

    cout << "Enriched convexes : " << enriched_elements << endl;
    cout << "Enriched dofs : " << enriched_dofs << endl;

    for (dal::bv_visitor i(mf.convex_index()); !i.finished(); ++i) {
      if (enriched_elements[i]) build_method_of_convex(i); else
	set_finite_element(i, mf.fem_of_element(i));
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                            */

