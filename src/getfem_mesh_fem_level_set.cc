/*===========================================================================
 
 Copyright (C) 2004-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

#include "getfem/getfem_mesh_fem_level_set.h"

namespace getfem {
  
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
    xfem_index = reserve_xfem_index();
    GMM_ASSERT1(mf.get_qdim() == 1, "base mesh_fem for mesh_fem_level_set has "
		"to be of qdim one for the moment ...");
    this->add_dependency(mls);
    is_adapted = false;
  }

  
  DAL_SIMPLE_KEY(special_mfls_key, pfem);

  void mesh_fem_level_set::build_method_of_convex(size_type cv) {
    pfem pf = new fem_level_set(gmm::index_ref_iterator
				(dof_enrichments.begin(),
				 mf.ind_basic_dof_of_element(cv).begin()) ,
				mf.fem_of_element(cv), mls, xfem_index);
    special_mfls_key *psm = new special_mfls_key(pf);
    dal::add_stored_object(psm, pf, pf->ref_convex(0), pf->node_tab(0));
    build_methods.push_back(pf);
    set_finite_element(cv, pf);
  }
  
  void mesh_fem_level_set::adapt(void) {
    context_check();
    clear();
    GMM_ASSERT1(!mf.is_reduced(), "Mesh fem level set not defined for reduced "
		"mesh fems (difficult or impossible to define it in the "
		"general case)");
    enriched_dofs.clear(); enriched_elements.clear();
    dof_enrichments.resize(0);
    dof_enrichments.resize(mf.nb_basic_dof(), 0);

    for (size_type i = 0; i < mf.nb_basic_dof(); ++i) {
      const mesh::ind_cv_ct &ct = mf.convex_to_basic_dof(i);
      bool touch_cut = false;
      for (mesh::ind_cv_ct::const_iterator it = ct.begin();
	   it != ct.end(); ++it)
	if (mls.is_convex_cut(*it)) { touch_cut = true; break; }
      

      if (touch_cut) {
	mesh_level_set::zoneset zones;
	
	for (mesh::ind_cv_ct::const_iterator it = ct.begin();
	     it != ct.end(); ++it) {
	  if (mls.is_convex_cut(*it)) {
	    mls.merge_zoneset(zones, mls.zoneset_of_convex(*it));
	  }
	  else {
	    mls.merge_zoneset(zones, mls.primary_zone_of_convex(*it));
	  }
	}
	
	if (zones.size() != 1) {
	  dof_enrichments[i] = &(*(enrichments.insert(zones).first));
	  enriched_dofs.add(i);
	  for (mesh::ind_cv_ct::const_iterator it = ct.begin();
	       it != ct.end(); ++it) enriched_elements.add(*it);
	}
      }
    }

    // cout << "Enriched convexes : " << enriched_elements << endl;
    // cout << "Enriched dofs : " << enriched_dofs << endl;

    for (dal::bv_visitor i(mf.convex_index()); !i.finished(); ++i) {
      if (enriched_elements[i]) build_method_of_convex(i); else
	set_finite_element(i, mf.fem_of_element(i));
    }
    is_adapted = true; touch();
  }


}  /* end of namespace getfem.                                            */

