// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mesh_adaptable_im.cc : adaptable integration methods
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


namespace getfem {
  
  void mesh_im_level_set::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_im_level_set::receipt(const MESH_DELETE &) { clear(); }
  void mesh_im_level_set::clear(void) { mesh_im::clear(); level_sets.clear(); }

  mesh_im_level_set::mesh_im_level_set(getfem_mesh &me,
				       pintegration_method reg,
				       pintegration_method sing):mesh_im(me) {
    regular_simplex_pim = reg;
    singular_simplex_pim = (sing == 0) ? reg : sing;
  }

  pintegration_method 
  mesh_im_level_set::int_method_of_element(size_type cv) const {
    if (cut_im.convex_index().is_in(cv)) 
      return cut_im.int_method_of_element(cv); 
    else return mesh_im::int_method_of_element(cv);
  }

  void mesh_adaptable_im::cut_element(size_type cv) {
    
    std::auto_ptr<mesher_signed_distance> ref_element;
    std::vector<mesher_level_set> mesher_level_sets;
    std::vector<const mesher_signed_distance *> signed_dits;
    bgeot::pgeometric_trans pgt = linked_mesh().trans_of_convex(cv);
    bool found = false;
    size_type n = pgt->structure()->dim();
    size_type nbp = pgt->basic_structure()->nb_points();
    
    /* Identifying simplexes.                                          */

    if (nbp == n+1 && pgt->basic_structure() == bgeot::simplex_structure(n)) {
	ref_element.reset(new mesher_simplex_ref(n)); found = true;
    }
    
    /* Identifying parallelepiped.                                     */

    if (!found && nbp == (size_type(1) << n) &&
	pgt->basic_structure() == bgeot::parallelepiped_structure(n)) {
      base_node rmin(n), rmax(n);
      std::fill(rmax.begin(), rmax.end(), scalar_type(1));
      ref_element.reset(new mesher_rectangle(rmin, rmax)); found = true;
    }

    /* Identifying prisms.                                             */
 
    if (!found && nbp == 2 * n &&
	pgt->basic_structure() == bgeot::prism_structure(n))
      { ref_element.reset(new mesher_prism_ref(n)); found = true; }
    
    if (!found) 
      DAL_THROW(to_be_done_error,
		"This element is not taken into account. Contact us");
    
    signed_dits.push_back(ref_element);

    for (std::set<plevel_set>::const_iterator it = level_sets.begin();
	 it != level_sets.end(); ++it) {
      pfem pf = it->get_mesh_fem().fem_of_element(cv);
      std::vector<scalar_type> coeff(pf->nb_dof(0));
      for (size_type i = 0; i < coeff.size(); ++i)
	coeff[i] = /* + - */ it->primary()[it->get_mesh_fem().ind_dof_of_element(cv, i)];
      mesher_level_sets.push_back(mesher_level_set(pf, coeff));
      signed_dits.push_back(&mesher_level_sets[mesher_level_sets.size()-1]);
    }

    mesher_intersection final_dist(signed_dits);
    getfem_mesh mesh;
    build_mesh(mesh, final_dist, 1.0); // ...
    
      
    

  }


  void mesh_im_level_set::adapt(void) {

    // compute the elements touched by each level set
    // for each element touched, compute the sub mesh
    //   then compute the adapted integration method
    
  }

#ifndef GETFEM_HAVE_QHULL_QHULL_H

    DAL_THROW(failure_error, "Qhull header files not installed. "
	      "This part of getfem++ needs Qhull."
	      "Install qhull library and reinstall Getfem++ library.");
    
#else
    
    


#endif

  }

  
}  /* end of namespace getfem.                                             */



