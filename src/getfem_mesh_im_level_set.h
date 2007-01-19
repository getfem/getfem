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

/**@file getfem_mesh_im_level_set.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>
   @date February 02, 2005.
   @brief a subclass of mesh_im which is conformal to a number of level sets.
*/

#ifndef GETFEM_MESH_IM_LEVEL_SET_H__
#define GETFEM_MESH_IM_LEVEL_SET_H__

#include <getfem_mesh_im.h>
#include <getfem_mesh_level_set.h>
#include <set>

namespace getfem {

  /** 
     Describe an adaptable integration method linked to a mesh cut by
     a level set.  It is possible to choose to integrate over the
     whole mesh, or to select integration on the "inside" (the
     intersection of the negative parts of the levelsets), the
     "outside" or just the levelset boundary.
  */

  class mesh_im_level_set : public mesh_im {
  protected :
    pintegration_method regular_simplex_pim;
    pintegration_method base_singular_pim;
    mesh_level_set &mls;

    mesh_im cut_im; /* stores an im only for convexes who are crossed
		       by a levelset */

    dal::bit_vector ignored_im; /* convex list whose integration method is
				   ignored (for instance because
				   INTEGRATE_INSIDE and the convex
				   is outside etc.) */
    std::vector<pintegration_method> build_methods;

    mutable bool is_adapted;
    int integrate_where; // INTEGRATE_INSIDE or INTEGRATE_OUTSIDE

    void clear_build_methods();
    void build_method_of_convex(size_type cv);
    bool is_point_in_selected_area
    (const std::vector<mesher_level_set> &mesherls0,
     const std::vector<mesher_level_set> &mesherls1, const base_node& P);

  public :
    enum { INTEGRATE_INSIDE = 1, INTEGRATE_OUTSIDE = 2, INTEGRATE_ALL = 2+1,
           INTEGRATE_BOUNDARY = 4};
    void update_from_context(void) const;
    
    /** Apply the adequate integration methods. */
    void adapt(void);
    void clear(void); // to be modified

    /** Set the specific integration methods. see the constructor
	documentation for more details. */
    void set_simplex_im(pintegration_method reg,
			pintegration_method sing = 0) {
      regular_simplex_pim = reg;
      base_singular_pim = sing;
    }
    
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    // void receipt(const MESH_ADD_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SUP_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SWAP_CONVEX &m) // to be modified ?
    
    size_type memsize() const {
      return mesh_im::memsize(); // + ... ;
    }

    /** 
	@param me the level-set.
	
	@param integrate_where : choose between INTEGRATE_ALL,
	INTEGRATE_BOUNDARY, INTEGRATE_INSIDE and INTEGRATE_OUTSIDE.

	@param reg the integration method (for simplices) that will be
	used on the sub-simplices of convexes crossed by the levelset.

	@param sing the (optional) integration method to use on the crack tips
	(i.e. when the levelset has a secondary level set), this is
	generally an IM_QUASI_POLAR method as it provides a good
	integration of singular XFEM functions.
    */
    mesh_im_level_set(mesh_level_set &me, 
		      int integrate_where_ = INTEGRATE_ALL,
		      pintegration_method reg = 0,
		      pintegration_method sing = 0);

    virtual pintegration_method int_method_of_element(size_type cv) 
      const;
    ~mesh_im_level_set() { clear_build_methods(); }
  private:
    mesh_im_level_set(const mesh_im_level_set &);
    mesh_im_level_set & operator=(const mesh_im_level_set &);
  };

  
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_IM_LEVEL_SET_H__  */
