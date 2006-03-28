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

  /** Describe an adaptable integration method linked to a mesh cutted by
   * a level set.
   */
  class mesh_im_level_set : public mesh_im {
  protected :
    pintegration_method regular_simplex_pim;
    pintegration_method base_singular_pim;
    mesh_level_set &mls;

    mesh_im cut_im; /* stores an im only for convexes who are crossed
		       by a levelset */
    std::vector<pintegration_method> build_methods;

    mutable bool is_adapted;
    int integrate_where; // INTEGRATE_INSIDE or INTEGRATE_OUTSIDE

    void clear_build_methods();
    void build_method_of_convex(size_type cv);

  public :
    enum { INTEGRATE_INSIDE = 1, INTEGRATE_OUTSIDE = 2, INTEGRATE_ALL = 2+1,
           INTEGRATE_BOUNDARY = 4};
    void update_from_context(void) const { is_adapted = false; }
    void adapt(void);
    void clear(void); // to be modified

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
