// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_im_level_set.h : adaptable Integration methods
//           on convex meshes.
//           
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


#ifndef GETFEM_MESH_IM_LEVEL_SET_H__
#define GETFEM_MESH_IM_LEVEL_SET_H__

#include <getfem_mesh_im.h>
#include <getfem_level_set.h>
#include <set>

namespace getfem {

  /// Describe an adaptable integration method linked to a mesh.
  class mesh_im_level_set : public mesh_im {
  protected :

    typedef level_set *plevel_set;

    std::set<plevel_set> level_sets; // set of level set
    pintegration_method regular_simplex_pim;
    pintegration_method singular_simplex_pim; // en 3D ?

    mesh_im cut_im; /* stores an im only for convexes who are crossed
		       by a levelset */
    std::vector<pintegration_method> build_methods;

  public :

    void adapt(void);

    void clear(void); // to be modified

    void add_level_set(level_set &ls) { level_sets.insert(&ls); }
    void sup_level_set(level_set &ls) { level_sets.erase(&ls); }
    void set_simplex_im(pintegration_method reg,
			pintegration_method sing = 0) {
      regular_simplex_pim = reg;
      singular_simplex_pim = (sing == 0) ? reg : sing;
    }
    
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    // void receipt(const MESH_ADD_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SUP_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SWAP_CONVEX &m) // to be modified ?
    
    size_type memsize() const {
      return mesh_im::memsize(); // + ... ;
    }
    
    mesh_im_level_set(getfem_mesh &me, pintegration_method reg = 0,
		      pintegration_method sing = 0);

    virtual pintegration_method int_method_of_element(size_type cv) 
      const;
  private:
    mesh_im_level_set(const mesh_im_level_set &);
    mesh_im_level_set & operator=(const mesh_im_level_set &);
    void cut_element(size_type cv, const dal::bit_vector &primary,
		const dal::bit_vector &secondary);    
    bool is_crossed_by(size_type c, plevel_set ls, unsigned lsnum);
    void find_crossing_level_set(size_type cv, dal::bit_vector &prim, dal::bit_vector &sec);
    void run_delaunay(std::vector<base_node> &fixed_points,
		      gmm::dense_matrix<size_type> &simplexes,
		      std::vector<dal::bit_vector> &fixed_points_constraints);

  };

  void getfem_mesh_im_level_set_noisy(void); // for debugging
  
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_IM_LEVEL_SET_H__  */
