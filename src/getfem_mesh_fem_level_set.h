// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_fem_level_set.h : definition of a finite element
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

/**@file getfem_mesh_fem_level_set.h
   @brief a subclass of mesh_fem which is conformal to a number of level sets.
*/

#ifndef GETFEM_MESH_FEM_LEVEL_SET_H__
#define GETFEM_MESH_FEM_LEVEL_SET_H__

#include <getfem_mesh_level_set.h>
#include <getfem_mesh_fem.h>
#include <getfem_fem_level_set.h>

namespace getfem {

  class mesh_fem_level_set : public mesh_fem, public boost::noncopyable {
  protected :
    const mesh_level_set &mls;
    const mesh_fem &mf;
    mutable std::vector<pfem> build_methods;
    mutable bool is_adapted;
    mutable dal::bit_vector enriched_dofs, enriched_elements;
    mutable std::set<std::set<const mesh_level_set::zone *> > enrichments;
    mutable std::vector< const std::set<const mesh_level_set::zone *> *> dof_enrichments;
    size_type xfem_index;
    void clear_build_methods();
    void build_method_of_convex(size_type cv);

  public :
    void update_from_context(void) const { is_adapted = false; }
    void adapt(void);
    void clear(void); // to be modified
    
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    // void receipt(const MESH_ADD_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SUP_CONVEX &m)  // to be modified ?
    // void receipt(const MESH_SWAP_CONVEX &m) // to be modified ?
    
    size_type memsize() const {
      return mesh_fem::memsize(); // + ... ;
    }
    
    mesh_fem_level_set(const mesh_level_set &me, const mesh_fem &mef);

    ~mesh_fem_level_set() { clear_build_methods(); }
  };



  
}  /* end of namespace getfem.                                            */

#endif
  
