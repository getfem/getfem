// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_mesh_adaptable_im.h : adaptable Integration methods
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


#ifndef GETFEM_MESH_IM_H__
#define GETFEM_MESH_IM_H__

#include <getfem_mesh_im.h>
#include <getfem_level_set.h>
#include <set>

namespace getfem {

  /// Describe an adaptable integration method linked to a mesh.
  class mesh_adaptable_im : public getfem_mesh_im {
  protected :

    typedef level_set *plevel_set;

    std::set<plevel_set> level_sets; // set of level set
    pintegration_method regular_simplex_pim;
    pintegration_method singular_simplex_pim;

    // + set of stored adapted integration methods
    
  public :

    add_level_set(level_set &ls) { level_sets.insert(&ls); }
    sup_level_set(level_set &ls) { level_sets.erase(&ls); }
    
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    void receipt(const MESH_ADD_CONVEX &m) {getfem_mesh_receiver::receipt(m); }
    void receipt(const MESH_SUP_CONVEX &m);
    void receipt(const MESH_SWAP_CONVEX &m);
    
    size_type memsize() const {
      return getfem_mesh_im::memsize() + ... ;
    }
    
    mesh_adaptable_im(getfem_mesh &me);
    
  private:
    mesh_adaptable_im(const mesh_adaptable_im &);
    mesh_adaptable_im & operator=(const meshadaptable__im &);
  };
  
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_MESH_IM_H__  */
