// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Elements Methods (getfem)
// File    : getfem_mesh_adaptable_im.cc : adaptable integration methods
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

#include <getfem_mesh_adaptable_im.h>


namespace getfem {
  
  void mesh_adaptable_im::receipt(const MESH_CLEAR &) { clear(); }
  void mesh_adaptable_im::receipt(const MESH_DELETE &) { clear(); }
  void mesh_adaptable_im::clear(void) { mesh_im::clear(); level_sets.clear(); }

  mesh_adaptable_im::mesh_adaptable_im(getfem_mesh &me) : mesh_im(me) { }


  
}  /* end of namespace getfem.                                             */



