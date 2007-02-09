// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2006 Yves Renard, Julien Pommier. Julien Pommier.
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

/**\file getfemint_mesh_fem.h
   \brief getfem::mesh_fem wrapper class
*/

#ifndef GETFEMINT_MESH_FEM_H__
#define GETFEMINT_MESH_FEM_H__

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_mesh_fem.h>

namespace getfemint
{

  class getfemint_mesh_fem : public getfemint::getfem_object {
  private:
    getfem::mesh_fem *mf;
    id_type linked_mesh_id_;
    getfemint_mesh_fem(getfem::mesh_fem *mf_, id_type idmesh);
  public:
    ~getfemint_mesh_fem();
    id_type class_id() const { return MESHFEM_CLASS_ID; }
    size_type memsize() const { return mf->memsize(); }
    
    void clear_before_deletion() { 
      if (!is_static()) mf->clear(); 
    }

    getfem::mesh_fem& mesh_fem() { return *mf; }
    const getfem::mesh_fem& mesh_fem() const { return *mf; }
    const getfem::mesh& linked_mesh() const { return mf->linked_mesh(); }
    id_type linked_mesh_id() const { return linked_mesh_id_;}

    static getfemint_mesh_fem* get_from(getfem::mesh_fem *mf, 
					int flags=0);
    static getfemint_mesh_fem* new_from(getfemint_mesh *mm,
					size_type qdim);
  };

  inline bool object_is_mesh_fem(getfem_object *o) {
    return (o->class_id() == MESHFEM_CLASS_ID);
  }

  inline getfemint_mesh_fem* object_to_mesh_fem(getfem_object *o) {
    if (object_is_mesh_fem(o)) return ((getfemint_mesh_fem*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_MESH_FEM_H__                                           */
