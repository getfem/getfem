// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2005-2006 Julien Pommier, Yves Renard.
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

/**\file getfemint_mesh.h
   \brief getfem::mesh interface.
*/

#ifndef GETFEMINT_MESH_H__
#define GETFEMINT_MESH_H__

#include <getfemint_std.h>
#include <getfemint_object.h>
#include <getfem/getfem_mesh.h>

namespace getfemint
{
  class getfemint_mesh : public getfem_object {
  private:
    getfem::mesh *m;
    getfemint_mesh(getfem::mesh *m);
  public:
    ~getfemint_mesh();
    id_type class_id() const { return MESH_CLASS_ID; }
    size_type memsize() const { return m->memsize(); }
    getfem::mesh& mesh() { return *m; }
    const getfem::mesh& mesh() const { return *m; }
    static getfemint_mesh* get_from(getfem::mesh *m, 
				    int flags = 0);
  };

  inline bool object_is_mesh(getfem_object *o) {
    return (o->class_id() == MESH_CLASS_ID);
  }

  inline getfemint_mesh* object_to_mesh(getfem_object *o) {
    if (object_is_mesh(o)) return ((getfemint_mesh*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_MESH_H__                                               */
