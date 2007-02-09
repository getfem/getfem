// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2001-2006 Yves Renard, Julien Pommier.
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


#ifndef GETFEMINT_MESH_SLICE_H__
#define GETFEMINT_MESH_SLICE_H__

#include <getfemint_std.h>
#include <getfemint_object.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_mesh_slice.h>

namespace getfemint
{
  class getfemint_mesh_slice : public getfem_object {
  private:
    getfem::stored_mesh_slice *sl;
    id_type linked_mesh_id_;
  public:
    getfemint_mesh_slice(getfemint_mesh& mi_m, getfem::stored_mesh_slice *sl_) {
      assert(workspace == 0);
      linked_mesh_id_ = mi_m.get_id();
      sl = sl_;
    }
    ~getfemint_mesh_slice() {
      delete sl;
    }
    id_type class_id() const { return SLICE_CLASS_ID; }
    size_type memsize() const { return sl->memsize(); }
    getfem::stored_mesh_slice& mesh_slice() { return *sl; }
    const getfem::stored_mesh_slice& mesh_slice() const { return *sl; }
    id_type linked_mesh_id() const { return linked_mesh_id_;}
  };

  inline bool object_is_mesh_slice(getfem_object *o) {
    return (o->class_id() == SLICE_CLASS_ID);
  }

  inline getfemint_mesh_slice* object_to_mesh_slice(getfem_object *o) {
    if (object_is_mesh_slice(o)) return ((getfemint_mesh_slice*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_MESH_SLICE_H__                                         */
