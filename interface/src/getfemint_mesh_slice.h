/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2001-2015 Yves Renard, Julien Pommier.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/


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
