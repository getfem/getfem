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

    static getfemint_mesh_fem* get_from(getfem::mesh_fem *mf_,
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
