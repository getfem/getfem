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

#ifndef GETFEMINT_LEVELSET_H__
#define GETFEMINT_LEVELSET_H__

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_level_set.h>

namespace getfemint
{
  class getfemint_levelset : public getfemint::getfem_object {
  private:
    getfem::level_set *ls;
    getfemint_levelset() {}
  public:
    /*getfemint_levelset(getfemint_mesh *m, dim_type degree, 
		       bool with_secondary) { 
      assert(workspace == 0);
      linked_mesh_id_ = m->get_id();
      ls = new getfem::level_set(m->mesh(), degree, with_secondary);
      }*/

    ~getfemint_levelset() {
      if (!is_static()) delete ls; ls = 0;
    }
    id_type class_id() const { return LEVELSET_CLASS_ID; }
    size_type memsize() const { return ls->memsize(); }
    
    getfem::level_set& levelset() { return *ls; }
    const getfem::mesh_fem& mesh_fem() const { return ls->get_mesh_fem(); }
    //id_type linked_mesh_id() const { return linked_mesh_id_;}

    static getfemint_levelset* get_from(getfem::level_set *ls, int flags=0);
    void values_from_poly(unsigned idx, const std::string &s);
  };

  inline bool object_is_levelset(getfem_object *o) {
    return (o->class_id() == LEVELSET_CLASS_ID);
  }

  inline getfemint_levelset* object_to_levelset(getfem_object *o) {
    if (object_is_levelset(o)) return ((getfemint_levelset*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_LEVELSET_H__                                           */
