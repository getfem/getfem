/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
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

#ifndef GETFEMINT_LEVELSET_H__
#define GETFEMINT_LEVELSET_H__

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_level_set.h>
#include <getfem/getfem_arch_config.h>

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
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H
    void values_from_func(unsigned idx, const std::string &s);
#endif
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
