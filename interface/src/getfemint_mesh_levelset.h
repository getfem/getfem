/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

#ifndef GETFEMINT_MESH_LEVELSET_H__
#define GETFEMINT_MESH_LEVELSET_H__

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_mesh_level_set.h>

namespace getfemint
{
  class getfemint_mesh_levelset : public getfemint::getfem_object {
  private:
    getfem::mesh_level_set *mls;
    getfemint_mesh_levelset();
  public:
    ~getfemint_mesh_levelset();
    id_type class_id() const { return MESH_LEVELSET_CLASS_ID; }
    size_type memsize() const { return mls->memsize(); }
    
    getfem::mesh_level_set& mesh_levelset() { return *mls; }
    static getfemint_mesh_levelset *get_from(getfem::mesh_level_set *mls,
					     int flag=0);
  };

  inline bool object_is_mesh_levelset(getfem_object *o) {
    return (o->class_id() == MESH_LEVELSET_CLASS_ID);
  }

  inline getfemint_mesh_levelset* object_to_mesh_levelset(getfem_object *o) {
    if (object_is_mesh_levelset(o)) 
      return ((getfemint_mesh_levelset*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_MESH_LEVELSET_H__                                           */
