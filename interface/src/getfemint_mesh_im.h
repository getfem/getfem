/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2005-2012 Julien pommier.
 
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


#ifndef GETFEMINT_MESH_IM_H__
#define GETFEMINT_MESH_IM_H__

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_mesh_im.h>

namespace getfemint
{

  class getfemint_mesh_im : public getfemint::getfem_object {
  private:
    getfem::mesh_im *mim;
    id_type linked_mesh_id_;
    getfemint_mesh_im(getfem::mesh_im *mim_, id_type idmesh);
  public:
    ~getfemint_mesh_im();
    id_type class_id() const { return MESHIM_CLASS_ID; }
    size_type memsize() const { return mim->memsize(); }
    
    void clear_before_deletion() { if (!is_static()) mim->clear(); }

    getfem::mesh_im& mesh_im() { return *mim; }
    const getfem::mesh_im& mesh_im() const { return *mim; }
    const getfem::mesh& linked_mesh() const { return mim->linked_mesh(); }
    id_type linked_mesh_id() const { return linked_mesh_id_;}

    static getfemint_mesh_im* get_from(getfem::mesh_im *mim, 
				       int flags=0);
    static getfemint_mesh_im* new_from(getfemint_mesh *mim);
  };

  inline bool object_is_mesh_im(getfem_object *o) {
    return (o->class_id() == MESHIM_CLASS_ID);
  }

  inline getfemint_mesh_im* object_to_mesh_im(getfem_object *o) {
    if (object_is_mesh_im(o)) return ((getfemint_mesh_im*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_MESH_IM_H__                                           */
