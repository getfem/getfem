/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2014-2015 Konstantinos Poulios.
 
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

/**\file getfemint_mesh_im_data.h
   \brief getfem::im_data wrapper class
*/

#ifndef GETFEMINT_MESH_IM_DATA_H__
#define GETFEMINT_MESH_IM_DATA_H__

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfemint_mesh_im.h>
#include <getfem/getfem_im_data.h>

namespace getfemint
{

  class getfemint_mesh_im_data : public getfemint::getfem_object {
  private:
    getfem::im_data *mimd;
    id_type linked_mesh_im_id_;
    getfemint_mesh_im_data(getfem::im_data *mimd_, id_type idmeshim);
  public:
    ~getfemint_mesh_im_data();
    id_type class_id() const { return MESHIMDATA_CLASS_ID; }
//    size_type memsize() const { return mimd->memsize(); }

//    void clear_before_deletion() {
//      if (!is_static()) mimd->clear();
//    }

    getfem::im_data& mesh_im_data() { return *mimd; }
    const getfem::im_data& mesh_im_data() const { return *mimd; }
    const getfem::mesh_im& linked_mesh_im() const { return mimd->linked_mesh_im(); }
    id_type linked_mesh_im_id() const { return linked_mesh_im_id_;}
    id_type linked_mesh_id() const;

    static getfemint_mesh_im_data* get_from(getfem::im_data *mimd_, int flags=0);
    static getfemint_mesh_im_data* new_from(getfemint_mesh_im *mim,
                                            size_type region,
                                            bgeot::multi_index tensor_size);
  };

  inline bool object_is_mesh_im_data(getfem_object *o) {
    return (o->class_id() == MESHIMDATA_CLASS_ID);
  }

  inline getfemint_mesh_im_data* object_to_mesh_im_data(getfem_object *o) {
    if (object_is_mesh_im_data(o)) return ((getfemint_mesh_im_data*)o);
    else THROW_INTERNAL_ERROR;
  }
}  /* end of namespace getfemint.                                          */

#endif /* GETFEMINT_MESH_IM_DATA_H__                                       */
