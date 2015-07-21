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
 
===========================================================================*/
#include <getfemint_mesh_im_data.h>
#include <getfemint_workspace.h>

namespace getfemint {

  getfemint_mesh_im_data::getfemint_mesh_im_data(getfem::im_data *mimd_,
                                                 id_type idmeshim) {
    linked_mesh_im_id_ = idmeshim;
    mimd = mimd_;
    ikey = getfem_object::internal_key_type(mimd);
  }

  getfemint_mesh_im_data::~getfemint_mesh_im_data() {
    if (!is_static()) delete mimd;
    mimd = 0;
  }

  id_type getfemint_mesh_im_data::linked_mesh_id() const {
      getfem_object *o = getfemint::workspace().object(linked_mesh_im_id(),
                                                       "gfMeshIm");
      return object_to_mesh_im(o)->linked_mesh_id();;
  }

  getfemint_mesh_im_data *
  getfemint_mesh_im_data::get_from(getfem::im_data *mimd_, int flags) {
    getfem_object *o =
      getfemint::workspace().object(getfem_object::internal_key_type(mimd_));
    getfemint_mesh_im_data *gmimd = 0;
    if (!o) {
      getfemint_mesh_im *gmim =
        getfemint_mesh_im::get_from(const_cast<getfem::mesh_im*>(&mimd_->linked_mesh_im()),
                                    flags);
      gmimd = new getfemint_mesh_im_data(mimd_, gmim->get_id());
      gmimd->set_flags(flags);
      getfemint::workspace().push_object(gmimd);
      getfemint::workspace().set_dependance(gmimd, gmim);

    } else gmimd = dynamic_cast<getfemint_mesh_im_data*>(o);
    assert(gmimd);
    return gmimd;
  }

  getfemint_mesh_im_data *
  getfemint_mesh_im_data::new_from(getfemint_mesh_im *mim,
                                   size_type region,
                                   bgeot::multi_index tensor_size) {
    //assert(::workspace == 0);
    getfem::im_data *mimd_ = new getfem::im_data(mim->mesh_im());
    mimd_->set_region(region);
    mimd_->set_tensor_size(tensor_size);
    getfemint_mesh_im_data *gmimd = get_from(mimd_);
    assert(gmimd->linked_mesh_im_id() == mim->get_id());
    return gmimd;
  }
}
