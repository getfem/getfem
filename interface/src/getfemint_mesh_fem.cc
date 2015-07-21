/*===========================================================================
 
 Copyright (C) 2007-2015 Yves Renard, Julien Pommier.
 
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
#include <getfemint_mesh_fem.h>
#include <getfemint_workspace.h>

namespace getfemint {

  getfemint_mesh_fem::getfemint_mesh_fem(getfem::mesh_fem *mf_,
                                         id_type idmesh) {
    linked_mesh_id_ = idmesh;
    mf = mf_;
    ikey = getfem_object::internal_key_type(mf);
  }

  getfemint_mesh_fem::~getfemint_mesh_fem() {
    if (!is_static()) delete mf;
    mf = 0;
  }

  getfemint_mesh_fem *
  getfemint_mesh_fem::get_from(getfem::mesh_fem *mf_,
                               int flags) {
    getfem_object *o =
      getfemint::workspace().object(getfem_object::internal_key_type(mf_));
    getfemint_mesh_fem *gmf = 0;
    if (!o) {
      getfemint_mesh *gm =
        getfemint_mesh::get_from(const_cast<getfem::mesh*>(&mf_->linked_mesh()),
                                 flags);
      gmf = new getfemint_mesh_fem(mf_, gm->get_id());
      gmf->set_flags(flags);
      getfemint::workspace().push_object(gmf);
      getfemint::workspace().set_dependance(gmf, gm);

    } else gmf = dynamic_cast<getfemint_mesh_fem*>(o);
    assert(gmf);
    return gmf;
  }

  getfemint_mesh_fem *
  getfemint_mesh_fem::new_from(getfemint_mesh *m,
                               size_type q_dim) {
    //assert(::workspace == 0);
    getfem::mesh_fem *mf_ = new getfem::mesh_fem(m->mesh());
    mf_->set_qdim(dim_type(q_dim));
    getfemint_mesh_fem *gmf = get_from(mf_);
    assert(gmf->linked_mesh_id() == m->get_id());
    return gmf;
  }
}
