// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================
#include <getfemint_mesh_fem.h>
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
  getfemint_mesh_fem::get_from(getfem::mesh_fem *mf, 
			       int flags) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(mf));
    getfemint_mesh_fem *gmf = 0;
    if (!o) {
      getfemint_mesh *gm = 
	getfemint_mesh::get_from(const_cast<getfem::mesh*>(&mf->linked_mesh()),
				 flags);
      gmf = new getfemint_mesh_fem(mf, gm->get_id());
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
    getfem::mesh_fem *mf = new getfem::mesh_fem(m->mesh());
    mf->set_qdim(q_dim);
    getfemint_mesh_fem *gmf = get_from(mf);
    assert(gmf->linked_mesh_id() == m->get_id());
    return gmf;
  }
}
