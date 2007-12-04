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
#include <getfemint_mesh_levelset.h>
#include <getfemint_mesh_levelset.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_mesh_levelset::getfemint_mesh_levelset() {
    mls = 0;
  }

  getfemint_mesh_levelset::~getfemint_mesh_levelset() {
    if (!is_static()) delete mls; 
    mls = 0;
  }

  getfemint_mesh_levelset *
  getfemint_mesh_levelset::get_from(getfem::mesh_level_set *mls,
				    int flags) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(mls));
    getfemint_mesh_levelset *gmls = 0;
    if (!o) {
      getfemint_mesh *gm = 
	getfemint_mesh::get_from(const_cast<getfem::mesh*>(&mls->linked_mesh()),
				 flags);
      gmls = new getfemint_mesh_levelset(); 
      gmls->mls = mls;
      gmls->ikey = getfem_object::internal_key_type(mls);
      gmls->set_flags(flags); 
      getfemint::workspace().push_object(gmls);
      getfemint::workspace().set_dependance(gmls, gm);
    } else gmls = dynamic_cast<getfemint_mesh_levelset*>(o);
    assert(gmls);
    return gmls;
  }
}
