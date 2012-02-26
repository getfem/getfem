/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
 
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
 
===========================================================================*/
#include <getfemint_mesh.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_mesh::getfemint_mesh(getfem::mesh *m_) {
    assert(workspace == 0);
    m = m_;
    ikey = getfem_object::internal_key_type(m);
  }

  getfemint_mesh::~getfemint_mesh() {
    if (!is_static()) {
      m->clear();
      delete m;
    }
  }

  getfemint_mesh *getfemint_mesh::get_from(getfem::mesh *m, int f) {
    getfem_object *o =
      getfemint::workspace().object(getfem_object::internal_key_type(m));
    getfemint_mesh *gm = 0;
    if (!o) {
      gm = new getfemint_mesh(m);
      gm->set_flags(f);
      getfemint::workspace().push_object(gm);
    } else gm = dynamic_cast<getfemint_mesh*>(o);
    assert(gm);
    return gm;
  }
}
