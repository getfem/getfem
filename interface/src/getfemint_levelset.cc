// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2009 Julien Pommier.
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
#include <getfemint_levelset.h>
#include <getfemint_levelset.h>
#include <getfemint_workspace.h>

namespace getfemint {
  getfemint_levelset* getfemint_levelset::get_from(getfem::level_set *ls,
						   int flags) {
    getfem_object *o = 
      getfemint::workspace().object(getfem_object::internal_key_type(ls));
    getfemint_levelset *gls = 0;
    if (!o) {
      const getfem::mesh &m = ls->get_mesh_fem().linked_mesh();
      getfemint_mesh *gm = 
	getfemint_mesh::get_from(const_cast<getfem::mesh*>(&m),
				 flags);
      gls = new getfemint_levelset(); 
      gls->ls = ls;
      gls->ikey = getfem_object::internal_key_type(ls);
      gls->set_flags(flags); 
      getfemint::workspace().push_object(gls);
      getfemint::workspace().set_dependance(gls, gm);
    } else gls = dynamic_cast<getfemint_levelset*>(o);
    assert(gls);
    return gls;

  }

  void getfemint_levelset::values_from_poly(unsigned idx, 
					    const std::string &s) {
    const getfem::mesh_fem &mf = levelset().get_mesh_fem();
    assert(!mf.is_reduced());
    bgeot::base_poly p = 
      bgeot::read_base_poly(mf.linked_mesh().dim(), s);
    ls->values(idx).resize(mf.nb_dof());
    for (unsigned i=0; i < mf.nb_dof(); ++i) {
      const getfem::base_node x = mf.point_of_basic_dof(i);
      ls->values(idx)[i] = p.eval(x.begin());
    }
  }
}
