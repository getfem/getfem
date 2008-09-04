// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Julien Pommier.
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

#include <getfemint.h>
#include <getfemint_levelset.h>
#include <getfemint_mesh_fem.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION I = gf_levelset_get(LS, ...)
    General function for querying information about LEVELSET objects.
    
MLABCOM*/

void gf_levelset_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_levelset *gls = in.pop().to_getfemint_levelset();
  getfem::level_set &ls = gls->levelset();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "values", in, out, 0, 1, 0, 1)) {
    size_type il = 0; 
    if (in.remaining()) il = in.pop().to_integer(0, 1);
    if (il != 0 && !ls.has_secondary())
      THROW_BADARG("The levelset has not secondary term");
    out.pop().from_dcvector(ls.values(unsigned(il)));
  } else if (check_cmd(cmd, "degree", in, out, 0, 0, 0, 1)) {
    out.pop().from_integer(ls.degree());
  } else if (check_cmd(cmd, "mf", in, out, 0, 0, 0, 1)) {
    
    getfem::mesh_fem *pmf = const_cast<getfem::mesh_fem*>(&ls.get_mesh_fem());
    getfemint_mesh_fem *gmf = 
      getfemint_mesh_fem::get_from(pmf, STATIC_OBJ | CONST_OBJ);
    /*cerr << "nb_dof = " << ls.get_mesh_fem().nb_dof() << "\n";
      cerr << "mf = " << &ls.get_mesh_fem() << " == " << &gmf->mesh_fem() << "\n";*/
    out.pop().from_object_id(gmf->get_id(), gmf->class_id());
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    out.pop().from_integer(int(ls.memsize()));
  } else bad_cmd(cmd);
}
