// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2005-2008 Julien Pommier.
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
#include <getfemint_mesh_levelset.h>
#include <getfemint_workspace.h>

using namespace getfemint;


void gf_mesh_levelset(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) THROW_BADARG("Wrong number of input arguments");
  if (!out.narg_in_range(1,1)) 
    THROW_BADARG("Wrong number of output arguments");
  getfemint_mesh *mm = in.pop().to_getfemint_mesh();
  getfem::mesh_level_set *mls = new getfem::mesh_level_set(mm->mesh());
  getfemint_mesh_levelset *gmls = getfemint_mesh_levelset::get_from(mls);
  out.pop().from_object_id(gmls->get_id(), MESH_LEVELSET_CLASS_ID);
}

