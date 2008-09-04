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
#include <getfemint_levelset.h>
#include <getfemint_workspace.h>

using namespace getfemint;

/*@TEXT LEVELSET:INIT('LEVELSET_init')
  General constructor for LevelSet objects.<Par>
  @*/
void gf_levelset(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) THROW_BADARG("Wrong number of input arguments");
  if (!out.narg_in_range(1,1)) 
    THROW_BADARG("Wrong number of output arguments");
  getfemint_mesh *mm = in.pop().to_getfemint_mesh();
  size_type degree = in.pop().to_integer(1, 20);
  bool with_secondary = false;
  std::string s1, s2;
  if (in.remaining() && in.front().is_string()) {
    s1 = in.pop().to_string();
    if (cmd_strmatch(s1, "with_secondary")) {
      with_secondary = true; s1 = "";
      if (in.remaining() && in.front().is_string())
	s1 = in.pop().to_string();
    } 
    if (in.remaining() && in.front().is_string()) {
      s2 = in.pop().to_string();
      with_secondary = true;
    }
  }

  if (in.remaining()) 
    THROW_BADARG("too many arguments");
    

  getfem::level_set *ls = 
    new getfem::level_set(mm->mesh(),dim_type(degree),with_secondary);
  getfemint_levelset *gls = getfemint_levelset::get_from(ls);
  if (s1.size()) gls->values_from_poly(0, s1);
  if (s2.size()) gls->values_from_poly(1, s2);

  out.pop().from_object_id(gls->get_id(), LEVELSET_CLASS_ID);
}

