/*===========================================================================
 
 Copyright (C) 2006-2015 Julien Pommier.
 
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

#include <getfemint.h>
#include <getfemint_levelset.h>
#include <getfemint_mesh_fem.h>

using namespace getfemint;

/*@GFDOC
  General function for modification of LEVELSET objects.
@*/

void gf_levelset_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_levelset *gls = in.pop().to_getfemint_levelset(true);
  getfem::level_set &ls = gls->levelset();
  std::string cmd = in.pop().to_string();
  if (check_cmd(cmd, "values", in, out, 1, 2, 0, 0)) {
    /*@SET ('values', {@mat v1|@str func_1}[, @mat v2|@str func_2])
    Set values of the vector of dof for the level-set functions.

    Set the primary function with the vector of dof `v1` (or the expression
    `func_1`) and the secondary function (if any) with  the vector of dof
    `v2` (or the expression `func_2`)@*/
    std::string s1, s2;
    darray v1, v2;
    if (in.front().is_string()) {
      s1 = in.pop().to_string();
    } else {
      v1 = in.pop().to_darray(int(ls.get_mesh_fem().nb_dof()));
    }
    if (in.remaining()) {
      if (!ls.has_secondary())
	THROW_BADARG("The levelset has not secondary term");
      if (in.front().is_string()) {
	s2 = in.pop().to_string();
      } else {
	v2 = in.pop().to_darray(int(ls.get_mesh_fem().nb_dof()));
      }
    }
    ls.values(0).resize(ls.get_mesh_fem().nb_dof());
    if (s1.size()) {
      gls->values_from_func(0, s1);
    } else {
      ls.values(0).assign(v1.begin(), v1.end());
    }
    if (ls.has_secondary()) {
      ls.values(1).resize(ls.get_mesh_fem().nb_dof());
      if (s2.size()) {
        gls->values_from_func(1, s2);
      } else {
	ls.values(1).assign(v2.begin(), v2.end());
      }
    }
  } else if (check_cmd(cmd, "simplify", in, out, 0, 1, 0, 0)) {
    /*@SET ('simplify'[, @scalar eps=0.01])
    Simplify dof of level-set optionally with the parameter `eps`.@*/
    if (in.remaining()==0) ls.simplify();
    else{
      ls.simplify(in.pop().to_scalar());
    }
  } else bad_cmd(cmd);
}
