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
The level-set object is represented by a primary level-set and optionally<par>
a secondary level-set used to represent fractures (if p(x) is the primary<par>
level-set function and s(x) is the secondary level-set, the crack is defined<par>
by p(x)=0 and s(x)<=0: the role of the secondary is stop the crack).<Par>

**IMPORTANT:**<Par>

All tools listed below needs the package qhull installed on your<par>
system. This package is widely available. It computes convex hull<par>
and delaunay triangulations in arbitrary dimension. Everything<par>
here is considered *work in progress*, it is still subject to<par>
major changes if needed.@*/

void gf_levelset(
getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  getfemint_levelset *gls = NULL;
  if (check_cmd("LevelSet", "LevelSet", in, out, 2, 4, 0, 1)) {
    /*@INIT LEVELSET:INIT('.mesh', @tmesh m, @int d[, 'ws'| @str poly1[, @str poly2| 'ws']])
    Create a @tls object on a @tmesh represented by a primary function
    (and optional secondary function, both) defined on a lagrange @tmf
    of degree `d`.<Par>

    If `ws` (with secondary) is set; this levelset is represented by a
    primary function and a secondary function. If `poly1` is set; the
    primary function is defined by that polynomial expression. If
    `poly2` is set; this levelset is represented by a primary function
    and a secondary function defined by these polynomials expressions.@*/
    getfemint_mesh *mm = in.pop().to_getfemint_mesh();
    size_type degree = in.pop().to_integer(1, 20);

    bool with_secondary = false;
    std::string s1="",s2="";
    if (in.remaining() && in.front().is_string()) s1 = in.pop().to_string();
    if (cmd_strmatch(s1, "ws") || cmd_strmatch(s1, "with_secondary")) {
      with_secondary = true; s1 = "";
    }else if (in.remaining() && in.front().is_string()) {
        s2 = in.pop().to_string();
        with_secondary = true;
        if (cmd_strmatch(s1, "ws") || cmd_strmatch(s2,"with_secondary")) s2 = "";
    }

    getfem::level_set *ls =
      new getfem::level_set(mm->mesh(),dim_type(degree),with_secondary);
    gls = getfemint_levelset::get_from(ls);

    if (s1.size()) gls->values_from_poly(0, s1);
    if (s2.size()) gls->values_from_poly(1, s2);
    workspace().set_dependance(gls, mm);
  }
  out.pop().from_object_id(gls->get_id(), LEVELSET_CLASS_ID);
}
