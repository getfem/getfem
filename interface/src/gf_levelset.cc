/*===========================================================================

 Copyright (C) 2005-2020 Julien Pommier.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
#include <getfemint_workspace.h>
#include <getfemint_levelset.h>
#include <getfem/getfem_arch_config.h>

using namespace getfemint;

/*@GFDOC

   The level-set object is represented by a primary level-set and optionally
   a secondary level-set used to represent fractures (if p(x) is the primary
   level-set function and s(x) is the secondary level-set, the crack is
   defined by :math:`p(x)=0` and :math:`s(x)\leq0` : the role of the secondary is to determine
   the crack front/tip).

   note:

      All tools listed below need the package qhull installed on your
      system. This package is widely available. It computes convex hull and
      delaunay triangulations in arbitrary dimension.

@*/

void gf_levelset(
getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (check_cmd("LevelSet", "LevelSet", in, out, 2, 4, 0, 1)) {
    /*@INIT LS = ('.mesh',@tmesh m, @int d[, @str 'ws'| @str f1[, @str f2 | @str 'ws']])
      Create a @tls object on a @tmesh represented by a primary function
      (and optional secondary function, both) defined on a lagrange @tmf of
      degree `d`.

      If `ws` (with secondary) is set; this levelset is represented by a
      primary function and a secondary function. If `f1` is set; the primary
      function is defined by that expression (with the syntax of the high
      generic assembly language). If `f2` is set; this levelset
      is represented by a primary function and a secondary function defined
      by these expressions. @*/
    getfem::mesh *mm = to_mesh_object(in.pop());
    size_type degree = in.pop().to_integer(1, 20);

    bool with_secondary = false;
    std::string s1="",s2="";
    if (in.remaining() && in.front().is_string()) s1 = in.pop().to_string();
    if (cmd_strmatch(s1, "ws") || cmd_strmatch(s1, "with_secondary")) {
      with_secondary = true; s1 = "";
    } else if (in.remaining() && in.front().is_string()) {
        s2 = in.pop().to_string();
        with_secondary = true;
        if (cmd_strmatch(s1, "ws") || cmd_strmatch(s2,"with_secondary")) s2 = "";
    }

    auto pls = std::make_shared<getfem::level_set>(*mm, dim_type(degree),
						   with_secondary);
    id_type id = store_levelset_object(pls);

    if (s1.size()) values_from_func(pls.get(), 0, s1);
    if (s2.size()) values_from_func(pls.get(), 1, s2);

    workspace().set_dependence(id, mm);
    out.pop().from_object_id(id, LEVELSET_CLASS_ID);
  }
}
