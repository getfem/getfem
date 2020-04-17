/*===========================================================================

 Copyright (C) 2009-2020 Yves Renard.

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

/**\file gf_model.cc
   \brief model construction wrapper.
*/

#include <getfem/getfem_mesh_im.h>
#include <getfem/getfem_models.h>
#include <getfemint.h>
#include <getfemint_workspace.h>


using namespace getfemint;


/*@GFDOC
  @tmodel variables store the variables and the state data and the
  description of a model. This includes the global tangent matrix, the right
  hand side and the constraints. There are two kinds of models, the `real`
  and the `complex` models.
@*/

void gf_model(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) THROW_BADARG( "Wrong number of input arguments");

  if (in.front().is_string()) {
    std::string cmd = in.pop().to_string();
    if (check_cmd(cmd, "real", in, out, 0, 0, 0, 1)) {
      /*@INIT MD = ('real')
      Build a model for real unknowns.@*/
      auto md = std::make_shared<getfem::model>(false);
      out.pop().from_object_id(store_model_object(md), MODEL_CLASS_ID);
    } else if (check_cmd(cmd, "complex", in, out, 0, 0, 0, 1)) {
      /*@INIT MD = ('complex')
      Build a model for complex unknowns.@*/
      auto md = std::make_shared<getfem::model>(true);
      out.pop().from_object_id(store_model_object(md), MODEL_CLASS_ID);
    } else bad_cmd(cmd);
  } else THROW_BADARG("expected a string");

  if (in.remaining()) THROW_BADARG("too many arguments");
}
