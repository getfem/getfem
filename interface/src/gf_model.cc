/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2009-2012 Yves Renard.
 
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

/**\file gf_model.cc
   \brief model construction wrapper.
*/

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_models.h>
#include <getfemint_mesh_im.h>
#include <getfemint_mesh_fem.h>


using namespace getfemint;


/*@GFDOC
  @tmodel variables store the variables and the state data and the
  description of a model. This includes the global tangent matrix, the right
  hand side and the constraints. There are two kinds of models, the `real`
  and the `complex` models.

  @tmodel object is the evolution for getfem++ 4.0 of the @tmdstate object.
@*/

void gf_model(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) THROW_BADARG( "Wrong number of input arguments");

  getfemint_model * md = new getfemint_model();
  out.pop().from_object_id(workspace().push_object(md), MODEL_CLASS_ID);

  if (in.front().is_string()) {
    std::string cmd = in.pop().to_string();
    if (check_cmd(cmd, "real", in, out, 0, 0, 0, 1)) {
      /*@INIT MD = ('real')
      Build a model for real unknowns.@*/
      md->set(new getfem::model(false));
    } else if (check_cmd(cmd, "complex", in, out, 0, 0, 0, 1)) {
      /*@INIT MD = ('complex')
      Build a model for complex unknowns.@*/
      md->set(new getfem::model(true));
    } else bad_cmd(cmd);
  } else THROW_BADARG("expected a string");

  if (in.remaining()) THROW_BADARG("too many arguments");
}
