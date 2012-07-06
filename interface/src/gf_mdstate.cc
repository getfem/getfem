/*===========================================================================
 
 Copyright (C) 2005-2012 Julien Pommier.
 
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

/**\file gf_mdbrick.cc
   \brief mdbrick construction wrapper.
*/

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_mdbrick.h>
#include <getfemint_mdstate.h>
#include <getfemint_mesh_im.h>
#include <getfemint_mesh_fem.h>


using namespace getfemint;


/*@GFDOC
  A model state is an object which store the state data for a chain of model
  bricks. This includes the global tangent matrix, the right hand side and
  the constraints.

  This object is now deprecated and replaced by the @tmodel object.

  There are two sorts of model states, the `real` and the `complex` models
  states.
@*/
void gf_mdstate(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }

  getfemint_mdstate * md = new getfemint_mdstate();
  out.pop().from_object_id(workspace().push_object(md),
			   MDSTATE_CLASS_ID);

  if (in.front().is_string()) {
    std::string cmd    = in.pop().to_string();
    if (check_cmd(cmd, "real", in, out, 0, 0, 0, 1)) {
      /*@INIT MDS = ('real')
        Build a model state for real unknowns.@*/
      md->set(new real_model_state);
    } else if (check_cmd(cmd, "complex", in, out, 0, 0, 0, 1)) {
      /*@INIT MDS = ('complex')
        Build a model state for complex unknowns.@*/
      md->set(new cplx_model_state);
    } else bad_cmd(cmd);
  } else if (in.front().is_mdbrick()) {
    /*@INIT MDS = ('.mdbrick', @tmdbrick B)
      Build a modelstate for the brick `B`.

      Selects the real or complex state from the complexity of `B`.@*/
    getfemint_mdbrick *b = in.pop().to_getfemint_mdbrick();
    if (!b->is_complex())
      md->set(new real_model_state(b->real_mdbrick()));
    else
      md->set(new cplx_model_state(b->cplx_mdbrick()));
  } else THROW_BADARG("expected a string or a mdbrick");

  if (in.remaining()) THROW_BADARG("too many arguments");
}
