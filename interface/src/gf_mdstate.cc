// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2005-2006 Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

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


/*MLABCOM

  FUNCTION M = gf_mdstate('complex' | 'real')
  General constructor for mdstate objects. Return a getfem handle to the newly
  created object.

  ``Model State'' variables store the state data for a set of model
  bricks. This includes the global tangent matrix, the right hand side
  and the constraints. There are two sorts of model states, the
  ``real'' and the ``complex'' models states.

  * MDS=gf_model_state(mdbrick B)
  Build a modelstate for the brick B (selects the real or
  complex state from the complexity of B).

  @INIT MDSTATE:INIT ('real')
  @INIT MDSTATE:INIT ('complex')
 
  $Id$
MLABCOM*/
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
      /*@INIT MDS=MDSTATE:INIT('real')
	Build a model state for real unknowns.
	@*/
      md->set(new real_model_state);
    } else if (check_cmd(cmd, "complex", in, out, 0, 0, 0, 1)) {
      /*@INIT MDS=MDSTATE:INIT('complex')
	Build a model state for complex unknowns.
	@*/
      md->set(new cplx_model_state);
    } else bad_cmd(cmd);
  } else if (in.front().is_mdbrick()) {
    getfemint_mdbrick *b = in.pop().to_getfemint_mdbrick();
    if (!b->is_complex()) 
      md->set(new real_model_state(b->real_mdbrick()));
    else 
      md->set(new cplx_model_state(b->cplx_mdbrick()));
  } else THROW_BADARG("expected a string or a mdbrick");

  if (in.remaining()) THROW_BADARG("too many arguments");
}
