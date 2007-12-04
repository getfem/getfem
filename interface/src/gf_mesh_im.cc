// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2006-2008 Yves Renard, Julien Pommier.
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

#include <getfemint_misc.h>
#include <getfemint_workspace.h>
#include <getfemint_mesh_im.h>
#include <getfemint_mesh.h>
#include <getfemint_mesh_levelset.h>
#include <getfem/getfem_mesh_im_level_set.h>

using namespace getfemint;

void gf_mesh_im_set_integ(getfem::mesh_im *mim, getfemint::mexargs_in& in);

/*MLABCOM
  FUNCTION MIM = gf_mesh_im(...)

  General constructor for @tmim object (integration methods on a mesh).

  * gf_mesh_im(mesh M [{integ IM|int IM_DEGREE}])

  Return a getfem handle to the newly created mesh_im object. For
  convenience, optional arguments (IM or IM_DEGREE) can be provided,
  in that case a call to gf_mesh_im_set(mim, 'integ', ..) is issued with these
  arguments.

  @INIT MESHIM:INIT('load')
  @INIT MESHIM:INIT('from string')
  @INIT MESHIM:INIT('clone')
  @INIT MESHIM:INIT('levelset')

  $Id$
MLABCOM*/

void gf_mesh_im(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) THROW_BADARG("Wrong number of input arguments");
  getfemint_mesh *mm = NULL;
  getfemint_mesh_im *mim = NULL;
  if (in.front().is_string()) {
    std::string cmd = in.pop().to_string();
    if (check_cmd(cmd, "load", in, out, 1, 2, 0, 1)) {
      /*@INIT MESHIM:INIT('load', fname[, @tmesh M]) 
	Load a @tmim from a file. 
	
	If the mesh M is not supplied (this kind of file does not
	store the mesh), then it is read from the file and its
	descriptor is returned as the second output argument. @*/
      std::string fname = in.pop().to_string();
      if (in.remaining()) mm = in.pop().to_getfemint_mesh();
      else {
	getfem::mesh *m = new getfem::mesh();
	m->read_from_file(fname);
	mm = getfemint_mesh::get_from(m);
      }
      mim = getfemint_mesh_im::new_from(mm);
      mim->mesh_im().read_from_file(fname);
    } else if (check_cmd(cmd, "from string", in, out, 1, 2, 0, 1)) {
      /*@INIT MESHIM:INIT('from string', str[, mesh M])
	Create a @tmim object from its string description.
	
	See also MESHIM:GET('char')
	@*/      
      std::stringstream ss(in.pop().to_string());
      if (in.remaining()) mm = in.pop().to_getfemint_mesh();
      else {
	getfem::mesh *m = new getfem::mesh();
	m->read_from_file(ss);
	mm = getfemint_mesh::get_from(m);
      }
      mim = getfemint_mesh_im::new_from(mm);
      mim->mesh_im().read_from_file(ss);
    } else if (check_cmd(cmd, "clone", in, out, 1, 1, 0, 1)) {
      /*@INIT MESHIM:INIT('clone', @tmim MIM2)
	Create a copy of a @tmim.
	@*/
      getfemint_mesh_im *mim2 = in.pop().to_getfemint_mesh_im();
      mm = object_to_mesh(workspace().object(mim2->linked_mesh_id()));
      mim = getfemint_mesh_im::new_from(mm);
      std::stringstream ss; /* not very elegant ! */
      mim2->mesh_im().write_to_file(ss);
      mim->mesh_im().read_from_file(ss);
    } else if (check_cmd(cmd, "levelset", in, out, 2, 4, 0, 1)) {
      /*@INIT MESHIM:INIT('levelset', @tmls MLS, @tstr WHERE, @integ IM [, @integ IMTIP])
	Build an integration method conformal to a partition defined
	implicitely by a levelset.

	The WHERE argument define the domain of integration with
	respect to the levelset, it has to be chosen among 'ALL',
	'INSIDE', 'OUTSIDE' and 'BOUNDARY'. 
	
	@*/
      getfemint_mesh_levelset *gmls = in.pop().to_getfemint_mesh_levelset();
      std::string swhere = in.pop().to_string();
      getfem::pintegration_method pim  = in.pop().to_integration_method();
      getfem::pintegration_method pim2 = 0;
      if (in.remaining()) pim2 = in.pop().to_integration_method();
      int where = 0;
      std::string csg_description;
      if (cmd_strmatch(swhere, "all")) 
	where = getfem::mesh_im_level_set::INTEGRATE_ALL;
      else {
	const char *slst[] = {"inside", "outside", "boundary", "all"};
	for (unsigned i=0; i < 4; ++i) {
	  if (cmd_strmatchn(swhere, slst[i], strlen(slst[i]))) {
	    csg_description.assign(swhere.begin() + strlen(slst[i]), swhere.end());
	    if (i == 0)      where = getfem::mesh_im_level_set::INTEGRATE_INSIDE;
	    else if (i == 1) where = getfem::mesh_im_level_set::INTEGRATE_OUTSIDE;
	    else if (i == 2) where = getfem::mesh_im_level_set::INTEGRATE_BOUNDARY;
	    else if (i == 3) where = getfem::mesh_im_level_set::INTEGRATE_ALL;
	  }
	}
      }
      if (where == 0) {
	THROW_BADARG("expecting 'inside', 'outside', 'boundary' or 'all'");
      }

      cerr << "csg_description: " << csg_description << "\n";

      getfem::mesh_im_level_set *mimls = 
	new getfem::mesh_im_level_set(gmls->mesh_levelset(),
				      where, pim, pim2);
      mimls->set_integration_method(mimls->linked_mesh().convex_index(), 1);
      if (csg_description.size()) {
	mimls->set_level_set_boolean_operations(csg_description);
      }
      mim = getfemint_mesh_im::get_from(mimls);
      workspace().set_dependance(mim, gmls);
      mimls->adapt();
    } else bad_cmd(cmd);
  } else {
    if (!out.narg_in_range(1, 1)) THROW_BADARG("Wrong number of output arguments");
    mm = in.pop().to_getfemint_mesh();
    mim = getfemint_mesh_im::new_from(mm);
    if (in.remaining()) {
      gf_mesh_im_set_integ(&mim->mesh_im(), in);
    }
    if (in.remaining()) THROW_BADARG("Wrong number of input arguments");
  }
  out.pop().from_object_id(mim->get_id(), MESHIM_CLASS_ID);
}

