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
#include <getfemint_mesh_fem.h>
#include <getfemint_mesh.h>
#include <getfem/getfem_mesh_fem_sum.h>
#include <getfemint_mesh_levelset.h>
#include <getfem/getfem_mesh_fem_level_set.h>
#include <getfem/getfem_partial_mesh_fem.h>

using namespace getfemint;

/*MLABCOM
  FUNCTION MF = gf_mesh_fem(...)

  General constructor for mesh_fem objects (Finite Element basis
  functions on a mesh).

  * gf_mesh_fem(mesh M [, int Qdim=1])

  Return a getfem handle to the newly created mesh_fem object.

  @INIT MESHFEM:INIT('load')
  @INIT MESHFEM:INIT('from string')
  @INIT MESHFEM:INIT('clone')
  @INIT MESHFEM:INIT('sum')
  @INIT MESHFEM:INIT('levelset')
  @INIT MESHFEM:INIT('partial')

  $Id$
MLABCOM*/

void gf_mesh_fem(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) THROW_BADARG("Wrong number of input arguments");
  getfemint_mesh *mm = NULL;
  getfemint_mesh_fem *mmf = NULL;
  unsigned q_dim = 1;
  if (in.front().is_string()) {
    std::string cmd = in.pop().to_string();
    if (check_cmd(cmd, "load", in, out, 1, 2, 0, 1)) {
      /*@INIT MESHFEM:INIT('load', @str fname[, @tmesh m])
      Load a @tmf from a file.

      If the mesh `m` is not supplied (this kind of file does not store
      the mesh), then it is read from the file `fname` and its descriptor
      is returned as the second output argument.@*/
      std::string fname = in.pop().to_string();
      if (in.remaining()) mm = in.pop().to_getfemint_mesh();
      else {
	getfem::mesh *m = new getfem::mesh();
	m->read_from_file(fname);
	mm = getfemint_mesh::get_from(m);
      }
      mmf = getfemint_mesh_fem::new_from(mm,q_dim);
      mmf->mesh_fem().read_from_file(fname);
    } else if (check_cmd(cmd, "from string", in, out, 1, 2, 0, 1)) {
      /*@INIT MESHFEM:INIT('from string', @str [, mesh M])
      Create a @tmf object from its string description.

      See also MESHFEM:GET('char')@*/
      std::stringstream ss(in.pop().to_string());
      if (in.remaining()) mm = in.pop().to_getfemint_mesh();
      else {
	getfem::mesh *m = new getfem::mesh();
	m->read_from_file(ss);
	mm = getfemint_mesh::get_from(m);
      }
      mmf = getfemint_mesh_fem::new_from(mm,q_dim);
      mmf->mesh_fem().read_from_file(ss);
    } else if (check_cmd(cmd, "clone", in, out, 1, 1, 0, 1)) {
      /*@INIT MESHFEM:INIT('clone', @tmf mf2)
      Create a copy of a @tmf.@*/
      getfemint_mesh_fem *mmf2 = in.pop().to_getfemint_mesh_fem();
      mm = object_to_mesh(workspace().object(mmf2->linked_mesh_id()));
      mmf = getfemint_mesh_fem::new_from(mm,q_dim);
      std::stringstream ss;
      mmf2->mesh_fem().write_to_file(ss);
      mmf->mesh_fem().read_from_file(ss);
    } else if (check_cmd(cmd, "sum", in, out, 1, -1, 0, 1)) {
      /*@INIT MESHFEM:INIT('sum', @tmf mf1, @tmf mf2[, @tmf mf3[, ...]])
      Create a @tmf that combines two (or more) @tmf's.

      All @tmf must share the same mesh (see FEM:INIT('interpolated_fem')
      to map a @tmf onto another).

      After that, you should not modify the FEM of `mf1`, `mf2` etc.@*/
      std::vector<const getfem::mesh_fem*> mftab;
      getfem::mesh_fem_sum *msum = 0;
      while (in.remaining()) {
	getfemint_mesh_fem *gfimf = in.pop().to_getfemint_mesh_fem();
	if (msum == 0) {
	  mm = object_to_mesh(workspace().object(gfimf->linked_mesh_id()));
	  msum = new getfem::mesh_fem_sum(gfimf->linked_mesh());

	  mmf = getfemint_mesh_fem::get_from(msum);
	}
	workspace().set_dependance(mmf, gfimf);
	mftab.push_back(&gfimf->mesh_fem());
      }
      msum->set_mesh_fems(mftab);
      msum->adapt();
    } else if (check_cmd(cmd, "levelset", in, out, 2, 2, 0, 1)) {
      /*@INIT MESHFEM:INIT('levelset', @tmf mf, @tls ls)
      Create a @tmf that is conformal to implicit surfaces defined by @tls's.@*/
      getfemint_mesh_fem *gmf = in.pop().to_getfemint_mesh_fem();
      getfemint_mesh_levelset *gmls = in.pop().to_getfemint_mesh_levelset();
      getfem::mesh_fem_level_set *mfls =
	new getfem::mesh_fem_level_set(gmls->mesh_levelset(),
				       gmf->mesh_fem());
      mmf = getfemint_mesh_fem::get_from(mfls);
      workspace().set_dependance(mmf, gmf);
      workspace().set_dependance(mmf, gmls);
      mfls->adapt();
    } else if (check_cmd(cmd, "partial", in, out, 2, 3, 0, 1)) {
      /*@INIT MESHFEM:INIT('partial', @tmf mf, @ivec DOFs[,@ivec RCVs])
      Build a restricted @tmf by keeping only a subset of the degrees of freedom of `mf`.

      If `RCVs` is given, no FEM will be put on the convexes listed
      in `RCVs`.@*/
      getfemint_mesh_fem *gmf = in.pop().to_getfemint_mesh_fem();
      dal::bit_vector doflst = in.pop().to_bit_vector();
      dal::bit_vector rcvlst;
      if (in.remaining()) rcvlst = in.pop().to_bit_vector();

      getfem::partial_mesh_fem *ppmf =
	new getfem::partial_mesh_fem(gmf->mesh_fem());
      ppmf->adapt(doflst, rcvlst);

      mmf = getfemint_mesh_fem::get_from(ppmf);
      workspace().set_dependance(mmf, gmf);
    } else bad_cmd(cmd);
    //out.pop().from_object_id(workspace().push_object(mmf), MESHFEM_CLASS_ID);
  } else {
    /*@INIT MESHFEM:INIT('.mesh', @tmesh m[, @int Qdim=1])
    Build a new @tmf object. The `Qdim` parameter is optional.@*/
    if (in.narg()  >  3) THROW_BADARG("Wrong number of input arguments");
    if (!out.narg_in_range(1, 1)) THROW_BADARG("Wrong number of output arguments");
    mm = in.pop().to_getfemint_mesh();
    if (in.remaining()) q_dim = in.pop().to_integer(0,256);
    if (in.remaining()) {
      unsigned q_dim2 = in.pop().to_integer(0,256);
      mmf = getfemint_mesh_fem::new_from(mm,q_dim * q_dim2);
      mmf->mesh_fem().set_qdim_mn(dim_type(q_dim), dim_type(q_dim2));
    } else
      mmf = getfemint_mesh_fem::new_from(mm,q_dim);
    //out.pop().from_object_id(workspace().push_object(mmf), MESHFEM_CLASS_ID);
  }
  out.pop().from_object_id(mmf->get_id(), MESHFEM_CLASS_ID);
}
