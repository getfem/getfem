// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Yves Renard.
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

/**\file gf_model_set.cc
   \brief getfemint_model setter.
*/

#include <getfemint.h>
#include <getfemint_models.h>
#include <getfemint_mesh_fem.h>

using namespace getfemint;


/*MLABCOM

  FUNCTION M = gf_model_set(cmd, [, args])
  Modify a model object.

  @SET MODEL:SET('state')
  @SET MODEL:SET('variable')
  @SET MODEL:SET('clear')
  @SET MODEL:SET('add fem variable')
  @SET MODEL:SET('add variable')
  @SET MODEL:SET('add fem data')
  @SET MODEL:SET('add data')

MLABCOM*/

void gf_model_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) THROW_BADARG( "Wrong number of input arguments");

  getfemint_model *md  = in.pop().to_getfemint_model(true);
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "state", in, out, 1, 1, 0, 0)) {
    /*@SET MODEL:SET('state',@vec U)
      Update the internal state with the vector `U`.@*/
    if (!md->is_complex()) {
      darray st = in.pop().to_darray(-1);
      GMM_ASSERT1(st.size() == md->model().real_state().size(),
		  "Bad size in assigment");
      md->model().real_state().assign(st.begin(), st.end());
    } else {
      carray st = in.pop().to_carray(-1);
      GMM_ASSERT1(st.size() == md->model().real_state().size(),
		  "Bad size in assigment");
      md->model().complex_state().assign(st.begin(), st.end());
    }
  } else if (check_cmd(cmd, "clear", in, out, 0, 0, 0, 1)) {
    /*@SET MODEL:SET('clear')
    Clear the model.@*/
    md->clear();
  } else if (check_cmd(cmd, "add fem variable", in, out, 2, 3, 0, 0)) {
    /*@SET MODEL:SET('add fem variable', @str name, @tmf mf[, @int niter])
    Add a variable to the model linked to a @tmf. `name` is the variable name
    and `niter` is the optional number of copy of the variable for time
    integration schemes. @*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fem_variable(name, gfi_mf->mesh_fem(), niter);
  } else if (check_cmd(cmd, "add variable", in, out, 2, 3, 0, 0)) {
    /*@SET MODEL:SET('add variable', @str name, @int size[, @int niter])
    Add a variable to the model of constant size. `name` is the variable name
    and `niter` is the optional number of copy of the variable for time
    integration schemes. @*/
    std::string name = in.pop().to_string();
    size_type s = in.pop().to_integer();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fixed_size_variable(name, s, niter);
  } else if (check_cmd(cmd, "add fem data", in, out, 2, 4, 0, 0)) {
    /*@SET MODEL:SET('add fem data', @str name, @tmf mf[, @int qdim, @int niter])
    Add a data to the model linked to a @tmf. `name` is the data name,
    `qdim` is the optional dimension of the data over the @tmf
    and `niter` is the optional number of copy of the data for time
    integration schemes. @*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    dim_type qdim = 1;
    if (in.remaining()) qdim = dim_type(in.pop().to_integer(1,255));
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fem_data(name, gfi_mf->mesh_fem(), qdim, niter);
  } else if (check_cmd(cmd, "add data", in, out, 2, 3, 0, 0)) {
    /*@SET MODEL:SET('add data', @str name, @int size[, @int niter])
    Add a data to the model of constant size. `name` is the data name
    and `niter` is the optional number of copy of the data for time
    integration schemes. @*/
    std::string name = in.pop().to_string();
    size_type s = in.pop().to_integer();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fixed_size_data(name, s, niter);
  } else if (check_cmd(cmd, "variable", in, out, 2, 3, 0, 0)) {
    /*@SET V = MODEL:SET('variable', @str name, @vec U[, @int niter])
    Set the value of a variable or data.@*/
    std::string name = in.pop().to_string();
    if (!md->is_complex()) {
      darray st = in.pop().to_darray(-1);
      size_type niter = 0;
      if (in.remaining()) niter = in.pop().to_integer(0,10);
      GMM_ASSERT1(st.size() == md->model().real_variable(name, niter).size(),
		  "Bad size in assigment");
      md->model().real_variable(name, niter).assign(st.begin(), st.end());
    } else {
      carray st = in.pop().to_carray(-1);
      size_type niter = 0;
      if (in.remaining())
	niter = in.pop().to_integer(0,10) - config::base_index();
      GMM_ASSERT1(st.size() == md->model().real_variable(name, niter).size(),
		  "Bad size in assigment");
      md->model().complex_variable(name, niter).assign(st.begin(), st.end());
    }
  } else bad_cmd(cmd);
}
