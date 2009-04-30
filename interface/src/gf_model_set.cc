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
#include <getfemint_workspace.h>
#include <getfemint_mesh_im.h>

using namespace getfemint;


/*MLABCOM

  FUNCTION M = gf_model_set(cmd, [, args])
  Modify a model object.

  @SET MODEL:SET('variable')
  @SET MODEL:SET('clear')
  @SET MODEL:SET('add fem variable')
  @SET MODEL:SET('add variable')
  @SET MODEL:SET('add mult on region')
  @SET MODEL:SET('add fem data')
  @SET MODEL:SET('add initialized fem data')
  @SET MODEL:SET('add data')
  @SET MODEL:SET('add initialized data')
  @SET MODEL:SET('add Laplacian brick')
  @SET MODEL:SET('add generic elliptic brick')
  @SET MODEL:SET('add source term brick')
  @SET MODEL:SET('add Dirichlet condition with multiplier')
  @SET MODEL:SET('add Dirichlet condition with penalization')
  @SET MODEL:SET('change penalization coeff Dirichlet')

MLABCOM*/

void gf_model_set(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 2) THROW_BADARG( "Wrong number of input arguments");

  getfemint_model *md  = in.pop().to_getfemint_model(true);
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "clear", in, out, 0, 0, 0, 1)) {
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
    workspace().set_dependance(md, gfi_mf);
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
  } else if (check_cmd(cmd, "add mult on region", in, out, 5, 6, 0, 0)) {
    /*@SET MODEL:SET('add mult on region', @str name, @tmf mf, @tmim mim, @str primalname,@int region[, @int niter])
    Add a particular variable linked to a fem beeing a multiplier with
    respect to a primal variable. The dof will be filtered with a mass
    matrix (computed with the integration method `mim`) to retain only
    linearly independant constraints on the primal
    variable. niter is the number of version of the data stored, for time
    integration schemes. @*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string primalname = in.pop().to_string();
    size_type region = in.pop().to_integer();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_mult_on_region(name, gfi_mf->mesh_fem(),
				   gfi_mim->mesh_im(), primalname,
				   region, niter);
  } else if (check_cmd(cmd, "add fem data", in, out, 2, 4, 0, 0)) {
    /*@SET MODEL:SET('add fem data', @str name, @tmf mf[, @int qdim, @int niter])
    Add a data to the model linked to a @tmf. `name` is the data name,
    `qdim` is the optional dimension of the data over the @tmf
    and `niter` is the optional number of copy of the data for time
    integration schemes. @*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    workspace().set_dependance(md, gfi_mf);
    dim_type qdim = 1;
    if (in.remaining()) qdim = dim_type(in.pop().to_integer(1,255));
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fem_data(name, gfi_mf->mesh_fem(), qdim, niter);
  } else if (check_cmd(cmd, "add initialized fem data", in, out, 3, 3, 0, 0)) {
    /*@SET MODEL:SET('add initialized fem data', @str name, @tmf mf, @vec V)
    Add a data to the model linked to a @tmf. `name` is the data name.
    The data is initiakized with `V`. The data can be a scalar or vector field.
    @*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    workspace().set_dependance(md, gfi_mf);
    if (!md->is_complex()) {
      darray st = in.pop().to_darray(-1);
      std::vector<double> V(st.begin(), st.end());
      md->model().add_initialized_fem_data(name, gfi_mf->mesh_fem(), V);
    } else {
      carray st = in.pop().to_carray(-1);
      std::vector<std::complex<double> > V(st.begin(), st.end());
      md->model().add_initialized_fem_data(name, gfi_mf->mesh_fem(), V);
    }
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
  } else if (check_cmd(cmd, "add initialized data", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('add initialized data', @str name, @vec V)
    Add a fixed size data to the model linked to a @tmf.
    `name` is the data name, `V` is the value of the data.
    @*/
    std::string name = in.pop().to_string();
    if (!md->is_complex()) {
      darray st = in.pop().to_darray(-1);
      std::vector<double> V(st.begin(), st.end());
      md->model().add_initialized_fixed_size_data(name, V);
    } else {
      carray st = in.pop().to_carray(-1);
      std::vector<std::complex<double> > V(st.begin(), st.end());
      md->model().add_initialized_fixed_size_data(name, V);
    }
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
      md->model().set_real_variable(name, niter).assign(st.begin(), st.end());
    } else {
      carray st = in.pop().to_carray(-1);
      size_type niter = 0;
      if (in.remaining())
	niter = in.pop().to_integer(0,10) - config::base_index();
      GMM_ASSERT1(st.size() == md->model().real_variable(name, niter).size(),
		  "Bad size in assigment");
      md->model().set_complex_variable(name, niter).assign(st.begin(), st.end());
    }
  } else if (check_cmd(cmd, "add Laplacian brick", in, out, 2, 3, 0, 1)) {
    /*@SET V = MODEL:SET('add Laplacian brick', @tmim mim, @str varname[, @int region])
    add a Laplacian term to the model relatively to the variable `varname`.
    If this is a vector valued variable, the Laplacian term is added
    componentwise. `region` is an optional mesh region on which the term is added. If it is not specified, it is added on the whole mesh. @*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    workspace().set_dependance(md, gfi_mim);
    getfem::mesh_im &mim = gfi_mim->mesh_im();
    std::string varname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_Laplacian_brick(md->model(), mim, varname, region)
      + config::base_index();
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add generic elliptic brick", in, out, 3, 4, 0, 1)) {
    /*@SET V = MODEL:SET('add generic elliptic brick', @tmim mim, @str varname, @str dataname[, @int region])
    add a generic elliptic term to the model relatively to the variable `varname`.
    The shape of the elliptic
    term depends both on the variable and the data. This corresponds to a
    term $-\text{div}(a\nabla u)$ where $a$ is the data and $u$ the variable.
    The data can be a scalar, a matrix or an order four tensor. The variable
    can be vector valued or not. If the data is a scalar or a matrix and
    the variable is vector valued then the term is added componentwise.
    An order four tensor data is allowed for vector valued variable only.
    The data can be constant or describbed on a fem. Of course, when
    the data is a tensor describe on a finite element method (a tensor
    field) the data can be a huge vector. The components of the
    matrix/tensor have to be stored with the fortran order (columnwise) in
    the data vector (compatibility with blas). The symmetry of the given
    matrix/tensor is not verified (but assumed).
    If this is a vector valued variable, the Laplacian term is added
    componentwise. `region` is an optional mesh region on which the term is 
    added. If it is not specified, it is added on the whole mesh. @*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    workspace().set_dependance(md, gfi_mim);
    getfem::mesh_im &mim = gfi_mim->mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_generic_elliptic_brick(md->model(), mim, varname,
					   dataname, region)
      + config::base_index();
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add source term brick", in, out, 3, 5, 0, 1)) {
    /*@SET V = MODEL:SET('add source term brick', @tmim mim, @str varname, @str dataname[, @int region])
    add a source term to the model relatively to the variable `varname`.
    The source term is represented by the data `dataname` which could be
    constant or described on a fem.  `region` is an optional mesh region
    on which the term is added. An additional optional data `directdataname`
    can be provided. The corresponding data vector will be directly added
    to the right hand side without assembly. @*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    workspace().set_dependance(md, gfi_mim);
    getfem::mesh_im &mim = gfi_mim->mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    std::string directdataname;
    if (in.remaining()) directdataname = in.pop().to_string();
    size_type ind
      = getfem::add_source_term_brick(md->model(), mim, varname,
				      dataname, region, directdataname)
      + config::base_index();
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add Dirichlet condition with multiplier", in, out, 4, 5, 0, 1)) {
    /*@SET V = MODEL:SET('add Dirichlet condition with multiplier', @tmim mim, @str varname, mult_description, @int region[, @str dataname])
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`. This region should be a boundary. The Dirichlet
    condition is prescribed with a multiplier variable described by
    `mult_description`. If `mult_description` is a string this is assumed
    to be the variable name correpsonding to the multiplier (which should be
    first declared as a multiplier  variable on the mesh region in the model).
    If it is a finite element method (mesh_fem object) then a multiplier
    variable will be added to the model and build on this finite element
    method (it will be restricted to the mesh region `region` and eventually
    some conflicting dofs with some other multiplier variables will be
    suppressed). If it is an integer, then a  multiplier variable will be
    added to the model and build on a classical finite element of degree
    that integer. `dataname` is the optional
    right hand side of  the Dirichlet condition. It could be constant or
    described on a fem; scalar or vector valued, depending on the variable
    on which the Dirichlet condition is prescribed. Return the brick index
    in the model. @*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    workspace().set_dependance(md, gfi_mim);
    getfem::mesh_im &mim = gfi_mim->mesh_im();
    std::string varname = in.pop().to_string();
    int version = 0;
    size_type degree = 0;
    std::string multname;
    getfem::mesh_fem *mf_mult = 0;
    mexarg_in argin = in.pop();

    if (argin.is_integer()) {
      degree = argin.to_integer();
      version = 1;
    } else if (argin.is_string()) {
      multname = argin.to_string();
      version = 2;
    } else {
      getfemint_mesh_fem *gfi_mf = argin.to_getfemint_mesh_fem();
      mf_mult = &(gfi_mf->mesh_fem());
      version = 3;
    }
    size_type region = in.pop().to_integer();
    std::string dataname;
    if (in.remaining()) dataname = in.pop().to_string();
    size_type ind = config::base_index();
    switch(version) {
    case 1:  ind += getfem::add_Dirichlet_condition_with_multipliers
	(md->model(), mim, varname, dim_type(degree), region, dataname);
      break;
    case 2:  ind += getfem::add_Dirichlet_condition_with_multipliers
	(md->model(), mim, varname, multname, region, dataname);
      break;
    case 3:  ind += getfem::add_Dirichlet_condition_with_multipliers
	(md->model(), mim, varname, *mf_mult, region, dataname);
      break;
    }
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add Dirichlet condition with penalization", in, out, 4, 5, 0, 1)) {
    /*@SET V = MODEL:SET('add Dirichlet condition with penalization', @tmim mim, @str varname, @int region, @scalar coeff[, @str dataname])
      Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`. This region should be a boundary. The Dirichlet
      condition is prescribed with penalization. The penalization coefficient
      is intially `coeff` and will be added to the data of
      the model. `dataname` is the optional
      right hand side of  the Dirichlet condition. It could be constant or
      described on a fem; scalar or vector valued, depending on the variable
      on which the Dirichlet condition is prescribed. Return the brick index
      in the model.
      @*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    workspace().set_dependance(md, gfi_mim);
    getfem::mesh_im &mim = gfi_mim->mesh_im();
    std::string varname = in.pop().to_string();
    size_type region = in.pop().to_integer();
    double coeff = in.pop().to_scalar();
    std::string dataname;
    if (in.remaining()) dataname = in.pop().to_string();
    size_type ind = config::base_index();
    ind += getfem::add_Dirichlet_condition_with_penalization
      (md->model(), mim, varname, coeff, region, dataname);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "change penalization coeff Dirichlet", in, out, 2, 2, 0, 0)) {
    /*@SET V = MODEL:SET('change penalization coeff Dirichlet', @int ind_brick, @scalar coeff)
      Change the penalization coefficient of a Dirichlet condition with
      penalization brick. If the brick is not of this kind,
      this function has an undefined behavior.
      @*/
    size_type ind_brick = in.pop().to_integer();
    double coeff = in.pop().to_scalar();
    getfem::change_penalization_coeff_Dirichlet(md->model(), ind_brick, coeff);
  } else bad_cmd(cmd);
}
