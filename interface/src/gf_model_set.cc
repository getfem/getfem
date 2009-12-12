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
#include <getfemint_gsparse.h>
#include <getfem/getfem_Coulomb_friction.h>

using namespace getfemint;


/*MLABCOM

  FUNCTION M = gf_model_set(cmd, [, args])
  Modify a model object.

  @SET MODEL:SET('variable')
  @SET MODEL:SET('clear')
  @SET MODEL:SET('add fem variable')
  @SET MODEL:SET('add variable')
  @SET MODEL:SET('resize variable')
  @SET MODEL:SET('add multiplier')
  @SET MODEL:SET('add fem data')
  @SET MODEL:SET('add initialized fem data')
  @SET MODEL:SET('add data')
  @SET MODEL:SET('add initialized data')
  @SET MODEL:SET('to variables')
  @SET MODEL:SET('add Laplacian brick')
  @SET MODEL:SET('add generic elliptic brick')
  @SET MODEL:SET('add source term brick')
  @SET MODEL:SET('add normal source term brick')
  @SET MODEL:SET('add Dirichlet condition with multipliers')
  @SET MODEL:SET('add Dirichlet condition with penalization')
  @SET MODEL:SET('add generalized Dirichlet condition with multipliers')
  @SET MODEL:SET('add generalized Dirichlet condition with penalization')
  @SET MODEL:SET('change penalization coeff')
  @SET MODEL:SET('add Helmholtz brick')
  @SET MODEL:SET('add Fourier Robin brick')
  @SET MODEL:SET('add constraint with multipliers')
  @SET MODEL:SET('add constraint with penalization')
  @SET MODEL:SET('add explicit matrix')
  @SET MODEL:SET('add explicit rhs')
  @SET MODEL:SET('set private matrix')
  @SET MODEL:SET('set private rhs')
  @SET MODEL:SET('disable bricks')
  @SET MODEL:SET('unable bricks')
  @SET MODEL:SET('add isotropic linearized elasticity brick')
  @SET MODEL:SET('add linear incompressibility brick')
  @SET MODEL:SET('add mass brick')
  @SET MODEL:SET('add basic d on dt brick')
  @SET MODEL:SET('add basic d2 on dt2 brick')
  @SET MODEL:SET('add theta method dispatcher')
  @SET MODEL:SET('velocity update for order two theta method')
  @SET MODEL:SET('add midpoint dispatcher')
  @SET MODEL:SET('velocity update for Newmark scheme')
  @SET MODEL:SET('first iter')
  @SET MODEL:SET('next iter')
  @SET MODEL:SET('add basic contact brick')
  @SET MODEL:SET('contact brick set BN')
  @SET MODEL:SET('contact brick set BT')
  @SET MODEL:SET('add contact with rigid obstacle brick')
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
    Add a variable to the model linked to a @tmf. `name` is the variable
    name and `niter` is the optional number of version of the data stored,
    for time integration schemes.@*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fem_variable(name, gfi_mf->mesh_fem(), niter);
    workspace().set_dependance(md, gfi_mf);
  } else if (check_cmd(cmd, "add variable", in, out, 2, 3, 0, 0)) {
    /*@SET MODEL:SET('add variable', @str name, @int size[, @int niter])
    Add a variable to the model of constant size. `name` is the variable
    name and `niter` is the optional number of version of the data stored,
    for time integration schemes. @*/
    std::string name = in.pop().to_string();
    size_type s = in.pop().to_integer();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fixed_size_variable(name, s, niter);
  } else if (check_cmd(cmd, "resize variable", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('resize variable', @str name, @int size)
    Resize a  constant size variable of the model. `name` is the variable
    name. @*/
    std::string name = in.pop().to_string();
    size_type s = in.pop().to_integer();
    md->model().resize_fixed_size_variable(name, s);
  } else if (check_cmd(cmd, "add multiplier", in, out, 3, 4, 0, 0)) {
    /*@SET MODEL:SET('add multiplier', @str name, @tmf mf, @str primalname[, @int niter])
    Add a particular variable linked to a fem being a multiplier with
    respect to a primal variable. The dof will be filtered with the
    ``gmm::range_basis`` function applied on the terms of the model
    which link the multiplier and the primal variable. This in order to
    retain only linearly independant constraints on the primal variable.
    Optimized for boundary multipliers. `niter` is the optional number
    of version of the data stored, for time integration schemes. @*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    std::string primalname = in.pop().to_string();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_multiplier(name, gfi_mf->mesh_fem(), primalname, niter);
    workspace().set_dependance(md, gfi_mf);
  } else if (check_cmd(cmd, "add fem data", in, out, 2, 4, 0, 0)) {
    /*@SET MODEL:SET('add fem data', @str name, @tmf mf[, @int qdim[, @int niter]])
    Add a data to the model linked to a @tmf. `name` is the data name,
    `qdim` is the optional dimension of the data over the @tmf and
    `niter` is the optional number of version of the data stored,
    for time integration schemes. @*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    dim_type qdim = 1;
    if (in.remaining()) qdim = dim_type(in.pop().to_integer(1,255));
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fem_data(name, gfi_mf->mesh_fem(), qdim, niter);
    workspace().set_dependance(md, gfi_mf);
  } else if (check_cmd(cmd, "add initialized fem data", in, out, 3, 3, 0, 0)) {
    /*@SET MODEL:SET('add initialized fem data', @str name, @tmf mf, @vec V)
    Add a data to the model linked to a @tmf. `name` is the data name.
    The data is initiakized with `V`. The data can be a scalar or vector
    field.@*/
    std::string name = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    if (!md->is_complex()) {
      darray st = in.pop().to_darray();
      std::vector<double> V(st.begin(), st.end());
      md->model().add_initialized_fem_data(name, gfi_mf->mesh_fem(), V);
    } else {
      carray st = in.pop().to_carray();
      std::vector<std::complex<double> > V(st.begin(), st.end());
      md->model().add_initialized_fem_data(name, gfi_mf->mesh_fem(), V);
    }
    workspace().set_dependance(md, gfi_mf);
  } else if (check_cmd(cmd, "add data", in, out, 2, 3, 0, 0)) {
    /*@SET MODEL:SET('add data', @str name, @int size[, @int niter])
    Add a data to the model of constant size. `name` is the data name
    and `niter` is the optional number of version of the data stored,
    for time integration schemes. @*/
    std::string name = in.pop().to_string();
    size_type s = in.pop().to_integer();
    size_type niter = 1;
    if (in.remaining()) niter = in.pop().to_integer(1,10);
    md->model().add_fixed_size_data(name, s, niter);
  } else if (check_cmd(cmd, "add initialized data", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('add initialized data', @str name, @vec V)
    Add a fixed size data to the model linked to a @tmf.
    `name` is the data name and `V` is the value of the data.@*/
    std::string name = in.pop().to_string();
    if (!md->is_complex()) {
      darray st = in.pop().to_darray();
      std::vector<double> V(st.begin(), st.end());
      md->model().add_initialized_fixed_size_data(name, V);
    } else {
      carray st = in.pop().to_carray();
      std::vector<std::complex<double> > V(st.begin(), st.end());
      md->model().add_initialized_fixed_size_data(name, V);
    }
  } else if (check_cmd(cmd, "variable", in, out, 2, 3, 0, 0)) {
    /*@SET MODEL:SET('variable', @str name, @vec V[, @int niter])
    Set the value of a variable or data. `name` is the data name
    and `niter` is the optional number of version of the data stored,
    for time integration schemes.@*/
    std::string name = in.pop().to_string();
    if (!md->is_complex()) {
      darray st = in.pop().to_darray();
      size_type niter = 0;
      if (in.remaining()) niter = in.pop().to_integer(0,10) - config::base_index();
      GMM_ASSERT1(st.size() == md->model().real_variable(name, niter).size(),
                  "Bad size in assigment");
      md->model().set_real_variable(name, niter).assign(st.begin(), st.end());
    } else {
      carray st = in.pop().to_carray();
      size_type niter = 0;
      if (in.remaining())
        niter = in.pop().to_integer(0,10) - config::base_index();
      GMM_ASSERT1(st.size() == md->model().complex_variable(name, niter).size(),
                  "Bad size in assigment");
      md->model().set_complex_variable(name, niter).assign(st.begin(), st.end());
    }
  } else if (check_cmd(cmd, "to variables", in, out, 1, 1, 0, 0)) {
    /*@SET MODEL:SET('to variables', @vec V)
    Set the value of the variables of the model with the vector `V`.
    Typically, the vector `V` results of the solve of the tangent
    linear system (usefull to solve your problem with you own solver).@*/
    if (!md->is_complex()) {
      darray st = in.pop().to_darray(-1);
      std::vector<double> V;
      V.assign(st.begin(), st.end());
      md->model().to_variables(V);
    } else {
      carray st = in.pop().to_carray(-1);
      std::vector<std::complex<double> > V;
      V.assign(st.begin(), st.end());
      md->model().to_variables(V);
    }
  } else if (check_cmd(cmd, "add Laplacian brick", in, out, 2, 3, 0, 1)) {
    /*@SET ind = MODEL:SET('add Laplacian brick', @tmim mim, @str varname[, @int region])
    Add a Laplacian term to the model relatively to the variable `varname`.
    If this is a vector valued variable, the Laplacian term is added
    componentwise. `region` is an optional mesh region on which the term
    is added. If it is not specified, it is added on the whole mesh. Return
    the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_Laplacian_brick(md->model(), gfi_mim->mesh_im(), varname, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add generic elliptic brick", in, out, 3, 4, 0, 1)) {
    /*@SET ind = MODEL:SET('add generic elliptic brick', @tmim mim, @str varname, @str dataname[, @int region])
    Add a generic elliptic term to the model relatively to the variable `varname`.
    The shape of the elliptic term depends both on the variable and the data.
    This corresponds to a term
    @PYTHON{:math:`-\\text{div}(a\\nabla u)`}@MATLAB{$-\text{div}(a\nabla u)$}
    where @PYTHON{:math:`a`}@MATLAB{$a$} is the data and
    @PYTHON{:math:`u`}@MATLAB{$u$} the variable. The data can be a scalar,
    a matrix or an order four tensor. The variable can be vector valued or
    not. If the data is a scalar or a matrix and the variable is vector
    valued then the term is added componentwise. An order four tensor data
    is allowed for vector valued variable only. The data can be constant or
    describbed on a fem. Of course, when the data is a tensor describe on a
    finite element method (a tensor field) the data can be a huge vector.
    The components of the matrix/tensor have to be stored with the fortran
    order (columnwise) in the data vector (compatibility with blas). The
    symmetry of the given matrix/tensor is not verified (but assumed). If
    this is a vector valued variable, the Laplacian term is added
    componentwise. `region` is an optional mesh region on which the term is
    added. If it is not specified, it is added on the whole mesh. Return the
    brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_generic_elliptic_brick(md->model(), gfi_mim->mesh_im(),
                                           varname, dataname, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add source term brick", in, out, 3, 5, 0, 1)) {
    /*@SET ind = MODEL:SET('add source term brick', @tmim mim, @str varname, @str dataname[, @int region[, @str directdataname]])
    Add a source term to the model relatively to the variable `varname`.
    The source term is represented by the data `dataname` which could be
    constant or described on a fem. `region` is an optional mesh region
    on which the term is added. An additional optional data `directdataname`
    can be provided. The corresponding data vector will be directly added
    to the right hand side without assembly. Return the brick index in the
    model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    std::string directdataname;
    if (in.remaining()) directdataname = in.pop().to_string();
    size_type ind
      = getfem::add_source_term_brick(md->model(), gfi_mim->mesh_im(),
                                      varname, dataname, region, directdataname)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add normal source term brick", in, out, 4, 4, 0, 1)) {
    /*@SET ind = MODEL:SET('add normal source term brick', @tmim mim, @str varname, @str dataname, @int region)
    Add a source term on the variable `varname` on a boundary `region`.
    This region should be a boundary. The source term is represented by the
    data `dataname` which could be constant or described on a fem. A scalar
    product with the outward normal unit vector to the boundary is performed.
    The main aim of this brick is to represent a Neumann condition with a
    vector data without performing the scalar product with the normal as a
    pre-processing. Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname = in.pop().to_string();
    size_type region = in.pop().to_integer();
    size_type ind
      = getfem::add_normal_source_term_brick(md->model(), gfi_mim->mesh_im(),
                                             varname, dataname, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add Dirichlet condition with multipliers", in, out, 4, 5, 0, 1)) {
    /*@SET ind = MODEL:SET('add Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region[, @str dataname])
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`. This region should be a boundary. The Dirichlet
    condition is prescribed with a multiplier variable described by
    `mult_description`. If `mult_description` is a string this is assumed
    to be the variable name correpsonding to the multiplier (which should be
    first declared as a multiplier variable on the mesh region in the model).
    If it is a finite element method (mesh_fem object) then a multiplier
    variable will be added to the model and build on this finite element
    method (it will be restricted to the mesh region `region` and eventually
    some conflicting dofs with some other multiplier variables will be
    suppressed). If it is an integer, then a  multiplier variable will be
    added to the model and build on a classical finite element of degree
    that integer. `dataname` is the optional right hand side of  the
    Dirichlet condition. It could be constant or described on a fem; scalar
    or vector valued, depending on the variable on which the Dirichlet
    condition is prescribed. Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    int version = 0;
    size_type degree = 0;
    std::string multname;
    getfemint_mesh_fem *gfi_mf = 0;
    mexarg_in argin = in.pop();
    if (argin.is_integer()) {
      degree = argin.to_integer();
      version = 1;
    } else if (argin.is_string()) {
      multname = argin.to_string();
      version = 2;
    } else {
      gfi_mf = argin.to_getfemint_mesh_fem();
      version = 3;
    }
    size_type region = in.pop().to_integer();
    std::string dataname;
    if (in.remaining()) dataname = in.pop().to_string();
    size_type ind = config::base_index();
    switch(version) {
    case 1:  ind += getfem::add_Dirichlet_condition_with_multipliers
        (md->model(), gfi_mim->mesh_im(), varname, dim_type(degree), region, dataname);
      break;
    case 2:  ind += getfem::add_Dirichlet_condition_with_multipliers
        (md->model(), gfi_mim->mesh_im(), varname, multname, region, dataname);
      break;
    case 3:  ind += getfem::add_Dirichlet_condition_with_multipliers
        (md->model(), gfi_mim->mesh_im(), varname, gfi_mf->mesh_fem(), region, dataname);
        workspace().set_dependance(md, gfi_mf);
      break;
    }
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add Dirichlet condition with penalization", in, out, 4, 5, 0, 1)) {
    /*@SET ind = MODEL:SET('add Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region[, @str dataname])
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`. This region should be a boundary. The Dirichlet
    condition is prescribed with penalization. The penalization coefficient
    is intially `coeff` and will be added to the data of the model.
    `dataname` is the optional right hand side of the Dirichlet condition.
    It could be constant or described on a fem; scalar or vector valued,
    depending on the variable on which the Dirichlet condition is prescribed.
    Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    double coeff = in.pop().to_scalar();
    size_type region = in.pop().to_integer();
    std::string dataname;
    if (in.remaining()) dataname = in.pop().to_string();
    size_type ind = config::base_index();
    ind += getfem::add_Dirichlet_condition_with_penalization
      (md->model(), gfi_mim->mesh_im(), varname, coeff, region, dataname);
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add generalized Dirichlet condition with multipliers", in, out, 6, 6, 0, 1)) {
    /*@SET ind = MODEL:SET('add generalized Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region, @str dataname, @str Hname)
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`.  This version is for vector field.
    It prescribes a condition @PYTHON{:math:`Hu = r`}@MATLAB{$Hu = r$}
    where `H` is a matrix field. The region should be a boundary. The Dirichlet
    condition is prescribed with a multiplier variable described by
    `mult_description`. If `mult_description` is a string this is assumed
    to be the variable name correpsonding to the multiplier (which should be
    first declared as a multiplier variable on the mesh region in the model).
    If it is a finite element method (mesh_fem object) then a multiplier
    variable will be added to the model and build on this finite element
    method (it will be restricted to the mesh region `region` and eventually
    some conflicting dofs with some other multiplier variables will be
    suppressed). If it is an integer, then a  multiplier variable will be
    added to the model and build on a classical finite element of degree
    that integer. `dataname` is the right hand side of  the
    Dirichlet condition. It could be constant or described on a fem; scalar
    or vector valued, depending on the variable on which the Dirichlet
    condition is prescribed. `Hname' is the data
    corresponding to the matrix field `H`.
    Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    int version = 0;
    size_type degree = 0;
    std::string multname;
    getfemint_mesh_fem *gfi_mf = 0;
    mexarg_in argin = in.pop();
    if (argin.is_integer()) {
      degree = argin.to_integer();
      version = 1;
    } else if (argin.is_string()) {
      multname = argin.to_string();
      version = 2;
    } else {
      gfi_mf = argin.to_getfemint_mesh_fem();
      version = 3;
    }
    size_type region = in.pop().to_integer();
    std::string dataname = in.pop().to_string();
    std::string Hname = in.pop().to_string();
    size_type ind = config::base_index();
    switch(version) {
    case 1:  ind += getfem::add_generalized_Dirichlet_condition_with_multipliers
        (md->model(), gfi_mim->mesh_im(), varname, dim_type(degree), region, dataname, Hname);
      break;
    case 2:  ind += getfem::add_generalized_Dirichlet_condition_with_multipliers
        (md->model(), gfi_mim->mesh_im(), varname, multname, region, dataname, Hname);
      break;
    case 3:  ind += getfem::add_generalized_Dirichlet_condition_with_multipliers
        (md->model(), gfi_mim->mesh_im(), varname, gfi_mf->mesh_fem(), region, dataname, Hname);
        workspace().set_dependance(md, gfi_mf);
      break;
    }
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add generalized Dirichlet condition with penalization", in, out, 6, 6, 0, 1)) {
    /*@SET ind = MODEL:SET('add generalized Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region, @str dataname, @str Hname)
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`. This version is for vector field.
    It prescribes a condition @PYTHON{:math:`Hu = r`}@MATLAB{$Hu = r$}
    where `H` is a matrix field.
    The region should be a boundary. The Dirichlet
    condition is prescribed with penalization. The penalization coefficient
    is intially `coeff` and will be added to the data of the model.
    `dataname` is the right hand side of the Dirichlet condition.
    It could be constant or described on a fem; scalar or vector valued,
    depending on the variable on which the Dirichlet condition is prescribed.
    `Hname' is the data
    corresponding to the matrix field `H`. It has to be a constant matrix
      or described on a scalar fem. 
    Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    double coeff = in.pop().to_scalar();
    size_type region = in.pop().to_integer();
    std::string dataname= in.pop().to_string();
    std::string Hname= in.pop().to_string();
    size_type ind = config::base_index();
    ind += getfem::add_generalized_Dirichlet_condition_with_penalization
      (md->model(), gfi_mim->mesh_im(), varname, coeff, region,
       dataname, Hname);
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "change penalization coeff", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('change penalization coeff', @int ind_brick, @scalar coeff)
    Change the penalization coefficient of a Dirichlet condition with
    penalization brick. If the brick is not of this kind, this
    function has an undefined behavior.@*/
    size_type ind_brick = in.pop().to_integer();
    double coeff = in.pop().to_scalar();
    getfem::change_penalization_coeff(md->model(), ind_brick, coeff);
  } else if (check_cmd(cmd, "add Helmholtz brick", in, out, 3, 4, 0, 1)) {
    /*@SET ind = MODEL:SET('add Helmholtz brick', @tmim mim, @str varname, @str dataname[, @int region])
    Add a Helmholtz term to the model relatively to the variable `varname`.
    `dataname` should contain the wave number. `region` is an optional mesh
    region on which the term is added. If it is not specified, it is added
    on the whole mesh. Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_Helmholtz_brick(md->model(), gfi_mim->mesh_im(),
                                    varname, dataname, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add Fourier Robin brick", in, out, 4, 4, 0, 1)) {
    /*@SET ind = MODEL:SET('add Fourier Robin brick', @tmim mim, @str varname, @str dataname, @int region)
    Add a Fourier-Robin term to the model relatively to the variable
    `varname`. This corresponds to a weak term of the form
    @PYTHON{:math:`\\int (qu).v`}@MATLAB{$\int (qu).v$}. `dataname`
    should contain the parameter @PYTHON{:math:`q`}@MATLAB{$q$} of
    the Fourier-Robin condition. `region` is the mesh region on which
    the term is added. Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_Fourier_Robin_brick(md->model(), gfi_mim->mesh_im(),
                                        varname,dataname, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add constraint with multipliers", in, out, 4, 4, 0, 1)) {
    /*@SET ind = MODEL:SET('add constraint with multipliers', @str varname, @str multname, @spmat B, @vec L)
    Add an additional explicit constraint on the variable `varname` thank to
    a multiplier `multname` peviously added to the model (should be a fixed
    size variable). The constraint is @PYTHON{:math:`BU=L`}@MATLAB{$BU=L$}
    with `B` being a rectangular sparse matrix. It is possible to change
    the constraint at any time whith the methods MODEL:SET('set private matrix')
    and MODEL:SET('set private rhs'). Return the brick index in the model.@*/
    std::string varname = in.pop().to_string();
    std::string multname = in.pop().to_string();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();
    if (B->is_complex() && !md->is_complex())
      THROW_BADARG("Complex constraint for a real model");
    if (!B->is_complex() && md->is_complex())
      THROW_BADARG("Real constraint for a complex model");

    size_type ind
      = getfem::add_constraint_with_multipliers(md->model(),varname,multname);

    if (md->is_complex()) {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_wsc());
      else
        THROW_BADARG("Constraint matrix should be a sparse matrix");
    } else {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_wsc());
      else
        THROW_BADARG("Constraint matrix should be a sparse matrix");
    }

    if (!md->is_complex()) {
      darray st = in.pop().to_darray();
      std::vector<double> V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    } else {
      carray st = in.pop().to_carray();
      std::vector<std::complex<double> > V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    }

    out.pop().from_integer(int(ind + config::base_index()));
  } else if (check_cmd(cmd, "add constraint with penalization", in, out, 4, 4, 0, 1)) {
    /*@SET ind = MODEL:SET('add constraint with penalization', @str varname, @scalar coeff, @spmat B, @vec L)
    Add an additional explicit penalized constraint on the variable `varname`.
    The constraint is $BU=L$ with `B` being a rectangular sparse matrix.
    Be aware that `B` should not contain a palin row, otherwise the whole
    tangent matrix will be plain. It is possible to change the constraint
    at any time whith the methods MODEL:SET('set private matrix')
    and MODEL:SET('set private rhs'). The method
    MODEL:SET('change penalization coeff') can be used. Return the brick
    index in the model.@*/
    std::string varname = in.pop().to_string();
    double coeff = in.pop().to_scalar();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();
    if (B->is_complex() && !md->is_complex())
      THROW_BADARG("Complex constraint for a real model");
    if (!B->is_complex() && md->is_complex())
      THROW_BADARG("Real constraint for a complex model");

    size_type ind
      = getfem::add_constraint_with_penalization(md->model(), varname, coeff);

    if (md->is_complex()) {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_wsc());
      else
        THROW_BADARG("Constraint matrix should be a sparse matrix");
    } else {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_wsc());
      else
        THROW_BADARG("Constraint matrix should be a sparse matrix");
    }

    if (!md->is_complex()) {
      darray st = in.pop().to_darray();
      std::vector<double> V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    } else {
      carray st = in.pop().to_carray();
      std::vector<std::complex<double> > V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    }

    out.pop().from_integer(int(ind + config::base_index()));
  } else if (check_cmd(cmd, "add explicit matrix", in, out, 3, 5, 0, 1)) {
    /*@SET ind = MODEL:SET('add explicit matrix', @str varname1, @str varname2, @spmat B[, @int issymmetric[, @int iscoercive]])
    Add a brick representing an explicit matrix to be added to the tangent
    linear system relatively to the variables 'varname1' and 'varname2'.
    The given matrix should have has many rows as the dimension of
    'varname1' and as many columns as the dimension of 'varname2'.
    If the two variables are different and if `issymmetric' is set to 1
    then the transpose of the matrix is also added to the tangent system
    (default is 0). Set `iscoercive` to 1 if the term does not affect the
    coercivity of the tangent system (default is 0). The matrix can be
    changed by the command MODEL:SET('set private matrix'). Return the
    brick index in the model.@*/
    std::string varname1 = in.pop().to_string();
    std::string varname2 = in.pop().to_string();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();
    bool issymmetric = false;
    bool iscoercive = false;
    if (in.remaining()) issymmetric = (in.pop().to_integer(0,1) != 1);
    if (!issymmetric && in.remaining())
      iscoercive = (in.pop().to_integer(0,1) != 1);

    size_type ind
      = getfem::add_explicit_matrix(md->model(), varname1, varname2,
                                    issymmetric, iscoercive);

    if (B->is_complex() && !md->is_complex())
      THROW_BADARG("Complex constraint for a real model");
    if (!B->is_complex() && md->is_complex())
      THROW_BADARG("Real constraint for a complex model");

    if (md->is_complex()) {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_wsc());
      else
      THROW_BADARG("Constraint matrix should be a sparse matrix");
    } else {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_wsc());
      else
      THROW_BADARG("Constraint matrix should be a sparse matrix");
    }

    out.pop().from_integer(int(ind + config::base_index()));
  } else if (check_cmd(cmd, "add explicit rhs", in, out, 2, 2, 0, 1)) {
    /*@SET ind = MODEL:SET('add explicit rhs', @str varname, @vec L)
    Add a brick representing an explicit right hand side to be added to
    the right hand side of the tangent linear system relatively to the
    variable 'varname'. The given rhs should have the same size than the
    dimension of 'varname'. The rhs can be changed by the command
    MODEL:SET('set private rhs'). Return the brick index in the model.@*/
    std::string varname = in.pop().to_string();
    size_type ind
      = getfem::add_explicit_rhs(md->model(), varname);

    if (!md->is_complex()) {
      darray st = in.pop().to_darray();
      std::vector<double> V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    } else {
      carray st = in.pop().to_carray();
      std::vector<std::complex<double> > V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    }

    out.pop().from_integer(int(ind + config::base_index()));
  } else if (check_cmd(cmd, "set private matrix", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('set private matrix', @int indbrick, @spmat B)
    For some specific bricks having an internal sparse matrix
    (explicit bricks: 'constraint brick' and 'explicit matrix brick'),
    set this matrix. @*/
    size_type ind = in.pop().to_integer();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();

    if (B->is_complex() && !md->is_complex())
      THROW_BADARG("Complex constraint for a real model");
    if (!B->is_complex() && md->is_complex())
      THROW_BADARG("Real constraint for a complex model");

    if (md->is_complex()) {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->cplx_wsc());
      else
      THROW_BADARG("Constraint matrix should be a sparse matrix");
    } else {
      if (B->storage()==gsparse::CSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_csc());
      else if (B->storage()==gsparse::WSCMAT)
        getfem::set_private_data_matrix(md->model(), ind, B->real_wsc());
      else
      THROW_BADARG("Constraint matrix should be a sparse matrix");
    }
  } else if (check_cmd(cmd, "set private rhs", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('set private rhs', @int indbrick, @vec B)
    For some specific bricks having an internal right hand side vector
    (explicit bricks: 'constraint brick' and 'explicit rhs brick'),
    set this rhs. @*/
    size_type ind = in.pop().to_integer();

    if (!md->is_complex()) {
      darray st = in.pop().to_darray();
      std::vector<double> V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    } else {
      carray st = in.pop().to_carray();
      std::vector<std::complex<double> > V(st.begin(), st.end());
      getfem::set_private_data_rhs(md->model(), ind, V);
    }

  } else if (check_cmd(cmd, "add isotropic linearized elasticity brick", in, out, 4, 5, 0, 1)) {
    /*@SET ind = MODEL:SET('add isotropic linearized elasticity brick', @tmim mim, @str varname, @str dataname_lambda, @str dataname_mu[, @int region])
    Add an isotropic linearized elasticity term to the model relatively to the
    variable `varname`. `dataname_lambda` and `dataname_mu` should contain the
    Lam\'e coefficients. `region` is an optional mesh region on which the term
    is added. If it is not specified, it is added on the whole mesh. Return the
    brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname_lambda = in.pop().to_string();
    std::string dataname_mu = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_isotropic_linearized_elasticity_brick
      (md->model(), gfi_mim->mesh_im(), varname, dataname_lambda, dataname_mu, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add linear incompressibility brick", in, out, 3, 5, 0, 1)) {
    /*@SET ind = MODEL:SET('add linear incompressibility brick', @tmim mim, @str varname, @str multname_pressure[, @int region[, @str dataname_coeff]])
    Add an linear incompressibility condition on `variable`. `multname_pressure`
    is a variable which represent the pressure. Be aware that an inf-sup
    condition between the finite element method describing the pressure and the
    primal variable has to be satisfied. `region` is an optional mesh region on
    which the term is added. If it is not specified, it is added on the whole mesh.
    `dataname_coeff` is an optional penalization coefficient for nearly
    incompressible elasticity for instance. In this case, it is the inverse
    of the Lam\'e coefficient $\lambda$. Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string multname = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    std::string dataname;
    if (in.remaining()) dataname = in.pop().to_string();
    size_type ind
      = getfem::add_linear_incompressibility
      (md->model(), gfi_mim->mesh_im(), varname, multname, region, dataname)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add mass brick", in, out, 2, 4, 0, 1)) {
    /*@SET ind = MODEL:SET('add mass brick', @tmim mim, @str varname[, @str dataname_rho[, @int region]])
    Add mass term to the model relatively to the variable `varname`. If specified,
    the data `dataname_rho` should contain the density (1 if omitted). `region` is
    an optional mesh region on which the term is added. If it is not specified, it
    is added on the whole mesh. Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname = in.pop().to_string();
    std::string dataname_rho;
    if (in.remaining()) dataname_rho = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_mass_brick
      (md->model(), gfi_mim->mesh_im(), varname, dataname_rho, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add basic d on dt brick", in, out, 3, 5, 0, 1)) {
    /*@SET ind = MODEL:SET('add basic d on dt brick', @tmim mim, @str varnameU, @str dataname_dt[, @str dataname_rho[, @int region]])
    Add the standard discretization of a first order time derivative on
    `varnameU`. The parameter `dataname_rho` is the density which could
    be omitted (the defaul value is 1). This brick should be used in
    addition to a time dispatcher for the other terms. Return the brick
    index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varnameU = in.pop().to_string();
    std::string varnamedt = in.pop().to_string();
    std::string dataname_rho;
    if (in.remaining()) dataname_rho = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_basic_d_on_dt_brick
      (md->model(), gfi_mim->mesh_im(), varnameU, varnamedt, dataname_rho, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add basic d2 on dt2 brick", in, out, 5, 7, 0,1)) {
    /*@SET ind = MODEL:SET('add basic d2 on dt2 brick', @tmim mim, @str varnameU,  @str datanameV, @str dataname_dt, @str dataname_alpha,[, @str dataname_rho[, @int region]])
    Add the standard discretization of a second order time derivative
    on `varnameU`. `datanameV` is a data represented on the same finite
    element method as U which represents the time derivative of U. The
    parameter `dataname_rho` is the density which could be omitted (the defaul
    value is 1). This brick should be used in addition to a time dispatcher for
    the other terms. The time derivative @PYTHON{:math:`v`}@MATLAB{$v$} of the
    variable @PYTHON{:math:`u`}@MATLAB{$u$} is preferably computed as a
    post-traitement which depends on each scheme. The parameter `dataname_alpha`
    depends on the time integration scheme. Return the brick index in the model.@*/
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varnameU = in.pop().to_string();
    std::string varnameV = in.pop().to_string();
    std::string varnamedt = in.pop().to_string();
    std::string varnamealpha = in.pop().to_string();
    std::string dataname_rho;
    if (in.remaining()) dataname_rho = in.pop().to_string();
    size_type region = size_type(-1);
    if (in.remaining()) region = in.pop().to_integer();
    size_type ind
      = getfem::add_basic_d2_on_dt2_brick
      (md->model(), gfi_mim->mesh_im(), varnameU,  varnameV, varnamedt, varnamealpha, dataname_rho, region)
      + config::base_index();
    workspace().set_dependance(md, gfi_mim);
    out.pop().from_integer(int(ind));
  } else if (check_cmd(cmd, "add theta method dispatcher", in, out, 2, 2, 0,0)) {
    /*@SET MODEL:SET('add theta method dispatcher', @ivec bricks_indices, @str theta)
    Add a theta-method time dispatcher to a list of bricks. For instance,
    a matrix term @PYTHON{:math:`K`}@MATLAB{$K$} will be replaced by
    @PYTHON{:math:`\\theta K U^{n+1} + (1-\\theta) K U^{n}`}
    @MATLAB{$\theta K U^{n+1} + (1-\theta) K U^{n}$}.@*/
    dal::bit_vector bv = in.pop().to_bit_vector();
    std::string datanametheta = in.pop().to_string();
    getfem::add_theta_method_dispatcher(md->model(), bv, datanametheta);
  } else if (check_cmd(cmd, "add midpoint dispatcher", in, out, 1, 1, 0,0)) {
    /*@SET MODEL:SET('add midpoint dispatcher', @ivec bricks_indices)
    Add a midpoint time dispatcher to a list of bricks. For instance, a
    nonlinear term @PYTHON{:math:`K(U)`}@MATLAB{$K(U)$} will be replaced by
    @PYTHON{:math:`K((U^{n+1} +  U^{n})/2)`}@MATLAB{$K((U^{n+1} +  U^{n})/2)$}.@*/
    dal::bit_vector bv = in.pop().to_bit_vector();
    getfem::add_midpoint_dispatcher(md->model(), bv);
  } else if (check_cmd(cmd, "velocity update for order two theta method", in, out, 4, 4, 0,0)) {
    /*@SET MODEL:SET('velocity update for order two theta method', @str varnameU,  @str datanameV, @str dataname_dt, @str dataname_theta)
    Function which udpate the velocity @PYTHON{:math:`v^{n+1}`}@MATLAB{$v^{n+1}$} after
    the computation of the displacement @PYTHON{:math:`u^{n+1}`}@MATLAB{$u^{n+1}$} and
    before the next iteration. Specific for theta-method and when the velocity is
    included in the data of the model. @*/
    std::string varnameU = in.pop().to_string();
    std::string varnameV = in.pop().to_string();
    std::string varnamedt = in.pop().to_string();
    std::string varnametheta = in.pop().to_string();
    velocity_update_for_order_two_theta_method
      (md->model(), varnameU, varnameV, varnamedt, varnametheta);
  } else if (check_cmd(cmd, "velocity update for Newmark scheme", in, out, 6, 6, 0,0)) {
    /*@SET MODEL:SET('velocity update for Newmark scheme', @int id2dt2_brick, @str varnameU,  @str datanameV, @str dataname_dt, @str dataname_twobeta, @str dataname_alpha)
    Function which udpate the velocity @PYTHON{:math:`v^{n+1}`}@MATLAB{$v^{n+1}$} after
    the computation of the displacement @PYTHON{:math:`u^{n+1}`}@MATLAB{$u^{n+1}$} and
    before the next iteration. Specific for Newmark scheme and when the velocity is
    included in the data of the model. This version inverts the mass matrix by a
    conjugate gradient.@*/
    size_type id2dt2 = in.pop().to_integer();
    std::string varnameU = in.pop().to_string();
    std::string varnameV = in.pop().to_string();
    std::string varnamedt = in.pop().to_string();
    std::string varnametwobeta = in.pop().to_string();
    std::string varnamegamma = in.pop().to_string();
    velocity_update_for_Newmark_scheme
      (md->model(), id2dt2, varnameU, varnameV, varnamedt,
       varnametwobeta, varnamegamma);
  } else if (check_cmd(cmd, "disable bricks", in, out, 1, 1, 0,0)) {
    /*@SET MODEL:SET('disable bricks', @ivec bricks_indices)
    Disable a brick (the brick will no longer participate to the
    building of the tangent linear system).@*/
    dal::bit_vector bv = in.pop().to_bit_vector();
    for (dal::bv_visitor ii(bv); !ii.finished(); ++ii)
      md->model().disable_brick(ii);
  } else if (check_cmd(cmd, "unable bricks", in, out, 1, 1, 0,0)) {
    /*@SET MODEL:SET('unable bricks', @ivec bricks_indices)
    Unable a disabled brick.@*/
    dal::bit_vector bv = in.pop().to_bit_vector();
    for (dal::bv_visitor ii(bv); !ii.finished(); ++ii)
      md->model().unable_brick(ii);
  } else if (check_cmd(cmd, "first iter", in, out, 0, 0, 0,0)) {
    /*@SET MODEL:SET('first iter')
    To be executed before the first iteration of a time integration scheme.@*/
    md->model().first_iter();
  } else if (check_cmd(cmd, "next iter", in, out, 0, 0, 0,0)) {
    /*@SET MODEL:SET('next iter')
    To be executed at the end of each iteration of a time integration scheme.@*/
    md->model().next_iter();
  } else if (check_cmd(cmd, "add basic contact brick", in, out, 4, 7, 0, 1)) {
    /*@SET ind = MODEL:SET('add basic contact brick', @str varname_u, @str multname_n, @str dataname_r, @spmat BN[, @str dataname_gap[, @str dataname_alpha[, @int symmetrized]])
    
    Add a contact without friction brick to the model. If U is the vector
    of degrees of freedom on which the unilateral contraint is applied,
    the matrix `BN` have to be such that this conctraint is defined by
    $B_N U \le 0$. The constraint is prescribed thank to a multiplier
    `multname_n` whose dimension should be equal to the number of lines of
    `BN`. The augmentation parameter `r` should be chosen in a range of
    acceptabe values (see Getfem user documentation). `dataname_gap` is an
    optional parameter representing the initial gap. It can be a single value
    or a vector of value. `dataname_alpha` is an optional homogenization
    parameter for the augmentation parameter
    (see Getfem user documentation). The parameter `symmetrized` indicates
    that the symmetry of the tangent matrix will be kept or not @*/
    
    std::string varname_u = in.pop().to_string();
    std::string multname_n = in.pop().to_string();
    std::string dataname_r = in.pop().to_string();
    dal::shared_ptr<gsparse> BN = in.pop().to_sparse();
    if (BN->is_complex()) THROW_BADARG("Complex matrix not allowed");
    std::string dataname_gap;
    if (in.remaining()) dataname_gap = in.pop().to_string();
    std::string dataname_alpha;
    if (in.remaining()) dataname_alpha = in.pop().to_string();
    bool symmetrized = false;
    if (in.remaining()) symmetrized = (in.pop().to_integer(0,1)) != 0;

    getfem::model_real_sparse_matrix BBN;
    if (BN->storage()==gsparse::CSCMAT) {
      gmm::resize(BBN, gmm::mat_nrows(BN->real_csc()),
                  gmm::mat_ncols(BN->real_csc()));
      gmm::copy(BN->real_csc(), BBN);
    }
    else if (BN->storage()==gsparse::WSCMAT) {
      gmm::resize(BBN, gmm::mat_nrows(BN->real_wsc()),
                  gmm::mat_ncols(BN->real_wsc()));
      gmm::copy(BN->real_wsc(), BBN);
    }
    else THROW_BADARG("Matrix BN should be a sparse matrix");

    size_type ind
      = getfem::add_basic_contact_brick(md->model(), varname_u, multname_n,
                                        dataname_r, BBN, dataname_gap,
                                        dataname_alpha, symmetrized);

    out.pop().from_integer(int(ind + config::base_index()));
  } else if (check_cmd(cmd, "contact brick set BN", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('contact brick set BN', @int indbrick, @spmat BN)
    Can be used to set the BN matrix of a basic contact/friction brick. @*/
    size_type ind = in.pop().to_integer();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();

    if (B->is_complex())
      THROW_BADARG("BN should be a real matrix");

    if (B->storage()==gsparse::CSCMAT)
      gmm::copy(B->real_csc(), getfem::contact_brick_set_BN(md->model(), ind));
    else if (B->storage()==gsparse::WSCMAT)
      gmm::copy(B->real_wsc(), getfem::contact_brick_set_BN(md->model(), ind));
    else
      THROW_BADARG("BN should be a sparse matrix");
  } else if (check_cmd(cmd, "contact brick set BT", in, out, 2, 2, 0, 0)) {
    /*@SET MODEL:SET('contact brick set BT', @int indbrick, @spmat BT)
    Can be used to set the BT matrix of a basic contact with friction brick. @*/
    size_type ind = in.pop().to_integer();
    dal::shared_ptr<gsparse> B = in.pop().to_sparse();

    if (B->is_complex())
      THROW_BADARG("BT should be a real matrix");

    if (B->storage()==gsparse::CSCMAT)
      gmm::copy(B->real_csc(), getfem::contact_brick_set_BT(md->model(), ind));
    else if (B->storage()==gsparse::WSCMAT)
      gmm::copy(B->real_wsc(), getfem::contact_brick_set_BT(md->model(), ind));
    else
      THROW_BADARG("BT should be a sparse matrix");

  } else if (check_cmd(cmd, "add contact with rigid obstacle brick", in, out, 6, 7, 0, 1)) {
    /*@SET ind = MODEL:SET('add contact with rigid obstacle brick',  @tmim mim, @str varname_u, @str multname_n, @str dataname_r, @int region, @str obstacle, [,  @int symmetrized])
    
    Add a contact without friction condition with a rigid obstacle
    to the model. The condition is applied on the variable `varname_u`
    on the boundary corresponding to `region`. The rigid obstacle should
    be described with the string `obstacle` being a signed distance to
    the obstacle. This string should be an expression where the coordinates
    are 'x', 'y' in 2D and 'x', 'y', 'z' in 3D. For instance, if the rigid
    obstacle corrspond to $z \le 0$, the corresponding signed distance will
    be simply "z". `multname` should be a fixed size variable whose size is
    the number of degrees of freedom on boundary `region`.
    The augmentation parameter `r` should be chosen in a
    range of acceptabe values (close to the Young modulus of the elastic
      body, see Getfem user documentation). The
    parameter `symmetrized` indicates that the symmetry of the tangent
    matrix will be kept or not. Basically, this brick compute the matrix BN
    and the vectors gap and alpha and call the basic contact brick. @*/
    
    getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
    std::string varname_u = in.pop().to_string();
    std::string multname_n = in.pop().to_string();
    std::string dataname_r = in.pop().to_string();
    size_type region = in.pop().to_integer();
    std::string obstacle = in.pop().to_string();
    bool symmetrized = false;
    if (in.remaining()) symmetrized = (in.pop().to_integer(0,1)) != 0;

    size_type ind
      = getfem::add_contact_with_rigid_obstacle_brick
      (md->model(), gfi_mim->mesh_im(), varname_u, multname_n,
       dataname_r, region, obstacle, symmetrized);

    out.pop().from_integer(int(ind + config::base_index()));
  } else bad_cmd(cmd);
}
