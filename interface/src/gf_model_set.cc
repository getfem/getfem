// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2010 Yves Renard.
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
#include <getfemint_misc.h>
#include <getfemint_models.h>
#include <getfemint_mesh_fem.h>
#include <getfemint_workspace.h>
#include <getfemint_mesh_im.h>
#include <getfemint_gsparse.h>
#include <getfem/getfem_Coulomb_friction.h>
#include <getfem/getfem_nonlinear_elasticity.h>
#include <getfem/getfem_plasticity.h>

using namespace getfemint;


/*@GFDOC
  Modifies a model object.
@*/


// Object for the declaration of a new sub-command.

struct sub_gf_md_set : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfemint_model *md) = 0;
};

typedef boost::intrusive_ptr<sub_gf_md_set> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_md_set {                                \
      virtual void run(getfemint::mexargs_in& in,                       \
                       getfemint::mexargs_out& out,                     \
                       getfemint_model *md)                             \
      { dummy_func(in); dummy_func(out); code }                         \
    };                                                                  \
    psub_command psubc = new subc;                                      \
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;         \
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;     \
    subc_tab[cmd_normalize(name)] = psubc;                              \
  }                           


void gf_model_set(getfemint::mexargs_in& m_in,
                  getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@SET ('clear')
      Clear the model.@*/
    sub_command
      ("clear", 0, 0, 0, 1,
       md->clear();
       );


    /*@SET ('add fem variable', @str name, @tmf mf[, @int niter])
      Add a variable to the model linked to a @tmf. `name` is the variable
      name and `niter` is the optional number of version of the data stored,
      for time integration schemes.@*/
    sub_command
      ("add fem variable", 2, 3, 0, 0,
       std::string name = in.pop().to_string();
       getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
       size_type niter = 1;
       if (in.remaining()) niter = in.pop().to_integer(1,10);
       md->model().add_fem_variable(name, gfi_mf->mesh_fem(), niter);
       workspace().set_dependance(md, gfi_mf);
       );


    /*@SET ('add variable', @str name, @int size[, @int niter])
      Add a variable to the model of constant size. `name` is the variable
      name and `niter` is the optional number of version of the data stored,
      for time integration schemes. @*/
    sub_command
      ("add variable", 2, 3, 0, 0,
       std::string name = in.pop().to_string();
       size_type s = in.pop().to_integer();
       size_type niter = 1;
       if (in.remaining()) niter = in.pop().to_integer(1,10);
       md->model().add_fixed_size_variable(name, s, niter);
       );


    /*@SET ('resize variable', @str name, @int size)
      Resize a  constant size variable of the model. `name` is the variable
      name. @*/
    sub_command
      ("resize variable", 2, 2, 0, 0,
       std::string name = in.pop().to_string();
       size_type s = in.pop().to_integer();
       md->model().resize_fixed_size_variable(name, s);
       );


    /*@SET ('add multiplier', @str name, @tmf mf, @str primalname[, @int niter])
    Add a particular variable linked to a fem being a multiplier with
    respect to a primal variable. The dof will be filtered with the
    ``gmm::range_basis`` function applied on the terms of the model
    which link the multiplier and the primal variable. This in order to
    retain only linearly independant constraints on the primal variable.
    Optimized for boundary multipliers. `niter` is the optional number
    of version of the data stored, for time integration schemes. @*/
    sub_command
      ("add multiplier", 3, 4, 0, 0,
       std::string name = in.pop().to_string();
       getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
       std::string primalname = in.pop().to_string();
       size_type niter = 1;
       if (in.remaining()) niter = in.pop().to_integer(1,10);
       md->model().add_multiplier(name, gfi_mf->mesh_fem(), primalname, niter);
       workspace().set_dependance(md, gfi_mf);
       );


    /*@SET ('add fem data', @str name, @tmf mf[, @int qdim[, @int niter]])
      Add a data to the model linked to a @tmf. `name` is the data name,
      `qdim` is the optional dimension of the data over the @tmf and
      `niter` is the optional number of version of the data stored,
      for time integration schemes. @*/
    sub_command
      ("add fem data", 2, 4, 0, 0,
       std::string name = in.pop().to_string();
       getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
       dim_type qdim = 1;
       if (in.remaining()) qdim = dim_type(in.pop().to_integer(1,255));
       size_type niter = 1;
       if (in.remaining()) niter = in.pop().to_integer(1,10);
       md->model().add_fem_data(name, gfi_mf->mesh_fem(), qdim, niter);
       workspace().set_dependance(md, gfi_mf);
       );


    /*@SET ('add initialized fem data', @str name, @tmf mf, @vec V)
      Add a data to the model linked to a @tmf. `name` is the data name.
      The data is initiakized with `V`. The data can be a scalar or vector
      field.@*/
    sub_command
      ("add initialized fem data", 3, 3, 0, 0,
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
       );


    /*@SET ('add data', @str name, @int size[, @int niter])
      Add a data to the model of constant size. `name` is the data name
      and `niter` is the optional number of version of the data stored,
      for time integration schemes. @*/
    sub_command
      ("add data", 2, 3, 0, 0,
       std::string name = in.pop().to_string();
       size_type s = in.pop().to_integer();
       size_type niter = 1;
       if (in.remaining()) niter = in.pop().to_integer(1,10);
       md->model().add_fixed_size_data(name, s, niter);
       );


    /*@SET ('add initialized data', @str name, @vec V)
      Add a fixed size data to the model linked to a @tmf.
      `name` is the data name and `V` is the value of the data.@*/
    sub_command
      ("add initialized data", 2, 2, 0, 0,
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
       );


    /*@SET ('variable', @str name, @vec V[, @int niter])
      Set the value of a variable or data. `name` is the data name
      and `niter` is the optional number of version of the data stored,
      for time integration schemes.@*/
    sub_command
      ("variable", 2, 3, 0, 0,
       std::string name = in.pop().to_string();
       if (!md->is_complex()) {
         darray st = in.pop().to_darray();
         size_type niter = 0;
         if (in.remaining())
           niter = in.pop().to_integer(0,10) - config::base_index();
         GMM_ASSERT1(st.size()==md->model().real_variable(name, niter).size(),
                     "Bad size in assigment");
         md->model().set_real_variable(name, niter).assign(st.begin(), st.end());
       } else {
         carray st = in.pop().to_carray();
         size_type niter = 0;
         if (in.remaining())
           niter = in.pop().to_integer(0,10) - config::base_index();
         GMM_ASSERT1(st.size() == md->model().complex_variable(name,
                                                               niter).size(),
                     "Bad size in assigment");
         md->model().set_complex_variable(name, niter).assign(st.begin(),
                                                              st.end());
       }
       );


    /*@SET ('to variables', @vec V)
      Set the value of the variables of the model with the vector `V`.
      Typically, the vector `V` results of the solve of the tangent
      linear system (usefull to solve your problem with you own solver).@*/
    sub_command
      ("to variables", 1, 1, 0, 0,
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
       );


    /*@SET ind = ('add Laplacian brick', @tmim mim, @str varname[, @int region])
    Add a Laplacian term to the model relatively to the variable `varname`.
    If this is a vector valued variable, the Laplacian term is added
    componentwise. `region` is an optional mesh region on which the term
    is added. If it is not specified, it is added on the whole mesh. Return
    the brick index in the model.@*/
    sub_command
      ("add Laplacian brick", 2, 3, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       std::string varname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_Laplacian_brick(md->model(), gfi_mim->mesh_im(), varname, region)
       + config::base_index();
       workspace().set_dependance(md, gfi_mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add generic elliptic brick', @tmim mim, @str varname, @str dataname[, @int region])
    Add a generic elliptic term to the model relatively to the variable `varname`.
    The shape of the elliptic term depends both on the variable and the data.
    This corresponds to a term
    :math:`-\text{div}(a\nabla u)`
    where :math:`a` is the data and :math:`u` the variable. The data can be a scalar,
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
    sub_command
      ("add generic elliptic brick", 3, 4, 0, 1,
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
       );

    
    /*@SET ind = ('add source term brick', @tmim mim, @str varname, @str dataname[, @int region[, @str directdataname]])
    Add a source term to the model relatively to the variable `varname`.
    The source term is represented by the data `dataname` which could be
    constant or described on a fem. `region` is an optional mesh region
    on which the term is added. An additional optional data `directdataname`
    can be provided. The corresponding data vector will be directly added
    to the right hand side without assembly. Return the brick index in the
    model.@*/
    sub_command
      ("add source term brick", 3, 5, 0, 1,
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
       );


    /*@SET ind = ('add normal source term brick', @tmim mim, @str varname, @str dataname, @int region)
      Add a source term on the variable `varname` on a boundary `region`.
      This region should be a boundary. The source term is represented by the
      data `dataname` which could be constant or described on a fem. A scalar
      product with the outward normal unit vector to the boundary is performed.
      The main aim of this brick is to represent a Neumann condition with a
      vector data without performing the scalar product with the normal as a
      pre-processing. Return the brick index in the model.@*/
    sub_command
      ("add normal source term brick", 4, 4, 0, 1,
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
       );


    /*@SET ind = ('add Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region[, @str dataname])
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
    sub_command
      ("add Dirichlet condition with multipliers", 4, 5, 0, 1,
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
       );


    /*@SET ind = ('add Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region[, @str dataname, @tmf mf_mult])
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`. This region should be a boundary. The Dirichlet
    condition is prescribed with penalization. The penalization coefficient
    is intially `coeff` and will be added to the data of the model.
    `dataname` is the optional right hand side of the Dirichlet condition.
    It could be constant or described on a fem; scalar or vector valued,
    depending on the variable on which the Dirichlet condition is prescribed.
    `mf_mult` is an optional parameter which allows to weaken the
      Dirichlet condition specifying a multiplier space.
    Return the brick index in the model.@*/
    sub_command
      ("add Dirichlet condition with penalization", 4, 6, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       size_type region = in.pop().to_integer();
       std::string dataname;
       if (in.remaining()) dataname = in.pop().to_string();
       const getfem::mesh_fem *mf_mult = 0;
       if (in.remaining()) mf_mult = in.pop().to_const_mesh_fem();
       size_type ind = config::base_index();
       ind += getfem::add_Dirichlet_condition_with_penalization
       (md->model(), gfi_mim->mesh_im(), varname, coeff, region,
	dataname, mf_mult);
       workspace().set_dependance(md, gfi_mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add generalized Dirichlet condition with multipliers', @tmim mim, @str varname, mult_description, @int region, @str dataname, @str Hname)
    Add a Dirichlet condition on the variable `varname` and the mesh
    region `region`.  This version is for vector field.
    It prescribes a condition :math:`Hu = r`
    where `H` is a matrix field. The region should be a boundary. The Dirichlet
    condition is prescribed with a multiplier variable described by
    `mult_description`. If `mult_description` is a string this is assumed
    to be the variable name corresponding to the multiplier (which should be
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
    sub_command
      ("add generalized Dirichlet condition with multipliers", 6, 6, 0, 1,
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
       );


    /*@SET ind = ('add generalized Dirichlet condition with penalization', @tmim mim, @str varname, @scalar coeff, @int region, @str dataname, @str Hname[, @tmf mf_mult])
      Add a Dirichlet condition on the variable `varname` and the mesh
      region `region`. This version is for vector field.
      It prescribes a condition :math:`Hu = r`
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
      `mf_mult` is an optional parameter which allows to weaken the
      Dirichlet condition specifying a multiplier space.
      Return the brick index in the model.@*/
    sub_command
      ("add generalized Dirichlet condition with penalization", 6, 7, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       std::string varname = in.pop().to_string();
       double coeff = in.pop().to_scalar();
       size_type region = in.pop().to_integer();
       std::string dataname= in.pop().to_string();
       std::string Hname= in.pop().to_string();
       const getfem::mesh_fem *mf_mult = 0;
       if (in.remaining()) mf_mult = in.pop().to_const_mesh_fem();
       size_type ind = config::base_index();
       ind += getfem::add_generalized_Dirichlet_condition_with_penalization
       (md->model(), gfi_mim->mesh_im(), varname, coeff, region,
        dataname, Hname, mf_mult);
       workspace().set_dependance(md, gfi_mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ('change penalization coeff', @int ind_brick, @scalar coeff)
    Change the penalization coefficient of a Dirichlet condition with
    penalization brick. If the brick is not of this kind, this
    function has an undefined behavior.@*/
    sub_command
      ("change penalization coeff", 2, 2, 0, 0,
       size_type ind_brick = in.pop().to_integer();
       double coeff = in.pop().to_scalar();
       getfem::change_penalization_coeff(md->model(), ind_brick, coeff);
       );


    /*@SET ind = ('add Helmholtz brick', @tmim mim, @str varname, @str dataname[, @int region])
      Add a Helmholtz term to the model relatively to the variable `varname`.
      `dataname` should contain the wave number. `region` is an optional mesh
      region on which the term is added. If it is not specified, it is added
      on the whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add Helmholtz brick", 3, 4, 0, 1,
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
       );


    /*@SET ind = ('add Fourier Robin brick', @tmim mim, @str varname, @str dataname, @int region)
    Add a Fourier-Robin term to the model relatively to the variable
    `varname`. This corresponds to a weak term of the form
    :math:`\int (qu).v`. `dataname`
    should contain the parameter :math:`q` of
    the Fourier-Robin condition. `region` is the mesh region on which
    the term is added. Return the brick index in the model.@*/
    sub_command
      ("add Fourier Robin brick", 4, 4, 0, 1,
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
       );


    /*@SET ind = ('add constraint with multipliers', @str varname, @str multname, @tspmat B, @vec L)
    Add an additional explicit constraint on the variable `varname` thank to
    a multiplier `multname` peviously added to the model (should be a fixed
    size variable). The constraint is :math:`BU=L`
    with `B` being a rectangular sparse matrix. It is possible to change
    the constraint at any time whith the methods MODEL:SET('set private matrix')
    and MODEL:SET('set private rhs'). Return the brick index in the model.@*/
    sub_command
      ("add constraint with multipliers", 4, 4, 0, 1,
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
       );


    /*@SET ind = ('add constraint with penalization', @str varname, @scalar coeff, @tspmat B, @vec L)
    Add an additional explicit penalized constraint on the variable `varname`.
    The constraint is :math`BU=L` with `B` being a rectangular sparse matrix.
    Be aware that `B` should not contain a palin row, otherwise the whole
    tangent matrix will be plain. It is possible to change the constraint
    at any time whith the methods MODEL:SET('set private matrix')
    and MODEL:SET('set private rhs'). The method
    MODEL:SET('change penalization coeff') can be used. Return the brick
    index in the model.@*/
    sub_command
      ("add constraint with penalization", 4, 4, 0, 1,
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
       );


    /*@SET ind = ('add explicit matrix', @str varname1, @str varname2, @tspmat B[, @int issymmetric[, @int iscoercive]])
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
    sub_command
      ("add explicit matrix", 3, 5, 0, 1,
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
       );


    /*@SET ind = ('add explicit rhs', @str varname, @vec L)
      Add a brick representing an explicit right hand side to be added to
      the right hand side of the tangent linear system relatively to the
      variable 'varname'. The given rhs should have the same size than the
      dimension of 'varname'. The rhs can be changed by the command
      MODEL:SET('set private rhs'). Return the brick index in the model.@*/
    sub_command
      ("add explicit rhs", 2, 2, 0, 1,
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
       );


    /*@SET ('set private matrix', @int indbrick, @tspmat B)
      For some specific bricks having an internal sparse matrix
      (explicit bricks: 'constraint brick' and 'explicit matrix brick'),
      set this matrix. @*/
    sub_command
      ("set private matrix", 2, 2, 0, 0,
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
       );


    /*@SET ('set private rhs', @int indbrick, @vec B)
      For some specific bricks having an internal right hand side vector
      (explicit bricks: 'constraint brick' and 'explicit rhs brick'),
      set this rhs. @*/
    sub_command
      ("set private rhs", 2, 2, 0, 0,
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
       );


    /*@SET ind = ('add isotropic linearized elasticity brick', @tmim mim, @str varname, @str dataname_lambda, @str dataname_mu[, @int region])
      Add an isotropic linearized elasticity term to the model relatively to
      the variable `varname`. `dataname_lambda` and `dataname_mu` should
      contain the Lam\'e coefficients. `region` is an optional mesh region
      on which the term is added. If it is not specified, it is added
      on the whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add isotropic linearized elasticity brick", 4, 5, 0, 1,
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
       );


    /*@SET ind = ('add linear incompressibility brick', @tmim mim, @str varname, @str multname_pressure[, @int region[, @str dataname_coeff]])
    Add an linear incompressibility condition on `variable`. `multname_pressure`
    is a variable which represent the pressure. Be aware that an inf-sup
    condition between the finite element method describing the pressure and the
    primal variable has to be satisfied. `region` is an optional mesh region on
    which the term is added. If it is not specified, it is added on the whole mesh.
    `dataname_coeff` is an optional penalization coefficient for nearly
    incompressible elasticity for instance. In this case, it is the inverse
    of the Lam\'e coefficient :math:`\lambda`. Return the brick index in the model.@*/
    sub_command
      ("add linear incompressibility brick", 3, 5, 0, 1,
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
       );


    /*@SET ind = ('add nonlinear elasticity brick', @tmim mim, @str varname, @str constitutive_law, @str dataname[, @int region])
    Add a nonlinear elasticity term to the model relatively to the
    variable `varname`. `lawname` is the constitutive law which
    could be 'SaintVenant Kirchhoff', 'Mooney Rivlin' or 'Ciarlet Geymonat'.
    `dataname` is a vector of parameters for the constitutive law. Its length
    depends on the law. It could be a short vector of constant values or a
    vector field described on a finite element method for variable
    coefficients. `region` is an optional mesh region on which the term
    is added. If it is not specified, it is added on the whole mesh. Return the
    brick index in the model.@*/
    sub_command
      ("add nonlinear elasticity brick", 4, 5, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       std::string varname = in.pop().to_string();
       std::string lawname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_nonlinear_elasticity_brick
       (md->model(), gfi_mim->mesh_im(), varname,
        abstract_hyperelastic_law_from_name(lawname), dataname, region);
       
       workspace().set_dependance(md, gfi_mim);
       out.pop().from_integer(int(ind));
       );

    /*@SET ind = ('add elastoplasticity brick', @tmim mim ,@str projname, @str varname, @str datalambda, @str datamu, @str datathreshold, @str datasigma[, @int region])
      Add a nonlinear elastoplastic term to the model relatively to the
      variable `varname`, in small deformations, for an isotropic material
      and for a quasistatic model. 'projname' is the type of projection that
      we want to use. For the moment, only the Von Mises projection is
      computing that we could entering 'VM' or 'Von Mises'.
      `datasigma` is the variable representing the constraints on the material.
      Be carefull that `varname` and `datasigma` are composed of two iterates
      for the time scheme needed for the Newton algorithm used.
      Moreover, the finite element method on which `varname` is described
      is an K ordered mesh_fem, the `datasigma` one have to be at least
      an K-1 ordered mesh_fem. 
      `datalambda` and `datamu` are the Lamé coefficients of the studied
      material.
      `datathreshold` is the plasticity threshold of the material.
      The three last variable could be constants or described on the
      same finite element method.
      `region` is an optional mesh region on which the term is added.
      If it is not specified, it is added on the whole mesh.
      Return the brick index in the model.@*/
    sub_command
      ("add elastoplasticity brick", 7, 8, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       std::string projname = in.pop().to_string();
       std::string varname = in.pop().to_string();
       std::string datalambda = in.pop().to_string();
       std::string datamu = in.pop().to_string();
       std::string datathreshold = in.pop().to_string();
       std::string datasigma = in.pop().to_string();
       size_type region = size_type(-1);

       // getfem::VM_projection proj(0);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind = config::base_index() +
       add_elastoplasticity_brick
       (md->model(), gfi_mim->mesh_im(),
        abstract_constraints_projection_from_name(projname), 
	varname, datalambda, datamu,
	datathreshold, datasigma,
	region);
       
       workspace().set_dependance(md, gfi_mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add nonlinear incompressibility brick', @tmim mim, @str varname, @str multname_pressure[, @int region])
    Add an nonlinear incompressibility condition on `variable` (for large
    strain elasticity). `multname_pressure`
    is a variable which represent the pressure. Be aware that an inf-sup
    condition between the finite element method describing the pressure and the
    primal variable has to be satisfied. `region` is an optional mesh region on
    which the term is added. If it is not specified, it is added on the
    whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add nonlinear incompressibility brick", 3, 4, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       std::string varname = in.pop().to_string();
       std::string multname = in.pop().to_string();
       size_type region = size_type(-1);
       if (in.remaining()) region = in.pop().to_integer();
       size_type ind
       = getfem::add_nonlinear_incompressibility_brick
       (md->model(), gfi_mim->mesh_im(), varname, multname, region)
       + config::base_index();
       workspace().set_dependance(md, gfi_mim);
       out.pop().from_integer(int(ind));
       );


    /*@SET ind = ('add mass brick', @tmim mim, @str varname[, @str dataname_rho[, @int region]])
      Add mass term to the model relatively to the variable `varname`.
      If specified, the data `dataname_rho` should contain the
      density (1 if omitted). `region` is an optional mesh region on
      which the term is added. If it is not specified, it
      is added on the whole mesh. Return the brick index in the model.@*/
    sub_command
      ("add mass brick", 2, 4, 0, 1,
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
       );


    /*@SET ind = ('add basic d on dt brick', @tmim mim, @str varnameU, @str dataname_dt[, @str dataname_rho[, @int region]])
    Add the standard discretization of a first order time derivative on
    `varnameU`. The parameter `dataname_rho` is the density which could
    be omitted (the defaul value is 1). This brick should be used in
    addition to a time dispatcher for the other terms. Return the brick
    index in the model.@*/
    sub_command
      ("add basic d on dt brick", 3, 5, 0, 1,
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
       );


    /*@SET ind = ('add basic d2 on dt2 brick', @tmim mim, @str varnameU,  @str datanameV, @str dataname_dt, @str dataname_alpha,[, @str dataname_rho[, @int region]])
    Add the standard discretization of a second order time derivative
    on `varnameU`. `datanameV` is a data represented on the same finite
    element method as U which represents the time derivative of U. The
    parameter `dataname_rho` is the density which could be omitted (the defaul
    value is 1). This brick should be used in addition to a time dispatcher for
    the other terms. The time derivative :math:`v` of the
    variable :math:`u` is preferably computed as a
    post-traitement which depends on each scheme. The parameter `dataname_alpha`
    depends on the time integration scheme. Return the brick index in the model.@*/
    sub_command
      ("add basic d2 on dt2 brick", 5, 7, 0, 1,
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
       );


    /*@SET ('add theta method dispatcher', @ivec bricks_indices, @str theta)
      Add a theta-method time dispatcher to a list of bricks. For instance,
      a matrix term :math:`K` will be replaced by
      :math:`\theta K U^{n+1} + (1-\theta) K U^{n}`.
      @*/
    sub_command
      ("add theta method dispatcher", 2, 2, 0,0,
       dal::bit_vector bv = in.pop().to_bit_vector();
       std::string datanametheta = in.pop().to_string();
       getfem::add_theta_method_dispatcher(md->model(), bv, datanametheta);
       );

    /*@SET ('add midpoint dispatcher', @ivec bricks_indices)
      Add a midpoint time dispatcher to a list of bricks. For instance, a
      nonlinear term :math:`K(U)` will be replaced by
      :math:`K((U^{n+1} +  U^{n})/2)`.@*/
    sub_command
      ("add midpoint dispatcher", 1, 1, 0,0,
       dal::bit_vector bv = in.pop().to_bit_vector();
       getfem::add_midpoint_dispatcher(md->model(), bv);
       );


    /*@SET ('velocity update for order two theta method', @str varnameU,  @str datanameV, @str dataname_dt, @str dataname_theta)
      Function which udpate the velocity :math:`v^{n+1}` after
      the computation of the displacement :math:`u^{n+1}` and
      before the next iteration. Specific for theta-method and when the velocity is
      included in the data of the model. @*/
    sub_command
      ("velocity update for order two theta method", 4, 4, 0,0,
       std::string varnameU = in.pop().to_string();
       std::string varnameV = in.pop().to_string();
       std::string varnamedt = in.pop().to_string();
       std::string varnametheta = in.pop().to_string();
       velocity_update_for_order_two_theta_method
       (md->model(), varnameU, varnameV, varnamedt, varnametheta);
       );


    /*@SET ('velocity update for Newmark scheme', @int id2dt2_brick, @str varnameU,  @str datanameV, @str dataname_dt, @str dataname_twobeta, @str dataname_alpha)
      Function which udpate the velocity
      :math:`v^{n+1}` after
      the computation of the displacement
      :math:`u^{n+1}` and
      before the next iteration. Specific for Newmark scheme
      and when the velocity is
      included in the data of the model.*
      This version inverts the mass matrix by a
      conjugate gradient.@*/
     sub_command
       ("velocity update for Newmark scheme", 6, 6, 0,0,
        size_type id2dt2 = in.pop().to_integer();
        std::string varnameU = in.pop().to_string();
        std::string varnameV = in.pop().to_string();
        std::string varnamedt = in.pop().to_string();
        std::string varnametwobeta = in.pop().to_string();
        std::string varnamegamma = in.pop().to_string();
        velocity_update_for_Newmark_scheme
        (md->model(), id2dt2, varnameU, varnameV, varnamedt,
         varnametwobeta, varnamegamma);
        );


     /*@SET ('disable bricks', @ivec bricks_indices)
       Disable a brick (the brick will no longer participate to the
       building of the tangent linear system).@*/
     sub_command
       ("disable bricks", 1, 1, 0, 0,
        dal::bit_vector bv = in.pop().to_bit_vector();
        for (dal::bv_visitor ii(bv); !ii.finished(); ++ii)
          md->model().disable_brick(ii);
        );


    /*@SET ('unable bricks', @ivec bricks_indices)
    Unable a disabled brick.@*/
     sub_command
       ("unable bricks", 1, 1, 0, 0,
        dal::bit_vector bv = in.pop().to_bit_vector();
        for (dal::bv_visitor ii(bv); !ii.finished(); ++ii)
          md->model().unable_brick(ii);
        );


    /*@SET ('first iter')
    To be executed before the first iteration of a time integration scheme.@*/
     sub_command
       ("first iter", 0, 0, 0, 0,
        md->model().first_iter();
        );


    /*@SET ('next iter')
      To be executed at the end of each iteration of a time
      integration scheme.@*/
     sub_command
       ("next iter", 0, 0, 0, 0,
        md->model().next_iter();
        );


    /*@SET ind = ('add basic contact brick', @str varname_u, @str multname_n[, @str multname_t], @str dataname_r, @tspmat BN[, @tspmat BT, @str dataname_friction_coeff][, @str dataname_gap[, @str dataname_alpha[, @int symmetrized]])
    
    Add a contact with  or without friction brick to the model.
    If U is the vector
    of degrees of freedom on which the unilateral constraint is applied,
    the matrix `BN` have to be such that this constraint is defined by
    :math:`B_N U \le 0`. A friction condition can be considered by adding the three
    parameters `multname_t`, `BT` and `dataname_friction_coeff`. In this case,
    the tangential displacement is :math:`B_T U` and the matrix `BT` should have as
    many rows as `BN` multiplied by :math:`d-1` where :math:`d` is the domain dimension.
    In this case also, `dataname_friction_coeff` is a data which represents the
    coefficient of friction. It can be a scalar or a vector representing a
    value on each contact condition.  The unilateral constraint is prescribed
    thank to a multiplier
    `multname_n` whose dimension should be equal to the number of rows of
    `BN`. If a friction condition is added, it is prescribed with a
    multiplier `multname_t` whose dimension should be equal to the number
    of rows of `BT`. The augmentation parameter `r` should be chosen in
    a range of
    acceptabe values (see Getfem user documentation). `dataname_gap` is an
    optional parameter representing the initial gap. It can be a single value
    or a vector of value. `dataname_alpha` is an optional homogenization
    parameter for the augmentation parameter
    (see Getfem user documentation). The parameter `symmetrized` indicates
    that the symmetry of the tangent matrix will be kept or not (except for
    the part representing the coupling between contact and friction which
    cannot be symmetrized). @*/
     sub_command
       ("add basic contact brick", 4, 10, 0, 1,
        
        bool friction = false;
        
        std::string varname_u = in.pop().to_string();
        std::string multname_n = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();
        std::string multname_t; std::string friction_coeff;
        
        mexarg_in argin = in.pop();
        if (argin.is_string()) {
          friction = true;
          multname_t = dataname_r;
          dataname_r = argin.to_string();
          argin = in.pop();
        }
        
        dal::shared_ptr<gsparse> BN = argin.to_sparse();
        if (BN->is_complex()) THROW_BADARG("Complex matrix not allowed");
        dal::shared_ptr<gsparse> BT;
        if (friction) {
          BT = in.pop().to_sparse();
          if (BT->is_complex()) THROW_BADARG("Complex matrix not allowed");
          friction_coeff = in.pop().to_string();
        }
        
        std::string dataname_gap;
        dataname_gap = in.pop().to_string();
        std::string dataname_alpha;
        if (in.remaining()) dataname_alpha = in.pop().to_string();
        bool symmetrized = false;
        if (in.remaining()) symmetrized = (in.pop().to_integer(0,1)) != 0;
        
        getfem::CONTACT_B_MATRIX BBN;
        getfem::CONTACT_B_MATRIX BBT;
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
        
        if (friction) {
          if (BT->storage()==gsparse::CSCMAT) {
            gmm::resize(BBT, gmm::mat_nrows(BT->real_csc()),
                        gmm::mat_ncols(BT->real_csc()));
            gmm::copy(BT->real_csc(), BBT);
          }
          else if (BT->storage()==gsparse::WSCMAT) {
            gmm::resize(BBT, gmm::mat_nrows(BT->real_wsc()),
                        gmm::mat_ncols(BT->real_wsc()));
            gmm::copy(BT->real_wsc(), BBT);
          }
          else THROW_BADARG("Matrix BT should be a sparse matrix");
        }
        
        size_type ind;
        if (friction) {
          ind = getfem::add_basic_contact_with_friction_brick
            (md->model(), varname_u, multname_n, multname_t, dataname_r, BBN, BBT,
             friction_coeff, dataname_gap, dataname_alpha, symmetrized);
        } else {
          ind = getfem::add_basic_contact_brick
            (md->model(), varname_u, multname_n, dataname_r, BBN, dataname_gap,
             dataname_alpha, symmetrized);
        }
        
        out.pop().from_integer(int(ind + config::base_index()));
        );


    /*@SET ('contact brick set BN', @int indbrick, @tspmat BN)
    Can be used to set the BN matrix of a basic contact/friction brick. @*/
     sub_command
       ("contact brick set BN", 2, 2, 0, 0,
        size_type ind = in.pop().to_integer();
        dal::shared_ptr<gsparse> B = in.pop().to_sparse();
        
        if (B->is_complex())
          THROW_BADARG("BN should be a real matrix");
        
        if (B->storage()==gsparse::CSCMAT)
          gmm::copy(B->real_csc(),
                    getfem::contact_brick_set_BN(md->model(), ind));
        else if (B->storage()==gsparse::WSCMAT)
          gmm::copy(B->real_wsc(),
                    getfem::contact_brick_set_BN(md->model(), ind));
        else
          THROW_BADARG("BN should be a sparse matrix");
        );


    /*@SET ('contact brick set BT', @int indbrick, @tspmat BT)
      Can be used to set the BT matrix of a basic contact with
      friction brick. @*/
     sub_command
       ("contact brick set BT", 2, 2, 0, 0,
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
        );


    /*@SET ind = ('add contact with rigid obstacle brick',  @tmim mim, @str varname_u, @str multname_n[, @str multname_t], @str dataname_r[, @str dataname_friction_coeff], @int region, @str obstacle[,  @int symmetrized])
    
    Add a contact with or without friction condition with a rigid obstacle
    to the model. The condition is applied on the variable `varname_u`
    on the boundary corresponding to `region`. The rigid obstacle should
    be described with the string `obstacle` being a signed distance to
    the obstacle. This string should be an expression where the coordinates
    are 'x', 'y' in 2D and 'x', 'y', 'z' in 3D. For instance, if the rigid
    obstacle correspond to :math:`z \le 0`, the corresponding signed distance will
    be simply "z". `multname_n` should be a fixed size variable whose size is
    the number of degrees of freedom on boundary `region`. It represent the
    contact equivalent nodal forces. In order to add a friction condition
    one has to add the `multname_t` and `dataname_friction_coeff` parameters.
    `multname_t` should be a fixed size variable whose size is
    the number of degrees of freedom on boundary `region` multiplied by :math:`d-1`
    where :math:`d` is the domain dimension. It represent the friction equivalent
    nodal forces.
    The augmentation parameter `r` should be chosen in a
    range of acceptabe values (close to the Young modulus of the elastic
    body, see Getfem user documentation).  `dataname_friction_coeff` is
    the friction coefficient. It could be a scalar or a vector of values
    representing the friction coefficient on each contact node. The
    parameter `symmetrized` indicates that the symmetry of the tangent
    matrix will be kept or not. Basically, this brick compute the matrix BN
    and the vectors gap and alpha and calls the basic contact brick. @*/
     sub_command
       ("add contact with rigid obstacle brick", 6, 9, 0, 1,
        
        bool friction = false;
        
        getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
        std::string varname_u = in.pop().to_string();
        std::string multname_n = in.pop().to_string();
        std::string dataname_r = in.pop().to_string();
        std::string multname_t;  std::string dataname_fr;
        
        mexarg_in argin = in.pop();
        if (argin.is_string()) {
          friction = true;
          multname_t = dataname_r;
          dataname_r = argin.to_string();
          dataname_fr = in.pop().to_string();
          argin = in.pop();
        }
        
        size_type region = argin.to_integer();
        std::string obstacle = in.pop().to_string();
        bool symmetrized = false;
        if (in.remaining()) symmetrized = (in.pop().to_integer(0,1)) != 0;
        
        size_type ind;
        
        if (friction)
          ind = getfem::add_contact_with_friction_with_rigid_obstacle_brick
            (md->model(), gfi_mim->mesh_im(), varname_u, multname_n,
             multname_t, dataname_r, dataname_fr, region, obstacle,
             symmetrized);
        else
          ind = getfem::add_contact_with_rigid_obstacle_brick
            (md->model(), gfi_mim->mesh_im(), varname_u, multname_n,
             dataname_r, region, obstacle, symmetrized);
        workspace().set_dependance(md, gfi_mim);
        out.pop().from_integer(int(ind + config::base_index()));
        );


    /*@SET ind = ('add nonmatching meshes contact brick',  @tmim mim1[, @tmim mim2], @str varname_u1[, @str varname_u2], @str multname_n[, @str multname_t], @str dataname_r[, @str dataname_fr], @int rg1, @int rg2[, @int slave1, @int slave2,  @int symmetrized])
    
    Add a contact with or without friction condition between two faces of
    one or two elastic bodies. The condition is applied on the variable
    `varname_u1` or the variables `varname_u1` and `varname_u2` depending
    if a single or two distinct displacement fields are given. Integers
    `rg1` and `rg2` represent the regions expected to come in contact with
    each other. In the single displacement variable case the regions defined
    in both `rg1` and `rg2` refer to the variable `varname_u1`. In the case
    of two displacement variables, `rg1` refers to `varname_u1` and `rg2`
    refers to `varname_u2`. `multname_n` should be a fixed size variable
    whose size is the number of degrees of freedom on those regions among
    the ones defined in `rg1` and `rg2` which are characterized as "slaves".
    It represents the contact equivalent nodal normal forces. `multname_t`
    should be a fixed size variable whose size corresponds to the size of
    `multname_n` multiplied by qdim - 1 . It represents the contact
    equivalent nodal tangent (frictional) forces. The augmentation parameter
    `r` should be chosen in a range of acceptabe values (close to the Young
    modulus of the elastic body, see Getfem user documentation). The
    friction coefficient stored in the parameter `fr` is either a single
    value or a vector of the same size as `multname_n`. The optional
    parameters `slave1` and `slave2` declare if the regions defined in `rg1`
    and `rg2` are correspondingly considered as "slaves". By default
    `slave1` is true and `slave2` is false, i.e. `rg1` contains the slave
    surfaces, while 'rg2' the master surfaces. Preferrably only one of
    `slave1` and `slave2` is set to true. The parameter `symmetrized`
    indicates that the symmetry of the tangent matrix will be kept or not.
    Basically, this brick computes the matrices BN and BT and the vectors
    gap and alpha and calls the basic contact brick. @*/
     sub_command
       ("add nonmatching meshes contact brick", 6, 13, 0, 1,

        bool two_variables = true;
        bool friction = false;

        getfemint_mesh_im *gfi_mim1;
        getfemint_mesh_im *gfi_mim2;
        std::string varname_u1;
        std::string varname_u2;
        bool slave1=true; bool slave2=false; bool symmetrized=false;

        gfi_mim1 = in.pop().to_getfemint_mesh_im();
        mexarg_in argin = in.pop();
        if (argin.is_string()) {
          two_variables = false;
          gfi_mim2 = gfi_mim1;
          varname_u1 = argin.to_string();
          varname_u2 = varname_u1;
        } else {
          gfi_mim2 = argin.to_getfemint_mesh_im();
          varname_u1 = in.pop().to_string();
          varname_u2 = in.pop().to_string();
        }
        std::string multname_n = in.pop().to_string();
        std::string multname_t;
        std::string dataname_r = in.pop().to_string();
        std::string dataname_fr;
        argin = in.pop();
        if (!argin.is_integer()) {
          friction = true;
          multname_t = dataname_r;
          dataname_r = in.pop().to_string();
          dataname_fr = in.pop().to_string();
          argin = in.pop();
        }
        std::vector<size_type> vrg1(1,argin.to_integer());
        std::vector<size_type> vrg2(1,in.pop().to_integer());
        if (in.remaining()) slave1 = (in.pop().to_integer(0,1)) != 0;
        if (in.remaining()) slave2 = (in.pop().to_integer(0,1)) != 0;
        if (in.remaining()) symmetrized = (in.pop().to_integer(0,1)) != 0;

        size_type ind;
        if (!friction)
          ind = getfem::add_nonmatching_meshes_contact_brick
            (md->model(), gfi_mim1->mesh_im(), gfi_mim2->mesh_im(),
             varname_u1, varname_u2, multname_n, dataname_r,
             vrg1, vrg2, slave1, slave2, symmetrized);
        else
          ind = getfem::add_nonmatching_meshes_contact_with_friction_brick
            (md->model(), gfi_mim1->mesh_im(), gfi_mim2->mesh_im(),
             varname_u1, varname_u2, multname_n, multname_t,
             dataname_r, dataname_fr,
             vrg1, vrg2, slave1, slave2, symmetrized);
        workspace().set_dependance(md, gfi_mim1);
        if (two_variables)
          workspace().set_dependance(md, gfi_mim2);
        out.pop().from_integer(int(ind + config::base_index()));
        );

  }
 
  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfemint_model *md  = m_in.pop().to_getfemint_model(true);
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, md);
  }
  else bad_cmd(init_cmd);

}
