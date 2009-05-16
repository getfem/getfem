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

/**\file gf_model_get.cc
   \brief getfemint_model getter.
*/

#include <getfemint.h>
#include <getfemint_models.h>
#include <getfem/getfem_model_solvers.h>
#include <getfemint_mdbrick.h>
#include <getfemint_mesh_fem.h>

using namespace getfemint;


#define RETURN_SPARSE(realmeth, cplxmeth)				\
  if (!md->is_complex()) {						\
    gf_real_sparse_by_col M(gmm::mat_nrows(md->model().realmeth()),	\
			    gmm::mat_ncols(md->model().realmeth()));	\
    gmm::copy(md->model().realmeth(), M);				\
    out.pop().from_sparse(M);						\
  } else {								\
    gf_cplx_sparse_by_col M(gmm::mat_nrows(md->model().cplxmeth()),	\
			    gmm::mat_ncols(md->model().cplxmeth()));	\
    gmm::copy(md->model().cplxmeth(), M);				\
    out.pop().from_sparse(M);						\
  }

#define RETURN_VECTOR(realmeth, cplxmeth)			\
  if (!md->is_complex()) {					\
    out.pop().from_dcvector(md->model().realmeth);		\
  } else {							\
    out.pop().from_dcvector(md->model().cplxmeth);		\
  }

/*MLABCOM

  FUNCTION M = gf_model_get(cmd, [, args])
  Get information from a model object.

  @RDATTR MODEL:GET('is_complex')
  @GET MODEL:GET('tangent_matrix')
  @GET MODEL:GET('rhs')
  @GET MODEL:GET('memsize')
  @GET MODEL:GET('listvar')
  @GET MODEL:GET('listbricks')
  @GET MODEL:GET('variable')
  @GET MODEL:GET('mult varname Dirichlet')
  @GET MODEL:GET('from variables')
  @GET MODEL:GET('solve')
  @GET MODEL:GET('compute isotropic linearized Von Mises or Tresca')

MLABCOM*/

void gf_model_get(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {

  if (in.narg() < 2) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  getfemint_model *md  = in.pop().to_getfemint_model();
  std::string cmd        = in.pop().to_string();
  if (check_cmd(cmd, "is_complex", in, out, 0, 0, 0, 1)) {
    /*@RDATTR b = MODEL:GET('is_complex')
    Return 0 is the model is real, 1 if it is complex.@*/
    out.pop().from_integer(md->is_complex());
  } else if (check_cmd(cmd, "tangent_matrix", in, out, 0, 0, 0, 1)) {
    /*@GET T = MODEL:GET('tangent_matrix')
    Return the tangent matrix stored in the model .@*/
    RETURN_SPARSE(real_tangent_matrix, complex_tangent_matrix);
  } else if (check_cmd(cmd, "rhs", in, out, 0, 0, 0, 1)) {
    /*@GET MODEL:GET('rhs')
    Return the right hand side of the tangent problem.@*/
    RETURN_VECTOR(real_rhs(), complex_rhs());
  } else if (check_cmd(cmd, "memsize", in, out, 0, 0, 0, 1)) {
    /*@GET z = MODEL:GET('memsize')
    Return a rough approximation of the amount of memory (in bytes) used by
    the model.@*/
    out.pop().from_integer(int(md->memsize()));
  } else if (check_cmd(cmd, "listvar", in, out, 0, 0, 0, 0)) {
    /*@GET MODEL:GET('listvar')
    print to the output the list of variables and constants of the model.@*/
    md->model().listvar(infomsg());
  } else if (check_cmd(cmd, "listbricks", in, out, 0, 0, 0, 0)) {
    /*@GET MODEL:GET('listbricks')
    print to the output the list of bricks of the model.@*/
    md->model().listbricks(infomsg());
  } else if (check_cmd(cmd, "variable", in, out, 1, 2, 0, 1)) {
    /*@GET V = MODEL:GET('variable', @str name[, @int niter])
    Gives the value of a variable or data.@*/
    std::string name = in.pop().to_string();
    size_type niter = 0;
    if (in.remaining())
      niter = in.pop().to_integer(0,10) - config::base_index();
    RETURN_VECTOR(real_variable(name, niter), complex_variable(name, niter));
  } else if (check_cmd(cmd, "mult varname Dirichlet", in, out, 1, 1, 0, 1)) {
    /*@GET name = MODEL:GET('mult varname Dirichlet', @int ind_brick)
    Gives the name of the multiplier variable for a Dirichlet brick.
    If the brick is not a Dirichlet condition with multiplier brick,
    this function has an undefined behavior@*/
    size_type ind_brick = in.pop().to_integer();
    out.pop().from_string(getfem::mult_varname_Dirichlet(md->model(),
							 ind_brick).c_str());
  } else if (check_cmd(cmd, "from variables", in, out, 0, 0, 0, 1)) {
    /*@GET V = MODEL:GET('from variables')
    Return the vector of all the degrees of freedom of the model consisting
    of the concatenation of the variables of the model (usefull
    to solve your problem with you own solver). @*/
    if (!md->is_complex()) {
      std::vector<double> V(md->model().nb_dof());
      md->model().from_variables(V);
      out.pop().from_dcvector(V);
    } else {
      std::vector<std::complex<double> > V(md->model().nb_dof());
      md->model().from_variables(V);
      out.pop().from_dcvector(V);
    }
  } else if (check_cmd(cmd, "assembly", in, out, 0, 1, 0, 0)) {
    /*@GET MODEL:GET('assembly'[, @str option])
    Assembly of the tangent system taking into account the terms
    from all bricks. `option`, if specified, should be 'build all',
    'build rhs' or 'build matrix'. The default is to build the whole
    tangent linear system (matrix and rhs). This function is usefull to solve
    your problem with you own solver. @*/
    std::string option = "build all";
    if (in.remaining()) option = in.pop().to_string();
    getfem::model::assembly_version version = getfem::model::BUILD_ALL;
    if (cmd_strmatch(option, "build all"))
      version = getfem::model::BUILD_ALL;
    else if (cmd_strmatch(option, "build rhs"))
      version = getfem::model::BUILD_RHS;
    else if (cmd_strmatch(option, "build matrix"))
      version = getfem::model::BUILD_MATRIX;
    else THROW_BADARG("bad option: " << option);
    md->model().assembly(version); 
  } else if (check_cmd(cmd, "solve", in, out, 0, -1, 0, 0)) {
    /*@GET MODEL:GET('solve'[,...])
    Run the standard getfem solver.

    Note that you should be able to use your own solver if you want
    (it is possible to obtain the tangent matrix and its right hand
    side with the MODEL:GET('tangent matrix') etc.).<Par>

    Various options can be specified:<Par>

    - 'noisy' or 'very noisy'<par>
       the solver will display some information showing the progress<par>
       (residual values etc.).<par>
    - 'max_iter', NIT<par>
       set the maximum iterations numbers.<par>
    - 'max_res', RES<par>
       set the target residual value.<par>
    - 'lsolver', SOLVERNAME<par>
       select explicitely the solver used for the linear systems (the<par>
       default value is 'auto', which lets getfem choose itself).<par>
       Possible values are 'superlu', 'mumps' (if supported),<par>
       'cg/ildlt', 'gmres/ilu' and 'gmres/ilut'.@*/

    getfemint::interruptible_iteration iter;
    std::string lsolver = "auto";
    while (in.remaining() && in.front().is_string()) {
      std::string opt = in.pop().to_string();
      if (cmd_strmatch(opt, "noisy")) iter.set_noisy(1);
      else if (cmd_strmatch(opt, "very noisy")) iter.set_noisy(2);
      else if (cmd_strmatch(opt, "max_iter")) {
	if (in.remaining()) iter.set_maxiter(in.pop().to_integer());
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "max_res")) {
	if (in.remaining()) iter.set_resmax(in.pop().to_scalar());
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "lsolver")) {
	if (in.remaining()) lsolver = in.pop().to_string();
	else THROW_BADARG("missing solver name for " << opt);
      } else THROW_BADARG("bad option: " << opt);
    }
    gmm::default_newton_line_search ls(size_t(-1), 5.0/3.0,
				       1.0/1000.0, 3.0/5.0, 1.6);

    if (!md->model().is_complex()) {
      getfem::standard_solve(md->model(), iter,
			     getfem::rselect_linear_solver(md->model(),
							   lsolver), ls);
    } else {
      getfem::standard_solve(md->model(), iter,
			     getfem::cselect_linear_solver(md->model(),
							   lsolver), ls);
    }
  } else if (check_cmd(cmd, "compute isotropic linearized Von Mises or Tresca", in, out, 4, 5, 0, 1)) {
    /*@GET V = MODEL:GET('compute isotropic linearized Von Mises or Tresca', @str varname, @str dataname_lambda, @str dataname_mu, @tmf mf_vm[, @str version])
      Compute the Von-Mises stress or the Tresca stress of a field
      (only valid for isotropic linearized elasticity in 3D).
      `version` should be  'Von Mises' or 'Tresca'
      ('Von Mises' is the default). @*/
    std::string varname = in.pop().to_string();
    std::string dataname_lambda = in.pop().to_string();
    std::string dataname_mu = in.pop().to_string();
    getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
    getfem::mesh_fem &mf_vm = gfi_mf->mesh_fem();
    std::string stresca = "Von Mises";
    if (in.remaining()) stresca = in.pop().to_string();
    bool tresca = false;
    if (cmd_strmatch(stresca, "Von Mises"))
      tresca = false;
    else if (cmd_strmatch(stresca, "Tresca"))
      tresca = true;
    else THROW_BADARG("bad option: " << stresca);

    getfem::model_real_plain_vector VMM(mf_vm.nb_dof());
    getfem::compute_isotropic_linearized_Von_Mises_or_Tresca
      (md->model(), varname, dataname_lambda, dataname_mu, mf_vm, VMM, tresca);
    out.pop().from_dcvector(VMM);
  } else bad_cmd(cmd);
}
