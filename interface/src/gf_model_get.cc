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

/**\file gf_model_get.cc
   \brief getfemint_model getter.
*/

#include <getfemint.h>
#include <getfemint_misc.h>
#include <getfemint_models.h>
#include <getfem/getfem_model_solvers.h>
#include <getfemint_mdbrick.h>
#include <getfemint_mesh_fem.h>

using namespace getfemint;


#define RETURN_SPARSE(realmeth, cplxmeth)                              \
  if (!md->is_complex()) {                                             \
    gf_real_sparse_by_col M(gmm::mat_nrows(md->model().realmeth),    \
                            gmm::mat_ncols(md->model().realmeth));   \
    gmm::copy(md->model().realmeth, M);                              \
    out.pop().from_sparse(M);                                          \
  } else {                                                             \
    gf_cplx_sparse_by_col M(gmm::mat_nrows(md->model().cplxmeth),    \
                            gmm::mat_ncols(md->model().cplxmeth));   \
    gmm::copy(md->model().cplxmeth, M);                              \
    out.pop().from_sparse(M);                                          \
  }

#define RETURN_VECTOR(realmeth, cplxmeth)                      \
  if (!md->is_complex()) {                                     \
    out.pop().from_dcvector(md->model().realmeth);             \
  } else {                                                     \
    out.pop().from_dcvector(md->model().cplxmeth);             \
  }

/*@GFDOC
  Get information from a model object.
@*/


// Object for the declaration of a new sub-command.

struct sub_gf_md_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfemint_model *md) = 0;
};

typedef boost::intrusive_ptr<sub_gf_md_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_md_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfemint_model *md)				\
      { dummy_func(in); dummy_func(out);  dummy_func(md); code }	\
    };									\
    psub_command psubc = new subc;					\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           


void gf_model_get(getfemint::mexargs_in& m_in,
		  getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@GET b = ('is_complex')
      Return 0 is the model is real, 1 if it is complex.@*/
    sub_command
      ("is_complex", 0, 0, 0, 1,
       out.pop().from_integer(md->is_complex());
       );


    /*@GET T = ('tangent_matrix')
      Return the tangent matrix stored in the model .@*/
    sub_command
      ("tangent_matrix", 0, 0, 0, 1,
       RETURN_SPARSE(real_tangent_matrix(), complex_tangent_matrix());
       );


    /*@GET ('rhs')
      Return the right hand side of the tangent problem.@*/
    sub_command
      ("rhs", 0, 0, 0, 1,
       RETURN_VECTOR(real_rhs(), complex_rhs());
       );


    /*@GET z = ('memsize')
      Return a rough approximation of the amount of memory (in bytes) used by
      the model.@*/
    sub_command
      ("memsize", 0, 0, 0, 1,
       out.pop().from_integer(int(md->memsize()));
       );


    /*@GET ('listvar')
      print to the output the list of variables and constants of the model.@*/
    sub_command
      ("listvar", 0, 0, 0, 0,
       md->model().listvar(infomsg());
       );


    /*@GET ('listbricks')
      print to the output the list of bricks of the model.@*/
    sub_command
      ("listbricks", 0, 0, 0, 0,
       md->model().listbricks(infomsg(), config::base_index());
       );


    /*@GET V = ('variable', @str name[, @int niter])
      Gives the value of a variable or data.@*/
    sub_command
      ("variable", 1, 2, 0, 1,
       std::string name = in.pop().to_string();
       size_type niter = 0;
       if (in.remaining())
	 niter = in.pop().to_integer(0,10) - config::base_index();
       RETURN_VECTOR(real_variable(name, niter),
		     complex_variable(name, niter));
       );


    /*@GET name = ('mult varname Dirichlet', @int ind_brick)
      Gives the name of the multiplier variable for a Dirichlet brick.
      If the brick is not a Dirichlet condition with multiplier brick,
      this function has an undefined behavior@*/
    sub_command
      ("mult varname Dirichlet", 1, 1, 0, 1,
       size_type ind_brick = in.pop().to_integer();
       out.pop().from_string
       (getfem::mult_varname_Dirichlet(md->model(), ind_brick).c_str());
       );


    /*@GET I = ('interval of variable', @str varname)
      Gives the interval of the variable `varname` in the linear system of
      the model.@*/
    sub_command
      ("interval of variable", 1, 1, 0, 1,
       std::string name = in.pop().to_string();
       const gmm::sub_interval I = md->model().interval_of_variable(name);
       iarray opids = out.pop().create_iarray_h(2);
       opids[0] = int(I.first() + config::base_index());
       opids[1] = int(I.size());
       );


    /*@GET V = ('from variables')
      Return the vector of all the degrees of freedom of the model consisting
      of the concatenation of the variables of the model (usefull
      to solve your problem with you own solver). @*/
    sub_command
      ("from variables", 0, 0, 0, 1,
       if (!md->is_complex()) {
	 std::vector<double> V(md->model().nb_dof());
	 md->model().from_variables(V);
	 out.pop().from_dcvector(V);
       } else {
	 std::vector<std::complex<double> > V(md->model().nb_dof());
	 md->model().from_variables(V);
	 out.pop().from_dcvector(V);
       }
       );


    /*@GET ('assembly'[, @str option])
      Assembly of the tangent system taking into account the terms
      from all bricks. `option`, if specified, should be 'build_all',
      'build_rhs' or 'build_matrix'. The default is to build the whole
      tangent linear system (matrix and rhs). This function is usefull
      to solve your problem with you own solver. @*/
    sub_command
      ("assembly", 0, 1, 0, 0,
       std::string option = "build_all";
       if (in.remaining()) option = in.pop().to_string();
       getfem::model::build_version version = getfem::model::BUILD_ALL;
       if (cmd_strmatch(option, "build all") ||
	   cmd_strmatch(option, "build_all"))
	 version = getfem::model::BUILD_ALL;
       else if (cmd_strmatch(option, "build rhs") ||
		cmd_strmatch(option, "build_rhs"))
	 version = getfem::model::BUILD_RHS;
       else if (cmd_strmatch(option, "build matrix") ||
		cmd_strmatch(option, "build_matrix"))
	 version = getfem::model::BUILD_MATRIX;
       else THROW_BADARG("bad option: " << option);
       md->model().assembly(version);
       );


    /*@GET ('solve'[, ...])
    Run the standard getfem solver.

    Note that you should be able to use your own solver if you want
    (it is possible to obtain the tangent matrix and its right hand
    side with the MODEL:GET('tangent matrix') etc.).

    Various options can be specified:

    - 'noisy' or 'very_noisy'
       the solver will display some information showing the progress
       (residual values etc.).
    - 'max_iter', @int NIT
       set the maximum iterations numbers.
    - 'max_res', @float RES
       set the target residual value.
    - 'lsolver', @str SOLVER_NAME
       select explicitely the solver used for the linear systems (the
       default value is 'auto', which lets getfem choose itself).
       Possible values are 'superlu', 'mumps' (if supported),
       'cg/ildlt', 'gmres/ilu' and 'gmres/ilut'.
    - 'with pseudo potential'
      for nonlinear problems, the criterion of the line search will
      be a pseudo potential instead of the residual. Still experimental since
      not all bricks define a pseudo potential. @*/
    sub_command
      ("solve", 0, 7, 0, 0,
       getfemint::interruptible_iteration iter;
       std::string lsolver = "auto";
       bool with_pseudo_pot = false;
       while (in.remaining() && in.front().is_string()) {
	 std::string opt = in.pop().to_string();
	 if (cmd_strmatch(opt, "noisy")) iter.set_noisy(1);
	 else if (cmd_strmatch(opt, "with pseudo potential"))
	   with_pseudo_pot = true;
	 else if (cmd_strmatch(opt, "very noisy") ||
		  cmd_strmatch(opt, "very_noisy")) iter.set_noisy(2);
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
       
       gmm::default_newton_line_search ls;
       
       if (!md->model().is_complex()) {
	 getfem::standard_solve(md->model(), iter,
				getfem::rselect_linear_solver(md->model(),
							      lsolver),
				ls, with_pseudo_pot);
       } else {
	 getfem::standard_solve(md->model(), iter,
				getfem::cselect_linear_solver(md->model(),
							      lsolver),
				ls, with_pseudo_pot);
       }
       );


    /*@GET V = ('compute isotropic linearized Von Mises or Tresca', @str varname, @str dataname_lambda, @str dataname_mu, @tmf mf_vm[, @str version])
      Compute the Von-Mises stress or the Tresca stress of a field (only
      valid for isotropic linearized elasticity in 3D). `version` should
      be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default). @*/
    sub_command
      ("compute isotropic linearized Von Mises or Tresca", 4, 5, 0, 1,
       std::string varname = in.pop().to_string();
       std::string dataname_lambda = in.pop().to_string();
       std::string dataname_mu = in.pop().to_string();
       getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
       std::string stresca = "Von Mises";
       if (in.remaining()) stresca = in.pop().to_string();
       bool tresca = false;
       if (cmd_strmatch(stresca, "Von Mises") ||
	   cmd_strmatch(stresca, "Von_Mises"))
	 tresca = false;
       else if (cmd_strmatch(stresca, "Tresca"))
	 tresca = true;
       else THROW_BADARG("bad option \'version\': " << stresca);
       
       getfem::model_real_plain_vector VMM((gfi_mf->mesh_fem()).nb_dof());
       getfem::compute_isotropic_linearized_Von_Mises_or_Tresca
       (md->model(), varname, dataname_lambda, dataname_mu, gfi_mf->mesh_fem(), VMM, tresca);
       out.pop().from_dcvector(VMM);
       );

    /*@GET V = ('compute Von Mises or Tresca', @str varname, @str lawname, @str dataname, @tmf mf_vm[, @str version])
      Compute on `mf_vm` the Von-Mises stress or the Tresca stress of a field
      for nonlinear elasticity in 3D. `lawname` is the constitutive law which
      could be 'SaintVenant Kirchhoff', 'Mooney Rivlin' or 'Ciarlet Geymonat'.
      `dataname` is a vector of parameters for the constitutive law. Its length
      depends on the law. It could be a short vector of constant values or a
      vector field described on a finite element method for variable coefficients.
      `version` should be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default). @*/
    sub_command
      ("compute Von Mises or Tresca", 4, 5, 0, 1,
       std::string varname = in.pop().to_string();
       std::string lawname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
       std::string stresca = "Von Mises";
       if (in.remaining()) stresca = in.pop().to_string();
       bool tresca = false;
       if (cmd_strmatch(stresca, "Von Mises") ||
	   cmd_strmatch(stresca, "Von_Mises"))
	 tresca = false;
       else if (cmd_strmatch(stresca, "Tresca"))
	 tresca = true;
       else THROW_BADARG("bad option \'version\': " << stresca);
       
       getfem::model_real_plain_vector VMM((gfi_mf->mesh_fem()).nb_dof());
       getfem::compute_Von_Mises_or_Tresca
       (md->model(), varname, abstract_hyperelastic_law_from_name(lawname),
	dataname, gfi_mf->mesh_fem(), VMM, tresca);
       out.pop().from_dcvector(VMM);
       );


    /*@GET M = ('matrix term', @int ind_brick, @int ind_term)
      Gives the matrix term ind_term of the brick ind_brick if it exists
      @*/
    sub_command
      ("matrix term", 2, 2, 0, 1,
       size_type ind_brick = in.pop().to_integer() - config::base_index();
       size_type ind_term = in.pop().to_integer() - config::base_index();
       RETURN_SPARSE(linear_real_matrix_term(ind_brick, ind_term),
		     linear_complex_matrix_term(ind_brick, ind_term));
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tmodel.

      This can be used to perform comparisons between two
      different @tmodel objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false, "Sorry, function to be done");
       // std::string s = ...;
       // out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
      displays a short summary for a @tmodel object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       if (md->is_complex()) infomsg() << "Complex "; else infomsg()<< "Real ";
       infomsg() << "gfModel object with " << md->model().nb_dof()
       << " degrees of freedom\n";
       );
    
  }

 
  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfemint_model *md  = m_in.pop().to_getfemint_model();
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
