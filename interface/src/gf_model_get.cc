// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2011 Yves Renard.
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
#include <getfemint_mesh_im.h>
#include <getfem/getfem_continuation.h>

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

    /*@GET ('brick term rhs', @int ind_brick[, @int ind_term, @int sym, @int ind_iter])
      Gives the access to the part of the right hand side of a term
      of a particular nonlinear brick. Does not account of the eventual
      time dispatcher. An assembly of the rhs has to be done first.
      `ind_brick` is the brick index. `ind_term` is the index of the
      term inside the brick (default value : @MATLAB{1}@PYTHON{0}@SCILAB{1}).
      `sym` is to access to the second right hand side of for symmetric
      terms acting on two different variables (default is 0).
      `ind_iter` is the iteration number when time dispatchers are
      used (default is @MATLAB{1}@PYTHON{0}@SCILAB{1}).
      @*/
    sub_command
      ("brick term rhs", 1, 4, 0, 1,

       size_type ind_brick = in.pop().to_integer() - config::base_index();
       size_type ind_term = 0;
       if (in.remaining())
	 ind_term = in.pop().to_integer() - config::base_index();
       bool sym = false;
       if (in.remaining())
	 sym = (in.pop().to_integer() != 0);
       size_type ind_iter = 0;
       if (in.remaining())
	 ind_iter = in.pop().to_integer() - config::base_index();

       if (md->model().is_complex())
	 out.pop().from_dcvector(md->model().complex_brick_term_rhs(ind_brick, ind_term, sym, ind_iter));
       else
	 out.pop().from_dcvector(md->model().real_brick_term_rhs(ind_brick, ind_term, sym, ind_iter));
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
       size_type ind_brick = in.pop().to_integer() - config::base_index();
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
      'build_rhs', 'build_matrix' or 'pseudo_potential' (in that case,
      the pseudo_potential is returned).
      The default is to build the whole
      tangent linear system (matrix and rhs). This function is usefull
      to solve your problem with you own solver. @*/
    sub_command
      ("assembly", 0, 1, 0, 1,
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
       else if (cmd_strmatch(option, "pseudo potential") ||
		cmd_strmatch(option, "pseudo_potential"))
	 version = getfem::model::BUILD_PSEUDO_POTENTIAL;
       else THROW_BADARG("bad option: " << option);
       md->model().assembly(version);
       if (version == getfem::model::BUILD_PSEUDO_POTENTIAL)
	 out.pop().from_scalar(md->model().pseudo_potential());
       );


    /*@GET @CELL{nbit, converged} = ('solve'[, ...])
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
    - 'diverged_res', @float RES
       set the threshold value of the residual beyond which the iterative
       method is considered to diverge (default is 1e200).
    - 'lsolver', @str SOLVER_NAME
       select explicitely the solver used for the linear systems (the
       default value is 'auto', which lets getfem choose itself).
       Possible values are 'superlu', 'mumps' (if supported),
       'cg/ildlt', 'gmres/ilu' and 'gmres/ilut'.
    - 'lsearch', @str LINE_SEARCH_NAME
       select explicitely the line search method used for the linear systems (the
       default value is 'default').
       Possible values are 'simplest', 'systematic', 'quadratic' or 'basic'.
    - 'with pseudo potential'
      for nonlinear problems, the criterion of the line search will
      be a pseudo potential instead of the residual. Still experimental since
      not all bricks define a pseudo potential.

      Return the number of iterations, if a iterative method is used. @*/
    sub_command
      ("solve", 0, 15, 0, 2,
       getfemint::interruptible_iteration iter;
       std::string lsolver = "auto";
       std::string lsearch = "default";
       bool with_pseudo_pot = false;
       scalar_type alpha_mult = 0.5;
       scalar_type alpha_min = 1.0/1000.0;
       scalar_type alpha_max_ratio = 6.0/5.0;
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
	 } else if (cmd_strmatch(opt, "diverged_res")) {
	   if (in.remaining()) iter.set_diverged_residual(in.pop().to_scalar());
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "lsolver")) {
	   if (in.remaining()) lsolver = in.pop().to_string();
	   else THROW_BADARG("missing solver name for " << opt);
	 } else if (cmd_strmatch(opt, "lsearch")) {
	   if (in.remaining()) lsearch = in.pop().to_string();
	   else THROW_BADARG("missing line search name for " << opt);
	 } else if (cmd_strmatch(opt, "alpha mult")) {
	   if (in.remaining()) alpha_mult = in.pop().to_scalar();
	   else THROW_BADARG("missing line search name for " << opt);
	 } else if (cmd_strmatch(opt, "alpha min")) {
	   if (in.remaining()) alpha_min = in.pop().to_scalar();
	   else THROW_BADARG("missing line search name for " << opt);
	 } else if (cmd_strmatch(opt, "alpha max ratio")) {
	   if (in.remaining()) alpha_max_ratio = in.pop().to_scalar();
	   else THROW_BADARG("missing line search name for " << opt);
	 } else THROW_BADARG("bad option: " << opt);
       }
       
       gmm::default_newton_line_search default_ls(size_type(-1), alpha_mult);
       gmm::simplest_newton_line_search simplest_ls(size_type(-1), alpha_max_ratio, alpha_min, alpha_mult);
       gmm::systematic_newton_line_search systematic_ls(size_type(-1), alpha_min, alpha_mult);
       gmm::basic_newton_line_search basic_ls(size_type(-1), alpha_min, alpha_mult);
       gmm::quadratic_newton_line_search quadratic_ls(size_type(-1));

       gmm::abstract_newton_line_search *ls = 0;

       if (lsearch == "default")
	 ls = &default_ls;
       else if (lsearch == "simplest")
	 ls = &simplest_ls;
       else if (lsearch == "basic")
	 ls = &basic_ls;
       else if (lsearch == "systematic")
	 ls = &systematic_ls;
       else if (lsearch == "quadratic")
	 ls = &quadratic_ls;
       else GMM_ASSERT1(false, "unknown line search");


       if (!md->model().is_complex()) {
	 getfem::standard_solve(md->model(), iter,
				getfem::rselect_linear_solver(md->model(),
							      lsolver),
				*ls, with_pseudo_pot);
       } else {
	 getfem::standard_solve(md->model(), iter,
				getfem::cselect_linear_solver(md->model(),
							      lsolver),
				*ls, with_pseudo_pot);
       }
       if (out.remaining()) out.pop().from_integer(int(iter.get_iteration()));
       if (out.remaining()) out.pop().from_integer(int(iter.converged()));
       );


    /*@FUNC E = ('init Moore-Penrose continuation', @str dataname_parameter[,@str dataname_init, @str dataname_final, @str dataname_current], @scalar init_dir[, ...])
    Initialize the Moore-Penrose continuation (for more details about the
    continuation see the Getfem++ user documentation): The variable
    `dataname_parameter` should parametrize the model. If the parametrization
    is done via some vector datum, `dataname_init` and `dataname_final`
    should store two given values of this datum determining the
    parametrization, and `dataname_current` serves for actual values of this
    datum. Return a tangent corresponding to the solution branch at the
    current solution and the current value of the parameter, and an initial
    step size for the continuation. Direction of the computed tangent with
    respect to the parameter is determined by the sign of `init_dir`.
    
    Additional options:
    
    - 'lsolver', @str SOLVER_NAME
       name of a solver to be used for the incorporated linear system (the
       default value is 'auto', which lets getfem choose itself); possible
       values are 'superlu', 'mumps' (if supported), 'cg/ildlt', 'gmres/ilu'
       and 'gmres/ilut';
    - 'noisy' or 'very_noisy'
       determines how detailed information has to be displayed during the
       computations (residual values etc.);
    - 'epsilon', @scalar EPS
       increment to be used to compute the incorporated finite
       difference (the default value is 1e-8);
    - 'max_res_solve', @scalar RES_SOLVE
       target residual value for the linear system to be solved (the default
       value is 1e-7);
    - 'h_init', @scalar HIN
       initial step size (the default value is 1e-2).@*/
    sub_command
      ("init Moore-Penrose continuation", 2, 14, 0, 3,

       bool with_parametrized_data = false;
       std::string dataname_parameter = in.pop().to_string();
       std::string dataname_init; std::string dataname_final;
       std::string dataname_current;

       mexarg_in argin = in.pop();
       if (argin.is_string()) {
	 with_parametrized_data = true;
	 dataname_init = argin.to_string();
	 dataname_final = in.pop().to_string();
	 dataname_current = in.pop().to_string();
	 argin = in.pop();
       }
       scalar_type t_gamma = argin.to_scalar();
       
       std::string lsolver = "auto";
       size_type maxit = 10;
       size_type thrit = 8;
       scalar_type maxres = 1.e-6;
       scalar_type maxdiff = 1.e-6;
       scalar_type minang = 0.9; 
       scalar_type h_init = 1.e-2;
       scalar_type h_max = 1.e-1;
       scalar_type h_min = 1.e-5;
       scalar_type h_inc = 1.3;
       scalar_type h_dec = 0.5;
       scalar_type epsilon = 1.e-8;
       scalar_type maxres_solve = 1.e-7;
       int noisy = 0;
       
       while (in.remaining() && in.front().is_string()) {
	 std::string opt = in.pop().to_string();
	 if (cmd_strmatch(opt, "lsolver"))  {
	   if (in.remaining()) lsolver = in.pop().to_string();
	   else THROW_BADARG("missing name for " << opt);
	 } else if (cmd_strmatch(opt, "noisy")) noisy = 1;
	 else if (cmd_strmatch(opt, "very noisy") ||
		  cmd_strmatch(opt, "very_noisy")) noisy = 2;
	 else if (cmd_strmatch(opt, "epsilon")) {
	   if (in.remaining()) epsilon = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "max_res_solve")) {
	   if (in.remaining()) maxres_solve = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "h_init")) {
	   if (in.remaining()) h_init = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else THROW_BADARG("bad option: " << opt);
       }

       if (md->model().is_complex())
	 THROW_BADARG("Sorry, Moore-Penrose continuation " 
		      "has only a real version.");
       
       getfem::S_getfem_model S;
       if (with_parametrized_data) {
	  getfem::S_getfem_model S1
	   (md->model(), dataname_parameter, dataname_init, dataname_final,
	    dataname_current,
	    getfem::rselect_linear_solver(md->model(), lsolver), maxit,
	    thrit, maxres, maxdiff, minang, h_init, h_max, h_min, h_inc,
	    h_dec, epsilon, maxres_solve, noisy);
	  S = S1;
       }  else {
	 getfem::S_getfem_model S1
	   (md->model(), dataname_parameter,
	    getfem::rselect_linear_solver(md->model(), lsolver),
	    maxit, thrit, maxres, maxdiff, minang, h_init, h_max, h_min,
	    h_inc, h_dec, epsilon, maxres_solve, noisy);
	 S = S1;
       }

       size_type nbdof = md->model().nb_dof();
       std::vector<double> yy(nbdof);
       md->model().from_variables(yy);
       const getfem::model_real_plain_vector &GAMMA
       = md->model().real_variable(dataname_parameter);
       GMM_ASSERT1(gmm::vect_size(GAMMA) == 1,
		   "The continuation parameter should be a real scalar!");
       scalar_type gamma = GAMMA[0];
       scalar_type h;
       std::vector<double> tt_y(nbdof);

       getfem::init_Moore_Penrose_continuation(S, yy, gamma,
					       tt_y, t_gamma, h);
       out.pop().from_dcvector(tt_y);
       out.pop().from_scalar(t_gamma);
       out.pop().from_scalar(h);
       );


    /*@FUNC E = ('Moore-Penrose continuation', @str dataname_parameter[, @str dataname_init, @str dataname_final, @str dataname_current], @vec tangent, @scalar tangent_parameter, @scalar h[, ...])
    Compute one step of the Moore-Penrose continuation (for more details
    about the continuation see the Getfem++ user documentation): The variable
    `dataname_parameter` should parametrize the model. If the parametrization
    is done via some vector datum, `dataname_init` and `dataname_final`
    should store two given values of this datum determining the
    parametrization, and `dataname_current` serves for actual values of this
    datum. Take the current solution, the current value of the parameter, the
    tangent given by `tangent` and `tangent_parameter`, and the step size
    `h`, return a new tangent and step size for the next step (and save a new
    point in the model). If the returned step size equals zero, the
    continuation has failed.

    Additional options:

    - 'lsolver', @str SOLVER_NAME
       name of a solver to be used for the incorporated linear systems (the
       default value is 'auto', which lets getfem choose itself); possible
       values are 'superlu', 'mumps' (if supported), 'cg/ildlt', 'gmres/ilu'
       and 'gmres/ilut';
    - 'noisy' or 'very_noisy'
       determines how detailed information has to be displayed during the
       process (residual values etc.);
    - 'max_iter', @int NIT
       maximum number of iterations allowed in the correction (the default
       value is 10);
    - 'thr_iter', @int TIT
       threshold number of iterations of the correction for enlarging the
       step size (the default value is 8);
    - 'max_res', @scalar RES
       target residual value of the new point (the default value is 1e-6);
    - 'max_diff', @scalar DIFF
       determines a criterion of convergence to the new tangent vector (the
       default value is 1e-6);
    - 'min_ang', @scalar ANG
       minimal value of the cosine of the angle between tangents to the
       solution curve at the old point and the new one (the default value is
       0.9);
    - 'h_init', @scalar HIN
       initial step size (the default value is 1e-2);
    - 'h_max', @scalar HMAX
       maximal step size (the default value is 1e-1);
    - 'h_min', @scalar HMIN
       minimal step size (the default value is 1e-5);
    - 'h_inc', @scalar HINC
       factor for enlarging the step size (the default value is 1.3);
    - 'h_dec', @scalar HDEC
       factor for diminishing the step size (the default value is 0.5);
    - 'epsilon', @scalar EPS
       increment to be used to compute the incorporated finite
       differences (the default value is 1e-8);
    - 'max_res_solve', @scalar RES_SOLVE
       target residual value for the linear systems to be solved (the default
       value is 1e-7).@*/
    sub_command
      ("Moore-Penrose continuation", 4, 34, 0, 3,

       bool with_parametrized_data = false;
       std::string dataname_parameter = in.pop().to_string();
       std::string dataname_init; std::string dataname_final;
       std::string dataname_current;

       mexarg_in argin = in.pop();
       if (argin.is_string()) {
	 with_parametrized_data = true;
	 dataname_init = argin.to_string();
	 dataname_final = in.pop().to_string();
	 dataname_current = in.pop().to_string();
	 argin = in.pop();
       }

       darray t_y = argin.to_darray();
       scalar_type t_gamma = in.pop().to_scalar();
       scalar_type h = in.pop().to_scalar();
       std::string lsolver = "auto";
       size_type maxit = 10;
       size_type thrit = 8;
       scalar_type maxres = 1.e-6;
       scalar_type maxdiff = 1.e-6;
       scalar_type minang = 0.9; 
       scalar_type h_init = 1.e-2;
       scalar_type h_max = 1.e-1;
       scalar_type h_min = 1.e-5;
       scalar_type h_inc = 1.3;
       scalar_type h_dec = 0.5;
       scalar_type epsilon = 1.e-8;
       scalar_type maxres_solve = 1.e-7;
       int noisy = 0;
       
       while (in.remaining() && in.front().is_string()) {
	 std::string opt = in.pop().to_string();
	 if (cmd_strmatch(opt, "lsolver"))  {
	   if (in.remaining()) lsolver = in.pop().to_string();
	   else THROW_BADARG("missing name for " << opt);
	 } else if (cmd_strmatch(opt, "max_iter")) {
	   if (in.remaining()) maxit = in.pop().to_integer();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "thr_iter")) {
	   if (in.remaining()) thrit = in.pop().to_integer();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "max_res")) {
	   if (in.remaining()) maxres = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "max_diff")) {
	   if (in.remaining()) maxdiff = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "min_ang")) {
	   if (in.remaining()) minang = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "h_init")) {
	   if (in.remaining()) h_init = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "h_max")) {
	   if (in.remaining()) h_max = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "h_min")) {
	   if (in.remaining()) h_min = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "h_inc")) {
	   if (in.remaining()) h_inc = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "h_dec")) {
	   if (in.remaining()) h_dec = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "epsilon")) {
	   if (in.remaining()) epsilon = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "max_res_solve")) {
	   if (in.remaining()) maxres_solve = in.pop().to_scalar();
	   else THROW_BADARG("missing value for " << opt);
	 } else if (cmd_strmatch(opt, "noisy")) noisy = 1;
	 else if (cmd_strmatch(opt, "very noisy") ||
		  cmd_strmatch(opt, "very_noisy")) noisy = 2;
	 else THROW_BADARG("bad option: " << opt);
       }
       
       if (md->model().is_complex())
	 THROW_BADARG("Sorry, Moore-Penrose continuation "
		      "has only a real version.");

       size_type nbdof = md->model().nb_dof();
       std::vector<double> yy(nbdof);
       md->model().from_variables(yy);
       const getfem::model_real_plain_vector &GAMMA
       = md->model().real_variable(dataname_parameter);
       GMM_ASSERT1(gmm::vect_size(GAMMA) == 1,
		   "The continuation parameter should be a real scalar!");
       scalar_type gamma = GAMMA[0];

       getfem::S_getfem_model S;
       if (with_parametrized_data) {
	 getfem::S_getfem_model S1
	   (md->model(), dataname_parameter, dataname_init, dataname_final,
	    dataname_current,
	    getfem::rselect_linear_solver(md->model(), lsolver), maxit,
	    thrit, maxres, maxdiff, minang, h_init, h_max, h_min, h_inc,
	    h_dec, epsilon, maxres_solve, noisy);
	 S = S1;
       }
       else {
	 getfem::S_getfem_model S1
	   (md->model(), dataname_parameter,
	    getfem::rselect_linear_solver(md->model(), lsolver),
	    maxit, thrit, maxres, maxdiff, minang, h_init, h_max, h_min,
	    h_inc, h_dec, epsilon, maxres_solve, noisy);
	 S = S1;
       }

       std::vector<double> tt_y(nbdof);
       gmm::copy(t_y, tt_y);
       getfem::Moore_Penrose_continuation(S, yy, gamma, tt_y, t_gamma, h);
       out.pop().from_dcvector(tt_y);
       out.pop().from_scalar(t_gamma);
       out.pop().from_scalar(h);
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
       (md->model(), varname, 
	abstract_hyperelastic_law_from_name
	(lawname, gfi_mf->mesh_fem().linked_mesh().dim()),
	dataname, gfi_mf->mesh_fem(), VMM, tresca);
       out.pop().from_dcvector(VMM);
       );



    /*@GET V = ('compute plasticity Von Mises or Tresca', @str datasigma, @tmf mf_vm[, @str version])
      Compute on `mf_vm` the Von-Mises or the Tresca stress of a field for plasticity and return it into the vector V.
      `datasigma` is a vector which contains the stress constraints values supported by the mesh.  
      `version` should be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default).@*/
    sub_command
      ("compute elastoplasticity Von Mises or Tresca", 2, 3, 0, 1,
       std::string datasigma = in.pop().to_string();
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
       getfem::compute_elastoplasticity_Von_Mises_or_Tresca
       (md->model(), datasigma, gfi_mf->mesh_fem(), VMM, tresca);
       out.pop().from_dcvector(VMM);
       );




        /*@GET ('compute plasticity constraints', @tmim mim, @str varname, @str projname, @str datalambda, @str datamu, @str datathreshold, @str datasigma)
      Compute and save the stress constraints sigma for other hypothetical iterations. 
      'mim' is the integration method to use for the computation.
      'varname' is the main variable of the problem.
      'projname' is the type of projection to use. For the moment it could only be 'Von Mises' or 'VM'.
      'datalambda' and 'datamu' are the Lame coefficients of the material.
      'datasigma' is a vector which will contains the new stress constraints values.@*/
    sub_command
      ("elastoplasticity next iter", 7, 7, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();;
       std::string varname = in.pop().to_string();
       std::string projname = in.pop().to_string();
       std::string datalambda = in.pop().to_string();
       std::string datamu = in.pop().to_string();
       std::string datathreshold = in.pop().to_string();
       std::string datasigma = in.pop().to_string();

      
       getfem::elastoplasticity_next_iter
       (md->model(), gfi_mim->mesh_im(), varname,
	  abstract_constraints_projection_from_name(projname), 
	  datalambda, datamu, datathreshold, datasigma);
       );



    
    /*@GET V = ('compute plastic part', @tmim mim, @tmf mf_pl, @str varname, @str projname, @str datalambda, @str datamu, @str datathreshold, @str datasigma)
      Compute on `mf_pl` the plastic part and return it into the vector V.
      `datasigma` is a vector which contains the stress constraints values supported by the mesh.@*/
     sub_command
      ("compute plastic part", 8, 8, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();;
       getfemint_mesh_fem *gfi_mf = in.pop().to_getfemint_mesh_fem();
       std::string varname = in.pop().to_string();
       std::string projname = in.pop().to_string();
       std::string datalambda = in.pop().to_string();
       std::string datamu = in.pop().to_string();
       std::string datathreshold = in.pop().to_string();
       std::string datasigma = in.pop().to_string();

       getfem::model_real_plain_vector plast((gfi_mf->mesh_fem()).nb_dof());
       getfem::compute_plastic_part
       (md->model(), gfi_mim->mesh_im(),  gfi_mf->mesh_fem(), varname,
	  abstract_constraints_projection_from_name(projname), 
	datalambda, datamu, datathreshold, datasigma, plast);
       out.pop().from_dcvector(plast);
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
