/*===========================================================================

 Copyright (C) 2009-2020 Yves Renard.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

/**\file gf_model_get.cc
   \brief getfemint_model getter.
*/

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfemint_misc.h>
#include <getfem/getfem_model_solvers.h>
#include <getfem/getfem_generic_assembly.h>
#include <getfem/getfem_nonlinear_elasticity.h>
#include <getfem/getfem_plasticity.h>
#include <getfem/getfem_contact_and_friction_large_sliding.h>
#include <getfem/getfem_im_data.h>

using namespace getfemint;


#define RETURN_SPARSE(realmeth, cplxmeth)                            \
  if (!md->is_complex()) {                                           \
    gf_real_sparse_by_col M(gmm::mat_nrows(md->realmeth),            \
                            gmm::mat_ncols(md->realmeth));           \
    gmm::copy(md->realmeth, M);                                      \
    out.pop().from_sparse(M);                                        \
  } else {                                                           \
    gf_cplx_sparse_by_col M(gmm::mat_nrows(md->cplxmeth),            \
                            gmm::mat_ncols(md->cplxmeth));           \
    gmm::copy(md->cplxmeth, M);                                      \
    out.pop().from_sparse(M);                                        \
  }

#define RETURN_VECTOR(realmeth, cplxmeth)                      \
  if (!md->is_complex()) {                                     \
    out.pop().from_dcvector(md->realmeth);                     \
  } else {                                                     \
    out.pop().from_dcvector(md->cplxmeth);                     \
  }

/*@GFDOC
  Get information from a model object.
@*/


// Object for the declaration of a new sub-command.

struct sub_gf_md_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfem::model *md) = 0;
};

typedef std::shared_ptr<sub_gf_md_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_md_get {                                \
      virtual void run(getfemint::mexargs_in& in,                       \
                       getfemint::mexargs_out& out,                     \
                       getfem::model *md)                               \
      { dummy_func(in); dummy_func(out);  dummy_func(md); code }        \
    };                                                                  \
    psub_command psubc = std::make_shared<subc>();                      \
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;         \
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;     \
    subc_tab[cmd_normalize(name)] = psubc;                              \
  }

static void filter_lawname(std::string &lawname) {
  for (auto &c : lawname)
    { if (c == ' ') c = '_'; if (c >= 'A' && c <= 'Z') c = char(c+'a'-'A'); }
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


    /*@GET T = ('nbdof')
      Return the total number of degrees of freedom of the model.@*/
    sub_command
      ("nbdof", 0, 0, 0, 1,
       out.pop().from_integer(int(md->nb_dof()));
       );

    /*@GET dt = ('get time step')
      Gives the value of the time step. @*/
    sub_command
      ("get time step", 0, 0, 0, 1,
       out.pop().from_scalar(md->get_time_step());
       );

    /*@GET t = ('get time')
      Give the value of the data `t` corresponding to the current time.
      @*/
    sub_command
      ("get time", 0, 0, 0, 1,
       out.pop().from_scalar(md->get_time());
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

       if (md->is_complex())
         out.pop().from_dcvector(md->complex_brick_term_rhs(ind_brick, ind_term, sym, ind_iter));
       else
         out.pop().from_dcvector(md->real_brick_term_rhs(ind_brick, ind_term, sym, ind_iter));
       );


    /*@GET z = ('memsize')
      Return a rough approximation of the amount of memory (in bytes) used by
      the model.@*/
    sub_command
      ("memsize", 0, 0, 0, 1,
       size_type msz = 0;
       /* Rough estimate. Quite optimistic with sparse matrices! */
       if (!md->is_complex()) {
         size_type szd = sizeof(double);
         size_type szs = sizeof(size_type);
         msz = gmm::nnz(md->real_tangent_matrix()) * (szd + szs)
           + gmm::vect_size(md->real_rhs()) *  szd * 3
           + sizeof(getfem::model);
       }
       else {
         size_type szc = sizeof(std::complex<double>);
         size_type szs = sizeof(size_type);
         msz = gmm::nnz(md->complex_tangent_matrix()) * (szc + szs)
           + gmm::vect_size(md->complex_rhs()) *  szc * 3
           + sizeof(getfem::model);
       }
       out.pop().from_integer(int(msz));
       );


    /*@GET ('variable list')
      print to the output the list of variables and constants of the model.@*/
    sub_command
      ("variable list", 0, 0, 0, 0,
       md->listvar(infomsg());
       );


    /*@GET ('brick list')
      print to the output the list of bricks of the model.@*/
    sub_command
      ("brick list", 0, 0, 0, 0,
       md->listbricks(infomsg(), config::base_index());
       );

    /*@GET ('list residuals')
      print to the output the residuals corresponding to all terms
      included in the model.@*/
    sub_command
      ("list residuals", 0, 0, 0, 0,
       md->listresiduals(cout);
       );


    /*@GET V = ('variable', @str name)
      Gives the value of a variable or data. @*/
    sub_command
      ("variable", 1, 1, 0, 1,
       std::string name = in.pop().to_string();
       RETURN_VECTOR(real_variable(name), complex_variable(name));
       );


    /*@GET V = ('interpolation', @str expr, {@tmf mf | @tmimd mimd | @vec pts,  @tmesh m}[, @int region[, @int extrapolation[, @int rg_source]]])
      Interpolate a certain expression with respect to the mesh_fem `mf`
      or the mesh_im_data `mimd` or the set of points `pts` on mesh `m`.
      The expression has to be valid according to the high-level generic
      assembly language possibly including references to the variables
      and data of the model.

      The options `extrapolation` and `rg_source` are specific to
      interpolations with respect to a set of points `pts`. @*/
    sub_command
      ("interpolation", 2, 6, 0, 1, // should be extended to complex models ...
       std::string expr = in.pop().to_string();
       getfem::base_vector result;
       if (is_meshfem_object(in.front())) {
         const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
         size_type rg = in.remaining() ? in.pop().to_integer() : size_type(-1);
         getfem::ga_interpolation_Lagrange_fem(*md, expr, *mf, result, rg);
       } else if (is_meshimdata_object(in.front())) {
         getfem::im_data *mimd = to_meshimdata_object(in.pop());
         size_type rg = in.remaining() ? in.pop().to_integer() : size_type(-1);
         getfem::ga_interpolation_im_data(*md, expr, *mimd, result, rg);
       } else {
         darray st = in.pop().to_darray();
         std::vector<double> PTS(st.begin(), st.end());
         const getfem::mesh *m = extract_mesh_object(in.pop());
         size_type N = m->dim();
         size_type nbpoints = gmm::vect_size(PTS) / N;
         getfem::base_node p(N);
         getfem::mesh_trans_inv mti(*m);
         for (size_type i = 0; i < nbpoints; ++i) {
           gmm::copy(gmm::sub_vector(PTS, gmm::sub_interval(i*N, N)), p);
           mti.add_point(p);
         }
         size_type rg = in.remaining() ? in.pop().to_integer()
                                       : size_type(-1);
         int extrapolation = in.remaining() ? in.pop().to_integer()
                                            : size_type(0);
         size_type rg_source = in.remaining() ? in.pop().to_integer()
                                              : size_type(-1);
         getfem::ga_interpolation_mti(*md, expr, mti,
                                      result, rg, extrapolation, rg_source);
       }
       out.pop().from_dcvector(result);
       );

    /*@GET V = ('local_projection', @tmim mim, @str expr, @tmf mf[, @int region])
      Make an elementwise L2 projection of an expression with respect
      to the mesh_fem `mf`. This mesh_fem has to be
      a discontinuous one.
      The expression has to be valid according to the high-level generic
      assembly language possibly including references to the variables
      and data of the model. @*/
    sub_command
      ("local_projection", 3, 4, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string expr = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       GMM_ASSERT1(!(mf->is_reduced()), "Sorry, cannot apply to reduced fems");
       size_type region = in.remaining() ? in.pop().to_integer():size_type(-1);

       getfem::base_vector result;
       getfem::ga_local_projection(*md, *mim, expr, *mf, result, region);
       out.pop().from_dcvector(result);
       );
    
    /*@GET mf = ('mesh fem of variable', @str name)
      Gives access to the `mesh_fem` of a variable or data.@*/
    sub_command
      ("mesh fem of variable", 1, 1, 0, 1,
       std::string name = in.pop().to_string();
       const getfem::mesh_fem &mf = md->mesh_fem_of_variable(name);
       id_type id = workspace().object(&mf);
       if (id == id_type(-1)) {
         id = store_meshfem_object(std::shared_ptr<getfem::mesh_fem>
                                   (std::shared_ptr<getfem::mesh_fem>(),
                                    const_cast<getfem::mesh_fem *>(&mf)));
         workspace().set_dependence(&mf, md);
       }
       out.pop().from_object_id(id, MESHFEM_CLASS_ID);
       );


    /*@GET name = ('mult varname Dirichlet', @int ind_brick)
      Gives the name of the multiplier variable for a Dirichlet brick.
      If the brick is not a Dirichlet condition with multiplier brick,
      this function has an undefined behavior@*/
    sub_command
      ("mult varname Dirichlet", 1, 1, 0, 1,
       size_type ind_brick = in.pop().to_integer() - config::base_index();
       out.pop().from_string
       (getfem::mult_varname_Dirichlet(*md, ind_brick).c_str());
       );


    /*@GET I = ('interval of variable', @str varname)
      Gives the interval of the variable `varname` in the linear system of
      the model.@*/
    sub_command
      ("interval of variable", 1, 1, 0, 1,
       std::string name = in.pop().to_string();
       const gmm::sub_interval I = md->interval_of_variable(name);
       iarray opids = out.pop().create_iarray_h(2);
       opids[0] = int(I.first() + config::base_index());
       opids[1] = int(I.size());
       );


    /*@GET V = ('from variables')
      Return the vector of all the degrees of freedom of the model consisting
      of the concatenation of the variables of the model (useful
      to solve your problem with you own solver). @*/
    sub_command
      ("from variables", 0, 0, 0, 1,
       if (!md->is_complex()) {
         std::vector<double> V(md->nb_dof());
         md->from_variables(V);
         out.pop().from_dcvector(V);
       } else {
         std::vector<std::complex<double> > V(md->nb_dof());
         md->from_variables(V);
         out.pop().from_dcvector(V);
       }
       );


    /*@GET ('assembly'[, @str option])
      Assembly of the tangent system taking into account the terms
      from all bricks. `option`, if specified, should be 'build_all',
      'build_rhs', 'build_matrix', 'build_rhs_with_internal',
      'build_matrix_condensed', 'build_all_condensed'.
      The default is to build the whole tangent linear system (matrix and rhs).
      This function is useful to solve your problem with you own solver. @*/
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
       else if (cmd_strmatch(option, "build rhs with internal") ||
                cmd_strmatch(option, "build_rhs_with_internal"))
         version = getfem::model::BUILD_RHS_WITH_INTERNAL;
       else if (cmd_strmatch(option, "build matrix condensed") ||
                cmd_strmatch(option, "build_matrix_condensed"))
         version = getfem::model::BUILD_MATRIX_CONDENSED;
       else if (cmd_strmatch(option, "build all condensed") ||
                cmd_strmatch(option, "build_all_condensed"))
         version = getfem::model::BUILD_ALL_CONDENSED;
       else THROW_BADARG("bad option: " << option);
       md->assembly(version);
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

      Return the number of iterations, if an iterative method is used.
      
      Note that it is possible to disable some variables
      (see MODEL:SET('disable variable') ) in order to
      solve the problem only with respect to a subset of variables (the
      disabled variables are then considered as data) for instance to
      replace the global Newton strategy with a fixed point one.

      @*/
    sub_command
      ("solve", 0, 17, 0, 2,
       getfemint::interruptible_iteration iter;
       std::string lsolver = "auto";
       std::string lsearch = "default";
       scalar_type alpha_max_ratio(-1);
       scalar_type alpha_min(-1);
       scalar_type alpha_mult(-1);
       scalar_type alpha_threshold_res(1e50);
       while (in.remaining() && in.front().is_string()) {
         std::string opt = in.pop().to_string();
         if (cmd_strmatch(opt, "noisy")) iter.set_noisy(1);
         else if (cmd_strmatch(opt, "very noisy") ||
                  cmd_strmatch(opt, "very_noisy")) iter.set_noisy(3);
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
           else THROW_BADARG("missing line search value for " << opt);
         } else if (cmd_strmatch(opt, "alpha min")) {
           if (in.remaining()) alpha_min = in.pop().to_scalar();
           else THROW_BADARG("missing line search value for " << opt);
         } else if (cmd_strmatch(opt, "alpha max ratio")) {
           if (in.remaining()) alpha_max_ratio = in.pop().to_scalar();
           else THROW_BADARG("missing line search value for " << opt);
         } else if (cmd_strmatch(opt, "alpha threshold res")) {
           if (in.remaining()) alpha_threshold_res = in.pop().to_scalar();
           else THROW_BADARG("missing line search value for " << opt);
         } else THROW_BADARG("bad option: " << opt);
       }

       // default values in sync with getfem_model_solvers.h
       if (alpha_max_ratio < scalar_type(0))
         alpha_max_ratio = (lsearch == "basic") ?  5.0/3.0 : 6.0/5.0;
       if (alpha_min < scalar_type(0))
         alpha_min = (lsearch == "systematic") ? 1.0/10000.0 : 1.0/1000.0;
       if (alpha_mult < scalar_type(0))
         alpha_mult = 3.0/5.0;

       // getfem::default_newton_line_search default_ls;
       getfem::newton_search_with_step_control default_ls;
       getfem::simplest_newton_line_search simplest_ls(size_type(-1), alpha_max_ratio,
                                                       alpha_min, alpha_mult, alpha_threshold_res);
       getfem::systematic_newton_line_search systematic_ls(size_type(-1), alpha_min, alpha_mult);
       getfem::basic_newton_line_search basic_ls(size_type(-1), alpha_max_ratio, alpha_min, alpha_mult);
       getfem::quadratic_newton_line_search quadratic_ls(size_type(-1));

       getfem::abstract_newton_line_search *ls = 0;

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


       if (!md->is_complex()) {
         getfem::standard_solve(*md, iter,
                                getfem::rselect_linear_solver(*md,
                                                              lsolver), *ls);
       } else {
         getfem::standard_solve(*md, iter,
                                getfem::cselect_linear_solver(*md,
                                                              lsolver), *ls);
       }
       if (out.remaining()) out.pop().from_integer(int(iter.get_iteration()));
       if (out.remaining()) out.pop().from_integer(int(iter.converged()));
       );


    /*@GET ('test tangent matrix'[, @scalar EPS[, @int NB[, @scalar scale]]])
      Test the consistency of the tangent matrix in some random positions
      and random directions (useful to test newly created bricks).
      `EPS` is the value of the small parameter for the finite difference
      computation of the derivative is the random direction (default is 1E-6).
      `NN` is the number of tests (default is 100). `scale` is a parameter
      for the random position (default is 1, 0 is an acceptable value) around
      the current position.
      Each dof of the random position is chosen in the range
      [current-scale, current+scale].
      @*/
    sub_command
      ("test tangent matrix", 0, 3, 0, 1,
       size_type nbdof = md->nb_dof();
       scalar_type EPS = 1E-6;
       if (in.remaining()) EPS = in.pop().to_scalar();
       scalar_type errmax = scalar_type(0);
       size_type NN = 100;
       if (in.remaining()) NN = in.pop().to_integer();
       scalar_type scalef = scalar_type(1);
       if (in.remaining()) scalef = in.pop().to_scalar();

       if (md->is_linear()) {
         GMM_WARNING1("Problem is linear, the test is not relevant");
       } else {
         if (md->is_complex()) {
           std::vector<complex_type> U(nbdof);
           std::vector<complex_type> dU(nbdof);
           std::vector<complex_type> DIR(nbdof);
           std::vector<complex_type> D1(nbdof);
           std::vector<complex_type> D2(nbdof);
           md->from_variables(U);
           for (size_type i = 0; i < NN; ++i) {
             gmm::fill_random(dU); gmm::scale(dU, complex_type(scalef));
             gmm::add(U, dU);
             gmm::fill_random(DIR);
             md->to_variables(dU);
             md->assembly(getfem::model::BUILD_ALL);
             gmm::copy(md->complex_rhs(), D2);
             gmm::mult(md->complex_tangent_matrix(), DIR, D1);
             gmm::add(gmm::scaled(DIR, complex_type(EPS)), dU);
             md->to_variables(dU);
             md->assembly(getfem::model::BUILD_RHS);
             gmm::add(gmm::scaled(md->complex_rhs(),
                                  -complex_type(1)), D2);
             gmm::scale(D2, complex_type(1)/complex_type(EPS));
             scalar_type err = gmm::vect_dist2(D1, D2);
             cout << "Error at step " << i << " : " << err << endl;
             errmax = std::max(err, errmax);
           }
           md->to_variables(U);
         } else {
           std::vector<scalar_type> U(nbdof);
           std::vector<scalar_type> dU(nbdof);
           std::vector<scalar_type> DIR(nbdof);
           std::vector<scalar_type> D1(nbdof);
           std::vector<scalar_type> D2(nbdof);
           md->from_variables(U);
           for (size_type i = 0; i < NN; ++i) {
             gmm::fill_random(dU); gmm::scale(dU, scalef);
             gmm::add(U, dU);
             gmm::fill_random(DIR);
             md->to_variables(dU);
             md->assembly(getfem::model::BUILD_ALL);
             gmm::copy(md->real_rhs(), D2);
             gmm::mult(md->real_tangent_matrix(), DIR, D1);
             gmm::add(gmm::scaled(DIR, EPS), dU);
             md->to_variables(dU);
             md->assembly(getfem::model::BUILD_RHS);
             gmm::add(gmm::scaled(md->real_rhs(),-scalar_type(1)), D2);
             gmm::scale(D2, scalar_type(1)/EPS);
             scalar_type err = gmm::vect_dist2(D1, D2);
             cout << "Error at step " << i << " : " << err << endl;
             errmax = std::max(err, errmax);
           }
           md->to_variables(U);
         }
       }
       out.pop().from_scalar(errmax);
       );


    /*@GET ('test tangent matrix term', @str varname1, @str varname2[, @scalar EPS[, @int NB[, @scalar scale]]])
      Test the consistency of a part of the tangent matrix in some
      random positions and random directions
      (useful to test newly created bricks).
      The increment is only made on variable `varname2` and tested on the
      part of the residual corresponding to `varname1`. This means that
      only the term (`varname1`, `varname2`) of the tangent matrix is tested.
      `EPS` is the value of the small parameter for the finite difference
      computation of the derivative is the random direction (default is 1E-6).
      `NN` is the number of tests (default is 100). `scale` is a parameter
      for the random position (default is 1, 0 is an acceptable value)
      around the current position.
      Each dof of the random position is chosen in the range
      [current-scale, current+scale].
      @*/
    sub_command
      ("test tangent matrix term", 2, 5, 0, 1,
       std::string varname1 = in.pop().to_string();
       std::string varname2 = in.pop().to_string();
       gmm::sub_interval I1 = md->interval_of_variable(varname1);
       gmm::sub_interval I2 = md->interval_of_variable(varname2);
       size_type nbdof1 = I1.size();
       size_type nbdof2 = I2.size();

       scalar_type EPS = 1E-6;
       if (in.remaining()) EPS = in.pop().to_scalar();
       scalar_type errmax = scalar_type(0);
       size_type NN = 100;
       if (in.remaining()) NN = in.pop().to_integer();
       scalar_type scalef = scalar_type(1);
       if (in.remaining()) scalef = in.pop().to_scalar();

       if (md->is_linear()) {
         GMM_WARNING1("Problem is linear, the test is not relevant");
       } else {
         if (md->is_complex()) {
           std::vector<complex_type> U2(nbdof2);
           std::vector<complex_type> dU2(nbdof2);
           std::vector<complex_type> DIR2(nbdof2);
           std::vector<complex_type> D1(nbdof1);
           std::vector<complex_type> D2(nbdof1);
           gmm::copy(md->complex_variable(varname2), U2);
           for (size_type i = 0; i < NN; ++i) {
             gmm::fill_random(dU2); gmm::scale(dU2, complex_type(scalef));
             gmm::add(U2, dU2);
             gmm::copy(dU2, md->set_complex_variable(varname2));
             gmm::fill_random(DIR2);
             md->assembly(getfem::model::BUILD_ALL);
             gmm::copy(gmm::sub_vector(md->complex_rhs(), I1), D2);
             gmm::mult(gmm::sub_matrix(md->complex_tangent_matrix(),
                                       I1, I2), DIR2, D1);
             gmm::add(gmm::scaled(DIR2, complex_type(EPS)), dU2);
             gmm::copy(dU2, md->set_complex_variable(varname2));
             md->assembly(getfem::model::BUILD_RHS);
             gmm::add(gmm::scaled(gmm::sub_vector(md->complex_rhs(),
                                                  I1), complex_type(-1)), D2);
             gmm::scale(D2, complex_type(1)/complex_type(EPS));
             scalar_type err = gmm::vect_dist2(D1, D2);
             cout << "Error at step " << i << " : " << err << endl;
             errmax = std::max(err, errmax);
           }
           gmm::copy(U2, md->set_complex_variable(varname2));
         } else {
           std::vector<scalar_type> U2(nbdof2);
           std::vector<scalar_type> dU2(nbdof2);
           std::vector<scalar_type> DIR2(nbdof2);
           std::vector<scalar_type> D1(nbdof1);
           std::vector<scalar_type> D2(nbdof1);
           gmm::copy(md->real_variable(varname2), U2);
           for (size_type i = 0; i < NN; ++i) {
             gmm::fill_random(dU2); gmm::scale(dU2, scalef);
             gmm::add(U2, dU2);
             gmm::copy(dU2, md->set_real_variable(varname2));
             gmm::fill_random(DIR2);
             md->assembly(getfem::model::BUILD_ALL);
             gmm::copy(gmm::sub_vector(md->real_rhs(), I1), D2);
             gmm::mult(gmm::sub_matrix(md->real_tangent_matrix(),
                                       I1, I2), DIR2, D1);
             gmm::add(gmm::scaled(DIR2, scalar_type(EPS)), dU2);
             gmm::copy(dU2, md->set_real_variable(varname2));
             md->assembly(getfem::model::BUILD_RHS);
             gmm::add(gmm::scaled(gmm::sub_vector(md->real_rhs(),
                                                  I1), scalar_type(-1)), D2);
             gmm::scale(D2, scalar_type(1)/EPS);
             scalar_type err = gmm::vect_dist2(D1, D2);
             cout << "Error at step " << i << " : " << err << endl;
             errmax = std::max(err, errmax);
           }
           gmm::copy(U2, md->set_real_variable(varname2));
         }
       }
       out.pop().from_scalar(errmax);
       );

    /*@GET expr = ('Neumann term', @str varname, @int region)
       Gives the assembly string corresponding to the Neumann term of
       the fem variable `varname` on `region`. It is deduced from the
       assembly string declared by the model bricks.
       `region` should be the index of a boundary region
       on the mesh where `varname` is defined. Care to call this function
       only after all the volumic bricks have been declared.
       Complains, if a brick
       omit to declare an assembly string. @*/
    sub_command
      ("Neumann term", 2, 2, 0, 1,
       std::string varname = in.pop().to_string();
       size_type rg =in.pop().to_integer();
       std::string expr = md->Neumann_term(varname, rg);
       out.pop().from_string(expr.c_str());
       );


    /*@GET V = ('compute isotropic linearized Von Mises or Tresca', @str varname, @str dataname_lambda, @str dataname_mu, @tmf mf_vm[, @str version])
      Compute the Von-Mises stress or the Tresca stress of a field (only
      valid for isotropic linearized elasticity in 3D). `version` should
      be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default).
      Parametrized by Lame coefficients.
      @*/
    sub_command
      ("compute isotropic linearized Von Mises or Tresca", 4, 5, 0, 1,
       std::string varname = in.pop().to_string();
       std::string dataname_lambda = in.pop().to_string();
       std::string dataname_mu = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       std::string stresca = "Von Mises";
       if (in.remaining()) stresca = in.pop().to_string();
       bool tresca = false;
       if (cmd_strmatch(stresca, "Von Mises") ||
           cmd_strmatch(stresca, "Von_Mises"))
         tresca = false;
       else if (cmd_strmatch(stresca, "Tresca"))
         tresca = true;
       else THROW_BADARG("bad option \'version\': " << stresca);

       getfem::model_real_plain_vector VMM(mf->nb_dof());
       getfem::compute_isotropic_linearized_Von_Mises_or_Tresca
       (*md, varname, dataname_lambda, dataname_mu, *mf, VMM, tresca);
       out.pop().from_dcvector(VMM);
       );

    /*@GET V = ('compute isotropic linearized Von Mises pstrain', @str varname, @str data_E, @str data_nu, @tmf mf_vm)
      Compute the Von-Mises stress  of a displacement field for isotropic
      linearized elasticity in 3D or in 2D with plane strain assumption.
      Parametrized by Young modulus and Poisson ratio.
      @*/
    sub_command
      ("compute isotropic linearized Von Mises pstrain", 4, 4, 0, 1,
       std::string varname = in.pop().to_string();
       std::string data_E = in.pop().to_string();
       std::string data_nu = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       getfem::model_real_plain_vector VMM(mf->nb_dof());
       getfem::compute_isotropic_linearized_Von_Mises_pstrain
       (*md, varname, data_E, data_nu, *mf, VMM);
       out.pop().from_dcvector(VMM);
       );

    /*@GET V = ('compute isotropic linearized Von Mises pstress', @str varname, @str data_E, @str data_nu, @tmf mf_vm)
      Compute the Von-Mises stress  of a displacement field for isotropic
      linearized elasticity in 3D or in 2D with plane stress assumption.
      Parametrized by Young modulus and Poisson ratio.
      @*/
    sub_command
      ("compute isotropic linearized Von Mises pstress", 4, 4, 0, 1,
       std::string varname = in.pop().to_string();
       std::string data_E = in.pop().to_string();
       std::string data_nu = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());

       getfem::model_real_plain_vector VMM(mf->nb_dof());
       getfem::compute_isotropic_linearized_Von_Mises_pstress
       (*md, varname, data_E, data_nu, *mf, VMM);
       out.pop().from_dcvector(VMM);
       );

    /*@GET V = ('compute Von Mises or Tresca', @str varname, @str lawname, @str dataname, @tmf mf_vm[, @str version])
      Compute on `mf_vm` the Von-Mises stress or the Tresca stress of a field
      for nonlinear elasticity in 3D. `lawname` is the constitutive law which
      could be 'SaintVenant Kirchhoff', 'Mooney Rivlin', 'neo Hookean' or
      'Ciarlet Geymonat'.
      `dataname` is a vector of parameters for the constitutive law. Its length
      depends on the law. It could be a short vector of constant values or a
      vector field described on a finite element method for variable coefficients.
      `version` should be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default).
      @*/
    sub_command
      ("compute Von Mises or Tresca", 4, 5, 0, 1,
       std::string varname = in.pop().to_string();
       std::string lawname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       std::string stresca = "Von Mises";
       if (in.remaining()) stresca = in.pop().to_string();
       bool tresca = false;
       if (cmd_strmatch(stresca, "Von Mises") ||
           cmd_strmatch(stresca, "Von_Mises"))
         tresca = false;
       else if (cmd_strmatch(stresca, "Tresca"))
         tresca = true;
       else THROW_BADARG("bad option \'version\': " << stresca);

       getfem::model_real_plain_vector VMM(mf->nb_dof());
       getfem::compute_Von_Mises_or_Tresca
       (*md, varname,
        abstract_hyperelastic_law_from_name
        (lawname, mf->linked_mesh().dim()), dataname, *mf, VMM, tresca);
       out.pop().from_dcvector(VMM);
       );

    /*@GET V = ('compute finite strain elasticity Von Mises',  @str lawname, @str varname, @str params, @tmf mf_vm[, @int region])
      Compute on `mf_vm` the Von-Mises stress of a field `varname`
      for nonlinear elasticity in 3D. `lawname` is the constitutive law which
      should be a valid name. `params` are the parameters law. It could be
      a short vector of constant values or may depend on data or variables
      of the model.
      Uses the high-level generic assembly.
      @*/
    sub_command
      ("compute finite strain elasticity Von Mises", 4, 5, 0, 1,
       std::string lawname = in.pop().to_string();
       std::string varname = in.pop().to_string();
       std::string params = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       size_type rg = size_type(-1);
       if (in.remaining()) rg = in.pop().to_integer();

       std::string ln = varname; // tolerance for the compatibility with 5.0
       filter_lawname(ln);
       if (ln.compare("saintvenant_kirchhoff") == 0 ||
           ln.compare("saint_venant_kirchhoff") == 0 ||
           ln.compare("generalized_blatz_ko") == 0 ||
           ln.compare("ciarlet_geymonat") == 0 ||
           ln.compare("incompressible_mooney_rivlin") == 0 ||
           ln.compare("compressible_mooney_rivlin") == 0 ||
           ln.compare("incompressible_neo_hookean") == 0 ||
           ln.compare("compressible_neo_hookean") == 0 ||
           ln.compare("compressible_neo_hookean_bonet") == 0 ||
           ln.compare("compressible_neo_hookean_ciarlet") == 0) {
         std::swap(lawname, varname);
       }

       getfem::model_real_plain_vector VMM(mf->nb_dof());
       getfem::compute_finite_strain_elasticity_Von_Mises
       (*md, lawname, varname, params, *mf, VMM, rg);
       out.pop().from_dcvector(VMM);
       );

    /*@GET V = ('compute second Piola Kirchhoff tensor', @str varname, @str lawname, @str dataname, @tmf mf_sigma)
      Compute on `mf_sigma` the second Piola Kirchhoff stress tensor of a field
      for nonlinear elasticity in 3D. `lawname` is the constitutive law which
      could be 'SaintVenant Kirchhoff', 'Mooney Rivlin', 'neo Hookean' or
      'Ciarlet Geymonat'.
      `dataname` is a vector of parameters for the constitutive law. Its length
      depends on the law. It could be a short vector of constant values or a
      vector field described on a finite element method for variable
      coefficients.
      @*/
    sub_command
      ("compute second Piola Kirchhoff tensor", 4, 4, 0, 1,
       std::string varname = in.pop().to_string();
       std::string lawname = in.pop().to_string();
       std::string dataname = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());

       size_type N = size_type(mf->linked_mesh().dim());
       size_type ratio = 1;
       if (mf->get_qdim() == 1) ratio = N*N;
       
       getfem::model_real_plain_vector VMM(ratio*mf->nb_dof());

       getfem::compute_sigmahathat
       (*md, varname,
        abstract_hyperelastic_law_from_name
        (lawname, mf->linked_mesh().dim()), dataname, *mf, VMM);
       out.pop().from_dcvector(VMM);
       );


    /*@GET ('elastoplasticity next iter', @tmim mim, @str varname, @str previous_dep_name, @str projname, @str datalambda, @str datamu, @str datathreshold, @str datasigma)
      Used with the old (obsolete) elastoplasticity brick to pass from an
      iteration to the next one.
      Compute and save the stress constraints sigma for the next iterations.
      'mim' is the integration method to use for the computation.
      'varname' is the main variable of the problem.
      'previous_dep_name' represents the displacement at the previous time step.
      'projname' is the type of projection to use. For the moment it could only be 'Von Mises' or 'VM'.
      'datalambda' and 'datamu' are the Lame coefficients of the material.
      'datasigma' is a vector which will contain the new stress constraints values.@*/
    sub_command
      ("elastoplasticity next iter", 8, 8, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string previous_dep = in.pop().to_string();
       std::string projname = in.pop().to_string();
       std::string datalambda = in.pop().to_string();
       std::string datamu = in.pop().to_string();
       std::string datathreshold = in.pop().to_string();
       std::string datasigma = in.pop().to_string();


       getfem::elastoplasticity_next_iter
       (*md, *mim, varname, previous_dep,
        abstract_constraints_projection_from_name(projname),
        datalambda, datamu, datathreshold, datasigma);
       );

    /*@GET ('small strain elastoplasticity next iter', @tmim mim,  @str lawname, @str unknowns_type [, @str varnames, ...] [, @str params, ...] [, @str theta = '1' [, @str dt = 'timestep']] [, @int region = -1]) 	 
      Function that allows to pass from a time step to another for the
      small strain plastic brick. The parameters have to be exactly the
      same than the one of `add_small_strain_elastoplasticity_brick`,
      so see the documentation of this function for the explanations.
      Basically, this brick computes the plastic strain
      and the plastic multiplier and stores them for the next step.
      Additionaly, it copies the computed displacement to the data
      that stores the displacement of the previous time step (typically
      'u' to 'Previous_u'). It has to be called before any use of
      `compute_small_strain_elastoplasticity_Von_Mises`.
      @*/
    sub_command
      ("small strain elastoplasticity next iter", 3, 15, 0, 0,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string lawname = in.pop().to_string();
       filter_lawname(lawname);
       size_type nb_var = 0; size_type nb_params = 0;
       if (lawname.compare("isotropic_perfect_plasticity") == 0 ||
	   lawname.compare("prandtl_reuss") == 0 ||
	   lawname.compare("plane_strain_isotropic_perfect_plasticity") == 0 ||
	   lawname.compare("plane_strain_prandtl_reuss") == 0) {
	 nb_var = nb_params = 3;
       } else if
	   (lawname.compare("isotropic_plasticity_linear_hardening") == 0 ||
	    lawname.compare("prandtl_reuss_linear_hardening") == 0 ||
	    lawname.compare("plane_strain_isotropic_plasticity_linear_hardening") == 0 ||
	    lawname.compare("plane_strain_prandtl_reuss_linear_hardening") == 0) {
	 nb_var = 4; nb_params = 5;
       } else
	 THROW_BADARG(lawname << " is not an implemented elastoplastic law");

       getfem::plasticity_unknowns_type unknowns_type(getfem::DISPLACEMENT_ONLY); 	 
       mexarg_in argin = in.pop();
       if (argin.is_string()) {
	 std::string opt = argin.to_string();
	 filter_lawname(opt);
	 if (opt.compare("displacement_only") == 0)
	   unknowns_type = getfem::DISPLACEMENT_ONLY;
	 else if (opt.compare("displacement_and_plastic_multiplier") == 0)
	   unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER;
	 else
	   THROW_BADARG("Wrong input");
       } else if (argin.is_integer())
	 unknowns_type = static_cast<getfem::plasticity_unknowns_type>
	   (argin.to_integer(0,1));

       std::vector<std::string> varnames;
       for (size_type i = 0; i < nb_var; ++i)
	 varnames.push_back(in.pop().to_string());
       
       std::vector<std::string> params;
       for (size_type i = 0; i < nb_params; ++i)
	 params.push_back(in.pop().to_string());
       
       std::string theta = "1";
       std::string dt = "timestep";
       size_type region = size_type(-1);
       for (size_type i=0; i < 3 && in.remaining(); ++i) {
	 argin = in.pop();
	 if (argin.is_string()) {
	   if (i==0)      theta = argin.to_string();
	   else if (i==1) dt = argin.to_string();
	   else           THROW_BADARG("Wrong input");
	 } else if (argin.is_integer()) {
	   region = argin.to_integer();
	   GMM_ASSERT1(!in.remaining(), "Wrong input");
	 }
       }
       params.push_back(theta);
       params.push_back(dt);
       
       getfem::small_strain_elastoplasticity_next_iter
       (*md, *mim, lawname, unknowns_type, varnames, params, region);
       workspace().set_dependence(md, mim);
       );
    
    /*@GET V = ('small strain elastoplasticity Von Mises', @tmim mim, @tmf mf_vm, @str lawname, @str unknowns_type [, @str varnames, ...] [, @str params, ...] [, @str theta = '1' [, @str dt = 'timestep']] [, @int region])
      This function computes the Von Mises stress field with respect to
      a small strain elastoplasticity term, approximated on `mf_vm`,
      and stores the result into `VM`.  All other parameters have to be
      exactly the same as for `add_small_strain_elastoplasticity_brick`.
      Remember that `small_strain_elastoplasticity_next_iter` has to be called
      before any call of this function.
      @*/
    sub_command 	 
      ("small strain elastoplasticity Von Mises", 4, 16, 0, 0,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       const getfem::mesh_fem *mf_vm = to_meshfem_object(in.pop());
       std::string lawname = in.pop().to_string();
       filter_lawname(lawname);
       size_type nb_var = 0; size_type nb_params = 0;
       if (lawname.compare("isotropic_perfect_plasticity") == 0 ||
	   lawname.compare("prandtl_reuss") == 0 ||
	   lawname.compare("plane_strain_isotropic_perfect_plasticity") == 0 ||
	   lawname.compare("plane_strain_prandtl_reuss") == 0) {
	 nb_var = nb_params = 3;
       } else if
	   (lawname.compare("isotropic_plasticity_linear_hardening") == 0 ||
	    lawname.compare("prandtl_reuss_linear_hardening") == 0 ||
	    lawname.compare("plane_strain_isotropic_plasticity_linear_hardening") == 0 ||
	    lawname.compare("plane_strain_prandtl_reuss_linear_hardening") == 0) {
	 nb_var = 4; nb_params = 5;
       } else
	 THROW_BADARG(lawname << " is not an implemented elastoplastic law");
       
       getfem::plasticity_unknowns_type unknowns_type(getfem::DISPLACEMENT_ONLY);
       mexarg_in argin = in.pop();
       if (argin.is_string()) {
	 std::string opt = argin.to_string();
	 filter_lawname(opt);
	 if (opt.compare("displacement_only") == 0)
	   unknowns_type = getfem::DISPLACEMENT_ONLY;
	 else if (opt.compare("displacement_and_plastic_multiplier") == 0)
	   unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER;
	 else
	   THROW_BADARG("Wrong input");
       } else if (argin.is_integer())
	 unknowns_type = static_cast<getfem::plasticity_unknowns_type>
	   (argin.to_integer(0,1));
       
       std::vector<std::string> varnames;
       for (size_type i = 0; i < nb_var; ++i)
	 varnames.push_back(in.pop().to_string());
      
       std::vector<std::string> params;
       for (size_type i = 0; i < nb_params; ++i)
	 params.push_back(in.pop().to_string());
       
       std::string theta = "1";
       std::string dt = "timestep";
       size_type region = size_type(-1);
       for (size_type i=0; i < 3 && in.remaining(); ++i) {
	 argin = in.pop();
	 if (argin.is_string()) {
	   if (i==0)      theta = argin.to_string();
	   else if (i==1) dt = argin.to_string();
	   else           THROW_BADARG("Wrong input");
	 } else if (argin.is_integer()) {
	   region = argin.to_integer();
	   GMM_ASSERT1(!in.remaining(), "Wrong input");
	 }
       }
       params.push_back(theta);
       params.push_back(dt);
       
       getfem::model_real_plain_vector VMM(mf_vm->nb_dof());
       getfem::compute_small_strain_elastoplasticity_Von_Mises
       (*md, *mim, lawname, unknowns_type, varnames, params, *mf_vm, VMM,
	region);
       out.pop().from_dcvector(VMM);
       );
    
    
    /*@GET V = ('compute elastoplasticity Von Mises or Tresca', @str datasigma, @tmf mf_vm[, @str version])
      Compute on `mf_vm` the Von-Mises or the Tresca stress of a field for plasticity and return it into the vector V.
      `datasigma` is a vector which contains the stress constraints values supported by the mesh.
      `version` should be  'Von_Mises' or 'Tresca' ('Von_Mises' is the default).@*/
    sub_command
      ("compute elastoplasticity Von Mises or Tresca", 2, 3, 0, 1,
       std::string datasigma = in.pop().to_string();
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       std::string stresca = "Von Mises";
       if (in.remaining()) stresca = in.pop().to_string();
       bool tresca = false;
       if (cmd_strmatch(stresca, "Von Mises") ||
           cmd_strmatch(stresca, "Von_Mises"))
         tresca = false;
       else if (cmd_strmatch(stresca, "Tresca"))
         tresca = true;
       else THROW_BADARG("bad option \'version\': " << stresca);

       getfem::model_real_plain_vector VMM(mf->nb_dof());
       getfem::compute_elastoplasticity_Von_Mises_or_Tresca
       (*md, datasigma, *mf, VMM, tresca);
       out.pop().from_dcvector(VMM);
       );

    /*@GET V = ('compute plastic part', @tmim mim, @tmf mf_pl, @str varname, @str previous_dep_name, @str projname, @str datalambda, @str datamu, @str datathreshold, @str datasigma)
      Compute on `mf_pl` the plastic part and return it into the vector V.
      `datasigma` is a vector which contains the stress constraints values supported by the mesh.@*/
     sub_command
      ("compute plastic part", 9, 9, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       const getfem::mesh_fem *mf = to_meshfem_object(in.pop());
       std::string varname = in.pop().to_string();
       std::string previous_dep = in.pop().to_string();
       std::string projname = in.pop().to_string();
       std::string datalambda = in.pop().to_string();
       std::string datamu = in.pop().to_string();
       std::string datathreshold = in.pop().to_string();
       std::string datasigma = in.pop().to_string();

       getfem::model_real_plain_vector plast(mf->nb_dof());
       getfem::compute_plastic_part
       (*md, *mim, *mf, varname, previous_dep,
        abstract_constraints_projection_from_name(projname),
        datalambda, datamu, datathreshold, datasigma, plast);
       out.pop().from_dcvector(plast);
       );


    /*@GET ('finite strain elastoplasticity next iter', @tmim mim, @str lawname, @str unknowns_type, [, @str varnames, ...] [, @str params, ...] [, @int region = -1])
      Function that allows to pass from a time step to another for the
      finite strain plastic brick. The parameters have to be exactly the
      same than the one of `add_finite_strain_elastoplasticity_brick`,
      so see the documentation of this function for the explanations.
      Basically, this brick computes the plastic strain
      and the plastic multiplier and stores them for the next step.
      For the Simo-Miehe law which is currently the only one implemented,
      this function updates the state variables defined in the last two
      entries of `varnames`, and resets the plastic multiplier field given
      as the second entry of `varnames`.
      @*/
    sub_command
      ("finite strain elastoplasticity next iter", 10, 11, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       std::string lawname = in.pop().to_string();
       filter_lawname(lawname);
       size_type nb_var = 0; size_type nb_params = 0;
       if (lawname.compare("simo_miehe") == 0 ||
           lawname.compare("eterovic_bathe") == 0) {
         nb_var = 4;
         nb_params = 3;
       } else
         THROW_BADARG(lawname << " is not an implemented finite strain"
                              << " elastoplastic law");

       getfem::plasticity_unknowns_type unknowns_type(getfem::DISPLACEMENT_ONLY);
       mexarg_in argin = in.pop();
       if (argin.is_string()) {
         std::string opt = argin.to_string();
         filter_lawname(opt);
         if (opt.compare("displacement_and_plastic_multiplier") == 0)
           unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER;
         else if (opt.compare("displacement_and_plastic_multiplier"
                              "_and_pressure") == 0)
           unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE;
         else
           THROW_BADARG("Wrong input");
       } else if (argin.is_integer()) {
         unknowns_type = static_cast<getfem::plasticity_unknowns_type>
                         (argin.to_integer());
         GMM_ASSERT1
           (unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER ||
            unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE,
            "Not valid input for unknowns_type");
       }
       if (unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE)
         nb_var += 1;

       std::vector<std::string> varnames;
       for (size_type i = 0; i < nb_var; ++i)
         varnames.push_back(in.pop().to_string());

       std::vector<std::string> params;
       for (size_type i = 0; i < nb_params; ++i)
         params.push_back(in.pop().to_string());

       size_type region(-1);
       if (in.remaining()) {
         argin = in.pop();
         if (!argin.is_integer())
           THROW_BADARG("Last optional argument must be an integer");
         region = argin.to_integer();
       }

       getfem::finite_strain_elastoplasticity_next_iter
         (*md, *mim, lawname, unknowns_type, varnames, params, region);
       );

    /*@GET V = ('compute finite strain elastoplasticity Von Mises', @tmim mim, @tmf mf_vm, @str lawname, @str unknowns_type, [, @str varnames, ...] [, @str params, ...] [, @int region = -1])
      Compute on `mf_vm` the Von-Mises or the Tresca stress of a field for plasticity and return it into the vector V.
      The first input parameters ar as in the function 'finite strain elastoplasticity next iter'.
      @*/
    sub_command
      ("compute finite strain elastoplasticity Von Mises", 11, 12, 0, 1,
       getfem::mesh_im *mim = to_meshim_object(in.pop());
       const getfem::mesh_fem *mf_vm = to_meshfem_object(in.pop());

       std::string lawname = in.pop().to_string();
       filter_lawname(lawname);
       size_type nb_var = 0; size_type nb_params = 0;
       if (lawname.compare("simo_miehe") == 0 ||
           lawname.compare("eterovic_bathe") == 0) {
         nb_var = 4;
         nb_params = 3;
       } else
         THROW_BADARG(lawname << " is not an implemented finite strain"
                              << " elastoplastic law");

       getfem::plasticity_unknowns_type unknowns_type(getfem::DISPLACEMENT_ONLY);
       mexarg_in argin = in.pop();
       if (argin.is_string()) {
         std::string opt = argin.to_string();
         filter_lawname(opt);
         if (opt.compare("displacement_and_plastic_multiplier") == 0)
           unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER;
         else if (opt.compare("displacement_and_plastic_multiplier"
                              "_and_pressure") == 0)
           unknowns_type = getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE;
         else
           THROW_BADARG("Wrong input");
       } else if (argin.is_integer()) {
         unknowns_type = static_cast<getfem::plasticity_unknowns_type>
                         (argin.to_integer());
         GMM_ASSERT1
           (unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER ||
            unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE,
            "Not valid input for unknowns_type");
       }
       if (unknowns_type == getfem::DISPLACEMENT_AND_PLASTIC_MULTIPLIER_AND_PRESSURE)
         nb_var += 1;

       std::vector<std::string> varnames;
       for (size_type i = 0; i < nb_var; ++i)
         varnames.push_back(in.pop().to_string());

       std::vector<std::string> params;
       for (size_type i = 0; i < nb_params; ++i)
         params.push_back(in.pop().to_string());

       size_type region(-1);
       if (in.remaining()) {
         argin = in.pop();
         if (!argin.is_integer())
           THROW_BADARG("Last optional argument must be an integer");
         region = argin.to_integer();
       }

       getfem::model_real_plain_vector VMM(mf_vm->nb_dof());
       getfem::compute_finite_strain_elastoplasticity_Von_Mises
         (*md, *mim, lawname, unknowns_type, varnames, params, *mf_vm, VMM);
       out.pop().from_dcvector(VMM);
       );


     /*@GET V = ('sliding data group name of large sliding contact brick', @int indbrick)
      Gives the name of the group of variables corresponding to the
      sliding data for an existing large sliding contact brick.@*/
      sub_command
      ("sliding data group name of large sliding contact brick", 1, 1, 0, 1,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string name
       = sliding_data_group_name_of_large_sliding_contact_brick(*md, ind);
       out.pop().from_string(name.c_str());
       );

     /*@GET V = ('displacement group name of large sliding contact brick', @int indbrick)
      Gives the name of the group of variables corresponding to the
      sliding data for an existing large sliding contact brick.@*/
      sub_command
      ("displacement group name of large sliding contact brick", 1, 1, 0, 1,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string name
       = displacement_group_name_of_large_sliding_contact_brick(*md, ind);
       out.pop().from_string(name.c_str());
       );

     /*@GET V = ('transformation name of large sliding contact brick', @int indbrick)
      Gives the name of the group of variables corresponding to the
      sliding data for an existing large sliding contact brick.@*/
      sub_command
      ("transformation name of large sliding contact brick", 1, 1, 0, 1,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string name
       = transformation_name_of_large_sliding_contact_brick(*md, ind);
       out.pop().from_string(name.c_str());
       );
  /*@GET V = ('sliding data group name of Nitsche large sliding contact brick', @int indbrick)
      Gives the name of the group of variables corresponding to the
      sliding data for an existing large sliding contact brick.@*/
      sub_command
      ("sliding data group name of Nitsche large sliding contact brick", 1, 1, 0, 1,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string name
       = sliding_data_group_name_of_Nitsche_large_sliding_contact_brick(*md, ind);
       out.pop().from_string(name.c_str());
       );

     /*@GET V = ('displacement group name of Nitsche large sliding contact brick', @int indbrick)
      Gives the name of the group of variables corresponding to the
      sliding data for an existing large sliding contact brick.@*/
      sub_command
      ("displacement group name of Nitsche large sliding contact brick", 1, 1, 0, 1,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string name
       = displacement_group_name_of_Nitsche_large_sliding_contact_brick(*md, ind);
       out.pop().from_string(name.c_str());
       );

     /*@GET V = ('transformation name of Nitsche large sliding contact brick', @int indbrick)
      Gives the name of the group of variables corresponding to the
      sliding data for an existing large sliding contact brick.@*/
      sub_command
      ("transformation name of Nitsche large sliding contact brick", 1, 1, 0, 1,
       size_type ind = in.pop().to_integer() - config::base_index();
       std::string name
       = transformation_name_of_Nitsche_large_sliding_contact_brick(*md, ind);
       out.pop().from_string(name.c_str());
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
       infomsg() << "gfModel object with " << md->nb_dof()
       << " degrees of freedom\n";
       );

  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::model *md  = to_model_object(m_in.pop());
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
