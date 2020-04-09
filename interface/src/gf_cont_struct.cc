/*===========================================================================

 Copyright (C) 2012-2020 Tomas Ligursky, Yves Renard, Konstantinos Poulios.

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

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <getfem/getfem_continuation.h>


using namespace getfemint;

/*@GFDOC
  This object serves for storing parameters and data used in numerical
  continuation of solution branches of models (for more details about
  continuation see the GetFEM user documentation).
@*/

void gf_cont_struct(getfemint::mexargs_in& in, getfemint::mexargs_out& out) {
  if (check_cmd("ContStruct", "ContStruct", in, out, 3, 43, 0, 1)) {

    /*@INIT S = ('.init', @tmodel md, @str dataname_parameter[,@str dataname_init, @str dataname_final, @str dataname_current], @scalar sc_fac[, ...])
    The variable `dataname_parameter` should parametrise the model given by
    `md`. If the parameterization is done via a vector datum, `dataname_init`
    and `dataname_final` should store two given values of this datum
    determining the parameterization, and `dataname_current` serves for actual
    values of this datum. `sc_fac` is a scale factor involved in the weighted
    norm used in the continuation.

    Additional options:

    - 'lsolver', @str SOLVER_NAME
       name of the solver to be used for the incorporated linear systems
       (the default value is 'auto', which lets getfem choose itself);
       possible values are 'superlu', 'mumps' (if supported), 'cg/ildlt',
       'gmres/ilu' and 'gmres/ilut';
    - 'h_init', @scalar HIN
       initial step size (the default value is 1e-2);
    - 'h_max', @scalar HMAX
       maximum step size (the default value is 1e-1);
    - 'h_min', @scalar HMIN
       minimum step size (the default value is 1e-5);
    - 'h_inc', @scalar HINC
       factor for enlarging the step size (the default value is 1.3);
    - 'h_dec', @scalar HDEC
       factor for diminishing the step size (the default value is 0.5);
    - 'max_iter', @int MIT
       maximum number of iterations allowed in the correction (the default
       value is 10);
    - 'thr_iter', @int TIT
       threshold number of iterations of the correction for enlarging the
       step size (the default value is 4);
    - 'max_res', @scalar RES
       target residual value of a new point on the solution curve (the
       default value is 1e-6);
    - 'max_diff', @scalar DIFF
       determines a convergence criterion for two consecutive points (the
       default value is 1e-6);
    - 'min_cos', @scalar MCOS
       minimal value of the cosine of the angle between tangents to the
       solution curve at an old point and a new one (the default value is
       0.9);
    - 'max_res_solve', @scalar RES_SOLVE
       target residual value for the linear systems to be solved (the
       default value is 1e-8);
    - 'singularities', @int SING
       activates tools for detection and treatment of singular points (1 for
       limit points, 2 for bifurcation points and points requiring special
       branching techniques);
    - 'non-smooth'
       determines that some special methods for non-smooth problems can be
       used;
    - 'delta_max', @scalar DMAX
       maximum size of division for evaluating the test function on the
       convex combination of two augmented Jacobians that belong to different
       smooth pieces (the default value is 0.005);
    - 'delta_min', @scalar DMIN
       minimum size of division for evaluating the test function on the
       convex combination (the default value is 0.00012);
    - 'thr_var', @scalar TVAR
       threshold variation for refining the division (the default value is
       0.02);
    - 'nb_dir', @int NDIR
       total number of the linear combinations of one couple of reference
       vectors when searching for new tangent predictions during location of
       new one-sided branches (the default value is 40);
    - 'nb_span', @int NSPAN
       total number of the couples of the reference vectors forming the
       linear combinations (the default value is 1);
    - 'noisy' or 'very_noisy'
       determines how detailed information has to be displayed during the
       continuation process (residual values etc.).@*/

    getfem::model *md = to_model_object(in.pop());
    std::string dataname_parameter = in.pop().to_string(),
      dataname_init, dataname_final, dataname_current;
    if (in.front().is_string()) {
      dataname_init = in.pop().to_string();
      dataname_final = in.pop().to_string();
      dataname_current = in.pop().to_string();
    }
    scalar_type scfac = in.pop().to_scalar();
    
    std::string lsolver("auto"), varname("");
    scalar_type h_init(1e-2), h_max(1e-1), h_min(1e-5), h_inc(1.3), h_dec(0.5),
      maxres(1e-6), maxdiff(1e-6), mincos(0.9), maxres_solve(1e-8),
      delta_max(0.005), delta_min(0.00012), thrvar(0.02);
    size_type maxit(10), thrit(4), nbdir(40), nbspan(1);
    int noisy(0), singularities(0);
    bool nonsmooth(false);
    
    while (in.remaining() && in.front().is_string()) {
      std::string opt = in.pop().to_string();
      if (cmd_strmatch(opt, "lsolver"))  {
	if (in.remaining()) lsolver = in.pop().to_string();
	else THROW_BADARG("missing name for " << opt);
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
      } else if (cmd_strmatch(opt, "min_cos")) {
	if (in.remaining()) mincos = in.pop().to_scalar();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "max_res_solve")) {
	if (in.remaining()) maxres_solve = in.pop().to_scalar();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "delta_max")) {
	if (in.remaining()) delta_max = in.pop().to_scalar();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "delta_min")) {
	if (in.remaining()) delta_min = in.pop().to_scalar();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "thr_var")) {
	if (in.remaining()) thrvar = in.pop().to_scalar();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "nb_dir")) {
	if (in.remaining()) nbdir = in.pop().to_integer();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "nb_span")) {
	if (in.remaining()) nbspan = in.pop().to_integer();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "singularities")) {
	if (in.remaining()) singularities = in.pop().to_integer();
	else THROW_BADARG("missing value for " << opt);
      } else if (cmd_strmatch(opt, "variable_name")) {
	if (in.remaining()) varname = in.pop().to_string();
	else THROW_BADARG("missing name for " << opt);
      } else if (cmd_strmatch(opt, "non-smooth")) nonsmooth = true;
      else if (cmd_strmatch(opt, "noisy")) noisy = 1;
      else if (cmd_strmatch(opt, "very noisy") ||
	       cmd_strmatch(opt, "very_noisy")) noisy = 2;
      else THROW_BADARG("bad option: " << opt);
    }
    
    auto shp = std::make_shared<getfem::cont_struct_getfem_model>
      (*md, dataname_parameter, scfac,
       getfem::rselect_linear_solver(*md, lsolver),
       h_init, h_max, h_min, h_inc, h_dec, maxit, thrit, maxres, maxdiff,
       mincos, maxres_solve, noisy, singularities, nonsmooth,
       delta_max, delta_min, thrvar, nbdir, nbspan);
    
    
    
    if (!dataname_current.empty())
      shp->set_parametrised_data_names(dataname_init, dataname_final,
				      dataname_current);
    if (!varname.empty())
      shp->set_interval_from_variable_name(varname);
    
    id_type id = store_cont_struct_object(shp);
    workspace().set_dependence(id, (const void *)(md));
    out.pop().from_object_id(id, CONT_STRUCT_CLASS_ID);
  }
  
}
