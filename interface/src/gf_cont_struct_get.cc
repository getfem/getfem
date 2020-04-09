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

#include <getfemint_misc.h>
#include <getfemint_workspace.h>
#include <getfemint.h>
#include <getfem/getfem_continuation.h>

using namespace getfemint;

// Object for the declaration of a new sub-command.

struct sub_gf_cont_struct_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
                   getfemint::mexargs_out& out,
                   getfem::cont_struct_getfem_model *ps) = 0;
};

typedef std::shared_ptr<sub_gf_cont_struct_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_cont_struct_get {                           \
      virtual void run(getfemint::mexargs_in& in,                           \
                       getfemint::mexargs_out& out,                         \
                       getfem::cont_struct_getfem_model *ps)                \
      { dummy_func(in); dummy_func(out); dummy_func(ps); code }             \
    };                                                                      \
    psub_command psubc = std::make_shared<subc>();			    \
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;             \
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;         \
    subc_tab[cmd_normalize(name)] = psubc;                                  \
  }


/*@GFDOC
  General function for querying information about cont_struct objects and for
  applying them to numerical continuation.
@*/

void gf_cont_struct_get(getfemint::mexargs_in& m_in,
                        getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@FUNC h = ('init step size')
      Return an initial step size for continuation.@*/
    sub_command
      ("init step size", 0, 0, 1, 1,
       out.pop().from_scalar(ps->h_init());
       );

    /*@FUNC h = ('min step size')
      Return the minimum step size for continuation.@*/
    sub_command
      ("min step size", 0, 0, 0, 1,
       out.pop().from_scalar(ps->h_min());
       );

    /*@FUNC h = ('max step size')
      Return the maximum step size for continuation.@*/
    sub_command
      ("max step size", 0, 0, 0, 1,
       out.pop().from_scalar(ps->h_max());
       );

    /*@FUNC h = ('step size decrement')
      Return the decrement ratio of the step size for continuation.@*/
    sub_command
      ("step size decrement", 0, 0, 0, 1,
       out.pop().from_scalar(ps->h_dec());
       );

    /*@FUNC h = ('step size increment')
      Return the increment ratio of the step size for continuation.@*/
    sub_command
      ("step size increment", 0, 0, 0, 1,
       out.pop().from_scalar(ps->h_inc());
       );

    /*@FUNC [@vec tangent_sol, @scalar tangent_par] = ('compute tangent', @vec solution, @scalar parameter, @vec tangent_sol, @scalar tangent_par)
      Compute and return an updated tangent.@*/
    sub_command
      ("compute tangent", 4, 4, 2, 2,
       size_type nbdof = ps->linked_model().nb_dof();
       darray x0 = in.pop().to_darray();
       scalar_type gamma = in.pop().to_scalar();
       darray tx0 = in.pop().to_darray();
       std::vector<double> x(nbdof); gmm::copy(x0, x);
       std::vector<double> tx(nbdof); gmm::copy(tx0, tx);
       scalar_type tgamma = in.pop().to_scalar();

       ps->compute_tangent(x, gamma, tx, tgamma);

       out.pop().from_dcvector(tx);
       out.pop().from_scalar(tgamma);
       );

    /*@FUNC E = ('init Moore-Penrose continuation', @vec solution, @scalar parameter, @scalar init_dir)
      Initialise the Moore-Penrose continuation: Return a unit tangent to
      the solution curve at the point given by `solution` and `parameter`,
      and an initial step size for the continuation. Orientation of the
      computed tangent with respect to the parameter is determined by the
      sign of `init_dir`.@*/
    sub_command
      ("init Moore-Penrose continuation", 3, 3, 3, 3,

       size_type nbdof = ps->linked_model().nb_dof();
       darray x_ = in.pop().to_darray();
       std::vector<double> x(nbdof); gmm::copy(x_, x);
       scalar_type gamma = in.pop().to_scalar();
       std::vector<double> tx(nbdof);
       scalar_type tgamma = in.pop().to_scalar();
       scalar_type h;

       ps->init_Moore_Penrose_continuation(x, gamma, tx, tgamma, h);
       out.pop().from_dcvector(tx);
       out.pop().from_scalar(tgamma);
       out.pop().from_scalar(h);
       );


    /*@FUNC E = ('Moore-Penrose continuation', @vec solution, @scalar parameter, @vec tangent_sol, @scalar tangent_par, @scalar h)
      Compute one step of the Moore-Penrose continuation: Take the point
      given by `solution` and `parameter`, the tangent given by `tangent_sol`
      and `tangent_par`, and the step size `h`. Return a new point on the
      solution curve, the corresponding tangent, a step size for the next
      step and optionally the current step size. If the returned step
      size equals zero, the continuation has failed. Optionally, return
      the type of any detected singular point.
      NOTE: The new point need not to be saved in the model in the end!@*/
    sub_command
      ("Moore-Penrose continuation", 5, 5, 5, 7,

       size_type nbdof = ps->linked_model().nb_dof();
       darray x_ = in.pop().to_darray();
       std::vector<double> x(nbdof); gmm::copy(x_, x);
       scalar_type gamma = in.pop().to_scalar();
       darray tx_ = in.pop().to_darray();
       std::vector<double> tx(nbdof); gmm::copy(tx_, tx);
       scalar_type tgamma = in.pop().to_scalar();
       scalar_type h = in.pop().to_scalar();
       scalar_type h0(0);

       ps->Moore_Penrose_continuation(x, gamma, tx, tgamma, h, h0);
       out.pop().from_dcvector(x);
       out.pop().from_scalar(gamma);
       out.pop().from_dcvector(tx);
       out.pop().from_scalar(tgamma);
       out.pop().from_scalar(h);
       if (out.remaining())
         out.pop().from_scalar(h0);
       if (out.remaining())
         out.pop().from_string(ps->get_sing_label().c_str());
       );


    /*@GET t = ('non-smooth bifurcation test', @vec solution1, @scalar parameter1, @vec tangent_sol1, @scalar tangent_par1, @vec solution2, @scalar parameter2, @vec tangent_sol2, @scalar tangent_par2)
      Test for a non-smooth bifurcation point between the point given by
      `solution1` and `parameter1` with the tangent given by `tangent_sol1`
      and `tangent_par1` and the point given by `solution2` and `parameter2`
      with the tangent given by `tangent_sol2` and `tangent_par2`.@*/
    sub_command
      ("non-smooth bifurcation test", 8, 8, 1, 1,

       size_type nbdof = ps->linked_model().nb_dof();
       darray x1_ = in.pop().to_darray();
       std::vector<double> x1(nbdof); gmm::copy(x1_, x1);
       scalar_type gamma1 = in.pop().to_scalar();
       darray tx1_ = in.pop().to_darray();
       std::vector<double> tx1(nbdof); gmm::copy(tx1_, tx1);
       scalar_type tgamma1 = in.pop().to_scalar();
       darray x2_ = in.pop().to_darray();
       std::vector<double> x2(nbdof); gmm::copy(x2_, x2);
       scalar_type gamma2 = in.pop().to_scalar();
       darray tx2_ = in.pop().to_darray();
       std::vector<double> tx2(nbdof); gmm::copy(tx2_, tx2);
       scalar_type tgamma2 = in.pop().to_scalar();
       ps->init_border(nbdof);
       ps->clear_tau_bp_currentstep();
       out.pop().from_integer(int(ps->test_nonsmooth_bifurcation
                                  (x1, gamma1, tx1, tgamma1,
                                   x2, gamma2, tx2, tgamma2)));
       );


    /*@GET t = ('bifurcation test function')
      Return the last value of the bifurcation test function and eventually
      the whole calculated graph when passing between different sub-domains
      of differentiability.@*/
    sub_command
      ("bifurcation test function", 0, 0, 1, 3,

       out.pop().from_scalar(ps->get_tau_bp_2());
       if (out.remaining()) out.pop().from_dcvector(ps->get_alpha_hist());
       if (out.remaining()) out.pop().from_dcvector(ps->get_tau_bp_hist());
       );


    /* @GET ('non-smooth branching', @vec solution, @scalar parameter, @vec tangent_sol, @scalar tangent_par)
       Approximate a non-smooth point close to the point given by `solution`
       and `parameter` and locate one-sided smooth solution branches
       emanating from there. Save the approximation of the non-smooth point
       as a singular point into the @tcs object together with the array of
       the tangents to the located solution branches that direct away from
       the non-smooth point. It is supposed that the point given by
       `solution` and `parameter` is a point on a smooth solution branch
       within the distance equal to the minimum step size from the end point
       of this branch, and the corresponding tangent given by `tangent_sol`
       and `tangent_par` is directed towards the end point.@*/
    sub_command
      ("non-smooth branching", 4, 4, 0, 0,

       size_type nbdof = ps->linked_model().nb_dof();
       darray x_ = in.pop().to_darray();
       std::vector<double> x(nbdof); gmm::copy(x_, x);
       scalar_type gamma = in.pop().to_scalar();
       darray tx_ = in.pop().to_darray();
       std::vector<double> tx(nbdof); gmm::copy(tx_, tx);
       scalar_type tgamma = in.pop().to_scalar();

       ps->clear_sing_data();
       ps->treat_nonsmooth_point(x, gamma, tx, tgamma, 0);
       );


    /*@GET @CELL{X, gamma, T_X, T_gamma} = ('sing_data')
      Return a singular point (`X`, `gamma`) stored in the @tcs object and a
      couple of arrays (`T_X`, `T_gamma`) of tangents to all located solution
      branches that emanate from there.@*/
    sub_command
      ("sing_data", 0, 0, 0, 4,

       out.pop().from_dcvector(ps->get_x_sing());
       out.pop().from_scalar(ps->get_gamma_sing());
       out.pop().from_vector_container(ps->get_tx_sing());
       out.pop().from_dcvector(ps->get_tgamma_sing());
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tcs.

      This can be used for performing comparisons between two
      different @tcs objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false, "Sorry, function to be done");
       // std::string s = ...;
       // out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
      Display a short summary for a @tcs object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfContStruct object\n";
       );

  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::cont_struct_getfem_model *ps = to_cont_struct_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
              it->second->arg_in_max, it->second->arg_out_min,
              it->second->arg_out_max);
    it->second->run(m_in, m_out, ps);
  }
  else bad_cmd(init_cmd);

}
