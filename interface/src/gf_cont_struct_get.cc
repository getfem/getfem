/*===========================================================================
 
 Copyright (C) 2012-2012 Tomas Ligursky, Yves Renard.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
#include <getfemint_cont_struct.h>

using namespace getfemint;

// Object for the declaration of a new sub-command.

struct sub_gf_cont_struct_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::cont_struct_getfem_model *ps) = 0;
};

typedef boost::intrusive_ptr<sub_gf_cont_struct_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_cont_struct_get {			\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::cont_struct_getfem_model *ps)		\
      { dummy_func(in); dummy_func(out); dummy_func(ps); code }		\
    };									\
    psub_command psubc = new subc;					\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
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
  

    /*@FUNC t = ('init test function', @vec tangent, @scalar tangent_parameter)
      Initialise the border of the bordered system that serves for
      calculating the test function. Return the value of the test function
      for the current point (i.e., the solution and the parameter saved in
      the corresponding model object) and the tangent given by `tangent` and
      `tangent_parameter`.@*/
    sub_command
      ("init test function", 2, 2, 0, 1,

       size_type nbdof = ps->linked_model().nb_dof();
       std::vector<double> yy(nbdof); ps->linked_model().from_variables(yy);
       const getfem::model_real_plain_vector &GAMMA =
       ps->linked_model().real_variable(ps->parameter_name());
       GMM_ASSERT1(gmm::vect_size(GAMMA) == 1,
                   "The continuation parameter should be a real scalar!");
       scalar_type gamma = GAMMA[0];
       darray t_y = in.pop().to_darray();
       std::vector<double> tt_y(nbdof); gmm::copy(t_y, tt_y);
       scalar_type t_gamma = in.pop().to_scalar();

       getfem::init_test_function(*ps, yy, gamma, tt_y, t_gamma);
       out.pop().from_scalar(ps->get_tau2());
       );
  

    /*@FUNC E = ('init Moore-Penrose continuation', @scalar init_dir)
      Initialise the Moore-Penrose continuation: Return a unit tangent
      corresponding to the solution branch at the solution and the
      value of the parameter saved in the corresponding model object,
      and an initial step size for the continuation. Direction of the
      computed tangent with respect to the parameter is determined by the
      sign of `init_dir`.@*/
    sub_command
      ("init Moore-Penrose continuation", 1, 1, 0, 3,

       size_type nbdof = ps->linked_model().nb_dof();
       std::vector<double> yy(nbdof); ps->linked_model().from_variables(yy);
       const getfem::model_real_plain_vector &GAMMA
       = ps->linked_model().real_variable(ps->parameter_name());
       GMM_ASSERT1(gmm::vect_size(GAMMA) == 1,
                   "The continuation parameter should be a real scalar!");
       scalar_type gamma = GAMMA[0];
       std::vector<double> tt_y(nbdof);
       scalar_type t_gamma = in.pop().to_scalar();
       scalar_type h;

       getfem::init_Moore_Penrose_continuation(*ps, yy, gamma,
					       tt_y, t_gamma, h);
       out.pop().from_dcvector(tt_y);
       out.pop().from_scalar(t_gamma);
       out.pop().from_scalar(h);
       );


    /*@FUNC E = ('Moore-Penrose continuation', @vec tangent, @scalar tangent_parameter, @scalar h)
      Compute one step of the Moore-Penrose continuation: Take the solution
      and the value of the parameter saved in the corresponding model object,
      the tangent given by `tangent` and `tangent_parameter`, and the step
      size `h`, save a new point on the solution curve into the model object,
      and return a new tangent and a step size for the next step. If the
      returned step size equals zero, the continuation has failed.
      Eventually, return a message.@*/
    sub_command
      ("Moore-Penrose continuation", 3, 3, 0, 4,

       size_type nbdof = ps->linked_model().nb_dof();
       std::vector<double> yy(nbdof); ps->linked_model().from_variables(yy);
       const getfem::model_real_plain_vector &GAMMA
       = ps->linked_model().real_variable(ps->parameter_name());
       GMM_ASSERT1(gmm::vect_size(GAMMA) == 1,
                   "The continuation parameter should be a real scalar!");
       scalar_type gamma = GAMMA[0];
       darray t_y = in.pop().to_darray();
       std::vector<double> tt_y(nbdof); gmm::copy(t_y, tt_y);
       scalar_type t_gamma = in.pop().to_scalar();
       scalar_type h = in.pop().to_scalar();

       getfem::Moore_Penrose_continuation(*ps, yy, gamma, tt_y, t_gamma, h);
       out.pop().from_dcvector(tt_y);
       out.pop().from_scalar(t_gamma);
       out.pop().from_scalar(h);
       if (out.remaining()) out.pop().from_string(ps->get_message().c_str());
       );


    /*@GET t = ('test function')
      Return the last value of the test function and eventaully all the
      points and the corresponding values calculated when passing through a
      boundary between different smooth pieces.@*/
    sub_command
      ("test function", 0, 0, 0, 3,
       out.pop().from_scalar(ps->get_tau2());
       if (out.remaining()) out.pop().from_dcvector(ps->get_alpha_hist());
       if (out.remaining()) out.pop().from_dcvector(ps->get_tau_hist());
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tcs.

      This can be used to perform comparisons between two
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

  getfem::cont_struct_getfem_model *ps = m_in.pop().to_cont_struct();
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
