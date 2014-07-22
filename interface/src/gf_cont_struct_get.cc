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
  
    
    /*@FUNC h = ('init step size')
      Return an initial step size for continuation.@*/
    sub_command
      ("init step size", 0, 0, 0, 1,
       
       out.pop().from_scalar(ps->h_init());
       );
  

    /*@FUNC E = ('init Moore-Penrose continuation', @vec solution, @scalar parameter, @scalar init_dir)
      Initialise the Moore-Penrose continuation: Return a unit tangent to
      the solution curve at the point given by `solution` and `parameter`,
      and an initial step size for the continuation. Orientation of the
      computed tangent with respect to the parameter is determined by the
      sign of `init_dir`.@*/
    sub_command
      ("init Moore-Penrose continuation", 3, 3, 0, 3,

       size_type nbdof = ps->linked_model().nb_dof();
       darray x = in.pop().to_darray();
       std::vector<double> xx(nbdof); gmm::copy(x, xx);
       scalar_type gamma = in.pop().to_scalar();
       std::vector<double> tt_x(nbdof);
       scalar_type t_gamma = in.pop().to_scalar();
       scalar_type h;

       getfem::init_Moore_Penrose_continuation(*ps, xx, gamma,
					       tt_x, t_gamma, h);
       out.pop().from_dcvector(tt_x);
       out.pop().from_scalar(t_gamma);
       out.pop().from_scalar(h);
       );


    /*@FUNC E = ('Moore-Penrose continuation', @vec solution, @scalar parameter, @vec tangent_sol, @scalar tangent_par, @scalar h)
      Compute one step of the Moore-Penrose continuation: Take the point
      given by `solution` and `parameter`, the tangent given by `tangent_sol`
      and `tangent_par`, and the step size `h`. Return a new point on the
      solution curve, the corresponding tangent and a step size for the next
      step. If the returned step size equals zero, the continuation has
      failed. Optionally, return the type of any detected singular point.
      NOTE: The new point need not to be saved in the model in the end!@*/
    sub_command
      ("Moore-Penrose continuation", 5, 5, 0, 6,

       size_type nbdof = ps->linked_model().nb_dof();
       darray x = in.pop().to_darray();
       std::vector<double> xx(nbdof); gmm::copy(x, xx);
       scalar_type gamma = in.pop().to_scalar();
       darray t_x = in.pop().to_darray();
       std::vector<double> tt_x(nbdof); gmm::copy(t_x, tt_x);
       scalar_type t_gamma = in.pop().to_scalar();
       scalar_type h = in.pop().to_scalar();

       getfem::Moore_Penrose_continuation(*ps, xx, gamma, tt_x, t_gamma, h);
       out.pop().from_dcvector(xx);
       out.pop().from_scalar(gamma);
       out.pop().from_dcvector(tt_x);
       out.pop().from_scalar(t_gamma);
       out.pop().from_scalar(h);
       if (out.remaining())
	 out.pop().from_string(ps->get_sing_label().c_str());
       );


    /*@GET t = ('bifurcation test function')
      Return the last value of the bifurcation test function and eventaully
      the whole calculated graph when passing between different sub-domains
      of differentiability.@*/
    sub_command
      ("bifurcation test function", 0, 0, 0, 3,
       out.pop().from_scalar(ps->get_tau_bp_2());
       if (out.remaining()) out.pop().from_dcvector(ps->get_alpha_hist());
       if (out.remaining()) out.pop().from_dcvector(ps->get_tau_bp_hist());
       );


    /*@GET @CELL{X, gamma, T_X, T_gamma} = ('sing_data')
      Return a singular point (`X`, `gamma`) encountered in the last
      continuation step (if any) and a couple of arrays (`T_X`, `T_gamma`) of
      tangents to all located solution branches, which emanate from there.@*/
    sub_command
      ("sing_data", 0, 0, 0, 4,
       out.pop().from_dcvector(ps->get_x_sing());
       out.pop().from_scalar(ps->get_gamma_sing());
       out.pop().from_vector_container(ps->get_t_x_sing());
       out.pop().from_dcvector(ps->get_t_gamma_sing());
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
