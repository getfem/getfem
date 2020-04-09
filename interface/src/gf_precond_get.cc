/*===========================================================================

 Copyright (C) 2006-2020 Yves Renard, Julien Pommier.

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

#include <getfemint_precond.h>

using namespace getfemint;

template <typename T> static void
mult_or_tmult(gprecond<T>& precond, mexargs_in& in, mexargs_out& out,
	      bool tmult) {
  garray<T> v = in.pop().to_garray(T());
  garray<T> w = out.pop().create_array(v.getm(), v.getn(), T());
  gmm::mult_or_transposed_mult(precond, v, w, tmult);
}


/*@GFDOC
  General function for querying information about @tprecond objects.
@*/




// Object for the declaration of a new sub-command.

struct sub_gf_precond_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   gprecond_base *precond) = 0;
};

typedef std::shared_ptr<sub_gf_precond_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_precond_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       gprecond_base *precond)				\
      { dummy_func(in); dummy_func(out); dummy_func(precond); code }	\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           



void gf_precond_get(getfemint::mexargs_in& m_in,
		    getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {


    /*@GET ('mult', @vec V)
    Apply the preconditioner to the supplied vector.@*/
    sub_command
      ("mult", 1, 1, 0, 1,
       gprecond<scalar_type> *rprecond
       = dynamic_cast<gprecond<scalar_type> *>(precond);
       gprecond<complex_type> *cprecond
       = dynamic_cast<gprecond<complex_type> *>(precond);
       if (rprecond) mult_or_tmult(*rprecond, in, out, false);
       else if (cprecond) mult_or_tmult(*cprecond, in, out, false);
       else THROW_INTERNAL_ERROR;
       );


    /*@GET ('tmult', @vec V)
      Apply the transposed preconditioner to the supplied vector.@*/
    sub_command
      ("tmult", 1, 1, 0, 1,
       gprecond<scalar_type> *rprecond
       = dynamic_cast<gprecond<scalar_type> *>(precond);
       gprecond<complex_type> *cprecond
       = dynamic_cast<gprecond<complex_type> *>(precond);
       if (rprecond) mult_or_tmult(*rprecond, in, out, false);
       else if (cprecond) mult_or_tmult(*cprecond, in, out, false);
       else THROW_INTERNAL_ERROR;
       );


    /*@GET ('type')
      Return a string describing the type of the preconditioner ('ilu', 'ildlt',..).@*/
    sub_command
      ("type", 0, 0, 0, 1,
       out.pop().from_string(precond->name());
       );


    /*@GET ('size')
      Return the dimensions of the preconditioner.@*/
    sub_command
      ("size", 0, 0, 0, 1,
       iarray sz = out.pop().create_iarray_h(2);
       sz[0] = int(precond->nrows());
       sz[1] = int(precond->ncols());
       );


    /*@GET ('is_complex')
      Return 1 if the preconditioner stores complex values.@*/
    sub_command
      ("is_complex", 0, 0, 0, 1,
       gprecond<scalar_type> *rprecond
       = dynamic_cast<gprecond<scalar_type> *>(precond);
       out.pop().from_integer(rprecond == 0);
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tprecond.

      This can be used to perform comparisons between two
      different @tprecond objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false, "Sorry, function to be done");
       // std::string s = ...;
       // out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
      displays a short summary for a @tprecond object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       gprecond<scalar_type> *rprecond
       = dynamic_cast<gprecond<scalar_type> *>(precond);
       infomsg() << "gfPrecond object with "
       << precond->nrows() << "x"
       << precond->ncols() << " "
       << ((rprecond == 0) ? "COMPLEX" : "REAL")
       << " " << precond->name() << " ["
       << precond->memsize() << " bytes]";
       );

  }


  if (m_in.narg() < 1)  THROW_BADARG( "Wrong number of input arguments");
 
  gprecond_base *precond = to_precond_object(m_in.pop());
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, precond);
  }
  else bad_cmd(init_cmd);

}
