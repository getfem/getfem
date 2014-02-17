/*===========================================================================
 
 Copyright (C) 2013-2013 Yves Renard.
 
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
#include <getfemint_multi_contact_frame.h>

using namespace getfemint;

// Object for the declaration of a new sub-command.

struct sub_gf_mcf_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::multi_contact_frame *ps) = 0;
};

typedef boost::intrusive_ptr<sub_gf_mcf_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mcf_get {       			\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::multi_contact_frame *ps)                 \
      { dummy_func(in); dummy_func(out); dummy_func(ps); code }		\
    };									\
    psub_command psubc = new subc;					\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }       


/*@GFDOC
  General function for querying information about multi contact frame objects.
@*/

void gf_multi_contact_frame_get(getfemint::mexargs_in& m_in,
                                getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {
  

    
    /*@GET s = ('compute pairs')
      Compute the contact pairs
      @*/
    sub_command
      ("compute pairs", 0, 0, 0, 0,
       ps->compute_contact_pairs();
       );

    /*@GET s = ('slave points')
      Get the slave points computed.
      @*/
    sub_command
      ("slave points", 0, 0, 0, 1,

       size_type nbp = ps->ct_pairs().size();
       size_type N = ps->dim();
       darray w1 = out.pop().create_darray((unsigned int)(N),
                                           (unsigned int)(nbp));

       for (size_type i = 0; i < nbp; ++i)
         for (size_type k = 0; k < N; ++k)
           w1(k, i) = ps->ct_pairs()[i].slave_point[k];

       );

    /*@GET s = ('master points')
      Get the master points computed.
      @*/
    sub_command
      ("master points", 0, 0, 0, 1,

       size_type nbp = ps->ct_pairs().size();
       size_type N = ps->dim();
       darray w1 = out.pop().create_darray((unsigned int)(N),
                                           (unsigned int)(nbp));

       for (size_type i = 0; i < nbp; ++i)
         for (size_type k = 0; k < N; ++k)
           w1(k, i) = ps->ct_pairs()[i].master_point[k];

       );

    /*@GET s = ('char')
      Output a (unique) string representation of the @tmcf.
      
      This can be used for performing comparisons between two
      different @tmcf objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false, "Sorry, function to be done");
       // std::string s = ...;
       // out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
      Display a short summary for a @tmcf object.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfMultiContactFrame object\n";
       );

  }


  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::multi_contact_frame *ps = m_in.pop().to_multi_contact_frame();
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
