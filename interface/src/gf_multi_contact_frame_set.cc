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
// $Id: gf_multi_contact_frame_set.cc 4114 2012-07-06 11:20:10Z renard $
#include <getfemint.h>
#include <getfemint_multi_contact_frame.h>
#include <getfemint_workspace.h>
#include <getfemint_models.h>
#include <getfemint_mesh_im.h>

using namespace getfemint;

/*@GFDOC
  General function for modification of @tmcf objects.
@*/




// Object for the declaration of a new sub-command.

struct sub_gf_mcf_set : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::multi_contact_frame *ps) = 0;
};

typedef boost::intrusive_ptr<sub_gf_mcf_set> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_mcf_set {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::multi_contact_frame *ps)			\
      { dummy_func(in); dummy_func(out);  dummy_func(ps); code }	\
    };									\
    psub_command psubc = new subc;					\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           




void gf_multi_contact_frame_set(getfemint::mexargs_in& m_in,
                                getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;
  
  if (subc_tab.size() == 0) {

    /*@SET ('add obstacle', @str obs)
    Add a rigid obstacle. The string `obs` is the expression of a
    function which should be closed to a signed distance to the obstacle.
    @*/
    sub_command
      ("add obstacle", 1, 1, 0, 1,
       std::string obs = in.pop().to_string();
       size_type ind = ps->add_obstacle(obs);
       out.pop().from_integer(int(ind + config::base_index()));
       );

    /*@SET ('add slave boundary', @tmim mim, @int region, @str varname [, @str multname])
    Add a slave contact bounary.
    @*/
    
    sub_command
      ("add slave boundary", 3, 4, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       size_type region = in.pop().to_integer();
       std::string varname = in.pop().to_string();
       std::string multname;
       if (in.remaining()) multname = in.pop().to_string();
       size_type ind = ps->add_slave_boundary(gfi_mim->mesh_im(), region,
                                              varname, multname);
       out.pop().from_integer(int(ind + config::base_index()));
       );

    /*@SET ('add master boundary', @tmim mim, @int region, @str varname [, @str multname])
    Add a master contact bounary.
    @*/
    
    sub_command
      ("add master boundary", 3, 4, 0, 1,
       getfemint_mesh_im *gfi_mim = in.pop().to_getfemint_mesh_im();
       size_type region = in.pop().to_integer();
       std::string varname = in.pop().to_string();
       std::string multname;
       if (in.remaining()) multname = in.pop().to_string();
       size_type ind = ps->add_master_boundary(gfi_mim->mesh_im(), region,
                                               varname, multname);
       out.pop().from_integer(int(ind + config::base_index()));
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
