/*===========================================================================

 Copyright (C) 2006-2020 Julien Pommier.

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

#include <getfemint_workspace.h>
#include <getfemint_levelset.h>
#include <getfem/getfem_mesh_fem.h>

using namespace getfemint;

/*@GFDOC
    General function for querying information about LEVELSET objects.
@*/



// Object for the declaration of a new sub-command.

struct sub_gf_ls_get : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out,
		   getfem::level_set &ls) = 0;
};

typedef std::shared_ptr<sub_gf_ls_get> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_ls_get {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out,			\
		       getfem::level_set &ls)				\
      { dummy_func(in); dummy_func(out); dummy_func(ls); code }		\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }                           







void gf_levelset_get(getfemint::mexargs_in& m_in,
		     getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;
  
  if (subc_tab.size() == 0) {
    

    /*@GET V = ('values', @int nls)
    Return the vector of dof for `nls` function.

    If `nls` is 0, the method return the vector of dof for the primary
    level-set function. If `nls` is 1, the method return the vector of
    dof for the secondary level-set function (if any).@*/
    sub_command
      ("values", 0, 1, 0, 1,
       size_type il = 0;
       if (in.remaining()) il = in.pop().to_integer(0, 1);
       if (il != 0 && !ls.has_secondary())
	 THROW_BADARG("The levelset has not secondary term");
       out.pop().from_dcvector(ls.values(unsigned(il)));
       );


    /*@RDATTR d = ('degree')
      Return the degree of lagrange representation.@*/
    sub_command
      ("degree", 0, 0, 0, 1,
       out.pop().from_integer(ls.degree());
       );


    /*@GET mf = ('mf')
    Return a reference on the @tmf object.@*/
    sub_command
      ("mf", 0, 0, 0, 1,
       getfem::mesh_fem *pmf=const_cast<getfem::mesh_fem*>(&ls.get_mesh_fem());
       id_type id = workspace().object(pmf);
       if (id == id_type(-1)) {
	 id = store_meshfem_object(std::shared_ptr<getfem::mesh_fem>
				   (std::shared_ptr<getfem::mesh_fem>(), pmf));
       }
       out.pop().from_object_id(id, MESHFEM_CLASS_ID);
       );


    /*@RDATTR z = ('memsize')
      Return the amount of memory (in bytes) used by the level-set.@*/
    sub_command
      ("memsize", 0, 0, 0, 1,
       out.pop().from_integer(int(ls.memsize()));
       );


    /*@GET s = ('char')
      Output a (unique) string representation of the @tls.

      This can be used to perform comparisons between two
      different @tls objects.
      This function is to be completed.
      @*/
    sub_command
      ("char", 0, 0, 0, 1,
       GMM_ASSERT1(false, "Sorry, function to be done");
       // std::string s = ...;
       // out.pop().from_string(s.c_str());
       );


    /*@GET ('display')
      displays a short summary for a @tls.@*/
    sub_command
      ("display", 0, 0, 0, 0,
       infomsg() << "gfLevelSet object\n";
       );


  }



  if (m_in.narg() < 2)  THROW_BADARG( "Wrong number of input arguments");

  getfem::level_set &ls = *(to_levelset_object(m_in.pop()));
  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out, ls);
  }
  else bad_cmd(init_cmd);

}
