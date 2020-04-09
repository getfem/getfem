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

#include <getfemint_workspace.h>
#include <algorithm>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <iomanip>

using namespace getfemint;
 
/*@GFDOC
    Getfem workspace management function. 

    Getfem uses its own workspaces in Matlab, independently of the
    matlab workspaces (this is due to some limitations in the memory
    management of matlab objects). By default, all getfem variables
    belong to the root getfem workspace. A function can create its own
    workspace by invoking gf_workspace('push') at its beginning. When
    exiting, this function MUST invoke gf_workspace('pop') (you can
    use matlab exceptions handling to do this cleanly when the
    function exits on an error).

 @*/



// Object for the declaration of a new sub-command.

struct sub_gf_workspace : virtual public dal::static_stored_object {
  int arg_in_min, arg_in_max, arg_out_min, arg_out_max;
  virtual void run(getfemint::mexargs_in& in,
		   getfemint::mexargs_out& out) = 0;
};

typedef std::shared_ptr<sub_gf_workspace> psub_command;

// Function to avoid warning in macro with unused arguments.
template <typename T> static inline void dummy_func(T &) {}

#define sub_command(name, arginmin, arginmax, argoutmin, argoutmax, code) { \
    struct subc : public sub_gf_workspace {				\
      virtual void run(getfemint::mexargs_in& in,			\
		       getfemint::mexargs_out& out)			\
      { dummy_func(in); dummy_func(out); code }				\
    };									\
    psub_command psubc = std::make_shared<subc>();			\
    psubc->arg_in_min = arginmin; psubc->arg_in_max = arginmax;		\
    psubc->arg_out_min = argoutmin; psubc->arg_out_max = argoutmax;	\
    subc_tab[cmd_normalize(name)] = psubc;				\
  }


void gf_workspace(getfemint::mexargs_in& m_in, getfemint::mexargs_out& m_out) {
  typedef std::map<std::string, psub_command > SUBC_TAB;
  static SUBC_TAB subc_tab;

  if (subc_tab.size() == 0) {

    /*@FUNC ('push')
      Create a new temporary workspace on the workspace stack. @*/
    sub_command
      ("push", 0, 1, 0, 0,
       std::string s = "unnamed";
       if (!in.remaining() == 0) s = in.pop().to_string();
       workspace().push_workspace(s);
       );

    /*@FUNC ('pop',  [,i,j, ...])
      Leave the current workspace, destroying all getfem objects
      belonging to it, except the one listed after 'pop', and the ones
      moved to parent workspace by ::WORKSPACE('keep'). @*/
    sub_command
      ("pop", 0, 256, 0, 0,
       if (workspace().get_current_workspace()
	   != workspace().get_base_workspace()) {
	 while (in.remaining()) {
	   workspace().send_object_to_parent_workspace
	     (in.pop().to_object_id());
	 }    
	 workspace().pop_workspace();
       } else THROW_ERROR("Can't pop main workspace");
       );


    /*@FUNC ('stat')
       Print informations about variables in current workspace. @*/
    sub_command
      ("stat", 0, 0, 0, 0,
       workspace().do_stats(infomsg(), workspace().get_current_workspace());
       infomsg() << endl;
       );


    /*@FUNC ('stats')
       Print informations about all getfem variables. @*/
    sub_command
      ("stats", 0, 0, 0, 0,
       workspace().do_stats(infomsg());
       infomsg() << endl;
       );


    /*@FUNC ('keep', i[,j,k...]) 
      prevent the listed variables from being deleted when
      ::WORKSPACE("pop") will be called by moving these variables in the
      parent workspace. @*/
    sub_command
      ("keep", 1, 256, 0, 0,
       while (in.remaining()) {
	 workspace().send_object_to_parent_workspace(in.pop().to_object_id());
       }
       );


    /*@FUNC ('keep all') 
      prevent all variables from being deleted when
      ::WORKSPACE("pop") will be called by moving the variables in the
      parent workspace. @*/
    sub_command
      ("keep all", 0, 0, 0, 0,
       workspace().send_all_objects_to_parent_workspace();
       );

    /*@FUNC ('clear')
      Clear the current workspace. @*/
    sub_command
      ("clear", 0, 0, 0, 0,
       workspace().clear_workspace();
       );

    /*@FUNC ('clear all')
      Clear every workspace, and returns to the main workspace (you
      should not need this command). @*/
    sub_command
      ("clear all", 0, 0, 0, 0,
       while (workspace().get_current_workspace()
	      != workspace().get_base_workspace()) {
	 workspace().pop_workspace();
	 //      mexPrintf("w <- %d\n", workspace().get_current_workspace());
       }
       workspace().clear_workspace();
       );


    /* Unofficial function */
#ifndef _MSC_VER
    sub_command
      ("chdir", 1, 1, 0, 0,
       if (::chdir(in.pop().to_string().c_str())) {}
       );
#endif

    /*@FUNC ('class name', i)
      Return the class name of object i (if I is a mesh handle, it 
      return gfMesh etc..) @*/
    sub_command
      ("class name", 0, 1, 0, 1,
       id_type id;  id_type cid;
       in.pop().to_object_id(&id, &cid);
       out.pop().from_string(name_of_getfemint_class_id(cid));
       );

    /* Unofficial function */
    sub_command
      ("connect", 0, -1, 0, -1,
       GMM_THROW(getfemint_error, "cannot connect: the toolbox was built "
		 "without rpc support");
       );


    /* Unofficial function */
    sub_command
      ("list static objects", 0, -1, 0, -1,
       dal::list_stored_objects(cout);
       );


    /* Unofficial function */
    sub_command
      ("nb static objects", 0, -1, 0, 1,
       out.pop().from_integer(int(dal::nb_stored_objects()));
       );

  }





  if (m_in.narg() < 1)  THROW_BADARG( "Wrong number of input arguments");

  std::string init_cmd   = m_in.pop().to_string();
  std::string cmd        = cmd_normalize(init_cmd);

  
  SUBC_TAB::iterator it = subc_tab.find(cmd);
  if (it != subc_tab.end()) {
    check_cmd(cmd, it->first.c_str(), m_in, m_out, it->second->arg_in_min,
	      it->second->arg_in_max, it->second->arg_out_min,
	      it->second->arg_out_max);
    it->second->run(m_in, m_out);
  }
  else bad_cmd(init_cmd);

}
