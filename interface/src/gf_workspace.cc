// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

#include <getfemint.h>
#include <getfemint_workspace.h>
#include <algorithm>
#include <unistd.h>
#include <iomanip>
#include <getfem/getfem_mat_elem.h>
#include <getfemint_mdbrick.h>

using namespace getfemint;

  class match_workspace {
    id_type wid;
  public:
    match_workspace(id_type wid_) { wid = wid_; }
    bool operator()(const getfem_object *o) {
      return o->get_workspace() == wid;
    }
  };

static void
do_stat(id_type wid) {
  const workspace_stack::obj_ct obj = workspace().get_obj_list();
  const workspace_stack::wrk_ct wrk = workspace().get_wrk_list();
  int cnt = 0;

  if (wid == workspace_stack::anonymous_workspace) {
    infomsg() << "Anonymous workspace (objects waiting for deletion)\n";
  } else {
    if (!wrk.index()[wid]) THROW_INTERNAL_ERROR;
    infomsg() << "Workspace " << wid << " [" << wrk[wid].get_name() << " -- " << 
      std::count_if(obj.tas_begin(), obj.tas_end(), match_workspace(wid)) << " objects]\n";
  }
  for (workspace_stack::obj_ct::const_tas_iterator it = obj.tas_begin(); it != obj.tas_end(); ++it, ++cnt) {
    if (!obj.index()[it.index()]) THROW_INTERNAL_ERROR;
    if (match_workspace(wid)(*it)) {
      std::string subclassname;
      if ((*it)->class_id() == MDBRICK_CLASS_ID)
	subclassname = "(" + dynamic_cast<getfemint_mdbrick*>(*it)->sub_class() + ")";
      infomsg() << " ID" << std::setw(4) << (*it)->get_id() << " " 
		<< std::setw(14) << name_of_getfemint_class_id((*it)->class_id())
	        << std::setw(30) << subclassname 
		<< "   " << std::setw(9) << (*it)->memsize() << " bytes";
      if ((*it)->is_static()) infomsg() << " * "; else infomsg() << "   ";
      if ((*it)->is_const()) infomsg() << "Const"; else infomsg() << "     ";
      const std::vector<id_type>& used_by = (*it)->get_used_by();
      if (used_by.size()) {
	infomsg() << " used by ";
	for (size_type i=0; i < used_by.size(); ++i) infomsg() << " ID" << used_by[i];
      }
      //      if ((*it)->is_anonymous()) mexPrintf("[ano]");
      infomsg() << endl;
    }
  }
}

static void
do_stats() {
  //infomsg() << "Memory used by elementary matrices structures : " << getfem::stored_mat_elem_memsize()/1024 << "KB\n";
  if (std::count_if(workspace().get_obj_list().tas_begin(), 
		    workspace().get_obj_list().tas_end(), 
		    match_workspace(workspace_stack::anonymous_workspace))) {
    do_stat(workspace_stack::anonymous_workspace);
  }
  for (dal::bv_visitor_c wid(workspace().get_wrk_list().index()); !wid.finished(); ++wid)
    do_stat(wid);
}

/*MLABCOM
  FUNCTION gf_workspace(operation)
    Getfem workspace management function. 

    Getfem uses its own workspaces in Matlab, independently of the
    matlab workspaces (this is due to some limitations in the memory
    management of matlab objects). By default, all getfem variables
    belong to the root getfem workspace. A function can create its own
    workspace by invoking gf_workspace('push') at its beginning. When
    exiting, this function MUST invoke gf_workspace('pop') (you can
    use matlab exceptions handling to do this cleanly when the
    function exits on an error).

    USAGE:
    * gf_workspace('push') 
    Create a new temporary workspace on the workspace stack.

    * gf_workspace('pop' [,i,j,..])  
    Leave the current workspace, destroying all getfem objects
    belonging to it, except the one listed after 'pop', and the ones
    moved to parent workspace by gf_workspace('keep').

    * gf_workspace('stat')
    Print informations about variables in current workspace.

    * gf_workspace('stats')
    Print informations about all getfem variables.

    * gf_workspace('keep', i[,j,k..]) 
    prevent the listed variables i from being deleted when
    gf_workspace("pop") will be called by moving this variable in the
    parent workspace.

    * gf_workspace('clear')
    Clear the current workspace.

    * gf_workspace('clear all') 
    Clear every workspace, and returns to the main workspace (you
    should not need this command).

    * S=gf_workspace('class name', i)
    Return the class name of object i (if I is a mesh handle, it 
    return gfMesh etc..)

 MLABCOM*/


void gf_workspace(getfemint::mexargs_in& in, getfemint::mexargs_out& out)
{
  if (in.narg() < 1) {
    THROW_BADARG( "Wrong number of input arguments");
  }
  mexarg_in arg = in.pop();

  std::string cmd = arg.to_string();

  if (check_cmd(cmd, "push", in, out, 0, 1, 0, 0)) {
    std::string s = "unnamed";
    if (!in.remaining() == 0) s = in.pop().to_string();
    workspace().push_workspace(s);
  } else if (check_cmd(cmd, "pop", in, out, 0, 256, 0, 0)) {
    if (workspace().get_current_workspace() != workspace().get_base_workspace()) {
      while (in.remaining()) {
	workspace().send_object_to_parent_workspace(in.pop().to_object_id());
      }    
      workspace().pop_workspace();
    } else THROW_ERROR("Can't pop main workspace");
  } else if (check_cmd(cmd, "stat", in, out, 0, 0, 0, 0)) {
    do_stat(workspace().get_current_workspace());
    infomsg() << endl;
  } else if (check_cmd(cmd, "stats", in, out, 0, 0, 0, 0)) {
    do_stats();
    infomsg() << endl;
  } else if (check_cmd(cmd, "keep", in, out, 1, 256, 0, 0)) {
    while (in.remaining()) {
      workspace().send_object_to_parent_workspace(in.pop().to_object_id());
    }
  } else if (check_cmd(cmd, "keep all", in, out, 0, 0, 0, 0)) {
    workspace().send_all_objects_to_parent_workspace();
  } else if (check_cmd(cmd, "clear", in, out, 0, 0, 0, 0)) {
    workspace().clear_workspace();
  } else if (check_cmd(cmd, "clear all", in, out, 0, 0, 0, 0)) {
    while (workspace().get_current_workspace() != workspace().get_base_workspace()) {
      workspace().pop_workspace();
      //      mexPrintf("w <- %d\n", workspace().get_current_workspace());
    }
    workspace().clear_workspace();
  } else if (check_cmd(cmd, "chdir", in, out, 1, 1, 0, 0)) {
    ::chdir(in.pop().to_string().c_str());
  } else if (check_cmd(cmd, "class name", in, out, 0, 1, 0, 1)) {
    id_type id, cid; in.pop().to_object_id(&id, &cid);
    out.pop().from_string(name_of_getfemint_class_id(cid));
  } else if (check_cmd(cmd, "connect", in, out, 0, -1, 0, -1)) {
    GMM_THROW(getfemint_error, "cannot connect: the toolbox was built without rpc support");
  } else bad_cmd(cmd);
}
