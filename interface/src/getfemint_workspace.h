// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Julien Pommier
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

#ifndef GETFEMINT_WORKSPACE_H__
#define GETFEMINT_WORKSPACE_H__

#include <getfem/dal_tas.h>
#include <getfemint_object.h>
namespace getfemint
{

  class workspace_data {
    std::string name;
    time_t creation_time;
  public:
    id_type parent_workspace;
    workspace_data() { name = "invalid"; creation_time = 0; parent_workspace = id_type(-2); }
    workspace_data(std::string n, id_type parent) : name(n), parent_workspace(parent) { 
      creation_time = ::time(NULL); 
    }
    const std::string& get_name() const { 
      return name; 
    }
    time_t get_creation_time() const { return creation_time; }    
  };

  class workspace_stack {
  public:
    typedef dal::dynamic_tas<getfem_object*>  obj_ct;
    typedef dal::dynamic_tas<workspace_data>  wrk_ct;
    static const id_type invalid_id = id_type(-1);
    static const id_type anonymous_workspace = getfem_object::anonymous_workspace;
  private:

    id_type current_workspace; /* stores the current workspace number */
    id_type base_workspace;    /* stores the initial workspace number */

    /* list of getfem object */
    obj_ct obj;
    /* list of used workspaces */
    wrk_ct wrk;

    std::map<getfem_object::internal_key_type, getfem_object*> kmap;

    std::vector<id_type> newly_created_objects;

    /* check if object 'id' can be deleted 
       (all objects which depend on it should be marked as anonymous)
     */
    void mark_deletable_objects(id_type id, dal::bit_vector& v) const;
  public:

    /* inserts a new object (and gives it an id) */
    id_type push_object(getfem_object *o);

    /* set the dependance of an object with respect to another */
    void set_dependance(getfem_object *user, getfem_object *used_by);
    void sup_dependance(getfem_object *user, getfem_object *used_by);

    /* at least mark the objet for future deletion (object becomes anonymous)
       and if possible, destroy the object (and all the objects which use this one
       if they are all anonymous) */
    void delete_object(id_type id);

    /* create a new workspace on top of the stack */
    void push_workspace(std::string n="unnamed");

    /* move the object in the parent workspace, in order to prevent
       the object from being deleted when the current workspace will
       be 'poped' */
    void send_object_to_parent_workspace(id_type obj_id);
    void send_all_objects_to_parent_workspace();

    /* delete every object in the workspace, but *does not* delete the workspace itself */
    void clear_workspace(id_type w);
    /* clears the current workspace */
    void clear_workspace() { clear_workspace(current_workspace); }

    /* deletes the current workspace and returns to the parent workspace */
    void pop_workspace(bool keep_all = false);

    /* throw an error if not found */
    getfem_object* object(id_type id, const char *expected_type="");

    /* return 0 is object not found */
    getfem_object* object(getfem_object::internal_key_type p);

    const obj_ct& get_obj_list() const { return obj; }
    const wrk_ct& get_wrk_list() const { return wrk; }
    id_type get_current_workspace() const { return current_workspace; }
    id_type get_base_workspace() const { return base_workspace; }
    workspace_stack() : current_workspace(anonymous_workspace) {
      push_workspace("main");
      base_workspace = current_workspace;
    }

    void commit_newly_created_objects();
    void destroy_newly_created_objects();
  };

  workspace_stack& workspace();
}
#endif
