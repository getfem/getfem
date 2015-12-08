/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2002-2015 Julien Pommier

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/
// $Id$
#ifndef GETFEMINT_WORKSPACE_H__
#define GETFEMINT_WORKSPACE_H__

#include <getfemint.h>
#include <getfemint_object.h>
#include <getfem/dal_bit_vector.h>
#include <getfem/dal_static_stored_objects.h>

namespace getfemint {

  class workspace_data {
    std::string name;
    time_t creation_time;
  public:
    workspace_data() { name = "invalid"; creation_time = 0; }
    workspace_data(std::string n) : name(n) { creation_time = ::time(NULL); }
    const std::string& get_name() const {  return name; }
    time_t get_creation_time() const { return creation_time; }    
  };

  class workspace_stack {
  public:
    typedef std::vector<getfem_object *>  obj_ct;
    typedef std::vector<workspace_data>  wrk_ct;
    static const id_type invalid_id = id_type(-1);
  private:

    id_type current_workspace; /* stores the current workspace number */
    id_type base_workspace;    /* stores the initial workspace number */

    /* list of getfem object */
    obj_ct obj;
    dal::bit_vector valid_objects;
    /* list of used workspaces */
    wrk_ct wrk;

    std::map<getfem_object::internal_key_type, getfem_object*> kmap;

    std::vector<id_type> newly_created_objects;

    /* check if object 'id' can be deleted 
       (all objects which depend on it should be marked as anonymous)
     */
    void mark_deletable_objects(id_type id, dal::bit_vector& v,
                                dal::bit_vector& g) const;
    void mark_deletable_objects(id_type id, dal::bit_vector& v) const
    { dal::bit_vector g; mark_deletable_objects(id, v, g); }
  public:

    /* inserts a new object (and gives it an id) */
    id_type push_object(getfem_object *o);

    /* set the dependance of an object with respect to another */
    void set_dependance(getfem_object *user, getfem_object *used_by);
    void sup_dependance(getfem_object *user, getfem_object *used_by);

    void undelete_object(id_type id);
    /* at least mark the objet for future deletion (object becomes anonymous)
       and if possible, destroy the object (and all the objects which use
       this one if they are all anonymous) */
    void delete_object(id_type id);

    /* create a new workspace on top of the stack */
    void push_workspace(std::string n="unnamed");

    /* move the object in the parent workspace, in order to prevent
       the object from being deleted when the current workspace will
       be 'poped' */
    void send_object_to_parent_workspace(id_type obj_id);
    void send_all_objects_to_parent_workspace();

    /* delete every object in the workspace, but *does not* delete the
       workspace itself */
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
    const dal::bit_vector& get_obj_index() const { return valid_objects; }
    const wrk_ct& get_wrk_list() const { return wrk; }
    id_type get_current_workspace() const { return current_workspace; }
    id_type get_base_workspace() const { return base_workspace; }
    workspace_stack() : current_workspace(anonymous_workspace)
    { push_workspace("main"); base_workspace = current_workspace; }

    void commit_newly_created_objects();
    void destroy_newly_created_objects();
  };

  workspace_stack& workspace();






  // The workspace_stack structure stores the various object assigned to
  // the script languages variables (Python, Scilab or Matlab).
  // This objects are organised in workspaces (for Matlab and Scilab only) in
  // order to be able to delete all the variables created locally in a
  // sub-program or a loop (however, this has to be managed by the user).
  // Additionnally, a variable may have dependances with respect to some
  // other variables, in the sense that is a deletion of a variable occurs
  // it will be delayed untill all the dependant variables are deleted
  // (implemented with shared pointers).

  class workspace2_stack {
    
    struct object_info {
      dal::pstatic_stored_object p;
      const void *raw_pointer;
      id_type workspace;
      getfemint_class_id class_id;
      std::vector<dal::pstatic_stored_object> used_by;

      object_info() : raw_pointer(0), class_id(GETFEMINT_NB_CLASS) {}
    };

    typedef std::vector<object_info>  obj_ct;
    typedef std::vector<std::string>  wrk_ct;
    static const id_type invalid_id = id_type(-1);

    obj_ct obj;                      // Array of getfem object.
    dal::bit_vector valid_objects;   // Indices of valid objects.
    wrk_ct wrk;                      // Stack of used workspaces.

    std::map<const void *, id_type> kmap;
    std::vector<id_type> newly_created_objects;

  public:

    // Creates a new workspace2 on top of the stack
    void push_workspace2(const std::string &n = "Unnamed") { wrk.push_back(n); }

    // Deletes the current workspace2 and returns to the parent workspace2
    void pop_workspace2(bool keep_all = false);

    // Inserts a new object (and gives it an id)
    id_type push_object(const dal::pstatic_stored_object &p,
			const void *raw_pointer, getfemint_class_id class_id);

    // Sets the dependance of an object with respect to another
    void set_dependance(id_type user, id_type used);
    void sup_dependance(id_type user, id_type used);

    /** At least mark the objet for future deletion (object becomes anonymous)
	and if possible, destroy the object (and all the objects which use
	this one if they are all anonymous).
    */
    void delete_object(id_type id);

    /* Move the object in the parent workspace2, in order to prevent
       the object from being deleted when the current workspace2 will
       be 'poped' */
    void send_object_to_parent_workspace2(id_type obj_id);
    void send_all_objects_to_parent_workspace2();

    id_type get_current_workspace2() const { return id_type(wrk.size()-1); }
    id_type get_base_workspace2() const { return id_type(0); }
    /* Delete every object in the workspace2, but *does not* delete the
       workspace2 itself */
    void clear_workspace2(id_type w);
    /* clears the current workspace2 */
    void clear_workspace2() { clear_workspace2(get_current_workspace2()); }

    
    /* Throw an error if not found */
    const void *object(id_type id, const char *expected_type="") const;

    /* Throw an error if not found */
    const dal::pstatic_stored_object &shared_pointer
    (id_type id, const char *expected_class="") const;

    /* Return id_type(-1) if not found */
    id_type object(const void *raw_pointer) const;
    
    workspace2_stack() { push_workspace2("main"); }

    void commit_newly_created_objects();
    void destroy_newly_created_objects();
  };

  workspace2_stack& workspace2();









}
#endif
