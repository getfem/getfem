/*===========================================================================

 Copyright (C) 2002-2015 Julien Pommier.

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

===========================================================================*/
// $Id$
#define GETFEMINT_WORKSPACE_C

#include <getfem/dal_singleton.h>
#include <getfem/bgeot_config.h>
#include <getfemint_workspace.h>

namespace getfemint {

  workspace_stack& workspace() {
    return dal::singleton<workspace_stack>::instance();
  }

  /* deletes the current workspace and returns to the parent workspace */
  void workspace_stack::pop_workspace(bool keep_all) {
    if (wrk.size() == 1) THROW_ERROR("You cannot pop the main workspace\n");
    if (keep_all) send_all_objects_to_parent_workspace();
    else clear_workspace();
    wrk.pop_back();
  }

  /* inserts a new object (and gives it an id) */
  id_type workspace_stack::push_object(const dal::pstatic_stored_object &p,
					const void *raw_pointer,
					getfemint_class_id class_id) {
    id_type id = id_type(valid_objects.first_false());
    valid_objects.add(id);
    if (id >= obj.size()) obj.push_back(object_info());
    
    object_info &o = obj[id];
    o.p = p;
    o.raw_pointer = raw_pointer;
    o.workspace = get_current_workspace();
    o.class_id = class_id;
    o.used_by.clear();

    kmap[raw_pointer] = id;
    newly_created_objects.push_back(id);
    return id;
  }

  void workspace_stack::set_dependance(id_type user, id_type used) {
    if (!(valid_objects.is_in(user)) || !(valid_objects.is_in(used)))
      THROW_ERROR("Invalid objects\n");
    auto &u = obj[used].used_by;
    auto &p = obj[user].p;
    for (auto it = u.begin(); it != u.end(); ++it)
      if (it->get() == p.get()) return;
    u.push_back(p);
  }

  void workspace_stack::sup_dependance(id_type user, id_type used) {
    if (!(valid_objects.is_in(user)) || !(valid_objects.is_in(used)))
      THROW_ERROR("Invalid objects\n");
    auto &u = obj[used].used_by;
    auto &p = obj[user].p;
    size_type i = 0, j = 0;
    for ( ; i < u.size(); ++i)
      { u[j] = u[i]; if (u[i].get() != p.get()) ++j; }
    u.resize(j);
  }

  void workspace_stack::delete_object(id_type id) {
    if (valid_objects[id]) {
      valid_objects.sup(id);
      kmap.erase(obj[id].raw_pointer);
      obj[id] = object_info();
    } else {
      std::stringstream s;
      s << "Object number " << id << " no longer exists : can't delete it";
      throw getfemint_error(s.str());
    }
  }

  void workspace_stack::send_object_to_parent_workspace(id_type id) {
    if (get_current_workspace() == 0) THROW_ERROR("Invalid operation\n");
    if (!(valid_objects.is_in(id))) THROW_ERROR("Invalid objects\n");
    auto &o = obj[id];
    o.workspace = id_type(get_current_workspace() - 1);
  }

  void workspace_stack::send_all_objects_to_parent_workspace() {
    id_type cw = get_current_workspace();
    for (dal::bv_visitor_c id(valid_objects); !id.finished(); ++id)
      if ((obj[id]).workspace == cw) obj[id].workspace = id_type(cw-1);
  }

  void workspace_stack::clear_workspace(id_type wid) {
    if (wid > get_current_workspace()) THROW_INTERNAL_ERROR;
    dal::bit_vector bv = valid_objects;
    for (dal::bv_visitor_c id(bv); !id.finished(); ++id) {
      if (valid_objects.is_in(id)) {
	id_type owid = obj[id].workspace;
	if (owid > get_current_workspace()) THROW_INTERNAL_ERROR;
	if (owid == wid) delete_object(id_type(id));
      }
    }
  }

  const void *workspace_stack::object(id_type id,
				       const char *expected_type) const {
    if (valid_objects[id] &&
        std::find(newly_created_objects.begin(),newly_created_objects.end(),id)
	== newly_created_objects.end()) {
      return obj[id].raw_pointer;
    } else {
      THROW_ERROR("object " << expected_type << " [id=" << id << "] not found");
    }
    return 0;
  }

  const dal::pstatic_stored_object &workspace_stack::shared_pointer
  (id_type id, const char *expected_type) const {
    if (valid_objects[id] &&
        std::find(newly_created_objects.begin(),newly_created_objects.end(),id)
	== newly_created_objects.end()) {
      return obj[id].p;
    } else {
      THROW_ERROR("object " << expected_type << " [id=" << id << "] not found");
    }
  }

  id_type workspace_stack::object(const void *raw_pointer) const {
    auto it = kmap.find(raw_pointer);
    if (it != kmap.end()) return it->second; else return id_type(-1);
  }

  void workspace_stack::commit_newly_created_objects() {
    newly_created_objects.resize(0);
  }

  void workspace_stack::destroy_newly_created_objects() {
    while (newly_created_objects.size()) {
      delete_object(newly_created_objects.back());
      newly_created_objects.pop_back();
    }
  }



}
