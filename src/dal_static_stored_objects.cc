// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_static_stored_objects.cc : object which should be stored.
//           
// Date    : February 19, 2005
// Authors : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2002-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================


#include <dal_static_stored_objects.h>
#include <map>
#include <list>


namespace dal {

  // Pointer to an object with the dependencies
  struct pstatic_stored_object {
    const static_stored_object *p;
    bool valid;
    int permanence; // 0 = not deletable object
                    // 1 = preferable not to delete it
                    // 2 = standard
                    // 3 = delete it if memory is necessary
    std::map<const static_stored_object *, bool> dependent_object;
    std::map<const static_stored_object *, bool> dependencies;
    pstatic_stored_object(const static_stored_object *o, int perma)
      : p(o), valid(true), permanence(perma) {}
    pstatic_stored_object(void) : p(0), valid(true), permanence(2) {}
    
  };
  
  // Pointer to a key with a coherent order
  struct pstatic_stored_object_key {
    const static_stored_object_key *p;
    bool operator < (const pstatic_stored_object_key &o) const
    { return (*p) < (*(o.p)); }
    pstatic_stored_object_key(const static_stored_object_key *o) : p(o) {}
  };

  // Storing array types
  typedef std::map<pstatic_stored_object_key, pstatic_stored_object>
  stored_object_tab;
  typedef std::map<const static_stored_object *,
		   const static_stored_object_key *> 
  stored_key_tab;
  
  // Storing array
  static stored_object_tab stored_objects;
  static stored_key_tab stored_keys;
  
  // Gives a pointer to a key of an object from its pointer
  static inline const static_stored_object_key *
  key_of_object(const static_stored_object *o) {
    stored_key_tab::iterator it = stored_keys.lower_bound(o);
    if (it != stored_keys.end()) return it->second;
    return 0;
  }

  // Gives a pointer to an object from a key pointer
  const static_stored_object *
  search_stored_object(const static_stored_object_key *k) {
    stored_object_tab::iterator it
      = stored_objects.lower_bound(pstatic_stored_object_key(k));
    if (it != stored_objects.end()) return it->second.p;
    return 0;
  }

  // Gives an iterator on stored object from a pointer object
  static inline stored_object_tab::iterator 
  iterator_of_object(const static_stored_object *o) {
    const static_stored_object_key *k = key_of_object(o);
    if (k) return stored_objects.lower_bound(pstatic_stored_object_key(k));
    return stored_objects.end();
  }

  // Add a dependency, object o1 will depend on object o2
  void add_dependency(const static_stored_object *o1,
		      const static_stored_object *o2) {
    stored_object_tab::iterator it1 = iterator_of_object(o1);
    stored_object_tab::iterator it2 = iterator_of_object(o2);
    if (it1 != stored_objects.end() && it2 != stored_objects.end()) {
      it2->second.dependent_object[o1] = true;
      it1->second.dependencies[o2] = true;
    }
  }

  // remove a dependency
  void del_dependency(const static_stored_object *o1,
		      const static_stored_object *o2) {
    stored_object_tab::iterator it1 = iterator_of_object(o1);
    stored_object_tab::iterator it2 = iterator_of_object(o2);
    if (it1 != stored_objects.end() && it2 != stored_objects.end()) {
      it2->second.dependent_object.erase(o1);
      it1->second.dependencies.erase(o2);
    }
  }

  // Add an object with two optional dependencies
  void add_stored_object(const static_stored_object_key *k,
			 const static_stored_object *o,
			 int permanence,
			 const static_stored_object *dep1,
			 const static_stored_object *dep2) {
    stored_keys[o] = k;
    stored_objects[pstatic_stored_object_key(k)]
      = pstatic_stored_object(o, permanence);
    if (dep1) add_dependency(o, dep1);
    if (dep2) add_dependency(o, dep2);
  }

  // Only delete the object but not the dependencies
  static void basic_delete(std::list<const static_stored_object *> &to_delete){
    std::list<const static_stored_object *>::iterator it;
    for (it = to_delete.begin(); it != to_delete.end(); ++it) {
      const static_stored_object_key *k = key_of_object(*it);
      stored_object_tab::iterator ito = stored_objects.lower_bound(k);
      if (k) stored_keys.erase(*it);
      if (ito != stored_objects.end()) {
	delete ito->first.p;
	delete ito->second.p;
	stored_objects.erase(ito);
      }
    }
  }

  
  // Delete a list of objects and their dependencies
  void del_stored_objects(std::list<const static_stored_object *> &to_delete) {
    std::list<const static_stored_object *>::iterator it;
    std::map<const static_stored_object *, bool>::iterator itd;
    for (it = to_delete.begin(); it != to_delete.end(); ++it) {
      if (*it) {
	stored_object_tab::iterator ito = iterator_of_object(*it);
	if (ito != stored_objects.end()) {
	  ito->second.valid = false;
	  for (itd = ito->second.dependencies.begin();
	       itd != ito->second.dependencies.end(); ++itd)
	    del_dependency(*it, itd->first);
	  for (itd = ito->second.dependent_object.begin();
	       itd != ito->second.dependent_object.end(); ++itd) {
	    stored_object_tab::iterator itod=iterator_of_object(itd->first);
	    if (itod != stored_objects.end()) {
	      if (itod->second.permanence == 0)
		DAL_THROW(failure_error,"Trying to delete a permanent object");
	      if (itod->second.valid) to_delete.push_back(itod->second.p);
	    }
	  }
	}
      }
    }
    basic_delete(to_delete);
  }

  // Delete an object and its dependencies
  void del_stored_object(const static_stored_object *o) {
    std::list<const static_stored_object *> to_delete;
    to_delete.push_back(o);
    del_stored_objects(to_delete);
  }
  
  // Delete all the object whose permanence is greater or equal to perm
  void del_stored_objects(int perm) {
    if (perm == 0) ++perm;
    std::list<const static_stored_object *> to_delete;
    stored_object_tab::iterator it;
    for (it = stored_objects.begin(); it != stored_objects.end(); ++it)
      if (it->second.permanence >= perm)
	to_delete.push_back(it->second.p);
    del_stored_objects(to_delete);
  }


}
