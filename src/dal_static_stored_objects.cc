/*===========================================================================

Copyright (C) 2002-2012 Yves Renard

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


#include "getfem/dal_static_stored_objects.h"
#include "getfem/dal_singleton.h"
#include <map>
#include <list>
#include <set>
#include <algorithm>
#include <deque>


namespace dal {


  // Gives a pointer to a key of an object from its pointer, while looking in the storage of
  // a specific thread
  pstatic_stored_object_key key_of_stored_object(pstatic_stored_object o, size_t thread) 
  {
    stored_object_tab::stored_key_tab& stored_keys 
      = dal::singleton<stored_object_tab>::instance(thread).stored_keys_;
    stored_object_tab::stored_key_tab::iterator it = stored_keys.find(o);
    if (it != stored_keys.end()) return it->second;
    return 0;
  }

  // gives a key of the stored object while looking in the storage of other threads
  pstatic_stored_object_key key_of_stored_object_other_threads(pstatic_stored_object o) 
  {
    for(size_t thread = 0; thread<getfem::num_threads();thread++)
    {
      if (thread == this_thread()) continue;
      pstatic_stored_object_key key = key_of_stored_object(o,thread);
      if (key) return key;
    }
    return 0;
  }

  pstatic_stored_object_key key_of_stored_object(pstatic_stored_object o) 
  {
    pstatic_stored_object_key key = key_of_stored_object(o,this_thread());
    if (key) return key;
    else return (num_threads() > 1) ? key_of_stored_object_other_threads(o) : 0;
  }

  bool exists_stored_object(pstatic_stored_object o) 
  {
    stored_object_tab::stored_key_tab& stored_keys 
       = dal::singleton<stored_object_tab>::instance().stored_keys_;
    return  (stored_keys.find(o) != stored_keys.end());
  }



/**
  STATIC_STORED_TAB -------------------------------------------------------
*/
  stored_object_tab::stored_object_tab() 
    : std::map<enr_static_stored_object_key, enr_static_stored_object>(),
    locks_(),
    stored_keys_()
  { }

  pstatic_stored_object stored_object_tab::
    search_stored_object(pstatic_stored_object_key k) const
  {
   getfem::local_guard guard = locks_.get_lock();
   stored_object_tab::const_iterator it = find(enr_static_stored_object_key(k));
   return (it != end()) ? it->second.p : 0;
  }

  bool stored_object_tab::add_dependency_(pstatic_stored_object o1,
    pstatic_stored_object o2)
  {
    getfem::local_guard guard = locks_.get_lock();
    stored_key_tab::const_iterator it = stored_keys_.find(o1);
    if (it == stored_keys_.end()) return false;
    iterator ito1 = find(it->second);
    GMM_ASSERT1(ito1 != end(), "Object has a key, but cannot be found");
    ito1->second.dependencies.insert(o2);
    return true;
  }
  void stored_object_tab::add_stored_object(pstatic_stored_object_key k, 
    pstatic_stored_object o,  permanence perm) 
  {
    getfem::local_guard guard = locks_.get_lock();
    GMM_ASSERT1(stored_keys_.find(o) == stored_keys_.end(),
      "This object has already been stored, possibly with another key");
    stored_keys_[o] = k;
    insert(std::make_pair(enr_static_stored_object_key(k),enr_static_stored_object(o, perm)));
    size_t t = this_thread();
    GMM_ASSERT2(stored_keys_.size() == size(), 
      "stored_keys are not consistent with stored_object tab");
  }

  bool stored_object_tab::add_dependent_(pstatic_stored_object o1,
    pstatic_stored_object o2)
  {
    getfem::local_guard guard = locks_.get_lock();
    stored_key_tab::const_iterator it = stored_keys_.find(o2);
    if (it == stored_keys_.end()) return false;
    iterator ito2 = find(it->second);
    GMM_ASSERT1(ito2 != end(), "Object has a key, but cannot be found");
    ito2->second.dependent_object.insert(o1);
    return true;
  }

  bool stored_object_tab::del_dependency_(pstatic_stored_object o1,
    pstatic_stored_object o2)
  {
    getfem::local_guard guard = locks_.get_lock();
    stored_key_tab::const_iterator it1 = stored_keys_.find(o1);
    if (it1 == stored_keys_.end()) return false;
    iterator ito1 = find(it1->second);
    GMM_ASSERT1(ito1 != end(), "Object has a key, but cannot be found");
    ito1->second.dependencies.erase(o2);
    return true;
  }

  stored_object_tab::iterator stored_object_tab::iterator_of_object_(pstatic_stored_object o)
  {
    stored_key_tab::const_iterator itk = stored_keys_.find(o);
    if (itk == stored_keys_.end()) return end();
    iterator ito = find(itk->second);
    GMM_ASSERT1(ito != end(), "Object has a key, but is not stored");
    return ito;
  }

  bool stored_object_tab::del_dependent_(pstatic_stored_object o1,
    pstatic_stored_object o2)
  {
    getfem::local_guard guard = locks_.get_lock();
    stored_key_tab::const_iterator it2 = stored_keys_.find(o2);
    if (it2 == stored_keys_.end()) return false;
    iterator ito2 = find(it2->second);
    GMM_ASSERT1(ito2 != end(), "Object has a key, but cannot be found");
    ito2->second.dependent_object.erase(o1);
    return true;
  }

  bool stored_object_tab::exists_stored_object(pstatic_stored_object o) const
  {
    getfem::local_guard guard = locks_.get_lock();
    return (stored_keys_.find(o) != stored_keys_.end());
  }

  bool stored_object_tab::has_dependent_objects(pstatic_stored_object o) const
  {
    getfem::local_guard guard = locks_.get_lock();
     stored_key_tab::const_iterator it = stored_keys_.find(o);
    GMM_ASSERT1(it != stored_keys_.end(), "Object is not stored");
    const_iterator ito = find(it->second);
    GMM_ASSERT1(ito != end(), "Object has a key, but cannot be found");
    return ito->second.dependent_object.empty();
  }



  void stored_object_tab::basic_delete_(std::list<pstatic_stored_object> &to_delete)
  {
    getfem::local_guard guard = locks_.get_lock();
    std::list<pstatic_stored_object>::iterator it;
    for (it = to_delete.begin(); it != to_delete.end(); ++it) 
    {
      stored_key_tab::iterator itk = stored_keys_.find(*it);
      stored_object_tab::iterator ito = end();
      if (itk != stored_keys_.end())
      {
          ito = find(itk->second);
          stored_keys_.erase(itk);
      }
      if (ito != end()) 
      {
        delete ito->first.p;
        erase(ito);
        it = to_delete.erase(it);
        --it;
      }
    }
  }


  pstatic_stored_object search_stored_object(pstatic_stored_object_key k) 
  {
    stored_object_tab& stored_objects
        = dal::singleton<stored_object_tab>::instance();
    pstatic_stored_object p = stored_objects.search_stored_object(k);
    if (p) return p;
    return 0;
  }

  stored_object_tab::iterator iterator_of_object(pstatic_stored_object o)
  {
    for(size_t thread=0; thread < num_threads(); ++thread)
    {
      stored_object_tab& stored_objects
        = dal::singleton<stored_object_tab>::instance(thread);
      stored_object_tab::iterator it = stored_objects.iterator_of_object_(o);
      if (it != stored_objects.end()) return it;
    }
    return dal::singleton<stored_object_tab>::instance().end();
  }


  void test_stored_objects(void) 
  {
    for(size_t thread = 0; thread < num_threads(); ++thread)
    {
      stored_object_tab& stored_objects 
        = dal::singleton<stored_object_tab>::instance(thread);
      stored_object_tab::stored_key_tab& stored_keys = stored_objects.stored_keys_;

      GMM_ASSERT1(stored_objects.size() == stored_keys.size(), 
        "keys and objects tables don't match");
      for (stored_object_tab::stored_key_tab::iterator it = stored_keys.begin();
        it != stored_keys.end(); ++it)
      {
        stored_object_tab::iterator ito = iterator_of_object(it->first);
        GMM_ASSERT1(ito !=dal::singleton<stored_object_tab>::instance().end(), 
          "Key without object is found");
      }
      for (stored_object_tab::iterator it = stored_objects.begin();
        it != stored_objects.end(); ++it)
        GMM_ASSERT1(iterator_of_object(it->second.p) != stored_objects.end(),
        "Object has key but cannot be found");
    }
  }

  void add_dependency(pstatic_stored_object o1,
    pstatic_stored_object o2) {
    bool dep_added = false;
    for(size_t thread=0; thread < num_threads(); ++thread)
    {
      stored_object_tab& stored_objects
          = dal::singleton<stored_object_tab>::instance(thread);
      if (dep_added = stored_objects.add_dependency_(o1,o2)) break;
    }
    GMM_ASSERT1(dep_added, "Failed to add dependency between " << o1 << " of type "
    << typeid(*o1).name() << " and " << o2 << " of type "  << typeid(*o2).name() << ". ");

    bool dependent_added = false;
    for(size_t thread=0; thread < num_threads(); ++thread)
    {
      stored_object_tab& stored_objects
          = dal::singleton<stored_object_tab>::instance(thread);
      if (dependent_added =stored_objects.add_dependent_(o1,o2)) break;
    }
    GMM_ASSERT1(dependent_added, "Failed to add dependent between " << o1 << " of type "
    << typeid(*o1).name() << " and " << o2 << " of type "  << typeid(*o2).name() << ". ");
  }



  /*remove a dependency (from storages of all threads). 
  Return true if o2 has no more dependent object. */
  bool del_dependency(pstatic_stored_object o1,
    pstatic_stored_object o2) 
  {
    bool dep_deleted = false;
    for(size_t thread=0; thread < num_threads(); ++thread)
    {
      stored_object_tab& stored_objects
          = dal::singleton<stored_object_tab>::instance(thread);
      if (dep_deleted = stored_objects.del_dependency_(o1,o2)) break;
    }
    GMM_ASSERT1(dep_deleted, "Failed to delete dependency between " << o1 << " of type "
    << typeid(*o1).name() << " and " << o2 << " of type "  << typeid(*o2).name() << ". ");

    bool dependent_deleted = false;
    bool dependent_empty = false;
    for(size_t thread=0; thread < num_threads(); ++thread)
    {
      stored_object_tab& stored_objects
          = dal::singleton<stored_object_tab>::instance(thread);
      if (dependent_deleted = stored_objects.del_dependent_(o1,o2))
      {
        dependent_empty = stored_objects.has_dependent_objects(o2);
        break;
      }
    }
    GMM_ASSERT1(dependent_deleted, "Failed to delete dependent between " << o1 << " of type "
    << typeid(*o1).name() << " and " << o2 << " of type "  << typeid(*o2).name() << ". ");

    return dependent_empty;
  }


  void add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
    permanence perm) 
  {
      dal::singleton<stored_object_tab>::instance().add_stored_object(k,o,perm);
  }

  void basic_delete(std::list<pstatic_stored_object> &to_delete)
  {
      stored_object_tab& stored_objects_this_thread
        = dal::singleton<stored_object_tab>::instance();
      stored_objects_this_thread.basic_delete_(to_delete);
    
      if (!to_delete.empty()) //need to delete from other threads
      {
        for(size_t thread=0; thread < num_threads(); ++thread)
        { 
          if (thread == this_thread()) continue;
          stored_object_tab& stored_objects
              = dal::singleton<stored_object_tab>::instance(thread);
          stored_objects.basic_delete_(to_delete);
          if (to_delete.empty()) break;
        }
      }
      if (me_is_multithreaded_now())
      {
        if (!to_delete.empty()) GMM_WARNING1("Not all objects were deleted");
      }
      else
      {
        GMM_ASSERT1(to_delete.empty(), "Could not delete objects");
      }
  }

  void del_stored_objects(std::list<pstatic_stored_object> &to_delete,
    bool ignore_unstored) 
  {
      getfem::omp_guard lock;
      stored_object_tab& stored_objects
        = dal::singleton<stored_object_tab>::instance();
      std::list<pstatic_stored_object>::iterator it, itnext;

      for (it = to_delete.begin(); it != to_delete.end(); it = itnext) 
      {
        itnext = it; itnext++;
        stored_object_tab::iterator ito = iterator_of_object(*it);
        if (ito == stored_objects.end()) 
        {
          if (ignore_unstored)
            to_delete.erase(it);
          else
            if (me_is_multithreaded_now())
            {
              GMM_WARNING1("This object is (already?) not stored : " << it->get()
              << " typename: " << typeid(*it->get()).name() 
              << "(which could happen in multithreaded code and is OK)");
            }
            else
            {
              GMM_ASSERT1(false, "This object is not stored : " << it->get()
              << " typename: " << typeid(*it->get()).name());
            }
        }
        else
          ito->second.valid = false;
      }

      std::set<pstatic_stored_object>::iterator itd;
      for (it = to_delete.begin(); it != to_delete.end(); ++it) 
      {
        if (*it) 
        {
          stored_object_tab::iterator ito = iterator_of_object(*it);
          GMM_ASSERT1(ito != stored_objects.end(), "An object disapeared !");
          ito->second.valid = false;
          std::set<pstatic_stored_object> dep = ito->second.dependencies;
          for (itd = dep.begin(); itd != dep.end(); ++itd) 
          {
            if (del_dependency(*it, *itd)) 
            {
              stored_object_tab::iterator itod=iterator_of_object(*itd);
              if (itod->second.perm == AUTODELETE_STATIC_OBJECT
                && itod->second.valid) 
              {
                  itod->second.valid = false;
                  to_delete.push_back(*itd);
              }
            }
          }
          for (itd = ito->second.dependent_object.begin();
                          itd != ito->second.dependent_object.end(); ++itd) 
          {
            stored_object_tab::iterator itod=iterator_of_object(*itd);
            if (itod != stored_objects.end()) 
            {
              GMM_ASSERT1(itod->second.perm != PERMANENT_STATIC_OBJECT,
              "Trying to delete a permanent object " << *itd);
              if (itod->second.valid) 
              {
                itod->second.valid = false;
                to_delete.push_back(itod->second.p);
              }
            }
          }
        }
      }
      basic_delete(to_delete);
  }


  void del_stored_object(pstatic_stored_object o, bool ignore_unstored) 
  {
    std::list<pstatic_stored_object> to_delete;
    to_delete.push_back(o);
    del_stored_objects(to_delete, ignore_unstored);
  }


  void del_stored_objects(permanence perm) 
  {
    std::list<pstatic_stored_object> to_delete;
    for(size_t thread=0; thread<getfem::num_threads();thread++)
    {
      stored_object_tab& stored_objects
        = dal::singleton<stored_object_tab>::instance(thread);
      if (perm == PERMANENT_STATIC_OBJECT) perm = STRONG_STATIC_OBJECT;
      stored_object_tab::iterator it;
      for (it = stored_objects.begin(); it != stored_objects.end(); ++it)
        if (it->second.perm >= perm)
          to_delete.push_back(it->second.p);
    }
    del_stored_objects(to_delete, false);
  }

  void list_stored_objects(std::ostream &ost) 
  {
    for(size_t thread=0; thread<getfem::num_threads();thread++)
    {
      stored_object_tab::stored_key_tab& stored_keys 
      = dal::singleton<stored_object_tab>::instance(thread).stored_keys_;

      if (stored_keys.begin() == stored_keys.end())
        ost << "No static stored objects" << endl;
      else ost << "Static stored objects" << endl;
      for (stored_object_tab::stored_key_tab::iterator it = stored_keys.begin();
        it != stored_keys.end(); ++it) 
      {
          ost << "Object: " << it->first << " typename: "
            << typeid(*it->first).name() << endl;
      }
    }
  }

  size_t nb_stored_objects(void) 
  {
    long num_objects=0;
    for(size_t thread=0;thread<getfem::num_threads(); ++thread)
    {
      stored_object_tab::stored_key_tab& stored_keys 
      = dal::singleton<stored_object_tab>::instance(thread).stored_keys_;
      num_objects+=stored_keys.size();
    }
    return num_objects;
  }

}
