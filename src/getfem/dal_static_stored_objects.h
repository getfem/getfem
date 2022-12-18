/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2002-2020 Yves Renard

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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/** @file dal_static_stored_objects.h
@author  Yves Renard <Yves.Renard@insa-lyon.fr>
@date February 19, 2005
@brief Stores interdependent getfem objects.

Stored object :

A type of object to be stored should derive from
dal::static_stored_object and a key should inherit from
static_stored_object_key with an overloaded "compare" method.

To store a new object, you have to test if the object is not
already stored and then call dal::add_stored_object:
@code
pstatic_stored_object_key p = make_shared<your_object_key>(parameters);
if (!search_stored_object(p)) {
add_stored_object(p, make_shared<your_object>(parameters));
}
@endcode
You can add a dependency of your new object with
@code
add_dependency(pointer_on_your_object,
pointer_on_the_object_object_from_which_it_depends);
@endcode
and then your object will be automatically deleted if the second object is
deleted.
The dependency can be added within the add_stored_object call:
@code
add_stored_object(new your_object_key(parameters),
new your_object(parameters),
dependency);
@endcode

std::shared_ptr are used.
*/
#ifndef DAL_STATIC_STORED_OBJECTS_H__
#define DAL_STATIC_STORED_OBJECTS_H__

#include "dal_config.h"
#include "getfem_omp.h"
#include <algorithm>
#include "dal_singleton.h"
#include <set>
#include <list>


#include "getfem/getfem_arch_config.h"

#include <atomic>

#define DAL_STORED_OBJECT_DEBUG 0

namespace dal {


#if DAL_STORED_OBJECT_DEBUG
// Little tool for debug : detects the deleted static stored objects for
// which the destructor is not called (i.e. a shared pointer is still stored
// somewhere. Not thread safe.
// Rule : Each potential stored object should call
// DAL_STORED_OBJECT_DEBUG_CREATED(this, "name") in the constructor and
// DAL_STORED_OBJECT_DEBUG_DESTROYED(this) in the destructor.
  class static_stored_object;
  void stored_debug_created(const static_stored_object *o,
                            const std::string &name);
  void stored_debug_added(const static_stored_object *o);
  void stored_debug_deleted(const static_stored_object *o);
  void stored_debug_destroyed(const static_stored_object *o,
                              const std::string &name);
# define DAL_STORED_OBJECT_DEBUG_CREATED(o, name) stored_debug_created(o, name)
# define DAL_STORED_OBJECT_DEBUG_ADDED(o)   stored_debug_added(o)
# define DAL_STORED_OBJECT_DEBUG_DELETED(o) stored_debug_deleted(o)
# define DAL_STORED_OBJECT_DEBUG_DESTROYED(o, name) \
                                            stored_debug_destroyed(o, name)
#else
# define DAL_STORED_OBJECT_DEBUG_CREATED(o, name)
# define DAL_STORED_OBJECT_DEBUG_ADDED(o)
# define DAL_STORED_OBJECT_DEBUG_DELETED(o)
# define DAL_STORED_OBJECT_DEBUG_DESTROYED(o, name)
#endif

  enum permanence { PERMANENT_STATIC_OBJECT = 0, // not deletable object
    STRONG_STATIC_OBJECT = 1,    // preferable not to delete it
    STANDARD_STATIC_OBJECT = 2,  // standard
    WEAK_STATIC_OBJECT = 3,      // delete it if necessary
    AUTODELETE_STATIC_OBJECT = 4 // automatically deleted
                                 // when the last dependent object is deleted
  };

  class static_stored_object_key {
  protected :
    virtual bool compare(const static_stored_object_key &) const = 0;
    virtual bool equal(const static_stored_object_key &) const = 0;

  public :
    bool operator < (const static_stored_object_key &o) const {
      // comparaison des noms d'objet
      if (typeid(*this).before(typeid(o))) return true;
      if (typeid(o).before(typeid(*this))) return false;
      return compare(o);
    }

    bool operator == (const static_stored_object_key &o) const {
      if (typeid(o)!=typeid(*this)) return false;
      return equal(o);
    }

    bool operator != (const static_stored_object_key &o) const {
      return !(*this == o);
    }

    virtual ~static_stored_object_key() {}
  };

  template <typename var_type>
  class simple_key : virtual public static_stored_object_key {
    var_type a;
  public :
     bool compare(const static_stored_object_key &oo) const override {
      auto &o = dynamic_cast<const simple_key &>(oo);
      return a < o.a;
    }

    bool equal(const static_stored_object_key &oo) const override {
      auto &o = dynamic_cast<const simple_key &>(oo);
      return a == o.a;
    }
    simple_key(var_type aa) : a(aa) {}
  };

#define DAL_SIMPLE_KEY(class_name, var_type)                            \
  struct class_name : public dal::simple_key<var_type> {                \
  class_name(var_type aa) : dal::simple_key<var_type>(aa) {}            \
  }

#define DAL_DOUBLE_KEY(class_name, var_type1, var_type2)                \
  struct class_name :                                                   \
  public dal::simple_key<std::pair<var_type1,var_type2> > {             \
  class_name(var_type1 aa, var_type2 bb) :                              \
  dal::simple_key<std::pair<var_type1,var_type2> >                      \
  (std::make_pair(aa,bb)) {}                                            \
  }

#define DAL_TRIPLE_KEY(class_name, var_type1, var_type2, var_type3)     \
  struct class_name :                                                   \
  public dal::simple_key<std::pair<var_type1,                           \
                                   std::pair<var_type2,var_type3> > > { \
  class_name(var_type1 aa, var_type2 bb, var_type3 cc) :                \
    dal::simple_key<std::pair<var_type1,                                \
                              std::pair<var_type2, var_type3> > >       \
    (std::make_pair(aa,std::make_pair(bb,cc))) {}                       \
  }

#define DAL_FOUR_KEY(class_name,var_type1,var_type2,var_type3,var_type4)\
  struct class_name : public                                            \
  dal::simple_key<std::pair                                             \
                  <var_type1, std::pair<var_type2, std::pair            \
                                        <var_type3,var_type4> > > > {   \
    class_name(var_type1 aa, var_type2 bb, var_type3 cc,var_type4 dd) : \
      dal::simple_key<std::pair                                         \
                      <var_type1, std::pair<var_type2,                  \
                                            std::pair<var_type3,        \
                                                      var_type4> > > >  \
      (std::make_pair(aa,std::make_pair(bb,std::make_pair(cc, dd)))) {} \
  }


  typedef std::shared_ptr<const static_stored_object_key>
    pstatic_stored_object_key;


  /**
  base class for static stored objects
  */
  class static_stored_object { public : virtual ~static_stored_object() {} };

  typedef std::shared_ptr<const static_stored_object> pstatic_stored_object;

  pstatic_stored_object_key key_of_stored_object(pstatic_stored_object o);

  /** Gives a pointer to an object from a key pointer. */
  pstatic_stored_object search_stored_object(pstatic_stored_object_key k);

  pstatic_stored_object search_stored_object_on_all_threads(pstatic_stored_object_key k);

  /** Test if an object is stored*/
  bool exists_stored_object(pstatic_stored_object o);

  /** Add a dependency, object o1 will depend on object o2. */
  void add_dependency(pstatic_stored_object o1, pstatic_stored_object o2);

  /** remove a dependency. Return true if o2 has no more dependent object. */
  bool del_dependency(pstatic_stored_object o1, pstatic_stored_object o2);

  /** Add an object with two optional dependencies. */
  void add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
                         permanence perm = STANDARD_STATIC_OBJECT);

  inline void
    add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
    pstatic_stored_object dep1,
    permanence perm = STANDARD_STATIC_OBJECT) {
      add_stored_object(k, o, perm);
      add_dependency(o, dep1);
  }

  inline void
  add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
    pstatic_stored_object dep1, pstatic_stored_object dep2,
    permanence perm = STANDARD_STATIC_OBJECT) {
      add_stored_object(k, o, perm);
      add_dependency(o, dep1);
      add_dependency(o, dep2);
  }

  inline void
  add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
                    pstatic_stored_object dep1, pstatic_stored_object dep2,
                    pstatic_stored_object dep3,
                    permanence perm = STANDARD_STATIC_OBJECT) {
      add_stored_object(k, o, perm);
      add_dependency(o, dep1);
      add_dependency(o, dep2);
      add_dependency(o, dep3);
  }

  inline void
  add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
                    pstatic_stored_object dep1, pstatic_stored_object dep2,
                    pstatic_stored_object dep3, pstatic_stored_object dep4,
                    permanence perm = STANDARD_STATIC_OBJECT) {
      add_stored_object(k, o, perm);
      add_dependency(o, dep1);
      add_dependency(o, dep2);
      add_dependency(o, dep3);
      add_dependency(o, dep4);
  }

  /** Delete an object and the object which depend on it. */
  void del_stored_object(const pstatic_stored_object &o,
                         bool ignore_unstored=false);

  /** Delete all the object whose permanence is greater or equal to perm. */
  void del_stored_objects(int perm);

  /** Show a list of stored objects (for debugging purpose). */
  void list_stored_objects(std::ostream &ost);

  /** Return the number of stored objects (for debugging purpose). */
  size_t nb_stored_objects(void);

  /** Delete a list of objects and their dependencies*/
  void del_stored_objects(std::list<pstatic_stored_object> &to_delete,
    bool ignore_unstored);

  /** Test the validity of the whole global storage */
  void test_stored_objects(void);


  /** Pointer to an object with the dependencies */
  struct enr_static_stored_object {
    pstatic_stored_object p;
    std::atomic_bool valid;
    const permanence perm;
    std::set<pstatic_stored_object> dependent_object;
    std::set<pstatic_stored_object> dependencies;
    enr_static_stored_object(pstatic_stored_object o, permanence perma)
      : p(o), perm(perma) {valid = true;}
    enr_static_stored_object()
      : perm(STANDARD_STATIC_OBJECT) {valid = true;}
    enr_static_stored_object(const enr_static_stored_object& enr_o)
      : p(enr_o.p), perm(enr_o.perm), dependent_object(enr_o.dependent_object),
      dependencies(enr_o.dependencies){valid = static_cast<bool>(enr_o.perm);}
  };



  /** Pointer to a key with a coherent order */
  struct enr_static_stored_object_key {
    pstatic_stored_object_key p;
    bool operator < (const enr_static_stored_object_key &o) const
    { return (*p) < (*(o.p)); }
    enr_static_stored_object_key(pstatic_stored_object_key o) : p(o) {}
  };



  /** Table of stored objects. Thread safe, uses thread specific mutexes. */
  struct stored_object_tab :
    public std::map<enr_static_stored_object_key, enr_static_stored_object> {

    typedef std::map<pstatic_stored_object,pstatic_stored_object_key>
      stored_key_tab;

    stored_object_tab();
    ~stored_object_tab();
    pstatic_stored_object
      search_stored_object(pstatic_stored_object_key k) const;
    bool has_dependent_objects(pstatic_stored_object o) const;
    bool exists_stored_object(pstatic_stored_object o) const;
    //adding the object to the storage on the current thread
    void add_stored_object(pstatic_stored_object_key k, pstatic_stored_object o,
    permanence perm);

    iterator iterator_of_object_(pstatic_stored_object o);
    //delete o2 from the dependency list of o1
    //true if successfull, false if o1 is not
    //on this thread
    bool del_dependency_(pstatic_stored_object o1,
    pstatic_stored_object o2);
    //delete o1 from the dependent list of o2
    //true if successfull, false if o1 is not
    //on this thread
    bool del_dependent_(pstatic_stored_object o1,
    pstatic_stored_object o2);
    //add o2 to the dependency list of o1
    //true if successfull, false if o1 is not
    //on this thread
    bool add_dependency_(pstatic_stored_object o1,
    pstatic_stored_object o2);
    //add o1 to the dependent list of o2
    //true if successfull, false if o1 is not
    //on this thread
    bool add_dependent_(pstatic_stored_object o1,
    pstatic_stored_object o2);
    void basic_delete_(std::list<pstatic_stored_object> &to_delete);

    getfem::lock_factory locks_;
    stored_key_tab stored_keys_;
  };



  /** delete all the specific type of stored objects*/
  template<typename OBJECT_TYPE>
  void delete_specific_type_stored_objects(bool all_threads = false){
    std::list<pstatic_stored_object> delete_object_list;

    auto filter_objects = [&](stored_object_tab &stored_objects){
      for(auto &&pair : stored_objects){
        auto p_object = std::dynamic_pointer_cast<const OBJECT_TYPE>(pair.second.p);
        if(p_object != nullptr) delete_object_list.push_back(pair.second.p);
      }
    };

    if (!all_threads){
      auto& stored_objects = singleton<stored_object_tab>::instance();
      filter_objects(stored_objects);
    }
    else{
      for(size_t thread = 0; thread < singleton<stored_object_tab>::num_threads(); ++thread)
      {
        auto& stored_objects = singleton<stored_object_tab>::instance(thread);
        filter_objects(stored_objects);
      }
    }
    del_stored_objects(delete_object_list, false);
  }
}

#endif /* DAL_STATIC_STORED_OBJECTS_H__ */
