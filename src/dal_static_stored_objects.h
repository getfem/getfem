// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_static_stored_objects.h : object which should be stored.
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


#ifndef DAL_STATIC_STORED_OBJECTS_H__
#define DAL_STATIC_STORED_OBJECTS_H__

#include <dal_except.h>

namespace dal {


  // Stored object :  
  // A type of object to be stored should derive from static_stored_object
  // and a key should derive from  static_stored_object_key with an overloaded
  // method compare.
  // 
  // To store a new object, you have to test if the object is not
  // already stored and then to call add_stored_object:
  // 
  // if (!search_stored_object(your_object_key(parameters))) {
  //   add_stored_object(new your_object_key(parameters),
  //                     new your_object(parameters));
  // }
  // You can add a dependency of your new object with
  // add_dependency(pointer_on_your_object,
  //                pointer_on_the_object_object_from_which_it_depends,
  //                perm)
  // your object will be automatically deleted if the second object is
  // deleted.
  // The dependency can be added when the addition is done by
  // add_stored_object(new your_object_key(parameters), 
  //                   new your_object(parameters),
  //                   dependence)


  class static_stored_object_key {
  protected :
    virtual bool compare(const static_stored_object_key &o) const {
      DAL_THROW(failure_error, "This method should not be called");
    }
    
  public :
    bool operator < (const static_stored_object_key &o) const {
      // comparaison des noms d'objet
      if (typeid(*this).before(typeid(o))) return true;
      if (typeid(o).before(typeid(*this))) return false;
      return compare(o);
    }
    
    
    virtual ~static_stored_object_key() {}
    
  };
  
  class static_stored_object {
    
    
  public :
    virtual ~static_stored_object() {}
  };

  /** Gives a pointer to an object from a key pointer. */
  const static_stored_object *
  search_stored_object(const static_stored_object_key *k);

  /** Gives a pointer to an object from a key reference. */
  inline const static_stored_object *
  search_stored_object(const static_stored_object_key &k)
  { return search_stored_object(&k); }

  /** Add a dependency, object o1 will depend on object o2. */
  void add_dependency(const static_stored_object *o1,
		      const static_stored_object *o2);

  /** remove a dependency. */
  void del_dependency(const static_stored_object *o1,
		      const static_stored_object *o2);

  /** Add an object with two optional dependencies. */
  void add_stored_object(const static_stored_object_key *k,
			 const static_stored_object *o,
			 int permanence = 2,
			 const static_stored_object *dep1 = 0,
			 const static_stored_object *dep2 = 0);

  /** Delete an object and its dependencies. */
  void del_stored_object(const static_stored_object *o);
  
  /** Delete all the object whose permanence is greater or equal to perm. */
  void del_stored_objects(int perm);

}

#endif /* DAL_STATIC_STORED_OBJECTS_H__ */
