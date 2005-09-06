// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_singleton.h : basic template class for singletons
//           
// Date    : May 2004.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Julien Pommier
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

/**@file dal_singleton.h
   @brief A simple singleton implementation

   Not thread safe, of course.
*/
#ifndef DAL_SINGLETON
#define DAL_SINGLETON

#include <vector>
#include <memory>

namespace dal {

  class singleton_instance_base {
  public:
    virtual ~singleton_instance_base() {}
    virtual int level() = 0;
  };

  class singletons_manager {
  protected:
    std::vector<singleton_instance_base *> lst;  
    static std::auto_ptr<singletons_manager> m;
  public:
    static void register_new_singleton(singleton_instance_base *p);
    ~singletons_manager();
  private:
    singletons_manager() {}
  };
  
  template <typename T, int LEV> class singleton_instance : public singleton_instance_base {
  public:
    static T *instance_;
    inline static T& instance() { 
      if (!instance_) {
	instance_ = new T(); 
	singletons_manager::register_new_singleton(new singleton_instance<T,LEV>());
      }
      return *instance_; 
    }
    int level() { return LEV; }
    singleton_instance() {}
    ~singleton_instance() { if (instance_) { delete instance_; instance_ = 0; } }
  };

  /** singleton class. 

     usage: 
     @code
     foo &f = singleton<foo>::instance();
     const foo &f = singleton<foo>::const_instance();
     @endcode
     the LEV template arguments allows one to choose the order of destruction
     of the singletons:
       lowest LEV will be destroyed first.
  */
  template <typename T, int LEV=1> class singleton {
  public:
    inline static T& instance() { 
      return singleton_instance<T,LEV>::instance();
    }
    inline static const T& const_instance() { return instance(); }
  protected:
    singleton() {}
    ~singleton() {}
  private:
    singleton(const singleton&);            
    singleton& operator=(const singleton&);
  };

  template <typename T, int LEV> T* singleton_instance<T,LEV>::instance_ = 0;
}

#endif
