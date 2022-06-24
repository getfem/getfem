/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2020 Julien Pommier

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

/**@file dal_singleton.h
@author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
@date May 2004.
@brief A simple singleton implementation

Singleton was made thread safe for OpenMP
However, now there is a singleton instance for every
thread (singleton is thread local). This replicates
the behaviour of singletons in distributed MPI-like
environment;
*/

#pragma once

#include <vector>
#include <memory>

#include "getfem_omp.h"


namespace dal {

  using bgeot::size_type;

  class singleton_instance_base {
  public:
    virtual ~singleton_instance_base() {};
    virtual int level() const = 0;
  };

  class singletons_manager {
    getfem::omp_distribute<std::vector<singleton_instance_base *>> lst;
    size_type nb_partitions;
    static singletons_manager& manager();

  public:
    static void register_new_singleton(singleton_instance_base *p);
    static void register_new_singleton(singleton_instance_base *p,
                                       size_t ithread);
    static void on_partitions_change();

    /**destroy singletons in increasing order*/
    ~singletons_manager();

  private:
    singletons_manager();
  };

  template <typename T, int LEV>
  class singleton_instance : public singleton_instance_base {

    static getfem::omp_distribute<T*>* initializing_pointer;

    static getfem::omp_distribute<T*>*& pointer() {
      static auto p = new getfem::omp_distribute<T*>{};
      return p;
    }

    static T*& instance_pointer() {
      return pointer()->thrd_cast();
    }

    static T*& instance_pointer(size_t ithread) {
      return (*pointer())(ithread);
    }

  public:

    /**Instance from thread ithread*/
    inline static T& instance(size_t ithread) {
      pointer()->on_thread_update();
      T*& tinstance_ = instance_pointer(ithread);
      if (!tinstance_) {
        tinstance_ = new T();
        singletons_manager::register_new_singleton(
          new singleton_instance<T,LEV>(), ithread);
      }
      return *instance_pointer(ithread);
    }

    /** Instance from the current thread*/
    inline static T& instance() {
      return instance(this_thread());
    }

    inline static size_type num_threads() {
      return pointer()->num_threads();
    }

    inline static size_type this_thread() {
      return pointer()->this_thread();
    }

    int level() const override {
      return LEV;
    }

    ~singleton_instance() {
      if (!pointer()) return;
      for(size_t i = 0; i != pointer()->num_threads(); ++i) {
        auto &p_singleton = (*pointer())(i);
        if(p_singleton){
          delete p_singleton;
          p_singleton = nullptr;
        }
      }
      delete pointer();
      pointer() = nullptr;
      if (initializing_pointer) initializing_pointer = nullptr;
    }
  };

  template<typename T, int LEV> getfem::omp_distribute<T*>*
  singleton_instance<T, LEV>::initializing_pointer = singleton_instance<T, LEV>::pointer();

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

    singleton(const singleton&) = delete;
    singleton& operator=(const singleton&) = delete;

    /** Instance from the current thread*/
    inline static T& instance() {
      return singleton_instance<T,LEV>::instance();
    }

    inline static const T& const_instance() {
      return instance();
    }

    inline static T& instance(size_t ithread) {
      return singleton_instance<T,LEV>::instance(ithread);
    }

    inline static const T& const_instance(size_t ithread){
      return instance(ithread);
    }

    /** number of threads this singleton is distributed on.*/
    inline static size_type num_threads(){
      return singleton_instance<T,LEV>::num_threads();
    }

    /** this thread number according to the threading policy of the singleton*/
    inline static size_type this_thread() {
      return singleton_instance<T, LEV>::this_thread();
    }

  protected:
    singleton() = default;
    ~singleton() = default;
  };

}/* end of namespace dal                                                             */