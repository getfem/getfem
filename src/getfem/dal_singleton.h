/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2015 Julien Pommier

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

/**@file dal_singleton.h
@author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
@date May 2004.
@brief A simple singleton implementation

Singleton was made thread safe for OpenMP
However, now there is a singleton instance for every 
thread (singleton is thread local). This replicates
the behaviour of singletons in distirbuted MPI-like 
environment;
*/
#ifndef DAL_SINGLETON
#define DAL_SINGLETON

#include <vector>
#include <memory>
#include "getfem_omp.h"


namespace dal {

  class singleton_instance_base {
  public:
    virtual ~singleton_instance_base() {}
    virtual int level() = 0;
  };


  class singletons_manager {
  protected:
    getfem::omp_distribute<std::vector<singleton_instance_base *> > lst;
    static singletons_manager& manager();
    static singletons_manager& m;

  public:
    static void register_new_singleton(singleton_instance_base *p);
    static void register_new_singleton(singleton_instance_base *p, size_t ithread);
    ~singletons_manager();
  private:
    singletons_manager();
  };




  template <typename T, int LEV> class singleton_instance : public singleton_instance_base 
  {
    static getfem::omp_distribute<T*>* instance_;
    static getfem::omp_distribute<T*>* omp_distro_pointer()
    {
      static getfem::omp_distribute<T*>* pointer = new getfem::omp_distribute<T*>( );
      return pointer;
    }
    static T*& instance_pointer() { return omp_distro_pointer( )->thrd_cast(); }
    static T*& instance_pointer(size_t ithread) { return (*omp_distro_pointer( ))(ithread);}

  public:

    singleton_instance() {}

    /** Instance from the current thread*/
    inline static T& instance() 
    { 
      T*& tinstance_ = instance_pointer();
      if (!tinstance_) {
        tinstance_ = new T();
        singletons_manager::register_new_singleton(new singleton_instance<T,LEV>());
      }
      return *tinstance_; 
    }

    /**Instance from thread ithread*/
    inline static T& instance(size_t ithread) 
    { 
      T*& tinstance_ = instance_pointer(ithread);
      if (!tinstance_) {
        tinstance_ = new T();
        singletons_manager::register_new_singleton(new singleton_instance<T,LEV>(),ithread);
      }
      return *tinstance_; 
    }

    int level() { return LEV; }

    ~singleton_instance() 
    {
      if (instance_) 
      {
        for(size_t i=0;i<getfem::num_threads();i++)
        {
          if((*instance_)(i)) 
          { 
            delete (*instance_)(i); 
            (*instance_)(i) = 0; 
          }
        } 
      }
      delete instance_; instance_=0;
    }
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

    /** Instance from the current thread*/
    inline static T& instance() { 
      return singleton_instance<T,LEV>::instance();
    }
    inline static const T& const_instance() { return instance(); }

    inline static T& instance(size_t ithread) { 
      return singleton_instance<T,LEV>::instance(ithread);
    }
    inline static const T& const_instance(size_t ithread) { return instance(ithread); }


  protected:
    singleton() {}
    ~singleton() {}
  private:
    singleton(const singleton&);            
    singleton& operator=(const singleton&);
  };

  template <typename T, int LEV> 
  getfem::omp_distribute<T*>* singleton_instance<T,LEV>::instance_
                        = singleton_instance<T,LEV>::omp_distro_pointer();
}

#endif


