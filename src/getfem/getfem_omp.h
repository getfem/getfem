/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

Copyright (C) 2000-2012 Yves Renard

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

As a special exception, you  may use  this file  as it is a part of a free
software  library  without  restriction.  Specifically,  if   other  files
instantiate  templates  or  use macros or inline functions from this file,
or  you compile this  file  and  link  it  with other files  to produce an
executable, this file  does  not  by itself cause the resulting executable
to be covered  by the GNU Lesser General Public License.  This   exception
does not  however  invalidate  any  other  reasons why the executable file
might be covered by the GNU Lesser General Public License.

===========================================================================*/
/**@file getfem_omp.h
@author  Andriy Andreykiv <andriy.andreykiv@gmail.com>
@date May 14th, 2013.
@brief Tools for multithreaded, OpenMP and Boost based parallelization.

This is the kernel of getfem.
*/
#pragma once
#ifndef GETFEM_OMP
#define GETFEM_OMP

#include <vector>
#include <algorithm>

#ifdef GETFEM_HAVE_OPENMP
#include <omp.h>
#include <boost/thread.hpp>
#endif

#ifdef _WIN32
#include <mbctype.h>
#endif

#include <locale.h>
#include <memory>
#include "dal_shared_ptr.h"
#include "gmm/gmm_std.h"
#include "bgeot_config.h"
#ifdef GETFEM_HAVE_OPENMP
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#endif


namespace getfem
{
  using bgeot::size_type;
  class mesh;
  class mesh_region;


#ifdef GETFEM_HAVE_OPENMP
  //declaring a thread lock, to protect multi-threaded accesses to
  //asserts, traces and warnings. Using a global mutex
  class omp_guard: public boost::lock_guard<boost::recursive_mutex>
  {
  public:
    omp_guard();
  private:
    static boost::recursive_mutex boost_mutex;
  };

  //like boost::lock_guard, but copyable
  class local_guard 
  {
  public:
    local_guard(boost::recursive_mutex&);
    local_guard(const local_guard&);

  private:
    boost::recursive_mutex& mutex_;
    boost::shared_ptr<boost::lock_guard<boost::recursive_mutex>> plock_;
  };

  //produces scoped lock on the 
  //mutex, held in this class
  class lock_factory
  {
  public:
    lock_factory();

    //get a lock object with RAII acuire/release semantics
    //on the mutex from this factory
    local_guard get_lock() const;
  private:
    mutable boost::recursive_mutex mutex_;
  };


#else 

  class omp_guard{};
  class local_guard{};
  struct lock_factory
  {
    inline local_guard get_lock() const {return local_guard();}
  };
#endif


#ifdef GETFEM_HAVE_OPENMP	
  /**number of OpenMP threads*/
  inline size_t num_threads(){return omp_get_max_threads();}

  /**set maximum number of OpenMP threads*/
  inline void set_num_threads(int n){omp_set_num_threads(n);}

  /**index of the current thread*/
  inline size_type this_thread() {return omp_get_thread_num();}
  /**is the program running in the parallel section*/
  inline bool me_is_multithreaded_now(){return static_cast<bool>(omp_in_parallel());}

#else
  inline size_type num_threads(){return size_type(1);}
  inline void set_num_threads(int n) { }
  inline size_type this_thread() {return size_type(0);}
  inline bool me_is_multithreaded_now(){return false;}
#endif



  /**use this template class for any object you want to 
  distribute to open_MP threads. The creation of this
  object should happen in serial, while accessing the individual
  thread local instances will take place in parallel. If 
  one needs creation of thread local object, use the macro
  DEFINE_STATIC_THREAD_LOCAL
  */
  template <typename T> class omp_distribute {
    std::vector<T> thread_values;
    friend struct all_values_proxy;
    struct all_values_proxy{
      omp_distribute& distro;
      all_values_proxy(omp_distribute& d): distro(d){}
      void operator=(const T& x){
        for(typename std::vector<T>::iterator it=distro.thread_values.begin();
          it!=distro.thread_values.end();it++) *it=x;
      }
    };
  public:
    omp_distribute() : thread_values(num_threads()) {}
    omp_distribute(const T& value) : 
      thread_values(num_threads(),value) {}
    operator T& (){return thread_values[this_thread()];}
    operator const T& () const {return thread_values[this_thread()];}
    T& thrd_cast(){return thread_values[this_thread()];}
    const T& thrd_cast() const {return thread_values[this_thread()];}
    T& operator()(size_type i) {
      return thread_values[i];
    }	
    const T& operator()(size_type i) const {
      return thread_values[i];
    }
    T& operator = (const T& x){ 
      return (thread_values[this_thread()]=x);
    }

    all_values_proxy all_threads(){return all_values_proxy(*this);}

    ~omp_distribute(){}
  };

  template <typename T> class omp_distribute<std::vector<T> > {
    std::vector<std::vector<T> > thread_values;
    friend struct all_values_proxy;
    struct all_values_proxy{
      omp_distribute& distro;
      all_values_proxy(omp_distribute& d): distro(d){}
      void operator=(const T& x){
        for(typename std::vector<T>::iterator it=distro.thread_values.begin();
          it!=distro.thread_values.end();it++) *it=x;
      }
    };

  public:
    typedef std::vector<T> VEC;
    omp_distribute() : thread_values(num_threads()) {}
    omp_distribute(size_t n, const T& value) : 
      thread_values(num_threads(), std::vector<T>(n,value)){}
    operator VEC& (){return thread_values[this_thread()];}
    operator const VEC& () const 
    {return thread_values[this_thread()];}
    VEC& operator()(size_type i) {return thread_values[i];}	
    const VEC& operator()(size_type i) const {return thread_values[i];}
    VEC& thrd_cast(){return thread_values[this_thread()];}
    const VEC& thrd_cast() const 
    {return thread_values[this_thread()];}
    T& operator[](size_type i) 
    {return thread_values[this_thread()][i];}	
    const T& operator[](size_type i) const 
    {return thread_values[this_thread()][i];}
    T& operator = (const T& value) {
      return (thread_values[this_thread()]=value);
    }
    all_values_proxy all_threads(){return all_values_proxy(*this);}
    ~omp_distribute(){}
  };

  /**specialization for bool, to circumvent the shortcommings
  of standards library specialization for std::vector<bool>*/
  template <> class omp_distribute<bool> {
    typedef int BOOL;
    std::vector<BOOL> thread_values;
    friend struct all_values_proxy;
    struct all_values_proxy{
      omp_distribute<bool>& distro;
      all_values_proxy(omp_distribute& d): distro(d){}
      void operator=(const bool& x);
    };


  public:

    omp_distribute() : thread_values(num_threads()) {}
    omp_distribute(const bool& value) : 
      thread_values(num_threads(),value) {}
    operator BOOL& (){return thread_values[this_thread()];}
    operator const BOOL& () const {return thread_values[this_thread()];}
    BOOL& thrd_cast(){return thread_values[this_thread()];}
    const BOOL& thrd_cast() const {return thread_values[this_thread()];}
    BOOL& operator()(size_type i) {
      return thread_values[i];
    }	
    const BOOL& operator()(size_type i) const {
      return thread_values[i];
    }
    BOOL& operator = (const BOOL& x){ 
      return (thread_values[this_thread()]=x);
    }
    all_values_proxy all_threads(){return all_values_proxy(*this);}
    ~omp_distribute(){}
  private:

  };

  /* Use these macros only in function local context to achieve
  the effect of thread local storage for any type of objects
  and their initialization (it's more general and portable
  then using __declspec(thread))*/
#ifdef GETFEM_HAVE_OPENMP
#define	DEFINE_STATIC_THREAD_LOCAL_INITIALIZED(Type,Var,initial) \
  static boost::thread_specific_ptr<Type> ptr_##Var; \
  if(!ptr_##Var.get()) {ptr_##Var.reset(new Type(initial));} \
  Type& Var=*ptr_##Var;

#define	DEFINE_STATIC_THREAD_LOCAL(Type,Var) \
  static boost::thread_specific_ptr<Type> ptr_##Var; \
  if(!ptr_##Var.get()) {ptr_##Var.reset(new Type());} \
  Type& Var=*ptr_##Var;

#define DEFINE_STATIC_THREAD_LOCAL_CONSTRUCTED(Type, Var, ...) \
  static boost::thread_specific_ptr<Type> ptr_##Var; \
  if(!ptr_##Var.get()) {ptr_##Var.reset(new Type(__VA_ARGS__));} \
  Type& Var=*ptr_##Var;

#else
#define	DEFINE_STATIC_THREAD_LOCAL_INITIALIZED(Type,Var,initial) \
  static Type Var(initial);

#define	DEFINE_STATIC_THREAD_LOCAL(Type,Var) \
  static Type Var;

#define DEFINE_STATIC_THREAD_LOCAL_CONSTRUCTED(Type, Var, ...) \
  static Type Var(__VA_ARGS__);

#endif

  class open_mp_is_running_properly{
    static omp_distribute<bool> answer;

  public:
    open_mp_is_running_properly();
    ~open_mp_is_running_properly();
    static bool is_it();
  };

#if defined _WIN32 && !defined (__GNUC__)
  /**parallelization function for a for loop*/
  template<class LOOP_BODY> 
  inline void open_mp_for(int begin, int end, 
    const LOOP_BODY& loop_body){
      _configthreadlocale(_ENABLE_PER_THREAD_LOCALE); 
      gmm::standard_locale locale;
      open_mp_is_running_properly check;
#pragma omp parallel default(shared) 
      { 
        _setmbcp(_MB_CP_ANSI);
#pragma omp for schedule(static)
        for(int i=begin;i<end;i++) loop_body(i);
      }
      _configthreadlocale(_DISABLE_PER_THREAD_LOCALE);
  }
#else /*LINUX*/
  /**parallelization function for a for loop*/
  template<class LOOP_BODY> 
  inline void open_mp_for(
    int begin, int end, const LOOP_BODY& loop_body){
      gmm::standard_locale locale;
      open_mp_is_running_properly check;
#pragma omp parallel default(shared) 
      { 
#pragma omp for schedule(static)
        for(int i=begin;i<end;i++) loop_body(i);
      }
  }
#endif








  /**parallelization macro of a for loop*/
#define OPEN_MP_FOR(begin,end,loop_counter,loop_body) \
  getfem::open_mp_for(begin,end,loop_body(loop_counter));

  /**used to partition a mesh region so that 
  each partition can be used on a different thread. Thread safe*/
  class region_partition {
    mesh* pparent_mesh;
    dal::shared_ptr<mesh_region> original_region;
    mutable std::vector<size_type> partitions;
  public:
    region_partition(mesh* mmesh=0,size_type id=-1);
    region_partition(const region_partition& rp);
    void operator=(const region_partition& rp);		
    size_type thread_local_partition() const;
  };


}

#endif //GETFEM_OMP

