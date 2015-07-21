/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2015 Julien Pommier
 
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

/**@file dal_shared_ptr.h
   @author  Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date June 2004.
   @brief A (very simplified and rough) version of boost::shared_ptr.
*/
#ifndef DAL_SHARED_PTR_H__
#define DAL_SHARED_PTR_H__

#include "dal_config.h"

namespace dal {
  /**
     a (very simplified and rough) version of boost::shared_ptr
     
     - not thread safe.
     - no fancy operators
     - possibility to handle non-refcounted pointers (i.e. normal pointers)

     caution: should not be used for data allocated with new[] 
  */
  template <typename T> class shared_ptr {
    T *p;
    unsigned long *refcnt;
  public:
    shared_ptr() : p(0), refcnt(0) {}
    explicit shared_ptr(T *q, bool refcounted = true) : p(q), refcnt(0) { if (refcounted) refcnt = new unsigned long(1); }
    shared_ptr(const shared_ptr<T> &other) : p(other.p), refcnt(other.refcnt) { if (refcnt) ++(*refcnt); }
    void reset(T *q = NULL, bool refcounted = true) { release(); if (q) {shared_ptr<T> tmp(q, refcounted); (*this).swap(tmp);} }
    void swap(shared_ptr<T> &other) { std::swap(p,other.p); std::swap(refcnt,other.refcnt); }
    void release() { if (refcnt && --(*refcnt) == 0) { if (p) delete p; delete refcnt; } p = 0; refcnt = 0; }
    ~shared_ptr() { release(); }
    shared_ptr<T> &operator=(const shared_ptr<T> &other) { shared_ptr<T> tmp(other); swap(tmp); return *this; }
    T *get() const { return p; }
    T& operator*() const { return *p; }
    T* operator->() const { return p; }
    bool counted() { return refcnt; }
    unsigned long use_count() const { return refcnt ? 0 : *refcnt; }
  };

  /**
     shared_array uses operator delete[] for destruction of data
   */
  template <typename T> class shared_array {
    T *p;
    unsigned long *refcnt;
  public:
    shared_array() : p(0), refcnt(0) {}
    explicit shared_array(T *q, bool refcounted = true) : p(q), refcnt(0) { if (refcounted) refcnt = new unsigned long(1); }
    shared_array(const shared_array<T> &other) : p(other.p), refcnt(other.refcnt) { if (refcnt) ++(*refcnt); }
    void reset(T *q, bool refcounted = true) { release(); shared_array<T> tmp(q,refcounted); (*this).swap(tmp); }
    void swap(shared_array<T> &other) { std::swap(p,other.p); std::swap(refcnt,other.refcnt); }
    void release() { if (refcnt && --(*refcnt) == 0) { if (p) delete[] p; delete refcnt; } p = 0; refcnt = 0; }
    ~shared_array() { release(); }
    shared_array<T> &operator=(const shared_array<T> &other) { shared_array<T> tmp(other); swap(tmp); return *this; }
    T *get() const { return p; }
    bool counted() { return refcnt; }
    unsigned long use_count() const { return refcnt ? 0 : *refcnt; }
    T& operator[](unsigned i) const { return p[i]; }
  };
}


#endif
