// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_shared_ptr.h : shared pointers
//           
// Date    : June 2004.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2006 Julien Pommier
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

/**@file dal_shared_ptr.h
   @brief A (very simplified and rough) version of boost::shared_ptr.
*/
#ifndef DAL_SHARED_PTR_H__
#define DAL_SHARED_PTR_H__

#include<dal_std.h>



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
    void reset(T *q, bool refcounted = true) { release(); shared_ptr<T> tmp(q, refcounted); (*this).swap(tmp); }
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
