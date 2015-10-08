/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2012-2015 Andriy Andreykiv

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
/**@file getfem_copyable_ptr.h
@author  Andriy Andreykiv <andriy.andreykiv@gmail.com>
@date October 6th, 2015.
@brief A smart pointer that copies the value it points to on copy operations
*/
#pragma once

#include <memory>

namespace getfem {

/**
  A wrapper around a unique_ptr that clones the value on copy
*/
template<class T> class copyable_ptr{
  std::unique_ptr<T> p_ = nullptr;
public:
  copyable_ptr() = default;

  copyable_ptr(std::unique_ptr<T> p) : p_(std::move(p)) {}

  copyable_ptr(const copyable_ptr<T> &x) : p_(x.p_ ? std::unique_ptr<T>(new T(*x.p_)) : nullptr)  {}

  copyable_ptr(copyable_ptr<T> &&x) : p_(std::move(*x.p_))  {}

  copyable_ptr<T> &operator=(const copyable_ptr<T> &x){
    if (x.p_) p_ = std::unique_ptr<T>(new T(*x.p_));
    return *this;
  }

  copyable_ptr<T> &operator=(copyable_ptr<T> &&x){
    p_ = std::move(x.p_);
    return *this;
  }

  operator bool() const {return p_.operator bool();}

  T& operator*() const {return *p_;}

  T& operator->() const {return *p_;}
};

} // namespace getfem