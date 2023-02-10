/*===========================================================================

 Copyright (C) 2004-2020 Yves Renard

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

===========================================================================*/



#include "getfem/getfem_context.h"

namespace getfem {

  context_dependencies::context_dependencies(const context_dependencies &cd)
    : state(cd.state),
      touched(static_cast<bool>(cd.touched)),
      dependencies(cd.dependencies),
      dependent(cd.dependent),
      locks_( )
  {}

  context_dependencies&
  context_dependencies::operator=(const context_dependencies &cd) {
    state = cd.state;
    touched = static_cast<bool>(cd.touched);
    dependencies = cd.dependencies;
    dependent = cd.dependent;
    return *this;
  }

  void context_dependencies::sup_dependent_
  (const context_dependencies &cd) const {
    getfem::local_guard lock = locks_.get_lock();
    size_type s = dependent.size();
    iterator_list it1 = dependent.begin(), it2 = it1, ite = dependent.end();
    for (; it1 != ite; ++it1) {
      *it2 = *it1;
      if (*it2 != &cd)
        ++it2;
      else
        --s;
    }
    dependent.resize(s);
  }
  
  void context_dependencies::sup_dependency_
  (const context_dependencies &cd) const {
    getfem::local_guard lock = locks_.get_lock();
    size_type s = dependencies.size();
    iterator_list it1=dependencies.begin(), it2=it1, ite=dependencies.end();
    for (; it1 != ite; ++it1) {
      *it2 = *it1;
      if (*it2 != &cd)
        ++it2;
      else
        --s;
    }
    dependencies.resize(s);
  }

  void context_dependencies::invalid_context() const {
    if (state != CONTEXT_INVALID) {
      for (auto &it : dependent)
        it->invalid_context();
      getfem::local_guard lock = locks_.get_lock();
      state = CONTEXT_INVALID;
    }
  }

  void context_dependencies::add_dependency(const context_dependencies &cd) {
    cd.context_check(); cd.touched = false;
    {
      getfem::local_guard lock = locks_.get_lock();
      for (auto &it : dependencies)
        if (it == &cd) return;
      dependencies.push_back(&cd);
    }
    getfem::local_guard lock = cd.locks_.get_lock();
    cd.dependent.push_back(this);
  }
  
  bool context_dependencies::go_check() const {
    if (state == CONTEXT_CHANGED) {
      for (auto &it : dependencies) {
        it->context_check(); 
        it->touched = false;
      }
      getfem::local_guard lock = locks_.get_lock();
      state = CONTEXT_NORMAL;
      update_from_context();
      return true;
    }
    GMM_ASSERT1(state != CONTEXT_INVALID, "Invalid context");
    return false;
  }
  
  void context_dependencies::touch() const {
    if (!touched) {
      for (auto &it : dependent)
        it->change_context();
      touched = true;
    }
  }

  void context_dependencies::clear_dependencies() {
    for (auto &it : dependencies)
      it->sup_dependent_(*this);
    dependencies.clear();
  }
 
  context_dependencies::~context_dependencies() {
    invalid_context();
    for (auto &it : dependencies) it->sup_dependent_(*this);
    for (auto &it : dependent)    it->sup_dependency_(*this);
  }
  
}
