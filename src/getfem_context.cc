// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2004-2007 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================



#include "getfem/getfem_context.h"

namespace getfem {

  static long current_id(0);
  
  long context_dependencies::new_ident(void)
  { return ++current_id; }
  
  void context_dependencies::sup_dependent_
  (const context_dependencies &cd) const {
    size_type s = dependent.size();
    iterator_list it1 = dependent.begin(), it2 = it1, ite = dependent.end();
    for (; it1 != ite; ++it1)
      { *it2 = *it1; if (*it2 != &cd) ++it2; else --s; }
    dependent.resize(s);
  }
  
  void context_dependencies::sup_dependency_
  (const context_dependencies &cd) const {
    size_type s = dependencies.size();
    iterator_list it1=dependencies.begin(), it2=it1, ite=dependencies.end();
    for (; it1 != ite; ++it1)
      { *it2 = *it1; if (*it2 != &cd) ++it2; else --s; }
    dependencies.resize(s);
  }

  void context_dependencies::invalid_context(void) const {
    if (state != CONTEXT_INVALID) {
      state = CONTEXT_INVALID;
      iterator_list it = dependent.begin(), ite = dependent.end();
      for (; it != ite; ++it) (*it)->invalid_context();
    }
  }

  void context_dependencies::add_dependency(const context_dependencies &cd) {
    // cout << "adding dep " << &cd << " à " << this << endl;
    cd.context_check(); cd.touched = false;
    iterator_list it = dependencies.begin(), ite = dependencies.end();
    for (; it != ite; ++it) if ((*it) == &cd) return;
    dependencies.push_back(&cd);
    cd.dependent.push_back(this);
  }
  
  bool context_dependencies::context_check(void) const {
    if (state == CONTEXT_CHANGED) {
      state = CONTEXT_NORMAL;
      iterator_list it = dependencies.begin(), ite = dependencies.end();
      for (; it != ite; ++it)
	{ (*it)->context_check(); (*it)->touched = false; }
      update_from_context();
      return true;
    }
    else if (state == CONTEXT_INVALID)
      GMM_ASSERT1(false, "Invalid context");
    return false;
  }
  
  void context_dependencies::touch(void) const {
    if (!touched) {
      touched = true;
      iterator_list it = dependent.begin(), ite = dependent.end();
      for (; it != ite; ++it)
	(*it)->change_context();
    }
  }
 
  context_dependencies::~context_dependencies() {
    //cout << "destruction de " << this <<  " " << typeid(*this).name() << endl;
    //cout << "state = " << state << endl;
    //cout << "nb dep = " << dependent.size() << endl;
    //cout << "nb depies = " << dependencies.size() << endl;
    invalid_context();
    iterator_list it = dependencies.begin(), ite = dependencies.end();
    for (; it != ite; ++it) (*it)->sup_dependent_(*this);
    it = dependent.begin(), ite = dependent.end();
    for (; it != ite; ++it)
      (*it)->sup_dependency_(*this);
  }
  
}
