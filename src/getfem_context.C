/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (Getfem)             */
/* File    :  getfem_context.C : Deals with version number and             */
/*            interdependencies of objects.                                */
/*     									   */
/* Date : June 17, 2004.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */



#include <getfem_context.h>

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
    size_type s = dependent.size();
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
  
  void context_dependencies::context_check(void) const {
    if (state == CONTEXT_CHANGED) {
      state = CONTEXT_NORMAL;
      iterator_list it = dependencies.begin(), ite = dependencies.end();
      for (; it != ite; ++it) (*it)->context_check();
      update_from_context();
    }
    else if (state == CONTEXT_INVALID)
      DAL_THROW(failure_error, "Invalid context");
  }
  
  void context_dependencies::touch(void) const {
    iterator_list it = dependent.begin(), ite = dependent.end();
    for (; it != ite; ++it) (*it)->change_context();
  }
  
  context_dependencies::~context_dependencies() {
    invalid_context();
    iterator_list it = dependencies.begin(), ite = dependencies.end();
    for (; it != ite; ++it) (*it)->sup_dependent_(*this);
    it = dependent.begin(), ite = dependent.end();
    for (; it != ite; ++it) (*it)->sup_dependency_(*this);
  }
  
}
