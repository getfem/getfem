/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (Getfem)             */
/* File    :  getfem_context.h : Deals with version number and             */
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


#ifndef GETFEM_CONTEXT_H__
#define GETFEM_CONTEXT_H__

#include <getfem_config.h>
#include <list>

namespace getfem {

  // An object can be in three different states :
  //   NORMAL  : no change is necessary
  //   CHANGED : something in the context has changed and the update function
  //             of the object has to be called.
  //   INVALID : one of the dependencies desappears, the object is invalid
  //             and should no longer be used.
  //
  // add_dependency(ct) : add a dependency to the dependency list.
  // touch()            : significate to the dependent objects that something
  //                      has change in the object. This make the dependent
  //                      objects to be in the CHANGED state
  // context_check()    : check if the object has to be updated. if it is the
  //                      case make first a check to the dependency list and
  //                      call the update function of the object.
  //                      (the update function of the dependencies are called
  //                      before the update function of the current object).
  // context_valid()    : says if the object has still a valid context
  //
  // Remarks :
  // - A protection against round dependencies exists. In this case, the order
  //   of call of the update functions can be arbitrary
  // - Detection of context changes is very fast (control the state). the touch
  //   operation can cover the whole tree of dependent object.
  //   But this is the case only for the first touch operation because once
  //   a dependent object is in the CHANGED state it will not be considered
  //   by next touch operations.


  class context_dependencies {

  protected :
    enum context_state { CONTEXT_NORMAL, CONTEXT_CHANGED, CONTEXT_INVALID };
    mutable context_state state;
    mutable bool touched;
    mutable std::vector<const context_dependencies *> dependencies;
    mutable std::vector<const context_dependencies *> dependent;
    typedef std::vector<const context_dependencies *>::iterator iterator_list;

    void sup_dependent_(const context_dependencies &cd) const;
    void sup_dependency_(const context_dependencies &cd) const;
    void change_context(void) const
    { if (state == CONTEXT_NORMAL) { state = CONTEXT_CHANGED; touch(); } }
    void invalid_context(void) const;

  public :
    
    // this function has to be defined and should update the object when
    // the context is modified.
    virtual void update_from_context(void) const = 0;

    void add_dependency(const context_dependencies &cd)
    { dependencies.push_back(&cd);cd.dependent.push_back(this);touched=false; }
    void sup_dependency(const context_dependencies &cd)
    { cd.sup_dependent_(*this); sup_dependency_(cd); }
    bool context_valid(void) const { return (state != CONTEXT_INVALID); }
    void context_check(void) const;
    void touch(void) const;
    static long new_ident(void); // outdated function
    virtual ~context_dependencies();
    context_dependencies() : state(CONTEXT_NORMAL), touched(true) {}

  };






}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTEXT_H__  */
