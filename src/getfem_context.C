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

  typedef context_dependencies::ident_type ident_type;

  static ident_type current_id(0);
  
  ident_type context_dependencies::new_ident(void)
  { return ++current_id; }

  ident_type context_dependencies::current_ident(void)
  { return current_id; }

  bool context_dependencies::context_changed(void) const {
    std::list<dependency>::iterator it = dependencies.begin(),
      ite = dependencies.end();
    if (c_ident == current_ident()) return false;
    bool b = false;
    for (; it != ite; ++it)
      if (it->first->context_changed()
	  || it->first->ident() != it->second) {
	b = true; it->second = it->first->ident();
      }
    if (b) touch();
    c_ident = current_ident();
    return b;
  }
  
}
