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

  class context_dependencies {
  public :
    typedef long ident_type;
    static ident_type new_ident(void);

  protected :
    typedef std::pair<context_dependencies *, long> dependency;
    ident_type ident_;
    std::list<dependency> dependencies;

  public :
    ident_type ident(void) const { return ident_; }

    void add_dependency(const context_dependencies &cd)
    { dependencies.push(dependency(&cd, cd.ident())); }
    
    bool context_changed(void) const;

  };

 

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONTEXT_H__  */
