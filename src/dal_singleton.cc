// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Dynamic Array Library (dal)
// File    : dal_singleton.cc : basic template class for singletons
//           
// Date    : May 2004.
// Author  : Julien Pommier <Julien.Pommier@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Julien Pommier
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

#include <dal_singleton.h>
#include <algorithm>
namespace dal {
  std::auto_ptr<singletons_manager> singletons_manager::m;

  void singletons_manager::register_new_singleton(singleton_instance_base *p) {
    if (!m.get()) m.reset(new singletons_manager());
    m->lst.push_back(p);
  }

  static int level_compare(singleton_instance_base *a, singleton_instance_base *b) {
    return a->level() < b->level();
  }

  singletons_manager::~singletons_manager() { 
    /* sort singletons in increasing levels,
       lowest levels will be destroyed first */
    std::sort(m->lst.begin(),m->lst.end(), level_compare);
    std::vector<singleton_instance_base *>::const_iterator 
      it = m->lst.begin(), ite = m->lst.end();
    for ( ; it != ite; ++it) { delete *it; }
  }
}
