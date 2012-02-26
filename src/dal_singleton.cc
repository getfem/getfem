/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2004-2012 Julien Pommier
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
===========================================================================*/

#include "getfem/dal_singleton.h"
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
