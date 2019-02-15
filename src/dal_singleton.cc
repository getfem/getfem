/*===========================================================================

 Copyright (C) 2004-2017 Julien Pommier

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

===========================================================================*/

#include "getfem/dal_singleton.h"
#include <algorithm>
#include "gmm/gmm.h"
#include "getfem/getfem_omp.h"


namespace dal {

  singletons_manager::singletons_manager()
    : lst{}, nb_partitions{lst.num_threads()}
  {}

  singletons_manager& singletons_manager::manager(){
    static singletons_manager x;
    return x;
  }

  void singletons_manager::register_new_singleton(singleton_instance_base *p){
    register_new_singleton(p, manager().lst.this_thread());
  }

  void singletons_manager::register_new_singleton(singleton_instance_base *p, size_t ithread){
    if (p) manager().lst(ithread).push_back(p);
  }

  void singletons_manager::on_partitions_change(){
    auto new_nb_partitions = manager().lst.num_threads();
    auto &nb_partitions = manager().nb_partitions;
    if (new_nb_partitions > nb_partitions){ //not allowing reducing nb. of partitions,
      manager().lst.on_thread_update();     //as it will invalidate global storage
      nb_partitions = new_nb_partitions;
    }
  }

  static int level_compare(singleton_instance_base *a, singleton_instance_base *b) {
    return a->level() < b->level();
  }

  singletons_manager::~singletons_manager() {
    for(size_type i = 0; i != nb_partitions; ++i){
      std::sort(lst(i).begin(), lst(i).end(), level_compare);
      for (auto &&p : lst(i)) delete p;
    }
  }

}/* end of namespace dal                                                             */