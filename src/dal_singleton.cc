/*===========================================================================
 
 Copyright (C) 2004-2015 Julien Pommier
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

 
  singletons_manager::singletons_manager() : lst() {}
  
  singletons_manager& singletons_manager::m = manager();

  singletons_manager& singletons_manager::manager()
  {
    static singletons_manager x;
    return x;
  }


	void singletons_manager::register_new_singleton(singleton_instance_base *p) 
  {  	
    manager().lst.thrd_cast().push_back(p);
	}

	void singletons_manager::register_new_singleton(singleton_instance_base *p, size_t ithread) 
  {  	
    manager().lst(ithread).push_back(p);
	}


	static int level_compare(singleton_instance_base *a,
		singleton_instance_base *b) 
	{
		return a->level() < b->level();
	}

  singletons_manager::~singletons_manager() { 
    GMM_ASSERT1(!getfem::me_is_multithreaded_now(), 
                "singletons_manager destructor should" 
                "not be running in parallel !!");
    //arrange distruction per thread
    for(size_t i=0;i<getfem::num_threads();i++)
    {
      /* sort singletons in increasing levels,
      lowest levels will be destroyed first */
      std::sort(lst(i).begin(),lst(i).end(), level_compare);
      std::vector<singleton_instance_base *>::const_iterator it  = lst(i).begin();
      std::vector<singleton_instance_base *>::const_iterator ite = lst(i).end();			
      for ( ; it != ite; ++it) { delete *it; }
    }
  }
}
