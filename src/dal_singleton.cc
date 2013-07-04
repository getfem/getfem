/*===========================================================================

Copyright (C) 2004-2012 Julien Pommier

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
namespace dal {
	shared_ptr<singletons_manager> singletons_manager::m(0);

	void singletons_manager::register_new_singleton(singleton_instance_base *p) {  	
		manager_pointer()->lst.thrd_cast().push_back(p);
	}

	void singletons_manager::register_new_singleton(singleton_instance_base *p, int ithread) {  	
		manager_pointer()->lst(ithread).push_back(p);
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
		for(size_t i=0;i<getfem::num_threads();i++){

			/* sort singletons in increasing levels,
			lowest levels will be destroyed first */
			std::sort(m->lst(i).begin(),m->lst(i).end(), level_compare);
			std::vector<singleton_instance_base *>::const_iterator 
				it = m->lst(i).begin(), 
                ite = m->lst(i).end();
			for ( ; it != ite; ++it) { delete *it; }
		}
	}
}
