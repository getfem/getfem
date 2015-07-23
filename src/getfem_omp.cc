/*===========================================================================

 Copyright (C) 2012-2015 Andriy Andreykiv.

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

#include "getfem/getfem_omp.h"
#include "getfem/getfem_omp.h"
#include "getfem/getfem_level_set_contact.h"

namespace getfem{ 

#ifdef GETFEM_HAVE_OPENMP

  boost::recursive_mutex omp_guard::boost_mutex;

  omp_guard::omp_guard() 
    : boost::lock_guard<boost::recursive_mutex>(boost_mutex) 
  {}

  local_guard::local_guard(boost::recursive_mutex& m) : 
    mutex_(m), 
    plock_(new boost::lock_guard<boost::recursive_mutex>(m))
  { }

  local_guard::local_guard(const local_guard& guard) 
    : mutex_(guard.mutex_), plock_(guard.plock_)
  { }

  lock_factory::lock_factory() : mutex_() {}
  local_guard lock_factory::get_lock() const
  {
    return local_guard(mutex_);
  }
#endif





  omp_distribute<bool> open_mp_is_running_properly::answer = false;
  open_mp_is_running_properly::open_mp_is_running_properly()
  {answer.all_threads()=true;}
  open_mp_is_running_properly::~open_mp_is_running_properly()
  {answer.all_threads()=false;}
  bool open_mp_is_running_properly::is_it(){return answer;}

  region_partition::region_partition(const region_partition& rp) : 
    pparent_mesh(rp.pparent_mesh),
    original_region(rp.original_region),
    partitions(rp.partitions)  {   }

  void region_partition::operator=(const region_partition& rp)
  {
    partitions.clear();

    if (!rp.pparent_mesh) return;
    pparent_mesh->copy_from(*rp.pparent_mesh);
    original_region = rp.original_region;
    partitions.resize(rp.partitions.size());
    gmm::copy(rp.partitions,partitions);
  }


  region_partition::region_partition(mesh* pm, size_type id) :
    pparent_mesh(pm),original_region(0),
    partitions(num_threads())
  {
    scalar_type time = gmm::uclock_sec();
    // in case of serial Getfem nothing to partition
    if (num_threads()==1) {partitions[0]=id; return;}

    //in case mesh is not provided, also don't do anything
    if (!pm) return;

    if (id == size_type(-1)) {
      original_region.reset(new mesh_region(pm->convex_index()));
      original_region->set_parent_mesh(pm);
    } else{
      GMM_ASSERT1(pm->has_region(id),"Improper region number");
      original_region.reset(new mesh_region(pm->region(id)));
    }
    if (me_is_multithreaded_now()) 
      GMM_WARNING0("building partitions inside parallel region");

    omp_guard scoped_lock;
    GMM_NOPERATION(scoped_lock);
    size_type Nelems = original_region->size();
    size_type psize = static_cast<size_type>
      (std::ceil(static_cast<scalar_type >(Nelems)/
       static_cast<scalar_type >(num_threads())));
    mr_visitor mr(*original_region);
    for(size_type thread = 0; thread<num_threads();thread++)
    {
      partitions[thread] = 
        getfem::mesh_region::free_region_id(*(original_region->get_parent_mesh()));
      mesh_region partition;
      for(size_type i=thread*psize;i<(thread+1)*psize && !mr.finished();i++,++mr)
      {
        if (mr.is_face()) partition.add(mr.cv(),mr.f());
        else partition.add(mr.cv());
      }
      pparent_mesh->region(partitions[thread]) = partition;
    }
    GMM_TRACE2("Partitioning time: "<<gmm::uclock_sec()-time<<" s.");
  }

  size_type region_partition::
    thread_local_partition() const {
      if (pparent_mesh==0 && num_threads() >1 ){
        GMM_WARNING1("partition is empty and cannot be used \
                     this means that the brick that created it should partition \
                     its domain by himself");
        return -10;
      }
      return partitions[this_thread()];
  }

  void omp_distribute<bool>::all_values_proxy::operator=(const bool& x)
  {
    for(std::vector<BOOL>::iterator it=distro.thread_values.begin();
      it!=distro.thread_values.end();it++) *it=x;

  }

}
