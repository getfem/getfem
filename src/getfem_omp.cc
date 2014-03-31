#include "getfem/getfem_omp.h"
#include "getfem/getfem_level_set_contact.h"


namespace getfem{ 

#ifdef GETFEM_HAVE_OPENMP
  boost::mutex omp_guard::boost_mutex;
  omp_guard::omp_guard() : boost::lock_guard<boost::mutex>(boost_mutex) {}
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


  region_partition::region_partition(mesh* pm,size_type id) :
    pparent_mesh(pm),original_region(0),
    partitions(num_threads())
  {
    scalar_type time = gmm::uclock_sec();
    // in case of serial Getfem nothing to partition
    if (num_threads()==1) {partitions[0]=id; return;}

    //in case mesh is not provided, also don't do anything
    if (!pm) return;

    if (id==-1) {
      original_region.reset(new mesh_region(pm->convex_index()));
      original_region->set_parent_mesh(pm);
    } else{
      GMM_ASSERT1(pm->has_region(id),"Improper region number");
      original_region.reset(new mesh_region(pm->region(id)));
    }
    if (me_is_multithreaded_now()) 
      GMM_WARNING0("building partitions inside parallel region");

    omp_guard scoped_lock;
    size_type Nelems = original_region->size();
    size_type psize = std::ceil(static_cast<scalar_type >(Nelems)/
      static_cast<scalar_type >(num_threads()));
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

