#include "getfem/getfem_omp.h"
#include "getfem/getfem_level_set_contact.h"


namespace getfem{ 
#ifdef _OPENMP
	omp_lock_t get_lock()
	{
		static omp_lock_t t;
		omp_init_lock(&t);
		return t;
	}
	omp_lock_t omp_guard::single_lock=get_lock();

	/** Construct guard object and acquire our lock */
	omp_guard::omp_guard (omp_lock_t &lock) : lock_ (&lock)
		, owner_ (false)
	{
		acquire ();
	}

	/** Explicitly set our lock */
	void omp_guard::acquire ()
	{
		if (me_is_multithreaded_now()){
			omp_set_lock (lock_);
		}
		owner_ = true;
	}

	/** Explicitly unset our lock.
	* Only unset it, though, if we are still the owner.
	*/
	void omp_guard::release ()
	{
		if (owner_ && me_is_multithreaded_now()) {
			owner_ = false;
			omp_unset_lock (lock_);
		}
	}

	/** Destruct guard object, release the lock */
	omp_guard::~omp_guard ()
	{
		release ();
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
	

	region_partition::region_partition(mesh* pm,size_type id) :
		pparent_mesh(pm),original_region(0),
		partitions(num_threads())
	{
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

		omp_guard local_lock;
        size_type Nelems = original_region->size();
		size_type psize = std::ceil(static_cast<scalar_type >(Nelems)/
            static_cast<scalar_type >(num_threads()));
		mr_visitor mr(*original_region);
        size_type dummy_=0;
		for(size_type thread = 0; thread<num_threads();thread++)
		{
			partitions[thread] = level_set_contact:: 
				free_region_num(*(original_region->get_parent_mesh()));
			mesh_region& partition = pparent_mesh->region(partitions[thread]);
			for(size_type i=thread*psize;i<(thread+1)*psize && !mr.finished();i++,++mr)
			{
                if(mr.is_face()) partition.add(mr.cv(),mr.f());
                else partition.add(mr.cv());
                dummy_=partition.size();
			}
		}
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

