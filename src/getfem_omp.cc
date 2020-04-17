/*===========================================================================

 Copyright (C) 2012-2020 Andriy Andreykiv.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
#include "getfem/getfem_locale.h"
#include "getfem/getfem_omp.h"

#ifdef GETFEM_HAS_OPENMP
  #include <thread>
  #include <omp.h>
#endif

using bgeot::scalar_type;

namespace getfem{

#ifdef GETFEM_HAS_OPENMP

  std::recursive_mutex omp_guard::mutex;

  omp_guard::omp_guard()
    : plock{me_is_multithreaded_now() ?
       std::make_unique<std::lock_guard<std::recursive_mutex>>(mutex)
      : nullptr}
  {}

  local_guard::local_guard(std::recursive_mutex& m) :
    mutex{m},
    plock{me_is_multithreaded_now() ?
      std::make_shared<std::lock_guard<std::recursive_mutex>>(m)
      : nullptr}
  {}

  local_guard lock_factory::get_lock() const{
    return local_guard{mutex};
  }

  size_type global_thread_policy::this_thread() {
    return partition_master::get().get_current_partition();
  }

  size_type global_thread_policy::num_threads(){
    return partition_master::get().get_nb_partitions();
  }

  size_type true_thread_policy::this_thread() {
    return omp_get_thread_num();
  }

  size_type true_thread_policy::num_threads(){
    return omp_get_max_threads();
  }

  void set_num_threads(int n){
    omp_set_num_threads(n);
    partition_master::get().check_threads();
  }

  bool me_is_multithreaded_now(){
    // serial region
    if(omp_get_num_threads() == 1 && omp_get_level() == 0) return false;
    // parallel region with one thread
    if(omp_get_num_threads() == 1 && omp_get_level() == 1) return true;
    // parallel region with more than one thread
    if(omp_in_parallel() == 1) return true;

    return false;
  }

  bool not_multithreaded(){
    return omp_get_max_threads() == 1;
  }

  size_type max_concurrency() {
    return std::thread::hardware_concurrency();
  }

#else

  size_type global_thread_policy::this_thread() {return 0;}

  size_type global_thread_policy::num_threads(){return 1;}

  size_type true_thread_policy::this_thread() {return 0;}

  size_type true_thread_policy::num_threads(){return 1;}

  bool me_is_multithreaded_now(){return false;}

  void set_num_threads(int /*n*/){}

  bool not_multithreaded(){return true;}

  size_type max_concurrency() {return 1;}

#endif

  /** Allows to re-throw exceptions, generated in OpemMP parallel section.
      Collects exceptions from all threads and on destruction re-throws
      the first one, so that
      it can be again caught in the master thread. */
  class thread_exception {
    std::vector<std::exception_ptr> exceptions;

    void captureException(){
      exceptions[true_thread_policy::this_thread()] = std::current_exception();
    }

  public:
    thread_exception()
      : exceptions(true_thread_policy::num_threads(), nullptr)
    {}

    template <typename function, typename... parameters>
    void run(function f, parameters... params){
        try {f(params...);} catch (...) {captureException();}
    }

    std::vector<std::exception_ptr> caughtExceptions() const{
      std::vector<std::exception_ptr> non_empty_exceptions;
      for (auto &&pException : exceptions){
        if (pException != nullptr) non_empty_exceptions.push_back(pException);
      }
      return non_empty_exceptions;
    }

    void rethrow() {
      for (auto &&pException : exceptions){
        if (pException != nullptr) std::rethrow_exception(pException);
      }
    }
  };

  partition_iterator::partition_iterator(
    partition_master &m, std::set<size_type>::const_iterator it_from_set)
    : master{m}, it{it_from_set}
  {}

  partition_iterator partition_iterator::operator++(){
    ++it;
    if (*this != master.end()) master.set_current_partition(*it);
    return *this;
  }

  bool partition_iterator::operator==(const partition_iterator &it1) const {
    return it == it1.it;
  }

  bool partition_iterator::operator!=(const partition_iterator &it1) const {
    return !(*this == it1);
  }

  size_type partition_iterator::operator*() const{
    return *it;
  }

  partition_master partition_master::instance;

  partition_master& partition_master::get(){
    return instance;
  }

  void partition_master::check_threads(){
    GLOBAL_OMP_GUARD
    auto must_update = false;
    if (nb_user_threads != true_thread_policy::num_threads()){
      nb_user_threads = true_thread_policy::num_threads();
      must_update = true;
    }
    if (nb_partitions < nb_user_threads && !partitions_set_by_user){
      nb_partitions = nb_user_threads;
      must_update = true;
    }
    if (must_update){
      update_partitions();
      dal::singletons_manager::on_partitions_change();
    }
  }

  void partition_master::set_nb_partitions(size_type n){
    GMM_ASSERT1 (!partitions_set_by_user,
                 "Number of partitions can be set only once.");
    if (n > nb_partitions){
      nb_partitions = n;
      nb_user_threads = true_thread_policy::num_threads();
      update_partitions();
      dal::singletons_manager::on_partitions_change();
    }
    else if (n < nb_partitions){
      GMM_WARNING1("Not reducing number of partitions from "
                   << nb_partitions <<" to " << n <<
                   " as it might invalidate global storage.");
    }
    partitions_set_by_user = true;
  }

  partition_iterator partition_master::begin(){
    GMM_ASSERT1(nb_user_threads == true_thread_policy::num_threads(),
                "The number of omp threads was changed outside partition_master."
                "Please use getfem::set_num_threads for this.");
    current_partition = *(std::begin(partitions.thrd_cast()));
    return partition_iterator{*this, std::begin(partitions.thrd_cast())};
  }

  partition_iterator partition_master::end(){
    return partition_iterator{*this, std::end(partitions.thrd_cast())};
  }

  void partition_master::set_behaviour(thread_behaviour b){
    if (b != behaviour){
      GMM_ASSERT1(!me_is_multithreaded_now(),
                  "Cannot change thread policy in parallel section.");
      behaviour = b;
      check_threads();
    }
  }

  partition_master::partition_master()
    : nb_user_threads{1}, nb_partitions{1} {
        partitions_updated = false;
        set_num_threads(1);
        update_partitions();
  }

  size_type partition_master::get_current_partition() const {
    GMM_ASSERT2(behaviour == thread_behaviour::partition_threads ?
                true_thread_policy::this_thread() < nb_partitions : true,
                "Requesting current partition for thread " <<
                true_thread_policy::this_thread() <<
                " while number of partitions is " << nb_partitions
                << ".");
    return behaviour == thread_behaviour::partition_threads ?
           current_partition : true_thread_policy::this_thread();
  }

  size_type partition_master::get_nb_partitions() const {
    return behaviour == thread_behaviour::partition_threads ?
           nb_partitions : true_thread_policy::num_threads();
  }

  void partition_master::set_current_partition(size_type p){
    if (behaviour == thread_behaviour::partition_threads){
      GMM_ASSERT2(partitions.thrd_cast().count(p) != 0, "Internal error: "
                  << p << " is not a valid partitions for thread "
                  << true_thread_policy::this_thread()
                  << ".");
      current_partition = p;
    }
  }

  void partition_master::rewind_partitions(){
    if (me_is_multithreaded_now()){
      current_partition = *(std::begin(partitions.thrd_cast()));
    }
    else{
      for (size_type t = 0; t != partitions.num_threads(); ++t){
        current_partition(t) = *(std::begin(partitions(t)));
      }
    }
  }

  void partition_master::update_partitions(){
    partitions_updated = false;

    GLOBAL_OMP_GUARD

    if (partitions_updated) return;

    partitions = decltype(partitions){};
    current_partition = decltype(current_partition){};

    auto n_threads = true_thread_policy::num_threads();
    if(n_threads > nb_partitions){
      GMM_WARNING0("Using " << n_threads <<
                   " threads which is above the maximum number of partitions :" <<
                   nb_partitions
                   << ".");
    }
    if (behaviour == thread_behaviour::partition_threads){
      for (size_type t = 0; t != n_threads; ++t){
        auto partition_size = static_cast<size_type>
                                (std::ceil(static_cast<scalar_type>(nb_partitions) /
                                           static_cast<scalar_type >(n_threads)));
        auto partition_begin = partition_size * t;
        if (partition_begin >= nb_partitions) break;
        auto partition_end = std::min(partition_size * (t + 1), nb_partitions);
        auto hint_it = std::begin(partitions(t));
        for (size_type i = partition_begin; i != partition_end; ++i){
          hint_it = partitions(t).insert(hint_it, i);
        }
        current_partition(t) = partition_begin;
      }
    }
    else{
      for (size_type t = 0; t != n_threads; ++t){
        partitions(t).insert(t);
        current_partition(t) = t;
      }
    }

    partitions_updated = true;
  }

  #if defined _WIN32 && !defined (__GNUC__)
    #define GETFEM_ON_WIN
  #endif

  parallel_boilerplate::
  parallel_boilerplate()
  : plocale{std::make_unique<standard_locale>()},
    pexception{std::make_unique<thread_exception>()} {
    #ifdef GETFEM_ON_WIN
      _configthreadlocale(_ENABLE_PER_THREAD_LOCALE);
    #endif
  }

  void parallel_boilerplate::run_lambda(std::function<void(void)> lambda){
    pexception->run(lambda);
  }

  parallel_boilerplate::~parallel_boilerplate(){
    #ifdef GETFEM_ON_WIN
    _configthreadlocale(_DISABLE_PER_THREAD_LOCALE);
    #endif
    pexception->rethrow();
  }

  void parallel_execution(std::function<void(void)> lambda,
                          bool iterate_over_partitions){
    if (me_is_multithreaded_now()) {
      lambda();
      return;
    }
    parallel_boilerplate boilerplate;
    auto &pm = partition_master::get();
    if (pm.get_nb_partitions() < true_thread_policy::num_threads()){
      pm.set_nb_partitions(true_thread_policy::num_threads());
    }
    #pragma omp parallel default(shared)
    {
      if (iterate_over_partitions) {
        for (auto &&partitions : partition_master::get()) {
          (void)partitions;
          boilerplate.run_lambda(lambda);
        }
      }
      else {
        boilerplate.run_lambda(lambda);
      }
    }
    if (iterate_over_partitions) partition_master::get().rewind_partitions();
  }

}  /* end of namespace getfem.                                             */