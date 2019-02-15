/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2018 Andriy Andreykiv

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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file getfem_accumulated_distro.h
@author  Andriy Andreykiv <andriy.andreykiv@gmail.com>
@date November 21, 2018.
@brief Distribution of assembly results (matrices/vectors) for parallel
       assembly.

This is the kernel of getfem.
*/
#pragma once

#include <gmm/gmm_def.h>

#include "getfem_omp.h"

namespace getfem
{

namespace detail {

  //T is a single vector
  template<class T>
  void add_spec(const T &a, T &b,
                gmm::abstract_null_type, gmm::abstract_vector){
    gmm::add(a, b);
  }

  //T is a single matrix
  template<class T>
  void add_spec(const T &a, T &b,
                gmm::abstract_null_type, gmm::abstract_matrix){
    gmm::add(a, b);
  }

  //T is a vector of either matrices or vectors
  template<class T>
  void add_list(const T &a, T &b){
    GMM_ASSERT2(a.size() == b.size(), "size mismatch");
    auto ita = begin(a);
    auto itb = begin(b);
    auto ita_end = end(a);
    for (;ita != ita_end; ++ita, ++itb) gmm::add(*ita, *itb);
  }

  //T is a vector of vectors
  template<class T>
  void add_spec(const T &a, T &b,
                gmm::abstract_vector, gmm::abstract_vector){
    add_list(a, b);
  }

  //T is a vector of matrices
  template<class T>
  void add_spec(const T &a, T &b,
                gmm::abstract_matrix, gmm::abstract_vector){
    add_list(a, b);
  }

  //T is a single vector
  template<class T>
  void equal_resize_spec(T &a, const T &b,
                         gmm::abstract_null_type, gmm::abstract_vector){
    gmm::resize(a, gmm::vect_size(b));
  }

  //T is a single matrix
  template<class T>
  void equal_resize_spec(T &a, const T &b,
                         gmm::abstract_null_type, gmm::abstract_matrix){
    gmm::resize(a, gmm::mat_nrows(b), gmm::mat_ncols(b));
  }

  //T is a vector of either matrices or vectors
  template<class T>
  void equal_resize_list(T &a, const T &b){
    GMM_ASSERT2(a.empty(), "the first list should be still empty");
    a.resize(b.size());
    auto ita = begin(a);
    auto itb = begin(b);
    auto ita_end = end(a);
    using Component = typename T::value_type;
    using AlgoC = typename gmm::linalg_traits<Component>::linalg_type;
    for (;ita != ita_end; ++ita, ++itb){
      equal_resize_spec(*ita, *itb, gmm::abstract_null_type{}, AlgoC{});
    }
  }

  //T is a vector of vectors
  template<class T>
  void equal_resize_spec(T &a, const T &b,
                         gmm::abstract_vector, gmm::abstract_vector){
    equal_resize_list(a, b);
  }

  //T is a vector of matrices
  template<class T>
  void equal_resize_spec(T &a, const T &b,
                         gmm::abstract_matrix, gmm::abstract_vector){
    equal_resize_list(a, b);
  }

} // namespace detail

  /**Generic addition for gmm types as well as vectors of gmm types*/
  template<class T>
  void gen_add(const T &a, T &b){
    using AlgoT = typename gmm::linalg_traits<T>::linalg_type;
    using Component = typename T::value_type;
    using AlgoC = typename gmm::linalg_traits<Component>::linalg_type;
    detail::add_spec(a, b, AlgoC{}, AlgoT{});
  }

  /**Resize 'a' to the same size as 'b'.
      a and b can be matrices, vectors, or lists of matrices or vectors
  */
  template<class T>
  void equal_resize(T &a, const T &b){
    using AlgoT = typename gmm::linalg_traits<T>::linalg_type;
    using Component = typename T::value_type;
    using AlgoC = typename gmm::linalg_traits<Component>::linalg_type;
    detail::equal_resize_spec(a, b, AlgoC{}, AlgoT{});
  }

  /**Takes a matrix or vector, or vector of matrices or vectors and
  creates an empty copy on each thread. When the
  thread computations are done (in the destructor), accumulates
  the assembled copies into the original*/
  template <class T>
  class accumulated_distro
  {
    T& original;
    omp_distribute<T> distributed;

  public:

    explicit accumulated_distro(T& l)
      : original{l}{
      if (distributed.num_threads() == 1) return;
      //intentionally skipping thread 0, as accumulated_distro will
      //use original for it
      for(size_type t = 1; t != distributed.num_threads(); ++t){
        equal_resize(distributed(t), original);
      }
    }

    operator T&(){
      if (distributed.num_threads() == 1 ||
          distributed.this_thread() == 0) return original;
      else return distributed;
    }

    T& operator = (const T &x){
      return distributed = x;
    }

    ~accumulated_distro(){
      if (distributed.num_threads() == 1) return;

      if (me_is_multithreaded_now()) {
        // GMM_ASSERT1 not convenient here
        cerr <<  "Accumulation distribution should not run in parallel";
        exit(1);
      }

      using namespace std;
      auto to_add = vector<T*>{};
      to_add.push_back(&original);
      for (size_type t = 1; t != distributed.num_threads(); ++t){
        to_add.push_back(&distributed(t));
      }

      //Accumulation in parallel.
      //Adding, for instance, elements 1 to 0, 2 to 3, 5 to 4 and 7 to 6
      //on separate 4 threads in case of parallelization of the assembly
      //on 8 threads.
      while (to_add.size() > 1){
        GETFEM_OMP_PARALLEL(
          auto i = distributed.this_thread() * 2;
          if (i + 1 < to_add.size()){
            auto &target = *to_add[i];
            auto &source = *to_add[i + 1];
            gen_add(source, target);
          }
        )
        //erase every second item , as it was already added
        for (auto it = begin(to_add);
             next(it) < end(to_add);
             it = to_add.erase(next(it)));
      }
    }
  };

} // namespace getfem
