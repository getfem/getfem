/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2003-2025 Yves Renard, Julien Pommier
               2025-2025 Konstantinos Poulios

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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file gmm_MUMPS_interface.h
   @author Yves Renard <Yves.Renard@insa-lyon.fr>,
   @author Julien Pommier <Julien.Pommier@insa-toulouse.fr>
   @date December 8, 2005.
   @brief Interface with MUMPS (LU direct solver for sparse matrices).
*/
#if defined(GMM_USES_MUMPS)

#ifndef GMM_MUMPS_INTERFACE_H
#define GMM_MUMPS_INTERFACE_H

#include "gmm_kernel.h"


extern "C" {

#include <smumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2
#include <dmumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2
#include <cmumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2
#include <zmumps_c.h>
#undef F_INT
#undef F_DOUBLE
#undef F_DOUBLE2

}

namespace gmm {

  template <typename T>
  struct ij_sparse_matrix {
    typedef typename number_traits<T>::magnitude_type R;
    std::vector<int> irn;
    std::vector<int> jcn;
    std::vector<T> a;

    static const std::vector<size_type> no_sel;

    // input 0-based, output 1-based
    void build_indices_vector(const std::vector<size_type> &inp,
                              std::vector<int> &out) {
      if (inp.empty())
        for (size_type i = 0; i < out.size(); ++i) out[i] = int(i+1);
      else
        for (size_type i = 0; i < inp.size(); ++i) out[inp[i]] = int(i+1);
    }

    // build an i,j,a matrix from matrix A, optionally performing row and
    // column selection/permutation on A, and optionally keeping only the
    // lower triangle part of the resulting matrix (after permutations)
    template <typename L>
    void build_from(const L& A, row_major,
                    bool lower_triangular=false,
                    const std::vector<size_type> &rows=no_sel,
                    const std::vector<size_type> &cols=no_sel) {
      std::vector<int> row_ind(mat_nrows(A),0), col_ind(mat_nrows(A),0);
      build_indices_vector(rows, row_ind);
      build_indices_vector(cols, col_ind);
      for (size_type i = 0; i < mat_nrows(A); ++i) {
        const int ir = row_ind[i];
        if (ir > 0) {
          auto row = mat_const_row(A, i);
          auto it = vect_const_begin(row), ite = vect_const_end(row);
          for (; it != ite; ++it) {
            const int jc = col_ind[it.index()];
            if (jc > 0 && (*it != T(0))
                       && (!lower_triangular || ir >= jc)) {
              irn.push_back(ir);
              jcn.push_back(jc);
              a.push_back(*it);
            }
          }
        }
      }
    }

    template <typename L>
    void build_from(const L& A, col_major,
                    bool lower_triangular=false,
                    const std::vector<size_type> &rows=no_sel,
                    const std::vector<size_type> &cols=no_sel) {
      std::vector<int> row_ind(mat_nrows(A),0), col_ind(mat_nrows(A),0);
      build_indices_vector(rows, row_ind);
      build_indices_vector(cols, col_ind);
      for (size_type j = 0; j < mat_ncols(A); ++j) {
        const int jc = col_ind[j];
        if (jc >0) {
          auto col = mat_const_col(A, j);
          auto it = vect_const_begin(col), ite = vect_const_end(col);
          for (; it != ite; ++it) {
            const int ir = row_ind[it.index()];
            if (ir > 0 && (*it != T(0))
                       && (!lower_triangular || ir >= jc)) {
              irn.push_back(ir);
              jcn.push_back(jc);
              a.push_back(*it);
            }
          }
        }
      }
    }

    template <typename L>
    ij_sparse_matrix(const L& A, bool lower_triangular=false,
                     const std::vector<size_type> &rows=no_sel,
                     const std::vector<size_type> &cols=no_sel)
    { // do not reserve nnz(A) entires in case only a sub-matrix is used
      //size_type nz = nnz(A);
      //irn.reserve(nz); jcn.reserve(nz); a.reserve(nz);
      build_from(A, typename principal_orientation_type
                             <typename linalg_traits<L>::sub_orientation>
                             ::potype(),
                 lower_triangular, rows, cols);
    }
  };

  template<typename T> const std::vector<size_type>
  ij_sparse_matrix<T>::no_sel= std::vector<size_type>();


  /* ********************************************************************* */
  /*   MUMPS solve interface                                               */
  /* ********************************************************************* */

  template <typename T> struct mumps_interf {};

  template <> struct mumps_interf<float> {
    typedef SMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef float value_type;

    static void mumps_c(MUMPS_STRUC_C &id) { smumps_c(&id); }
  };

  template <> struct mumps_interf<double> {
    typedef DMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef double value_type;
    static void mumps_c(MUMPS_STRUC_C &id) { dmumps_c(&id); }
  };

  template <> struct mumps_interf<std::complex<float> > {
    typedef CMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef mumps_complex value_type;
    static void mumps_c(MUMPS_STRUC_C &id) { cmumps_c(&id); }
  };

  template <> struct mumps_interf<std::complex<double> > {
    typedef ZMUMPS_STRUC_C  MUMPS_STRUC_C;
    typedef mumps_double_complex value_type;
    static void mumps_c(MUMPS_STRUC_C &id) { zmumps_c(&id); }
  };


  static inline bool
  mumps_error_check(int INFO1, int INFO2) {
    if (INFO1 < 0) {
      switch (INFO1) {
        case -2:
          GMM_ASSERT1(false, "Solve with MUMPS failed: NZ = " << INFO2
                      << " is out of range");
          break;
        case -6 : case -10 :
          GMM_WARNING1("Solve with MUMPS failed: matrix is singular");
          return false;
        case -9:
          GMM_ASSERT1(false, "Solve with MUMPS failed: error "
                      << INFO1 << ", increase ICNTL(14)");
          break;
        case -13 :
          GMM_ASSERT1(false, "Solve with MUMPS failed: not enough memory");
          break;
        default :
          GMM_ASSERT1(false, "Solve with MUMPS failed with error "
                      << INFO1);
          break;
      }
    }
    return true;
  }

  template <typename MUMPS_STRUCT>
  [[deprecated("Use updated error handling mechanism")]]
  static inline bool mumps_error_check(MUMPS_STRUCT &id) {
    return mumps_error_check(id.info[0], id.info[1]);
  }


  /**
   * The MUMPS interface is controlled by the JOB parameter.
   *
   * The following values are possible:
   * - JOB= 1: performs the analysis phase.
   * - JOB= 2: performs the factorization phase.
   * - JOB= 3: computes the solution.
   * - JOB= 4: combines the actions of JOB= 1 with those of JOB= 2.
   * - JOB= 5: combines the actions of JOB=2 and JOB= 3.
   * - JOB= 6: combines the actions of calls with JOB= 1, JOB= 2, and JOB= 3.
   * - JOB= 7: save / restore feature: saves MUMPS internal data to disk.
   * - JOB= 8: save / restore feature: restores MUMPS internal data from disk.
   * - JOB= 9: computes before the solution phase a possible distribution for
   *           the right-hand sides.
   * - JOB=-1: initializes an instance of the package.
   * - JOB=-2: terminates an instance of the package.
   * - JOB=-3: save / restore feature: removes data saved to disk.
   * - JOB=-4: after factorization or solve phases, frees all MUMPS internal
   *           data structures except the ones from analysis.
   *
   * MUMPS structure parameter and information vector sizes
   * ICNTL(60) CNTL(15) INFO(80) INFOG(80) RINFO(20) RINFOG(20)
   */

   /** MUMPS context

   This class encapsulates the context of a MUMPS computation. It allocates a
   a copy of the system matrix and vector (used both for right-hand sides and
   solutions) and stores the MUMPS internal data structure.

   The symmetry option, 0 for unsymmetric (default), 1 for symmetric positive
   definite, and 2 for general symmetric, is passed to the constructor and
   cannot be changed later.

   To solve a linear system in one step, the "analyze_factorize_and_solve"
   function has to be used. Calling just "solve" will not run the analysis
   and factorization phases.

   */
  template <typename T>
  class mumps_context {
    typedef typename mumps_interf<T>::value_type MUMPS_T;
    typedef typename number_traits<T>::magnitude_type MUMPS_R;
    typedef typename mumps_interf<T>::MUMPS_STRUC_C MUMPS_STRUC;
    MUMPS_STRUC id;
    std::unique_ptr<ij_sparse_matrix<T> > pK;
    std::vector<T> rhs_or_sol;
    int rank; // MPI rank
    int nrows_;

  public:

    inline void run_job(int job) {
      id.job = job;
      mumps_interf<T>::mumps_c(id);
    }

    mumps_context(int sym=0) : id(), rank(0), nrows_(0) {
#ifdef GMM_USES_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
      id.job = -1;
      id.par = 1;
      id.sym = sym;
      id.comm_fortran = -987654; // USE_COMM_WORLD
      run_job(-1); // JOB_INIT
    }

    ~mumps_context() { run_job(-2); } // JOB_END

    inline int nrows() const { return nrows_; }
    inline int &ICNTL(int I) { return id.icntl[I-1]; }
    inline MUMPS_R &CNTL(int I) { return id.cntl[I-1]; }
    inline const int &INFO(int I) { return id.info[I-1]; }
    inline const int &INFOG(int I) { return id.infog[I-1]; }
    inline const MUMPS_R &RINFO(int I) { return id.rinfo[I-1]; }
    inline const MUMPS_R &RINFOG(int I) { return id.rinfog[I-1]; }
    //inline int MPI_rank() const { return rank; }

    // Row/column index vectors are expected to be 0-based, the conversion to
    // 1-based indexing, required by MUMPS, is done inside ij_sparse_matrix
    template <typename MAT>
    inline void set_matrix(const MAT &K, bool distributed,
                           const std::vector<size_type> &
                             rows=ij_sparse_matrix<T>::no_sel,
                           const std::vector<size_type> &
                             cols=ij_sparse_matrix<T>::no_sel) {
      static_assert(std::is_same<typename linalg_traits<MAT>::value_type,
                                 T>::value,
                    "value_type of MAT and T must be the same");
      GMM_ASSERT2(gmm::mat_nrows(K) == gmm::mat_ncols(K), "Non-square matrix");
      nrows_ = int(gmm::mat_nrows(K));
      if (!distributed && rank != 0)
        return;
      id.n = nrows_;
      pK = std::make_unique< ij_sparse_matrix<T> >(K, id.sym > 0, rows, cols);
      if (distributed) {
        id.nz_loc = int(pK->irn.size());
        id.irn_loc = &(pK->irn[0]);
        id.jcn_loc = &(pK->jcn[0]);
        id.a_loc = (MUMPS_T*)(&(pK->a[0]));
      } else {
        id.nz = int(pK->irn.size());
        id.irn = &(pK->irn[0]);
        id.jcn = &(pK->jcn[0]);
        id.a = (MUMPS_T *)(&(pK->a[0]));
      }
    }

    template <typename VEC>
    inline void set_vector(const VEC &rhs) {
      static_assert(std::is_same<typename linalg_traits<VEC>::value_type,
                                 T>::value,
                    "value_type of MAT and T must be the same");
      GMM_ASSERT2(nrows() > 0,
                  "System size not defined, need to call set_matrix first.");
      const int nrhs = int(gmm::vect_size(rhs)/nrows());
      GMM_ASSERT2(size_type(nrhs*nrows()) == gmm::vect_size(rhs),
                  "Size of rhs (" << gmm::vect_size(rhs) << ") must be an "
                  "integer multiple of the matrix size (" << id.n << ")");
      rhs_or_sol.resize(gmm::vect_size(rhs));
      gmm::copy(rhs, rhs_or_sol);
      if (rank == 0) {
        id.nrhs = nrhs;
        id.rhs = (MUMPS_T*)(&rhs_or_sol[0]);
        id.lrhs = id.n;
      }
    }

    const std::vector<T> &vector() const { return rhs_or_sol; }

    inline void analyze() { run_job(1); }
    inline void factorize() { run_job(2); }
    inline void analyze_and_factorize() { run_job(4); }
    inline void solve() { run_job(3); }
    inline void factorize_and_solve() { run_job(5); }
    inline void analyze_factorize_and_solve() { run_job(6); }
    inline bool error_check() { return mumps_error_check(INFO(1), INFO(2)); }

    inline void mpi_broadcast() {
#ifdef GMM_USES_MPI
      MPI_Bcast(&(rhs_or_sol[0]), int(rhs_or_sol.size()),
                gmm::mpi_type(T()), 0, MPI_COMM_WORLD);
#endif
    }
  };


  /** MUMPS solve interface
   *  Works only with sparse or skyline matrices
   */
  template <typename MAT, typename VECTX, typename VECTB>
  bool MUMPS_solve(const MAT &A, VECTX &X, const VECTB &B,
                   bool sym = false, bool distributed = false) {

    typedef typename linalg_traits<MAT>::value_type T;

    bool ok=false;
    {
      mumps_context<T> mumps_ctx(sym ? 2 : 0); // General symmetric (2) or unsymmetric (0)
      mumps_ctx.set_matrix(A, distributed);
      mumps_ctx.set_vector(B);
      mumps_ctx.ICNTL(1) = -1;   // output stream for error messages
      mumps_ctx.ICNTL(2) = -1;   // output stream for other messages
      mumps_ctx.ICNTL(3) = -1;   // output stream for global information
      mumps_ctx.ICNTL(4) = 0;    // verbosity level
      if (distributed)
        mumps_ctx.ICNTL(5) = 0;  // assembled input matrix (default)
      mumps_ctx.ICNTL(14) += 80; // small boost to the workspace size as we have encountered
                                 // some problems that did not fit in the default settings of
                                 // mumps... ICNTL(14) = 15 or 20
      if (distributed)
        mumps_ctx.ICNTL(18) = 3; // strategy for distributed input matrix
        //mumps_ctx.ICNTL(22) = 1; // enables out-of-core support
      mumps_ctx.analyze_factorize_and_solve();
      ok = mumps_ctx.error_check();
      mumps_ctx.mpi_broadcast();
      gmm::copy(mumps_ctx.vector(), X);
    } // end scope of mumps_ctx, destructor calls mumps job=-2
    return ok;
  }


  /** MUMPS solve interface for distributed matrices
   *  Works only with sparse or skyline matrices
   */
  template <typename MAT, typename VECTX, typename VECTB>
  bool MUMPS_distributed_matrix_solve(const MAT &A, VECTX &X_,
                                      const VECTB &B, bool sym=false) {
    return MUMPS_solve(A, X_, B, sym, true);
  }


  template<typename T>
  inline T real_or_complex(std::complex<T> a) { return a.real(); }
  template<typename T>
  inline T real_or_complex(T &a) { return a; }


  /** Evaluate matrix determinant with MUMPS
   *  Works only with sparse or skyline matrices
   */
  template <typename MAT, typename T=typename linalg_traits<MAT>::value_type>
  T MUMPS_determinant(const MAT &A, int &exponent,
                      bool sym=false, bool distributed=false) {
    exponent = 0;
    typedef typename number_traits<T>::magnitude_type R;

    mumps_context<T> mumps_ctx(sym ?  2 : 0); // General symmetric (2)
    mumps_ctx.set_matrix(A, distributed);     // or unsymmetric (0)
    mumps_ctx.ICNTL(4) = 0;    // verbosity level
    if (distributed)
      mumps_ctx.ICNTL(5) = 0;  // assembled input matrix (default)
    mumps_ctx.ICNTL(14) += 80; // small boost to the workspace size
    if (distributed)
      mumps_ctx.ICNTL(18) = 3; // strategy for distributed input matrix
    mumps_ctx.ICNTL(31) = 1;   // only factorization, no solution to follow
    mumps_ctx.ICNTL(33) = 1;   // request determinant calculation

    mumps_ctx.analyze_and_factorize();

    T det = real_or_complex(std::complex<R>(mumps_ctx.RINFOG(12),
                                            mumps_ctx.RINFOG(13)));
    exponent = mumps_ctx.INFOG(34);
    return det;
  }

}

#endif // GMM_MUMPS_INTERFACE_H

#endif // GMM_USES_MUMPS
