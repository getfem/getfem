/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2000-2015 Yves Renard

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

/*! @mainpage GetFEM++ reference documentation.

  <center><img src="http://download.gna.org/getfem/html/homepage/_static/logo_getfem_small.png"></center>
 
  @section intro Introduction
 
  This documentation has been generated from the GetFEM++ sources, this is not a tutorial.
 
  @section Terminology
 
  This is just a short summary of the terms employed in getfem.

  The @link mesh mesh@endlink is composed of @link convexes
  convexes@endlink. What we call convexes can be simple line segments,
  prisms, tetrahedrons, curved triangles, of even something which is
  not convex (in the geometrical sense). They all have an associated
  @link refconvexes reference convex@endlink: for segments, this will
  be the @f$[0,1]@f$ segment, for triangles this will be the canonical
  triangle @f$(0,0)-(0,1)-(1,0)@f$ etc...  All convexes of the mesh
  are constructed from the reference convex through a @link
  bgeot_geometric_trans.h geometric transformation@endlink. In
  simple cases (when the convexes are simplices for example), this
  transformation will be linear (hence it is easily inverted, which
  can be a great advantage). In order to define the geometric
  transformation, one defines <i>geometrical nodes</i> on the
  reference convex. The geometrical transformation maps these nodes to
  the <i>mesh nodes</i>.

  On the mesh, one defines a set a <i>basis functions</i>: the @link pfem FEM@endlink. A FEM
  should be associated with each convex. The basis functions are also
  attached to some geometrical points (which can be arbitrarily
  chosen). These points are similar to the mesh nodes, but <b>they
  don't have to be the same</b> (this only happens on very simple cases,
  such as a classical P1 fem on a triangular mesh). The set of all
  basis functions on the mesh forms the basis of a vector space, on
  which the PDE will be solved. These basis functions (and their
  associated geometrical point) are the <b>degrees of freedom
  (dof)</b>. The FEM is said to be <i>Lagrangian</i> when each of its basis
  functions is equal to one at its attached geometrical point, and is
  null at the geometrical points of others basis functions. This is an
  important property as it is very easy to <i>interpolate</i> an arbitrary
  function on the finite elements space.

  The finite elements methods involves evaluation of integrals of
  these basis functions (or product of basis functions etc...) on 
  convexes (and faces of convexes). In simple cases (polynomial basis
  functions and linear geometrical transformation), one can evaluate
  analytically these integrals. In other cases, one has to approximate
  it, using <i>quadrature formulas</i>. Hence, at each convex is attached
  an @link getfem_integration.h integration method@endlink along with the FEM.  If you have to use an
  approximate integration method, always choose carefully its
  order(i.e. highest degree of the polynomials who are exactly
  integrated with the method) : the degree of the FEM, of the
  polynomial degree of the geometrical transformation, and the nature
  of the elementary matrix have to be taken into account. If you are
  unsure about the appropriate degree, always prefer a high order
  integration method (which will slow down the assembly) to a low
  order one which will produce a useless linear-system.

  The process of construction of a global linear system from integrals
  of basis functions on each convex is the @link asm assembly@endlink.

  A mesh, with a set of FEM attached to its convexes is called a @link getfem_mesh_fem.h
  mesh_fem@endlink object in GetFEM++.

  A mesh, with a set of integration methods attached to its convexes
  is called a @link getfem_mesh_im.h mesh_im@endlink object in GetFEM++.

  A @c mesh_fem can be used to approximate scalar fields (heat, pression,
  ..), or vector fields (displacement, electric field, ..). A @c mesh_im
  will be used to perform numerical integrations on these fields.
  Although GetFEM++ supports vector FEMs, intrinsic vector FEMs are
  essentially used in mixed methods (for instance Raviart-Thomas element).
  Most of the FEM currently available are scalar. Of course,
  these scalar FEMs can be used to approximate each component of a
  vector field. This is done by setting the Qdim of the mesh_fem to
  the dimension of the vector field (i.e. Qdim=1 for scalar field,
  Qdim=2 for 2D vector field etc...).

  When solving a PDE, one often has to use more than one FEM. The most
  important one will be of course the one on which is defined the
  solution of the PDE. But most PDEs involve various coefficients, for
  example: @f[ \nabla.(\lambda(x)\nabla u) = f(x).  @f] Hence one has
  to define a FEM for the main unknown @f$u@f$, but also for the data
  @f$\lambda(x)@f$ and @f$f(x)@f$ if they are not constant. In order
  to interpolate easily these coefficients in their finite element
  space, one often choose a Lagrangian FEM.

  The convexes, mesh nodes, and dof are all numbered. We sometimes
  refer to the number associated to a convex as its convex id or
  convex number. Mesh node numbers are also called point id or point
  number.  Faces of convexes do not have a global numbering, but only
  a local number in each convex.

  While the @c dof are always numbered consecutively, <b>this is not
  always the case for point ids and convex ids</b>, especially if you
  have removed points or convexes from the mesh. To ensure that they
  form a continuous sequence (starting from 0), you have to use the
  getfem::getfem_mesh::optimize_structure member function.

  Most of the example programs of getfem now uses @link bricks model bricks@endlink. 
  Connecting these basic blocks one to each other, solvers for many PDE problems can be built.

  @section Examples
  
  - @link laplacian.cc tests/laplacian.cc@endlink: solve the laplace equation.
  - @link elastostatic.cc tests/elastostatic.cc@endlink: solve a static linear elasticity problem.
  - @link helmholtz.cc tests/helmholtz.cc@endlink: scalar wave equation, with robin boundary conditions.
  - @link stokes.cc tests/stokes.cc@endlink: the Stokes equation (incompressible viscuous fluid).
  - @link nonlinear_elastostatic.cc tests/nonlinear_elastostatic.cc@endlink: large strain elastostatic problem (torsion of a bar).
  - @link icare.cc contrib/icare/icare.cc@endlink: Navier-Stokes equation (fluid flow around an obstacle).
*/

/**@file getfem_config.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date November 19, 2000.
   @brief defines and typedefs for namespace getfem
*/


#ifndef GETFEM_CONFIG_H__
#define GETFEM_CONFIG_H__

#include "bgeot_config.h"

// Parallelisation options

// GETFEM_PARA_LEVEL is the parallelisation level of Getfem
//    0 - Sequential
//    1 - Only the resolution of linear systems are parallelized
//    2 - Assembly procedures are also parallelized
#ifndef GETFEM_PARA_LEVEL
# define GETFEM_PARA_LEVEL 0
#endif

#define GETFEM_MPI_INIT(argc, argv) {GMM_TRACE1("Running sequential Getfem");}
#define GETFEM_MPI_FINALIZE {}

#if defined(GETFEM_HAVE_DMUMPS_C_H)
# ifndef GMM_USES_MUMPS
#   define GMM_USES_MUMPS
# endif
#endif


#if GMM_USES_MPI > 0
# include <mpi.h>

# undef GMM_TRACE_MSG_MPI
# define GMM_TRACE_MSG_MPI					         \
  int mip_rk__; MPI_Comm_rank(MPI_COMM_WORLD, &mip_rk__);	         \
  if (mip_rk__ == 0) 

# undef GETFEM_MPI_INIT
# define GETFEM_MPI_INIT(argc, argv) {                                   \
  MPI_Init(&argc, &argv);                                                \
  GMM_TRACE1("Running parallelized Getfem level " << GETFEM_PARA_LEVEL); \
  }
# undef GETFEM_MPI_FINALIZE 
# define GETFEM_MPI_FINALIZE { MPI_Finalize(); }

// GETFEM_PARA_SOLVER is the parallelisation solver used
//    MUMPS       : use direct parallel solver MUMPS
//    SCHWARZADD  : use a Schwarz additive method
#define MUMPS_PARA_SOLVER 1
#define SCHWARZADD_PARA_SOLVER 2

# ifndef GETFEM_PARA_SOLVER
#   define GETFEM_PARA_SOLVER MUMPS_PARA_SOLVER
# endif

# if GETFEM_PARA_SOLVER == MUMPS_PARA_SOLVER
#  ifndef GMM_USES_MUMPS
#    define GMM_USES_MUMPS
#  endif
# endif

#endif


#include "bgeot_tensor.h"
#include "bgeot_poly.h"
#include "getfem_superlu.h"

/// GEneric Tool for Finite Element Methods.
namespace getfem {

  using std::endl; using std::cout; using std::cerr;
  using std::ends; using std::cin;


#if GETFEM_PARA_LEVEL > 1
  template <typename T> inline T MPI_SUM_SCALAR(T a) {
    T b; MPI_Allreduce(&a,&b,1,gmm::mpi_type(a), MPI_SUM, MPI_COMM_WORLD);
    return b;
  }

  template <typename VECT> inline void MPI_SUM_VECTOR(const VECT &VV) {
    VECT &V = const_cast<VECT &>(VV);
    typedef typename gmm::linalg_traits<VECT>::value_type T;
    std::vector<T> W(gmm::vect_size(V));
    MPI_Allreduce((void *)(&(V[0])), &(W[0]), gmm::vect_size(V),
                  gmm::mpi_type(T()), MPI_SUM, MPI_COMM_WORLD);
    gmm::copy(W, V);
  }

  template <typename VECT> inline void MPI_MAX_VECTOR(const VECT &VV) {
    VECT &V = const_cast<VECT &>(VV);
    typedef typename gmm::linalg_traits<VECT>::value_type T;
    std::vector<T> W(gmm::vect_size(V));
    MPI_Allreduce((void *)(&(V[0])), &(W[0]), gmm::vect_size(V),
                  gmm::mpi_type(T()), MPI_MAX, MPI_COMM_WORLD);
    gmm::copy(W, V);
  }

  template <typename T> void MPI_BCAST0_SCALAR(T &a) {
    MPI_Bcast((void *)(&a), 1, gmm::mpi_type(a), 0, MPI_COMM_WORLD);
  }

  template <typename VECT> inline void MPI_BCAST0_VECTOR(const VECT &VV) {
    VECT &V = const_cast<VECT &>(VV);
    typedef typename gmm::linalg_traits<VECT>::value_type T;
    MPI_Bcast((void *)(&(V[0])), gmm::vect_size(V), gmm::mpi_type(T()), 0,
              MPI_COMM_WORLD);
  }

  template <typename VECT1, typename VECT2>
  inline void MPI_SUM_VECTOR(const VECT1 &VV, const VECT2 &WW) {
    VECT1 &V = const_cast<VECT1 &>(VV);
    VECT2 &W = const_cast<VECT2 &>(WW);
    typedef typename gmm::linalg_traits<VECT1>::value_type T;
    MPI_Allreduce((void *)(&(V[0])), &(W[0]),
		  gmm::vect_size(V), gmm::mpi_type(T()),
		  MPI_SUM, MPI_COMM_WORLD);
  }

  inline bool MPI_IS_MASTER(void)
  { int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk); return !rk; }

  template <typename T> inline
  void MPI_SUM_SPARSE_MATRIX(const gmm::col_matrix< gmm::rsvector<T> > &M) {
    typedef typename gmm::col_matrix< gmm::rsvector<T> > MAT;
    MAT &MM = const_cast<MAT &>(M);
    int rk, nbp;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &nbp);
    if (nbp < 2) return;
    size_type nr = gmm::mat_nrows(MM), nc = gmm::mat_ncols(MM);

    gmm::dense_matrix<int> all_nnz(nc, nbp);
    std::vector<int> my_nnz(nc), owner(nc);
    gmm::rsvector<T> V(nr);
    
    // Step 1 : each process fill the corresponding nnz line
    for (size_type i = 0; i < nc; ++i)
      my_nnz[i] = gmm::nnz(gmm::mat_col(MM, i));

    // Step 2 : All gather : each proc has all the information
    MPI_Allgather((void *)(&(my_nnz[0])), nc, MPI_INT,
                  (void *)(&(all_nnz(0,0))), nc, MPI_INT, MPI_COMM_WORLD);   
    
    // Step 3 : Scan each column and message send/receive
    std::vector<int> contributors(nc);
    for (int i = 0; i < nc; ++i) {
      if (my_nnz[i]) {
        int column_master = -1, max_nnz = 0;
        contributors.resize(0);
        
        // Determine who is the master for the column
        for (int j = nbp-1; j >= 0; --j) {
          if (all_nnz(i, j)) {
            if (rk != j) contributors.push_back(j);
            if (column_master == -1 || all_nnz(i, j) > max_nnz)
              { column_master = j; max_nnz = all_nnz(i, j); }
          }
        }

        for (int j = 0; j < int(contributors.size()); ++j) {
          if (column_master == rk) { // receive a column and store
            typename gmm::rsvector<T>::base_type_ &VV = V;
            int si = all_nnz(i, contributors[j]);
            VV.resize(si);
            MPI_Recv((void *)(&(VV[0])),
                     si*sizeof(gmm::elt_rsvector_<T>),
                     MPI_BYTE, contributors[j], 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            gmm::add(V, gmm::mat_col(MM, i));
          } else { // send the column to the column master
            typename gmm::rsvector<T>::base_type_ &VV = MM.col(i);
            MPI_Send((void *)(&(VV[0])),
                     my_nnz[i]*sizeof(gmm::elt_rsvector_<T>),
                     MPI_BYTE, column_master, 0,
                     MPI_COMM_WORLD);
            MM.col(i).clear();
          }
        }
      }
    }

    // Step 3 : gather the total nnz
    for (size_type i = 0; i < nc; ++i) {
      my_nnz[i] = gmm::nnz(gmm::mat_col(MM, i));
      owner[i] = (my_nnz[i]) ? rk : 0;
    }
    MPI_SUM_VECTOR(my_nnz);
    MPI_SUM_VECTOR(owner);
    
    // Step 4 : distribute each column to all the processes
    // Would it be more efficient to perform a global MPI_SUM on a compressed
    // storage ?
    for (size_type i = 0; i < nc; ++i) {
      if (my_nnz[i]) {
        typename gmm::rsvector<T>::base_type_ &VV = MM.col(i);
        if (owner[i] != rk) VV.resize(my_nnz[i]);
        MPI_Bcast((void *)(&(VV[0])), my_nnz[i]*sizeof(gmm::elt_rsvector_<T>),
                  MPI_BYTE, owner[i], MPI_COMM_WORLD);
      }
    }
  }
  
  template <typename MAT> inline void MPI_SUM_SPARSE_MATRIX(const MAT &M) {
    typedef typename gmm::linalg_traits<MAT>::value_type T;
    int nbp; MPI_Comm_size(MPI_COMM_WORLD, &nbp);
    if (nbp < 2) return;
    size_type nr = gmm::mat_nrows(M), nc = gmm::mat_ncols(M);
    gmm::col_matrix< gmm::rsvector<T> > MM(nr, nc);
    GMM_WARNING2("MPI_SUM_SPARSE_MATRIX: A copy of the matrix is done. "
               "To avoid it, adapt MPI_SUM_SPARSE_MATRIX to "
               "your matrix type.");
    gmm::copy(M, MM);
    MPI_SUM_SPARSE_MATRIX(MM);
    gmm::copy(MM, const_cast<MAT &>(M));
  }
#else
  template <typename T> inline T MPI_SUM_SCALAR(T a) { return a; }
  template <typename VECT> inline void MPI_SUM_VECTOR(const VECT &) {}
  template <typename VECT> inline void MPI_MAX_VECTOR(const VECT &) {}
  template <typename T> void MPI_BCAST0_SCALAR(T &a) {}
  template <typename VECT> inline void MPI_BCAST0_VECTOR(const VECT &) {}
  template <typename MAT> inline void MPI_SUM_SPARSE_MATRIX(const MAT &) {}
  template <typename VECT1, typename VECT2>
  inline void MPI_SUM_VECTOR(const VECT1 &V, const VECT2 &WW)
  {
    VECT2 &W = const_cast<VECT2 &>(WW);
    gmm::copy(V, W);
  }
  inline bool MPI_IS_MASTER(void) { return true; }
#endif

  using bgeot::ST_NIL;
  using bgeot::size_type;
  using bgeot::dim_type;
  using bgeot::short_type;
  using bgeot::short_type;
  using bgeot::scalar_type;
  using bgeot::complex_type;
  using bgeot::long_scalar_type;
  using bgeot::opt_long_scalar_type;
  
  using bgeot::base_small_vector;
  using bgeot::base_vector;
  using bgeot::base_complex_vector;
  using bgeot::base_matrix;
  using bgeot::base_complex_matrix;
  using bgeot::base_tensor;
  using bgeot::base_complex_tensor;
  using bgeot::base_poly;
  using bgeot::base_node;
  
  // For compatibility with Getfem 2.0

  using std::invalid_argument;
  using gmm::dimension_error;
  using gmm::file_not_found_error;
  using gmm::internal_error;
  using gmm::to_be_done_error;
  using gmm::failure_error;

#if defined(__GNUC__)
  using std::isnan;
#else
  inline bool isnan(scalar_type x) { return x != x; } 
#endif

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONFIG_H__  */
