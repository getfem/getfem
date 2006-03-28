// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2000-2006 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================

/*! @mainpage Getfem++ reference documentation.

  <center><img src="../logo_getfem_small.png"></center>
 
  @section intro Introduction
 
  This documentation has been generated from the Getfem++ sources, this is not a tutorial.
 
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
  mesh_fem@endlink object in Getfem++.

  A mesh, with a set of integration methods attached to its convexes
  is called a @link getfem_mesh_im.h mesh_im@endlink object in Getfem++.

  A @c mesh_fem can be used to approximate scalar fields (heat, pression,
  ..), or vector fields (displacement, electric field, ..). A @c mesh_im
  will be used to perform numerical integrations on these fields.
  Although Getfem++ supports vector FEMs, all of the FEM currently
  available are scalar (hopefully this will change soon). Of course,
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
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>
   @date November 19, 2000.
   @brief defines and typedefs for namespace getfem
*/


#ifndef GETFEM_CONFIG_H__
#define GETFEM_CONFIG_H__

#include <bgeot_config.h>

// Parallelisation options

// GETFEM_PARA_LEVEL is the parallelisation level of Getfem
//    0 - Sequential
//    1 - Only the resolution of linear systems are parallelized
//    2 - Assembly procedures are also parallelized
#ifndef GETFEM_PARA_LEVEL
# define GETFEM_PARA_LEVEL 0
#endif

#define GETFEM_MPI_INIT(argc, argv) {DAL_TRACE1("Running sequential Getfem");}
#define GETFEM_MPI_FINALIZE {}

#if GETFEM_PARA_LEVEL > 0

# include <mpi.h>

# undef DAL_TRACE_MSG_MPI
# define DAL_TRACE_MSG_MPI					         \
  int mip_rk__; MPI_Comm_rank(MPI_COMM_WORLD, &mip_rk__);	         \
  if (mip_rk__ == 0) 

# undef GETFEM_MPI_INIT
# define GETFEM_MPI_INIT(argc, argv) {                                   \
  MPI_Init(&argc, &argv);                                                \
  DAL_TRACE1("Running parallelized Getfem level " << GETFEM_PARA_LEVEL); \
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

# if GETFEM_PARA_LEVEL > 1
extern "C" void METIS_PartMeshNodal(int *, int *, int *, int *, int *, int *,
				    int *, int *, int *);
extern "C" void METIS_PartGraphKway(int *, int *, int *, int *, int *, int *,
				    int *, int *, int *, int *, int *);
# endif

#endif

#include <bgeot_tensor.h>
#include <bgeot_poly.h>
#include <getfem_superlu.h>

/// GEneric Tool for Finite Element Methods.
namespace getfem {

#if GETFEM_PARA_LEVEL > 1
  template <typename T> inline T MPI_SUM_SCALAR(T a)
  { T b; MPI_Allreduce(&a,&b,1,gmm::mpi_type(a),MPI_SUM,MPI_COMM_WORLD); return b; }
  template <typename VECT> inline void MPI_SUM_VECTOR(VECT V) {
    typedef typename gmm::linalg_traits<VECT>::value_type T;
    std::vector<T> W(gmm::vect_size(V)); gmm::copy(V, W);
    MPI_Allreduce(&(V[0]), &(W[0]), gmm::vect_size(V), gmm::mpi_type(T()),
		  MPI_SUM, MPI_COMM_WORLD);
  }
  inline bool MPI_IS_MASTER(void)
  { int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk); return !rk; }
#else
  template <typename T> inline T MPI_SUM_SCALAR(T a) { return a; }
  template <typename VECT> inline void MPI_SUM_VECTOR(VECT) {}
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
  using bgeot::base_matrix;
  using bgeot::base_tensor;
  using bgeot::base_poly;
  using bgeot::base_node;

  using std::invalid_argument;
  using dal::dimension_error;
  using dal::file_not_found_error;
  using dal::internal_error;
  using dal::to_be_done_error;
  using dal::failure_error;

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_CONFIG_H__  */
