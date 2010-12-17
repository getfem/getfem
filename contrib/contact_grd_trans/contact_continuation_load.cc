// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2009 Yves Renard, Julien  Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

/**
   @file contact_continuation_load.cc
   @brief Nonlinear Problem with Friction(large strain).

   A rubber bar is submitted to a large torsion. If the standard
   solver (Newton) fails, one tries to find a more suitable initial
   approximation by numerical continuation with respect to the
   Dirichlet condition and then restart the standard solver.
   NOTE: The continuation is proposed only for the quasi-static case.
   
   This program is used to check that getfem++ is working. This is
   also a good example of use of Getfem++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_nonlinear_elasticity.h"
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_superlu.h"
#include "gmm/gmm.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector;
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types.  These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
  structure for the frictional problem
*/
struct friction_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, FRICTION_BOUNDARY_NUM = 1, NEUMANN_BOUNDARY_NUM = 2};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim;       /* the integration methods */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure.                   */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  getfem::mesh_fem mf_vm;    /* mesh_fem used for the VonMises stress        */
  scalar_type p1, p2, p3;    /* elastic coefficients.                        */
  scalar_type LX, LY, LZ;    /* system dimensions                            */
  scalar_type lambda, mu;    /* Lame coefficients.                           */
  scalar_type residual;      /* max residual for the iterative solvers       */
  scalar_type rho;           /* density.                                     */
  size_type nocontact_mass;
  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  friction_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh), mf_rhs(mesh), mf_coef(mesh), mf_vm(mesh) {}
};


/* Read parameters from the .param file, build the mesh, set finite element
   and integration methods and selects the boundaries.

   (this is boilerplate code, not very interesting)
 */
void friction_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name for the pressure");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  
  LX = PARAM.real_value("LX", "Length along X axis");
  LY = PARAM.real_value("LY", "Length along Y axis");
  LZ = PARAM.real_value("LZ", "Length along Z axis");
  std::string meshfilename = PARAM.string_value("MESHFILENAME");
  int nb_refine = PARAM.int_value("NBREFINE") ? PARAM.int_value("NBREFINE") : 0;
  scalar_type layerx = PARAM.real_value("LAYERX",
					"Thickness of refinement along the X-axis");
  scalar_type layery = PARAM.real_value("LAYERY",
					"Thickness of refinement along the Y-axis");
  scalar_type layerx_fact = PARAM.real_value("LAYERX_FACT", "Multiplicative factor");
  scalar_type layery_fact = PARAM.real_value("LAYERY_FACT", "Multiplicative factor");
  lambda = PARAM.real_value("LAMBDA", "Lame coefficient lambda");
  mu = PARAM.real_value("MU", "Lame coefficient mu");

  if (meshfilename.size() > 0)
    mesh.read_from_file(meshfilename);
  else {
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(),
	      PARAM.int_value("NX", "Number of space steps "));
    nsubdiv[1] = PARAM.int_value("NY") ? PARAM.int_value("NY") : nsubdiv[0];
    if (N>2) nsubdiv[2] = PARAM.int_value("NZ") ? PARAM.int_value("NZ") : nsubdiv[0];
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			      PARAM.int_value("MESH_NOISED") != 0);

    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    if (N>1) { M(0,1) = PARAM.real_value("INCLINE") * PARAM.real_value("LY"); }
    
    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);
  }

  if (nb_refine > 0) { /* refinement in the right-bottom corner */
    dal::bit_vector cvref;
    for (int i = 0; i < nb_refine; ++i) {
      cvref.clear();
      for (dal::bv_visitor j(mesh.convex_index()); !j.finished(); ++j) {
	if ((mesh.points_of_convex(j)[0][0] > LX * (1 - layerx))
	    && (mesh.points_of_convex(j)[1][0] > LX * (1 - layerx))
	    && (mesh.points_of_convex(j)[2][0] > LX * (1 - layerx))
	    && (mesh.points_of_convex(j)[0][N-1] < LY * layery)
	    &&  (mesh.points_of_convex(j)[1][N-1] < LY * layery)
	    &&  (mesh.points_of_convex(j)[2][N-1] < LY * layery))
	  cvref.add(j);
      }
      mesh.Bank_refine(cvref);
      layerx *= layerx_fact;
      layery *= layery_fact;
    }
  }

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;

  rho = PARAM.real_value("RHO", "Density");
  nocontact_mass = PARAM.int_value("NOCONTACT_MASS", "Suppress the mass "
				   "of contact nodes");

  p1 = PARAM.real_value("P1", "First Elastic coefficient");
  p2 = PARAM.real_value("P2", "Second Elastic coefficient");
  p3 = PARAM.real_value("P3", "Third Elastic coefficient");
  
  mf_u.set_qdim(bgeot::dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mim.set_integration_method(ppi);
  mf_u.set_finite_element(pf_u);

  mf_p.set_finite_element(getfem::fem_descriptor(FEM_TYPE_P));

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM"
		". In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));

  mf_vm.set_classical_discontinuous_finite_element(1);

  /* set boundary conditions
   * (Dirichlet on the upper face, contact on the bottom face 
   *  and on parts of the adjacent faces, Neumann elsewhere) */
      cout << "Selecting the Dirichlet, Neumann and contact boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
    assert(it.is_face());
    base_node un = mesh.normal_of_face_of_convex(it.cv(), it.f());
    un /= gmm::vect_norm2(un);
    if (gmm::abs(un[N-1] - 1.0) < 1.0E-7) { 
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
    } else if (gmm::abs(un[N-1] + 1.0) < 1.0E-7) {
      mesh.region(FRICTION_BOUNDARY_NUM).add(it.cv(),it.f());
    } else if (mesh.points_of_convex(it.cv())[0][N-1] < LZ * 0.2)
      mesh.region(FRICTION_BOUNDARY_NUM).add(it.cv(),it.f());
  }
}

/* The Iteration object for corrections */
class iteration_corr {
  size_type maxiter; /* Max. number of iterations.                       */
  int noise;         /* if noise > 0 iterations are printed.             */
  double resmax;     /* maximum residual.                                */
  double diffmax;    /* maximum difference.                              */
  size_type nit;     /* iteration number.                                */
  double res;        /* last computed residual.                          */
  double diff;       /* last computed difference.                        */
public:
  void init(void) { 
    nit = 0; res = 0.0; diff = 0.0; 
    }

  iteration_corr(double r = 1.0E-8, double d = 1.0E-8, int noi = 0,
		 size_type mit = size_type(-1))
    : maxiter(mit), noise(noi), resmax(r), diffmax(d) { init(); }

  void  operator ++(int) {  nit++; }

  int get_noisy(void) const { return noise; }
  void set_res(double r) { res = r; }
  void set_diff(double d) { diff = d; }
  size_type get_iteration(void) const { return nit; }

  bool converged(void) { return ((res <= resmax) && (diff <= diffmax)); }

  bool finished(void) {
    if (noise > 0)
      cout << "iter " << nit << " residual " << gmm::abs(res)
	   << " difference " << gmm::abs(diff) << endl;
    return ((nit >= maxiter) || converged());
  }
};

/* The object for adaptation of the continuation to non-smoothness of
   the problem by means of active functions */
class active_functions {
  size_type row, nba;
  gmm::dense_matrix<bool> CH_M;
  gmm::dense_matrix<size_type> ACT_M;
  plain_vector limit_val, tangent;
public:
  void clear(void) { row = nba = 0; gmm::resize(tangent, 0); }
  active_functions(void) { clear(); }
  bool empty(void) { return (nba == 0); }

  template<typename MODEL_STATE, typename MATRIX, typename VECT>
  void compute_Jacobians
  (getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const VECT &contact_nodes,
   const plain_vector &UN, const plain_vector &UT, const plain_vector &UT0,
   const plain_vector &LN, const plain_vector &LT, const MATRIX &CH_, scalar_type limit,
   int noisy);

  template<typename MODEL_STATE>
  bool compute_new_tangent
  (MODEL_STATE &MS, getfem::mdbrick_abstract<MODEL_STATE> &final_model,
   getfem::mdbrick_Dirichlet<MODEL_STATE> &DIRICHLET,
   getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, plain_vector grad_XI,
   const plain_vector &F2_init, const plain_vector &deltaF2, const plain_vector &X,
   const plain_vector &T, int noisy);

  template<typename MODEL_STATE>
  void tangent_orientation
  (getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const plain_vector &UN,
   const plain_vector &UT, const plain_vector &UT0, const plain_vector &LN,
   const plain_vector &LT, plain_vector &T);

  plain_vector &get_tangent(void) { return tangent; }
  
};

template<typename MAT>
void compute_condition_number(const MAT &M) {
/* computes the smallest and the largest eigenvalue of the matrix M
   and consequently its condition number */
  
  size_type N = gmm::mat_nrows(M);
  plain_vector V(N), U(N), W(N);
  plain_vector V2(N), U2(N);  
  gmm::fill_random(V);
  gmm::fill_random(V2);
  scalar_type lambda = 0., lambda2 = 0.;
  gmm::SuperLU_factor<scalar_type> SLU;
  SLU.build_with(M);
  
  for(size_type i = 0; i < 100; ++i) {
    gmm::scale(V, 1./gmm::vect_norm2(V));
    gmm::scale(V2, 1./gmm::vect_norm2(V2));
    SLU.solve(U, V);
    gmm::mult(M, V2, U2);
    
    lambda = gmm::vect_sp(U, V) / gmm::vect_sp(V, V);
    lambda2 = gmm::vect_sp(U2, V2) / gmm::vect_sp(V2, V2);
    
    gmm::add(V, gmm::scaled(U, -1./lambda), W);
//     if (i > 90)
//       cout << "iteration " << i << " norm(W) = " <<  gmm::vect_norm2(W)
// 	   << " lambda = " << lambda 
// 	   << " lambda2 = " << lambda2 << endl;
    if (gmm::vect_norm2(W) < 1E-2) break;
    
    gmm::copy(U, V);
    gmm::copy(U2, V2);
  }
  
  // cout << "Smallest eigenvalue : " << 1.0/lambda << endl;
  // cout << "Largest eigenvalue : " << lambda2 << endl;
  cout << "Condition number: " << lambda2 * lambda << endl;
  
}

template<typename MAT, typename VECT1, typename VECT2>
void compute_tangent(const MAT &A, const VECT1 &B, VECT2 &T, int noisy) {
/* computes a unit tangent at the point X to the curve given implicitly by the equation H = 0
   with the Jacobian grad H(X) = (A|B) which is positively oriented w.r.t. the vector T
   and saves it into T, i.e.:
   solves (A|B)T_+ = 0 & T.T_+ = 1;
   sets T = T_+ / ||T_+||; */

  bool singular_A = false;
  size_type N = gmm::mat_nrows(A);
  double rcond;
		
  try { /* whether A is singular; if not, the calculation can be accelerated provided that 
	   A is sparse  */
    scalar_type sp;
    plain_vector Y(N);
    if (noisy > 1) cout << "starting computing tangent" << endl;

    gmm::SuperLU_solve(A, Y, B, rcond);
    sp = gmm::vect_sp(gmm::sub_vector(T, gmm::sub_interval(0,N)), Y);
		  
    T[N] = 1.0 / (T[N] - sp);
    gmm::copy(gmm::scaled(Y, -1.0*T[N]),
	      gmm::sub_vector(T, gmm::sub_interval(0, N)));
  }
  catch (...) {
    singular_A = true;
  }
  
  if (singular_A) {
    sparse_matrix C(N+1, N+1);
    plain_vector R(N+1);

    gmm::copy(A, gmm::sub_matrix(C, gmm::sub_interval(0, N), gmm::sub_interval(0, N)));
    gmm::copy(gmm::col_vector(B),
	      gmm::sub_matrix(C, gmm::sub_interval(0, N), gmm::sub_interval(N, 1)));
    gmm::copy(gmm::row_vector(T),
	      gmm::sub_matrix(C, gmm::sub_interval(N, 1), gmm::sub_interval(0, N+1)));

    gmm::clear(R);
    R[N] = 1.0;		  
    
    gmm::SuperLU_solve(C, T, R, rcond);
  }

  gmm::scale(T, 1.0/gmm::vect_norm2(T));

}

template<typename MAT, typename VECT1, typename VECT2>
void correction_solver
(const MAT &A, const VECT1 &B, VECT2 &X, VECT2 &T, const VECT1 &P, iteration_corr &iter) {
/* computes one correction step of the Moore-Penrose method for continuing the implicit curve
   H = 0 with P = H(X), (A|B) = grad H(X), X and T being the approximations of the point
   and the tangent, respectively, i.e.:
   solves (A|B)T_+ = 0 & T.T_+ = 1;
   solves (A|B)dX = P; T.dX = 0;
   sets X = X - dX;
   sets T = T_+ / ||T_+||;  */

  bool singular_A = false;
  size_type N = gmm::mat_nrows(A);
  plain_vector dX(N+1);

  try { /* whether A is singular; if not, the calculation can be accelerated provided that 
	   A is sparse */
    scalar_type sp1, sp2;
    plain_vector Y1(N), Y2(N);

    compute_condition_number(A);

    gmm::SuperLU_factor<scalar_type> SLU;
    SLU.build_with(A);
    SLU.solve(Y1, B);
    SLU.solve(Y2, P);

    sp1 = gmm::vect_sp(gmm::sub_vector(T, gmm::sub_interval(0,N)), Y1);
    sp2 = gmm::vect_sp(gmm::sub_vector(T, gmm::sub_interval(0,N)), Y2);
		
    dX[N] = sp2 / (sp1 - T[N]);
    gmm::add(gmm::scaled(Y1, -1.0*dX[N]), Y2);
    gmm::copy(Y2, gmm::sub_vector(dX, gmm::sub_interval(0, N)));
		
    T[N] = 1.0 / (T[N] - sp1);
    gmm::copy(gmm::scaled(Y1, -1.0*T[N]),
	      gmm::sub_vector(T, gmm::sub_interval(0, N)) );
  }
  catch (...) {
    singular_A = true;
  }


  try {
    if (singular_A) {
      plain_vector Q(N+1), R(N+1);
      sparse_matrix C(N+1, N+1);

      gmm::copy(A, gmm::sub_matrix(C, gmm::sub_interval(0, N), gmm::sub_interval(0, N)));
      gmm::copy(gmm::col_vector(B),
		gmm::sub_matrix(C, gmm::sub_interval(0, N), gmm::sub_interval(N, 1)));
      gmm::copy(gmm::row_vector(T),
		gmm::sub_matrix(C, gmm::sub_interval(N, 1), gmm::sub_interval(0, N+1)));
		
      gmm::copy(P, gmm::sub_vector(Q, gmm::sub_interval(0, N)));
      Q[N] = 0.0;

      gmm::clear(R);
      R[N] = 1.0;	
		    
      gmm::SuperLU_factor<scalar_type> SLU;
      SLU.build_with(C);
      SLU.solve(dX, Q);
      SLU.solve(T, R);
    }
  }
  catch (...) {
    cout << "SuperLU failed!" << endl;
  }
    
  gmm::add(gmm::scaled(dX, -1.0), X);
  gmm::scale(T, 1.0/gmm::vect_norm2(T));

  iter.set_diff(gmm::vect_norm2(dX));

}

template<typename MODEL_STATE>
void standard_solver
(MODEL_STATE &MS, getfem::mdbrick_abstract<MODEL_STATE> &final_model,
 getfem::mdbrick_Dirichlet<MODEL_STATE> &DIRICHLET, const plain_vector &F2,
gmm::iteration &iter) {
  /* standard solver (Newton) */

  DIRICHLET.rhs().set(F2);
// cout << "F2 = " << F2 << endl;
// cout << "|U| = " << gmm::vect_norm2(MS.state()) << "\n";

  gmm::basic_newton_line_search ls(size_type(-1), 5.0/3.0, 0.1, 0.5, 3.0);
//   gmm::default_newton_line_search ls;
  getfem::standard_solve(MS, final_model, iter,
			 getfem::default_linear_solver(final_model), ls);
}

template<typename MODEL_STATE>
void Newton_correction
(MODEL_STATE &MS, getfem::mdbrick_abstract<MODEL_STATE> &final_model,
 getfem::mdbrick_Dirichlet<MODEL_STATE> &DIRICHLET, const plain_vector &grad_XI,
 const plain_vector &F2_init, const plain_vector &deltaF2, plain_vector &X, plain_vector &T,
 iteration_corr &iter) {

  size_type stot = gmm::mat_ncols(MS.tangent_matrix());
  plain_vector F2 = F2_init;

  gmm::copy(gmm::sub_vector(X, gmm::sub_interval(0, stot)), MS.state());
  gmm::add(gmm::scaled(deltaF2, X[stot]), F2);
  DIRICHLET.rhs().set(F2);

  final_model.compute_residual(MS);

  do {
    if (iter.get_noisy() > 1)
      cout << "starting computing tangent matrix" << endl;
    final_model.compute_tangent_matrix(MS);

    if (iter.get_noisy() > 1)
      cout << "starting linear solver" << endl;
    correction_solver(MS.tangent_matrix(), grad_XI, X, T, MS.residual(), iter);
    if (iter.get_noisy() > 1) 
      cout << "linear solver done" << endl;

    gmm::copy(gmm::sub_vector(X, gmm::sub_interval(0, stot)), MS.state());
    gmm::copy(F2_init, F2);
    gmm::add(gmm::scaled(deltaF2, X[stot]), F2);
    DIRICHLET.rhs().set(F2);

    final_model.compute_residual(MS);
    iter.set_res(gmm::vect_norm2(MS.residual()));

    iter++;
  } while (!iter.finished());
}

template<typename MODEL_STATE>
void compute_displacement
(MODEL_STATE &MS, getfem::mdbrick_abstract<MODEL_STATE> &final_model,
 getfem::abstract_hyperelastic_law &l,
 getfem::mdbrick_nonlinear_elasticity<MODEL_STATE> &ELAS_nonlin,
 getfem::mdbrick_isotropic_linearized_elasticity<MODEL_STATE> &ELAS_lin, plain_vector &U,
 size_type law_num) {

  if (law_num < 4) {
    l.reset_unvalid_flag();
    final_model.compute_residual(MS);
    if (l.get_unvalid_flag()) 
      GMM_WARNING1("The solution is not completely valid, the determinant "
		   "of the transformation is negative on "
		   << l.get_unvalid_flag() << " Gauss points");
    
    gmm::copy(ELAS_nonlin.get_solution(MS), U);
  } else
    gmm::copy(ELAS_lin.get_solution(MS), U);
}

template<typename MODEL_STATE>
void compute_Von_Mises
(MODEL_STATE &MS, getfem::mdbrick_nonlinear_elasticity<MODEL_STATE> &ELAS_nonlin,
 getfem::mdbrick_isotropic_linearized_elasticity<MODEL_STATE> &ELAS_lin,
 getfem::mesh_fem &mf_vm, plain_vector &VM, size_type law_num) {
  /* computes the Von Mises stress from the solution */
  if (law_num < 4)
    ELAS_nonlin.compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
  else 
    ELAS_lin.compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
}

template<typename MODEL_STATE, typename VECT>
void write_charact_solution
(getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const VECT &contact_nodes,
 const plain_vector &UN, const plain_vector &UT, const plain_vector &UT0,
 const plain_vector &LN, const plain_vector &LT) {
  
  scalar_type r = FRICTION.get_r();
  scalar_type alpha = FRICTION.get_alpha(), beta = FRICTION.get_beta();
  plain_vector gap(LN.size()), friction_coef(LN.size());
    
  gmm::copy(FRICTION.get_gap(), gap);
  gmm::copy(FRICTION.get_friction_coef(), friction_coef);

  cout << "characters of the solution:" << endl;

  for (size_type i = 0; i < LN.size(); ++i) {
    cout << "node " << i << ": " << contact_nodes[i] << " " 
	 << (LN[i] - r * alpha * (UN[i] - gap[i]) <= 0)
	 << (- friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]) >= 0)
	 << (friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]) <= 0) << endl;

// 	cout << "LN - r * alpha * (UN - gap) = " << LN[i] - r * alpha * (UN[i] - gap[i])
// 	     << ", - friction_coef * LN + LT -  r * beta * (UT - UT0) = "
// 	     << -friction_coef[i] * LN[i] + LT[i] -  r * beta * (UT[i]-UT0[i])
// 	     << ", friction_coef * LN + LT -  r * beta * (UT - UT0) = "
// 	     << friction_coef[i] * LN[i] + LT[i] -  r * beta * (UT[i]-UT0[i]) << endl;
// 	cout << "LN = " << LN[i] 
// 	     << ", r * alpha * (UN - gap) = " << r * alpha * (UN[i] - gap[i])
// 	     << ", - friction_coef * LN + LT = "
// 	     << -friction_coef[i] * LN[i] + LT[i]
// 	     << ", friction_coef * LN + LT = "
// 	     << friction_coef[i] * LN[i] + LT[i] 
// 	     << ", r * beta * (UT - UT0) = "
// 	     << r * beta * (UT[i]-UT0[i]) << endl;
  }
}

template<typename MODEL_STATE, typename MATRIX, typename VECT>
void compute_charact_solution
(getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const VECT &contact_nodes,
 const plain_vector &UN, const plain_vector &UT, const plain_vector &UT0,
 const plain_vector &LN, const plain_vector &LT, MATRIX &CH, int noisy) {
  
  scalar_type r = FRICTION.get_r();
  scalar_type alpha = FRICTION.get_alpha(), beta = FRICTION.get_beta();
  plain_vector gap(LN.size()), friction_coef(LN.size());
  base_vector test_fcn(3);
    
  gmm::copy(FRICTION.get_gap(), gap);
  gmm::copy(FRICTION.get_friction_coef(), friction_coef);

  if (noisy > 1)
    write_charact_solution(FRICTION, contact_nodes, UN, UT, UT0, LN, LT);

  for (size_type i = 0; i < LN.size(); ++i) {
    CH(i, 0) = (LN[i] - r * alpha * (UN[i] - gap[i]) <= 0);
    CH(i, 1) = (- friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]) >= 0);
    CH(i, 2) = (friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]) <= 0);
  }

}

template<typename MODEL_STATE, typename MATRIX, typename VECT>
size_type compute_charact_solution
(getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const VECT &contact_nodes,
 const plain_vector &UN, const plain_vector &UT, const plain_vector &UT0,
 const plain_vector &LN, const plain_vector &LT, MATRIX &CH, const MATRIX &CH0, int noisy) {
  
  size_type nb_change = 0;
  bool change;

  compute_charact_solution(FRICTION, contact_nodes, UN, UT, UT0, LN, LT, CH, noisy);

  for (size_type i = 0; i < LN.size(); ++i) {
    change = false;
    for (size_type j = 0; j < 3; ++j) {
      if ((CH(i, j) != CH0(i, j)) && (j == 0 || (j > 0 && CH(i, 0))))
	change = true;
    }

    if (change) {
      if (nb_change == 0)
	cout << "changes of characters:" << endl;
      cout << "node " << i << ": " << contact_nodes[i] << " "
	   << CH0(i, 0) << CH0(i, 1) << CH0(i, 2) << " -> "
	   << CH(i, 0) << CH(i, 1) << CH(i, 2) << endl;
      ++nb_change;
    }
  }
  
  return nb_change;
}

template<typename MATRIX>
void straight_insertion
(MATRIX &M, plain_vector &key, size_type i) {
/* places the i-th row of M according to the absolut values of the corresponding components
   in key; increasing order is wanted in the first i components of key */

  bool found;
  size_type j = i;
  std::vector<size_type> X(2);
  scalar_type X_key;

  X[0] = M(i, 0); X[1] = M(i, 1);
  X_key = key[i];
  if (j == 0) found = true;
  else found = (gmm::abs(X_key) >= gmm::abs(key[j - 1]));
 
  while (!found) { // seeking the appropriate  position
    key[j] = key[j - 1];
    M(j, 0) = M(j - 1, 0);
    M(j, 1) = M(j - 1, 1);
    --j;
    if (j == 0) found = true;
    else found = (gmm::abs(X_key) >= gmm::abs(key[j - 1]));
  }

  M(j, 0) = X[0]; M(j, 1) = X[1];
  key[j] = X_key;
}

template<typename MODEL_STATE, typename MATRIX, typename VECT>
void active_functions::compute_Jacobians
(getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const VECT &contact_nodes,
 const plain_vector &UN, const plain_vector &UT, const plain_vector &UT0,
 const plain_vector &LN, const plain_vector &LT, const MATRIX &CH_, scalar_type limit,
 int noisy) {

  scalar_type r = FRICTION.get_r();
  scalar_type alpha = FRICTION.get_alpha(), beta = FRICTION.get_beta();
  plain_vector gap(LN.size()), friction_coef(LN.size());
  base_vector test_fcn(3);
    
  gmm::copy(FRICTION.get_gap(), gap);
  gmm::copy(FRICTION.get_friction_coef(), friction_coef);

  if (noisy > 1)
    write_charact_solution(FRICTION, contact_nodes, UN, UT, UT0, LN, LT);

  if (gmm::mat_nrows(CH_M) == 0) {
    gmm::resize(CH_M, gmm::mat_nrows(CH_), gmm::mat_ncols(CH_));
    gmm::resize(ACT_M, 3 * gmm::mat_nrows(CH_), 2);
    gmm::resize(limit_val, 3 * gmm::mat_nrows(CH_));
  }
  
  gmm::copy(CH_, CH_M);

  nba = 0;

  for (size_type i = 0; i < LN.size(); ++i) {
    test_fcn[0] = LN[i] - r * alpha * (UN[i] - gap[i]);
    test_fcn[1] = - friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);
    test_fcn[2] = friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);
    
    for (size_type j = 0; j < 3; ++j)
      if (gmm::abs(test_fcn[j]) <= limit) {
	if (nba == 0)
	  cout << "test functions with values close to zero:" << endl;
	cout << "test_fcn(" << i << ", " << j << ") = " << test_fcn[j]
	     << " (node " << i << ", " << contact_nodes[i] << " " 
// 	     << "LN = " << LN[i] 
// 	     << ", r * alpha * (UN - gap) = " << r * alpha * (UN[i] - gap[i])
// 	     << ", - friction_coef * LN + LT = "
// 	     << -friction_coef[i] * LN[i] + LT[i]
// 	     << ", friction_coef * LN + LT = "
// 	     << friction_coef[i] * LN[i] + LT[i] 
// 	     << ", r * beta * (UT - UT0) = "
// 	     << r * beta * (UT[i]-UT0[i])
	     << CH_M(i, 0) << CH_M(i, 1) << CH_M(i, 2) << ")" << endl;
	
	ACT_M(nba, 0) = i; ACT_M(nba, 1) = j;
	limit_val[nba] = test_fcn[j];
	straight_insertion(ACT_M, limit_val, nba);
	++nba;
      }
  }

  if (nba == 0)
    cout << "no test functions with values close to zero were founded " << endl;
}

template<typename MODEL_STATE>
bool active_functions::compute_new_tangent
(MODEL_STATE &MS, getfem::mdbrick_abstract<MODEL_STATE> &final_model,
 getfem::mdbrick_Dirichlet<MODEL_STATE> &DIRICHLET,
 getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, plain_vector grad_XI,
 const plain_vector &F2_init, const plain_vector &deltaF2, const plain_vector &X,
 const plain_vector &T, int noisy) {
  /* cames through the proposed Jacobians with a change exactly at one character */

  bool new_tangent = false;
  size_type stot = gmm::mat_ncols(MS.tangent_matrix());
  plain_vector F2 = F2_init;
  
  if (row < nba) {
    gmm::copy(gmm::sub_vector(X, gmm::sub_interval(0, stot)), MS.state());
    gmm::add(gmm::scaled(deltaF2, X[stot]), F2);
    DIRICHLET.rhs().set(F2);
    
    while(!new_tangent && row < nba) {
      size_type i = ACT_M(row, 0), j = ACT_M(row, 1);
      
      CH_M(i, j) = !CH_M(i, j);
      
      if((CH_M(i, 0) || j == 0) && (CH_M(i, 1) || CH_M(i, 2))) {
	/* the character can be accepted */
	cout << "trying the following character of node " << i << ": "
	     << CH_M(i, 0) << CH_M(i, 1) << CH_M(i, 2) << endl;
	
	FRICTION.set_character_matrix(CH_M);
	final_model.compute_tangent_matrix(MS);
	
	if (gmm::vect_size(tangent) == 0)
	  gmm::resize(tangent, stot + 1);
	gmm::copy(T, tangent);
	compute_tangent(MS.tangent_matrix(), grad_XI, tangent, noisy);
	cout << "T0.tangent = " << gmm::vect_sp(T, tangent) << ", ";
	new_tangent = true;
      }
      
      CH_M(i, j) = !CH_M(i, j);
      ++row;
    }
    
    FRICTION.clear_character_matrix();
  }

  return new_tangent;

}

template<typename MODEL_STATE>
void active_functions::tangent_orientation
(getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const plain_vector &UN,
 const plain_vector &UT, const plain_vector &UT0, const plain_vector &LN,
 const plain_vector &LT, plain_vector &T) {

  size_type i = ACT_M(row - 1, 0), j = ACT_M(row - 1, 1);
  scalar_type r = FRICTION.get_r(), alpha = FRICTION.get_alpha(), beta = FRICTION.get_beta(),
    test_fcn;
  plain_vector gap(LN.size()), friction_coef(LN.size());
    
  gmm::copy(FRICTION.get_gap(), gap);
  gmm::copy(FRICTION.get_friction_coef(), friction_coef);

  switch (j) {
  case 0: test_fcn = LN[i] - r * alpha * (UN[i] - gap[i]); break;
  case 1: test_fcn = - friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]); break;
  default: test_fcn = friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]); break;
  }

  cout << "test_fcn(i, j)(tangent) = " << test_fcn - limit_val[row - 1] << ", ";
//   cout << " LN = " << LN[i] 
//        << ", r * alpha * (UN - gap) = " << r * alpha * (UN[i] - gap[i]) << ", ";
// 	     << ", - friction_coef * LN + LT = "
// 	     << -friction_coef[i] * LN[i] + LT[i]
// 	     << ", friction_coef * LN + LT = "
// 	     << friction_coef[i] * LN[i] + LT[i] 
// 	     << ", r * beta * (UT - UT0) = "
// 	     << r * beta * (UT[i]-UT0[i])
  if (limit_val[row - 1] * (test_fcn - limit_val[row - 1])  > 0) {
    gmm::scale(T, -1.0);
    cout << "T0 := -tangent" << endl;
  } else
    cout << "T0 := tangent" << endl;
}
  
/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool friction_problem::solve(plain_vector &U) {

  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  size_type law_num = PARAM.int_value("LAW");
  // Linearized elasticity brick.
  base_vector p(3); p[0] = p1; p[1] = p2; p[2] = p3;
 
  /* choose the material law */
  getfem::abstract_hyperelastic_law *pl = 0;
  switch (law_num) {
    case 0:
    case 1: pl = new getfem::SaintVenant_Kirchhoff_hyperelastic_law(); break;
    case 2: pl = new getfem::Ciarlet_Geymonat_hyperelastic_law(); break;
    case 3: pl = new getfem::Mooney_Rivlin_hyperelastic_law(); break;
    case 4 : case 5 : break;
    default: GMM_ASSERT1(false, "no such law");
  }
  
  getfem::mdbrick_abstract<> *pELAS = 0;
  getfem::mdbrick_nonlinear_elasticity<> *pELAS_nonlin = 0;
  getfem::mdbrick_isotropic_linearized_elasticity<> *pELAS_lin = 0;
  
  if (law_num < 4) {
    p.resize(pl->nb_params());
    pELAS = pELAS_nonlin
      = new getfem::mdbrick_nonlinear_elasticity<>(*pl, mim, mf_u, p);
  } else {
    pELAS = pELAS_lin = new getfem::mdbrick_isotropic_linearized_elasticity<>
      (mim, mf_u, law_num == 5 ? 0.0 : lambda, mu);
  }
  
  getfem::mdbrick_abstract<> *pINCOMP = pELAS;
  switch (law_num) {
    case 1: 
    case 3: pINCOMP = new getfem::mdbrick_nonlinear_incomp<>(*pELAS, mf_p);
      break;
    case 5: {
      getfem::mdbrick_linear_incomp<> *pb = 
        new getfem::mdbrick_linear_incomp<>(*pELAS, mf_p);
      pINCOMP = pb;
      pb->penalization_coeff().set(1.0/lambda);
    }
  }

  size_type nb_step = int(PARAM.int_value("NBSTEP"));
  scalar_type deltat = PARAM.real_value("DELTAT", "Time step");
 
  // contact condition for the Lagrange elements
  dal::bit_vector cn = mf_u.basic_dof_on_region(FRICTION_BOUNDARY_NUM);
//   cout << "cn = " << cn << endl;
//   cout << "cn.card()/N = " << cn.card()/N << endl;
  sparse_matrix BN(cn.card()/N, mf_u.nb_dof());
  sparse_matrix BT((N-1)*cn.card()/N, mf_u.nb_dof());
  std::vector<base_node> contact_nodes;
  plain_vector gap(cn.card()/N);
  size_type jj = 0;
  for (dal::bv_visitor i(cn); !i.finished(); ++i)
    if (i % N == 0) {
//       cout << "point de contact " << mf_u.point_of_dof(i) << endl;
      BN(jj, i+N-1) = -1.;
      gap[jj] = mf_u.point_of_basic_dof(i)[N-1];
      contact_nodes.push_back(mf_u.point_of_basic_dof(i));
      for (size_type k = 0; k < N-1; ++k) BT((N-1)*jj+k, i+k) = 1.;
      ++jj;
    }
  
    // creating force density vectors
  int nbc = int(jj);
//   sparse_matrix MMBN(nbc, nbc), MMBT(nbc*(N-1), nbc*(N-1));
  plain_vector UN(nbc), UT(nbc*(N-1)), LN1(nbc), LT1(nbc*(N-1));
//   {
//     sparse_matrix BB(mf_u.nb_dof(), mf_u.nb_dof());
//     getfem::asm_mass_matrix(BB, mim, mf_u, mf_u, FRICTION_BOUNDARY_NUM);
//     std::vector<size_type> indN, indT;
//     for (dal::bv_visitor i(cn); !i.finished(); ++i)
//       if ((i%N) == N-1) indN.push_back(i); else indT.push_back(i);
//     gmm::sub_index SUBI(indN);
//     gmm::copy(gmm::sub_matrix(BB, SUBI, SUBI), MMBN);
//     gmm::sub_index SUBJ(indT);
//     gmm::copy(gmm::sub_matrix(BB, SUBJ, SUBJ), MMBT);    
//   }



  scalar_type friction_coef = PARAM.real_value("FRICTION_COEFF",
					       "Friction cefficient");
  scalar_type r = PARAM.real_value("R", "Augmentation parameter");
  scalar_type alpha = PARAM.real_value("ALPHA") ? PARAM.real_value("ALPHA") : 1.0;
  

  getfem::mdbrick_Coulomb_friction<> FRICTION(*pINCOMP, BN, gap,
					      friction_coef, BT);
  FRICTION.set_r(r);
  FRICTION.set_alpha(alpha);
  FRICTION.set_beta(1./deltat);
 
  // Defining the volumic source term.
  base_vector f(N);
  f[0] = PARAM.real_value("FORCEX","Amplitude of the gravity");
  f[1] = PARAM.real_value("FORCEY","Amplitude of the gravity");
  if (N>2)
    f[2] = PARAM.real_value("FORCEZ","Amplitude of the gravity");
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i) {
    gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  }

  getfem::mdbrick_source_term<> VOL_F(FRICTION, mf_rhs, F);

  // Dirichlet condition
  plain_vector F2(nb_dof_rhs * N);
  getfem::mdbrick_Dirichlet<> DIRICHLET(VOL_F, DIRICHLET_BOUNDARY_NUM);
  DIRICHLET.rhs().set(mf_rhs, F2);
  DIRICHLET.set_constraints_type(getfem::constraints_type
				   (PARAM.int_value("DIRICHLET_VERSION")));

  getfem::mdbrick_abstract<> *pfinal_model = &DIRICHLET;
  
  getfem::mdbrick_dynamic<> *pDYNAMIC = 0;

  // Eventual dynamic brick
  if (rho > 0) {
    pfinal_model = pDYNAMIC = new getfem::mdbrick_dynamic<>(DIRICHLET, rho);
    pDYNAMIC->set_dynamic_coeff(1./(deltat*deltat), 1.);
    if (nocontact_mass)
      pDYNAMIC->no_mass_on_boundary( FRICTION_BOUNDARY_NUM, nocontact_mass == 2);
  }

  // Generic solver.
  getfem::standard_model_state MS(*pfinal_model);
  size_type step0 = PARAM.int_value("STEP0") ? PARAM.int_value("STEP0") : 0;
  size_type step0_cont = PARAM.int_value("STEP0_CONT") ?
    PARAM.int_value("STEP0_CONT") : 0;
  bool X0_given = PARAM.int_value("X0_GIVEN") ? (PARAM.int_value("X0_GIVEN") != 0) : false;
  size_type maxit = PARAM.int_value("MAXITER"); 
  gmm::iteration iter;

  scalar_type dy = PARAM.real_value("DIRICHLET_Y",
				    "Prescribed displacement in y");
  scalar_type dxv = PARAM.real_value("DIRICHLET_X_SPEED",
				     "Prescribed velocity in x");
  scalar_type limit = PARAM.real_value("LIMIT", "limit for the test functions");
  int noisy = PARAM.int_value("NOISY");


  plain_vector F2_init = F2;
  plain_vector U0 = U, U_init = U, UT0 = UT, UT_init = UT,
    V0(gmm::vect_size(U)), DF(gmm::vect_size(U));
  short convergence = -1;
//  plain_vector SumN(nb_step), SumT(nb_step),Pressure(nb_step);
//  plain_vector maxT(nb_step);
//  plain_vector maxN(nb_step);
//  plain_vector max(nb_step);
//  double Fmax =0;

  if (step0 > 0){ /* load the foregoing values */
    char s[100]; sprintf(s, "step%d", step0);
    gmm::vecload(datafilename + s + ".U", U0);
    gmm::vecload(datafilename + s + ".Y", MS.state());
    if (rho > 0)
      gmm::vecload(datafilename + s + ".V", V0);
  }
 

  for (size_type step = step0; step < nb_step; ++step) {
    cout << "beginning of step " << step+1 
	     << ", number of variables : " << pfinal_model->nb_dof() << endl ;

    if (rho > 0) {
      gmm::mult(pDYNAMIC->get_M(), gmm::scaled(U0, 1./(deltat*deltat)), DF);
      gmm::mult_add(pDYNAMIC->get_M(), gmm::scaled(V0, 1./deltat), DF);
      pDYNAMIC->set_DF(DF);
    }

    convergence = -1;

    do { /* seek the solution */

      /* in the first eleven iterations, the body is compressed,
	 afterwards it us pulling to the right while subject to a constant pressure */
      for (size_type i = 0; i < nb_dof_rhs; ++i) {
	F2[i*N+N-1] = (step < 10) ? (dy * step/10.0) : dy;
	//  F2[i*N+N-1] = dy;
	F2[i*N+N-2] = (step < 10) ? 0.0 : (step-10)*deltat*dxv;
	// F2[i*N+N-1] = dy+dy*3*step/nb_step;
      }

      FRICTION.set_WT(gmm::scaled(U0, -1.0));

      if (step0_cont ==  0 && !X0_given) { /* the standard solver */
	iter = gmm::iteration(residual, noisy, maxit ? maxit : 40000);
	standard_solver(MS, *pfinal_model, DIRICHLET, F2, iter);

	if (iter.converged())
	  convergence = 1;
	else {
	  cout << "the standard solver has failed, ";
	  if (convergence == 0) /* the initial approximation was given by continuation */
	    convergence = - 2;
	}
	  
      } // the standard solver

      if (convergence == -1) {
	/* the solution has not been found by the standard solver; 
	   proceed with continuation --- the continuation parameter XI is such that 
	   the actual Dirichlet condition (load) is F2_init + XI * (F2 - F2_init) */

	cout << "starting to continue";
	GMM_ASSERT1(PARAM.int_value("DIRICHLET_VERSION") == 0, "The continuation only "
		    "implemented for the Dirichlet condition with multipliers");
	GMM_ASSERT1(rho <= 0, "The continuation is proposed only for the quasi-static case");

	size_type step_dec = PARAM.int_value("STEP_DEC") ? PARAM.int_value("STEP_DEC") : 0;
	size_type step_init = step - step_dec;
	cout << ", initial step for continuation = " << step_init << endl;
	scalar_type difference = PARAM.real_value("DIFFERENCE");
	if (difference == 0.) difference = 1e-10;
	scalar_type minangle = PARAM.real_value("ANGLE");
	scalar_type XI_end = PARAM.real_value("XI_END");
	scalar_type maxdist = PARAM.real_value("DISTANCE");
	size_type maxit_corr = PARAM.int_value("MAXITER_CORR");
	size_type thr_corr = PARAM.int_value("THRESHOLD_CORR");
	size_type nb_step_cont = PARAM.int_value("NBSTEP_CONT");

	size_type n = pfinal_model->nb_dof();
	scalar_type h_init = PARAM.real_value("H_INIT", "Initial step length") * n;
	scalar_type h_max = PARAM.real_value("H_MAX", "Maximal step length") * n;
	scalar_type h_min = PARAM.real_value("H_MIN", "Minimal step length") * n;
	scalar_type h_inc = PARAM.real_value("H_INC");
	scalar_type h_dec = PARAM.real_value("H_DEC");
	scalar_type h_change = PARAM.real_value("H_CHANGE") * n;
	scalar_type h = PARAM.real_value("H") ? PARAM.real_value("H") * n : h_init;


	/* set the initial values and compute the gradient with respect to XI */
	scalar_type XI0;
	size_type stot = gmm::mat_ncols(MS.tangent_matrix());
	plain_vector Y0(stot);
	plain_vector X0(stot+1), X(stot+1), T0(stot+1), T(stot+1);
	gmm::dense_matrix<bool> CH(nbc, 3), CH0(nbc, 3);
	plain_vector deltaF2 = F2;
	plain_vector grad_XI(stot);

	if (step_init > 0) {
	  char s[100]; sprintf(s, "step%d", step_init);
	  gmm::vecload(datafilename + s + ".F2", F2_init);
	  gmm::vecload(datafilename + s + ".U", U_init);
	  FRICTION.set_WT(gmm::scaled(U_init, -1.0));
	}

	size_type ind =  DIRICHLET.first_ind();
	size_type sc = gmm::vect_size(DIRICHLET.get_CRHS());
	gmm::add(gmm::scaled(F2_init, -1.0), deltaF2);
	DIRICHLET.rhs().set(deltaF2);
	pfinal_model->compute_tangent_matrix(MS);
	gmm::copy(DIRICHLET.get_CRHS(),
		  gmm::sub_vector(grad_XI, gmm::sub_interval(ind, sc)));
	gmm::scale(grad_XI, -1.0);

	if (X0_given || step0_cont > 0) {
	  char s1[100]; sprintf(s1, "step%d", step + 1);
	  char s2[100]; sprintf(s2, "%d", step0_cont);
	  gmm::vecload(datafilename + s1 + "_" + s2 + ".X", X0);
	  gmm::copy(gmm::sub_vector(X0, gmm::sub_interval(0, stot)), MS.state());
	  XI0 = X0[stot];
	  if (step0_cont > 0)
	    gmm::vecload(datafilename + s1 + "_" + s2 + ".T", T0);
	}

	if (step0_cont == 0){

	  if (!X0_given) {
	    cout << "starting computing an initial point" << endl;
	    if (step_init > 0) {
	      char s[100]; sprintf(s, "step%d", step_init);
	      gmm::vecload(datafilename + s + ".Y", Y0);
	    }
	    gmm::copy(Y0, MS.state()); XI0 = maxdist;
	  } else {
	    cout << "correcting the given initial point" << endl;
	    XI0 /= h_dec;
	  }

	  do {
	    XI0 *= h_dec; 
	    cout << "XI0 = " << XI0 << endl;
	    gmm::copy(F2_init, F2);
	    gmm::add(gmm::scaled(deltaF2, XI0), F2);
	    iter = gmm::iteration(residual, noisy, maxit_corr ? maxit_corr : 40000);
	    standard_solver(MS, *pfinal_model, DIRICHLET, F2, iter);
	  } while (!iter.converged());
	  
	  gmm::copy(MS.state(), gmm::sub_vector(X0, gmm::sub_interval(0, stot)));
	  X0[stot] = XI0;
	  
	  pfinal_model->compute_tangent_matrix(MS); 
	  gmm::fill_random(T0);
	  compute_tangent(MS.tangent_matrix(), grad_XI, T0, noisy);
	  if (T0[stot] < 0)
	    gmm::scale(T0, -1.0);
	}

	compute_displacement(MS, *pfinal_model, *pl, *pELAS_nonlin, *pELAS_lin, U, law_num);
	gmm::mult(BN, U, UN); gmm::mult(BT, U, UT); gmm::mult(BT, U_init, UT_init);
	gmm::copy(FRICTION.get_LN(MS), LN1); gmm::copy(FRICTION.get_LT(MS), LT1);
	compute_charact_solution(FRICTION, contact_nodes, UN, UT, UT_init, LN1, LT1, CH0,
				 noisy + 1);

	if (step0_cont == 0){
	  char s1[100]; sprintf(s1, "step%d", step + 1);
	  gmm::vecsave(datafilename + s1 + "_0.U", U);
	  gmm::vecsave(datafilename + s1 + "_0.UN", UN);
	  gmm::vecsave(datafilename + s1 + "_0.UT", UT);
	  gmm::vecsave(datafilename + s1 + "_0.LN", LN1);
	  gmm::vecsave(datafilename + s1 + "_0.LT", LT1);
	  gmm::vecsave(datafilename + s1 + "_0.X", X0);
	    
	  plain_vector VM(mf_vm.nb_dof());
	  compute_Von_Mises(MS, *pELAS_nonlin, *pELAS_lin, mf_vm, VM, law_num);
	  gmm::vecsave(datafilename + s1 + "_0.VM", VM);
	}

	size_type nb_dec;
	active_functions act_fcn;
	short new_point = -1;
	iteration_corr iter_corr;
	size_type step_cont = step0_cont;
	while ((XI0 < XI_end - maxdist) && (step_cont < nb_step_cont)) {
	  cout << "beginning of step " << step + 1 << "_" << step_cont + 1
	     << ", number of variables : " << stot + 1<< endl;
	  nb_dec = 0;

	  do { /* seek a new point */
	    new_point = -1;
	    cout << "XI0 = " << XI0 << ", h = " << h / n <<  ", deltaXI = " << h * T0[stot]
		 << endl;

	    gmm::copy(X0, X); gmm::copy(T0, T);
	    gmm::add(gmm::scaled(T, h), X);

	    iter_corr = iteration_corr(residual, difference, noisy,
				       maxit_corr ? maxit_corr : 40000);
	    Newton_correction(MS, *pfinal_model, DIRICHLET, grad_XI, F2_init, deltaF2, X, T,
			      iter_corr);

	    if (iter_corr.converged()) {
	      compute_displacement(MS, *pfinal_model, *pl, *pELAS_nonlin, *pELAS_lin, U,
				   law_num);
	      gmm::mult(BN, U, UN); gmm::mult(BT, U, UT);
	      gmm::copy(FRICTION.get_LN(MS), LN1); gmm::copy(FRICTION.get_LT(MS), LT1);
	      size_type nb_change =
		compute_charact_solution(FRICTION, contact_nodes, UN, UT, UT_init, LN1, LT1,
					 CH, CH0, noisy);

	      scalar_type XI = X[stot], angle =  gmm::vect_sp(T0, T);
	      cout << "XI = " << XI << ", XI - XI0 = " << XI - XI0 << ", T0.T = " << angle;
	      if ((angle >= minangle || (nb_change == 1 && h <= h_change))
		  && (XI <= XI_end + maxdist)) {
		XI0 = XI; new_point = 1;
	      }
	      else
		cout << " - the point rejected";
	      cout << endl;
	    }

	    if (new_point <= 0) {
	      if (h > h_min) {
		h = (h_dec * h > h_min) ? h_dec * h : h_min;
		++nb_dec;
		new_point = 0;
	      } else { /* try to change the Jacobian */
		if (act_fcn.empty()) {
		  cout << "classical continuation has broken down, "
		       << "starting searching for a new Jacobian" << endl;
		  char s1[100]; sprintf(s1, "step%d", step + 1);
		  char s2[100]; sprintf(s2, "%d", step_cont);
		  gmm::vecload(datafilename + s1 + "_" + s2 + ".UN", UN);
		  gmm::vecload(datafilename + s1 + "_" + s2 + ".UT", UT);
		  gmm::vecload(datafilename + s1 + "_" + s2 + ".LN", LN1);
		  gmm::vecload(datafilename + s1 + "_" + s2 + ".LT", LT1);
		  act_fcn.compute_Jacobians(FRICTION, contact_nodes, UN, UT, UT_init, LN1,
					    LT1, CH0, limit, noisy);
		}
		if (act_fcn.compute_new_tangent(MS, *pfinal_model, DIRICHLET, FRICTION,
						grad_XI, F2_init, deltaF2, X0, T0, noisy)) {
		  gmm::copy(act_fcn.get_tangent(), T0);

		  cout << "determining direction of the new tangent" << endl;
		  gmm::add(gmm::sub_vector(T0, gmm::sub_interval(0, stot)),
			   gmm::sub_vector(X0, gmm::sub_interval(0, stot)), MS.state());
		  compute_displacement(MS, *pfinal_model, *pl, *pELAS_nonlin, *pELAS_lin, U,
				       law_num);
		  gmm::mult(BN, U, UN); gmm::mult(BT, U, UT);
		  gmm::copy(FRICTION.get_LN(MS), LN1); gmm::copy(FRICTION.get_LT(MS), LT1);
		  act_fcn.tangent_orientation(FRICTION, UN, UT, UT_init, LN1, LT1, T0);

		  h = h_init; nb_dec = 0;
		  new_point = 0;
		}
	      }
	    }
	  } while (new_point == 0);

	  if (new_point == 1) {

	    char s1[100]; sprintf(s1, "step%d", step + 1);
	    char s2[100]; sprintf(s2, "%d", step_cont+1);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".U", U);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".UN", UN);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".UT", UT);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".LN", LN1);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".LT", LT1);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".X", X);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".T", T);
	    
	    plain_vector VM(mf_vm.nb_dof());
	    compute_Von_Mises(MS, *pELAS_nonlin, *pELAS_lin, mf_vm, VM, law_num);
	    gmm::vecsave(datafilename + s1 + "_" + s2 + ".VM", VM);

	    gmm::copy(X, X0); gmm::copy(T, T0); gmm::copy(CH, CH0);
	    act_fcn.clear();
	    if ((nb_dec == 0) && (iter_corr.get_iteration() <= thr_corr))
	      h = (h_inc * h < h_max) ? h_inc * h : h_max; 
	    cout << "end of Step n° " << step + 1 << "_" << step_cont + 1 << endl;

	  } else
	    break;

	  ++step_cont;
	} // the main loop of continuation

	if (gmm::abs(XI_end - XI0) <= maxdist) {
	  X0_given = false; step0_cont = 0; 
	  convergence = 0;
	  cout << "stop continuing, restarting the standard solver"
	       << " with a new initial approximation" << endl;
	}
      } // continuation

    } while (convergence == 0);
    
    if (convergence == 1) { /* solution in the actual step has been found */

      compute_displacement(MS, *pfinal_model, *pl, *pELAS_nonlin, *pELAS_lin, U, law_num);
      if (noisy > 1) {
	gmm::mult(BN, U, UN); gmm::mult(BT, U, UT); gmm::mult(BT, U0, UT0);
	gmm::copy(FRICTION.get_LN(MS), LN1); gmm::copy(FRICTION.get_LT(MS), LT1);
	write_charact_solution(FRICTION, contact_nodes, UN, UT, UT0, LN1, LT1);
      }
      
      // gmm::copy(FRICTION.get_LN(MS), LN1);
      // gmm::copy(FRICTION.get_LT(MS), LT1);
      
      //     {
      //       gmm::iteration itercg(1e-12, 1);
      //       plain_vector LLN(nbc), LLT(nbc*(N-1));
      //       gmm::cg(MMBN, LLN, LN1, gmm::identity_matrix(), itercg);
      //       itercg.init();
      //       gmm::cg(MMBT, LLT, LT1, gmm::identity_matrix(), itercg);
      //     }
      /*    	cout << "LT1 = " << LT1 << endl;
	    cout << "LN1 = " << LN1 << endl;
      */
      
      
      //    plain_vector LT2(nbc);
      // 	maxT[step]=-LT1[0];
      // 	maxN[step]=-LN1[0];
      // Calculating forces on contact nodes and total contact pressure
      //    for (int y = 0; y < nbc; ++y) {
      //     
      //     if (N>2) 
      //         LT2[y] = sqrt(LT1[N*y]*LT1[N*y] + LT1[N*y +1]*LT1[N*y +1]);
      // 	
      // 	
      //     if (N<3) 
      //     	LT2[y]=LT1[y];
      // 	
      // // 	if (LT2[y]> maxT[step])
      // // 		maxT[step]=-LT2[y];
      // // 	
      // // 	if (LN1[y]> maxN[step])
      // // 		maxN[step]=-LN1[y];
      // 	
      // 	
      //     	SumT[step]+= -LT2[y];
      // 	SumN[step]+= -LN1[y]; 
      // 	Pressure[step]+= SumN[step]/(LX*LY);
      //     }
      //     cout << "Tu vas marcher BORDEL 3 \n";
      // // 	max[step]=sqrt(maxN[step]*maxN[step]+maxT[step]*maxT[step]);
      //    
      //     
      //     cout << "\n  \n \n normal contact pressure " << Pressure << endl;
      //     cout << "\n Total force along normal " << SumN << endl;
      //     cout << "\n Total force along tangential plane " << SumT << endl;
      // 
      // 	gmm::vecsave(datafilename + ".SumT",SumT);
      // 	gmm::vecsave(datafilename + ".SumN",SumN);
      // // 	cout << "LT2 = " << LT2 << endl;
      // 	cout << "LN2 = " << LN1 << endl;
      
      //    cout << "CN = " << FRICTION.get_LN(MS) << endl;

      //   plain_vector UN(gmm::mat_nrows(BN));
      //   gmm::mult(BN, U, UN);
      //    cout << "UN = " << UN << endl;
      
      char s[100]; sprintf(s, "step%d", step+1);
      gmm::vecsave(datafilename + s + ".F2", F2);
      gmm::vecsave(datafilename + s + ".U", U);
      gmm::vecsave(datafilename + s + ".Y", MS.state());
      
      
      plain_vector VM(mf_vm.nb_dof());
      compute_Von_Mises(MS, *pELAS_nonlin, *pELAS_lin, mf_vm, VM, law_num);
      gmm::vecsave(datafilename + s + ".VM", VM);
      
      if (rho > 0) {
	gmm::add(U, gmm::scaled(U0, -1.0), V0);
	gmm::scale(V0, 1./deltat);
	gmm::vecsave(datafilename + s + ".V", V0);
      }

      gmm::copy(U, U0);
      cout << "end of Step n° " << step+1 << " / " << nb_step << endl;
    
    // 	if (max[step]> Fmax)
    // 		Fmax = max[step];

    } else /* solution has not been found */
      break;
    
  } // the main cycle
  
//  cout << "end of all steps \n" ;
    if (law_num == 5 || law_num == 3 || law_num == 1) delete pINCOMP;
 /*  plain_vector VM(mf_u.nb_dof());
     cout << "calcul van mises\n";
   calcul_von_mises(mf_u, U0, mf_vm, VM, mu);
     cout << "Fin calcul von mises\n";
 */    

  return (convergence == 1);
  
}

  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

   // try {    
    friction_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    p.mf_u.write_to_file(p.datafilename + ".mf", true);
    p.mf_rhs.write_to_file(p.datafilename + ".mfd", true);
    p.mf_vm.write_to_file(p.datafilename + ".mfvm", true);
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) cerr << "Solve has failed\n";
    if (p.PARAM.int_value("VTK_EXPORT")) {
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "elastostatic_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << ".vtk -f ExtractVectorNorm -f "
	"WarpVector -m BandedSurfaceMap -m Outline &\n";
    }
  // }   GMM_STANDARD_CATCH_ERROR;

  return 0;
  
}
