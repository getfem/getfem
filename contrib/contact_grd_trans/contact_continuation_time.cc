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
   @file nonlinear_elastostatic.cc
   @brief Nonlinear Elastostatic problem (large strain).

   A rubber bar is submitted to a large torsion. An attempt to compute
   the solution by numerical continuation w.r.t. time. Not too
   satisfactory because the structure of retrograde branches seems to
   be very complicated.

   NOTE: In the case of neg_deltat = 0 the absolut value of deltat is
   set in the friction condition, in the opposite case deltat is used
   instead; if fact, neg_deltat = 0 should correspond to the case when
   the term deltat does not occur in the denominator of the friction
   condition whereas neg_deltat = 1 should correspond to the case when
   it occurs there; however, in both cases the derivative with respect
   to time is zero in the rows corresponding to the friction
   conditions and both cases are actually some mixture of the
   possibilities described above.
   
   This program is used to check that getfem++ is working. This is
   also a good example of use of Getfem++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_nonlinear_elasticity.h"
#include "my_getfem_Coulomb_friction.h"
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
  structure for the elastostatic problem
*/
struct elastostatic_problem {

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
  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  elastostatic_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh), mf_rhs(mesh), mf_coef(mesh), mf_vm(mesh) {}
};


/* Read parameters from the .param file, build the mesh, set finite element
   and integration methods and selects the boundaries.

   (this is boilerplate code, not very interesting)
 */
void elastostatic_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name for the pressure");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */
  int nb_refine = PARAM.int_value("NBREFINE");
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Number of space steps "));
  nsubdiv[1] = PARAM.int_value("NY") ? PARAM.int_value("NY") : nsubdiv[0];
  if (N>2) nsubdiv[2] = PARAM.int_value("NZ") ? PARAM.int_value("NZ") : nsubdiv[0];
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);

  LX = PARAM.real_value("LX", "Length along X axis");
  LY = PARAM.real_value("LY", "Length along Y axis");
  LZ = PARAM.real_value("LZ", "Length along Z axis");
  lambda = PARAM.real_value("LAMBDA", "Lame coefficient lambda");
  mu = PARAM.real_value("MU", "Lame coefficient mu");
  bgeot::base_matrix M(N,N);
  for (size_type i=0; i < N; ++i) {
    static const char *t[] = {"LX","LY","LZ"};
    M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
  }
  if (N>1) { M(0,1) = PARAM.real_value("INCLINE") * PARAM.real_value("LY"); }
  
  /* scale the unit mesh to [LX,LY,..] and incline it */
  mesh.transformation(M);

  if (nb_refine > 0) { /* refinement in the right-bottom corner */
    scalar_type layerx, layerz;
    dal::bit_vector cvref;
    layerx = 0.201; layerz = 1.5 * layerx;
    for (int i=0; i < nb_refine; ++i) {
      cvref.clear();
      for (dal::bv_visitor j(mesh.convex_index()); !j.finished(); ++j) {
	if ((mesh.points_of_convex(j)[0][0] > LX*(1-layerx))
	    && (mesh.points_of_convex(j)[1][0] > LX*(1-layerx))
	    && (mesh.points_of_convex(j)[2][0] > LX*(1-layerx))
	    && (mesh.points_of_convex(j)[0][N-1] < LZ*layerz)
	    &&  (mesh.points_of_convex(j)[1][N-1] < LZ*layerz)
	    &&  (mesh.points_of_convex(j)[2][N-1] < LZ*layerz))
	  cvref.add(j);
      }
      mesh.Bank_refine(cvref);
      layerx = 0.5 * layerx;
      if (i < 9) layerz = 0.5 * layerz;
    }
  }

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;

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

  /* set the finite element on mf_rhs (same as mf_u if DATA_FEM_TYPE is
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
   * since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));
  
  mf_vm.set_classical_discontinuous_finite_element(1);

  /* set boundary conditions
   * (Dirichlet on the upper face, contact on the bottom face 
   *  and on parts of the adjacent faces, Neumann elsewhere) */
  cout << "Selecting Dirichlet, Neumann and contact boundaries\n";
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
    } else if (mesh.points_of_convex(it.cv())[0][N-1] < LZ*0.3) {
      mesh.region(FRICTION_BOUNDARY_NUM).add(it.cv(),it.f());
    }
  }
}

/* The Iteration object for corrections; to accelerate the
   computation, heuristics is used if anglemin is used */
class iteration_corr {
  size_type maxiter; /* Max. number of iterations.                       */
  int noise;         /* if noise > 0 iterations are printed.             */
  double resmax;     /* maximum residual.                                */
  double diffmax;    /* maximum difference.                              */
  double anglemin;   /* minimum angle.                                   */
  size_type nit;     /* iteration number.                                */
  double res;        /* last computed residual.                          */
  double diff;       /* last computed difference.                        */
  double angle;      /* last computed angle.                             */
public:
  void init(void) { 
    nit = 0; res = 0.0; diff = 0.0; angle = 0.0;
    }

  iteration_corr(double r = 1.0E-8, double d = 1.0E-8, int noi = 0,
		 size_type mit = size_type(-1), double a = - 2.0)
    : maxiter(mit), noise(noi), resmax(r), diffmax(d), anglemin(a) { init(); }

  void  operator ++(int) {  nit++; }

  int get_noisy(void) const { return noise; }
  void set_res(double r) { res = r; }
  void set_diff(double d) { diff = d; }
  void set_angle(double a) { angle = a; }
  size_type get_iteration(void) const { return nit; }
  void set_iteration(size_type i) { nit = i; }

  bool converged(void) { return ((res <= resmax) && (diff <= diffmax)); }

  bool finished(void) {
    if (noise > 0) {
      cout << "iter " << nit << " residual " << gmm::abs(res)
	   << " difference " << gmm::abs(diff);
      if (anglemin >= -1.0)
	cout << " angle = " << angle;
      cout << endl;
    }
    return ((nit >= maxiter) || converged() || ((nit == 1) && (angle < anglemin)));
  }
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
  plain_vector dX(N+1), T0 = T;

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
  iter.set_angle(gmm::vect_sp(T, T0));

}

template<typename MODEL_STATE>
void Newton_correction
(MODEL_STATE &MS, getfem::mdbrick_Dirichlet<MODEL_STATE> &problem,
 getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const plain_vector &grad_t,
 const plain_vector &F2_init, const plain_vector &deltaF2, plain_vector &X, plain_vector &T,
 scalar_type t0, bool neg_deltat, iteration_corr &iter) {

  size_type stot = gmm::mat_ncols(MS.tangent_matrix());
  scalar_type deltat = X[stot] - t0, t;
  plain_vector F2 = F2_init;

  gmm::copy(gmm::sub_vector(X, gmm::sub_interval(0, stot)), MS.state());
  t = X[stot];

  deltat = t - t0; if ((deltat < 0) && (!neg_deltat)) deltat = -deltat;
  FRICTION.set_beta(1./deltat); 
  // (deltat0 < 0) ? FRICTION.set_beta(-1.0) : FRICTION.set_beta(1.0);
    
  gmm::add(gmm::scaled(deltaF2, t), F2);
  problem.rhs().set(F2);

  problem.compute_residual(MS);

  do {
    if (iter.get_noisy() > 1)
      cout << "starting computing tangent matrix" << endl;
    problem.compute_tangent_matrix(MS);

    if (iter.get_noisy() > 1)
      cout << "starting linear solver" << endl;
    correction_solver(MS.tangent_matrix(), grad_t, X, T, MS.residual(), iter);
    if (iter.get_noisy() > 1) 
      cout << "linear solver done" << endl;

    gmm::copy(gmm::sub_vector(X, gmm::sub_interval(0, stot)), MS.state());
    t = X[stot];
    
    deltat = t - t0; if ((deltat < 0) && (!neg_deltat)) deltat = -deltat;
    FRICTION.set_beta(1./deltat); 
    // (deltat0 < 0) ? FRICTION.set_beta(-1.0) : FRICTION.set_beta(1.0);
    
    gmm::copy(F2_init, F2);
    gmm::add(gmm::scaled(deltaF2, t), F2);
    problem.rhs().set(F2);

    problem.compute_residual(MS);
    iter.set_res(gmm::vect_norm2(MS.residual()));

    iter++;
  } while (!iter.finished());
}

template<typename MODEL_STATE>
void compute_displacement
(MODEL_STATE &MS, getfem::mdbrick_Dirichlet<MODEL_STATE> &problem,
 getfem::abstract_hyperelastic_law *pl,
 getfem::mdbrick_nonlinear_elasticity<MODEL_STATE> *pELAS_nonlin,
 getfem::mdbrick_isotropic_linearized_elasticity<MODEL_STATE> *pELAS_lin, plain_vector &U,
 size_type law_num) {

  if (law_num < 4) {
    pl->reset_unvalid_flag();
    problem.compute_residual(MS);
    if (pl->get_unvalid_flag()) 
      GMM_WARNING1("The solution is not completely valid, the determinant "
		   "of the transformation is negative on "
		   << pl->get_unvalid_flag() << " Gauss points");
    
    gmm::copy(pELAS_nonlin->get_solution(MS), U);
  } else
    gmm::copy(pELAS_lin->get_solution(MS), U);
}

template<typename MODEL_STATE>
void compute_Von_Mises
(MODEL_STATE &MS, getfem::mdbrick_nonlinear_elasticity<MODEL_STATE> *pELAS_nonlin,
 getfem::mdbrick_isotropic_linearized_elasticity<MODEL_STATE> *pELAS_lin,
 getfem::mesh_fem &mf_vm, plain_vector &VM, size_type law_num) {
  /* computes the Von Mises stress from the solution */
  if (law_num < 4)
    pELAS_nonlin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
  else 
    pELAS_lin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
}

template<typename MODEL_STATE, typename MATRIX, typename VECT>
void test_functions
(MODEL_STATE &MS, getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION,
 const VECT &contact_nodes, const sparse_matrix &BN, const sparse_matrix &BT,
 const plain_vector &U, const plain_vector &U0, MATRIX &TST, scalar_type classboundary,
 int noisy) {
  /* saves the values of the test functions at the solution given by
     MS.state() and U into TST */

  size_type nb_dof_N = gmm::mat_nrows(BN), nb_dof_T = gmm::mat_nrows(BT);
  scalar_type r = FRICTION.get_r(), alpha = FRICTION.get_alpha(), beta = FRICTION.get_beta();
  plain_vector gap(nb_dof_N), friction_coef(nb_dof_N);
  plain_vector UN(nb_dof_N), UT(nb_dof_T), UT0(nb_dof_T), LN(nb_dof_N), LT(nb_dof_T);
    
  gmm::copy(FRICTION.get_gap(), gap); gmm::copy(FRICTION.get_friction_coef(), friction_coef);
  gmm::mult(BN, U, UN); gmm::mult(BT, U, UT); gmm::mult(BT, U0, UT0);
  gmm::copy(FRICTION.get_LN(MS), LN); gmm::copy(FRICTION.get_LT(MS), LT);

  if (noisy > 1)
    cout << "characters, eventually limit values of the test functions:" << endl;  

  for (size_type i = 0; i < nb_dof_N; ++i) {
    TST(i, 0) = (LN[i] - r * alpha * (UN[i] - gap[i]));
    TST(i, 1) = - friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);
    TST(i, 2) = friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);

    if (noisy > 1) {
      cout << "node " << i << ": " << contact_nodes[i] << " "
	   << (TST(i, 0) <= 0) << (TST(i, 1) >= 0) << (TST(i, 2) <= 0) << " ";

      for (size_type j = 0; j < 3; ++j){
	if ((TST(i,j) >= -classboundary && TST(i,j) <= classboundary)
		 && ((j == 0) || (j > 0 && TST(i,0) <= classboundary)))
	  cout << " ("<< j+1 << "-th component close to the limit: " << TST(i,j) << ") ";
      }

      cout << endl;
      if (noisy > 2)
	cout << "LN[i] = " << LN[i] 
	     << ", r * alpha * (UN[i] - gap[i]) = " << r * alpha * (UN[i] - gap[i])
	     << ", - friction_coef[i] * LN[i] + LT[i] = "
	     << -friction_coef[i] * LN[i] + LT[i]
	     << ", friction_coef[i] * LN[i] + LT[i] = "
	     << friction_coef[i] * LN[i] + LT[i] 
	     << ", r * beta * (UT[i]-UT0[i]) = "
	     << r * beta * (UT[i]-UT0[i]) << endl;
    }
  }
}

template<typename MODEL_STATE, typename MATRIX, typename VECT>
size_type test_functions
(MODEL_STATE &MS, getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION,
 const VECT &contact_nodes, const sparse_matrix &BN, const sparse_matrix &BT,
 const plain_vector &U, const plain_vector &U0, MATRIX &TST, MATRIX &TST0,
 scalar_type classboundary, int noisy) {
  /* saves the values of the test functions at the solution given by
     MS.state() and U into TST; moreover, returns the number of
     changes in comparison with TST0 */

  size_type nb_dof_N = gmm::mat_nrows(BN), nb_dof_T = gmm::mat_nrows(BT);
  size_type nb_change = 0;
  scalar_type r = FRICTION.get_r(), alpha = FRICTION.get_alpha(), beta = FRICTION.get_beta();
  plain_vector gap(nb_dof_N), friction_coef(nb_dof_N);
  plain_vector UN(nb_dof_N), UT(nb_dof_T), UT0(nb_dof_T), LN(nb_dof_N), LT(nb_dof_T);
    
  gmm::copy(FRICTION.get_gap(), gap); gmm::copy(FRICTION.get_friction_coef(), friction_coef);
  gmm::mult(BN, U, UN); gmm::mult(BT, U, UT); gmm::mult(BT, U0, UT0);
  gmm::copy(FRICTION.get_LN(MS), LN); gmm::copy(FRICTION.get_LT(MS), LT);

  if (noisy > 1)
    cout << "changes of signs or limit values of the test functions (if any):" << endl;  

  for (size_type i = 0; i < nb_dof_N; ++i) {
    TST(i, 0) = (LN[i] - r * alpha * (UN[i] - gap[i]));
    TST(i, 1) = - friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);
    TST(i, 2) = friction_coef[i] * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);

    for (size_type j = 0; j < 3; ++j) {
      if (((TST0(i,j) < 0 && TST(i,j) > 0) || (TST0(i,j) > 0 && TST(i,j) < 0))
	  && (j == 0 || (j > 0 && TST(i,0) < 0))) {
	++nb_change;
	if (noisy > 1)
	  cout << "node " << i << ": " << contact_nodes[i]
	       << (TST(i, 0) <= 0) << (TST(i, 1) >= 0) << (TST(i, 2) <= 0)
	       << " - change of the sign of the " << j+1 << "-th component: "
	       << TST0(i,j) << " -> " << TST(i,j) << endl;
      } else if ((TST(i,j) >= -classboundary && TST(i,j) <= classboundary)
		 && ((j == 0) || (j > 0 && TST(i,0) <= classboundary)))
	cout << "node " << i << ": " << contact_nodes[i] 
	     << j+1 << "-th component close to the limit: " << TST(i,j) << endl;

      if (noisy > 2)
	cout << "LN[i] = " << LN[i] 
	     << ", r * alpha * (UN[i] - gap[i]) = " << r * alpha * (UN[i] - gap[i])
	     << ", - friction_coef[i] * LN[i] + LT[i] = "
	     << -friction_coef[i] * LN[i] + LT[i]
	     << ", friction_coef[i] * LN[i] + LT[i] = "
	     << friction_coef[i] * LN[i] + LT[i] 
	     << ", r * beta * (UT[i]-UT0[i]) = "
	     << r * beta * (UT[i]-UT0[i]) << endl;
    }
  }
  
  return nb_change;
}

template<typename MATRIX>
void straight_insertion
(MATRIX &ACT, plain_vector &KEY, size_type i){
/* places the i-th row of ACT according to the corresponding values in KEY;
   increasing order is wanted in the first i components of KEY */

  bool found;
  size_type j = i;
  scalar_type X_key;
  std::vector<size_type> X(3);

  X_key = KEY[i];
  X[0] = ACT(i, 0);
  X[1] = ACT(i, 1);
  X[2] = ACT(i, 2);
  if (j == 0) found = true;
  else found = (X_key >= KEY[j - 1]);
 
  while (!found) { // seeking the appropriate  position
    KEY[j] = KEY[j - 1];
    ACT(j, 0) = ACT(j - 1, 0);
    ACT(j, 1) = ACT(j - 1, 1);
    ACT(j, 2) = ACT(j - 1, 2);
    --j;
    if (j == 0) found = true;
    else found = (X_key >= KEY[j - 1]);
  }

  KEY[j] = X_key;
  ACT(j, 0) = X[0];
  ACT(j, 1) = X[1];
  ACT(j, 2) = X[2];
}
	    
template<typename MATRIX, typename T_MATRIX, typename VECT>
void compute_Jacobians
(const VECT &contact_nodes, MATRIX &JAC, MATRIX &ACT, const T_MATRIX &TST,
 scalar_type classboundary) {
  /* computes the arrays JAC and ACT determining the Jacobian and the
     active selection functions, respectively; the items in ACT are
     ordered so that the corresponding values in KEY are increasing (a
     heuristics is used!), the values of the test functions are saved
     in TST */
  
  size_type nbc = gmm::mat_nrows(TST);
  plain_vector KEY(3*nbc);
  
  size_type pos = 0;
  for (size_type i = 0; i < nbc; ++i) {
    JAC(i, 0) = i;

    JAC(i, 1) = (TST(i, 0) <= 0);
    if ((TST(i, 0) > -classboundary) && (TST(i, 0) < classboundary)) {
      ACT(pos, 0) = i;
      ACT(pos, 1) = 1;
      ACT(pos, 2) = JAC(i, 1);
      KEY[pos] = gmm::abs(TST(i, 0));
      straight_insertion(ACT, KEY, pos);
      ++pos;
    }

    if (TST(i, 1) < 0) {                  // positive slip
      JAC(i, 2) = 0; JAC(i, 3) = 1;
      if (TST(i, 0) < classboundary) {    // contact is possible
	ACT(pos, 0) = i;
	ACT(pos, 1) = 2;
	ACT(pos, 2) = JAC(i, 2);
	KEY[pos] = (-TST(i, 1) < 1./contact_nodes[i][0]) ?
	  -TST(i, 1) : 1./contact_nodes[i][0];                  // (!)
	straight_insertion(ACT, KEY, pos);
	++pos;
      }
    }
    else if (TST(i,2) > 0) {              // negative slip
      JAC(i, 2) = 1; JAC(i, 3) = 0;
      if (TST(i, 0) < classboundary) {    // contact is possible
	ACT(pos, 0) = i;
	ACT(pos, 1) = 3;
	ACT(pos, 2) = JAC(i, 3);
	KEY[pos] = (TST(i, 2) < 1./contact_nodes[i][0]) ?
	  TST(i, 2) : 1./contact_nodes[i][0];                  // (!)
	straight_insertion(ACT, KEY, pos);
	++pos;
      }
    }
    else {                                 // stick
      JAC(i, 2) = 1; JAC(i, 3) = 1;
      if (TST(i, 0) < classboundary) {     // contact is possible
	if (TST(i, 1) < classboundary) {
	  ACT(pos, 0) = i;
	  ACT(pos, 1) = 2;
	  ACT(pos, 2) = JAC(i, 2);
	  KEY[pos] = TST(i, 1);                                // (!)
	  straight_insertion(ACT, KEY, pos);
	  ++pos;
	}

	if (TST(i, 2) > -classboundary) {
	  ACT(pos, 0) = i;
	  ACT(pos, 1) = 3;
	  ACT(pos, 2) = JAC(i, 3);
	  KEY[pos] = -TST(i, 2);                               // (!)
	  straight_insertion(ACT, KEY, pos);
	  ++pos;
	}
      }
    }
  }

  JAC(nbc, 0) = size_type(-1);
  ACT(pos, 0) = size_type(-1);
}

template<typename MODEL_STATE, typename MATRIX, typename VECT>
bool test_Jacobian
(MODEL_STATE &MS, getfem::mdbrick_Dirichlet<MODEL_STATE> &problem,
 getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const VECT &contact_nodes,
 const plain_vector &grad_t, const plain_vector &F2_init, const plain_vector &deltaF2,
 const plain_vector &U00, const plain_vector &U0, scalar_type deltat0, plain_vector &T0,
 const plain_vector &X0, const MATRIX &JAC, scalar_type h, scalar_type minangle_back,
 bool neg_deltat, iteration_corr &iter) {
  /* tries to continuate according to the Jacobian given by JAC;
     returns `true' in the case of success, `false' otherwise */

  bool success = false;
  size_type stot = gmm::mat_ncols(MS.tangent_matrix());
  plain_vector F2 = F2_init, T_test = T0, T(stot+1), X = X0;

  if (iter.get_noisy() > 1)
    cout << "testing the following Jacobian:" << endl;

  for (size_type i = 0; i < gmm::mat_nrows(JAC) - 1; ++i) {
    if (iter.get_noisy() > 1)
      cout << "node " << JAC(i, 0) << ": " << contact_nodes[i] << " "
	   << JAC(i, 1) << JAC(i, 2) << JAC(i, 3) << endl;
    if (((JAC(i, 2) == 0) && (JAC(i, 3) == 0))
	|| ((JAC(i,1) == 0) && ((JAC(i, 2) == 1) && (JAC(i, 3) == 1))))  /* not wanted */
      return success; 
    }
	   
  gmm::copy(gmm::sub_vector(X0, gmm::sub_interval(0, stot)), MS.state());
  
  scalar_type t0 = X0[stot];
  gmm::add(gmm::scaled(deltaF2, t0), F2);
  problem.rhs().set(F2);

  FRICTION.set_WT(gmm::scaled(U00, -1.0));
  FRICTION.set_beta(1./deltat0);
  
  problem.compute_tangent_matrix(MS);

  FRICTION.set_JAC(JAC);
  FRICTION.do_compute_tangent_matrix(MS, 0, 0);
  compute_tangent(MS.tangent_matrix(), grad_t, T_test, iter.get_noisy());
  
  FRICTION.clear_JAC();
  FRICTION.set_WT(gmm::scaled(U0, -1.0));
  gmm::copy(T_test, T);

  gmm::add(gmm::scaled(T, h), X);
  iter.set_iteration(0);
  Newton_correction(MS, problem, FRICTION, grad_t, F2_init, deltaF2, X, T, t0, neg_deltat,
		    iter);
  
  if (iter.converged()) {
    scalar_type angle = gmm::vect_sp(T, T0);
    if (iter.get_noisy()) cout << "angle = " << angle << endl;

    if (angle >= minangle_back) { // tests whether we have not arrived at the incoming branch
      cout << ", success? (1/0): ";  // for confirmation
      cin >> success;
    }
  }

  if (!success) { // try the opposite tangent
    gmm::copy(X0, X);
    gmm::scale(T_test, -1.0);
    gmm::copy(T_test, T);
    
    gmm::add(gmm::scaled(T, h), X);
    iter.set_iteration(0);
    Newton_correction(MS, problem, FRICTION, grad_t, F2_init, deltaF2, X, T, t0, neg_deltat,
		      iter);
  
    if (iter.converged()) {
      scalar_type angle = gmm::vect_sp(T, T0);
      if (iter.get_noisy()) cout << "angle = " << angle;

      if (angle >= minangle_back) { // tests whether we have not
				    // arrived at the incoming branch
	cout << ", success? (1/0): ";  // for confirmation
	cin >> success;
      }
    }
  }

  if (success) gmm::copy(T_test, T0);
  return success;
}

template<typename MODEL_STATE, typename MATRIX, typename VECT>
bool systematic_test_Jacobians
(MODEL_STATE &MS, getfem::mdbrick_Dirichlet<MODEL_STATE> &problem,
 getfem::mdbrick_Coulomb_friction<MODEL_STATE> &FRICTION, const VECT &contact_nodes,
 const plain_vector &grad_t, const plain_vector &F2_init, const plain_vector &deltaF2,
 const plain_vector &U00, const plain_vector &U0, scalar_type deltat0, plain_vector &T0,
 const plain_vector &X0, MATRIX &JAC, const MATRIX &ACT, size_type level, scalar_type h,
 scalar_type minangle_back, bool neg_deltat, iteration_corr &iter) {
  /* a recursive procedure for systematic testing of possible Jacobians, starts with 
     the changes corresponding to the test functions with the closest values to zero;
     returns true in the case success, false otherwise */

  bool success = false;

  if (level == 0)
    success =
      test_Jacobian(MS, problem, FRICTION, contact_nodes, grad_t, F2_init, deltaF2, U00, U0,
		    deltat0, T0, X0, JAC, h, minangle_back, neg_deltat, iter);
  else
    success = 
      systematic_test_Jacobians(MS, problem, FRICTION, contact_nodes, grad_t, F2_init,
				deltaF2, U00, U0, deltat0, T0, X0, JAC, ACT, level - 1, h,
				minangle_back, neg_deltat, iter);

  if (!success){
    JAC(ACT(level, 0), ACT(level, 1)) = (ACT(level, 2) == 0);
    if (level == 0)
      success = 
	test_Jacobian(MS, problem, FRICTION, contact_nodes, grad_t, F2_init, deltaF2, U00, U0,
		      deltat0, T0, X0, JAC, h, minangle_back, neg_deltat, iter);
    else
      success =
	systematic_test_Jacobians(MS, problem, FRICTION, contact_nodes, grad_t, F2_init,
				  deltaF2, U00, U0, deltat0, T0, X0, JAC, ACT, level - 1, h,
				  minangle_back, neg_deltat, iter);
    JAC(ACT(level, 0), ACT(level, 1)) = ACT(level, 2);
  }

  return success;
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool elastostatic_problem::solve(plain_vector &U) {
  
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
  bool neg_deltat = (PARAM.int_value("NEGATIVE_DELTAT",
				     "Negative deltat is permitted or not") != 0);
  scalar_type deltat = PARAM.real_value("DELTAT", "Time step");
  if ((deltat < 0) && (!neg_deltat)) deltat = -deltat;
 
  // contact condition for Lagrange elements
  dal::bit_vector cn = mf_u.basic_dof_on_region(FRICTION_BOUNDARY_NUM);
  sparse_matrix BN(cn.card()/N, mf_u.nb_dof());
  sparse_matrix BT((N-1)*cn.card()/N, mf_u.nb_dof());
  std::vector<base_node> contact_nodes;
  plain_vector gap(cn.card()/N);
  size_type jj = 0;
  for (dal::bv_visitor i(cn); !i.finished(); ++i)
    if (i % N == 0) {
      BN(jj, i+N-1) = -1.;
      gap[jj] = mf_u.point_of_basic_dof(i)[N-1];
      contact_nodes.push_back(mf_u.point_of_basic_dof(i));
      for (size_type k = 0; k < N-1; ++k) BT((N-1)*jj+k, i+k) = 1.;
      ++jj;
    }
  
  // creating force density vectors
  size_type nbc = int(jj);
  sparse_matrix MMBN(nbc, nbc), MMBT(nbc*(N-1), nbc*(N-1));
//   plain_vector LN1(nbc), LT1(nbc*(N-1));
  {
    sparse_matrix BB(mf_u.nb_dof(), mf_u.nb_dof());
    getfem::asm_mass_matrix(BB, mim, mf_u, mf_u, FRICTION_BOUNDARY_NUM);
    std::vector<size_type> indN, indT;
    for (dal::bv_visitor i(cn); !i.finished(); ++i)
      if ((i%N) == N-1) indN.push_back(i); else indT.push_back(i);
    gmm::sub_index SUBI(indN);
    gmm::copy(gmm::sub_matrix(BB, SUBI, SUBI), MMBN);
    gmm::sub_index SUBJ(indT);
    gmm::copy(gmm::sub_matrix(BB, SUBJ, SUBJ), MMBT);    
  }



  scalar_type friction_coef = PARAM.real_value("FRICTION_COEFF",
					       "Friction cefficient");
  scalar_type r = PARAM.real_value("R", "Augmentation parameter");
  scalar_type alpha = PARAM.real_value("ALPHA", "Augmentation parameter");
  gmm::dense_matrix<size_type> JAC(nbc+1, 4);
  JAC(0, 0) = size_type(-1);
  

  getfem::mdbrick_Coulomb_friction<> FRICTION(*pINCOMP, BN, gap,
					      friction_coef, BT);
  FRICTION.set_r(r);
  FRICTION.set_alpha(alpha);
  FRICTION.set_beta(1./deltat);
  FRICTION.set_JAC(JAC);
 
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
  getfem::mdbrick_Dirichlet<> final_model(VOL_F, DIRICHLET_BOUNDARY_NUM);
  final_model.rhs().set(mf_rhs, F2);
  final_model.set_constraints_type(getfem::constraints_type
				   (PARAM.int_value("DIRICHLET_VERSION")));

  // Generic solver.
  getfem::standard_model_state MS(final_model);
  size_type step0 = PARAM.int_value("STEP0") ? PARAM.int_value("STEP0") : 0;
  size_type maxit = PARAM.int_value("MAXITER");
  gmm::iteration iter;

  scalar_type dz = PARAM.real_value("DIRICHLET_Z",
				    "Prescribed displacement in z");
  scalar_type dyv = PARAM.real_value("DIRICHLET_Y_SPEED",
				     "Prescribed velocity in y"); 
  int noisy = PARAM.int_value("NOISY");

  plain_vector U0 = U;
  size_type stot = gmm::mat_ncols(MS.tangent_matrix());
  plain_vector X0(stot+1);

  if (step0 > 0) { /* load the foregoing values */
    char s[100]; sprintf(s, "step%d", step0);
    gmm::vecload(datafilename + s + ".U", U0);
    gmm::vecload(datafilename + s + ".X", X0);
    gmm::copy(gmm::sub_vector(X0, gmm::sub_interval(0, stot)), MS.state());
  }

  /* the standard solver is used for the first eleven steps when the body
     is compressed, for the next iterations, when the body is pulling
     to the right while subject to a constant pressure, numerical
     continuation is used */

  bool converged = true;
  size_type step = step0;
  while ((step < nb_step) && (step < 11) && converged) { /* standard solver */
    cout << "beginning of step " << step+1
	 << ", number of variables : " << final_model.nb_dof() << endl;

    for (size_type i = 0; i < nb_dof_rhs; ++i) {
      F2[i*N+N-1] = dz * step / 10.0;
      F2[i*N+N-2] = 0.0;
    }
    final_model.rhs().set(F2);
    
    FRICTION.set_WT(gmm::scaled(U0, -1.0));
    
    
    iter = gmm::iteration(residual, noisy, maxit ? maxit : 40000);
//     cout << "|U| = " << gmm::vect_norm2(MS.state()) << "\n";
    
    
    
    /* let the default non-linear solve (Newton) do its job */
    gmm::default_newton_line_search ls; // (size_type(-1), 5.0/3.0, 0.1, 0.5, 3.0)
    getfem::standard_solve(MS, final_model, iter,
			   getfem::default_linear_solver(final_model), ls);
    
    converged = iter.converged();

    if (converged) {
    
      compute_displacement(MS, final_model, pl, pELAS_nonlin, pELAS_lin, U, law_num);
    
      char s[100]; sprintf(s, "step%d", step+1);
      gmm::vecsave(datafilename + s + ".U",U);
      
      gmm::copy(MS.state(), gmm::sub_vector(X0, gmm::sub_interval(0, stot))); X0[stot] = 0.0;
      gmm::vecsave(datafilename + s + ".X", X0);
    
      plain_vector VM(mf_vm.nb_dof());
      compute_Von_Mises(MS, pELAS_nonlin, pELAS_lin, mf_vm, VM, law_num);
      gmm::vecsave(datafilename + s + ".VM", VM);
    
      gmm::copy(U, U0);
      cout << "end of Step n° " << step+1 << " / " << nb_step << endl;
    }
     
    ++step;
  } // standard solver
 
  if ((nb_step > 11) && converged) { /* continuation */
    GMM_ASSERT1(PARAM.int_value("DIRICHLET_VERSION") == 0,
		"Continuations only implemented for Dirichlet "
		"condition with multipliers");

    scalar_type difference = PARAM.real_value("DIFFERENCE");
    if (difference == 0.) difference = 1e-10;
    scalar_type minangle = PARAM.real_value("ANGLE");
    scalar_type minangle_back = PARAM.real_value("ANGLE_BACK");
    size_type maxit_corr = PARAM.int_value("MAXITER_CORR"); 
    size_type thr_corr = PARAM.int_value("THRESHOLD_CORR");
    scalar_type classboundary = PARAM.real_value
      ("CLASSBOUNDARY", "Parameter for classification of the boundary");
 
    scalar_type h_init = PARAM.real_value("H_INIT", "Initial step length");
    scalar_type h_max = PARAM.real_value("H_MAX", "Maximal step length");
    scalar_type h_min = PARAM.real_value("H_MIN", "Minimal step length");
    scalar_type h_inc = PARAM.real_value("H_INC");
    scalar_type h_dec = PARAM.real_value("H_DEC");
    scalar_type h = PARAM.real_value("H") ? PARAM.real_value("H") : h_init;

    scalar_type deltat0, t0, t;
    plain_vector deltaF2(nb_dof_rhs * N);
    plain_vector U00 = U0;
    plain_vector X(stot+1), T0(stot+1), T(stot+1);
    gmm::dense_matrix<scalar_type> TST(nbc,3), TST0(nbc,3);
    size_type ind =  final_model.first_ind(), sc = gmm::vect_size(final_model.get_CRHS());
    plain_vector grad_t(stot);

    if (step0 >= 11)
      for (size_type i = 0; i < nb_dof_rhs; ++i) {
	F2[i*N+N-1] = dz;
	F2[i*N+N-2] = 0.0;
      }
    for (size_type i = 0; i < nb_dof_rhs; ++i) {
      deltaF2[i*N+N-1] = 0.0;
      deltaF2[i*N+N-2] = dyv;
    }    
    final_model.rhs().set(deltaF2);  
    final_model.compute_tangent_matrix(MS);
    gmm::copy(final_model.get_CRHS(),
	      gmm::sub_vector(grad_t, gmm::sub_interval(ind, sc)));
    gmm::scale(grad_t, -1.0); 

    char s[100]; sprintf(s, "step%d", step - 1);
    plain_vector X00(stot+1);
    gmm::vecload(datafilename + s + ".X",X00);
    gmm::vecload(datafilename + s + ".U",U00); 

    sprintf(s, "step%d", step);
    gmm::vecload(datafilename + s + ".U",U0);
    t0 = X0[stot];
    deltat0 = (step > 11) ? t0 - X00[stot] : deltat;
    if ((deltat0 < 0) && (!neg_deltat)) deltat0 = -deltat0;

    gmm::copy(gmm::sub_vector(X0, gmm::sub_interval(0, stot)), MS.state());
    FRICTION.set_beta(1./deltat0);
    // (deltat0 < 0) ? FRICTION.set_beta(-1.0) : FRICTION.set_beta(1.0);
    FRICTION.set_WT(gmm::scaled(U00, -1.0));

    test_functions(MS, FRICTION, contact_nodes, BN, BT, U0, U00, TST0, classboundary, noisy);

      if (step == 11) {
	gmm::fill_random(T0);
	compute_tangent(MS.tangent_matrix(), grad_t, T0, noisy);
	if (T0[stot] < 0)
	  gmm::scale(T0, -1.0);
      }
      else
	gmm::vecload(datafilename + s + ".T",T0);


      size_type nb_dec;
      bool new_Jacobian = false;
      iteration_corr iter_corr;
      short new_point = 1;

      while ((new_point == 1) && (step < nb_step)) {
	cout << "beginning of step " << step+1
	     << ", number of variables : " << final_model.nb_dof() << endl;
	
	FRICTION.set_WT(gmm::scaled(U0, -1.0));
	nb_dec = 0;
      
// 	cout << " |U| = " << gmm::vect_norm2(MS.state()) << "\n";
	
	do { /* seek a new point */
	   new_point = -1;
	   cout << "t0 = " << t0 <<  ", h = " << h << ", deltat = " << h * T0[stot] << endl;
     
	   gmm::copy(X0, X); gmm::copy(T0, T);
	   gmm::add(gmm::scaled(T, h), X);
	   
	   iter_corr = iteration_corr(residual, difference, noisy,
				      maxit_corr ? maxit_corr : 40000);
	   Newton_correction(MS, final_model, FRICTION, grad_t, F2, deltaF2, X, T, t0,
			     neg_deltat, iter_corr);
	   
	   if (iter_corr.converged()) {
	     compute_displacement(MS, final_model, pl, pELAS_nonlin, pELAS_lin, U, law_num);
	     size_type nb_change = test_functions(MS, FRICTION, contact_nodes, BN, BT, U,
						  U0, TST, TST0, classboundary, noisy);
	     
	     t = X[stot]; scalar_type angle = gmm::vect_sp(T, T0);
	     cout << "t = " << t << ", deltat = " << t - t0 << ", T0.T = " << angle;
	     if ((angle > 0.9) || ((angle > 0.7) && (nb_change < 2)) || new_Jacobian) {
	       new_point = 1;
	       new_Jacobian = false;
	     } else
// 	     { cout << "Do you want to accept this increment? (1, 0, -1) ";
// 	      int sign;
// 	      cin >> sign;
// 	      if (sign != 0) {
// // 		if (sign < 0) gmm::scale(T, -1.0);
// // 		h = h_init;
// // 		++nb_dec;
// 		break;
// 	      }
// 	    } 
	       cout << " - the point rejected";
	     cout << endl;
	   }
	   
	   if (new_point <= 0) {
	     if (h > h_min) {
	       h = (h_dec * h > h_min) ? h_dec * h : h_min;
	       ++nb_dec;
	       new_point = 0;
	     } else { /* try to change the Jacobian manually */
	       
	       cout << "classical continuation has broken down"
		    << ", starting to change the Jacobian manually" << endl;
	       
	       new_Jacobian = false;
	       iteration_corr iter_test(residual, difference, noisy,
					maxit_corr ? maxit_corr : 40000, minangle);
	       gmm::dense_matrix<size_type> ACT(3*nbc+1, 3);
	       
	       compute_Jacobians(contact_nodes, JAC, ACT, TST0, classboundary);
	       
	       for (size_type level = 0; level < 3*nbc+1; ++level) {
		 if (ACT(level, 0) == size_type(-1)) break;
		 
		 JAC(ACT(level, 0), ACT(level, 1)) = (ACT(level, 2) == 0);
		 cout << "change for node " << ACT(level, 0) << endl;
		 
		 if (level > 0) 
		   new_Jacobian = systematic_test_Jacobians
		     (MS, final_model, FRICTION, contact_nodes, grad_t, F2, deltaF2, U00, U0,
		      deltat0, T0, X0, JAC, ACT, level - 1, h_min, minangle_back, neg_deltat,
		      iter_test);
		 else
		   new_Jacobian = test_Jacobian
		     (MS, final_model, FRICTION, contact_nodes, grad_t, F2, deltaF2, U00, U0,
		      deltat0, T0, X0, JAC, h, minangle_back, neg_deltat, iter_test);
		 
		 if (new_Jacobian)
		   break;
		 else
		   JAC(ACT(level, 0), ACT(level, 1)) = ACT(level, 2);
	       }
	       
	       if (new_Jacobian) {
		 cout << "new Jacobian found, restarting classical continuation" << endl;
		 h = h_init;
		 nb_dec = 0;
		 new_point = 0;
	       }
	     }
	   }
	} while (new_point == 0);
	  
	if (new_point == 1) {

	  gmm::copy(U0, U00); gmm::copy(U, U0); 
	  gmm::copy(X, X0); gmm::copy(T, T0);
	  gmm::copy(TST, TST0);
	  deltat0 = deltat; t0 = X[stot];
	  
	  sprintf(s, "step%d", step+1);
	  gmm::vecsave(datafilename + s + ".U",U);
	  gmm::vecsave(datafilename + s + ".X", X);
	  gmm::vecsave(datafilename + s + ".T", T);
	  
	  plain_vector VM(mf_vm.nb_dof());
	  compute_Von_Mises(MS, pELAS_nonlin, pELAS_lin, mf_vm, VM, law_num);
	  gmm::vecsave(datafilename + s + ".VM", VM);
	  
	  if ((h < h_max) && (nb_dec == 0) && (iter_corr.get_iteration() <= thr_corr))
	    h = (h*h_inc < h_max) ? h*h_inc : h_max;
	  
	  cout << "end of Step n° " << step+1 << " / " << nb_step << endl;
	  ++step;
	}
      } // the main loop of continuation

      converged = (new_point == 1);
  } // continuation
  

  if (law_num == 5 || law_num == 3 || law_num == 1) delete pINCOMP;    
  
  return converged;
  
}

  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

   // try {    
    elastostatic_problem p;
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
