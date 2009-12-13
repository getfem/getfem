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

   A rubber bar is submitted to a large torsion.
   
   This program is used to check that getfem++ is working. This is also 
   a good example of use of Getfem++.
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

template<typename MAT> void smallest_eigenvalue(const MAT &M);

template<class VEC, class VECT, class MAT>
void compute_tangent(const VEC &T, const MAT &A, VECT &V, int noisy);

template<class VEC, class VECT, class MAT>
int compute_Newton_iteration(const VEC &T, const MAT &A, VECT &V, VECT &X, const VEC &P,
			     scalar_type tolerance, scalar_type residual, int noisy);

/*
  structure for the elastostatic  problem
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
  scalar_type LX, LY, LZ;    /* system dimensions   */
  // scalar_type lambda, mu;    /* LamÈ coefficients.                           */
  scalar_type residual;        /* max residual for the iterative solvers         */
  std::string datafilename;
  bgeot::md_param PARAM;


  template<class VEC, class STDVEC, class ALLOCVEC, class MAT>
  int compute_test_functions(const VEC &UN, const VEC &UT, const VEC &UT0, const VEC &LN,
			     const VEC &LT, MAT &TST, const MAT &TST0,
			     const STDVEC &contact_nodes, const ALLOCVEC &gap,
			     scalar_type deltat, scalar_type friction_coef);
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
  int refine = PARAM.int_value("REFINE");
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
  // mu = PARAM.real_value("MU", "Lam√© coefficient mu");
  // lambda = PARAM.real_value("LAMBDA", "Lam√© coefficient lambda");
  bgeot::base_matrix M(N,N);
  for (size_type i=0; i < N; ++i) {
    static const char *t[] = {"LX","LY","LZ"};
    M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
  }
  if (N>1) { M(0,1) = PARAM.real_value("INCLINE") * PARAM.real_value("LY"); }
  
  /* scale the unit mesh to [LX,LY,..] and incline it */
  mesh.transformation(M);

  if (refine > 0) { /* manual refinement */
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
   * (Dirichlet on the upper face, contact on the bottom face) */
      cout << "Selecting Dirichlet and contact boundaries\n";
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

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

template<typename MAT>
void smallest_eigenvalue(const MAT &M) {
  
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
    if (0) // (i > 90)
      cout << "iteration " << i << " norm(W) = " <<  gmm::vect_norm2(W)
	   << " lambda = " << lambda 
	   << " lambda2 = " << lambda2 << endl;
    if (gmm::vect_norm2(W) < 1E-2) break;
    
    gmm::copy(U, V);
    gmm::copy(U2, V2);
  }
  
  // cout << "Smallest eigenvalue : " << 1.0/lambda << endl;
  // cout << "Largest eigenvalue : " << lambda2 << endl;
  cout << "Condition number: " << lambda2 * lambda << endl;
  
}


/* computes initial tangent for the Moore-Penrose method */
template<class VEC, class VECT, class MAT>
void compute_tangent(const VEC &T, const MAT &A, VECT &V, int noisy) {

  bool singular_tangent_matrix = false;
  size_type N = gmm::mat_nrows(A);
  double rcond;
  scalar_type sp;
  plain_vector Y(N);
		
  try { /* whether A is singular */
    if (noisy > 1) cout << "starting computing tangent" << endl;

    gmm::SuperLU_solve(A, Y, T, rcond);;
    sp = gmm::vect_sp(gmm::sub_vector(V, gmm::sub_interval(1,N)), Y);
		  
    V[0] = 1.0 / (V[0] - sp);
    gmm::copy(gmm::scaled(Y, -1.0*V[0]),
	      gmm::sub_vector(V, gmm::sub_interval(1, N)));
  }
  catch (...) {
    singular_tangent_matrix = true;
    cout << "Tangent matrix is close to be singular!" << endl;
  }
  
  if (singular_tangent_matrix) {
    plain_vector R(1+N);
    sparse_matrix B(1+N, 1+N);

    gmm::copy(gmm::col_vector(T),
	      gmm::sub_matrix(B, gmm::sub_interval(0, N), gmm::sub_interval(0, 1)));
    gmm::copy(A, gmm::sub_matrix(B, gmm::sub_interval(0, N), gmm::sub_interval(1, N)));
    gmm::copy(gmm::row_vector(V),
	      gmm::sub_matrix(B, gmm::sub_interval(N, 1), gmm::sub_interval(0, 1+N)));

    gmm::clear(R);
    R[N] = 1.0;		  
    
    gmm::SuperLU_solve(B, V, R, rcond);
  }

  gmm::scale(V, 1.0/gmm::vect_norm2(V));

  if (noisy > 1) cout << "tangent computed" << endl;

} // compute_tangent


/* Newton iteration for the Moore-Penrose method */
template<class VEC, class VECT, class MAT>
int compute_Newton_iteration(const VEC &T, const MAT &A, VECT &V, VECT &X, const VEC &P,
			     scalar_type tolerance, scalar_type residual, int noisy) {

  bool singular_tangent_matrix = false;
  int convergence = 1;
  size_type N = gmm::mat_nrows(A);
  scalar_type sp1, sp2, normdX, normP;
  plain_vector Y1(N), Y2(N);
  plain_vector dX(1+N);

  try { /* whether A is singular */
    if (noisy > 1) {
      cout << "starting linear solvers" << endl;
      smallest_eigenvalue(A);
    }

    gmm::SuperLU_factor<scalar_type> SLU;
    SLU.build_with(A);
    SLU.solve(Y1, T);
    SLU.solve(Y2, P);

    sp1 = gmm::vect_sp(gmm::sub_vector(V, gmm::sub_interval(1,N)), Y1);
    sp2 = gmm::vect_sp(gmm::sub_vector(V, gmm::sub_interval(1,N)), Y2);
		
    dX[0] = sp2 / (sp1 - V[0]);
    gmm::add(gmm::scaled(Y1, -1.0*dX[0]), Y2);
    gmm::copy(Y2, gmm::sub_vector(dX, gmm::sub_interval(1, N)));
		
    V[0] = 1.0 / (V[0] - sp1);
    gmm::copy(gmm::scaled(Y1, -1.0*V[0]),
	      gmm::sub_vector(V, gmm::sub_interval(1, N)) );
  }
  catch (...) {
    singular_tangent_matrix = true;
    cout << "Tangent matrix is close to be singular!" << endl;
  }


  try {	      
    if (singular_tangent_matrix) {
      plain_vector Q(1+N), R(1+N);
      sparse_matrix B(1+N, 1+N);

      gmm::copy(gmm::col_vector(T),
		gmm::sub_matrix(B, gmm::sub_interval(0, N), gmm::sub_interval(0, 1)));
      gmm::copy(A, gmm::sub_matrix(B, gmm::sub_interval(0, N), gmm::sub_interval(1, N)));
      gmm::copy(gmm::row_vector(V),
		gmm::sub_matrix(B, gmm::sub_interval(N, 1), gmm::sub_interval(0, 1+N)));
		
      gmm::copy(P, gmm::sub_vector(Q, gmm::sub_interval(0, N)));
      Q[N] = 0.0;

      gmm::clear(R);
      R[N] = 1.0;	
		    
      gmm::SuperLU_factor<scalar_type> SLU;
      SLU.build_with(B);
      SLU.solve(dX, Q);
      SLU.solve(V, R);
    }
    
    gmm::add(gmm::scaled(dX, -1.0), X);
    gmm::scale(V, 1.0/gmm::vect_norm2(V));

    if (noisy > 1) cout << "linear solver done" << endl;
  
    normdX = gmm::vect_norm2(dX);
    normP = gmm::vect_norm2(P);
    if (noisy) cout << "difference " << normdX << " residual " << normP << endl;
	      
    if ((normdX < tolerance) && (normP < residual)) convergence = 2;
    else if (normP > 1E10) convergence = 0;
  }
  catch (...) {
    convergence = 0;
    cout << "SuperLU failed!" << endl;
  }

  return convergence;

} // compute_Newton_iteration
	  

template<class VEC, class STDVEC, class ALLOCVEC, class MAT>
int elastostatic_problem::compute_test_functions
(const VEC &UN, const VEC &UT, const VEC &UT0, const VEC &LN, const VEC &LT,
 MAT &TST, const MAT &TST0,  const STDVEC &contact_nodes, const ALLOCVEC &gap, 
 scalar_type deltat, scalar_type friction_coef) {

  scalar_type r = PARAM.real_value("R", "Augmentation parameter");
  scalar_type alpha = PARAM.real_value("ALPHA", "Augmentation parameter");
  scalar_type beta = PARAM.real_value("BETA_T", "Augmentation parameter");
  scalar_type classboundary = PARAM.real_value
    ("CLASSBOUNDARY", "parameter for classification of the boundary"); 
  int noisy = PARAM.int_value("NOISY");
  size_type nbchange = 0, nblimit = 0;
//   size_type N = mesh.dim();
  
  for (size_type i = 0; i < LN.size(); ++i) {
//     if (LN[i]  > -1E-12)
//       cout << "node " << i << " : " << contact_nodes[i] << " no contact ";
//     else {
//       scalar_type lt =
// 	gmm::vect_norm2(gmm::sub_vector(LT, gmm::sub_interval(i*(N-1), N-1)));
//       if (lt < - friction_coef * LN[i] - 1E-12)
// 	cout << "node " << i << " : " << contact_nodes[i] << " sticking ";
//       else
// 	cout << "node " << i << " : " << contact_nodes[i] << " sliding ";
//     }
    
    TST(i,0) = (LN[i] - r * alpha * (UN[i] - gap[i]));
    if (deltat >= 0) {
      TST(i,1) = -friction_coef * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);
      TST(i,2) = friction_coef * LN[i] + LT[i] - r * beta * (UT[i] - UT0[i]);
    } else {
      TST(i,1) = -friction_coef * LN[i] + LT[i] + r * beta * (UT[i] - UT0[i]);
      TST(i,2) = friction_coef * LN[i] + LT[i] + r * beta * (UT[i] - UT0[i]);
    }

    if (noisy > 1) {
      cout << "node " << i << " : " << contact_nodes[i] << " ";
      cout << ((LN[i] - r * alpha * (UN[i] - gap[i])) <= 0)
	   << (-friction_coef * LN[i] + LT[i] - r * (UT[i] - UT0[i]) / deltat >= 0)
	   << (friction_coef * LN[i] + LT[i] - r * (UT[i] - UT0[i]) / deltat <= 0) << " ";
      for (size_type j = 0; j < 3; ++j){
	if (((TST0(i,j) < 0 && TST(i,j) > classboundary)
	     || (TST0(i,j) > 0 && TST(i,j) < -classboundary))
	    && (j == 0 || (j != 0 && TST(i,0) <= classboundary))) {
	  cout << " change of the sign: " << TST0(i,j) << " -> ";
	  ++nbchange;
	}
	else if (((TST(i,j) >= -classboundary && TST(i,j) <= classboundary)
	     || (TST(i,j) >= -classboundary && TST(i,j) <= classboundary))
	    && (j == 0 || (j != 0 && TST(i,0) <= classboundary))) {
	  cout << " close to a limit: ";
	  ++nblimit;
	}
	cout << TST(i,j) << " ";
      }
      cout << endl;

      if (noisy > 2)
	cout << LN[i] << " " << r * alpha * (UN[i] - gap[i]) << " "
	     << -friction_coef * LN[i] + LT[i] << " "
	     << friction_coef * LN[i] + LT[i] << " " << r * beta * (UT[i]-UT0[i]) << endl;
    }
  }

  if (noisy > 1) {
    cout << "Number of changes of characters = " << nbchange << endl;
    cout << "Number of test functions on a limit = " << nblimit << endl;
  }

  return nbchange + nblimit;
}


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
      (mim, mf_u, law_num == 5 ? 0.0 : p[0], p[1]);
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
      pb->penalization_coeff().set(1.0/p[0]);
    }
  }
  
  int nb_step = int(PARAM.int_value("NBSTEP"));
  scalar_type deltat = PARAM.real_value("DELTAT", "Time step");
 
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
  int nbc = int(jj);
  sparse_matrix MMBN(nbc, nbc), MMBT(nbc*(N-1), nbc*(N-1));
  plain_vector LN1(nbc), LT1(nbc*(N-1));
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
  gmm::dense_matrix<size_type> CH(nbc+1, 4);
  CH(0,0) = size_type(-1);
  

  getfem::mdbrick_Coulomb_friction<> FRICTION(*pINCOMP, BN, gap,
					      friction_coef, BT);
  FRICTION.set_r(r);
  FRICTION.set_alpha(alpha);
  // if (deltat < 0) FRICTION.set_beta(-1.0); 
  FRICTION.set_beta(1./deltat);
  FRICTION.set_CH(CH);
 
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
  size_type maxit = PARAM.int_value("MAXITER"); 
  int noisy = PARAM.int_value("NOISY");
  gmm::iteration iter;

  size_type continuation = PARAM.int_value("CONTINUATION", "Continuation method");
  int step_init = PARAM.int_value("STEP_INIT");

  scalar_type dz = PARAM.real_value("DIRICHLET_Z",
				    "Prescribed displacement in z");
  scalar_type dyv = PARAM.real_value("DIRICHLET_Y_SPEED",
				     "Prescribed velocity in y");


  bool converged = true;
  char strstep[100];
  size_type stot = gmm::mat_ncols(MS.tangent_matrix());
  scalar_type t = 0, deltat0 = deltat;
  plain_vector UN(gmm::mat_nrows(BN)), UT(gmm::mat_nrows(BT)), UT0(gmm::mat_nrows(BT));
  plain_vector U0 = U, U00 = U0;
  plain_vector X0(1+stot);
  bgeot::base_matrix TST0(nbc,3), TST(nbc,3);

  if (step_init > 0) { // load initial values
    sprintf(strstep, "step%d", step_init);
    gmm::vecload(datafilename + strstep + ".U",U0);
    gmm::vecload(datafilename + strstep + ".X",X0);
    t = X0[0];
    gmm::copy(gmm::sub_vector(X0, gmm::sub_interval(1, stot)), MS.state());
    Harwell_Boeing_load(datafilename + strstep + ".TST", TST0);
  }

  for (int step = step_init; step < nb_step; ++step) { /* Newton */
    
    if ((step > 2) && (continuation > 0)) break;  /* proceed to continuations */
    
    cout << "beginning of step " << step+1
	 << ", number of variables : " << final_model.nb_dof() << endl;
    
    t += (step < 2) ? 0.0 : deltat;
    if (noisy) cout << "t = " << t << " (deltat = " << deltat << ")" << endl;

    for (size_type i = 0; i < nb_dof_rhs; ++i) {
      F2[i*N+N-1] = (step < 2) ? dz*(step+1)/2.0 : dz;
      F2[i*N+N-2] = t*dyv;
    }
    final_model.rhs().set(F2);
    
    FRICTION.set_WT(gmm::scaled(U0, -1.0));
    
    
    iter = gmm::iteration(residual, noisy, maxit ? maxit : 40000);
    cout << "|U| = " << gmm::vect_norm2(MS.state()) << "\n";
    
    
    
    /* let the default non-linear solve (Newton) do its job */
    gmm::default_newton_line_search ls;
    getfem::standard_solve(MS, final_model, iter,
			   getfem::default_linear_solver(final_model), ls);
    
    

    if (law_num < 4) {
      pl->reset_unvalid_flag();
      final_model.compute_residual(MS);
      if (pl->get_unvalid_flag()) 
	GMM_WARNING1("The solution is not completely valid, the determinant "
		     "of the transformation is negative on "
		     << pl->get_unvalid_flag() << " gauss points");
    }
    
    if (law_num < 4)
      gmm::copy(pELAS_nonlin->get_solution(MS), U);
    else
      gmm::copy(pELAS_lin->get_solution(MS), U);
    
    gmm::copy(FRICTION.get_LN(MS), LN1);
    gmm::copy(FRICTION.get_LT(MS), LT1);
    
    gmm::mult(BN, U, UN);
    gmm::mult(BT, U, UT);
    gmm::mult(BT, U0, UT0);

    compute_test_functions(UN, UT, UT0, LN1, LT1, TST, TST0, contact_nodes,
			   gap, deltat, friction_coef);
    
    
    converged = iter.converged();
    if (!converged) break;
    
    gmm::copy(U0, U00);
    gmm::copy(U, U0);
    gmm::copy(TST, TST0);
    X0[0] = t;
    gmm::copy(MS.state(), gmm::sub_vector(X0, gmm::sub_interval(1, stot)));
    
    char s[100]; sprintf(s, "step%d", step+1);
    gmm::vecsave(datafilename + s + ".U",U);
    
    plain_vector VM(mf_vm.nb_dof());
    if (law_num < 4)
      pELAS_nonlin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
    else
      pELAS_lin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
    gmm::vecsave(datafilename + s + ".VM", VM);
    
    gmm::vecsave(datafilename + s + ".X", X0);
    Harwell_Boeing_save(datafilename + s + ".TST", TST);
    
    cout << "end of Step n∞ " << step+1 << " / " << nb_step << endl;
    
  } // Newton
  

  if (continuation > 0){
    GMM_ASSERT1(PARAM.int_value("DIRICHLET_VERSION") == 0,
		"Continuations only implemented for Dirichlet "
		"condition with multipliers");

    size_type maxcorrit = PARAM.int_value("MAXCORRITER");
    size_type ind =  final_model.first_ind(), sc = gmm::vect_size(final_model.get_CRHS());
    plain_vector T(stot);

    /* load initial values if necessary and compute derivative with respect to time */
    if (step_init < 3) step_init = 3;
    else {
      sprintf(strstep, "step%d", step_init-1);
      plain_vector X00(1+stot);
      gmm::vecload(datafilename + strstep + ".X",X00);
      deltat0 = X0[0] - X00[0];
      // (deltat0 < 0) ? FRICTION.set_beta(-1.0) : FRICTION.set_beta(1.0);
      FRICTION.set_beta(1./deltat0);
      gmm::vecload(datafilename + strstep + ".U",U00);  
      FRICTION.set_WT(gmm::scaled(U00, -1.0));
      
      for (size_type i = 0; i < nb_dof_rhs; ++i) {
	F2[i*N+N-1] = dz;
	F2[i*N+N-2] = t*dyv;
      }    
      final_model.rhs().set(F2);	    
      
      final_model.compute_tangent_matrix(MS);
    }
   
    gmm::copy(final_model.get_CRHS(),
	      gmm::sub_vector(T, gmm::sub_interval(ind, sc)));
    for (size_type i = 0; i < sc/N; ++i) 
      T[ind+i*N + N-1] = 0.;
    gmm::scale(T, -1.0/t);
    
    sprintf(strstep, "step%d", step_init);
    gmm::vecload(datafilename + strstep + ".U",U0);
    

    if (continuation < 2){ /* Euler-Newton */
      scalar_type deltat_max = PARAM.real_value("DELTAT_MAX", "Maximal step length");
      scalar_type deltat_min = PARAM.real_value("DELTAT_MIN", "Minimal step length");
      scalar_type deltat_inc_fac = PARAM.real_value("DELTAT_INC_FAC", "Increment factor");
      scalar_type deltat_dec_fac = PARAM.real_value("DELTAT_DEC_FAC", "Decrement factor");

      int corrections;
      double rcond;
      plain_vector V0(stot), V(stot);

      gmm::SuperLU_solve(MS.tangent_matrix(), V0, T, rcond); // compute initial tangent

      // (deltat < 0) ? FRICTION.set_beta(-1.0) : FRICTION.set_beta(1.0);
      FRICTION.set_beta(1./deltat);
      
      for (int step = step_init; step < nb_step; ++step) {
	cout << "beginning of step " << step+1 
	     << ", number of variables : " << final_model.nb_dof() << endl;
	  
	FRICTION.set_WT(gmm::scaled(U0, -1.0));
	corrections = 1;

	while (1) {
	  /* predict next point */
	  gmm::add(gmm::scaled(V0, -deltat), MS.state()); 

	  t = X0[0] + deltat;
	  if (noisy) cout << "t = " << t << " (deltat = " << deltat << ")" << endl;
	  
	  for (size_type i = 0; i < nb_dof_rhs; ++i) {
	    F2[i*N+N-1] = dz;
	    F2[i*N+N-2] = t*dyv;
	  }
	  final_model.rhs().set(F2);
	  
	  
	  iter = gmm::iteration(residual, noisy, maxcorrit ? maxcorrit : 40000);
	  cout << "|U| = " << gmm::vect_norm2(MS.state()) << "\n";
	  
	  
	  
	  /* let the default non-linear solve (Newton) do its job */
	  gmm::default_newton_line_search ls;
	  getfem::standard_solve(MS, final_model, iter,
				 getfem::default_linear_solver(final_model), ls);


	  
	  if (law_num < 4) {
	    pl->reset_unvalid_flag();
	    final_model.compute_residual(MS);
	    if (pl->get_unvalid_flag()) 
	      GMM_WARNING1("The solution is not completely valid, the determinant "
			   "of the transformation is negative on "
			   << pl->get_unvalid_flag() << " gauss points");
	  }
	  
	  if (law_num < 4)
	    gmm::copy(pELAS_nonlin->get_solution(MS), U);
	  else
	    gmm::copy(pELAS_lin->get_solution(MS), U);
	  
	  gmm::copy(FRICTION.get_LN(MS), LN1);
	  gmm::copy(FRICTION.get_LT(MS), LT1);
	  
	  gmm::mult(BN, U, UN);
	  gmm::mult(BT, U, UT);
	  gmm::mult(BT, U0, UT0);
	  
	  int nb_nodes = compute_test_functions(UN, UT, UT0, LN1, LT1, TST,TST0, 
				       	contact_nodes, gap, deltat, friction_coef);
	  
	  
	  gmm::SuperLU_solve(MS.tangent_matrix(), V, T, rcond);
	  
	  scalar_type SP = (gmm::vect_sp(V, V0) / (gmm::vect_norm2(V) * gmm::vect_norm2(V0)));
	  cout << "SP = " << SP << endl;
	  
	  if (iter.converged() && ((SP > 0.9) || ((SP <= 0.9) && (nb_nodes < 2))))
	    break;
	  
	  if (deltat > deltat_min) {
	    deltat = (deltat*deltat_dec_fac > deltat_min)
	      ? deltat*deltat_dec_fac : deltat_min;
	    // if (deltat < 0) FRICTION.set_beta(-1.0); 
	    FRICTION.set_beta(1./deltat);
	    ++corrections;
	    gmm::copy(gmm::sub_vector(X0, gmm::sub_interval(1, stot)), MS.state());
	  }
	} // while (1)
   
	converged = iter.converged();
	if (!converged) break;
	
	deltat0 = deltat;
	gmm::copy(U0, U00);
	gmm::copy(U, U0); 
	gmm::copy(V, V0);
	gmm::copy(TST, TST0);
	X0[0] = t;
	gmm::copy(MS.state(), gmm::sub_vector(X0, gmm::sub_interval(1, stot)));

	char s[100]; sprintf(s, "step%d", step+1);
	gmm::vecsave(datafilename + s + ".U",U);
   
	plain_vector VM(mf_vm.nb_dof());
	if (law_num < 4)
	  pELAS_nonlin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
	else
	  pELAS_lin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
	gmm::vecsave(datafilename + s + ".VM", VM);

	gmm::vecsave(datafilename + s + ".X", X0);
	Harwell_Boeing_save(datafilename + s + ".TST", TST);

	if ((deltat < deltat_max) && (corrections == 1) && (iter.converged())) {
	  deltat = (deltat*deltat_inc_fac < deltat_max) ? deltat*deltat_inc_fac : deltat_max;
	  // if (deltat < 0) FRICTION.set_beta(-1.);
	  FRICTION.set_beta(1./deltat);
	}

	cout << "end of Step n∞ " << step+1 << " / " << nb_step << endl;	
      } // for (int step = 3; step < nb_step; ++step)

    } // Euler-Newton
  

    else { /* Moore-Penrose */
      scalar_type tolerance = PARAM.real_value("TOLERANCE");
      if (tolerance == 0.) tolerance = 1e-10;
      scalar_type h = PARAM.real_value("H", "Actual step length");
      scalar_type h_init = PARAM.real_value("H_INIT", "Initial step length");
      scalar_type h_max = PARAM.real_value("H_MAX", "Maximal step length");
      scalar_type h_min = PARAM.real_value("H_MIN", "Minimal step length");
      scalar_type h_inc_fac = PARAM.real_value("H_INC_FAC", "Increment factor");
      scalar_type h_dec_fac = PARAM.real_value("H_DEC_FAC", "Decrement factor");

      int convergence = 1, corrections;
      plain_vector V0(1+stot), V(1+stot), X(1+stot);
       gmm::dense_matrix<size_type> CH0(nbc+1, 4);
       CH0(0,0) = size_type(-1);

      if (step_init < 4) {
	gmm::clear(V0);
	V0[0] = 1.0;
	compute_tangent(T, MS.tangent_matrix(), V0, noisy);
      }
      else {
	gmm::vecload(datafilename + strstep + ".V",V0);
      }

      gmm::copy(X0, X);
      gmm::copy(V0, V);
      
      for (int step = step_init; step < nb_step; ++step) {
	cout << "beginning of step " << step+1
	     << ", number of variables : " << final_model.nb_dof() << endl;
	
	FRICTION.set_WT(gmm::scaled(U0, -1.0));
	corrections = 1;
      
	cout << " |U| = " << gmm::vect_norm2(MS.state()) << "\n";
	
	while (1) {
	  
	  /* predict next point */
	  gmm::add(gmm::scaled(V0, h), X);
	  cout << "h = " << h << endl;
	  
	  for (size_type j = 0; j < maxcorrit; ++j) {
	    if (noisy) cout << "iter " << j << endl;
	    if (noisy > 1) cout << "starting computing matrices" << endl;
	    
	    gmm::copy(gmm::sub_vector(X, gmm::sub_interval(1, stot)), MS.state());

	    deltat = X[0]-X0[0];
	    // (deltat < 0) ? FRICTION.set_beta(-1.0) : FRICTION.set_beta(1.0); 
	    FRICTION.set_beta(1./deltat);
	  	    
	    t = X[0];
	    for (size_type i = 0; i < nb_dof_rhs; ++i) {
	      F2[i*N+N-1] = dz;
	      F2[i*N+N-2] = t*dyv;
	    }
	    final_model.rhs().set(F2);	    
	    
	    final_model.compute_tangent_matrix(MS);
	    final_model.compute_residual(MS);

	    if (noisy) cout << "t = " << t << " (deltat = " << deltat << ")" << endl;

	    if (noisy > 2) {
	      if (law_num < 4)
		gmm::copy(pELAS_nonlin->get_solution(MS), U);
	      else
		gmm::copy(pELAS_lin->get_solution(MS), U);
	      
	      gmm::copy(FRICTION.get_LN(MS), LN1);
	      gmm::copy(FRICTION.get_LT(MS), LT1);
	      
	      gmm::mult(BN, U, UN);
	      gmm::mult(BT, U, UT);
	      gmm::mult(BT, U0, UT0);
	      
 	      compute_test_functions(UN, UT, UT0, LN1, LT1, TST, TST0,
				     contact_nodes, gap, deltat, friction_coef);
	    }
	    
	    convergence = compute_Newton_iteration(T, MS.tangent_matrix(), V, X,
						  MS.residual(), tolerance, residual, noisy);
	    cout << "SP = " << gmm::vect_sp(V, V0) << endl;
	  
	    if (convergence < 1) break; /* an error occured */

	    if (convergence > 1) {/* convergence */
	      
	      gmm::copy(gmm::sub_vector(X, gmm::sub_interval(1, stot)), MS.state());
	   
	      CH(0,0) = size_type(-1);
	      FRICTION.set_CH(CH);
	      
	      deltat = X[0]-X0[0];
	      // (deltat < 0.) ? FRICTION.set_beta(-1.) : FRICTION.set_beta(1.); 
	      FRICTION.set_beta(1./deltat);
	      
	      t = X[0];
	      for (size_type i = 0; i < nb_dof_rhs; ++i) {
		F2[i*N+N-1] = dz;
		F2[i*N+N-2] = t*dyv;
	      }
	      
	      final_model.rhs().set(F2);	    
	      
	      final_model.compute_tangent_matrix(MS);
	      final_model.compute_residual(MS);

	      scalar_type res = gmm::vect_norm2(MS.residual());
	      cout << "final residual = " << res << endl;
	      if (res < residual)
		compute_tangent(T, MS.tangent_matrix(), V, noisy);
	      else convergence = 1;
		
	      break;

	    } // if (convergence > 1)
	  
	  } // for (size_type j = 1; j < maxcorrit; ++j)
	
	  if (convergence > 1) { /* convergence */
	    if (law_num < 4) {
	      pl->reset_unvalid_flag();
	      final_model.compute_residual(MS);
	      if (pl->get_unvalid_flag()) 
		GMM_WARNING1("The solution is not completely valid, the determinant "
			     "of the transformation is negative on "
			     << pl->get_unvalid_flag() << " gauss points");
	    }
	    
	    if (law_num < 4)
	      gmm::copy(pELAS_nonlin->get_solution(MS), U);
	    else
	      gmm::copy(pELAS_lin->get_solution(MS), U);
	    
	    deltat = X[0]-X0[0];

	    if (deltat < 0.0) cout << "Warning: deltat is negative!" << endl;
	    
	    gmm::copy(FRICTION.get_LN(MS), LN1);
	    gmm::copy(FRICTION.get_LT(MS), LT1);

	    gmm::mult(BN, U, UN);
	    gmm::mult(BT, U, UT);
	    gmm::mult(BT, U0, UT0);
	     
	    int nb_nodes = compute_test_functions(UN, UT, UT0, LN1, LT1, TST, TST0,
					      contact_nodes, gap, deltat, friction_coef);

	    
	    scalar_type SP = gmm::vect_sp(V, V0);
	    cout << "SP = " << SP << endl;
	    
	    if ((SP > 0.9) || ((SP > 0.7) && (nb_nodes < 2)))
	      break;
// 	    else {
// 	      cout << "Do you want to accept this increment? (1, 0, -1) ";
// 	      int sign;
// 	      cin >> sign;
// 	      if (sign != 0) {
// 		if (sign < 0) gmm::scale(V, -1.0);
// 		h = h_init;
// 		++corrections;
// 		break;
// 	      }
// 	    }
	    
	  } // if (convergence > 1)
	  
	  if (h > h_min) {
	    FRICTION.set_CH(CH0);
	    gmm::copy(X0, X);
	    gmm::copy(V0, V);
	    h = (h*h_dec_fac > h_min) ? h*h_dec_fac : h_min;
	    ++corrections;
	  }
	  else { /* try to change the Jacobian */
	      
	    cout << "Set changes: ";
	    for (int i = 0; i < nbc+1; ++i) {
	      cin >> CH(i,0) >> CH(i,1) >> CH(i,2) >> CH(i,3);
	      if (CH(i,0) == size_type(-1)) break;
	    }
	    
	    for (int i = 0; i < nbc; ++i) {
	      CH0(i,0) = i;
	      CH0(i,1) = (TST0(i,0) <= 0) ? 1 : 0;
	      if (TST0(i,1) < 0) { CH0(i,2) = 0; CH0(i,3) = 1; }
	      else if (TST0(i,2) > 0) { CH0(i,2) = 1; CH0(i,3) = 0; }
	      else { CH0(i,2) = 1; CH0(i,3) = 1; }
	    }
	    for (int i = 0; i < nbc; ++i) {
	      if (CH(i,0) == size_type(-1)) break;
	      CH0(CH(i,0), 1) = CH(i,1);
	      CH0(CH(i,0), 2) = CH(i,2);
	      CH0(CH(i,0), 3) = CH(i,3);
	    }
	    CH0(nbc,0) = size_type(-1);
	    
	    gmm::copy(gmm::sub_vector(X0, gmm::sub_interval(1, stot)), MS.state());
	    gmm::copy(V0, V);
	    
	    t = X0[0];
	    for (size_type i = 0; i < nb_dof_rhs; ++i) {
	      F2[i*N+N-1] = dz;
	      F2[i*N+N-2] = t*dyv;
	    }
	    
	    final_model.rhs().set(F2);
	    FRICTION.set_WT(gmm::scaled(U00, -1.0));
	    // (deltat0 < 0.0) ? FRICTION.set_beta(-1.0) : FRICTION.set_beta(1.0);
	    FRICTION.set_beta(1./deltat0);
	    
	    final_model.compute_tangent_matrix(MS);
	   
	    FRICTION.set_CH(CH0);
	    FRICTION.do_compute_tangent_matrix(MS, 0, 0);
	    compute_tangent(T, MS.tangent_matrix(), V0, noisy);
	   
	    CH0(0,0) = size_type(-1);
	    FRICTION.set_CH(CH0);
	    FRICTION.set_WT(gmm::scaled(U0, -1.0));
	      
	    cout << "Choose a sign for the tangent (1, 0, -1): ";
	    int sign;
	    cin >> sign;
	    if (sign == 0) {
	      convergence = 1; /* no convergence */
	      break;
	    }
	    if (sign < 0) gmm::scale(V0, -1.0);
	    cout << "SP = " << gmm::vect_sp(V, V0) << endl;
	    
	    gmm::copy(X0, X);
	    gmm::copy(V0, V);
	    h = h_init;

	  } // try to change the Jacobian
	
	} // while (1)
	    
	converged = (convergence > 1);
	if (!converged) break;

	CH0(0,0) = size_type(-1);
	deltat0 = deltat;
	gmm::copy(U0, U00);
	gmm::copy(U, U0); 
	gmm::copy(X, X0);
	gmm::copy(V, V0);
	gmm::copy(TST, TST0);

	char s[100]; sprintf(s, "step%d", step+1);
	gmm::vecsave(datafilename + s + ".U",U);
	
	plain_vector VM(mf_vm.nb_dof());
	if (law_num < 4)
	  pELAS_nonlin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
	else
	  pELAS_lin->compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
	gmm::vecsave(datafilename + s + ".VM", VM);

	gmm::vecsave(datafilename + s + ".X", X);
	gmm::vecsave(datafilename + s + ".V", V);
	Harwell_Boeing_save(datafilename + s + ".TST", TST);
	
	if ((h < h_max) && (corrections == 1) && (convergence > 1)) {
	  h = (h*h_inc_fac < h_max) ? h*h_inc_fac : h_max;
	}

	cout << "end of Step n∞ " << step+1 << " / " << nb_step << endl;
      
      } // for (int step = 3; step < nb_step; ++step)
    
    } // Moore-Penrose

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
