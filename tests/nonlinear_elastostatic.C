/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004 Yves Renard, Julien Pommier.                    */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

/**
 * Nonlinear Elastostatic problem (large strain).
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <getfem_nonlinear_elasticity.h>
#include <gmm.h>

/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector;
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
  structure for the elastostatic problem
*/
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type p1, p2, p3;    /* elastic coefficients.                        */

  scalar_type residu;        /* max residu for the iterative solvers         */

  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  elastostatic_problem(void) : mf_u(mesh), mf_rhs(mesh), mf_coef(mesh) {}
};

void test(scalar_type a, scalar_type b, scalar_type c) {
  getfem::Hooke_hyperelastic_law AHL;
  base_matrix L(3,3), L2(3,3), gradU(3,3), gradU2(3,3), dgradU(3,3), Sigma(3,3), Sigma2(3,3), NL(3,3), NL2(3,3),A(3,3),B(3,3);
  getfem::base_tensor dNL(3,3,3,3), tt(3,3,3,3);
  base_vector params(2);

  params[0] = a; params[1] = b;
  gmm::fill_random(gradU); gmm::scale(gradU, c); 
  gmm::fill_random(dgradU); gmm::scale(dgradU,0.000001);
  
  //gmm::clear(dgradU); dgradU(0,0) = 0.0001;

  gmm::copy(gradU,gradU2); gmm::add(dgradU,gradU2);

  gmm::mult(gmm::transposed(gradU), gradU, L);
  gmm::add(gradU, L);
  gmm::add(gmm::transposed(gradU), L);
  gmm::scale(L, scalar_type(0.5));
  gmm::add(gmm::identity_matrix(), gradU);

  gmm::mult(gmm::transposed(gradU2), gradU2, L2);
  gmm::add(gradU2, L2);
  gmm::add(gmm::transposed(gradU2), L2);
  gmm::scale(L2, scalar_type(0.5));
  gmm::add(gmm::identity_matrix(), gradU2);

  AHL.sigma(L, Sigma, params);
  AHL.sigma(L2, Sigma2, params);

  AHL.grad_sigma(L, tt, params);
  size_type N=3;
  for (size_type i = 0; i < N; ++i)
    for (size_type j = 0; j < N; ++j) {
      scalar_type aux(0), aux2(0);
      for (size_type k = 0; k < N; ++k) {
	aux += gradU(i, k) * Sigma(k, j);
	aux2 += gradU2(i, k) * Sigma2(k, j);
      }
      NL(i,j) = aux;
      NL2(i,j) = aux2;
    }


  for (size_type n = 0; n < N; ++n)
    for (size_type m = 0; m < N; ++m)
      for (size_type l = 0; l < N; ++l)
	for (size_type k = 0; k < N; ++k) {
	  //scalar_type aux = (k == l) ? Sigma(m, l) : 0.0;
	  //scalar_type aux = (m == l) ? Sigma(k,n) : 0.0;
	  scalar_type aux = (k == n) ? Sigma(m,l) : 0.0;
	  for (size_type j = 0; j < N; ++j)
	    for (size_type i = 0; i < N; ++i) {
	      aux += gradU(n ,j) * gradU(k, i) * tt(j, m, i, l);
	    }
	  dNL(n, m, k, l) = aux;
	}

  gmm::mult(gmm::transposed(gradU),dgradU,A);
  for (size_type n = 0; n < N; ++n)
    for (size_type m = 0; m < N; ++m) {
      B(n,m)=0;
      for (size_type j = 0; j < N; ++j)
	for (size_type i = 0; i < N; ++i) 
	  B(n,m) += tt(n,m,i,j)*A(i,j);
    }
  gmm::clear(A); gmm::mult(gradU,B,A);
  gmm::clear(B); gmm::mult(dgradU,Sigma,B);
  gmm::add(B,A);

  for (size_type n = 0; n < N; ++n)
    for (size_type m = 0; m < N; ++m) {
      scalar_type aux = 0;
      for (size_type j = 0; j < N; ++j)
	for (size_type i = 0; i < N; ++i) {
	  aux += dNL(i,j,n,m)*dgradU(i,j);
	}
      printf("%d,%d: %12.7f %12.7f %12.9f %12.9f %12.9f\n", n,m,NL(n,m),NL2(n,m),NL2(n,m)-NL(n,m),aux,A(n,m));
    }
  printf("\n");
}

namespace getfem {
template <typename MODEL_STATE> void
  nl_solve(MODEL_STATE &MS, mdbrick_abstract<MODEL_STATE> &problem,
	   gmm::iteration &iter) {

    typedef typename MODEL_STATE::vector_type VECTOR;
    typedef typename MODEL_STATE::value_type value_type;
    typedef typename MODEL_STATE::tangent_matrix_type T_MATRIX;
    typedef typename MODEL_STATE::constraints_matrix_type C_MATRIX;
    typedef typename gmm::number_traits<value_type>::magnitude_type mtype;

    // TODO : take iter into account for the Newton. compute a consistent 
    //        max residu.
    //        detect the presence of multipliers before using a preconditioner

    size_type ndof = problem.nb_dof();
    bool is_linear = problem.is_linear();
    //mtype alpha, alpha_min=mtype(1)/mtype(32), alpha_mult=mtype(3)/mtype(4);
    mtype alpha, alpha_min=mtype(1)/mtype(8), alpha_mult=mtype(2)/mtype(3);
    mtype alpha_max_ratio(1);
    dal::bit_vector mixvar;
    gmm::iteration iter_linsolv0 = iter;
    iter_linsolv0.set_maxiter(10000);
    if (!is_linear) { iter_linsolv0.reduce_noisy(); iter_linsolv0.set_resmax(iter.get_resmax()/10000.0); }

    MS.adapt_sizes(problem);
    /*if (!is_linear) gmm::fill_random(MS.state()); 
      else */
    //gmm::clear(MS.state());
    //gmm::fill_random(MS.state()); 
    problem.compute_residu(MS);
    problem.compute_tangent_matrix(MS);
    MS.compute_reduced_system();
    mtype act_res = MS.reduced_residu_norm(), act_res_new(0);
    cout << "residu initial: " << gmm::vect_norm2(MS.residu()) << ", " 
	 << gmm::vect_norm2(MS.reduced_residu()) << ", |U0|" << gmm::vect_norm2(MS.state()) << "\n";
    while (!iter.finished(act_res)) {
      gmm::iteration iter_linsolv = iter_linsolv0;
      VECTOR d(ndof), dr(gmm::vect_size(MS.reduced_residu()));

      if (!(iter.first())) {
	problem.compute_tangent_matrix(MS);
	MS.compute_reduced_system();
      }

      gmm::ildlt_precond<T_MATRIX> P(MS.reduced_tangent_matrix());
      gmm::cg(MS.reduced_tangent_matrix(), dr, 
	      gmm::scaled(MS.reduced_residu(), value_type(-1)),
	      P, iter_linsolv);
      if (!iter_linsolv.converged()) DAL_WARNING(2,"cg did not converge!");
      MS.unreduced_solution(dr,d);

      VECTOR stateinit(ndof);
      gmm::copy(MS.state(), stateinit);
      
      {
	problem.compute_residu(MS);
	MS.compute_reduced_system();
	scalar_type r0 = MS.reduced_residu_norm();
	cout << "R0 === " << r0 << "\n";
	VECTOR W(gmm::vect_size(dr));
	gmm::mult(MS.reduced_tangent_matrix(), dr,W);
	scalar_type w = 2*gmm::vect_sp(MS.reduced_residu(), W);

	for (alpha = mtype(1); alpha >= alpha_min/100; alpha *= alpha_mult) {
	  gmm::add(stateinit, gmm::scaled(d, -alpha), MS.state());
	  problem.compute_residu(MS);
	  MS.compute_reduced_system(); // The whole reduced system do not
	  // have to be computed, only the RHS. To be adapted.
	  act_res_new = dal::sqr(MS.reduced_residu_norm());
	  printf("%+12.5g  %+12.5g  %+12.5g\n", -alpha, act_res_new, w*(-alpha)+dal::sqr(r0));
	}
	for (alpha = mtype(1); alpha >= alpha_min/100; alpha *= alpha_mult) {
	  gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	  problem.compute_residu(MS);
	  MS.compute_reduced_system(); // The whole reduced system do not
	  // have to be computed, only the RHS. To be adapted.
	  act_res_new = dal::sqr(MS.reduced_residu_norm());
	  printf("%+12.5g  %+12.5g  %+12.5g\n", alpha, act_res_new, w*(alpha)+dal::sqr(r0));
	}
	gmm::copy(stateinit, MS.state());
	problem.compute_residu(MS);
	MS.compute_reduced_system();
	r0 = MS.reduced_residu_norm();
	cout << "R0 === " << r0 << "\n";
      }



      for (alpha = mtype(1); alpha >= alpha_min; alpha *= alpha_mult) {
	gmm::add(stateinit, gmm::scaled(d, alpha), MS.state());
	problem.compute_residu(MS);
	MS.compute_reduced_system(); // The whole reduced system do not
	// have to be computed, only the RHS. To be adapted.
	act_res_new = MS.reduced_residu_norm();
	if (act_res_new <= act_res * alpha_max_ratio) break;

	gmm::add(stateinit, gmm::scaled(d, -alpha), MS.state());
	problem.compute_residu(MS);
	MS.compute_reduced_system(); // The whole reduced system do not
	// have to be computed, only the RHS. To be adapted.
	act_res_new = MS.reduced_residu_norm();
	if (act_res_new <= act_res * alpha_max_ratio) {
	  cout << "WWWWWRONG DIRECTION ALPHA < 0 is BETTER THAN ALPHA > 0 !!!!\n";
	  break;
	}
      }
      cout << "alpha = " << alpha << ", |U| = " << gmm::vect_norm2(MS.state()) << ", residu = " << act_res_new << "\n";
      act_res = act_res_new; ++iter;
    }
  }
}

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void elastostatic_problem::init(void) {
  const char *MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  const char *FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
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

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;

  p1 = PARAM.real_value("P1", "First Elastic coefficient");
  p2 = PARAM.real_value("P2", "Second Elastic coefficient");
  p3 = PARAM.real_value("P3", "Third Elastic coefficient");
  
  mf_u.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);

  mf_u.set_finite_element(mesh.convex_index(), pf_u, ppi);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
    if (!pf_u->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM "
		<< data_fem_name << ". In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u, ppi);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name), ppi);
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0), ppi);

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::convex_face_ct border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
       it != border_faces.end(); ++it) {
    assert(it->f != size_type(-1));
    base_node un = mesh.normal_of_face_of_convex(it->cv, it->f);
    un /= gmm::vect_norm2(un);
    if (dal::abs(un[N-1] - 1.0) < 1.0E-7) { // new Neumann face
      mf_u.add_boundary_elt(NEUMANN_BOUNDARY_NUM, it->cv, it->f);
    } else if (dal::abs(un[N-1] + 1.0) < 1.0E-7) {
      mf_u.add_boundary_elt(DIRICHLET_BOUNDARY_NUM, it->cv, it->f);
    }
  }
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool elastostatic_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  // Linearized elasticity brick.
  getfem::Hooke_hyperelastic_law l;

  if (0) {
    cout << "test derivees Hooke_hyperelastic_law\n";
    base_vector param(2); param[0] = 1.; param[1] = .7423;
    for (size_type itest = 0; itest < 100; ++itest) {
      base_matrix L(3,3), L2(3,3); 
      gmm::fill_random(L); //gmm::copy(L2,L); gmm::add(transposed(L2),L);
      base_matrix DL(3,3); 
      gmm::fill_random(DL);
      //gmm::fill_random(L2); gmm::copy(L2,DL); gmm::add(transposed(L2),DL); gmm::scale(DL,0.1);
      base_matrix sigma1(3,3), sigma2(3,3);
      getfem::base_tensor tdsigma(3,3,3,3);
      base_matrix dsigma(3,3);
      gmm::copy(L,L2);
      gmm::add(DL,L2);
      l.sigma(L, sigma1, param);l.sigma(L2, sigma2, param);
      l.grad_sigma(L,tdsigma,param);
      for (size_type i=0; i < 3; ++i) {
	for (size_type j=0; j < 3; ++j) {
	  dsigma(i,j) = 0;
	  for (size_type k=0; k < 3; ++k) {
	    for (size_type m=0; m < 3; ++m) {
	      dsigma(i,j) += tdsigma(i,j,k,m)*DL(k,m);
	    }
	  }
	  sigma2(i,j) -= sigma1(i,j);
	  if (dal::abs(dsigma(i,j) - sigma2(i,j)) > 1e-13) {
	    cout << "erreur derivees i=" << i << ", j=" << j << ", dsigma=" << dsigma(i,j) << ", var sigma = " << sigma2(i,j) << "\n";
	  }
	}
      }
    }
  }

  //getfem::Mooney_Rivlin_hyperelastic_law l;
  getfem::mdbrick_nonlinear_elasticity<>  ELAS(l, mf_u, mf_coef, p1, p2);

  // Defining the volumic source term.
  base_vector f(N); 
  f[0] = PARAM.real_value("FORCEX","Amplitude of the force");
  f[1] = PARAM.real_value("FORCEY","Amplitude of the force");
  if (N>2)
    f[2] = PARAM.real_value("FORCEZ","Amplitude of the force");
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  
  // Volumic source term brick.
  int nb_step = PARAM.int_value("NBSTEP");

  // Generic solve.
  getfem::standard_model_state MS;
  size_type maxit = PARAM.int_value("MAXITER"); 
  gmm::iteration iter;
  gmm::clear(MS.state());
  for (int step = 0; step < nb_step; ++step) {
    plain_vector DF(F);
    gmm::copy(gmm::scaled(F, (step+1.)/(scalar_type)nb_step), DF);
    getfem::mdbrick_source_term<> VOL_F(ELAS, mf_rhs, DF);
    gmm::clear(DF);
    getfem::mdbrick_Dirichlet<> final_model(VOL_F, mf_rhs,
					    DF, DIRICHLET_BOUNDARY_NUM,
					    PARAM.int_value("USE_MULTIPLIERS"));
    cout << "step " << step << ", number of variables : " << final_model.nb_dof() << endl;
    iter = gmm::iteration(residu, 1, maxit ? maxit : 40000);
    cout << "|U0| = " << gmm::vect_norm2(MS.state()) << "\n";
    getfem::nl_solve(MS, final_model, iter);
  }

  // Solution extraction
  ELAS.get_solution(MS, U);
  
  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb); // to debug ...

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  try {    
    elastostatic_problem p;
    test(1,0, 0);
    test(0,1, 0);
    test(1,0, 1);
    test(0,1, 1);
    test(.35,.54, 1.);

    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    p.mf_u.write_to_file(p.datafilename + ".mf", true);
    p.mf_rhs.write_to_file(p.datafilename + ".mfd", true);
    plain_vector U(p.mf_u.nb_dof());
    if (p.PARAM.int_value("VTK_EXPORT")) {
      if (!p.solve(U)) 
	//DAL_THROW(dal::failure_error,"Solve has failed");
	cerr << "Solve has failed\n";
      cout << "export to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk",
			     p.PARAM.int_value("VTK_EXPORT")==1);
      exp.exporting(p.mf_u); 
      exp.write_point_data(p.mf_u, U, "elastostatic_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << ".vtk -f ExtractVectorNorm -f "
	"WarpVector -m BandedSurfaceMap -m Outline\n";
    }
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
