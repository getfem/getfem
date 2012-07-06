/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien  Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
#include "getfem/getfem_Coulomb_friction.h"
#include "getfem/getfem_superlu.h"
#include "gmm/gmm.h"
#include "getfem/getfem_import.h"  //ajout des fonctionnalités d'importation de maillage (pfe)

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
  structure for the elastostatic  problem
*/
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, FRICTION_BOUNDARY_NUM = 1};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim;       /* the integration methods */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure.                   */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  getfem::mesh_fem mf_vm;    /* mesh_fem used for the VonMises stress        */
  scalar_type p1, p2, p3;    /* elastic coefficients.                        */
  scalar_type LX, LY, LZ;    /* system dimensions   */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type residual;        /* max residual for the iterative solvers         */
  std::string datafilename;
  bgeot::md_param PARAM;


  bool solve(plain_vector &U);
  void init(void);
  elastostatic_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh), mf_rhs(mesh), mf_coef(mesh), mf_vm(mesh) {}
//   elastostatic_problem(void) : mim(mesh), mf_u(mesh), mf_p(mesh), mf_rhs(mesh), mf_coef(mesh) {}
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
  std::string mesh_filename = PARAM.string_value("MESH_FILE_NAME");
  size_type MESH_CATEGORY = PARAM.int_value("MESH_CATEGORY");
  size_type REFINE_MESH = PARAM.int_value("REFINE_MESH");
  size_type N = 0;

  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  /* First step : build the mesh */

    bgeot::pgeometric_trans pgt =
    bgeot::geometric_trans_descriptor(MESH_TYPE);

  // According to the value of MESH_CATEGORY (user parameter), getfem builds the mesh or imports the mesh defined in the file MESH_FILE_NAME

  if (MESH_CATEGORY == 0) {

    N = pgt->dim();
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(),
  	    PARAM.int_value("NX", "Nomber of space steps "));
    nsubdiv[1] = PARAM.int_value("NY") ? PARAM.int_value("NY") : nsubdiv[0];
    if (N>2) nsubdiv[2] = PARAM.int_value("NZ") ? PARAM.int_value("NZ") : nsubdiv[0];
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
  			    PARAM.int_value("MESH_NOISED") != 0);

    LX = PARAM.real_value("LX", "Length along X axis");
    LY = PARAM.real_value("LY", "Length along Y axis");
    LZ = PARAM.real_value("LZ", "Length along Z axis");
    mu = PARAM.real_value("MU", "Lamé coefficient mu");
    lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    if (N>1) { M(0,1) = PARAM.real_value("INCLINE") * PARAM.real_value("LY"); }

    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);

    }

  else if (MESH_CATEGORY == 1) {

   // CHARGEMENT DU MAILLAGE ISSU DE MATLAB PDE TOOL

   getfem::import_mesh(mesh_filename,mesh);
   N = mesh.dim();

   LX = PARAM.real_value("LX", "Length along X axis");
   LY = PARAM.real_value("LY", "Length along Y axis");
   LZ = PARAM.real_value("LZ", "Length along Z axis");
   mu = PARAM.real_value("MU", "Lamé coefficient mu");
   lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");

   }



  // Refine the mesh "REFINE_MESH" times  only if REFINE_MESH (user parameter) > 0

  if (REFINE_MESH > 0 ) {

    scalar_type e = 0.26; // parameter defining the refining area
    for (size_type nb_refine = 0; nb_refine < REFINE_MESH; ++nb_refine) {

       dal::bit_vector bv;
        for (dal::bv_visitor iter1(mesh.convex_index()); !iter1.finished(); ++iter1)
        { if (mesh.points_of_convex(iter1)[0][N-2]>(LX - LX*e) && mesh.points_of_convex(iter1)[0][N-1]<LY*e)
	    { bv.add(iter1);}
        }

       cout << "nb elt refined : " << bv.card() << endl;
       mesh.Bank_refine(bv);

      e = e / 2.;
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
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
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
    } else if (mesh.points_of_convex(it.cv())[0][N-1]<0.003 && mesh.points_of_convex(it.cv())[1][N-1]<0.003 && mesh.points_of_convex(it.cv())[2][N-1]<0.003) {
      mesh.region(FRICTION_BOUNDARY_NUM).add(it.cv(),it.f());
    }
  }
}



//  template <typename VEC1, typename VEC2>
// void calcul_von_mises(const getfem::mesh_fem &mf_u, const VEC1 &U,
// 		      const getfem::mesh_fem &mf_vm, VEC2 &VM,
// 		      scalar_type mu=1) {
//   /* DU=gf_compute(mfu,U,'gradient',mf_vm);
//
//   from the derivative, we compute the von mises stress
//   VM=zeros(1,gf_mesh_fem_get(mf_vm,'nbdof'));
//   N=gf_mesh_get(m,'dim');
//   for i=1:size(DU,3),
//   t=DU(:,:,i);
//   E=(t+t')/2;
//   VM(i) = sum(E(:).^2) - (1./N)*sum(diag(E))^2;
//   end;
//   VM = 4*pde.mu{1}^2*VM;
//   */
//   assert(mf_vm.get_qdim() == 1);
//   unsigned N = mf_u.linked_mesh().dim(); assert(N == mf_u.get_qdim());
//   std::vector<scalar_type> DU(mf_vm.nb_dof() * N * N);
//
//   getfem::compute_gradient(mf_u, mf_vm, U, DU);
//
//   gmm::resize(VM, mf_vm.nb_dof());
//   scalar_type vm_min, vm_max;
//   for (size_type i=0; i < mf_vm.nb_dof(); ++i) {
//     VM[i] = 0;
//     scalar_type sdiag = 0.;
//     for (unsigned j=0; j < N; ++j) {
//       sdiag += DU[i*N*N + j*N + j];
//       for (unsigned k=0; k < N; ++k) {
// 	scalar_type e = .5*(DU[i*N*N + j*N + k] + DU[i*N*N + k*N + j]);
// 	VM[i] += e*e;
//       }
//     }
//     VM[i] -= 1./N * sdiag * sdiag;
//     vm_min = (i == 0 ? VM[0] : std::min(vm_min, VM[i]));
//     vm_max = (i == 0 ? VM[0] : std::max(vm_max, VM[i]));
//     assert(VM[i] > -1e-6);
//   }
//   cout << "Von Mises : min=" << 4*mu*mu*vm_min << ", max=" << 4*mu*mu*vm_max << "\n";
//   gmm::scale(VM, 4*mu*mu);
//   gmm::vecsave("nonlinear_elastostatic.VM",VM);
// }

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
    default: GMM_ASSERT1(false, "no such law");
  }

  p.resize(pl->nb_params());
  getfem::mdbrick_nonlinear_elasticity<>  ELAS(*pl, mim, mf_u, p);

  getfem::mdbrick_dynamic<> MASS_MATRIX(ELAS, 0.0);

  getfem::mdbrick_abstract<> *pINCOMP = &MASS_MATRIX;
  switch (law_num) {
    case 1:
    case 3: pINCOMP = new getfem::mdbrick_nonlinear_incomp<>(ELAS, mf_p);
  }

  size_type nb_step = size_type(PARAM.int_value("NBSTEP"));
  scalar_type deltat = PARAM.real_value("DELTAT", "Time step");

  // contact condition for Lagrange elements
  dal::bit_vector cn = mf_u.basic_dof_on_region(FRICTION_BOUNDARY_NUM);
//   cout << "cn = " << cn << endl;
//   cout << "cn.card()/N = " << cn.card()/N << endl;
  std::vector<size_type> mouchard(cn.card()/N);
  plain_vector coord_contacts(cn.card());
  int compteur=0;
//  for (int ii=0  ; ii < cn.card(); ++ii) {
//  mouchard[ii]=cn[ii];
//  }
//  gmm::vecsave(datafilename + ".friction_index",mouchard);

  sparse_matrix BN(cn.card()/N, mf_u.nb_dof());
  sparse_matrix BT((N-1)*cn.card()/N, mf_u.nb_dof());
  plain_vector gap(cn.card()/N);
  size_type jj = 0;




  for (dal::bv_visitor i(cn); !i.finished(); ++i)

    if (i % N == 0) {
      cout << "point de contact " << mf_u.point_of_basic_dof(i) << endl;
      coord_contacts[compteur]= mf_u.point_of_basic_dof(i)[0];
      coord_contacts[compteur+1]= mf_u.point_of_basic_dof(i)[1];
      compteur=compteur+2;
      BN(jj, i+N-1) = -1.;
      gap[jj] = mf_u.point_of_basic_dof(i)[N-1];
      mouchard[jj]=i;
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


  getfem::mdbrick_Coulomb_friction<> FRICTION(*pINCOMP, BN, gap,
					      friction_coef, BT);
  FRICTION.set_r(r);
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
  getfem::mdbrick_Dirichlet<> final_model(VOL_F, DIRICHLET_BOUNDARY_NUM);
  final_model.rhs().set(mf_rhs, F2);
  final_model.set_constraints_type(getfem::constraints_type
				   (PARAM.int_value("DIRICHLET_VERSION")));
  // Generic solver.
  getfem::standard_model_state MS(final_model);
  size_type maxit = PARAM.int_value("MAXITER");
  gmm::iteration iter;

  scalar_type dz = PARAM.real_value("DIRICHLET_Z",
				    "Prescribed displacement in z");
  scalar_type dyv = PARAM.real_value("DIRICHLET_Y_SPEED",
				     "Prescribed velocity in y");


  plain_vector U0 = U;
  plain_vector SumN(nb_step), SumT(nb_step),Pressure(nb_step);
//  plain_vector maxT(nb_step);
//  plain_vector maxN(nb_step);
//  plain_vector max(nb_step);
//  double Fmax =0;


 gmm::vecsave(datafilename + ".contact_index",mouchard);
 gmm::vecsave(datafilename + ".contact_coord",coord_contacts);


 for (size_type step = 0; step < nb_step; ++step) {
	cout << "beginning of step " << step+1 << ", number of variables : " << final_model.nb_dof() << endl ;

    for (size_type i = 0; i < nb_dof_rhs; ++i) {
      F2[i*N+N-1] = dz;
     // F2[i*N+N-1] = (step < 10) ? (dz * step/10.0) : dz;
     // F2[i*N+N-1] = dz+dz*3*step/nb_step;
      F2[i*N+N-2] = (step < 10) ? 0.0 : scalar_type(step-10)*deltat*dyv;
    }
    final_model.rhs().set(F2);

    // cout << "F2 = " << F2 << endl;

    FRICTION.set_WT(gmm::scaled(U0, -1.0));

    //Introduction of a mass matrix for one or several steps starting at step DYNAMIC_STEP_NUMBER
    size_type DYNAMIC_STEP_NUMBER = PARAM.int_value("DYNAMIC_STEP_NUMBER") ;

    if (DYNAMIC_STEP_NUMBER > 0  &&  step >= DYNAMIC_STEP_NUMBER) {
      double rho = step > DYNAMIC_STEP_NUMBER ? 0 : 1.0;
      MASS_MATRIX.rho().set(rho/deltat);
    }

    plain_vector MU0(mf_u.nb_dof());
    gmm::mult(MASS_MATRIX.get_M(), U0, MU0);
    VOL_F.set_auxF(MU0);





    iter = gmm::iteration(residual, int(PARAM.int_value("NOISY")),
			  maxit ? maxit : 40000);
    cout << "|U| = " << gmm::vect_norm2(MS.state()) << "\n";



    /* let the default non-linear solve (Newton) do its job */
    getfem::default_newton_line_search ls;
    getfem::standard_solve(MS, final_model, iter,
			   getfem::default_linear_solver(final_model), ls);



    pl->reset_unvalid_flag();
    final_model.compute_residual(MS);
    if (pl->get_unvalid_flag())
      GMM_WARNING1("The solution is not completely valid, the determinant "
		   "of the transformation is negative on "
		   << pl->get_unvalid_flag() << " gauss points");

    gmm::copy(ELAS.get_solution(MS), U);

    gmm::copy(FRICTION.get_LN(MS), LN1);

    double force_totale = 0.;
    int k=0;

    for (k = 0; k < nbc; ++k) {
      force_totale += LN1[k];
	}
      cout << "Force totale = " << force_totale << endl;

    gmm::copy(FRICTION.get_LT(MS), LT1);


    plain_vector LLN_save(nbc), LLT_save(nbc*(N-1));


    {
      gmm::iteration itercg(1e-12, 1);
      plain_vector LLN(nbc), LLT(nbc*(N-1));
      gmm::cg(MMBN, LLN, LN1, gmm::identity_matrix(), itercg);
      itercg.init();
      gmm::cg(MMBT, LLT, LT1, gmm::identity_matrix(), itercg);

      LLN_save = LLN;
      LLT_save = LLT;
    }

    cout << "Pression de contact : " << LN1 << endl;

/*    	cout << "LT1 = " << LT1 << endl;
	cout << "LN1 = " << LN1 << endl;
	*/
	

	plain_vector LT2(nbc);
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
   plain_vector UN(gmm::mat_nrows(BN));
   gmm::mult(BN, U, UN);
//    cout << "UN = " << UN << endl;

    char s[100]; sprintf(s, "step%ld", step+1);
    gmm::vecsave(datafilename + s + ".U",U);


   plain_vector VM(mf_vm.nb_dof());
   ELAS.compute_Von_Mises_or_Tresca(MS, mf_vm, VM, false);
   gmm::vecsave(datafilename + s + ".VM",VM);


   //Sauvegarde des efforts de contact normaux et tangentiels : variables LLN et LLT
   gmm::vecsave(datafilename + s + ".LN",LLN_save);
   gmm::vecsave(datafilename + s + ".LT",LLT_save);

     gmm::copy(U, U0);
        cout << "end of Step nº " << step+1 << " / " << nb_step << endl;


// 	if (max[step]> Fmax)
// 		Fmax = max[step];

}

//  cout << "end of all steps \n" ;
  if (law_num == 3 || law_num == 1) delete pINCOMP;
 /*  plain_vector VM(mf_u.nb_dof());
     cout << "calcul van mises\n";
   calcul_von_mises(mf_u, U0, mf_vm, VM, mu);
     cout << "Fin calcul von mises\n";
 */

  return (iter.converged());

}


/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //  try {
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
    // }  GMM_STANDARD_CATCH_ERROR;

  return 0;

}
