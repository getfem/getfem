/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
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
/* *********************************************************************** */
/*                                                                         */
/*  Laplacian (Poisson) problem.                                           */
/*                                                                         */
/* *********************************************************************** */

#include <getfem_assembling.h>
#include <getfem_export.h>
#include <getfem_norm.h>
#include <getfem_regular_meshes.h>
#include <gmm.h>


using bgeot::base_vector;
using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::size_type;
using bgeot::dim_type;

typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;
typedef std::vector<scalar_type> linalg_vector;

/**************************************************************************/
/*  exact solution                                                        */
/**************************************************************************/

base_vector sol_K;

scalar_type sol_u(const base_node &x)
{ return sin(bgeot::vect_sp(sol_K, base_vector(x))); }

scalar_type sol_f(const base_vector &x)
{ return bgeot::vect_sp(sol_K, sol_K) * sin(bgeot::vect_sp(sol_K, x)); }

base_vector sol_grad(const base_vector &x)
{
  base_vector res = sol_K;
  res *= cos(bgeot::vect_sp(sol_K, x));
  return res;
}

/**************************************************************************/
/*  structure representing the problem.                                   */
/**************************************************************************/

struct lap_pb
{
  getfem::getfem_mesh mesh;
  getfem::mesh_fem mef;
  getfem::mesh_fem mef_data;
  getfem::mesh_fem mef_data2;

  scalar_type LX, LY, LZ, incline, residu;
  size_type N;
  int NX, K, fem_type, KI, gen_dirichlet;
  size_type nb_dof;

  sparse_matrix_type SM;   /* stiffness matrix.                           */
  linalg_vector U, B; /* inconnue et second membre.                       */

  linalg_vector Ud;  /* reduced solution for generic Dirichlet condition. */
  col_sparse_matrix_type NN;

 
  int integration, mesh_type;

  std::string datafilename;
  ftool::md_param PARAM;
  bool failed;

  void assemble(void);
  void solve(void);
  void init(void);
  lap_pb(void) : mef(mesh), mef_data(mesh), mef_data2(mesh) {}
};

void lap_pb::init(void)
{
  dal::bit_vector nn;

  /***********************************************************************/
  /* READING PARAMETER FILE                                              */
  /***********************************************************************/
  
  N = PARAM.int_value("N", "Domaine dimension");
  LX = PARAM.real_value("LX", "Size in X");
  LY = PARAM.real_value("LY", "Size in Y");
  LZ = PARAM.real_value("LZ", "Size in Y");
  incline = PARAM.real_value("INCLINE", "incline of the mesh");
  NX = PARAM.int_value("NX", "Nomber of sace steps ");
  integration = PARAM.int_value("INTEGRATION", "integration method");
  mesh_type = PARAM.int_value("MESH_TYPE", "Mesh type ");
  residu = PARAM.real_value("RESIDU", "Residu for c.g.");
  K = PARAM.int_value("K", "Finite element degree");
  KI = PARAM.int_value("KI", "Integration degree");
  gen_dirichlet = PARAM.int_value("GENERIC_DIRICHLET",
				  "Generic Dirichlet condtion");
  fem_type = PARAM.int_value("FEM_TYPE", "Finite element method");
  datafilename = std::string( PARAM.string_value("ROOTFILENAME",
			     "File name for saving"));

  scalar_type FT = PARAM.real_value("FT", "parameter for exact solution");

  sol_K = base_vector(N);
  for (dim_type j = 0; j < N; j++)
    sol_K[j] = ((j & 1) == 0) ? FT : -FT;

  /***********************************************************************/
  /*  BUILD MESH.                                                        */
  /***********************************************************************/

  cout << "Mesh generation\n";

  base_node org(N); org.fill(0.0);
  std::vector<base_vector> vtab(N);
  std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
  for (dim_type i = 0; i < N; i++)
  { 
    vtab[i] = base_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
  }
  if (N > 1) vtab[N-1][0] = incline * LX / scalar_type(NX);

  switch (mesh_type) {
  case 0 : getfem::parallelepiped_regular_simplex_mesh
		   (mesh, N, org, vtab.begin(), ref.begin()); break;
  case 1 : getfem::parallelepiped_regular_mesh
		   (mesh, N, org, vtab.begin(), ref.begin()); break;
  case 2 : getfem::parallelepiped_regular_prism_mesh
		   (mesh, N, org, vtab.begin(), ref.begin()); break;
  default : DAL_THROW(dal::internal_error, "Unknown type of mesh");
  }

  mesh.optimize_structure();

  if (mesh_type == 2 && N <= 1) mesh_type = 0;

  cout << "Selecting finite element method.\n";

  switch(fem_type) {
  case 0 : break;
  case 1 :
    if (N != 1 || mesh_type != 0)
      DAL_THROW(dal::internal_error,
		"This element is only defined on segments");
     K = 3;
    break;
  case 2 : 
    if (mesh_type != 0)
      DAL_THROW(dal::internal_error,
		"This element is only defined on simplexes");
    break;
  case 3 : 
    if (mesh_type != 0)
      DAL_THROW(dal::internal_error,
		"This element is only defined on simplexes");
    break;
  default : DAL_THROW(dal::internal_error, "Unknown finite element method");
  }

  getfem::pintegration_method ppi;
  char meth[500];
  nn = mesh.convex_index(N);
  switch (integration) {
  case 0 :
    switch (mesh_type) { 
    case 0 : sprintf(meth, "IM_EXACT_SIMPLEX(%d)", int(N)); break;
    case 1 : sprintf(meth, "IM_EXACT_PARALLELEPIPED(%d)", int(N)); break;
    default : DAL_THROW(dal::internal_error, 
    "Exact integration not allowed in this context");
    }
    break;
  case 1 :
    switch (mesh_type) { 
    case 0 : 
      sprintf(meth, "IM_NC(%d,%d)", int(N), int(2*K));
      break;
    case 1 : 
      sprintf(meth, "IM_NC_PARALLELEPIPED(%d,%d)", int(N), int(2*K));
      break;
    case 2 :
      sprintf(meth, "IM_NC_PRISM(%d,%d)", int(N), int(2*K));
      break;
    }
    break;
  case 2 :
    if (mesh_type == 1 || N == 1)
      sprintf(meth, "IM_GAUSS_PARALLELEPIPED(%d,%d)", int(N), int(2*K));
    else
      DAL_THROW(dal::internal_error,
		"Product of 1D Gauss only for parallelepipeds");
    break;
  case 3 :
    if (mesh_type == 0) {
      if (N == 1)
	sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_GAUSS1D(%d), %d)",2,int(KI));
      else if (N == 2)
	sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(%d), %d)",
		2,int(KI));
      else
	sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_NC(%d, %d), %d)",
		int(N), 2, int(KI));
    }
    else
      DAL_THROW(dal::internal_error,
		"Composite integration only for simplexes");
    break;
  case 11 : sprintf(meth, "IM_TRIANGLE(1)"); break;
  case 12 : sprintf(meth, "IM_TRIANGLE(2)"); break;
  case 13 : sprintf(meth, "IM_TRIANGLE(3)"); break;
  case 14 : sprintf(meth, "IM_TRIANGLE(4)"); break;
  case 15 : sprintf(meth, "IM_TRIANGLE(5)"); break;
  case 16 : sprintf(meth, "IM_TRIANGLE(6)"); break;
  case 17 : sprintf(meth, "IM_TRIANGLE(7)"); break;
  case 21 : sprintf(meth, "IM_TETRAHEDRON(1)"); break;
  case 22 : sprintf(meth, "IM_TETRAHEDRON(2)"); break;
  case 23 : sprintf(meth, "IM_TETRAHEDRON(3)"); break;
  case 25 : sprintf(meth, "IM_TETRAHEDRON(5)"); break;
  case 32 : sprintf(meth, "IM_QUAD(2)"); break;
  case 33 : sprintf(meth, "IM_QUAD(3)"); break;
  case 35 : sprintf(meth, "IM_QUAD(5)"); break;
  default : DAL_THROW(std::logic_error, "Undefined integration method");
  }
  ppi = getfem::int_method_descriptor(meth);
  getfem::pfem pfprinc = 0;
  switch (mesh_type) {
  case 0 :
    sprintf(meth, "FEM_PK(%d,%d)", int(N), int(K));
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    mef_data.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    sprintf(meth, "FEM_PK(%d,%d)", int(N), 0);
    mef_data2.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    break;
  case 1 :
    sprintf(meth, "FEM_QK(%d,%d)", int(N), K);
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, getfem::fem_descriptor(meth), ppi); 
    mef_data.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    sprintf(meth, "FEM_QK(%d,%d)", int(N), 0);
    mef_data2.set_finite_element(nn, getfem::fem_descriptor(meth),  ppi);
    break;
  case 2 :
    sprintf(meth, "FEM_PK_PRISM(%d,%d)", int(N), K);
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    mef_data.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    sprintf(meth, "FEM_PK_PRISM(%d,%d)", int(N), 0);
    mef_data2.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    break;
  }

  switch(fem_type) {

  case 0 : break;

  case 1 :
    sprintf(meth, "FEM_HERMITE_SEGMENT");
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    break;
    
  case 2 :
    sprintf(meth, "FEM_PK_HIERARCHICAL(%d, %d)", int(N), int(K));
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    break;

  case 3 :
    sprintf(meth, "FEM_PK_HIERARCHICAL_COMPOSITE(%d,%d,%d)", int(N), 1, int(K));
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
    break;
  
  }
  
  cout << "Name of principal finite element method : "
       << getfem::name_of_fem(pfprinc) << endl;
  cout << "Name of principal integration method : "
       << getfem::name_of_int_method(ppi) << endl;
//   if (pfprinc->is_polynomial()) {
//     cout << "basis of the principal finite element method : " << endl;
//     for (size_type l = 0; l < pfprinc->nb_dof(); ++l) {
//       cout << "base " << l << " : " 
//            << (((getfem::ppolyfem)(pfprinc))->base())[l] << endl;
//     }
//   }
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  nn = mesh.convex_index(N);
  base_vector un;
  size_type j;
  for (j << nn; j != size_type(-1); j << nn) {
    size_type k = mesh.structure_of_convex(j)->nb_faces();
    for (size_type i = 0; i < k; i++) {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty()) {
        un = mesh.normal_of_face_of_convex(j, i, 0);
	un /= bgeot::vect_norm2(un);

	// if (false) {
	if (dal::abs(un[N-1] - 1.0) < 1.0E-7) {
	  mef.add_boundary_elt(0, j, i);
	  // cout << "ajout bord Dirichlet, cv\t" << j << "\tf " << i << endl;
	}
	else {
	  size_type nb;
	  for(nb = 0; nb < N; ++nb) {
	    if (dal::abs(un[nb] + 1.0) < 1.0E-7) { nb = 1 + nb * 2; break; }
	    if (dal::abs(un[nb] - 1.0) < 1.0E-7) { nb = 2 + nb * 2; break; }
	  }
	  mef.add_boundary_elt(nb, j, i);
	  // cout << "ajout bord Neumann, cv\t" << j << "\tf " << i << endl;
	}
      }
    }
  }
}

void lap_pb::assemble(void)
{
  nb_dof = mef.nb_dof();
  size_type nb_dof_data = mef_data.nb_dof();
  size_type nb_dof_data2 = mef_data2.nb_dof();
  B = linalg_vector(nb_dof); gmm::clear(B);
  U = linalg_vector(nb_dof); gmm::clear(U); 
  SM = sparse_matrix_type(nb_dof, nb_dof);
  linalg_vector ST;
  
  cout << "Number of dof : " << nb_dof << endl;
  // cout << "dof number 1 : " << nb_dof_data << endl;
  // cout << "dof number 2 : " << nb_dof_data2 << endl;

  cout << "Assembling stiffness matrix" << endl;
  ST = linalg_vector(nb_dof_data2);
  std::fill(ST.begin(), ST.end(), 1.0);
  getfem::asm_stiffness_matrix_for_laplacian(SM, mef, mef_data2, ST);
  
  cout << "Assembling source term" << endl;
  ST = linalg_vector(nb_dof_data);
  for (size_type i = 0; i < nb_dof_data; ++i)
    ST[i] = sol_f(mef_data.point_of_dof(i));
  getfem::asm_source_term(B, mef, mef_data, ST);

  cout << "Assembling Neumann condition" << endl;
  ST = linalg_vector(nb_dof_data);
  getfem::base_node pt(N);

  for (size_type nb = 1; nb <= 2*N; ++nb) {
    dal::bit_vector nn = mesh.convex_index(N);
    base_vector un;
    size_type j;
    for (j << nn; j != size_type(-1); j << nn) {
      getfem::pfem pf = mef_data.fem_of_element(j);
      dal::bit_vector nf = mef.faces_of_convex_on_boundary(j, nb);
      size_type i;
      for (i << nf; i != size_type(-1); i << nf) {
	for (size_type l = 0; l<pf->structure()->nb_points_of_face(i); ++l) {
	  size_type n = pf->structure()->ind_points_of_face(i)[l];
	  un = mesh.normal_of_face_of_convex(j, i, pf->node_of_dof(n));
	  un /= bgeot::vect_norm2(un);
	  size_type dof = mef_data.ind_dof_of_element(j)[n];
	  ST[dof] = bgeot::vect_sp(sol_grad(mef_data.point_of_dof(dof)), un);
	}
      }
    }
    getfem::asm_source_term(B, mef, mef_data, ST, nb);
  }
  
  cout << "take Dirichlet condition into account" << endl;
  
  if (!gen_dirichlet) {
    ST = linalg_vector(nb_dof);
    for (size_type i = 0; i < nb_dof; ++i)
      ST[i] = sol_u(mef.point_of_dof(i));
    getfem::assembling_Dirichlet_condition(SM, B, mef, 0, ST);
  }
  else {

    ST = linalg_vector(nb_dof_data);
    for (size_type i = 0; i < nb_dof_data; ++i)
      ST[i] = sol_u(mef_data.point_of_dof(i));
    
    Ud = linalg_vector(nb_dof);
    NN = col_sparse_matrix_type(nb_dof, nb_dof);
    col_sparse_matrix_type HH(nb_dof, nb_dof);
    linalg_vector RR(nb_dof), RHaux(nb_dof);

    getfem::asm_dirichlet_constraints(HH, RR, mef, mef_data, ST, 0);
    
    // cout << "HH = " << HH  << endl;
    gmm::clean(HH, 1e-15);

    int nbcols = getfem::Dirichlet_nullspace(HH, NN, RR, Ud);
    // cerr << "Number of unknowns : " << nbcols << endl;
    NN.resize(nbcols);

    gmm::mult(SM, Ud, gmm::scaled(B, -1.0), RHaux);
    B = linalg_vector(nbcols);
    U = linalg_vector(nbcols);
    gmm::mult(gmm::transposed(NN), gmm::scaled(RHaux, -1.0), B);
    sparse_matrix_type SMaux(nbcols, nb_dof);
    gmm::mult(gmm::transposed(NN), SM, SMaux);
    SM = sparse_matrix_type(nbcols, nbcols);
    sparse_matrix_type NNaux(nb_dof, nbcols);
    gmm::copy(NN, NNaux);
    gmm::mult(SMaux, NNaux, SM);

  }
}

void lap_pb::solve(void) {

  cout << "Compute preconditionner\n";
  double time = ftool::uclock_sec();
  gmm::iteration iter(residu, 1, 40000);
  gmm::identity_matrix P;
  // gmm::diagonal_precond<sparse_matrix_type> P(SM);
  // gmm::mr_approx_inverse_precond<sparse_matrix_type> P(SM, 10, 10E-17);
  // gmm::cholesky_precond<sparse_matrix_type> P(SM);
  // gmm::ilu_precond<sparse_matrix_type> P(SM);
  // gmm::ilut_precond<sparse_matrix_type> P(SM, 5, 1E-6);
  // gmm::choleskyt_precond<sparse_matrix_type> P(SM, 5, 1E-6);

  cout << "Time to compute preconditionner : "
       << ftool::uclock_sec() - time << " seconds\n";

  cout << "Solve\n";
  // gmm::gmres(SM, U, B, P, 50, iter);
  gmm::cg(SM, U, B, P, iter);

  if (gen_dirichlet) {
    linalg_vector Uaux(nb_dof);
    gmm::mult(NN, U, Ud, Uaux);
    U = linalg_vector(nb_dof);
    gmm::copy(Uaux, U);
  }

  failed = !(iter.converged());
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[])
{
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb);

  try {
    
    lap_pb p;
    scalar_type exectime = ftool::uclock_sec(), total_time = 0.0;
    
    cout << "initialisation ...\n";
    p.PARAM.read_command_line(argc, argv);
    p.init();
    
    total_time += ftool::uclock_sec() - exectime;
    
    p.mesh.write_to_file(p.datafilename + ".mesh" + char(0));
    
    exectime = ftool::uclock_sec();
    int nb_dof = p.mef.nb_dof();
    cout << "Nb dofs : " << nb_dof << endl;
    
    total_time += ftool::uclock_sec() - exectime;
    
    cout << "Assembling \n";
    exectime = ftool::uclock_sec();
    p.assemble();
    
    total_time += ftool::uclock_sec() - exectime;
    
    //   cout << "Stifness matrix\n" << p.SM << endl;
    
    exectime = ftool::uclock_sec();
    p.solve();
    cout << "solve time : " << ftool::uclock_sec() - exectime << "seconds" << endl;
    if (p.failed) { cerr << "Solve procedure has failed\n"; }
    
    total_time += ftool::uclock_sec() - exectime;
    exectime = ftool::uclock_sec();
    
    size_type nbdof = p.mef_data.nb_dof();
    linalg_vector V(nbdof), W(nbdof);
    scalar_type linfnorm = 0.0;
    getfem::interpolation_solution_same_mesh(p.mef, p.mef_data, p.U, V, 1);
    for (size_type i = 0; i < nbdof; ++i) {
      W[i] = sol_u(p.mef_data.point_of_dof(i)); V[i] -= W[i];
      linfnorm = std::max(linfnorm, dal::abs(V[i]));
    }
    
    scalar_type l2norm = getfem::asm_L2_norm(p.mef_data, V);
    total_time += ftool::uclock_sec() - exectime;
    exectime = ftool::uclock_sec();
    scalar_type h1norm = getfem::asm_H1_norm(p.mef_data, V);
    total_time += ftool::uclock_sec() - exectime;

    cout.precision(40);
    cout << "L2 error = " << l2norm << endl
	 << "H1 error = " << h1norm << endl
	 << "Linfty error = " << linfnorm << endl;
     
    getfem::save_solution(p.datafilename + ".dataelt", p.mef, p.U, p.K);
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
}
