
/**************************************************************************/
/*                                                                        */
/*  Laplacian problem.                                                    */
/*                                                                        */
/**************************************************************************/

#include <getfem_assembling.h>
#include <getfem_export.h>
#include <getfem_norm.h>
#include <getfem_regular_meshes.h>
#include <bgeot_smatrix.h>


using bgeot::base_vector;
using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::size_type;

typedef bgeot::smatrix<scalar_type> sparse_matrix_type;
typedef bgeot::vsvector<scalar_type> linalg_vector;

/**************************************************************************/
/*  exact solution                                                        */
/**************************************************************************/

base_vector sol_K;

scalar_type sol_u(const base_node &x)
{ return sin(bgeot::vect_sp(sol_K, base_vector(x))); }

scalar_type sol_f(const base_node &x)
{ return bgeot::vect_sp(sol_K, sol_K) * sin(bgeot::vect_sp(sol_K, x)); }

base_vector sol_grad(const base_node &x)
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
  int NX, N, K, fem_type;

  sparse_matrix_type RM;   /* rigidity matrix.                            */
  linalg_vector U, B; /* inconnue et second membre.                       */
 
  int integration, mesh_type;

  std::string datafilename;
  ftool::md_param PARAM;

  void assemble(void);
  void solve(void);
  void init(void);
  lap_pb(void) : mef(mesh), mef_data(mesh), mef_data2(mesh) {}
};

void lap_pb::init(void)
{
  dal::bit_vector nn;
  size_type i, j, k;

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
  fem_type = PARAM.int_value("FEM_TYPE", "Finite element method");
  datafilename = std::string( PARAM.string_value("ROOTFILENAME",
			     "File name for saving"));

  scalar_type FT = PARAM.real_value("FT", "parameter for exact solution");

  sol_K = base_vector(N);
  for (j = 0; j < N; j++)
    sol_K[j] = ((j & 1) == 0) ? FT : -FT;

  /***********************************************************************/
  /*  BUILD MESH.                                                        */
  /***********************************************************************/

  cout << "Mesh generation\n";

  base_node org(N); org.fill(0.0);
  std::vector<base_vector> vtab(N);
  std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
  for (i = 0; i < N; i++)
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
  default : DAL_THROW(dal::internal_error, "Unknown finite element method");
  }

  getfem::pintegration_method ppi;
  nn = mesh.convex_index(N);
  switch (integration) {
  case 0 :
    switch (mesh_type) { 
    case 0 : ppi = bgeot::simplex_poly_integration(N); break;
    case 1 : ppi = bgeot::parallelepiped_poly_integration(N); break;
    default : DAL_THROW(dal::internal_error, 
    "Exact integration not allowed in this context");
    }
    break;
  case 1 :
    switch (mesh_type) { 
    case 0 : 
      ppi = bgeot::Newton_Cotes_approx_integration(N,2*K);
      break;
    case 1 : 
      ppi = bgeot::parallelepiped_Newton_Cotes_approx_integration(N, 2*K);
      break;
    case 2 :
      ppi = bgeot::prism_Newton_Cotes_approx_integration(N, 2*K);
      break;
    }
    break;
  case 2 :
    if (mesh_type == 1)
      ppi = bgeot::parallelepiped_Gauss_approx_integration(N, K+1);
    else
      DAL_THROW(dal::internal_error,
		"Product of 1D Gauss only for parallelepipeds");
    break;
  case 11 : ppi = bgeot::triangle1_approx_integration(); break;
  case 12 : ppi = bgeot::triangle2_approx_integration(); break;
  case 13 : ppi = bgeot::triangle3_approx_integration(); break;
  case 14 : ppi = bgeot::triangle4_approx_integration(); break;
  case 15 : ppi = bgeot::triangle5_approx_integration(); break;
  case 16 : ppi = bgeot::triangle6_approx_integration(); break;
  case 17 : ppi = bgeot::triangle7_approx_integration(); break;
  case 21 : ppi = bgeot::tetrahedron1_approx_integration(); break;
  case 22 : ppi = bgeot::tetrahedron2_approx_integration(); break;
  case 23 : ppi = bgeot::tetrahedron3_approx_integration(); break;
  case 25 : ppi = bgeot::tetrahedron5_approx_integration(); break;
  case 32 : ppi = bgeot::quad2_approx_integration(); break;
  case 33 : ppi = bgeot::quad3_approx_integration(); break;
  case 35 : ppi = bgeot::quad5_approx_integration(); break;
  default : DAL_THROW(std::logic_error, "Undefined integration method");
  }

  switch (mesh_type) {
  case 0 : 
    mef.set_finite_element(nn, getfem::PK_fem(N, K), ppi);
    mef_data.set_finite_element(nn, getfem::PK_fem(N, K),
				bgeot::simplex_poly_integration(N));
    mef_data2.set_finite_element(nn, getfem::PK_fem(N, 0),
				 bgeot::simplex_poly_integration(N));
    break;
  case 1 : 
    mef.set_finite_element(nn, getfem::QK_fem(N, K), ppi); 
    mef_data.set_finite_element(nn, getfem::QK_fem(N, K), ppi);
    mef_data2.set_finite_element(nn, getfem::QK_fem(N, 0),  ppi);
    break;
  case 2 : 
    mef.set_finite_element(nn, getfem::product_fem(getfem::PK_fem(N-1, K),
						   getfem::PK_fem(1, K)), ppi);
    mef_data.set_finite_element(nn,
				getfem::product_fem(getfem::PK_fem(N-1, K),
						 getfem::PK_fem(1, K)), ppi);
    mef_data2.set_finite_element(nn, 
				 getfem::product_fem(getfem::PK_fem(N-1, 0),
						  getfem::PK_fem(1, 0)), ppi);
    break;
  }

  switch(fem_type) {

  case 0 : break;

  case 1 :
    mef.set_finite_element(nn, getfem::segment_Hermite_fem(), ppi);
    break;
  
  }
  
  
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  nn = mesh.convex_index(N);
  base_vector un;
  for (j << nn; j != size_type(-1); j << nn) {
    k = mesh.structure_of_convex(j)->nb_faces();
    for (i = 0; i < k; i++) {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty()) {
        un = mesh.normal_of_face_of_convex(j, i, 0);
	un /= bgeot::vect_norm2(un);

	// if (true) {
	if (dal::abs(un[N-1] - 1.0) < 1.0E-7) {
	  mef.add_boundary_elt(0, j, i);
	  // cout << "ajout bord Dirichlet, cv\t" << j << "\tf " << i << endl;
	}
	else {
	  mef.add_boundary_elt(1, j, i);
	  // cout << "ajout bord Neumann, cv\t" << j << "\tf " << i << endl;
	}
      }
    }
  }
}

void lap_pb::assemble(void)
{
  int nb_dof = mef.nb_dof(), nb_dof_data = mef_data.nb_dof();
  int nb_dof_data2 = mef_data2.nb_dof();
  B = linalg_vector(nb_dof); B.fill(0.0);
  U = linalg_vector(nb_dof); U.fill(0.0); 
  RM = sparse_matrix_type(nb_dof, nb_dof);
  linalg_vector ST;
  
  cout << "Number of dof : " << nb_dof << endl;
  // cout << "dof number 1 : " << nb_dof_data << endl;
  // cout << "dof number 2 : " << nb_dof_data2 << endl;

  cout << "Assembling rigidity matrix" << endl;
  ST = linalg_vector(nb_dof_data2);
  std::fill(ST.begin(), ST.end(), 1.0);
  getfem::assembling_rigidity_matrix_for_laplacian(RM, mef, mef_data2, ST);
  
  cout << "Assembling source term" << endl;
  ST = linalg_vector(nb_dof_data);
  for (size_type i = 0; i < nb_dof_data; ++i)
    ST[i] = sol_f(mef_data.point_of_dof(i));
  getfem::assembling_volumic_source_term(B, mef, mef_data, ST, 1);

  cout << "Assembling Neumann condition" << endl;
  ST = linalg_vector(nb_dof_data);
  getfem::base_node pt(N);

  dal::bit_vector nn = mesh.convex_index(N);
  base_vector un;
  size_type j;
  for (j << nn; j != size_type(-1); j << nn)
  { // Ne tiens pas compte des coins ...
    size_type k = mesh.structure_of_convex(j)->nb_faces();
    getfem::pfem pf = mef_data.fem_of_element(j);
    for (size_type i = 0; i < k; ++i) {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty()) {
	
	for (size_type l = 0; l < pf->structure()->nb_points_of_face(i); ++l) {
	  size_type n = pf->structure()->ind_points_of_face(i)[l];
	  un = mesh.normal_of_face_of_convex(j, i, pf->node_of_dof(n));
	  un /= bgeot::vect_norm2(un);
	  size_type dof = mef_data.ind_dof_of_element(j)[n];
	  ST[dof] = bgeot::vect_sp(sol_grad(mef_data.point_of_dof(dof)), un);
	}
      }
    }
  }
  getfem::assembling_Neumann_condition(B, mef, 1, mef_data, ST, 1);
  
  cout << "take Dirichlet condition into account" << endl;
  // nn = mef.dof_on_boundary(0);
  // cout << "Number of Dirichlet nodes : " << nn.card() << endl;
  // cout << "Dirichlet nodes : " << nn << endl;
  ST = linalg_vector(nb_dof);
  for (size_type i = 0; i < nb_dof; ++i)
    ST[i] = sol_u(mef.point_of_dof(i));
  getfem::assembling_Dirichlet_condition(RM, B, mef, 0, ST, 1);
}

void lap_pb::solve(void)
{
  bgeot::cg(RM, U, B, 20000, residu, false);
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[])
{
  try
    {
    
    lap_pb p;
    scalar_type exectime = ftool::uclock_sec(), total_time = 0.0;
    
    // cout << "initialisation ...\n";
    p.PARAM.read_command_line(argc, argv);
    p.init();
    // cout << "Initialisation terminee\n";
    
    std::ofstream cres((p.datafilename + ".res").c_str());
    cres << p.N << "\t" <<  p.K << "\t" << p.NX << "\t";
    cres << ftool::uclock_sec() - exectime << "  ";
    
    total_time += ftool::uclock_sec() - exectime;
    
    // p.mesh.write_to_file(cout);
    // p.mesh.stat();
    
    p.mesh.write_to_file(p.datafilename + ".mesh" + char(0));
    
    exectime = ftool::uclock_sec();
    int nb_dof = p.mef.nb_dof();
    
    total_time += ftool::uclock_sec() - exectime;
    
    cres << nb_dof << "\t" <<  ftool::uclock_sec() - exectime << "\t";
    
    cout << "Assembling \n";
    exectime = ftool::uclock_sec();
    p.assemble();
    
    cres << ftool::uclock_sec() - exectime << "\t";
    total_time += ftool::uclock_sec() - exectime;
    
    //   cout << "Matrice de rigidite\n";
    //   for (i = 0; i < p.RM.nrows(); i++)
    //   { 
    //     cout << "ligne " << i << " [ ";
    //     for (l = 0; l < p.RM.nrows(); l++)
    //       if (p.RM(i, l) != 0.0)
    // 	cout << "(" << l << "," << p.RM(i, l) << ")  ";
    //     cout << "]" << endl;
    //   }
    //   cout << endl << endl;
    
    cout << "Solving the system\n";
    exectime = ftool::uclock_sec();
    p.solve();
    
    cres << ftool::uclock_sec() - exectime << "\t";
    total_time += ftool::uclock_sec() - exectime;
    exectime = ftool::uclock_sec();
    
    int nbdof = p.mef_data.nb_dof();
    linalg_vector V(nbdof);
    interpolation_solution_same_mesh(p.mef, p.mef_data, p.U, V, 1);
    for (int i = 0; i < nbdof; ++i)
      V[i] -= sol_u(p.mef_data.point_of_dof(i));
    
    scalar_type l2norm = getfem::L2_norm(p.mef_data, V, 1);
    cres << l2norm << "\t";
    cres << ftool::uclock_sec() - exectime << "\t";
    total_time += ftool::uclock_sec() - exectime;
    exectime = ftool::uclock_sec();
    scalar_type h1norm = getfem::H1_norm(p.mef_data, V, 1);
    cres << h1norm << "\t";
    cres << ftool::uclock_sec() - exectime << "\t";
    total_time += ftool::uclock_sec() - exectime;
    cres << total_time << endl;

    cout.precision(16);
    cout << "L2 error = " << l2norm << endl
	 << "H1 error = " << h1norm << endl;

    getfem::save_solution(p.datafilename + ".dataelt", p.mef, p.U, 1, p.K);
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
}
