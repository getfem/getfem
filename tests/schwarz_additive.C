/**************************************************************************/
/*                                                                        */
/*  Schwarz additive test program on an elastostatic problem with         */
/*  an optional coarse mesh.                                              */
/*                                                                        */
/**************************************************************************/

#include <getfem_assembling.h>
#include <getfem_norm.h>
#include <getfem_regular_meshes.h>
#include <getfem_export.h>
#include <gmm.h>

using bgeot::base_vector;
using bgeot::base_node;
using bgeot::size_type;
using bgeot::scalar_type;

typedef gmm::row_matrix<gmm::rsvector<scalar_type> > general_sparse_matrix;
typedef std::vector<scalar_type>                     linalg_vector;


struct pb_data {
  getfem::getfem_mesh mesh;
  getfem::getfem_mesh mesh_coarse;

  getfem::mesh_fem mef;
  getfem::mesh_fem mef_data;
  getfem::mesh_fem mef_coarse;

  double mu, lambda, rho, gravity;
  double LX, LY, LZ, residu, overlap, subdomsize;
  int NX, N, NXCOARSE, USECOARSE, K;
  base_vector D;

  general_sparse_matrix RM;   /* stifness matrix.                         */
  linalg_vector U, F;         /* Unknown and right hand side.             */
  int solver;

  void assemble(void);
  void init(ftool::md_param &params);

  int solve_cg(void);
  int solve_cg2(void);
  int solve_schwarz(int);

  int solve(void) {
    cout << "solving" << endl;
    switch (solver) {
    case 0 : return solve_cg();
    case 1 : case 2 : return solve_schwarz(solver);
    case 3 : return solve_cg2();
    }
    return 0;
  }

  base_vector vol_force(const base_node &x)
  { base_vector res(x.size()); res[N-1] = -rho*gravity; return res; }

  pb_data(void) : mef(mesh), mef_data(mesh), mef_coarse(mesh_coarse)  {}
};

ftool::md_param PBSTFR_PARAM;

void pb_data::init(ftool::md_param &params) {

  /***********************************************************************/
  /*  READING PARAMETER FILE.                                            */
  /***********************************************************************/
  
  /* parametres physiques */
  N = params.int_value("N", "Dimension");
  mu = params.real_value("MU", "Stiffness parameter mu");
  gravity = params.real_value("PG", "G");
  rho = params.real_value("RHO", "RHO");
  lambda = params.real_value("LAMBDA", "lambda");
  D.resize(N); D.fill(0.0);
  D[N-1] = params.real_value("D", "Dirichlet condition");
  
  /* parametres numeriques */
  LX = params.real_value("LX", "Size in X");
  LY = params.real_value("LY", "Size in Y");
  LZ = params.real_value("LZ", "Size in Y");
  NX = params.int_value("NX", "Nomber of space step ");
  NXCOARSE = params.int_value("NXCOARSE", "Nombre of space step ");
  USECOARSE = params.int_value("USECOARSE", "Coarser mesh or not");
  residu = params.real_value("RESIDU", "residu");
  overlap = params.real_value("OVERLAP", "overlap");
  K = params.int_value("K", "Degree");
  solver = params.int_value("SOLVER", "solver");
  subdomsize = params.real_value("SUBDOMSIZE", "sub-domains size");  
  std::string meshname(params.string_value("MESHNAME",
			     "mesh file name"));
  std::cout << "\n\n";

  /***********************************************************************/
  /*  BUILDING MESH.                                                     */
  /***********************************************************************/

  std::cout << "building mesh\n";

  if (meshname.size() > 0) {
    mesh.read_from_file(meshname);
  }
  else {
    base_node org(N); org.fill(0.0);
    std::vector<bgeot::base_small_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
    for (int i = 0; i < N; i++) { 
      vtab[i] = bgeot::base_small_vector(N); vtab[i].fill(0.0);
      (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
    }
    getfem::parallelepiped_regular_simplex_mesh(mesh, N, org,
						&(vtab[0]), &(ref[0]));
  }

  if (USECOARSE) { // coarse mesh
    base_node org(N); org.fill(0.0);
    std::vector<bgeot::base_small_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NXCOARSE);
    for (int i = 0; i < N; i++) { 
      vtab[i] = bgeot::base_small_vector(N); vtab[i].fill(0.0);
      (vtab[i])[i] = 
	((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NXCOARSE);
    }
    getfem::parallelepiped_regular_simplex_mesh(mesh_coarse, N, org,
						&(vtab[0]), &(ref[0]));
  }
  
  mesh.trans_of_convex(0);
  mesh.optimize_structure();

  dal::bit_vector nn = mesh.convex_index(N);
  char method[500];
  sprintf(method, "IM_EXACT_SIMPLEX(%d)", N);
  getfem::pintegration_method ppi = getfem::int_method_descriptor(method);
  
  sprintf(method, "FEM_PK(%d, %d)", N, K);
  mef.set_finite_element(nn, getfem::fem_descriptor(method), ppi);
  mef_coarse.set_finite_element(mesh_coarse.convex_index(N),
				getfem::fem_descriptor(method), ppi);
  mef_data.set_finite_element(nn, getfem::fem_descriptor(method), ppi);
  mef.set_qdim(N);
  mef_coarse.set_qdim(N);

  nn = mesh.convex_index(N);
  base_vector un(N);
  for (int j = nn.take_first(); j >= 0; j << nn) {
    int k = mesh.structure_of_convex(j)->nb_faces();
    for (int i = 0; i < k; i++) {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty()) {
	gmm::copy(mesh.normal_of_face_of_convex(j, i, 0), un);
	un /= bgeot::vect_norm2(un);	
	if (dal::abs(un[N-1] - 1.0) < 1.0E-3) mef.add_boundary_elt(0, j, i);
      }
    }
  }
}

void pb_data::assemble(void) {
  size_type nb_dof = mef.nb_dof();
  std::cout << "number of dof : "<< nb_dof << endl;
  size_type nb_dof_data = mef_data.nb_dof();
  dal::bit_vector ddlD = mef.dof_on_boundary(0);
 
  F.resize(nb_dof); gmm::clear(F);
  U.resize(nb_dof); gmm::clear(U);
  gmm::resize(RM, nb_dof, nb_dof);

  std::cout << "Assembly of stiffness matrix" << endl;

  linalg_vector STLA(nb_dof_data), STG(nb_dof_data);
  std::fill(STLA.begin(), STLA.end(), lambda);
  std::fill(STG.begin(), STG.end(), mu);
  getfem::asm_stiffness_matrix_for_linear_elasticity(RM,mef,mef_data,STLA,STG);

  std::cout << "Assembly of source term" << endl;
  linalg_vector STF(N * nb_dof_data);
  for (size_type j = 0; j < nb_dof_data; j++)
    for (int k = 0; k < N; k++)
      STF[j*N + k] = (vol_force(mef_data.point_of_dof(j)))[k];
  getfem::asm_source_term(F, mef, mef_data, STF);
  
  linalg_vector UD(nb_dof);
  for (size_type j = 0; j < nb_dof/N; j++)
    for (size_type k = 0; k < size_type(N); k++) UD[j*N + k] = D[k];
  getfem::assembling_Dirichlet_condition(RM, F, mef, 0, UD);
}

int pb_data::solve_cg(void) {
  gmm::iteration iter(residu, 1, 1000000);
  gmm::ildlt_precond<general_sparse_matrix> P(RM);
  gmm::cg(RM, U, F, gmm::identity_matrix(), P, iter);
  return iter.get_iteration();
}

int pb_data::solve_cg2(void) {
  gmm::iteration iter(residu, 1, 1000000);
  gmm::cg(RM, U, F, gmm::identity_matrix(), gmm::identity_matrix(), iter);
  return iter.get_iteration();
}

int pb_data::solve_schwarz(int version) {

  size_type nb_dof = mef.nb_dof();
  std::vector<base_node> pts(nb_dof);
  for (size_type i = 0; i < nb_dof; ++i) pts[i] = mef.point_of_dof(i);

  std::vector<general_sparse_matrix> vB;
  gmm::rudimentary_regular_decomposition(pts, subdomsize, overlap, vB);

  size_type nsd = vB.size();

  cout << "Nomber of sub-domains = " << nsd + (USECOARSE != 0) << endl;
  
  if (USECOARSE) {
    vB.resize(nsd+1);
    cout << "interpolation coarse mesh\n";
    size_type nb_dof_coarse = mef_coarse.nb_dof();
    gmm::resize(vB[nsd], nb_dof, nb_dof_coarse);
    getfem::interpolation(mef_coarse, mef, vB[nsd], true);
    ++nsd;
  }
  
  gmm::iteration iter(residu, 1, 1000000);
  switch (version) {
  case 1 : gmm::sequential_additive_schwarz(RM, U, F,
	      gmm::ildlt_precond<general_sparse_matrix>(), vB, iter,
	      gmm::using_cg(), gmm::using_cg());
    break;
  case 2 : gmm::sequential_additive_schwarz(RM, U, F,
	      gmm::ilu_precond<general_sparse_matrix>(), vB, iter,
	      gmm::using_gmres(), gmm::using_gmres());
    break;
  }
  return 0;
}

  
struct exception_cb : public dal::exception_callback  {
   virtual void callback(const std::string& msg)
   { cerr << msg << endl; *(int *)(0) = 0; } 
};

int main(int argc, char *argv[]) {
   exception_cb cb;
   dal::exception_callback::set_exception_callback(&cb);

  try {
    ftool::md_param params;
    pb_data p;
    
    std::cout << "initialization ...\n";
    params.read_command_line(argc, argv);
    p.init(params);
    p.mesh.stat();
    
    p.assemble();

    double rutime = ftool::uclock_sec();
    int itebilan = p.solve();
    std::cout << "resolution time : " << ftool::uclock_sec() - rutime << endl;
    cout << "itebilan = " << itebilan << endl;

    gmm::mult(p.RM, gmm::scaled(p.U, -1.0), p.F, p.F);
    cout << "final residu : " << gmm::vect_norm2(p.F) << endl;
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0;
}
