/*===========================================================================
 
 Copyright (C) 2007-2012 Yves Renard, Julien Pommier.
 
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
#define GMM_USES_SUPERLU
/**************************************************************************/
/*                                                                        */
/*  Schwarz additive test program on an elastostatic problem with         */
/*  an optional coarse mesh.                                              */
/*                                                                        */
/**************************************************************************/

#define GMM_USES_SUPERLU

#include <getfem/getfem_assembling.h>
#include <getfem/getfem_norm.h>
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_export.h>
#include <gmm/gmm.h>
#ifdef GMM_USES_MPI
#include <mpi++.h>
#endif

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;
using bgeot::base_vector;
using bgeot::base_node;
using bgeot::size_type;
using bgeot::short_type;
using bgeot::dim_type;
using bgeot::scalar_type;

typedef gmm::row_matrix<gmm::rsvector<scalar_type> > general_sparse_matrix;
typedef std::vector<scalar_type>                     linalg_vector;


struct pb_data {
  getfem::mesh mesh;
  getfem::mesh mesh_coarse;

  getfem::mesh_im  mim;
  getfem::mesh_fem mef;
  getfem::mesh_fem mef_data;
  getfem::mesh_fem mef_coarse;

  double mu, lambda, rho, gravity;
  double LX, LY, LZ, residual, overlap, subdomsize;
  int NX, N, NXCOARSE, USECOARSE, K;
  base_vector D;

  general_sparse_matrix RM;   /* stifness matrix.                         */
  linalg_vector U, F;         /* Unknown and right hand side.             */
  int solver;

  void assemble(void);
  void init(ftool::md_param &params);

  int solve_cg(void);
  int solve_cg2(void);
  int solve_superlu(void);
  int solve_schwarz(int);

  int solve(void) {
    cout << "solving" << endl;
    switch (solver) {
    case 0 : return solve_cg();
    case 1 : return solve_cg2();
    case 2 : return solve_superlu();
    default : return solve_schwarz(solver);
    }
    return 0;
  }

  base_vector vol_force(const base_node &x)
  { base_vector res(x.size()); res[N-1] = -rho*gravity; return res; }

  pb_data(void) : mim(mesh), mef(mesh), mef_data(mesh), mef_coarse(mesh_coarse)  {}
};

ftool::md_param PBSTFR_PARAM;

void pb_data::init(ftool::md_param &params) {

  /***********************************************************************/
  /*  READING PARAMETER FILE.                                            */
  /***********************************************************************/
  
  /* parametres physiques */
  N = int(params.int_value("N", "Dimension"));
  mu = params.real_value("MU", "Stiffness parameter mu");
  gravity = params.real_value("PG", "G");
  rho = params.real_value("RHO", "RHO");
  lambda = params.real_value("LAMBDA", "lambda");
  D.resize(N); gmm::clear(D);
  D[N-1] = params.real_value("D", "Dirichlet condition");
  
  /* parametres numeriques */
  LX = params.real_value("LX", "Size in X");
  LY = params.real_value("LY", "Size in Y");
  LZ = params.real_value("LZ", "Size in Y");
  NX = int(params.int_value("NX", "Nomber of space step "));
  NXCOARSE = int(params.int_value("NXCOARSE", "Nombre of space step "));
  USECOARSE = int(params.int_value("USECOARSE", "Coarser mesh or not"));
  residual = params.real_value("RESIDUAL", "residual");
  overlap = params.real_value("OVERLAP", "overlap");
  K = int(params.int_value("K", "Degree"));
  solver = int(params.int_value("SOLVER", "solver"));
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
    base_node org(N); gmm::clear(org);
    std::vector<bgeot::base_small_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
    for (int i = 0; i < N; i++) { 
      vtab[i] = bgeot::base_small_vector(N); gmm::clear(vtab[i]);
      (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
    }
    getfem::parallelepiped_regular_simplex_mesh(mesh, dim_type(N), org,
						&(vtab[0]), &(ref[0]));
  }

  if (USECOARSE) { // coarse mesh
    base_node org(N); gmm::clear(org);
    std::vector<bgeot::base_small_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NXCOARSE);
    for (int i = 0; i < N; i++) { 
      vtab[i] = bgeot::base_small_vector(N); gmm::clear(vtab[i]);
      (vtab[i])[i] = 
	((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NXCOARSE);
    }
    getfem::parallelepiped_regular_simplex_mesh(mesh_coarse, dim_type(N), org,
						&(vtab[0]), &(ref[0]));
  }
  
  mesh.trans_of_convex(0);
  mesh.optimize_structure();

  dal::bit_vector nn = mesh.convex_index(dim_type(N));
  char method[500];
  sprintf(method, "IM_EXACT_SIMPLEX(%d)", N);
  getfem::pintegration_method ppi = getfem::int_method_descriptor(method);
  
  sprintf(method, "FEM_PK(%d, %d)", N, K);
  mim.set_integration_method(nn, ppi);
  mef.set_finite_element(nn, getfem::fem_descriptor(method));
  mef_coarse.set_finite_element(mesh_coarse.convex_index(dim_type(N)),
				getfem::fem_descriptor(method));
  mef_data.set_finite_element(nn, getfem::fem_descriptor(method));
  mef.set_qdim(dim_type(N));
  mef_coarse.set_qdim(dim_type(N));

  nn = mesh.convex_index(dim_type(N));
  base_vector un(N);
  for (int j = nn.take_first(); j >= 0; j << nn) {
    int k = mesh.structure_of_convex(j)->nb_faces();
    for (short_type i = 0; i < k; i++) {
      if (mesh.is_convex_having_neighbour(j, i)) {
	gmm::copy(mesh.normal_of_face_of_convex(j, i, 0), un);
	gmm::scale(un, 1/gmm::vect_norm2(un));
	if (gmm::abs(un[N-1] - 1.0) < 1.0E-3) mesh.region(0).add(j, i);
      }
    }
  }
}

void pb_data::assemble(void) {
  size_type nb_dof = mef.nb_dof();
  std::cout << "number of dof : "<< nb_dof << endl;
  size_type nb_dof_data = mef_data.nb_dof();
  dal::bit_vector ddlD = mef.basic_dof_on_region(0);
 
  F.resize(nb_dof); gmm::clear(F);
  U.resize(nb_dof); gmm::clear(U);
  gmm::resize(RM, nb_dof, nb_dof);

  std::cout << "Assembly of stiffness matrix" << endl;

  linalg_vector STLA(nb_dof_data), STG(nb_dof_data);
  std::fill(STLA.begin(), STLA.end(), lambda);
  std::fill(STG.begin(), STG.end(), mu);
  getfem::asm_stiffness_matrix_for_linear_elasticity(RM,mim,mef,mef_data,STLA,STG);

  std::cout << "Assembly of source term" << endl;
  linalg_vector STF(N * nb_dof_data);
  for (size_type j = 0; j < nb_dof_data; j++)
    for (int k = 0; k < N; k++)
      STF[j*N + k] = (vol_force(mef_data.point_of_basic_dof(j)))[k];
  getfem::asm_source_term(F, mim, mef, mef_data, STF);
  
  linalg_vector UD(nb_dof);
  for (size_type j = 0; j < nb_dof/N; j++)
    for (size_type k = 0; k < size_type(N); k++) UD[j*N + k] = D[k];
  getfem::assembling_Dirichlet_condition(RM, F, mef, 0, UD);
}

int pb_data::solve_cg(void) {
  gmm::iteration iter(residual, 1, 1000000);
  gmm::ildlt_precond<general_sparse_matrix> P(RM);
  gmm::cg(RM, U, F, gmm::identity_matrix(), P, iter);
  return int(iter.get_iteration());
}

int pb_data::solve_superlu(void) {
  double rcond;
  SuperLU_solve(RM, U, F, rcond);
  return 1;
}

int pb_data::solve_cg2(void) {
  gmm::iteration iter(residual, 1, 1000000);
  gmm::cg(RM, U, F, gmm::identity_matrix(), gmm::identity_matrix(), iter);
  return int(iter.get_iteration());
}

int pb_data::solve_schwarz(int version) {

  size_type nb_dof = mef.nb_dof();
  std::vector<base_node> pts(nb_dof);
  for (size_type i = 0; i < nb_dof; ++i) pts[i] = mef.point_of_basic_dof(i);

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
  
  gmm::iteration iter(residual, 1, 1000000);
  switch (version) {
  case 3 : gmm::additive_schwarz(RM, U, F,
	      gmm::ildlt_precond<general_sparse_matrix>(), vB, iter,
	      gmm::using_cg(), gmm::using_cg());
    break;
  case 4 : gmm::additive_schwarz(RM, U, F,
	      gmm::ilu_precond<general_sparse_matrix>(), vB, iter,
	      gmm::using_gmres(), gmm::using_gmres());
    break;
  case 5 : gmm::additive_schwarz(RM, U, F,
	      gmm::ilu_precond<general_sparse_matrix>(), vB, iter,
	      gmm::using_superlu(), gmm::using_cg());
    break;
  }
  return 0;
}

  
int main(int argc, char *argv[]) {
#ifdef GMM_USES_MPI
    MPI_Init(&argc,&argv);
#endif
 
  try {
    ftool::md_param params;
    pb_data p;
    
    std::cout << "initialization ...\n";
    params.read_command_line(argc, argv);
    p.init(params);
    p.mesh.stat();
    
    p.assemble();

    double rutime = dal::uclock_sec();
    int itebilan = p.solve();
    std::cout << "resolution time : " << dal::uclock_sec() - rutime << endl;
    cout << "itebilan = " << itebilan << endl;

    gmm::mult(p.RM, gmm::scaled(p.U, -1.0), p.F, p.F);
    cout << "final residual : " << gmm::vect_norm2(p.F) << endl;
  }
  GMM_STANDARD_CATCH_ERROR;
#ifdef GMM_USES_MPI
   MPI_Finalize();
#endif
  return 0;
}
