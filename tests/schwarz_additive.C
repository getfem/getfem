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

typedef gmm::rsvector<scalar_type>     sparse_vector;
typedef gmm::row_matrix<sparse_vector> general_sparse_matrix;
typedef gmm::row_matrix<sparse_vector> symmetric_sparse_matrix;
typedef std::vector<scalar_type>       linalg_vector;



struct pb_data {
  getfem::getfem_mesh mesh;
  getfem::getfem_mesh mesh_coarse;
  getfem::mesh_fem mef_coarse;

  getfem::mesh_fem mef;
  getfem::mesh_fem mef_data;

  double PG, mu, lambda, rho;
  double LX, LY, LZ, residu, overlap;
  int NX, N, NXCOARSE, USECOARSE;
  int K;     /* finite element degree.                                    */
  base_vector D;

  std::vector<size_type> nsdm;

  symmetric_sparse_matrix RM; /* stifness matrix.                         */
  linalg_vector U, F;         /* Unknown and right hand side.             */
  int solver;

  std::string datafilename;

  void assemble(void);
  void init(void);

  int solve_cg(void);
  int solve_schwarz(int);

  int solve(void) {
    switch (solver) {
    case 0 : return solve_cg();  
    case 1 : return solve_schwarz(1); 
    case 2 : return solve_schwarz(2); 
    }
    return -1;
  }

  pb_data(void) : mef_coarse(mesh_coarse), mef(mesh), mef_data(mesh) {}
};

ftool::md_param PBSTFR_PARAM;
double gravity = 9.81, prho = 1.0;


base_vector second_membre(const base_node &x) {
  int N = x.size();
  base_vector res(N);
  double z = 0; for (int i = 0; i < N-1; i++) z += x[i];
  res.fill(0.0); res[N-1] = -prho*gravity;
  return res;
}

void pb_data::init(void) {
  dal::bit_vector nn;
  int i, j, k;

  /***********************************************************************/
  /*  READING PARAMETER FILE.                                            */
  /***********************************************************************/
  
  /* parametres physiques */
  N = PBSTFR_PARAM.int_value("N", "Dimension");
  mu = PBSTFR_PARAM.real_value("MU", "Stiffness parameter mu");
  PG = PBSTFR_PARAM.real_value("PG", "G");
  prho = rho = PBSTFR_PARAM.real_value("RHO", "RHO");
  gravity = PG * rho;
  lambda = PBSTFR_PARAM.real_value("LAMBDA", "lambda");
  D.resize(N); D.fill(0.0);
  D[N-1] = PBSTFR_PARAM.real_value("D", "Dirichlet condition");
  
  /* parametres numeriques */
  LX = PBSTFR_PARAM.real_value("LX", "Size in X");
  LY = PBSTFR_PARAM.real_value("LY", "Size in Y");
  LZ = PBSTFR_PARAM.real_value("LZ", "Size in Y");
  NX = PBSTFR_PARAM.int_value("NX", "Nomber of space step ");
  NXCOARSE = PBSTFR_PARAM.int_value("NXCOARSE", "Nombre of space step ");
  USECOARSE = PBSTFR_PARAM.int_value("USECOARSE", "Coarser mesh or not");
  residu = PBSTFR_PARAM.real_value("RESIDU", "residu");
  overlap = PBSTFR_PARAM.real_value("OVERLAP", "overlap");
  K = PBSTFR_PARAM.int_value("K", "Degree");
  solver = PBSTFR_PARAM.int_value("SOLVER", "solver");
  nsdm.resize(std::max(3, N));
  nsdm[0] = PBSTFR_PARAM.int_value("NSDMX", "Nomber of sub-domains");
  nsdm[1] = PBSTFR_PARAM.int_value("NSDMY", "Nombre of sub-domains");
  nsdm[2] = PBSTFR_PARAM.int_value("NSDMZ", "Nombre of sub-domains");
  for (int i = 3; i < N; ++i) nsdm[i] = nsdm[2];
  
  datafilename = std::string(PBSTFR_PARAM.string_value("ROOTFILENAME",
			     "file name"));
  std::string meshname(PBSTFR_PARAM.string_value("MESHNAME",
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
    for (i = 0; i < N; i++)
    { 
      vtab[i] = bgeot::base_small_vector(N); vtab[i].fill(0.0);
      (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
    }
    getfem::parallelepiped_regular_simplex_mesh(mesh, N, org,
						&(vtab[0]), &(ref[0]));
  }



  {
    base_node org(N); org.fill(0.0);
    std::vector<bgeot::base_small_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NXCOARSE);
    for (i = 0; i < N; i++)
    { 
      vtab[i] = bgeot::base_small_vector(N); vtab[i].fill(0.0);
      (vtab[i])[i] = 
	((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NXCOARSE);
    }
    getfem::parallelepiped_regular_simplex_mesh(mesh_coarse, N, org,
						&(vtab[0]), &(ref[0]));
  }
  
  mesh.trans_of_convex(0);
  mesh.optimize_structure();

  nn = mesh.convex_index(N);
  char method[500];
  sprintf(method, "IM_EXACT_SIMPLEX(%d)", N);
  getfem::pintegration_method ppi = getfem::int_method_descriptor(method);
  
  sprintf(method, "FEM_PK(%d, %d)", N, K);
  mef.set_finite_element(nn, getfem::fem_descriptor(method), ppi);
  mef_coarse.set_finite_element(mesh_coarse.convex_index(N),
				getfem::fem_descriptor(method), ppi);
  sprintf(method, "FEM_PK(%d, %d)", N, K);
  mef_data.set_finite_element(nn, getfem::fem_descriptor(method), ppi);
  mef.set_qdim(N);
  mef_coarse.set_qdim(N);

  nn = mesh.convex_index(N);
  base_vector un(N);
  for (j << nn; j >= 0; j << nn) {
    k = mesh.structure_of_convex(j)->nb_faces();
    for (i = 0; i < k; i++) {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty()) {
	gmm::copy(mesh.normal_of_face_of_convex(j, i, 0), un);
	un /= bgeot::vect_norm2(un);	
	if (dal::abs(un[N-1] - 1.0) < 1.0E-3)
	  mef.add_boundary_elt(0, j, i);
      }
    }
  }
}

void pb_data::assemble(void) {
  size_type nb_dof = mef.nb_dof();
  std::cout << "number of dof : "<< nb_dof << endl;
  size_type nb_dof_data = mef_data.nb_dof();
  dal::bit_vector ddlD = mef.dof_on_boundary(0);
 
  F = linalg_vector(nb_dof);
  gmm::clear(F);
  U = linalg_vector(nb_dof);
  gmm::clear(U);
  RM = symmetric_sparse_matrix(nb_dof, nb_dof);

  std::cout << "Assembly of stiffness matrix" << endl;

  linalg_vector STLA(nb_dof_data), STG(nb_dof_data);
  std::fill(STLA.begin(), STLA.end(), lambda);
  std::fill(STG.begin(), STG.end(), mu);
  getfem::asm_stiffness_matrix_for_linear_elasticity(RM,mef,mef_data,STLA,STG);

  std::cout << "Assembly of source term" << endl;
  linalg_vector STF(N * nb_dof_data);
  for (size_type j = 0; j < nb_dof_data; j++)
    for (int k = 0; k < N; k++)
      STF[j*N + k] = (second_membre(mef_data.point_of_dof(j)))[k];
  getfem::asm_source_term(F, mef, mef_data, STF);
  
  linalg_vector UD(nb_dof);
  for (size_type j = 0; j < nb_dof/N; j++)
    for (size_type k = 0; k < size_type(N); k++) UD[j*N + k] = D[k];

  getfem::assembling_Dirichlet_condition(RM, F, mef, 0, UD);
}

std::vector<size_type> extract_sub_domain(const getfem::getfem_mesh &mesh,
					  getfem::mesh_fem &mef,
					  const base_vector &min,
					  const base_vector &max) {
  size_type i, j;
  const base_node *pt;
  dal::bit_vector mm = mesh.convex_index(), nn;
  for (i << mm; i != size_type(-1); i << mm) {
    bool in = false;
    for (j = 0; j < mesh.nb_points_of_convex(i) && !in; ++j) {
      pt = &(mesh.points_of_convex(i)[j]);
      in = true;
      for (size_type k = 0; k < pt->size(); ++k)
	if ((*pt)[k] <  min[k] || (*pt)[k] > max[k])
	  { in = false; break; }
    }
    if (in) {
      for (j = 0; j < mef.nb_dof_of_element(i); ++j)
	nn.add(mef.ind_dof_of_element(i)[j]);
    }
  }
  std::vector<size_type> res(nn.card());
  for (i << nn, j = 0; i != size_type(-1); i << nn, ++j) res[j] = i;
  return res;
}

int pb_data::solve_cg(void) {
  gmm::iteration iter(residu, 1, 1000000);
  gmm::ildltt_precond<general_sparse_matrix> P(RM, 10, 1E-7);
  gmm::cg(RM, U, F, gmm::identity_matrix(), P, iter);
  return iter.get_iteration();
}

int pb_data::solve_schwarz(int version) {
  size_type nsd = 1, nb_dof = gmm::mat_nrows(RM);
  for (int i = 0; i < N; ++i) nsd *= nsdm[i];
  std::vector<scalar_type> L(N);
  L[0] = LX; if (N >= 1) L[1] = LY;
  for (int i = 2; i < N; ++i) L[i] = LZ;

  cout << "Nomber of sub-domains = " << nsd + (USECOARSE != 0) << endl;

  std::vector< gmm::sub_index > index_tab(nsd);
  std::vector<general_sparse_matrix> vB(nsd+1);

  std::vector<size_type> ind(N);
  std::fill(ind.begin(), ind.end(), 0);
  for ( size_type n = 0; n < nsd; ++n) {

    base_vector min(N), max(N);
    for (int i = 0; i < N; ++i) {
      min[i] = (L[i] / nsdm[i]) * (ind[i] - overlap); 
      max[i] = (L[i] / nsdm[i]) * (ind[i] + 1.0 + overlap);
    }
    std::vector<size_type> sd = extract_sub_domain(mesh, mef, min, max);
    gmm::resize(vB[n], sd.size(), nb_dof);

    for (size_type i = 0; i < sd.size(); ++i)
      vB[n](i, sd[i]) = 1.0;
  
    for (int i = 0; i < N; ++i) {
      (ind[i])++;
      if (ind[i] == nsdm[i]) { ind[i] = 0; } else break;
    }
  }
  
  if (USECOARSE) {
    cout << "interpolation coarse mesh\n";
    size_type nb_dof_coarse = mef_coarse.nb_dof();
    gmm::resize(vB[nsd], nb_dof_coarse, nb_dof);
    general_sparse_matrix aux(nb_dof, nb_dof_coarse);
    getfem::interpolation_solution(mef_coarse, mef, aux);
    gmm::copy(gmm::transposed(aux), vB[nsd]);
  }
  else resize(vB, nsd);
  
  gmm::iteration iter(residu, 1, 1000000);
  switch (version) {
  case 1 :
    return gmm::sequential_additive_schwarz(RM, U, F,
	      gmm::ildltt_precond<general_sparse_matrix>(10, 1E-7), vB, iter,
					gmm::using_cg(), gmm::using_cg());
  case 2 :
    return gmm::sequential_additive_schwarz(RM, U, F,
	      gmm::ilut_precond<general_sparse_matrix>(10, 1E-7), vB, iter,
				     gmm::using_gmres(), gmm::using_gmres());
  }

  return 0;
}

  
class exception_cb : public dal::exception_callback  {
   public:
   virtual void callback(const std::string& msg)
   { cerr << msg << endl; *(int *)(0) = 0; } 
};

int main(int argc, char *argv[])
{
   exception_cb cb;
   dal::exception_callback::set_exception_callback(&cb);

  try {
    pb_data p;
    
    std::cout << "initialization ...\n";
    PBSTFR_PARAM.read_command_line(argc, argv);
    p.init();
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
