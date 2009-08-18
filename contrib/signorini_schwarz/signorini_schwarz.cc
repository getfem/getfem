// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2009 Yves Renard.
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

/**************************************************************************/
/*                                                                        */
/*  Calcul d'une structure lineairement elastique en contact unilateral   */
/*  sans frottement avec une fondation rigide en mouvement uniforme.      */
/*                                                                        */
/**************************************************************************/

#include <getfem/getfem_assembling.h>
#include <getfem/getfem_norm.h>
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_export.h>
#include <gmm/gmm.h>

using std::cout;

using bgeot::base_vector;
using bgeot::base_node;
using bgeot::size_type;
using bgeot::scalar_type;

typedef gmm::rsvector<scalar_type> sparse_vector;
typedef gmm::row_matrix<sparse_vector> general_sparse_matrix;
typedef gmm::row_matrix<sparse_vector> symmetric_sparse_matrix;
typedef std::vector<scalar_type>       linalg_vector;



struct pb_data
{
  getfem::mesh mesh;
  getfem::mesh mesh_coarse;
  getfem::mesh_fem mef_coarse;

  getfem::mesh_fem mef;
  getfem::mesh_fem mef_P1;     /* for Taylor-Hood version.                */
  getfem::mesh_fem mef_data;   /* for Taylor-Hood version.                */
  getfem::mesh_im mim;

  double PG, mu, lambda, rho;
  double LX, LY, LZ, residu, cUsawa, cgradient, overlap, subdomsize;
  bgeot::dim_type N;
  int NX, NXCOARSE, USECOARSE;
  int K;     /* finite element degree.                                    */
  base_vector D;

  double critip;

  general_sparse_matrix B;    /* for Taylor-Hood version.                 */
  symmetric_sparse_matrix M;  /* for Taylor-Hood version.                 */
  symmetric_sparse_matrix RM; /* matrice de rigidite.                     */
  linalg_vector U, F, UD;     /* inconnue et second membre.               */
  int option;
  int solver;

  int NBLEVEL;
  bool draw;
  std::string datafilename;

  double muldep;

  void assemble(void);
  void init(void);
  void visu(void);

  general_sparse_matrix CO;
  linalg_vector S;

  int solve_cg(void);  
  int solve_cgc(void);  
  int solve_schwarz(int);
  int solve_superlu(void);

  int solve(void) {
    switch (solver) {
    case 0 : return solve_cg();
    case 1 : return solve_cgc();
    case 2 : return solve_superlu();  
    case 3 : return solve_schwarz(0); 
    case 4 : return solve_schwarz(1); 
    case 5 : return solve_schwarz(2); 
    }
    return -1;
  }

  pb_data(void) : mef_coarse(mesh_coarse), mef(mesh), mef_P1(mesh),
		  mef_data(mesh), mim(mesh) {}
};

ftool::md_param PBSTFR_PARAM;
double gravity = 9.81, prho = 1.0;
double perts = 10.0;


base_vector second_membre(const base_node &x)
{
  int N = x.size();
  base_vector res(N);
  double z = 0; for (int i = 0; i < N-1; i++) z += x[i];
  res[N-1] = prho*gravity*(perts*cos(2.0*M_PI*z) - 1.0);
  return res;
}

void pb_data::init(void)
{
  dal::bit_vector nn;
  size_type j, k;

  /***********************************************************************/
  /*  LECTURE DES PARAMETRES SUR FICHIER.                                */
  /***********************************************************************/
  
  /* parametres physiques */
  N = bgeot::dim_type(PBSTFR_PARAM.int_value("N", "Dimension de l'espace"));
  mu = PBSTFR_PARAM.real_value("MU", "Raideur de cisaillement");
  PG = PBSTFR_PARAM.real_value("PG", "Acceleration de gravitation");
  prho = rho = PBSTFR_PARAM.real_value("RHO", "Densite massique");
  gravity = PG * rho;
  lambda = PBSTFR_PARAM.real_value("LAMBDA", "Coeff elastique lambda");
  D.resize(N); gmm::fill(D, 0.0);
  D[N-1] = PBSTFR_PARAM.real_value("D", "parametre condition de Dirichlet");
  perts = PBSTFR_PARAM.real_value("PERTS", "perturbation terme source");
  critip = PBSTFR_PARAM.real_value("CRITIP",
				   "Critere arret pour points interieurs");

  /* parametres numeriques */
  LX = PBSTFR_PARAM.real_value("LX", "Taille en X");
  LY = PBSTFR_PARAM.real_value("LY", "Taille en Y");
  LZ = PBSTFR_PARAM.real_value("LZ", "Taille en Y");
  NX = PBSTFR_PARAM.int_value("NX", "Nombre de pas d'espace ");
  NXCOARSE = PBSTFR_PARAM.int_value("NXCOARSE", "Nombre de pas d'espace ");
  USECOARSE = PBSTFR_PARAM.int_value("USECOARSE", "Coarser mesh or not");
  residu = PBSTFR_PARAM.real_value("RESIDU", "Valeur pour test d'arret");
  muldep = PBSTFR_PARAM.real_value("MULDEP", "Coeff de grossissement");
  overlap = PBSTFR_PARAM.real_value("OVERLAP",
				    "taille des recouvrement");
  subdomsize = PBSTFR_PARAM.real_value("SUBDOMSIZE",
				       "taille des sous domaines");
  cUsawa = PBSTFR_PARAM.real_value("CUSAWA", "Coeff pour algo USAWA");
  cgradient = PBSTFR_PARAM.real_value("CGRADIENT", "Coeff pour algo gradient");
  NBLEVEL = PBSTFR_PARAM.int_value("NBLEVEL", "Nombre de courbes de niveau");
  K = PBSTFR_PARAM.int_value("K", "Degre de l'element fini de Lagrange");
  option = PBSTFR_PARAM.int_value("OPTION", "Option de resolution");
  solver = PBSTFR_PARAM.int_value("SOLVER", "Type de solver");
  datafilename = std::string(PBSTFR_PARAM.string_value("ROOTFILENAME",
			     "Nom du fichier de sauvegarde sans extension"));
  std::string meshname(PBSTFR_PARAM.string_value("MESHNAME",
			     "Nom du fichier de maillage"));
  std::string dessin = PBSTFR_PARAM.string_value("GRAPHICS", "Dessin ? ");
  const char *dds = dessin.c_str();
  if (!strcmp("N", dds) || !strcmp("n", dds))
   { draw = false; }
  else
   { draw = true; }

  std::cout << "\n\n";

  /***********************************************************************/
  /*  CONSTRUCTION DU MAILLAGE.                                          */
  /***********************************************************************/

  std::cout << "construction du maillage\n";

  if (meshname.size() > 0)
  {
    mesh.read_from_file(meshname);
  }
  else
  {
    base_node org(N); org.fill(0.0);
    std::vector<bgeot::base_small_vector> vtab(N);
    std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
    for (bgeot::dim_type i = 0; i < N; i++)
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
    for (bgeot::dim_type i = 0; i < N; i++)
    { 
      vtab[i] = bgeot::base_small_vector(N); vtab[i].fill(0.0);
      (vtab[i])[i] = 
	((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NXCOARSE);
    }
    getfem::parallelepiped_regular_simplex_mesh
      (mesh_coarse, N, org, &(vtab[0]), &(ref[0]));
  }
  

  mesh.trans_of_convex(0);

  std::cout << "Optimisation de la structure\n";
  mesh.optimize_structure();

  std::cout << "Choix de l'element fini.\n";
  nn = mesh.convex_index(N);
  std::cout << "Nombre d'elements : " << nn.card() << endl;
  char method[500];
  sprintf(method, "IM_EXACT_SIMPLEX(%d)", N);
  getfem::pintegration_method ppi = getfem::int_method_descriptor(method);
  mim.set_integration_method(nn, ppi);
  if (option & 1)
  {
    K = 2;
    sprintf(method, "FEM_PK(%d, %d)", N, 1);
    mef_P1.set_finite_element(nn, getfem::fem_descriptor(method));
    sprintf(method, "FEM_PK(%d, %d)", N, 2);
    mef.set_finite_element(nn, getfem::fem_descriptor(method));
    mef_coarse.set_finite_element(mesh_coarse.convex_index(N),
				  getfem::fem_descriptor(method));
  }
  else if (option & 16)
  {
    K = 1; assert(N == 2);
    // mef.set_finite_element(nn, getfem::PK_fem(N, 1));
    sprintf(method, "FEM_P1_NONCONFORMING");
    mef.set_finite_element(nn, getfem::fem_descriptor(method));
    mef_coarse.set_finite_element(mesh_coarse.convex_index(N),
				  getfem::fem_descriptor(method));
  }
  else
  {
    sprintf(method, "FEM_PK(%d, %d)", N, K);
    mef.set_finite_element(nn, getfem::fem_descriptor(method));
    mef_coarse.set_finite_element(mesh_coarse.convex_index(N),
				  getfem::fem_descriptor(method));
  }
  sprintf(method, "FEM_PK(%d, %d)", N, K);
  mef_data.set_finite_element(nn, getfem::fem_descriptor(method));
  mef.set_qdim(N);
  mef_coarse.set_qdim(N);
 
  std::cout << "Reperage des bords de contact et Dirichlet\n";
  nn = mesh.convex_index(N);
  bgeot::base_small_vector un(N);
  for (j << nn; j != size_type(-1); j << nn)
  {
    k = mesh.structure_of_convex(j)->nb_faces();
    for (bgeot::short_type i = 0; i < k; i++) {
      if (!mesh.is_convex_having_neighbour(j, i)) {
	gmm::copy(mesh.normal_of_face_of_convex(j, i, 0), un);
	un /= gmm::vect_norm2(un);	
	
        //if (true)
	if (gmm::abs(un[N-1] + 1.0) < 1.0E-3)
	  mesh.region(1).add(j, i);
	else  if (gmm::abs(un[N-1] - 1.0) < 1.0E-3)
	  mesh.region(0).add(j, i);
      }
    }
  }
}

void pb_data::assemble(void) {
  GMM_ASSERT1(!mef.is_reduced(), "To be adapted");
  size_type nb_dof = mef.nb_dof();
  std::cout << "nombre de ddl pour l'elasticite lineaire : "<< nb_dof << endl;
  size_type nb_dof_data = mef_data.nb_dof();
  dal::bit_vector ddlD = mef.basic_dof_on_region(0);
  dal::bit_vector ddlC = mef.basic_dof_on_region(1);
 
  F = linalg_vector(nb_dof);
  gmm::clear(F);
  U = linalg_vector(nb_dof);
  gmm::clear(U);
  RM = symmetric_sparse_matrix(nb_dof, nb_dof);

  // std::cout << "Dirichlet nodes : " << ddlD << endl;
  // std::cout << "Contact nodes : " << ddlC << endl;

  std::cout << "Assemblage de la matrice de rigidite" << endl;

  linalg_vector STLA(nb_dof_data), STG(nb_dof_data);
  std::fill(STLA.begin(), STLA.end(), lambda);
  std::fill(STG.begin(), STG.end(), mu);
  getfem::asm_stiffness_matrix_for_linear_elasticity(RM,mim,mef,mef_data,STLA,STG);

  std::cout << "Assemblage du terme source" << endl;
  linalg_vector STF(N * nb_dof_data);
  getfem::interpolation_function(mef_data, STF, second_membre);
  getfem::asm_source_term(F, mim, mef, mef_data, STF);
  
  UD = linalg_vector(nb_dof);
  for (size_type j = 0; j < nb_dof/N; j++)
    for (size_type k = 0; k < size_type(N); k++) UD[j*N + k] = D[k];
  std::cout << "Nettoyage de la matrice RM" << endl;
  getfem::assembling_Dirichlet_condition(RM, F, mef, 0, UD);
  
  if (option & 2) ddlC.clear();  // probleme lineaire.
  // std::cout << "Contact nodes : " << ddlC << endl;
  CO = general_sparse_matrix(ddlC.card() / N, U.size());
  S = linalg_vector(ddlC.card() / N);
  gmm::clear(S);
  size_type j, i = 0;
  bgeot::dim_type k = 0;
  for (j << ddlC; j != size_type(-1); j << ddlC,
	 k = bgeot::dim_type((k + 1) % N))
    if (k == N-1) { CO(i, j) = -1.0; ++i; }

  // cout << "Matrice des contraintes : " << CO << endl;

}

void pb_data::visu(void) {
  dal::dynamic_array<base_node> ptab;
  dal::bit_vector nn = mef.convex_index();
  size_type c, i;
  for (i = 0; i <= size_type(mesh.points().index().last()); ++i)
    { ptab[i] = base_node(N); ptab[i].fill(0.0); }
  
  for (c << nn; c != size_type(-1); c << nn) {
    for (i = 0; i < mef.nb_basic_dof_of_element(c); ++i) {
      size_type l = mef.ind_basic_dof_of_element(c)[i];
      size_type k = mesh.search_point(mef.point_of_basic_dof(l));
      if (k != size_type(-1)) ptab[k][i % N] = U[l] * muldep;
    }
  }
  // simple_mesh_visu(mesh, "Configuration deformee", NX <= 5, ptab);
}

int pb_data::solve_cg(void) {
  gmm::iteration iter(residu, 1, 1000000);
  gmm::ildltt_precond<general_sparse_matrix> P(RM, 10, 1E-7);
  gmm::cg(RM, U, F, gmm::identity_matrix(), P, iter);
  return iter.get_iteration();
} 

int pb_data::solve_superlu(void) {
  scalar_type rcond;
  gmm::SuperLU_solve(RM, U, F, rcond);
  cout << "SuperLU condest = " << 1.0/rcond << "\n";
  return 1;
} 

int pb_data::solve_cgc(void) {
  gmm::iteration iter(residu, 1, 1000000);
  std::vector<double> COF(gmm::mat_nrows(CO));
  gmm::clear(COF);
  gmm::constrained_cg(RM, CO, U, F, COF, gmm::identity_matrix(),
			     gmm::identity_matrix(), iter);
  return iter.get_iteration();
}  

int pb_data::solve_schwarz(int version) {
  size_type nb_dof = gmm::mat_nrows(RM);

  std::vector<base_node> pts(nb_dof);
  for (size_type i = 0; i < nb_dof; ++i) pts[i] = mef.point_of_basic_dof(i);

  std::vector<general_sparse_matrix> vB;
  gmm::rudimentary_regular_decomposition(pts, subdomsize, overlap, vB);

  size_type nsd = vB.size();

  cout << "Nombre de sous domaines = " << nsd + (USECOARSE != 0) << endl;
  
  if (USECOARSE) {
    cout << "interpolation coarse mesh\n";
    vB.resize(nsd+1);
    size_type nb_dof_coarse = mef_coarse.nb_dof();
    gmm::resize(vB[nsd], nb_dof, nb_dof_coarse);
    // general_sparse_matrix aux(nb_dof, nb_dof_coarse);
    getfem::interpolation(mef_coarse, mef, vB[nsd]);
    // gmm::copy(aux, vB[nsd]);
    ++nsd;
  }
  
  gmm::iteration iter(residu, 1, 1000000);
  switch (version) {
  case 0 :
    gmm::additive_schwarz(RM, U, F, 
                                   gmm::identity_matrix(), vB, iter,
				   gmm::using_superlu(), gmm::using_cg());
    break;
  case 1 :
    gmm::additive_schwarz(RM, U, F,
	      gmm::ildltt_precond<general_sparse_matrix>(10, 1E-7), vB, iter,
					gmm::using_cg(), gmm::using_cg());
    break;
  case 2 :
    gmm::additive_schwarz(RM, U, F,
	      gmm::ilut_precond<general_sparse_matrix>(10, 1E-7), vB, iter,
				     gmm::using_gmres(), gmm::using_gmres());
    break;
  }

  return 0;
}

  
int main(int argc, char *argv[])
{

  try {
    pb_data p;
    dal::bit_vector nn;
    double utime;
    
    utime = dal::uclock_sec();
    
    std::cout << "initialisation ...\n";
    PBSTFR_PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.stat();
    // p.mesh.write_to_file(p.datafilename + ".mesh" + char(0));
    
    std::cout << "Debut de l'assemblage\n";
    p.assemble();

    // std::cout << "Stiffness matrix ="; gmm::write(p.RM, cout);
    
    std::cout << "Resolution, il y a " << p.U.size() << " ddl\n";
    double rutime = dal::uclock_sec();
    int itebilan = p.solve();
    std::cout << "temps de resolution = " << dal::uclock_sec() - rutime 
	      << endl;

    cout << "itebilan = " << itebilan << endl;

    gmm::mult(p.RM, gmm::scaled(p.U, -1.0), p.F, p.F);
    cout << "residu final : " << gmm::vect_norm2(p.F) << endl;

    // cout << "U = " << p.U << endl;
    
    if (p.draw) p.visu();
    
    // getfem::save_solution(p.datafilename + ".dataelt", p.mef, p.U,p.K);
    
    std::cout << "temps d'execution : " << dal::uclock_sec() - utime << endl;
    

  }
  GMM_STANDARD_CATCH_ERROR;
  return 0;
}
