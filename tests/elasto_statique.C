/**************************************************************************/
/*                                                                        */
/*  Calcul sur une structure lineairement elastique                       */
/*                                                                        */
/**************************************************************************/

#include <getfem_assembling.h>
#include <getfem_regular_meshes.h>

#include <linear_systems.h>
#include <model_param.h>
#include <graph_visu_3D.h>
#include <geo_mesh_visu_3D.h>
#include <clock.h>


/**************************************************************************/
/*  Definition de la solution test.                                       */
/**************************************************************************/

dal::dynamic_array<base_vector> sol_K;
scalar_type sol_lambda, sol_G;

base_vector sol_u(const base_node &x)
{
  int N = x.size(); base_vector res(N);
  for (int i = 0; i < N; ++i) res[i] = sin(bgeot::vect_sp(sol_K[i], x));
  return res;
}

base_vector sol_f(const base_node &x)
{
  int N = x.size();
  base_vector res(N);
  for (int i = 0; i < N; i++)
  {
    res[i] = ( sol_G * bgeot::vect_sp(sol_K[i], sol_K[i]) )
                  * sin(bgeot::vect_sp(sol_K[i], x));
    for (int j = 0; j < N; j++)
      res[i] += ( (sol_lambda + sol_G) * sol_K[j][j] * sol_K[j][i])
	          * sin(bgeot::vect_sp(sol_K[j], x));
  }
  return res;
}

base_matrix sol_sigma(const base_node &x)
{
  int N = x.size();
  base_matrix res(N,N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j <= i; j++)
    {
      res(j,i) = res(i,j) = sol_G *
	( sol_K[i][j] * cos(bgeot::vect_sp(sol_K[i], x))
       +  sol_K[j][i] * cos(bgeot::vect_sp(sol_K[j], x))
	);
      if (i == j)
	for (int k = 0; k < N; k++)
	  res(i,j) += sol_lambda * sol_K[k][k]
	                         * cos(bgeot::vect_sp(sol_K[k], x));
    }
  return res;
}

/**************************************************************************/
/*  Structure definissant le probleme.                                    */
/**************************************************************************/

struct pb_data
{
  getfem::getfem_mesh mesh;
  getfem::mesh_fem mef;
  getfem::mesh_fem mef_data;
  getfem::mesh_fem mef_data2;

  scalar_type G, lambda;
  scalar_type LX, LY, LZ, residu;
  int NX, N, K;

  symmetric_sparse_matrix RM; /* matrice de rigidite.                     */
  linalg_vector U, B; /* inconnue et second membre.                       */

  int iteimpl, integration;

  md_param PBSTFR_PARAM;

  int NBLEVEL;
  bool draw;
  string datafilename;
  
  scalar_type muldep;

  void init();
  void assemble(void);
  void solve(void);
  void visu(void);

  pb_data(void) : mef(mesh), mef_data(mesh), mef_data2(mesh) {}
};

void pb_data::init(void)
{
  dal::bit_vector nn;
  int i, j, k, h;

  /***********************************************************************/
  /*  LECTURE DES PARAMETRES SUR FICHIER.                                */
  /***********************************************************************/
  
  /* parametres physiques */
  N = PBSTFR_PARAM.int_value("N", "Dimension de l'espace");
  G = PBSTFR_PARAM.real_value("G", "Raideur de cisaillement");
  lambda = PBSTFR_PARAM.real_value("LAMBDA", "Coeff elastique lambda");
  
  /* parametres numeriques */
  LX = PBSTFR_PARAM.real_value("LX", "Taille en X");
  LY = PBSTFR_PARAM.real_value("LY", "Taille en Y");
  LZ = PBSTFR_PARAM.real_value("LZ", "Taille en Y");
  NX = PBSTFR_PARAM.int_value("NX", "Nombre de pas d'espace ");
  iteimpl = PBSTFR_PARAM.int_value("ITEIMPL", "Nombre d'iteration pt fixe");
  integration = PBSTFR_PARAM.int_value("INTEGRATION", "Type d'integration");
  residu = PBSTFR_PARAM.real_value("RESIDU", "Valeur pour test d'arret");
  muldep = PBSTFR_PARAM.real_value("MULDEP", "Coeff de grossissement");
  NBLEVEL = PBSTFR_PARAM.int_value("NBLEVEL", "Nombre de courbes de niveau");
  K = PBSTFR_PARAM.int_value("K", "Degre de l'element fini de Lagrange");

  datafilename = string(PBSTFR_PARAM.string_value("ROOTFILENAME",
			     "Nom du fichier de sauvegarde sans extension"));
  const char *dds = PBSTFR_PARAM.string_value("GRAPHICS", "Dessin ? ");
  draw = (strcmp("N", dds) && strcmp("n", dds));
 
  scalar_type FT = PBSTFR_PARAM.real_value("FT", "parametre pour solution exacte");
  for (i = 0; i < N; i++)
  {
    sol_K[i] = base_vector(N);
    for (j = 0; j < N; j++)  sol_K[i][j] = (i == j) ? FT : -FT;
  }

  sol_lambda = lambda; sol_G = G;

  /***********************************************************************/
  /*  CONSTRUCTION DU MAILLAGE.                                          */
  /***********************************************************************/

  cout << "construction du maillage\n";

  base_node org(N); org.fill(0.0);
  std::vector<base_vector> vtab(N);
  std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
  for (i = 0; i < N; i++)
  { 
    vtab[i] = base_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
  }
  getfem::parallelepiped_regular_simplex_mesh(mesh, N, org,
					     vtab.begin(), ref.begin());

  cout << "Optimisation de la structure\n";
  mesh.optimize_structure();

  cout << "Choix de l'element fini.\n";
  getfem::pintegration_method ppi;
  nn = mesh.convex_index(N);
  if (integration == 0) ppi = bgeot::simplex_poly_integration(N);
  else ppi = bgeot::Newton_Cotes_approx_integration(N,2*K);
  mef.set_finite_element(nn, getfem::PK_fem(N, K), ppi);

  mef_data.set_finite_element(nn, getfem::PK_fem(N, K), ppi);
  mef_data2.set_finite_element(nn, getfem::PK_fem(N, 0), ppi);

  if (draw) simple_mesh_visu(mesh, "Maillage", NX <= 5);

  cout << "Reperage des bord Neumann et Dirichlet\n";
  nn = mesh.convex_index(N);
  base_vector un;
  for (j << nn; j >= 0; j << nn)
  {
    k = mesh.structure_of_convex(j)->nb_faces();
    for (i = 0; i < k; i++)
    {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty())
      {
	un = mesh.convex(j).unit_norm_of_face(i);
	
	// if (true)
	if (abs(un[N-1] - 1.0) < 1.0E-3)
	  mef.add_boundary_elt(0, j, i);
	else
	  mef.add_boundary_elt(1, j, i);
      }
    }
  }
}

void pb_data::assemble(void)
{
  int nb_dof = mef.nb_dof(), nb_dof_data = mef_data.nb_dof();
  int nb_dof_data2 = mef_data2.nb_dof();
  B = linalg_vector(N*nb_dof, 0.0);
  U = linalg_vector(N*nb_dof, 0.0);
  RM = symmetric_sparse_matrix(N*nb_dof, N*nb_dof);
  linalg_vector ST1, ST2;

  cout << "nombre de ddl de l'element fini : " << nb_dof << endl;
  cout << "nombre de ddl pour l'elasticite lineaire : " << nb_dof * N << endl;

  cout << "Assemblage de la matrice de rigidite" << endl;
  ST1 = linalg_vector(nb_dof_data2); ST2 = linalg_vector(nb_dof_data2);
  std::fill(ST1.begin(), ST1.end(), lambda);
  std::fill(ST2.begin(), ST2.end(), G);
  assembling_rigidity_matrix_for_linear_elasticity(RM,mef,mef_data2,ST1,ST2);

  cout << "Assemblage du terme source" << endl;
  ST1 = linalg_vector(nb_dof_data * N);
  for (size_type i = 0; i < nb_dof_data; ++i)
    for (size_type j = 0; j < N; ++j) 
      ST1[i*N+j] = sol_f(mef_data.point_of_dof(i))[j];
  getfem::assembling_volumic_source_term(B, mef, mef_data, ST1, N);

  cout << "Assemblage de la condition de Neumann" << endl;
  ST1 = linalg_vector(nb_dof_data * N);
  getfem::base_node pt(N); getfem::base_vector n(N), v;
  for (size_type i = 0; i < nb_dof_data; ++i)
  {
    pt = mef_data.point_of_dof(i);
    if (dal::abs(pt[0]-LX) < 10E-6) n[0] = 1.0; // pas terrible ... !!
      else if (dal::abs(pt[0]) < 10E-6) n[0] = -1.0; else n[0] = 0.0;
    if (N > 1) {
      if (dal::abs(pt[1]-LY) < 10E-6) n[1] = 1.0;
      else if (dal::abs(pt[1]) < 10E-6) n[1] = -1.0; else n[1] = 0.0;
      for (int k = 2; k < N; ++k)
	if (dal::abs(pt[k]-LZ) < 10E-6) n[k] = 1.0;
	else if (dal::abs(pt[k]) < 10E-6) n[k] = -1.0; else n[k] = 0.0;
    }
    v = sol_sigma(pt) *  n;
    for (size_type j = 0; j < N; ++j)
      ST1[i*N+j] = v[j];
  }
  getfem::assembling_Neumann_condition(B, mef, 1, mef_data, ST1, N);

  cout << "Prise en compte de la condition de Dirichlet" << endl;
  ST1 = linalg_vector(nb_dof * N);
  for (size_type i = 0; i < nb_dof; ++i)
    for (size_type j = 0; j < N; ++j)
      ST1[i*N+j] = sol_u(mef.point_of_dof(i))[j];
  getfem::assembling_Dirichlet_condition(RM, B, mef, 0, ST1, N);
}

void pb_data::solve(void)
{ precond_cg(RM, U, B, 20000, residu); }

void pb_data::visu(void)
{
  dal::dynamic_array<base_node> ptab;
  dal::bit_vector nn = mef.convex_index();
  int c;
  for (int i = 0; i < mesh.points().index().last()+1; ++i)
  { ptab[i] = base_node(N); ptab[i].fill(0.0); }

  for (c << nn; c >= 0; c << nn)
  {
    for (int i = 0; i < mef.nb_dof_of_element(c); ++i)
    {
      int l = mef.ind_dof_of_element(c)[i];
      int k = mesh.search_point(mef.point_of_dof(l));
      if (k >= 0)
      {
	for (int j = 0; j < N; ++j) ptab[k][j] = U[l * N + j];
	ptab[k] *= muldep;
      }
    }
  }
  simple_mesh_visu(mesh, "Configuration deformee", NX <= 5, ptab);
}

int main(int argc, char *argv[])
{
  pb_data p;

  cout << "initialisation ...\n";
  p.PBSTFR_PARAM.read_command_line(argc, argv);
  p.init();
  cout << "Initialisation terminee\n";
  p.mesh.stat();
  p.mesh.write_to_file(p.datafilename + ".mesh" + char(0));
 
  cout << "Debut de l'assemblage\n";
  p.assemble();
  cout << "Assemblage termine\n";

//   cout << "Matrice de rigidite\n";
//   for (int i = 0; i < p.RM.nrows(); i++)
//   { 
//     cout << "ligne " << i << " [ ";
//     for (l = 0; l < p.RM.nrows(); l++)
//       if (p.RM(i, l) != 0.0)
// 	cout << "(" << l << "," << p.RM(i, l) << ")  ";
//     cout << "]" << endl;
//   }
//   cout << endl << endl;

  cout << "Resolution\n";
  p.solve();

  cout << "Calcul de l'erreur\n";

  int nbdof = p.mef.nb_dof();
  linalg_vector V(nbdof*p.N); mtl::copy(p.U, V);
  base_vector S;

  for (int i = 0; i < nbdof; ++i)
  {
    S = sol_u(p.mef.point_of_dof(i));
    for (int k = 0; k < p.N; ++k) V[i*p.N + k] -= S[k];
  }

  cout <<  "Error L^2 " << getfem::L2_norm(p.mef, V, p.N) << endl;
  cout <<  "Error H^1 " << getfem::H1_norm(p.mef, V, p.N) << endl;

  if (p.draw) p.visu();

  return 0;
}
