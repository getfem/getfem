
/**************************************************************************/
/*                                                                        */
/*  Probleme du laplacien.                                                */
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

  scalar_type LX, LY, LZ, residu;
  int NX, N, K;

  symmetric_sparse_matrix RM; /* matrice de rigidite.                     */
  linalg_vector U, B; /* inconnue et second membre.                       */
 
  bool mixte;
  int iteimpl, integration;

  int NBLEVEL;
  bool draw;
  std::string datafilename;
  scalar_type muldep;
  md_param PARAM;

  void assemble(void);
  void solve(void);
  void init(void);
  lap_pb(void) : mef(mesh), mef_data(mesh), mef_data2(mesh) {}
};

void lap_pb::init(void)
{
  dal::bit_vector nn;
  int i, j, k, h;

  /***********************************************************************/
  /*  LECTURE DES PARAMETRES SUR FICHIER.                                */
  /***********************************************************************/
  
  /* parametres physiques */
  N = PARAM.int_value("N", "Dimension de l'espace");
  /* parametres numeriques */
  LX = PARAM.real_value("LX", "Taille en X");
  LY = PARAM.real_value("LY", "Taille en Y");
  LZ = PARAM.real_value("LZ", "Taille en Y");
  NX = PARAM.int_value("NX", "Nombre de pas d'espace ");
  iteimpl = PARAM.int_value("ITEIMPL", "Nombre d'iteration pt fixe");
  integration = PARAM.int_value("INTEGRATION", "Type d'integration");
  residu = PARAM.real_value("RESIDU", "Valeur pour test d'arret");
  muldep = PARAM.real_value("MULDEP", "Coeff de grossissement");
  NBLEVEL = PARAM.int_value("NBLEVEL", "Nombre de courbes de niveau");
  K = PARAM.int_value("K", "Degre de l'element fini de Lagrange");
  datafilename = std::string( PARAM.string_value("ROOTFILENAME",
			     "Nom du fichier de sauvegarde sans extension"));
  char *dds = PARAM.string_value("GRAPHICS", "Dessin ? ");
  draw = (strcmp("N", dds) && strcmp("n", dds));

  scalar_type FT = PARAM.real_value("FT", "parametre pour la solution exacte");

  dds = PARAM.string_value("MIXTE", "Methode de Leila ? ");
  mixte = (strcmp("N", dds) && strcmp("n", dds));

  sol_K = base_vector(N);
  for (j = 0; j < N; j++)
    sol_K[j] = ((j & 1) == 0) ? FT : -FT;

  cout << "\n\n";

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

  // mesh.write_to_file(cout);

  // if (draw) simple_mesh_visu(mesh, "Maillage", true);

  cout << "Choix de l'element fini.\n";
  getfem::pintegration_method ppi;
  nn = mesh.convex_index(N);
  if (mixte) {
    if (integration == 0) ppi = bgeot::simplex_poly_integration(N);
    else ppi = bgeot::Newton_Cotes_approx_integration(N,2);
    mef.set_finite_element(nn, getfem::P1_nonconforming_fem(), ppi);
    assert(N == 2);
  }
  else {
    if (integration == 0) ppi = bgeot::simplex_poly_integration(N);
    else ppi = bgeot::Newton_Cotes_approx_integration(N,2*K);
    mef.set_finite_element(nn, getfem::PK_fem(N, K), ppi);
  }

  cout << "Nombre d'elements : " << nn.card() << endl;
  mef_data.set_finite_element(nn, getfem::PK_fem(N, K),
			      bgeot::simplex_poly_integration(N));
  mef_data2.set_finite_element(nn, getfem::PK_fem(N, 0),
			       bgeot::simplex_poly_integration(N));

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
	 if (std::abs(un[N-1] - 1.0) < 1.0E-3)
	{
	  mef.add_boundary_elt(0, j, i);
// 	  cout << "ajout d'un bord Dirichlet : convexe\t" << j << "\tface "
// 	       << i << endl;

	}
	else
	{
	  mef.add_boundary_elt(1, j, i);
//	  cout << "ajout d'un bord Neumann   : convexe\t" << j << "\tface "
// 	       << i << endl;
	}
// 	  }
      }
    }
  }
}

void lap_pb::assemble(void)
{
  int nb_dof = mef.nb_dof(), nb_dof_data = mef_data.nb_dof();
  int nb_dof_data2 = mef_data2.nb_dof();
  B = linalg_vector(nb_dof, 0.0);
  U = linalg_vector(nb_dof, 0.0);
  RM = symmetric_sparse_matrix(nb_dof, nb_dof);
  linalg_vector ST;
  
  cout << "nombre de ddl de l'element fini : " << nb_dof << endl;
  cout << "nombre de ddl de l'element pour les données : " << nb_dof_data
       << endl;

  cout << "Assemblage de la matrice de rigidite" << endl;
  ST = linalg_vector(nb_dof_data2);
  std::fill(ST.begin(), ST.end(), 1.0);
  getfem::assembling_rigidity_matrix_for_laplacian(RM, mef, mef_data2, ST);
  
  cout << "Assemblage du terme source" << endl;
  ST = linalg_vector(nb_dof_data);
  for (size_type i = 0; i < nb_dof_data; ++i)
    ST[i] = sol_f(mef_data.point_of_dof(i));
  getfem::assembling_volumic_source_term(B, mef, mef_data, ST, 1);

  cout << "Assemblage de la condition de Neumann" << endl;
  ST = linalg_vector(nb_dof_data);
  getfem::base_node pt(N); getfem::base_vector n(N);
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
    ST[i] = bgeot::vect_sp(sol_grad(pt), n);
  }
  getfem::assembling_Neumann_condition(B, mef, 1, mef_data, ST, 1);
  
  cout << "Prise en compte de la condition de Dirichlet" << endl;
  dal::bit_vector nn = mef.dof_on_boundary(0);
  cout << "Number of Dirichlet nodes : " << nn.card() << endl;
  cout << "Dirichlet nodes : " << nn << endl;
  ST = linalg_vector(nb_dof);
  for (size_type i = 0; i < nb_dof; ++i)
    ST[i] = sol_u(mef.point_of_dof(i));
  getfem::assembling_Dirichlet_condition(RM, B, mef, 0, ST, 1);
}

void lap_pb::solve(void)
{
  precond_cg(RM, U, B, 20000, residu);
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[])
{
  cout << "debut du programme" << endl;

  lap_pb p;
  dal::bit_vector nn;
  int l, i;
  scalar_type exectime, total_time = 0.0;

  exectime = uclock_sec();

  cout << "initialisation ...\n";
  p.PARAM.read_command_line(argc, argv);
  p.init();
  cout << "Initialisation terminee\n";

  ofstream cres((p.datafilename + ".res" + char(0)).begin());

  cres << p.N << "\t" <<  p.K << "\t" << p.NX << "\t";
  cres << uclock_sec() - exectime << "  ";

  total_time += uclock_sec() - exectime;

  // p.mesh.write_to_file(cout);
  // p.mesh.stat();
  
  // p.mesh.write_to_file(p.datafilename + ".mesh" + char(0));
  if (p.draw) simple_mesh_visu(p.mesh, "Maillage", true);
 
  exectime = uclock_sec();
  int nb_dof = p.mef.nb_dof();

  total_time += uclock_sec() - exectime;

  cres << nb_dof << "\t" <<  uclock_sec() - exectime << "\t";

  cout << "Debut de l'assemblage\n";
  exectime = uclock_sec();
  p.assemble();

  cres << uclock_sec() - exectime << "\t";
  total_time += uclock_sec() - exectime;

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
  
  cout << "Assemblage termine\n";
  // mtl::print_row(p.RM);
  
//    cout << "Second membre : " << endl;
//    mtl::print_vector(p.B);

  cout << "Resolution\n";
  exectime = uclock_sec();
  p.solve();

  cres << uclock_sec() - exectime << "\t";
  total_time += uclock_sec() - exectime;
  exectime = uclock_sec();

  int nbdof = p.mef.nb_dof();
  linalg_vector V(nbdof); mtl::copy(p.U, V);
  for (int i = 0; i < nbdof; ++i)
    V[i] -= sol_u(p.mef.point_of_dof(i));

  scalar_type l2norm = getfem::L2_norm(p.mef, V, 1);
  cres << l2norm << "\t";
  cres << uclock_sec() - exectime << "\t";
  total_time += uclock_sec() - exectime;
  exectime = uclock_sec();
  scalar_type h1norm = getfem::H1_norm(p.mef, V, 1);
  cres << h1norm << "\t";
  cres << uclock_sec() - exectime << "\t";
  total_time += uclock_sec() - exectime;
  cres << total_time << endl;

  cout << "Erreur L2 : " << l2norm << endl << "Erreur H1 : " << h1norm << endl;

  cout << "calcul termine" << endl; exit(0);
  return 0; 
}
