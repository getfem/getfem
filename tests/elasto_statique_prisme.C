/**************************************************************************/
/*                                                                        */
/*  Calcul sur une structure lineairement elastique                       */
/*                                                                        */
/**************************************************************************/

#include <getfem_applied_fem.h>
#include <getfem_regular_meshes.h>

#include <linear_systems.h>
#include <model_param.h>
#include <graph_visu_3D.h>
#include <geo_mesh_visu_3D.h>
#include <clock.h>


typedef mtl::matrix<double, mtl::symmetric<mtl::upper>, 
                            mtl::array< mtl::compressed<> >,
                            mtl::row_major >::type sparse_matrix;

typedef mtl::dense1D<double> linalg_vector;

typedef bgeot::vsvector<double> plain_vector;
typedef bgeot::vsmatrix<double> plain_matrix;
typedef bgeot::PT< bgeot::vsvector<double> > node;

md_param PBSTFR_PARAM;

/**************************************************************************/
/*  Definition de la solution test.                                       */
/**************************************************************************/

dal::dynamic_array<plain_vector> sol_K;
double sol_lambda, sol_G;

plain_vector sol_u(const node &x)
{
  int N = x.size();plain_vector res(N);
  for (int i = 0; i < N; ++i) res[i] = sin(bgeot::vect_sp(sol_K[i], x));
  return res;
}

struct Dirichlet_term : public getfem::vector_function_by_node<plain_vector>
{
  const getfem::mesh_fem *mf;
  inline plain_vector operator() (int i) const
  { return sol_u(mf->point_of_dof(i)); }
  Dirichlet_term(const getfem::mesh_fem &m) { mf = &m; }
};

plain_vector sol_f(const node &x)
{
  int N = x.size();
  plain_vector res(N);
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

struct source_term : public getfem::vector_function_by_node<plain_vector>
{
  const getfem::mesh_fem *mf;
  inline plain_vector operator() (int i) const
  { return sol_f(mf->point_of_dof(i)); }
  source_term(const getfem::mesh_fem &m) { mf = &m; }
};

plain_matrix sol_sigma(const node &x)
{
  int N = x.size();
  plain_matrix res(N,N);
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

struct Neumann_term_n : public getfem::matrix_function_by_node<plain_matrix>
{
  const getfem::mesh_fem *mf;
  inline plain_matrix operator() (int i) const
  { return sol_sigma(mf->point_of_dof(i)); }
  Neumann_term_n(const getfem::mesh_fem &m) { mf = &m; }
};

struct Neumann_term_e : public getfem::matrix_function_by_element<plain_matrix>
{
  const getfem::mesh_fem *mf;
  inline plain_matrix operator() (int i) const
  {
    plain_matrix m = sol_sigma(mf->linked_mesh().points_of_convex(i)[0]);
    int n = mf->linked_mesh().nb_points_of_convex(i);
    for (int k = 1; k < n; ++k)
      m += sol_sigma(mf->linked_mesh().points_of_convex(i)[k]);
    m *= (1.0 / double(n));
    return m;
  }
  Neumann_term_e(const getfem::mesh_fem &m) { mf = &m; }
};

struct data_pb
{
  getfem::getfem_mesh mesh;
  getfem::mesh_fem mesh_fem;

  double G, lambda;
  double LX, LY, LZ, errgauss;
  int NX, N;
  int K;     /* finite element degree.                                    */


  sparse_matrix RM; /* matrice de rigidite.                               */
  linalg_vector U, B; /* inconnue et second membre.                       */
  dal::bit_vector ddlD, ddlN;  /* ddl de Dirichlet et de Neumann.            */

  linalg_vector cDi;
  linalg_vector cNe;
 

  int iteimpl;

  void assemble(void);
  void solve(void);

  data_pb(void) : mesh_fem(mesh) {}
};

void data_pb::assemble(void)
{
  int nb_dof = mesh_fem.nb_dof();
  B = linalg_vector(N*nb_dof, 0.0);
  U = linalg_vector(N*nb_dof, 0.0);
  RM = sparse_matrix(N*nb_dof, N*nb_dof);
  
  cout << "nombre de ddl de l'element fini : " << nb_dof << endl;
  cout << "nombre de ddl pour l'elasticite lineaire : " << nb_dof * N << endl;

  cout << "Assemblage de la matrice de rigidite" << endl;
  getfem::rigidity_matrix_for_linear_elasticity(RM, lambda, G, mesh_fem);

  cout << "Assemblage du terme source" << endl;
  getfem::volumic_source_term(B, source_term(mesh_fem), mesh_fem);
  
  cout << "Assemblage de la condition de Neumann" << endl;
  getfem::Neumann_condition(B, Neumann_term_n(mesh_fem), mesh_fem, 1);

  cout << "Prise en compte de la condition de Dirichlet" << endl;
  getfem::Dirichlet_condition(RM, B, Dirichlet_term(mesh_fem), mesh_fem, 0);
  
}

void data_pb::solve(void)
{
  itl::cholesky<sparse_matrix> precond(RM);
  // itl::SSOR<sparse_matrix> precond(RM);
  itl::noisy_iteration<double> iter(B, 100000, errgauss);
  itl::cg(RM, U, B, precond(), iter);
}


int NBLEVEL;
bool draw_solutions;
string datafilename;

double muldep;

void init_pbstfr(data_pb &p)
{
  dal::bit_vector nn;
  int i, j, k, h;

  /***********************************************************************/
  /*  LECTURE DES PARAMETRES SUR FICHIER.                                */
  /***********************************************************************/
  
  /* parametres physiques */
  p.N = PBSTFR_PARAM.int_value("N", "Dimension de l'espace");
  p.G = PBSTFR_PARAM.real_value("G", "Raideur de cisaillement");
  p.lambda = PBSTFR_PARAM.real_value("LAMBDA", "Coeff elastique lambda");
 
  /* parametres numeriques */
  p.LX = PBSTFR_PARAM.real_value("LX", "Taille en X");
  p.LY = PBSTFR_PARAM.real_value("LY", "Taille en Y");
  p.LZ = PBSTFR_PARAM.real_value("LZ", "Taille en Y");
  p.NX = PBSTFR_PARAM.int_value("NX", "Nombre de pas d'espace ");
  p.iteimpl = PBSTFR_PARAM.int_value("ITEIMPL", "Nombre d'iteration pt fixe");
  p.errgauss = PBSTFR_PARAM.real_value("RESIDU", "Valeur pour test d'arret");
  muldep = PBSTFR_PARAM.real_value("MULDEP", "Coeff de grossissement");
  
  NBLEVEL = PBSTFR_PARAM.int_value("NBLEVEL", "Nombre de courbes de niveau");
  p.K = PBSTFR_PARAM.int_value("K", "Degre de l'element fini de Lagrange");
  datafilename = string(PBSTFR_PARAM.string_value("ROOTFILENAME",
			     "Nom du fichier de sauvegarde sans extension"));
  char *dds = PBSTFR_PARAM.string_value("GRAPHICS", "Dessin ? ");
  draw_solutions = !(!strcmp("N", dds) || !strcmp("n", dds));

  double FT;
  FT = PBSTFR_PARAM.real_value("FT", "parametre pour la solution exacte");

  for (i = 0; i < p.N; i++)
  {
    sol_K[i] = plain_vector(p.N);
    for (j = 0; j < p.N; j++)
      sol_K[i][j] = (i == j) ? FT : -FT;
  }

  sol_lambda = p.lambda;
  sol_G = p.G;

  /***********************************************************************/
  /*  CONSTRUCTION DU MAILLAGE.                                          */
  /***********************************************************************/

  cout << "construction du maillage\n";

  node org(p.N-1); org.fill(0.0);

  plain_vector vtab[10];
  int ref[30];
  assert(p.N < 10);

  for (i = 0; i < p.N-1; i++)
  { 
    double a = p.LZ / double(p.NX); int n = p.NX;
    if (i == 0) { a = p.LX / double(p.NX); n = p.NX; }
    if (i == 1) { a = p.LY / double(p.NX); n = p.NX; }
    vtab[i] = plain_vector(p.N-1);
    vtab[i].fill(0.0); (vtab[i])[i] = a; ref[i] = n;
  }

  getfem::getfem_mesh little_mesh;

  parallelepiped_regular_simplex_mesh(little_mesh, p.N - 1, org, vtab, ref);

  if (draw_solutions)
    simple_mesh_visu(little_mesh, "Petit Maillage", p.NX <= 5);

  nn = little_mesh.convex_index();
  for (j << nn; j >= 0; j << nn)
  {
    for (int i = 0; i < p.NX; ++i)
    {
      for (int l = 0; l < p.N; ++l)
      {
	node pt(p.N);
	for (int k = 0; k < p.N-1; ++k)
	  pt[k] = (little_mesh.points_of_convex(j)[l])[k];
	pt[p.N-1] = i * p.LZ / double(p.NX);
	ref[l] = p.mesh.add_point(pt);
	pt[p.N-1] = (i+1) * p.LZ / double(p.NX);
	ref[l+p.N] = p.mesh.add_point(pt);
      }  
      p.mesh.add_prism(p.N, ref);
    }
  }

  cout << "Optimisation de la structure des points\n";
  p.mesh.optimize_structure();

  cout << "Choix de l'element fini.\n";
  nn = p.mesh.convex_index(p.N);
  cout << "Nombre d'elements : " << nn.card() << endl;
  p.mesh_fem.set_finite_element(nn,
        getfem::fem_exact_integration(
	 bgeot::product_interpolation( bgeot::PK_fem_interpolation(p.N-1, p.K),
	                               bgeot::PK_fem_interpolation(1, p.K)),
	 bgeot::linear_product_trans(  bgeot::simplex_trans(p.N-1, 1),
				       bgeot::simplex_trans(1, 1)),
	 bgeot::prism_poly_integration(p.N)));

  if (draw_solutions) simple_mesh_visu(p.mesh, "Grand maillage", (p.NX <= 5));

  cout << "nbdof : " << p.mesh_fem.nb_dof() << endl;

  cout << "Reperage des bord Neumann et Dirichlet\n";
  nn = p.mesh.convex_index(p.N);
  base_vector un;
  for (j << nn; j >= 0; j << nn)
  {
    k = p.mesh.structure_of_convex(j)->nb_faces();
    for (i = 0; i < k; i++)
    {
      if (bgeot::neighbour_of_convex(p.mesh, j, i).empty())
      {
	un = p.mesh.convex(j).unit_norm_of_face(i);
	
	// if (true)
	if (abs(un[p.N-1] - 1.0) < 1.0E-3)
	  p.mesh_fem.add_boundary_elt(0, j, i);
	else
	  p.mesh_fem.add_boundary_elt(1, j, i);
      }
    }
  }
  
}

void visu_result(data_pb &p)
{
  dal::dynamic_array<node> ptab;
  for (int i = 0; i < p.mesh.points().index().last(); ++i)
  { ptab[i] = node(p.N); ptab[i].fill(0.0); }

  dal::bit_vector nn = p.mesh_fem.convex_index();
  int c;
  for (c << nn; c >= 0; c << nn)
  {
    for (int i = 0; i < p.mesh_fem.nb_dof_of_element(c); ++i)
    {
      int k = p.mesh.search_point(
	      p.mesh_fem.point_of_dof(p.mesh_fem.ind_dof_of_element(c)[i]));
      if (k >= 0)
      {
	for (int j = 0; j < p.N; ++j)
	  ptab[k][j] = p.U[p.mesh_fem.ind_dof_of_element(c)[i] * p.N + j];
	ptab[k] *= muldep;
      }
    }
  }

  simple_mesh_visu(p.mesh, "Configuration deformee", (p.NX <= 5), ptab);
}

int main(int argc, char *argv[])
{

  data_pb p;
  dal::bit_vector nn;
  int l, i;


  cout << "initialisation ...\n";
  PBSTFR_PARAM.read_command_line(argc, argv);
  init_pbstfr(p);
  cout << "Initialisation terminee\n";
  
  p.mesh.stat();
  cout << "nbdof : " << p.mesh_fem.nb_dof() << endl;
  // p.mesh.write_to_file(cout);
  string meshout = datafilename; meshout += ".mesh"; meshout += char(0);
  p.mesh.write_to_file(meshout);
 
  cout << "Debut de l'assemblage\n";
  p.assemble();

  cout << "Matrice de rigidite\n";
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
  
  // cout << "Second membre : " << endl;
  // mtl::print_vector(p.B);

  cout << "Resolution\n";
  p.solve();

//   cout << "Vecteur solution calculee : " << endl;
//   mtl::print_vector(p.U);

//   int nbd = p.mesh_fem.nb_dof();
//   cout << "Vecteur solution exacte : " << endl;
//   for (i = 0; i < nbd; i++)
//   {
//     plain_vector V = sol_u(p.mesh_fem.point_of_dof(i));
//     for (int k = 0; k < p.N; k++) cout << V[k] << " ";
//   }
//   cout << endl;


 cout << "Calcul de l'erreur\n";

  int nbdof = p.mesh_fem.nb_dof();
  linalg_vector V(nbdof*p.N); mtl::copy(p.U, V);
  plain_vector S;

  for (int i = 0; i < nbdof; ++i)
  {
    S = sol_u(p.mesh_fem.point_of_dof(i));
    for (int k = 0; k < p.N; ++k) V[i*p.N + k] -= S[k];
  }

  cout <<  "Erreur  L^2 " << getfem::L2_norm(p.mesh_fem, V, p.N) << endl;
  cout <<  "Erreur  H^1 " << getfem::H1_norm(p.mesh_fem, V, p.N) << endl;



  if (draw_solutions) visu_result(p);
  
  return 0;
}
