#include <getfem_export.h>
#include <getfem_assembling.h>

using std::cout;
using std::cerr;

typedef getfem::size_type size_type;

void err_msg(void)
{
  cerr << "Bad format for arguments of command interpole_solution\n";
  cerr << "The right format is\n";
  cerr << "interpole_solution filename1 filename2 filename3\n";
  cerr << "where filename1 and filename2 are two files containing some\n";
  cerr << "solutions exported by GETFEM++. filename3 is the name of the result\n";
  exit(1);
}

void bad_format(const std::string &fi)
{
  cerr << "Data file " << fi << " is corrupted or has a bad format\n";
  exit(1);
}

void read_file(const std::string &fi, getfem::getfem_mesh &mesh,
	       getfem::mesh_fem &mef, std::vector<getfem::scalar_type> &U, 
	       size_type &N, size_type &P, size_type &K)
{
  std::ifstream ist((fi + char(0)).data());
  size_type i;

  cout << "reading file " << fi << endl;

  if (ist == NULL)
  { cerr << "At least one of the file doesn't exist\n"; exit(1); }

  char tmp[100], c;

  if (!(ftool::read_untill(ist, "DATA BY ELEMENT"))) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  if (strcmp(tmp, "N")) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  if (strcmp(tmp, "=")) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  N = atoi(tmp);
  if (N > 255 || N == 0) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  if (strcmp(tmp, "P")) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  if (strcmp(tmp, "=")) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  P = atoi(tmp);
  if (P > 255 || P == 0) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  if (strcmp(tmp, "K")) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  if (strcmp(tmp, "=")) bad_format(fi);
  ftool::get_token(ist, tmp, 99);
  K = atoi(tmp);
  if (K > 255) bad_format(fi);

  cout << "N = " << N <<  "  P = " << P << "  K = " << K << endl;
  
  dal::dynamic_array<getfem::base_node> ptab;
  dal::dynamic_array<getfem::scalar_type> vtab;
  size_type nbvtab = 0;

  size_type Np = 1; for (i = 0; i < N; ++i) Np *= K + 1;
 
  // prise en compte des solutions sur des simplexes de dim inferieure au maillage ...

  while(!(ist.eof()))
  {
    size_type nbpt = 0;
    while(!(ist.eof()) && nbpt < 100)
    {
      if (ptab[nbpt].size() != N) ptab[nbpt].resize(N);
      for (i = 0; i < N; ++i) ist >> (ptab[nbpt])[i];
      if (ist.eof()) break;
      for (i = 0; i < P; ++i) ist >> vtab[nbvtab * P + i];
      nbpt++; nbvtab++;
      // cout << "nbvtab = " << nbvtab << " memsize = " << vtab.memsize() << endl;
      while(ist.get(c))
	{ if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
      ist.get(c); if (c != 10 && !(ist.eof())) bad_format(fi);
      while(ist.get(c))
	{ if (!(isspace(c)) || c == 10) { ist.putback(c); break; } }
      ist.get(c); if (c == 10 || ist.eof()) break; else ist.putback(c);
    }
    if (nbpt == 100) bad_format(fi);
    if (nbpt == 0) break;
    
    if (nbpt == bgeot::alpha(N, K))
    {
      std::vector<getfem::base_node> pnode(N+1);
      for (i = 0; i <= N; ++i) pnode[i] = ptab[bgeot::alpha(i, K) - 1];
      i = mesh.add_simplex_by_points(N, pnode.begin());
      mef.set_finite_element(i, getfem::PK_fem(N, K),
			     bgeot::simplex_poly_integration(N));
    }
    else if (nbpt == Np) 
    { // a verifier
      bgeot::pgeometric_trans pgt = bgeot::simplex_trans(1, K);
      for (i = 1; i < N; ++i) // à simplifier !!
        pgt = bgeot::convex_product_trans(pgt, bgeot::simplex_trans(1, K));
      
      std::vector<getfem::base_node> pnode(1 << N);
      size_type j;
      for (i = 0, j = 0; i <= (1 << N); ++i)
      {
	pnode[i] = ptab[j];
	size_type k = i + 1, l = K-1;
	while (!(k & 1)) { l *= K+1; k >>= 1; }
	j += l + 1;
      }
      i = mesh.add_parallelepiped_by_points(N, pnode.begin());
      mef.set_finite_element(i, getfem::QK_fem(N, K),
			     bgeot::parallelepiped_poly_integration(N));
    }
    else if (nbpt == bgeot::alpha(N - 1, K) * bgeot::alpha(1, K))
    { // a verifier
      std::vector<getfem::base_node> pnode(2*N);
      for (i = 0; i < N; ++i)
	pnode[i] = ptab[bgeot::alpha(i, K) - 1];
      for (i = 0; i < N; ++i)
	pnode[i+N] = ptab[bgeot::alpha(i, K) - 1 + bgeot::alpha(N-1, K)];
      i = mesh.add_prism_by_points(N, pnode.begin());
      mef.set_finite_element(i, getfem::PK_prism_fem(N, K),
			     bgeot::prism_poly_integration(N));
    }
    else 
    { cerr << "Unknown element in data file\n"; exit(1); }
  }

  ist.close();

  cout << "Data repartition\n";
  
  dal::bit_vector nn = mesh.convex_index();
  U.resize(mef.nb_dof() * P);
  std::fill(U.begin(), U.end(), 0.0);
  std::vector<int> cp(mef.nb_dof());
  std::fill(cp.begin(), cp.end(), 0);
  size_type l = 0;
  for (i << nn; i != size_type(-1); i << nn)
  {
    size_type nbd = mef.nb_dof_of_element(i);
    for (int j = 0; j < nbd; ++j, ++l)
    {
      size_type dof = mef.ind_dof_of_element(i)[j];
      for (int k = 0; k < P; ++k) U[dof*P +k] += vtab[l*P+k];
      (cp[dof])++;
    }
  }
  for (i = 0; i < mef.nb_dof(); ++i)
  {
    if (cp[i] == 0) { cerr << "Internal error\n"; exit(1); }
    for (int k = 0; k < P; ++k) U[i*P +k] /= getfem::scalar_type(cp[i]);
  }
   
}


int main(int argc, char *argv[])
{
  try {
  int found = 0;
  std::string fi1, fi2, fi3;

  for (int aa = 1; aa < argc; aa++)
  {
    if (argv[aa][0] != '-')
    {
      switch(found)
      {
        case 0  : fi1 = std::string(argv[aa]); found++; break;
        case 1  : fi2 = std::string(argv[aa]); found++; break;
        case 2  : fi3 = std::string(argv[aa]); found++; break;
        default : err_msg();
      }
    }
    else
      err_msg();
  }
  if (found != 3) err_msg();

  cout.precision(14);

  getfem::getfem_mesh mesh1;
  getfem::mesh_fem mef1(mesh1);
  size_type N1, P1, K1;
  std::vector<getfem::scalar_type> U1;
  read_file(fi1, mesh1, mef1, U1, N1, P1, K1);
  
  getfem::getfem_mesh mesh2;
  getfem::mesh_fem mef2(mesh2);
  size_type N2, P2, K2;
  std::vector<getfem::scalar_type> U2, U3;
  read_file(fi2, mesh2, mef2, U2, N2, P2, K2);

  if (N1 != N2)
  { cerr << "Dimensions of the two meshes don't match\n"; exit(1); }

  if (P1 != P2)
  { cerr << "Dimensions of the two solutions don't match\n"; exit(1); }


  cout << "interpolation of the solution in " << fi1
       << " on the mesh of " << fi2 << endl;
  U3.resize(mef2.nb_dof() * P1);
  getfem::scalar_type errin = 0.0;
  getfem::interpolation_solution(mef1, mef2, U1, U3, P1);
  save_solution(fi3, mef2, U3, P1, K2);
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
}
