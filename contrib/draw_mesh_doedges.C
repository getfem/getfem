#include <getfem_export.h>

using bgeot::size_type;

void err_msg(void)
{
  cerr << "Bad format for arguments of command draw_mesh. The right format is"
       << " 'draw_mesh filename1' where filename1 is a file containing a mesh"
       << " exported by GETFEM++.\n"; exit(1);
}

int main(int argc, char *argv[])
{
  try {
    int found = 0;
    std::string fi1;
    for (int aa = 1; aa < argc; aa++) {
      if (argv[aa][0] != '-') { fi1 = std::string(argv[aa]); found++; }
      else err_msg();
    }
    if (found != 1) err_msg();
    cout.precision(16);
    getfem::getfem_mesh mesh;
    mesh.read_from_file(fi1);
    int N = mesh.dim();
    cout << "N = " << N << endl;
    if (N <= 0 || N >= 4)
      DAL_THROW(dal::internal_error,
	     "Do not draw meshes in dimension lower than 1 or greater than 3");
    getfem::edge_list el;
    getfem::mesh_edges_list(mesh, el);
    dal::bit_vector nn = el.index();
    size_type i;
    for (i << nn; i != size_type(-1); i << nn) {
      for (int k = 0; k < N; ++k) cout << mesh.points()[el[i].i][k] << "\t";
      if (N == 1) cout << "0.0"; cout << endl;
      for (int k = 0; k < N; ++k) cout << mesh.points()[el[i].j][k] << "\t";
      if (N == 1) cout << "0.0"; cout << endl << endl;
    }
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
}
