
#include <getfem_export.h>

using std::cout;
using std::cerr;

using bgeot::size_type;

void err_msg(void)
{
  cerr << "Bad format for arguments of command draw_mesh\n";
  cerr << "The right format is\n";
  cerr << "compare_solutions filename1\n";
  cerr << "where filename1 is files containing a mesh\n";
  cerr << "exported by GETFEM++.\n";
  exit(1);
}


int main(int argc, char *argv[])
{
  try {
    int found = 0;
    std::string fi1;

    for (int aa = 1; aa < argc; aa++) {
      if (argv[aa][0] != '-') {
	switch(found) {
        case 0  : fi1 = std::string(argv[aa]); found++; break;
        default : err_msg();
	}
      }
      else
	err_msg();
    }
    if (found != 3) err_msg();
    
    cout.precision(16);
    
    getfem::getfem_mesh mesh;
    mesh.read_from_file(fi1);
    getfem::edge_list el;
    getfem::mesh_edges_list(mesh, el);
    dal::bit_vector nn = el.index();
    size_type i;
    for (i << nn; i != size_type(-1); i << nn) {
      cout << "edge : " << el[i].i << " : " << el[i].j << endl;
    }

  
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
}
