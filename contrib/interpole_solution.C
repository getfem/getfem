#include <getfem_export.h>
#include <getfem_assembling.h>

using bgeot::size_type;
using bgeot::dim_type;
using bgeot::short_type;

void err_msg(void)
{
  cerr << "Bad format for arguments of command interpole_solution\n";
  cerr << "The right format is\n";
  cerr << "interpole_solution filename1 filename2 filename3\n";
  cerr << "where filename1 and filename2 are two files containing some\n";
  cerr << "solutions exported by GETFEM++. filename3 is the name of the ";
  cerr << "result\n";
  exit(1);
}


int main(int argc, char *argv[])
{
  try {
    int found = 0;
    std::string fi1, fi2, fi3;
    
    for (int aa = 1; aa < argc; aa++) {
      if (argv[aa][0] != '-') {
	switch(found) {
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
    
    cout << "Reading file " << fi1 << endl;
    getfem::getfem_mesh mesh1;
    getfem::mesh_fem mef1(mesh1);
    dim_type N1, P1; short_type K1;
    std::vector<getfem::scalar_type> U1;
    getfem::load_solution(fi1, mesh1, mef1, U1, K1);
    N1 = mesh1.dim();
    P1 = mef1.get_qdim();
    
    cout << "Reading file " << fi2 << endl;
    getfem::getfem_mesh mesh2;
    getfem::mesh_fem mef2(mesh2);
    dim_type N2; short_type K2;
    std::vector<getfem::scalar_type> U2, U3;
    getfem::load_solution(fi2, mesh2, mef2, U2, K2);
    N2 = mesh1.dim();
  
    if (N1 != N2) DAL_THROW(std::invalid_argument,
			    "Dimensions of the two meshes mismatch\n");
    

    cout << "interpolation of the solution in " << fi1
	 << " on the mesh of " << fi2 << endl;
    U3.resize(mef2.nb_dof() * P1);
    getfem::interpolation_solution(mef1, mef2, U1, U3, P1);
    getfem::save_solution(fi3, mef2, U3, K2);
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
}
