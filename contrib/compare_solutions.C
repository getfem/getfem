#include <getfem_export.h>
#include <getfem_assembling.h>

using bgeot::size_type;
using bgeot::dim_type;
using bgeot::short_type;

/* ugly remplacement for "pre-qdim era" interpolation
 * ( contrib/compare_solutions uses that.. ) */
template<class VECT>
void interpolation_solution(getfem::mesh_fem &mf, const getfem::mesh_fem &mf_target,
			    const VECT &U, VECT &V, dim_type P) {
  if (mf.get_qdim() == 1) {
    mf.set_qdim(P);
    getfem::interpolation_solution(mf,mf_target,U,V); mf.set_qdim(1);
  }
  else getfem::interpolation_solution(mf,mf_target,U,V);
}


void err_msg(void)
{
  cerr << "Bad format for arguments of command compare_solutions\n";
  cerr << "The right format is\n";
  cerr << "compare_solutions filename1 filename2\n";
  cerr << "where filename1 and filename2 are two files containing some\n";
  cerr << "solutions exported by GETFEM++.\n";
  exit(1);
}

int main(int argc, char *argv[])
{
  try {

    int found = 0;
    std::string fi1, fi2;

    for (int aa = 1; aa < argc; aa++) {
      if (argv[aa][0] != '-') {
	switch(found) {
        case 0  : fi1 = std::string(argv[aa]); found++; break;
        case 1  : fi2 = std::string(argv[aa]); found++; break;
        default : err_msg();
	}
      }
      else
	err_msg();
    }
    if (found != 2) err_msg();

    cout.precision(16);

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
    dim_type N2, P2; short_type K2;
    std::vector<getfem::scalar_type> U2, U3;
    getfem::load_solution(fi2, mesh2, mef2, U2, K2);
    N2 = mesh1.dim();
    P2 = mef2.get_qdim();
    
    if (N1 != N2) DAL_THROW(std::invalid_argument,
			    "Dimensions of the two meshes mismatch\n");
    
    if (P1 != P2) DAL_THROW(std::invalid_argument,
			    "Dimensions of the two solutions mismatch\n");
    
    // test d'égalité des maillages
    
    if (mesh1.nb_points() == mesh2.nb_points()
	&& mesh1.nb_convex() == mesh2.nb_convex()
	&& mef1.nb_dof() == mef2.nb_dof()) {
      // + test complet : pour chaque point verifier qu'il existe dans l'autre
      // maillage puis pour chaque convexe aussi ...
      // faire un tableau de correspondance des points.
      // ou une fonction mesh_egal ... ?
      
      if (false) {
	getfem::scalar_type errin = 0.0;
	std::vector<getfem::scalar_type>::iterator it=U1.begin(), ite=U1.end();
	std::vector<getfem::scalar_type>::iterator it2 = U2.begin();
	for ( ; it != ite; ++it, ++it2)
	  { *it -= *it2; errin = std::max(errin, dal::abs(*it)); }
	cout <<  "Computation of norms\n";
	getfem::scalar_type errl2 = getfem::asm_L2_norm(mef2, U1);
	getfem::scalar_type errh1 = getfem::asm_H1_semi_norm(mef2, U1);
	cout <<  "L^2     Error : " << errl2 << endl;
	cout <<  "H^1     Error : " << sqrt(errh1*errh1 + errl2*errl2) << endl;
	cout <<  "L^infty Error : " << errin << endl;
	return 0;
      }
    }
  
    cout << "Begin comparaison\n";
    if (mef1.nb_dof() < mef2.nb_dof()) {
      cout << "interpolation of the solution in " << fi1
	   << " on the mesh of " << fi2 << endl;
      U3.resize(mef2.nb_dof() * P1);
      getfem::scalar_type errin = 0.0;
      interpolation_solution(mef1, mef2, U1, U3, P1);
      std::vector<getfem::scalar_type>::iterator it=U3.begin(), ite = U3.end();
      std::vector<getfem::scalar_type>::iterator it2 = U2.begin();
      for ( ; it != ite; ++it, ++it2)
	{ *it -= *it2; errin = std::max(errin, dal::abs(*it)); }
      
      getfem::scalar_type errl2 = getfem::asm_L2_norm(mef2, U3);
      getfem::scalar_type errh1 = getfem::asm_H1_semi_norm(mef2, U3);
      cout <<  "L^2     Error : " << errl2 << endl;
      cout <<  "H^1     Error : " << sqrt(errh1*errh1 + errl2*errl2) << endl;
      cout <<  "L^infty Error : " << errin << endl;
    }
    else {
      cout << "interpolation of the solution in " << fi2
	   << " on the mesh of " << fi1 << endl;
      U3.resize(mef1.nb_dof()*P1);
      getfem::scalar_type errin = 0.0;
      interpolation_solution(mef2, mef1, U2, U3, P1);
      std::vector<getfem::scalar_type>::iterator it = U3.begin(), ite=U3.end();
      std::vector<getfem::scalar_type>::iterator it2 = U1.begin();
      for ( ; it != ite; ++it, ++it2)
	{ *it -= *it2; errin = std::max(errin, dal::abs(*it)); }
      
      getfem::scalar_type errl2 = getfem::asm_L2_norm(mef1, U3);
      getfem::scalar_type errh1 = getfem::asm_H1_semi_norm(mef1, U3);
      cout <<  "L^2     Error : " << errl2 << endl;
      cout <<  "H^1     Error : " << sqrt(errh1*errh1 + errl2*errl2) << endl;
      cout <<  "L^infty Error : " << errin << endl;
    }
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
  
}
