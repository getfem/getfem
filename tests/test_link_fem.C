/**************************************************************************/
/*                                                                        */
/*  Test the virtual_link_fem method.                                     */
/*                                                                        */
/**************************************************************************/

#include <getfem_assembling.h>
#include <getfem_export.h>
#include <getfem_norm.h>
#include <getfem_regular_meshes.h>
#include <bgeot_smatrix.h>

using bgeot::base_vector;
using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::size_type;

typedef bgeot::smatrix<scalar_type> sparse_matrix_type;
typedef bgeot::vsvector<scalar_type> linalg_vector;

/**************************************************************************/
/*  structure representing the problem.                                   */
/**************************************************************************/

struct lap_pb
{
  getfem::getfem_mesh mesh1, mesh2;
  getfem::mesh_fem     mef1,  mef2, meflink;

  scalar_type LX, LY, LZ, residu;
  int NX1, NX2, N, K;

  int integration;

  ftool::md_param PARAM;

  void assemble(void);
  void init(void);
  lap_pb(void) : mef1(mesh1), mef2(mesh2), meflink(mesh1) {}
};

void lap_pb::init(void)
{
  dal::bit_vector nn;
  size_type i, j, k;

  /***********************************************************************/
  /* READING PARAMETER FILE                                              */
  /***********************************************************************/
  
  N = PARAM.int_value("N", "Domaine dimension");
  LX = PARAM.real_value("LX", "Size in X");
  LY = PARAM.real_value("LY", "Size in Y");
  LZ = PARAM.real_value("LZ", "Size in Y");
  NX1 = PARAM.int_value("NX1", "Nomber of sace steps ");
  NX2 = PARAM.int_value("NX2", "Nomber of sace steps ");
  integration = PARAM.int_value("INTEGRATION", "integration method");
  residu = PARAM.real_value("RESIDU", "Residu for c.g.");
  K = PARAM.int_value("K", "Finite element degree");
  
  /***********************************************************************/
  /*  BUILD MESH.                                                        */
  /***********************************************************************/

  cout << "Mesh generation\n";

  base_node org(N); org.fill(0.0);
  std::vector<base_vector> vtab(N);
  std::vector<size_type> ref(N);
  std::fill(ref.begin(), ref.end(), NX1);
  for (i = 0; i < N; i++)
  { 
    vtab[i] = base_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX1);
  }
  getfem::parallelepiped_regular_simplex_mesh(mesh1, N, org,
					      vtab.begin(), ref.begin());
  mesh1.optimize_structure();

  std::fill(ref.begin(), ref.end(), NX2);
  for (i = 0; i < N; i++)
  { 
    vtab[i] = base_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX2);
  }
  getfem::parallelepiped_regular_simplex_mesh(mesh2, N, org,
					      vtab.begin(), ref.begin());
  mesh2.optimize_structure();
  cout << "Selecting finite element method.\n";

  getfem::pintegration_method ppi;
  switch (integration) {
  case 0 : ppi = bgeot::simplex_poly_integration(N); break;
  case 1 : ppi = bgeot::Newton_Cotes_approx_integration(N,2*K); break;
  case 11 : ppi = bgeot::triangle1_approx_integration(); break;
  case 12 : ppi = bgeot::triangle2_approx_integration(); break;
  case 13 : ppi = bgeot::triangle3_approx_integration(); break;
  case 14 : ppi = bgeot::triangle4_approx_integration(); break;
  case 15 : ppi = bgeot::triangle5_approx_integration(); break;
  case 16 : ppi = bgeot::triangle6_approx_integration(); break;
  case 17 : ppi = bgeot::triangle7_approx_integration(); break;
  case 21 : ppi = bgeot::tetrahedron1_approx_integration(); break;
  case 22 : ppi = bgeot::tetrahedron2_approx_integration(); break;
  case 23 : ppi = bgeot::tetrahedron3_approx_integration(); break;
  case 25 : ppi = bgeot::tetrahedron5_approx_integration(); break;
  default : DAL_THROW(std::logic_error, "Undefined integration method");
  }
  nn = mesh1.convex_index(N);
  mef1.set_finite_element(nn, getfem::PK_fem(N, K), ppi);
  nn = mesh2.convex_index(N);
  mef2.set_finite_element(nn, getfem::PK_fem(N, K), ppi);
  nn = mesh1.convex_index(N);
  meflink.set_finite_element(nn, getfem::virtual_link_fem(mef2, meflink, ppi),
			     ppi);
 
}

void lap_pb::assemble(void)
{
  int nb_dof1 = mef1.nb_dof(), nb_dof2 = mef2.nb_dof();
  sparse_matrix_type RM(nb_dof1, nb_dof2);
  
  cout << "Number of dof : " << nb_dof1 << " : " << nb_dof2 << endl;

  cout << "Number of dof of interpolated method: " << meflink.nb_dof() << endl;
 
  cout << "Assembling mass matrix" << endl;
  getfem::mass_matrix(RM, mef1, meflink, 1);

  cout << "Matrice de rigidite\n";
  for (int i = 0; i < RM.nrows(); i++) { 
    cout << "ligne " << i << " [ ";
    for (int l = 0; l < RM.nrows(); l++)
      if (RM(i, l) != 0.0)
	cout << "(" << l << "," << RM(i, l) << ")  ";
    cout << "]" << endl;
  }
  cout << endl << endl;
  
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[])
{
  try {
    
    lap_pb p;
    
    cout << "initialisation ...\n";
    p.PARAM.read_command_line(argc, argv);
    p.init();
    
    cout << "Assembling \n";
    p.assemble();
    
  }
  DAL_STANDARD_CATCH_ERROR;
  return 0; 
}
