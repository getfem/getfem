/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */
/**************************************************************************/
/*                                                                        */
/*  Test the virtual_link_fem method.                                     */
/*                                                                        */
/**************************************************************************/

#include <getfem_assembling.h>
#include <getfem_export.h>
#include <getfem_norm.h>
#include <getfem_regular_meshes.h>
#include <gmm.h>

using bgeot::base_vector;
using bgeot::base_small_vector;
using bgeot::base_node;
using bgeot::scalar_type;
using bgeot::size_type;
using bgeot::dim_type;

typedef gmm::wsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef std::vector<scalar_type> linalg_vector;

/**************************************************************************/
/*  structure representing the problem.                                   */
/**************************************************************************/

struct lap_pb
{
  getfem::getfem_mesh mesh1, mesh2;
  getfem::mesh_fem     mef1,  mef2, meflink;

  scalar_type LX, LY, LZ;
  int NX1, NX2, N, K, KI, integration;

  ftool::md_param PARAM;

  void assemble(void);
  void init(void);
  lap_pb(void) : mef1(mesh1), mef2(mesh2), meflink(mesh1) {}
};

void lap_pb::init(void)
{
  dal::bit_vector nn;

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
  K = PARAM.int_value("K", "Finite element degree");
  KI = PARAM.int_value("KI", "Integration degree");
  
  /***********************************************************************/
  /*  BUILD MESH.                                                        */
  /***********************************************************************/

  cout << "Mesh generation\n";

  base_node org(N); org.fill(0.0);
  std::vector<base_small_vector> vtab(N);
  std::vector<size_type> ref(N);
  std::fill(ref.begin(), ref.end(), NX1);
  for (dim_type i = 0; i < N; i++)
  { 
    vtab[i] = base_small_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX1);
  }
  getfem::parallelepiped_regular_simplex_mesh(mesh1, N, org,
					      vtab.begin(), ref.begin());
  mesh1.optimize_structure();

  std::fill(ref.begin(), ref.end(), NX2);
  for (dim_type i = 0; i < N; i++)
  { 
    vtab[i] = base_small_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX2);
  }
  getfem::parallelepiped_regular_simplex_mesh(mesh2, N, org,
					      vtab.begin(), ref.begin());
  mesh2.optimize_structure();
  cout << "Selecting finite element method.\n";
  char meth[500];
  getfem::pintegration_method ppi;
  switch (integration) {
  case 0  : sprintf(meth, "IM_EXACT_SIMPLEX(%d)", int(N)); break;
  case 1  : sprintf(meth, "IM_NC(%d, %d)", int(N), int(KI)); break;
  case 2  : sprintf(meth, "IM_GAUSS1D(%d)", int(KI)); break;
  case 3  : sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_NC(%d, %d), %d)",
		    int(N), int(2*K), int(KI)); break;
  case 11 : sprintf(meth, "IM_TRIANGLE(1)"); break;
  case 12 : sprintf(meth, "IM_TRIANGLE(2)"); break;
  case 13 : sprintf(meth, "IM_TRIANGLE(3)"); break;
  case 14 : sprintf(meth, "IM_TRIANGLE(4)"); break;
  case 15 : sprintf(meth, "IM_TRIANGLE(5)"); break;
  case 16 : sprintf(meth, "IM_TRIANGLE(6)"); break;
  case 17 : sprintf(meth, "IM_TRIANGLE(7)"); break;
  case 21 : sprintf(meth, "IM_TETRAHEDRON(1)"); break;
  case 22 : sprintf(meth, "IM_TETRAHEDRON(2)"); break;
  case 23 : sprintf(meth, "IM_TETRAHEDRON(3)"); break;
  case 25 : sprintf(meth, "IM_TETRAHEDRON(5)"); break;
  default : DAL_THROW(std::logic_error, "Undefined integration method");
  }
  ppi = getfem::int_method_descriptor(meth);
  
  sprintf(meth, "FEM_PK(%d,%d)", int(N), int(K));
  nn = mesh1.convex_index(N);
  mef1.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
  nn = mesh2.convex_index(N);
  mef2.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
  nn = mesh1.convex_index(N);
  meflink.set_finite_element(nn, getfem::virtual_link_fem(mef2, meflink, ppi),
			     ppi);
 
}

void lap_pb::assemble(void)
{
  int nb_dof1 = mef1.nb_dof(), nb_dof2 = mef2.nb_dof();
  sparse_matrix_type RM1(nb_dof2, nb_dof2);
  double sum, diff;
  
  cout << "Number of dof : " << nb_dof1 << " : " << nb_dof2 << endl;

  cout << "Number of dof of interpolated method: " << meflink.nb_dof() << endl;
 
  cout << "Assembling interpolated mass matrix" << endl;
  getfem::asm_mass_matrix(RM1, meflink, meflink);

  cout << "Matrice de masse\n";
  sum = 0.0;
  for (size_type i = 0; i < RM1.nrows(); i++) { 
    cout << "ligne " << i << " [ ";
    scalar_type slig = 0;
    for (size_type l = 0; l < RM1.nrows(); l++)
      if (RM1(i, l) != 0.0) {
	cout << "(" << l << "," << RM1(i, l) << ")  ";
	slig += RM1(i, l);
      }
    sum += slig;
    cout << "] -> sum(line)=" << slig << endl;
  }
  cout << endl << " sum: " << sum << endl << endl;

  sparse_matrix_type RM2 = sparse_matrix_type(nb_dof2, nb_dof2);
  cout << "Assembling normal mass matrix" << endl;
  getfem::asm_mass_matrix(RM2, mef2);
  
  cout << "Matrice de masse\n";
  sum = 0.0; diff = 0.0;
  for (size_type i = 0; i < RM2.nrows(); i++) { 
    cout << "ligne " << i << " [ ";
    scalar_type slig = 0;
    for (size_type l = 0; l < RM2.nrows(); l++) {
      diff += dal::abs(RM2(i, l) - RM1(i, l));
      if (RM2(i, l) != 0.0) {
	cout << "(" << l << "," << RM2(i, l) << ")  ";
	slig = slig + RM2(i, l);
      }
    }
    sum += slig;
    cout << "] -> sum(line)=" << slig << endl;
  }
  cout << endl << " sum: " << sum << endl << endl;
  cout << endl << " diff: " << diff << "max_norm=" << gmm::mat_maxnorm(RM1) << endl << endl;

  
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
