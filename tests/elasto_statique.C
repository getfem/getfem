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
/* *********************************************************************** */
/*                                                                         */
/*  Calcul sur une structure lineairement elastique                        */
/*                                                                         */
/* *********************************************************************** */

#include <getfem_assembling.h>
#include <getfem_norm.h>
#include <getfem_regular_meshes.h>
#include <gmm.h>

using bgeot::base_vector;
using bgeot::base_node;
using bgeot::base_matrix;
using bgeot::scalar_type;
using bgeot::size_type;
using bgeot::dim_type;

typedef gmm::wsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef std::vector<scalar_type> linalg_vector;



/**************************************************************************/
/*  Definition de la solution test.                                       */
/**************************************************************************/

dal::dynamic_array<base_vector> sol_K;
scalar_type sol_lambda, sol_G;

base_vector sol_u(const base_vector &x) {
  int N = x.size(); base_vector res(N);
  for (int i = 0; i < N; ++i) res[i] = sin(bgeot::vect_sp(sol_K[i], x));
  return res;
}

base_vector sol_f(const base_vector &x) {
  int N = x.size();
  base_vector res(N);
  for (int i = 0; i < N; i++) {
    res[i] = ( sol_G * bgeot::vect_sp(sol_K[i], sol_K[i]) )
                  * sin(bgeot::vect_sp(sol_K[i], x));
    for (int j = 0; j < N; j++)
      res[i] += ( (sol_lambda + sol_G) * sol_K[j][j] * sol_K[j][i])
	          * sin(bgeot::vect_sp(sol_K[j], x));
  }
  return res;
}

base_matrix sol_sigma(const base_vector &x) {
  int N = x.size();
  base_matrix res(N,N);
  for (int i = 0; i < N; i++)
    for (int j = 0; j <= i; j++) {
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

struct pb_data {
  getfem::getfem_mesh mesh;
  getfem::mesh_fem mef;
  getfem::mesh_fem mef_data;
  getfem::mesh_fem mef_data2;

  scalar_type G, lambda;
  scalar_type LX, LY, LZ, residu;
  size_type NX, N, K, KI, mesh_type;

  sparse_matrix_type SM; /* matrice de rigidite.                     */
  linalg_vector U, B; /* inconnue et second membre.                       */

  int integration;

  ftool::md_param PBSTFR_PARAM;

  std::string datafilename;

  void init();
  void assemble(void);
  void solve(void);

  pb_data(void) : mef(mesh), mef_data(mesh), mef_data2(mesh) {}
};

void pb_data::init(void) {
  dal::bit_vector nn;
  size_type i, j, k;

  /***********************************************************************/
  /*  LECTURE DES PARAMETRES SUR FICHIER.                                */
  /***********************************************************************/
  
  /* parametres physiques */
  N = PBSTFR_PARAM.int_value("N", "Domain dimension");
  G = PBSTFR_PARAM.real_value("G", "Lamé coefficient G");
  mef.set_qdim(N);
  lambda = PBSTFR_PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  
  /* parametres numeriques */
  LX = PBSTFR_PARAM.real_value("LX", "Size in X");
  LY = PBSTFR_PARAM.real_value("LY", "Size in Y");
  LZ = PBSTFR_PARAM.real_value("LZ", "Size in Y");
  mesh_type = PBSTFR_PARAM.int_value("MESH_TYPE", "Mesh type ");
  KI = PBSTFR_PARAM.int_value("KI", "Integration degree");
  NX = PBSTFR_PARAM.int_value("NX", "Number of space steps ");
  integration = PBSTFR_PARAM.int_value("INTEGRATION", "integration method");
  residu = PBSTFR_PARAM.real_value("RESIDU", "residu for c.g.");
  K = PBSTFR_PARAM.int_value("K", "Finite element degree");

  datafilename = std::string(PBSTFR_PARAM.string_value("ROOTFILENAME",
						       "Filename for saving"));
  scalar_type FT = PBSTFR_PARAM.real_value("FT", 
					   "parameter for exact solution");
  for (i = 0; i < N; i++) {
    sol_K[i] = base_vector(N);
    for (j = 0; j < N; j++)  sol_K[i][j] = (i == j) ? FT : -FT;
  }

  sol_lambda = lambda; sol_G = G;

  /***********************************************************************/
  /*  CONSTRUCTION DU MAILLAGE.                                          */
  /***********************************************************************/

  cout << "Mesh generation\n";

  base_node org(N); org.fill(0.0);
  std::vector<base_vector> vtab(N);
  std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
  for (i = 0; i < N; i++) { 
    vtab[i] = base_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
  }

  if (mesh_type == 2 && N <= 1) mesh_type = 0;
  
  switch (mesh_type) {
  case 0 : getfem::parallelepiped_regular_simplex_mesh
		   (mesh, N, org, vtab.begin(), ref.begin()); break;
  case 1 : getfem::parallelepiped_regular_mesh
		   (mesh, N, org, vtab.begin(), ref.begin()); break;
  case 2 : getfem::parallelepiped_regular_prism_mesh
		   (mesh, N, org, vtab.begin(), ref.begin()); break;
  default : DAL_THROW(dal::internal_error, "Unknown type of mesh");
  }

  mesh.optimize_structure();

  cout << "Selecting finite element method.\n";
  char meth[500];
  getfem::pintegration_method ppi;
  nn = mesh.convex_index(N);

  switch (integration) {
  case 0 :
    switch (mesh_type) { 
    case 0 : sprintf(meth, "IM_EXACT_SIMPLEX(%d)", int(N)); break;
    case 1 : sprintf(meth, "IM_EXACT_PARALLELEPIPED(%d)", int(N)); break;
    default : DAL_THROW(dal::internal_error, 
    "Exact integration not allowed in this context");
    }
    break;
  case 1 :
    switch (mesh_type) { 
    case 0 : 
      sprintf(meth, "IM_NC(%d,%d)", int(N), int(2*K));
      break;
    case 1 : 
      sprintf(meth, "IM_NC_PARALLELEPIPED(%d,%d)", int(N), int(2*K));
      break;
    case 2 :
      sprintf(meth, "IM_NC_PRISM(%d,%d)", int(N), int(2*K));
      break;
    }
    break;
  case 2 :
    if (mesh_type == 1)
      sprintf(meth, "IM_GAUSS_PARALLELEPIPED(%d,%d)", int(N), int(2*K));
    else
      DAL_THROW(dal::internal_error,
		"Product of 1D Gauss only for parallelepipeds");
    break;
  case 3 :
    if (mesh_type == 0) {
      if (N == 1)
	sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_GAUSS1D(%d), %d)",2,int(KI));
      else if (N == 2)
	sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(%d), %d)",
		2,int(KI));
      else
	sprintf(meth, "IM_STRUCTURED_COMPOSITE(IM_NC(%d, %d), %d)",
		int(N), int(2*K), int(KI));
    }
    else
      DAL_THROW(dal::internal_error,
		"Composite integration only for simplexes");
    break;
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
  case 32 : sprintf(meth, "IM_QUAD(2)"); break;
  case 33 : sprintf(meth, "IM_QUAD(3)"); break;
  case 35 : sprintf(meth, "IM_QUAD(5)"); break;
  default : DAL_THROW(std::logic_error, "Undefined integration method");
  }
  ppi = getfem::int_method_descriptor(meth);
  getfem::pfem pfprinc = 0;
  switch (mesh_type) {
  case 0 :
    sprintf(meth, "FEM_PK(%d,%d)", int(N), int(K));
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, pfprinc, ppi);
    mef_data.set_finite_element(nn, pfprinc,
				getfem::exact_simplex_im(N));
    sprintf(meth, "FEM_PK(%d,%d)", int(N), 0);
    mef_data2.set_finite_element(nn, pfprinc,
				 getfem::exact_simplex_im(N));
    break;
  case 1 :
    sprintf(meth, "FEM_QK(%d,%d)", int(N), K);
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, pfprinc, ppi); 
    mef_data.set_finite_element(nn, pfprinc, ppi);
    sprintf(meth, "FEM_QK(%d,%d)", int(N), 0);
    mef_data2.set_finite_element(nn, pfprinc,  ppi);
    break;
  case 2 :
    sprintf(meth, "FEM_PK_PRISM(%d,%d)", int(N), K);
    pfprinc = getfem::fem_descriptor(meth);
    mef.set_finite_element(nn, pfprinc, ppi);
    mef_data.set_finite_element(nn, pfprinc, ppi);
    sprintf(meth, "FEM_PK_PRISM(%d,%d)", int(N), 0);
    mef_data2.set_finite_element(nn, pfprinc, ppi);
    break;
  }

//   if (integration == 0) sprintf(meth, "IM_EXACT_SIMPLEX(%d)", int(N));
//   else sprintf(meth, "IM_NC(%d, %d)", int(N), int(2*K));
//   ppi = getfem::int_method_descriptor(meth);
//   sprintf(meth, "FEM_PK(%d,%d)", int(N), int(K));
//   mef.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
//   mef_data.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);
//   sprintf(meth, "FEM_PK(%d,%d)", int(N), 0);
//   mef_data2.set_finite_element(nn, getfem::fem_descriptor(meth), ppi);

  cout << "Selecting Neumann and Dirichlet boundaries\n";
  nn = mesh.convex_index(N);
  base_vector un;
  for (j << nn; j != size_type(-1); j << nn) {
    k = mesh.structure_of_convex(j)->nb_faces();
    for (i = 0; i < k; i++) {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty()) {
	un = mesh.normal_of_face_of_convex(j, i, 0);
	un /= bgeot::vect_norm2(un);
	
	// if (true)
	if (dal::abs(un[N-1] - 1.0) < 1.0E-3)
	  mef.add_boundary_elt(0, j, i);
	else
	  mef.add_boundary_elt(1, j, i);
      }
    }
  }
}

void pb_data::assemble(void)
{
  size_type nb_dof = mef.nb_dof(), nb_dof_data = mef_data.nb_dof();
  size_type nb_dof_data2 = mef_data2.nb_dof();
  B = linalg_vector(nb_dof); gmm::clear(B);
  U = linalg_vector(nb_dof); gmm::clear(U);
  SM = sparse_matrix_type(nb_dof, nb_dof);
  linalg_vector ST1, ST2;

  cout << "dof number fo linear elasticity : " << nb_dof << ", data nb_dof=" << mef_data.nb_dof() << endl;
  
  // cout << "Assemblage de la matrice de rigidite" << endl;
  ST1 = linalg_vector(nb_dof_data2); ST2 = linalg_vector(nb_dof_data2);
  std::fill(ST1.begin(), ST1.end(), lambda);
  std::fill(ST2.begin(), ST2.end(), G);
  getfem::asm_stiffness_matrix_for_linear_elasticity(SM,mef,mef_data2,ST1,ST2);
  // cout << "Assemblage du terme source" << endl;
  ST1 = linalg_vector(nb_dof_data * N);
  for (size_type i = 0; i < nb_dof_data; ++i)
    for (size_type j = 0; j < N; ++j) 
      ST1[i*N+j] = sol_f(mef_data.point_of_dof(i))[j];
  //  getfem::assembling_volumic_source_term(B, mef, mef_data, ST1, N);
  getfem::asm_source_term(B, mef, mef_data, ST1);

  // cout << "Assemblage de la condition de Neumann" << endl;
  ST1 = linalg_vector(nb_dof_data * N);
  getfem::base_node pt(N); getfem::base_vector n(N), v;
  for (size_type i = 0; i < nb_dof_data; ++i) {
    pt = mef_data.point_of_dof(i);
    if (dal::abs(pt[0]-LX) < 10E-6) n[0] = 1.0; // pas terrible ... !!
      else if (dal::abs(pt[0]) < 10E-6) n[0] = -1.0; else n[0] = 0.0;
    if (N > 1) {
      if (dal::abs(pt[1]-LY) < 10E-6) n[1] = 1.0;
      else if (dal::abs(pt[1]) < 10E-6) n[1] = -1.0; else n[1] = 0.0;
      for (dim_type k = 2; k < N; ++k)
	if (dal::abs(pt[k]-LZ) < 10E-6) n[k] = 1.0;
	else if (dal::abs(pt[k]) < 10E-6) n[k] = -1.0; else n[k] = 0.0;
    }
    v = sol_sigma(pt) *  n;
    for (size_type j = 0; j < N; ++j)
      ST1[i*N+j] = v[j];
  }
  getfem::asm_source_term(B, mef, mef_data, ST1, 1);

  // cout << "Prise en compte de la condition de Dirichlet" << endl;
  ST1 = linalg_vector(nb_dof);
  for (size_type i = 0; i < nb_dof/N; ++i)
    for (size_type j = 0; j < N; ++j)
      ST1[i*N+j] = sol_u(mef.point_of_dof(i*N))[j];
  getfem::assembling_Dirichlet_condition(SM, B, mef, 0, ST1);
}

void pb_data::solve(void) {
  gmm::iteration iter(residu);
  gmm::cg(SM, U, B, gmm::identity_matrix(), gmm::identity_matrix(), iter);
}

int main(int argc, char *argv[]) {
  try {
    pb_data p;
    
    // cout << "initialisation ...\n";
    p.PBSTFR_PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.stat();
    p.mesh.write_to_file(p.datafilename + ".mesh" + char(0));
    
    cout << "Assembling\n";
    p.assemble();
    
    // cout << "Matrice de rigidite\n";
    // cout << p.SM << endl;
    
    cout << "Solving linear system\n";
    p.solve();
    
    cout << "Error computation\n";
    
    int nbdof = p.mef.nb_dof();
    linalg_vector V(nbdof); V = p.U;
    base_vector S;
    
    for (int i = 0; i < int(nbdof/p.N); ++i) {
      S = sol_u(p.mef.point_of_dof(i*p.N));
      for (dim_type k = 0; k < p.N; ++k) V[i*p.N + k] -= S[k];
    }
    
    cout <<  "L2 Error = " << getfem::asm_L2_norm(p.mef, V) << endl;
    cout <<  "H1 Error = " << getfem::asm_H1_norm(p.mef, V) << endl;
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
