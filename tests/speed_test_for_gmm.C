/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2003  Yves Renard.                                   */
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
/*  Speed tests for GMM, matrix coming from an elastostatic computation.   */
/*                                                                         */
/* *********************************************************************** */

#include <getfem_assembling.h>
#include <getfem_norm.h>
#include <getfem_regular_meshes.h>
#include <gmm.h>

using bgeot::base_small_vector;
using bgeot::base_vector;
using bgeot::base_node;
using bgeot::base_matrix;
using bgeot::scalar_type;
using bgeot::size_type;
using bgeot::dim_type;

typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef std::vector<scalar_type> linalg_vector;

/**************************************************************************/
/*  Definition de la solution test.                                       */
/**************************************************************************/

dal::dynamic_array<base_small_vector> sol_K;
scalar_type sol_lambda, sol_G;

base_small_vector sol_u(const base_node &x) {
  int N = x.size(); base_small_vector res(N);
  for (int i = 0; i < N; ++i) res[i] = sin(bgeot::vect_sp(sol_K[i], x));
  return res;
}

base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N);
  for (int i = 0; i < N; i++) {
    res[i] = ( sol_G * bgeot::vect_sp(sol_K[i], sol_K[i]) )
                  * sin(bgeot::vect_sp(sol_K[i], x));
    for (int j = 0; j < N; j++)
      res[i] += ( (sol_lambda + sol_G) * sol_K[j][j] * sol_K[j][i])
	          * sin(bgeot::vect_sp(sol_K[j], x));
  }
  return res;
}

base_matrix sol_sigma(const base_node &x) {
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
  size_type NX, N, K, KI           ;

  sparse_matrix_type SM; /* matrice de rigidite.                     */
  linalg_vector U, B; /* inconnue et second membre.                       */

  ftool::md_param PBSTFR_PARAM;

  std::string datafilename;

  void init(size_type, size_type);
  void assemble(void);
  void solve(void);

  pb_data(void) : mef(mesh), mef_data(mesh), mef_data2(mesh) {}
};

void pb_data::init(size_type NNX, size_type KK) {
  dal::bit_vector nn;
  size_type i, j, k;

  /***********************************************************************/
  /*  LECTURE DES PARAMETRES SUR FICHIER.                                */
  /***********************************************************************/
  
  /* parametres physiques */
  N = 2;
  G = 1.0;
  mef.set_qdim(N);
  lambda = 1.0;
  
  /* parametres numeriques */
  LX = LY = LZ = 1.0;
  KI = 1;
  NX = NNX;
  residu = 10E-9;
  K = KK;
  for (i = 0; i < N; i++) {
    sol_K[i] = base_small_vector(N);
    for (j = 0; j < N; j++)  sol_K[i][j] = (i == j) ? 0.1 : -0.1;
  }
  sol_lambda = lambda; sol_G = G;

  /***********************************************************************/
  /*  CONSTRUCTION DU MAILLAGE.                                          */
  /***********************************************************************/

  // cout << "Mesh generation\n";
  base_node org(N); org.fill(0.0);
  std::vector<base_small_vector> vtab(N);
  std::vector<size_type> ref(N); std::fill(ref.begin(), ref.end(), NX);
  for (i = 0; i < N; i++) { 
    vtab[i] = base_small_vector(N); vtab[i].fill(0.0);
    (vtab[i])[i] = ((i == 0) ? LX : ((i == 1) ? LY : LZ)) / scalar_type(NX);
  }
  
  getfem::parallelepiped_regular_simplex_mesh(mesh, N, 
					      org, vtab.begin(), ref.begin());
  mesh.optimize_structure();

  // cout << "Selecting finite element method.\n";
  char meth[500];
  getfem::pintegration_method ppi;
  nn = mesh.convex_index(N);
  sprintf(meth, "IM_EXACT_SIMPLEX(%d)", int(N)); 
  ppi = getfem::int_method_descriptor(meth);
  getfem::pfem pfprinc = 0;
  sprintf(meth, "FEM_PK(%d,%d)", int(N), int(K));
  pfprinc = getfem::fem_descriptor(meth);
  mef.set_finite_element(nn, pfprinc, ppi);
  mef_data.set_finite_element(nn, pfprinc,
			      getfem::exact_simplex_im(N));
  sprintf(meth, "FEM_PK(%d,%d)", int(N), 0);
  mef_data2.set_finite_element(nn, pfprinc,
			       getfem::exact_simplex_im(N));
  
  // cout << "Selecting Neumann and Dirichlet boundaries\n";
  nn = mesh.convex_index(N);
  base_small_vector un;
  for (j << nn; j != size_type(-1); j << nn) {
    k = mesh.structure_of_convex(j)->nb_faces();
    for (i = 0; i < k; i++) {
      if (bgeot::neighbour_of_convex(mesh, j, i).empty()) {
	un = mesh.normal_of_face_of_convex(j, i, 0);
	un /= bgeot::vect_norm2(un);
	
	// if (true)
	if (dal::abs(un[N-1] - 1.0) < 1.0E-3)
	  mef.add_boundary_elt(0, j, i);
	else {
	  size_type nb;
	  for(nb = 0; nb < N; ++nb) {
	    if (dal::abs(un[nb] + 1.0) < 1.0E-7) { nb = 1 + nb * 2; break; }
	    if (dal::abs(un[nb] - 1.0) < 1.0E-7) { nb = 2 + nb * 2; break; }
	  }
	  mef.add_boundary_elt(nb, j, i);
	  // cout << "ajout bord Neumann, cv\t" << j << "\tf " << i << endl;
	}
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

  cout << "dof number for linear elasticity : " << nb_dof << endl;
  
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


  ST1 = linalg_vector(nb_dof_data * N);
  getfem::base_node pt(N);
  for (size_type nb = 1; nb <= 2*N; ++nb) {
    dal::bit_vector nn = mesh.convex_index(N);
    getfem::base_small_vector un, v(N);
    size_type j;
     for (j << nn; j != size_type(-1); j << nn) {
      getfem::pfem pf = mef_data.fem_of_element(j);
      dal::bit_vector nf = mef.faces_of_convex_on_boundary(j, nb);
      size_type i;
      for (i << nf; i != size_type(-1); i << nf) {
	for (size_type l = 0; l<pf->structure()->nb_points_of_face(i); ++l) {
	  size_type n = pf->structure()->ind_points_of_face(i)[l];
	  un = mesh.normal_of_face_of_convex(j, i, pf->node_of_dof(n));
	  un /= bgeot::vect_norm2(un);
	  size_type dof = mef_data.ind_dof_of_element(j)[n];
	  gmm::mult(sol_sigma(mef_data.point_of_dof(dof)), un, v);
	  for (size_type k = 0; k < N; ++k) ST1[dof*N+k] = v[k];
	}
      }
    }
    getfem::asm_source_term(B, mef, mef_data, ST1, nb);
  }

  ST1 = linalg_vector(nb_dof);
  for (size_type i = 0; i < nb_dof/N; ++i)
    for (size_type j = 0; j < N; ++j)
      ST1[i*N+j] = sol_u(mef.point_of_dof(i*N))[j];
  getfem::assembling_Dirichlet_condition(SM, B, mef, 0, ST1);
}

void pb_data::solve(void) {
  gmm::iteration iter(residu, 1);
  gmm::choleskyt_precond<sparse_matrix_type> P(SM, 10, 1E-7);
  // gmm::cholesky_precond<sparse_matrix_type> P(SM);
  gmm::cg(SM, U, B, P, iter);
}

int main(void) {
  try {

    cout << "***********************************************************\n";
    cout << "*      Dense matrix-vector mult.                          *\n";  
    cout << "***********************************************************\n";
    {
      cout << "\nbuilding the stiffness matrix\n";
      scalar_type exectime = ftool::uclock_sec();
      pb_data p;
      p.init(5 /* NX */ , 3 /* K */);
      p.assemble();
      cout << "Stiffness matrix computation : "
	   << ftool::uclock_sec() - exectime << " seconds\n";

      gmm::dense_matrix<double> M(gmm::mat_nrows(p.SM), gmm::mat_nrows(p.SM));
      gmm::copy(p.SM, M);

      exectime = ftool::uclock_sec();
      for (size_type i = 0; i < 10; ++i) {
	gmm::mult(M, p.B, p.U);
      }
      cout << "10 x Dense matrix-vector mult : "
	   << ftool::uclock_sec() - exectime << " seconds\n";
    
      cout << "***********************************************************\n";
      cout << "*      Dense LU decomposition.                            *\n";
      cout << "***********************************************************\n";
      
      exectime = ftool::uclock_sec();
      for (size_type i = 0; i < 1; ++i)
	gmm::lu_solve(M, p.U, p.B);
      cout << "\nLu solve : "
	   << ftool::uclock_sec() - exectime << " seconds\n";
      
    }
    /*

    cout << "Solving linear system\n";
    exectime = ftool::uclock_sec();
    p.solve();
    cout << "Time for solving : "
	 << ftool::uclock_sec() - exectime << " seconds\n";

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

    */
  }
  DAL_STANDARD_CATCH_ERROR;

  return 0;
}
