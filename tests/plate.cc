// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2002-2008 Yves Renard, Michel Salaün.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

/**
   @file plate.cc
   @brief Linear Elastostatic plate problem.
   
   This program is used to check that getfem++ is working. This is
   also a good example of use of Getfem++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_linearized_plates.h"
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/*
  structure for the elastostatic problem
*/
struct plate_problem {

  enum { SIMPLY_FIXED_BOUNDARY_NUM = 0 };
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim, mim_subint;
  getfem::mesh_fem mf_ut;
  getfem::mesh_fem mf_u3;
  getfem::mesh_fem mf_theta;
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type epsilon;       /* thickness of the plate.                      */
  scalar_type pressure;
  scalar_type residual;        /* max residual for the iterative solvers         */
  scalar_type LX , LY ;       // default : LX = LY = 1
  bool mixed, symmetrized;
  bool mitc;
  int sol_ref;               // sol_ref = 0 : simple support on the vertical edges
                             // sol_ref = 1 : homogeneous on the vertical edges
                             // sol_ref = 2 : homogeneous on the 4 vertical
                             //       edges with solution u3 = sin²(x)*sin²(y)
  scalar_type eta;           // usefull only if sol_ref == 2 :
                             // eta = 0 => Kirchoff-Love
			     // eta = small => Mindlin 
  size_type N_Four ;
  base_matrix theta1_Four, theta2_Four, u3_Four ;
  
  int study_flag;            // if studyflag = 1, then the loadings applied are chosen
                             // in order to have a maximal vertical displacement equal to one.
			     // Nothing is done if study_flag has another value.   

  std::string datafilename;
  bgeot::md_param PARAM;

  base_small_vector theta_exact(base_node P);
  scalar_type u3_exact(base_node P);

  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);
  plate_problem(void) : mim(mesh), mim_subint(mesh), mf_ut(mesh), mf_u3(mesh),
			mf_theta(mesh), mf_rhs(mesh), mf_coef(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void plate_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE_UT  = PARAM.string_value("FEM_TYPE_UT","FEM name");
  std::string FEM_TYPE_U3  = PARAM.string_value("FEM_TYPE_U3","FEM name");
  std::string FEM_TYPE_THETA = PARAM.string_value("FEM_TYPE_THETA","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string INTEGRATION_CT = PARAM.string_value("INTEGRATION_CT",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE_UT="  << FEM_TYPE_UT << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  cout << "INTEGRATION_CT=" << INTEGRATION_CT << "\n";

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  GMM_ASSERT1(N == 2, "For a plate problem, N should be 2");
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Number of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  mixed = (PARAM.int_value("MIXED", "Mixed version ?") != 0);
  mitc = (PARAM.int_value("MITC", "Mitc version ?") != 0);
  symmetrized = (PARAM.int_value("SYMMETRIZED",
				 "Mixed symmetrized version ?") != 0);

  cout << "MITC = " ;
  if (mitc) cout << "true \n" ; else cout << "false \n" ;
  sol_ref = (PARAM.int_value("SOL_REF") ) ;
  study_flag = (PARAM.int_value("STUDY_FLAG") ) ;
  eta = (PARAM.real_value("ETA") );
  N_Four = (PARAM.int_value("N_Four") ) ;
  
    
  LX = PARAM.real_value("LX");
  LY = PARAM.real_value("LY");
  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  epsilon = PARAM.real_value("EPSILON", "thickness of the plate");
  pressure = PARAM.real_value("PRESSURE",
			      "pressure on the top surface of the plate.");
			      
			      
  cout << "SOL_REF = " ;
  if (sol_ref==0) cout << "appui simple aux 2 bords verticaux\n" ;
  if (sol_ref==1) cout << "encastrement aux 2 bords verticaux\n" ;
  if (sol_ref==2) {
     cout << "encastrement aux 4 bords verticaux, solution en sin(x)^2*sin(y)^2\n" ;
     cout << "eta = " << eta <<"\n";
  }
  if (sol_ref==4) {
     cout << "bord en appuis simple\n" ; 
     cout << "nombre de terme pour calcul sol exacte : " << N_Four << " \n" ;
     // Calcul des coeeficients de Fourier de la solution exacte :
     // Cas où le chargement est seulement vertical (pas de moment appliqué)
     gmm::resize( theta1_Four, N_Four, N_Four) ;
     gmm::resize( theta2_Four, N_Four, N_Four) ;
     gmm::resize( u3_Four, N_Four, N_Four) ;
     base_matrix Jmn(3, 3) ; 
     base_small_vector Bmn(3), Xmn(3) ; 
     scalar_type /*det_Jmn, */E, nu, A, B, e2, Pmn ;
     E = 4.*mu*(mu+lambda) / (2. * mu + lambda);
     nu = lambda / (2. * mu + lambda);
     e2 = epsilon * epsilon ;
     for(size_type i = 0 ; i < N_Four ; i++) {
        for(size_type j = 0 ; j < N_Four ; j++) {
	   A = (j + 1) * M_PI / LX ; 
	   B = (i + 1) * M_PI / LY ; 
	   Jmn(0, 0) = 2. * A * A / (1. - nu) + B * B + 3. / e2   ;
	   Jmn(0, 1) = A * B * (1. +nu) / (1. - nu) ;
	   Jmn(0, 2) = A * 3. / e2 ;
	   Jmn(1, 0) = A * B * (1. +nu) / (1. - nu) ;
	   Jmn(1, 1) = 2. * B * B / (1. - nu) + A * A + 3. / e2  ;
	   Jmn(1, 2) = B * 3. / e2 ;
	   Jmn(2, 0) = - A ;
	   Jmn(2, 1) = - B ;
	   Jmn(2, 2) = A * A + B * B ;
	   gmm::scale(Jmn,  - E*epsilon / (1. + nu) ) ;
	   
	   // calcul du développement de Fourrier du chargement :
	   if ( ( (i + 1) % 2 == 1 ) && ( (j + 1) % 2 == 1) ) {
	      Pmn =  16. * pressure / ( (i + 1) * (j + 1) * M_PI * M_PI) ; }
	   else {
	      Pmn = 0. ; }	      
	   Bmn[0] = 0. ;
	   Bmn[1] = 0. ;
	   Bmn[2] = Pmn ;
	   gmm::lu_solve(Jmn, Xmn, Bmn) ;
// 	   det_Jmn = pow( E * epsilon / (nu + 1.), 3) * (
// 	           ( A * A * 2. / (1. - nu)  + B * B + 3. / e2 ) * 
// 		   ( B * B * 2. / (1. - nu)  + A * A + 3. / e2 ) * ( A * A + B * B) 
// 		 + ( 6. * (1. + nu) * A * A * B * B / ( e2 * (1. - nu) ) )
// 		 - ( A * A * 2. / (1. - nu)  + B * B + 3. / e2 ) * (3. * A * A ) / e2 
// 		 - ( B * B * 2. / (1. - nu)  + A * A + 3. / e2 ) * (3. * B * B ) / e2  
// 		 - A * A * B * B * ( A * A + B * B ) * (1. + nu) * (1. + nu) / ( (1. - nu) * (1. - nu) )
// 		   )  ;
// 	   Jmn(2, 0) =   ( E * epsilon / (nu + 1.) ) * ( E * epsilon / (nu + 1.) ) * (
// 	                 ( 3. * (1. + nu) * A * B * B / ( e2 * (1. - nu) ) )
// 		       - ( B * B * 2. / (1. - nu)  + A * A + 3. / e2 ) * 3. * A / e2 
// 		         ) / det_Jmn ;
// 	   Jmn(2, 1) = - ( E * epsilon / (nu + 1.) ) * ( E * epsilon / (nu + 1.) ) * (
// 	                 ( A * A * 2. / (1. - nu)  + B * B + 3. / e2 ) * 3. * A / e2 
// 		       - ( 3. * (1. + nu) * A * A * B / ( e2 * (1. - nu) ) )
// 		         ) / det_Jmn ;
// 	   Jmn(2, 2) =   ( E * epsilon / (nu + 1.) ) * ( E * epsilon / (nu + 1.) ) * (
//                          ( A * A * 2. / (1. - nu)  + B * B + 3. /  e2  ) *
// 		         ( B * B * 2. / (1. - nu)  + A * A + 3. /  e2  ) 
// 		       - ( A * B * (1. + nu) / (1. - nu)) * ( A * B * (1. + nu) / (1. - nu) ) 
// 		         ) / det_Jmn ;
/*      
	   theta1_Four(i, j) = Jmn(2, 0) * Pmn ;
	   theta1_Four(i, j) = Jmn(2, 1) * Pmn ;
               u3_Four(i, j) = Jmn(2, 2) * Pmn ;*/
	   theta1_Four(i, j) = Xmn[0] ;
	   theta1_Four(i, j) = Xmn[1] ;
               u3_Four(i, j) = Xmn[2] ;	       
	   }
        }
     }

  mf_ut.set_qdim(N);
  mf_theta.set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_ut = getfem::fem_descriptor(FEM_TYPE_UT);
  getfem::pfem pf_u3 = getfem::fem_descriptor(FEM_TYPE_U3);
  getfem::pfem pf_theta = getfem::fem_descriptor(FEM_TYPE_THETA);
  
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method ppi_ct = 
    getfem::int_method_descriptor(INTEGRATION_CT);

  mim.set_integration_method(mesh.convex_index(), ppi);
  mim_subint.set_integration_method(mesh.convex_index(), ppi_ct);
  mf_ut.set_finite_element(mesh.convex_index(), pf_ut);
  mf_u3.set_finite_element(mesh.convex_index(), pf_u3);
  mf_theta.set_finite_element(mesh.convex_index(), pf_theta);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_ut->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_ut);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));

  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    assert(i.is_face());
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    switch(sol_ref){
       case 0 :
            if (gmm::abs(un[1]) <= 1.0E-7)  // new Neumann face
              mesh.region(SIMPLY_FIXED_BOUNDARY_NUM).add(i.cv(), i.f());
            break ;
       case 1 :
            if (gmm::abs(un[1]) <= 1.0E-7)  // new Neumann face
              mesh.region(SIMPLY_FIXED_BOUNDARY_NUM).add(i.cv(), i.f());
            break ;
       case 2 :
            if ( (gmm::abs(un[0]) <= 1.0E-7) || (gmm::abs(un[1]) <= 1.0E-7) ) 
              mesh.region(SIMPLY_FIXED_BOUNDARY_NUM).add(i.cv(), i.f());
            break ;
       case 3 :
            if (un[0] <= (- 1. + 1.0E-7))  // new Neumann face
              mesh.region(SIMPLY_FIXED_BOUNDARY_NUM).add(i.cv(), i.f());
            break ;
       case 4 :
            if ( (gmm::abs(un[0]) <= 1.0E-7) || (gmm::abs(un[1]) <= 1.0E-7) ) 
              mesh.region(SIMPLY_FIXED_BOUNDARY_NUM).add(i.cv(), i.f());
            break ;
       default :
	    GMM_ASSERT1(false, "SOL_REF parameter is undefined");
	    break ;
       } 
  }
}

base_small_vector plate_problem::theta_exact(base_node P) {
  base_small_vector theta(2);
  if (sol_ref == 0) { // appui simple aux 2 bords
     theta[0] = (-pressure / (32. * mu * epsilon * epsilon * epsilon))
              * (4. * pow(P[0] - .5, 3.) - 3 * (P[0] - .5));
     theta[1] = 0.;
   }
  if (sol_ref == 1) { // encastrement aux 2 bords
     theta[0] = (-pressure / (16. * mu * epsilon * epsilon * epsilon))
              * P[0] * ( 2.*P[0]*P[0] - 3.* P[0] + 1. ) ;
     theta[1] = 0.;
     }
  if (sol_ref == 2) { // encastrement aux 4 bords et sols vraiment 2D, non polynomiale
     theta[0] = - (1. + eta) * M_PI * sin(2.*M_PI*P[0]) * sin(M_PI*P[1])* sin(M_PI*P[1]) ;
     theta[1] = - (1. + eta) * M_PI * sin(2.*M_PI*P[1]) * sin(M_PI*P[0])* sin(M_PI*P[0]) ;
     }
  if (sol_ref == 3) { // plaque cantilever
     theta[0] = (- 3. * pressure / (8. * mu * epsilon * epsilon * epsilon ))
              * P[0] * ( 0.25 * P[0] * P[0] - P[0] + 1. ) ;
     theta[1] = 0.;      
     }
  if (sol_ref == 4) { // bord entier en appui simple
     theta[0] = 0. ;
     theta[1] = 0. ;
     for(size_type i = 0 ; i < N_Four ; i ++) {
        for(size_type j = 0 ; j < N_Four ; j ++) {
	   theta[0] += theta1_Four(i, j) * cos( (j + 1) * M_PI * P[0] / LX ) * sin( (i + 1) * M_PI * P[1] / LY ) ;
	   theta[0] += theta2_Four(i, j) * sin( (j + 1) * M_PI * P[0] / LX ) * cos( (i + 1) * M_PI * P[1] / LY ) ;
           }
        }
     }
  return theta;
}

scalar_type plate_problem::u3_exact(base_node P) {
  switch(sol_ref) {
  case 0 : return (pressure / (32. * mu * epsilon * epsilon * epsilon))
       * P[0] * (P[0] - 1.)
       * (gmm::sqr(P[0] - .5) -1.25-(mixed ? 0 : 8.*epsilon*epsilon));
       break ;
  case 1 : return (pressure /(32.* mu * epsilon * epsilon * epsilon))
       * P[0] * (P[0] - 1.)
       * ( P[0] * P[0] - P[0] - 8. * epsilon *epsilon) ;
       break ;
  case 2 : return  gmm::sqr(sin(M_PI*P[0])) * gmm::sqr(sin(M_PI*P[1]));
       break ;
  case 3 : return (3. * pressure / (4. * mu * epsilon * epsilon * epsilon ))
       * P[0] * ( P[0] * P[0] * P[0] / 24. - P[0] * P[0] / 6. + P[0] / 4. 
       - epsilon * epsilon * P[0] / 3. + 2. * epsilon * epsilon / 3.) ;
      break ;
  case 4 : 
       scalar_type u3_local ;
       u3_local = 0. ;
       for(size_type i = 0 ; i < N_Four ; i ++) {
          for(size_type j = 0 ; j < N_Four ; j ++) 
	     u3_local += u3_Four(i, j) * sin( (j + 1) * M_PI * P[0] / LX ) * sin( (i + 1) * M_PI * P[1] / LY ) ;
	     }
       return (u3_local) ;
       break ; 
  default : GMM_ASSERT1(false, "indice de solution de référence incorrect");
  }
}


/* compute the error with respect to the exact solution */
void plate_problem::compute_error(plain_vector &U) {
  cout.precision(16);

  std::vector<scalar_type> V(mf_rhs.nb_dof()*2);
  
  size_type i1 = mf_ut.nb_dof();
  size_type i2 = mf_u3.nb_dof();
  size_type i3 = mf_theta.nb_dof();
  
  getfem::interpolation(mf_ut, mf_rhs,
			gmm::sub_vector(U, gmm::sub_interval(0, i1)), V);
  mf_rhs.set_qdim(2);
  scalar_type l2 = gmm::sqr(getfem::asm_L2_norm(mim, mf_rhs, V));
  scalar_type h1 = gmm::sqr(getfem::asm_H1_norm(mim, mf_rhs, V));
  scalar_type linf = gmm::vect_norminf(V);
  mf_rhs.set_qdim(1);
  cout << "L2 error = " << sqrt(l2) << endl
       << "H1 error = " << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;

  getfem::interpolation(mf_theta, mf_rhs,
			gmm::sub_vector(U, gmm::sub_interval(i1+i2, i3)), V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i) {
    gmm::add(gmm::scaled(theta_exact(mf_rhs.point_of_dof(i)), -1.0),
	     gmm::sub_vector(V, gmm::sub_interval(i*2, 2)));
  }
  mf_rhs.set_qdim(2);
  l2 += gmm::sqr(getfem::asm_L2_norm(mim, mf_rhs, V));
  h1 += gmm::sqr(getfem::asm_H1_norm(mim, mf_rhs, V));
  linf = std::max(linf, gmm::vect_norminf(V));
  mf_rhs.set_qdim(1);
  cout << "L2 error = " << sqrt(l2) << endl
       << "H1 error = " << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;

  gmm::resize(V, mf_rhs.nb_dof());
  getfem::interpolation(mf_u3, mf_rhs,
			gmm::sub_vector(U, gmm::sub_interval(i1, i2)), V);

  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= u3_exact(mf_rhs.point_of_dof(i));

  l2 += gmm::sqr(getfem::asm_L2_norm(mim, mf_rhs, V));
  h1 += gmm::sqr(getfem::asm_H1_norm(mim, mf_rhs, V));
  linf = std::max(linf, gmm::vect_norminf(V));

  cout.precision(16);
  cout << "L2 error = " << sqrt(l2) << endl
       << "H1 error = " << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;
  
       // stockage de l'erreur H1 :
   if (PARAM.int_value("SAUV")){
      std::ofstream f_out("errH1.don");
      if (!f_out) throw std :: runtime_error("Impossible to open file") ;
      f_out << sqrt(h1) <<  "\n" ;
   }
       
  
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool plate_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  
  cout << "Number of dof for ut: " << mf_ut.nb_dof() << endl;
  cout << "Number of dof for u3: " << mf_u3.nb_dof() << endl;
  cout << "Number of dof for theta: " << mf_theta.nb_dof() << endl;

  getfem::mdbrick_abstract<> *ELAS, *SIMPLE(0);

  // Linearized plate brick.
  getfem::mdbrick_isotropic_linearized_plate<>
    ELAS1(mim, mim_subint, mf_ut, mf_u3, mf_theta, lambda,
	  mu, epsilon);

  if (mitc) ELAS1.set_mitc();
  
  getfem::mdbrick_mixed_isotropic_linearized_plate<>
    ELAS2(mim, mf_ut, mf_u3, mf_theta, lambda, mu, epsilon,
	  symmetrized);

  if (mixed) ELAS = &ELAS2; else ELAS = &ELAS1;

  // Defining the surface source term.
  if (study_flag == 1 ){
     cout << "Attention : l'intensité de la pression verticale " ;
     cout << "a été choisie pour que le déplacement maximal soit unitaire." ;
     cout << "Pour annuler cette option, faire STUDY_FLAG = 0\n" ;
     switch(sol_ref) {
           case 0 : 
                pressure  = 128. * mu * epsilon * epsilon * epsilon ;
                pressure /= 1.25 + 8. * epsilon * epsilon ; 
                break ;
           case 1 :  
                pressure  = 128. * mu * epsilon * epsilon * epsilon ;
                pressure /= 0.25 + 8. * epsilon * epsilon ; 
                break ;
	   case 3 :
	        pressure  = 32. * mu * epsilon * epsilon * epsilon ;
		pressure /= 3.  + 8. * epsilon * epsilon ;
           default : 
	        break ;	
	}
  }
  plain_vector F(nb_dof_rhs * 3); 
  plain_vector M(nb_dof_rhs * 2);
  if (sol_ref == 2) {
    base_small_vector P(2) ;
     scalar_type E, nu, sx, sy, cx, cy, s2x, s2y, c2x, c2y ;
     E = 4.*mu*(mu+lambda) / (2. * mu + lambda);
     nu = lambda / (2. * mu + lambda);
     for (size_type i = 0; i < nb_dof_rhs; ++i) {
       P   = mf_rhs.point_of_dof(i);
       sx  = sin(M_PI*P[0]) ;
       cx  = cos(M_PI*P[0]) ;
       sy  = sin(M_PI*P[1]) ;
       cy  = cos(M_PI*P[1]) ;
       c2x = cos(2.*M_PI*P[0]) ;
       c2y = cos(2.*M_PI*P[1]) ;
       s2x = sin(2.*M_PI*P[0]) ;
       s2y = sin(2.*M_PI*P[1]) ;
       F[3*i+2] = 2. * epsilon * E * M_PI * M_PI * eta *
                  ( sy * sy * c2x + sx * sx * c2y ) / ( 1. + nu ) ;
       M[2*i]   = (epsilon * epsilon * epsilon * E * M_PI * s2x / 3. / (1. + nu))
	 * ( (4. * M_PI * M_PI * (1. + eta) * (2. * c2y - 1.) / (1.- nu))  
	      - 3. * eta  * sy * sy / epsilon / epsilon ) ;
       M[2*i+1] = (epsilon * epsilon * epsilon * E * M_PI * s2y / 3. / (1. + nu))
	 * (  (4. * M_PI * M_PI * (1. + eta) * (2. * c2x - 1.) / (1.- nu))  
	       - 3. * eta  * sx * sx / epsilon / epsilon ) ;
     }
  }
  else  // sol_ref = 0 ou 1 ou 3 ou 4: pression verticale uniforme
     for (size_type i = 0; i < nb_dof_rhs; ++i) {
       F[3*i+2] = pressure;
  }
  
  getfem::mdbrick_plate_source_term<> VOL_F(*ELAS, mf_rhs, F, M);
  
  getfem::mdbrick_plate_simple_support<> SIMPLE0
    (VOL_F, SIMPLY_FIXED_BOUNDARY_NUM, 0, getfem::AUGMENTED_CONSTRAINTS);
  getfem::mdbrick_plate_clamped_support<> SIMPLE1
    (VOL_F, SIMPLY_FIXED_BOUNDARY_NUM, 0, getfem::AUGMENTED_CONSTRAINTS);
    
  if (sol_ref == 0) SIMPLE = &SIMPLE0 ;
  if (sol_ref == 1) SIMPLE = &SIMPLE1 ;
  if (sol_ref == 2) SIMPLE = &SIMPLE1 ;
  if (sol_ref == 3) SIMPLE = &SIMPLE1 ;
  if (sol_ref == 4) SIMPLE = &SIMPLE0 ;
  
  
  getfem::mdbrick_plate_closing<> final_model(*SIMPLE, 0, 1);
  

  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);


  gmm::resize(U, mf_u3.nb_dof());
  if (mixed) gmm::copy(ELAS2.get_u3(MS), U);
  else gmm::copy(ELAS1.get_u3(MS), U);
  
  if (PARAM.int_value("VTK_EXPORT")) {
    cout << "export to " << datafilename + ".vtk" << "..\n";
    getfem::vtk_export exp(datafilename + ".vtk",
			   PARAM.int_value("VTK_EXPORT")==1);
    exp.exporting(mf_u3); 
    exp.write_point_data(mf_u3, U, "plate_normal_displacement");
    cout << "export done, you can view the data file with (for example)\n"
       "mayavi -d " << datafilename << ".vtk -f "
       "WarpScalar -m BandedSurfaceMap -m Outline\n";
//     cout << "export done, you can view the data file with (for example)\n"
//       "mayavi -d " << datafilename << ".vtk -f ExtractVectorNorm -f "
//       "WarpVector -m BandedSurfaceMap -m Outline\n";
  }
    if (PARAM.int_value("DX_EXPORT")) {
    cout << "export to " << datafilename + ".dx" << "..\n";
    getfem::dx_export exp(datafilename + ".dx",
			   PARAM.int_value("DX_EXPORT")==1);
    exp.exporting(mf_u3); 
    exp.write_point_data(mf_u3, U, "plate_normal_displacement");
  }
  // Solution extraction
  gmm::resize(U, ELAS->nb_dof());
  if (mixed) gmm::copy(ELAS2.get_solution(MS), U);
  else gmm::copy(ELAS1.get_solution(MS), U);

  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.  

  try {    
    plate_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    if ((p.study_flag != 1)&&((p.sol_ref == 0) || (p.sol_ref ==1)))
    p.pressure *= p.epsilon * p.epsilon * p.epsilon;
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U;
    if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");
    p.compute_error(U);
    
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
