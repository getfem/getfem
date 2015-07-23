/*===========================================================================

 Copyright (C) 2002-2015 Yves Renard, Michel Salaün.

 This file is a part of GetFEM++

 GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

/**
   @file plate.cc
   @brief Linear Elastostatic plate problem.
   
   This program is used to check that getfem++ is working. This is
   also a good example of use of GetFEM++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_linearized_plates.h"
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some GetFEM++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type;
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
  scalar_type E, nu;    /* Lamé coefficients.                           */
  scalar_type epsilon;       /* thickness of the plate.                      */
  scalar_type pressure;
  scalar_type residual;      /* max residual for the iterative solvers         */
  scalar_type LX , LY ;      // default : LX = LY = 1
  bool mitc;
  int sol_ref;               // sol_ref = 0 : simple support on the vertical edges
                             // sol_ref = 1 : homogeneous on the vertical edges
                             // sol_ref = 2 : homogeneous on the 4 vertical
                             //       edges with solution u3 = sin²(x)*sin²(y)
  scalar_type eta;           // useful only if sol_ref == 2 :
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

  bool solve(plain_vector &Ut, plain_vector &U3, plain_vector &THETA);
  void init(void);
  void compute_error(plain_vector &Ut, plain_vector &U3, plain_vector &THETA);
  plate_problem(void) : mim(mesh), mim_subint(mesh), mf_ut(mesh), mf_u3(mesh),
			mf_theta(mesh), mf_rhs(mesh), mf_coef(mesh) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void plate_problem::init(void) {
  std::string MESH_FILE = PARAM.string_value("MESH_FILE");
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
  size_type N;
    bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  if (!MESH_FILE.empty()) {
    cout << "MESH_FILE=" << MESH_FILE << "\n";
    mesh.read_from_file(MESH_FILE);
    MESH_TYPE = bgeot::name_of_geometric_trans
      (mesh.trans_of_convex(mesh.convex_index().first_true()));
    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    N = mesh.dim();
  } else {
  N = pgt->dim();
  GMM_ASSERT1(N == 2, "For a plate problem, N should be 2");
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Number of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  }

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  mitc = (PARAM.int_value("MITC", "Mitc version ?") != 0);
  

  cout << "MITC = " ;
  if (mitc) cout << "true \n" ; else cout << "false \n" ;
  sol_ref = int(PARAM.int_value("SOL_REF") ) ;
  study_flag = int(PARAM.int_value("STUDY_FLAG") ) ;
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
     scalar_type /*det_Jmn, */ A, B, e2, Pmn ;
     E = 4.*mu*(mu+lambda) / (2. * mu + lambda);
     nu = lambda / (2. * mu + lambda);
     e2 = epsilon * epsilon / 4. ;
     for(size_type i = 0 ; i < N_Four ; i++) {
        for(size_type j = 0 ; j < N_Four ; j++) {
	   A = scalar_type(j + 1) * M_PI / LX ; 
	   B = scalar_type(i + 1) * M_PI / LY ; 
	   Jmn(0, 0) = 2. * A * A / (1. - nu) + B * B + 3. / e2   ;
	   Jmn(0, 1) = A * B * (1. +nu) / (1. - nu) ;
	   Jmn(0, 2) = A * 3. / e2 ;
	   Jmn(1, 0) = A * B * (1. +nu) / (1. - nu) ;
	   Jmn(1, 1) = 2. * B * B / (1. - nu) + A * A + 3. / e2  ;
	   Jmn(1, 2) = B * 3. / e2 ;
	   Jmn(2, 0) = - A ;
	   Jmn(2, 1) = - B ;
	   Jmn(2, 2) = A * A + B * B ;
	   gmm::scale(Jmn,  - E*(epsilon/2.) / (1. + nu) ) ;
	   
	   // calcul du développement de Fourrier du chargement :
	   if ( ( (i + 1) % 2 == 1 ) && ( (j + 1) % 2 == 1) ) {
	      Pmn =  16. * pressure / ( scalar_type(i + 1) * scalar_type(j + 1) * M_PI * M_PI) ; }
	   else {
	      Pmn = 0. ; }	      
	   Bmn[0] = 0. ;
	   Bmn[1] = 0. ;
	   Bmn[2] = Pmn ;
	   gmm::lu_solve(Jmn, Xmn, Bmn) ;
	   theta1_Four(i, j) = Xmn[0] ;
	   theta1_Four(i, j) = Xmn[1] ;
               u3_Four(i, j) = Xmn[2] ;	       
	   }
        }
     }

  mf_ut.set_qdim(dim_type(N));
  mf_theta.set_qdim(dim_type(N));

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
     theta[0] = - (-pressure / (32. * mu * epsilon * epsilon * epsilon / 8.))
              * (4. * pow(P[0] - .5, 3.) - 3 * (P[0] - .5));
     theta[1] = 0.;
   }
  if (sol_ref == 1) { // encastrement aux 2 bords
     theta[0] = - (-pressure / (16. * mu * epsilon * epsilon * epsilon / 8.))
              * P[0] * ( 2.*P[0]*P[0] - 3.* P[0] + 1. ) ;
     theta[1] = 0.;
     }
  if (sol_ref == 2) { // encastrement aux 4 bords et sols vraiment 2D, non polynomiale
     theta[0] = (1. + eta) * M_PI * sin(2.*M_PI*P[0]) * sin(M_PI*P[1])* sin(M_PI*P[1]) ;
     theta[1] = (1. + eta) * M_PI * sin(2.*M_PI*P[1]) * sin(M_PI*P[0])* sin(M_PI*P[0]) ;
     }
  if (sol_ref == 3) { // plaque cantilever
     theta[0] = - (- 3. * pressure / (8. * mu * epsilon * epsilon * epsilon / 8.))
              * P[0] * ( 0.25 * P[0] * P[0] - P[0] + 1. ) ;
     theta[1] = 0.;      
     }
  if (sol_ref == 4) { // bord entier en appui simple
     theta[0] = 0. ;
     theta[1] = 0. ;
     for(size_type i = 0 ; i < N_Four ; i ++) {
        for(size_type j = 0 ; j < N_Four ; j ++) {
	   theta[0] -= theta1_Four(i, j) * cos( scalar_type(j + 1) * M_PI * P[0] / LX ) * sin( scalar_type(i + 1) * M_PI * P[1] / LY ) ;
	   theta[0] -= theta2_Four(i, j) * sin( scalar_type(j + 1) * M_PI * P[0] / LX ) * cos( scalar_type(i + 1) * M_PI * P[1] / LY ) ;
           }
        }
     }
  return theta;
}

scalar_type plate_problem::u3_exact(base_node P) {
  switch(sol_ref) {
  case 0 : return (pressure / (32. * mu * epsilon * epsilon * epsilon / 8.))
       * P[0] * (P[0] - 1.)
       * (gmm::sqr(P[0] - .5) -1.25-(8.* epsilon*epsilon / 4.));
       break ;
  case 1 : return (pressure /(32.* mu * epsilon * epsilon * epsilon / 8.))
       * P[0] * (P[0] - 1.)
       * ( P[0] * P[0] - P[0] - 8. * epsilon *epsilon / 4.) ;
       break ;
  case 2 : return  gmm::sqr(sin(M_PI*P[0])) * gmm::sqr(sin(M_PI*P[1]));
       break ;
  case 3 : return (3. * pressure / (4. * mu * epsilon * epsilon * epsilon / 8. ))
       * P[0] * ( P[0] * P[0] * P[0] / 24. - P[0] * P[0] / 6. + P[0] / 4. 
                  - (epsilon * epsilon / 4.) * P[0] / 3.
                  + 2. * (epsilon * epsilon / 4.) / 3.) ;
      break ;
  case 4 : 
       scalar_type u3_local ;
       u3_local = 0. ;
       for(size_type i = 0 ; i < N_Four ; i ++) {
          for(size_type j = 0 ; j < N_Four ; j ++) 
	     u3_local += u3_Four(i, j) * sin( scalar_type(j + 1) * M_PI * P[0] / LX ) * sin( scalar_type(i + 1) * M_PI * P[1] / LY ) ;
	     }
       return (u3_local) ;
       break ; 
  default : GMM_ASSERT1(false, "indice de solution de référence incorrect");
  }
}


/* compute the error with respect to the exact solution */
void plate_problem::compute_error(plain_vector &Ut, plain_vector &U3, plain_vector &THETA) {
  cout.precision(16);
  if (PARAM.int_value("SOL_EXACTE") == 1) {
    gmm::clear(Ut); gmm::clear(U3); gmm::clear(THETA); 
  }

  std::vector<scalar_type> V(mf_rhs.nb_dof()*2);
  
  getfem::interpolation(mf_ut, mf_rhs, Ut, V);
  mf_rhs.set_qdim(2);
  scalar_type l2 = gmm::sqr(getfem::asm_L2_norm(mim, mf_rhs, V));
  scalar_type h1 = gmm::sqr(getfem::asm_H1_norm(mim, mf_rhs, V));
  scalar_type linf = gmm::vect_norminf(V);
  mf_rhs.set_qdim(1);
  cout << "L2 error = " << sqrt(l2) << endl
       << "H1 error = " << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;

  getfem::interpolation(mf_theta, mf_rhs, THETA, V);
  GMM_ASSERT1(!mf_rhs.is_reduced(),
	      "To be adapted, use interpolation_function");
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i) {
    gmm::add(gmm::scaled(theta_exact(mf_rhs.point_of_basic_dof(i)), -1.0),
	     gmm::sub_vector(V, gmm::sub_interval(i*2, 2)));
  }
  mf_rhs.set_qdim(2);
  l2 += gmm::sqr(getfem::asm_L2_norm(mim, mf_rhs, V));
  h1 += gmm::sqr(getfem::asm_H1_semi_norm(mim, mf_rhs, V));
  linf = std::max(linf, gmm::vect_norminf(V));
  mf_rhs.set_qdim(1);
  cout << "L2 error theta:" << sqrt(l2) << endl
       << "H1 error theta:" << sqrt(h1) << endl
       << "Linfty error = " << linf << endl;

  gmm::resize(V, mf_rhs.nb_dof());
  getfem::interpolation(mf_u3, mf_rhs, U3, V);

  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= u3_exact(mf_rhs.point_of_basic_dof(i));

  l2 = gmm::sqr(getfem::asm_L2_norm(mim, mf_rhs, V));
  h1 = gmm::sqr(getfem::asm_H1_semi_norm(mim, mf_rhs, V));
  linf = std::max(linf, gmm::vect_norminf(V));

  cout.precision(16);
  cout << "L2 error u3:" << sqrt(l2) << endl
       << "H1 error u3:" << sqrt(h1) << endl
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

bool plate_problem::solve(plain_vector &Ut, plain_vector &U3, plain_vector &THETA) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  
  cout << "Number of dof for ut: " << mf_ut.nb_dof() << endl;
  cout << "Number of dof for u3: " << mf_u3.nb_dof() << endl;
  cout << "Number of dof for theta: " << mf_theta.nb_dof() << endl;

  E = 4.*mu*(mu+lambda) / (2. * mu + lambda);
  nu = lambda / (2. * mu + lambda);
  scalar_type kappa = 5./6.;

  getfem::model md;
  md.add_fem_variable("ut", mf_ut);
  md.add_fem_variable("u3", mf_u3);
  md.add_fem_variable("theta", mf_theta);

  // Linearized plate brick.
  md.add_initialized_scalar_data("E", E);
  md.add_initialized_scalar_data("nu", nu);
  md.add_initialized_scalar_data("lambda", lambda);
  md.add_initialized_scalar_data("mu", mu);
  md.add_initialized_scalar_data("epsilon", epsilon);
  md.add_initialized_scalar_data("kappa", kappa);
  getfem::add_Mindlin_Reissner_plate_brick(md, mim, mim_subint, "u3", "theta",
                                           "E", "nu", "epsilon", "kappa",
                                           (mitc) ? 2 : 1);
  getfem::add_isotropic_linearized_elasticity_brick(md, mim, "ut", "lambda", "mu");

  // Defining the surface source term.
  if (study_flag == 1 ){
     cout << "Attention : l'intensité de la pression verticale " ;
     cout << "a été choisie pour que le déplacement maximal soit unitaire." ;
     cout << "Pour annuler cette option, faire STUDY_FLAG = 0\n" ;
     switch(sol_ref) {
           case 0 : 
                pressure  = 128. * mu * epsilon * epsilon * epsilon / 8. ;
                pressure /= 1.25 + 8. * epsilon * epsilon / 4. ; 
                break ;
           case 1 :  
                pressure  = 128. * mu * epsilon * epsilon * epsilon / 8. ;
                pressure /= 0.25 + 8. * epsilon * epsilon / 4. ; 
                break ;
	   case 3 :
	        pressure  = 32. * mu * epsilon * epsilon * epsilon / 8.;
		pressure /= 3.  + 8. * epsilon * epsilon / 4.;
           default : 
	        break ;	
	}
  }
  plain_vector F(nb_dof_rhs); 
  plain_vector M(nb_dof_rhs * 2);
  if (sol_ref == 2) {
    base_small_vector P(2) ;
     scalar_type sx, sy, cx, cy, s2x, s2y, c2x, c2y ;
     E = 4.*mu*(mu+lambda) / (2. * mu + lambda);
     nu = lambda / (2. * mu + lambda);
     for (size_type i = 0; i < nb_dof_rhs; ++i) {
       P   = mf_rhs.point_of_basic_dof(i);
       sx  = sin(M_PI*P[0]) ;
       cx  = cos(M_PI*P[0]) ;
       sy  = sin(M_PI*P[1]) ;
       cy  = cos(M_PI*P[1]) ;
       c2x = cos(2.*M_PI*P[0]) ;
       c2y = cos(2.*M_PI*P[1]) ;
       s2x = sin(2.*M_PI*P[0]) ;
       s2y = sin(2.*M_PI*P[1]) ;
       F[i] = 2. * (epsilon / 2.) * E * M_PI * M_PI * eta *
                  ( sy * sy * c2x + sx * sx * c2y ) / ( 1. + nu ) ;
       M[2*i]   = -((epsilon * epsilon * epsilon / 8.) * E * M_PI * s2x / 3. / (1. + nu))
	 * ( (4. * M_PI * M_PI * (1. + eta) * (2. * c2y - 1.) / (1.- nu))  
             - 3. * eta  * sy * sy / (epsilon/2.) / (epsilon/2.) ) ;
       M[2*i+1] = -((epsilon * epsilon * epsilon/8.) * E * M_PI * s2y / 3. / (1. + nu))
	 * (  (4. * M_PI * M_PI * (1. + eta) * (2. * c2x - 1.) / (1.- nu))  
              - 3. * eta  * sx * sx / (epsilon/2.) / (epsilon/2.) ) ;
     }
  }
  else  // sol_ref = 0 ou 1 ou 3 ou 4: pression verticale uniforme
     for (size_type i = 0; i < nb_dof_rhs; ++i) {
       F[i] = pressure;
  }
  
  md.add_initialized_fem_data("VF", mf_rhs, F);
  getfem::add_source_term_brick(md, mim, "u3", "VF");
  md.add_initialized_fem_data("VM", mf_rhs, M);
  getfem::add_source_term_brick(md, mim, "theta", "VM");
  
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "u3", mf_u3, SIMPLY_FIXED_BOUNDARY_NUM);
  
  if (sol_ref == 1 || sol_ref == 2 || sol_ref == 3)
    getfem::add_Dirichlet_condition_with_multipliers
      (md, mim, "theta", mf_ut, SIMPLY_FIXED_BOUNDARY_NUM);
  

  // Generic solve.
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(md, iter);


  gmm::resize(Ut, mf_ut.nb_dof());
  gmm::copy(md.real_variable("ut"), Ut);
  gmm::resize(U3, mf_u3.nb_dof());
  gmm::copy(md.real_variable("u3"), U3);
  gmm::resize(THETA, mf_theta.nb_dof());
  gmm::copy(md.real_variable("theta"), THETA);
  
  if (PARAM.int_value("VTK_EXPORT")) {
    cout << "export to " << datafilename + ".vtk" << "..\n";
    getfem::vtk_export exp(datafilename + ".vtk",
			   PARAM.int_value("VTK_EXPORT")==1);
    exp.exporting(mf_u3); 
    exp.write_point_data(mf_u3, U3, "plate_normal_displacement");
    cout << "export done, you can view the data file with (for example)\n"
       "mayavi2 -d " << datafilename << ".vtk -f "
       "WarpScalar -m Surface -m Outline\n";
//     cout << "export done, you can view the data file with (for example)\n"
//       "mayavi -d " << datafilename << ".vtk -f ExtractVectorNorm -f "
//       "WarpVector -m BandedSurfaceMap -m Outline\n";
  }
    if (PARAM.int_value("DX_EXPORT")) {
    cout << "export to " << datafilename + ".dx" << ".\n";
    getfem::dx_export exp(datafilename + ".dx",
			   PARAM.int_value("DX_EXPORT")==1);
    exp.exporting(mf_u3); 
    exp.write_point_data(mf_u3, U3, "plate_normal_displacement");
  }
  

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
    p.pressure *= p.epsilon * p.epsilon * p.epsilon / 8.;
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector Ut, U3, THETA;
    if (!p.solve(Ut, U3, THETA)) GMM_ASSERT1(false, "Solve has failed");
    p.compute_error(Ut, U3, THETA);
    
  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
