// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2002-2007 Yves Renard, Julien Pommier.
//
// This file is a part of GETFEM++
//
// Getfem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301,
// USA.
//
//========================================================================


/**
 * Linear Elastostatic plate problem with a through crack.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_linearized_plates.h"
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_mesh_fem_sum.h"

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


// Parameters for the exact solution  ---------------------------


scalar_type nu = 0.3 ;

scalar_type A0 = 1. ;

scalar_type A1 = 0.1 ;
scalar_type B1 = 0.0 ;
scalar_type C1 = 3.* A1 ;   
scalar_type D1 = 0.0 ;  

scalar_type A2 = 0.1 ;
scalar_type B2 = 0.0 ;
scalar_type C2 = 3.* A2 + nu * (D1 - B1) ; 
scalar_type D2 = B2 ;  
 

// Singular functions for Mindlin ------------------------------

struct mindlin_singular_functions : public getfem::global_function, public getfem::context_dependencies {
  size_type l;             // singular function number
  const getfem::level_set &ls;
  mutable getfem::mesher_level_set mls0, mls1;
  mutable size_type cv;
  
  void update_mls(size_type cv_) const { 
    if (cv_ != cv) { 
      cv=cv_; 
      mls0=ls.mls_of_convex(cv, 0); 
      mls1=ls.mls_of_convex(cv, 1); 
    }
  }


  scalar_type sing_function(scalar_type x, scalar_type y) const {
    scalar_type r = sqrt(x*x + y*y);
    scalar_type theta = atan2(y,x); 
 
    switch (l) {
    case 0: {
      return sqrt(r)*cos(3.0 * theta/2);
    } break;
    case 1: {
      return sqrt(r)*sin(3.0 * theta/2);
    } break;
    case 2: {
      return sqrt(r)*cos(theta/2);
    } break;
    case 3: {
      return sqrt(r)*sin(theta/2);
    } break;
    case 4: {   // usefull for isotropic_linear_elasticity_2D
      return sqrt(r)*sin(theta/2)*cos(theta);
    } break;
    case 5: {   // usefull for isotropic_linear_elasticity_2D
      return sqrt(r)*cos(theta/2)*cos(theta);
    } break;
    default: assert(0); 
    }
    return 0;
  }

  /* derivative in the levelset coordinates */
  void sing_function_grad(scalar_type x, scalar_type y, base_small_vector &g) const {
    g.resize(2);
    scalar_type r = sqrt(x*x + y*y);
    scalar_type theta = atan2(y,x);

    switch (l) {
    case 0: {
      g[0] = (   cos(theta/2.0) - cos(5.0/2.0*theta) / 2.0) / sqrt(r) ; 
      g[1] = ( - sin(theta/2.0) - sin(5.0/2.0*theta) / 2.0) / sqrt(r) ; 
    } break;
    case 1: {
      g[0] = (   sin(theta/2.0) - sin(5.0/2.0*theta) / 2.0) / sqrt(r) ; 
      g[1] = (   cos(theta/2.0) + cos(5.0/2.0*theta) / 2.0) / sqrt(r) ;
    } break;
    case 2: {
      g[0] =   cos(theta/2.0)  / 2.0   / sqrt(r) ; 
      g[1] =   sin(theta/2.0)  / 2.0   / sqrt(r) ; 
    } break;
    case 3: {
      g[0] = - sin(theta/2.0)  / 2.0   / sqrt(r) ; 
      g[1] =   cos(theta/2.0)  / 2.0   / sqrt(r) ; 
    } break;
    case 4: {  // usefull for isotropic_linear_elasticity_2D   -> expression to check with maple !
      g[0] =   (   sin(theta/2.0) - sin(5.0/2.0*theta) )  / 4.0   / sqrt(r) ; 
      g[1] =   (   cos(theta/2.0) + cos(5.0/2.0*theta) )  / 4.0   / sqrt(r) ; 
    } break;
    case 5: {  // usefull for isotropic_linear_elasticity_2D   -> expression to check with maple !
      g[0] = (  3./2.* cos(theta/2.0) - cos(5.0/2.0*theta) / 2.0)   / 2.0   / sqrt(r) ; 
      g[1] = (       - sin(theta/2.0) - sin(5.0/2.0*theta)      )   / 4.0   / sqrt(r) ; 
    } break;
    default: assert(0); 
    }
  }
   
 
   
  
  virtual scalar_type val(const getfem::fem_interpolation_context& c) const {

    assert(ls.get_mesh_fem().convex_index().is_in(c.convex_num()));
    update_mls(c.convex_num());
    scalar_type x = mls1(c.xref()), y = mls0(c.xref());
    scalar_type v = sing_function(x, y);
    return v;
  }

  virtual void grad(const getfem::fem_interpolation_context& c,
		    base_small_vector &v) const {
    update_mls(c.convex_num());
    size_type P = c.xref().size();
    base_small_vector dx(P), dy(P), dfr(2);
    scalar_type x = mls1.grad(c.xref(), dx), y = mls0.grad(c.xref(), dy);
    if (x*x + y*y < 1e-20) {
      cerr << "Warning, point very close to the singularity. xreal = "
	   << c.xreal() << ", x_crack = " << x << ", y_crack=" << y << "\n";
    } 
    sing_function_grad(x, y, dfr); // gradient in the coord of the levelset

    base_small_vector dref = dfr[0]*dx + dfr[1]*dy; // gradient in the coordinate the ref convex.
    gmm::mult(c.B(), dref, v); // gradient in the coords of the real element.
  }
  

    
  void update_from_context(void) const { cv =  size_type(-1); }

  mindlin_singular_functions(size_type l_, const getfem::level_set &ls_) : l(l_), ls(ls_) {
    cv = size_type(-1);
    this->add_dependency(ls);
  }

};

getfem::pglobal_function mindlin_crack_singular(size_type i, const getfem::level_set &ls){ 
  return new mindlin_singular_functions(i, ls);
}


struct exact_solution {
  getfem::mesh_fem_global_function mf_ut, mf_u3, mf_theta;
  getfem::base_vector UT ;
  getfem::base_vector U3 ; 
  getfem::base_vector THETA;
  
  exact_solution(getfem::mesh &me) : mf_ut(me), mf_u3(me), mf_theta(me) {}
  
  void init(getfem::level_set &ls) {
    std::vector<getfem::pglobal_function> cfun_ut(4), cfun_u3(1), cfun_theta(4) ;
    for (unsigned j=0; j < 4; ++j) {
      cfun_ut[j] = mindlin_crack_singular(j+2, ls) ;
      cfun_theta[j] = mindlin_crack_singular(j, ls) ;
    }
    cfun_u3[0] = mindlin_crack_singular(3,ls) ;
    
    // Initialising  ut
    mf_ut.set_qdim(2) ;
    mf_ut.set_functions(cfun_ut);   
    UT.resize(8); assert(mf_ut.nb_dof() == 8);
    for( unsigned i = 0; i< 8 ; ++i)
        UT[i] = 0. ;
      
    // Initialising  theta
    mf_theta.set_qdim(2) ;
    mf_theta.set_functions(cfun_theta);   
    THETA.resize(8); assert(mf_theta.nb_dof() == 8);
    THETA[0] = A1 ;
    THETA[1] = A2 ;
    THETA[2] = B1 ;
    THETA[3] = B2 ;
    THETA[4] = C1 ;
    THETA[5] = C2 ;
    THETA[6] = D1 ;
    THETA[7] = D2 ;
    
    // Initialising u3
    mf_u3.set_functions(cfun_u3) ;
    U3.resize(1) ;
    U3[0] = A0 ;  
  }
  
};

/*                                                                            */
/*   -------------  Struct for the (cracked) plate problem ---------------    */ 
/*                                                                            */

struct crack_mindlin_problem{

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_level_set mls;       /* the integration methods.              */
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_ut, mf_pre_u3, mf_pre_theta ; 
  getfem::mesh_fem_level_set mfls_ut, mfls_u3, mfls_theta ; 
  getfem::mesh_fem_global_function mf_sing_ut, mf_sing_u3, mf_sing_theta ;
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_ut_product, mf_u3_product, mf_theta_product;
  getfem::mesh_fem_sum mf_ut_sum, mf_u3_sum, mf_theta_sum;

  getfem::mesh_fem& mf_ut() { return mf_ut_sum; }
  getfem::mesh_fem& mf_u3() { return mf_u3_sum; }
  getfem::mesh_fem& mf_theta() { return mf_theta_sum; }

  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_mult_ut, mf_mult_u3, mf_mult_theta ;  /* mesh_fem for the md_brick_dirichlet */ 
  getfem::constraints_type dirichlet_version;
  getfem::level_set ls;      /* The two level sets defining the crack.       */
 
  scalar_type residual;       /* max residual for the iterative solvers        */
  scalar_type cutoff_radius, enr_area_radius;
  int enrichment_option;
  std::string datafilename;
  bgeot::md_param PARAM;
  size_type mitc ;
  size_type sol_ref ;          /* sol_ref = 0 : plate subject to vertical pressure, 
                                                positive or negative on each side of the crack.
		                  sol_ref = 1 : same as sol_ref = 2 but a side is clamped on a 
				                side and the other is subjected to a 2 M0 
						horizontal moment.   
			          sol_ref = 2 : plate subject to horizontal moment M0 on the
				                2 vertical bounds of the plate      
				  sol_ref = 3 : non-homogeneous Dirichlet conditions on the
				                4 boundaries of the plate. The exact solution
						is equal to the singularities.   */

  scalar_type h_crack_length ;  // half-crack-length

  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type epsilon, pressure ;
  
  exact_solution exact_sol;
  
  // methods
  
  bool solve(plain_vector &UT, plain_vector &U3, plain_vector &THETA);
  void init(void);
  void compute_error(plain_vector &UT, plain_vector &U3, plain_vector &THETA);
  crack_mindlin_problem(void) : mls(mesh), mim(mls),
  			mf_pre_ut(mesh), mf_pre_u3(mesh), mf_pre_theta(mesh), 
			mfls_ut(mls, mf_pre_ut), mfls_u3(mls, mf_pre_u3), mfls_theta(mls, mf_pre_theta),
			mf_sing_ut(mesh),mf_sing_u3(mesh),mf_sing_theta(mesh),
			mf_partition_of_unity(mesh),
			mf_ut_product(mf_partition_of_unity, mf_sing_ut),
			mf_u3_product(mf_partition_of_unity, mf_sing_u3),
			mf_theta_product(mf_partition_of_unity, mf_sing_theta),
			mf_ut_sum(mesh), mf_u3_sum(mesh), mf_theta_sum(mesh),
			mf_rhs(mesh), mf_mult_ut(mesh), mf_mult_u3(mesh), mf_mult_theta(mesh), 
			ls(mesh, 1, true), exact_sol(mesh) {} 

};

void crack_mindlin_problem::init(void) {
  
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE_UT  = PARAM.string_value("FEM_TYPE_UT","FEM name");
  std::string FEM_TYPE_U3  = PARAM.string_value("FEM_TYPE_U3","FEM name");
  std::string FEM_TYPE_THETA = PARAM.string_value("FEM_TYPE_THETA","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
				      "Enrichment option");

  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE_UT="  << FEM_TYPE_UT << "\n";
  cout << "FEM_TYPE_U3="  << FEM_TYPE_U3 << "\n";
  cout << "FEM_TYPE_THETA="  << FEM_TYPE_THETA << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  
  // build the mesh :
  bgeot::pgeometric_trans pgt = 
  bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  GMM_ASSERT1(N == 2, "For a plate problem, N should be 2");
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Number of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);

  base_small_vector tt(N); 
  if (sol_ref > 0) tt[0] = -0.5 ;
  tt[1] = -0.5;
  mesh.translation(tt); 
  
  // setting parameters :
  epsilon = PARAM.real_value("EPSILON", "thickness") ; 
  mitc = (PARAM.int_value("MITC") != 0);
  pressure = PARAM.real_value("PRESSURE") ;
  int dv = PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version");
  dirichlet_version = getfem::constraints_type(dv);
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  cutoff_radius = PARAM.real_value("CUTOFF", "Cutoff");
  h_crack_length = PARAM.real_value("HALF_CRACK_LENGTH") ;
  sol_ref = PARAM.int_value("SOL_REF") ;
  
  mf_ut().set_qdim(N);
  mf_theta().set_qdim(N);
  
  // set the finite element method on ut, u3, theta
  getfem::pfem pf_ut = getfem::fem_descriptor(FEM_TYPE_UT);
  getfem::pfem pf_u3 = getfem::fem_descriptor(FEM_TYPE_U3);
  getfem::pfem pf_theta = getfem::fem_descriptor(FEM_TYPE_THETA);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method sppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
    
  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  mim.set_simplex_im(sppi);
  mf_pre_ut.set_finite_element(mesh.convex_index(), pf_ut);
  mf_pre_u3.set_finite_element(mesh.convex_index(), pf_u3);
  mf_pre_theta.set_finite_element(mesh.convex_index(), pf_theta);
  mf_partition_of_unity.set_classical_finite_element(1);

  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  mf_mult_ut.set_qdim(2) ;
  mf_mult_theta.set_qdim(2) ;
  if (dirichlet_fem_name.size() == 0) {
    mf_mult_u3.set_finite_element(mesh.convex_index(), pf_u3) ;
    mf_mult_theta.set_finite_element(mesh.convex_index(), pf_theta) ; }
  else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult_ut.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
    mf_mult_u3.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
    mf_mult_theta.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
  }
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

  /* set boundary conditions : Dirichlet on the right face */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor it(border_faces); !it.finished(); ++it) {
    assert(it.is_face());
    base_node un = mesh.normal_of_face_of_convex(it.cv(), it.f());
    un /= gmm::vect_norm2(un);
    switch (sol_ref) {
    case 0 : {
	if ( un[0] >= 1. - 1.0E-7) { // new Dirichlet
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
	} else {
	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(), it.f());
	}
      }
    break ;
    case 1 : {
	if ( un[0] >= 1. - 1.0E-7) { // new Dirichlet
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());
	} else {
	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(), it.f());
	}
      }
    break ;
//    case 2 : {
// 	if ( gmm::abs(un[1]) <  1.0E-7) { // new Neumann
// 	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(), it.f());
// 	} else {
// 	base_small_vector index ;
// 	index = mesh.ind_points_of_convex(it.cv()) ;
// 	size_type i, j, found ;
// 	i = 0 ; j = 0 ; found = 0 ;
// 	for( i = 0 ; i < 3 ; i++ ) { 
// 	  for( j = i + 1 ; j < 4 ; j++ ) { 
// 	    if (  gmm::abs( mesh.points()[index[i]][1] ) >= 0.5 - 1.0E-7  && 
// 	          gmm::abs( mesh.points()[index[j]][1] ) >= 0.5 - 1.0E-7  && 
// 		  mesh.points()[index[i]][0] * mesh.points()[index[j]][0] <=  1.0E-7 ) {
// 	       mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f());       
// 	       found = 1 ;
// 	    }
// 	  }
// 	}
// 	if ( found == 0)  
// 	mesh.region(NEUMANN_BOUNDARY_NUM).add(it.cv(), it.f());
// 	}
//       }
//    break ;
    case 3 : // Dirichlet on the hole boundary of the domain
        mesh.region(DIRICHLET_BOUNDARY_NUM).add(it.cv(),it.f()); 
      break ;
    default : {
        cout << "wrong parameter for sol_ref. One should set sol_ref to 0, 1 or 3 \n" ;
      }
    break ;
    }
  }
  
  exact_sol.init(ls);
  
}

bool crack_mindlin_problem::solve(plain_vector &UT, plain_vector &U3, plain_vector &THETA) {


  // Setting the enrichment --------------

  size_type nb_dof_rhs = mf_rhs.nb_dof();
  //size_type N = mesh.dim();

  ls.reinit();  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
  if (sol_ref == 0 || sol_ref == 3) {
    ls.values(0)[d] = (ls.get_mesh_fem().point_of_dof(d))[1];
    ls.values(1)[d] =  (ls.get_mesh_fem().point_of_dof(d))[0];}
  if (sol_ref == 1){
    ls.values(0)[d] = (ls.get_mesh_fem().point_of_dof(d))[0];
    ls.values(1)[d] = gmm::abs( (ls.get_mesh_fem().point_of_dof(d))[1] ) - h_crack_length ;}
    }
  
  ls.touch();
  
  mls.adapt();
  mim.adapt();
  mfls_ut.adapt();
  mfls_u3.adapt();
  mfls_theta.adapt();
  
  cout << "setting singularities... \n" ;
  std::vector<getfem::pglobal_function> utfunc(4) ;
  std::vector<getfem::pglobal_function> u3func(1) ;
  std::vector<getfem::pglobal_function> theta_func(4) ;
  

  for (size_type i = 0 ; i < 4 ; ++i){
        utfunc[i] = mindlin_crack_singular(i+2, ls); 
	theta_func[i] = mindlin_crack_singular(i, ls);
  }  
  u3func[0] = mindlin_crack_singular(3, ls) ;
  
  mf_sing_ut.set_functions(utfunc);
  mf_sing_u3.set_functions(u3func);
  mf_sing_theta.set_functions(theta_func);
  
  switch (enrichment_option) {
  case 1 : 
    {
      mf_ut_sum.set_mesh_fems(mf_sing_ut, mfls_ut);
      mf_u3_sum.set_mesh_fems(mf_sing_u3, mfls_u3);
      mf_theta_sum.set_mesh_fems(mf_sing_theta, mfls_theta);
      }
    break;
  case 2 :
    {
      dal::bit_vector enriched_dofs;
      plain_vector X(mf_partition_of_unity.nb_dof());
      plain_vector Y(mf_partition_of_unity.nb_dof());
      getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			    ls.values(1), X);
      getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			    ls.values(0), Y);
      for (size_type j = 0; j < mf_partition_of_unity.nb_dof(); ++j) {
	if (gmm::sqr(X[j]) + gmm::sqr(Y[j]) <= gmm::sqr(enr_area_radius))
	  enriched_dofs.add(j);
      }
      if (enriched_dofs.card() < 3)
	GMM_WARNING0("There is " << enriched_dofs.card() <<
		     " enriched dofs for the crack tip");
      mf_ut_product.set_enrichment(enriched_dofs);
      mf_u3_product.set_enrichment(enriched_dofs);
      mf_theta_product.set_enrichment(enriched_dofs);
      mf_ut_sum.set_mesh_fems(mf_ut_product, mfls_ut);
      mf_u3_sum.set_mesh_fems(mf_u3_product, mfls_u3);
      mf_theta_sum.set_mesh_fems(mf_theta_product, mfls_theta);
    }
    break;
  default : {
  mf_ut_sum.set_mesh_fems(mfls_ut);
  mf_u3_sum.set_mesh_fems(mfls_u3);
  mf_theta_sum.set_mesh_fems(mfls_theta);}
   break;
  }

  
  cout << " Solving ---------------------- \n" ;
  
  
  getfem::mdbrick_abstract<> *ELAS, *SIMPLE, *LAST(0);

  // Linearized plate brick.
  getfem::mdbrick_isotropic_linearized_plate<>
    ELAS1(mim, mim, mf_ut(), mf_u3(), mf_theta(), lambda,
	  mu, epsilon);
  if (mitc) ELAS1.set_mitc();  
  
  ELAS = &ELAS1;

  // Defining the surface source term.
  plain_vector F(nb_dof_rhs * 3); 
  plain_vector M(nb_dof_rhs * 2);
  scalar_type r, theta, E, MapleGenVar1, MapleGenVar2, MapleGenVar3, MapleGenVar4, MapleGenVar5, MapleGenVar6, MapleGenVar7 ;
  E = - 4.*mu*(mu+lambda) / (2. * mu + lambda) ; 
  for (size_type i = 0; i < nb_dof_rhs; ++i){
  if (sol_ref == 0){
    if (mf_rhs.point_of_dof(i)[1] > 0 )
      F[3*i+2] = pressure;
    else F[3*i+2] = -pressure;
    }
  if (sol_ref == 1) {
    if (mf_rhs.point_of_dof(i)[0] <  1E-7 - 0.5 )
      M[2*i] =  - 2. * epsilon * epsilon  / (3. * gmm::sqrt(h_crack_length) ) ;
    }
  
  if (sol_ref == 3) {
     r = sqrt ( mf_rhs.point_of_dof(i)[0] * mf_rhs.point_of_dof(i)[0] + mf_rhs.point_of_dof(i)[1] * mf_rhs.point_of_dof(i)[1] ) ;
     theta = atan2( mf_rhs.point_of_dof(i)[1], mf_rhs.point_of_dof(i)[0] ) ;
     F[3*i+2] = - epsilon * E / (1.+nu) / sqrt(r) * ( 
                   (A1 + C1/2. + B2 + D2/2.) * cos(theta/2.)
                 + (B1 - D1/2. - A2 + C2/2.) * sin(theta/2.)
                 + (- A1/2. + B2/2. ) * cos(5. * theta/2.) 
                 + (- B1/2. - A2/2. ) * sin(5. * theta/2.)     ) ;  
      
      MapleGenVar1 = E ;
      MapleGenVar3 = epsilon*epsilon*epsilon;
      MapleGenVar5 = 1/(3.0+3.0*nu);
      MapleGenVar7 = 2.0/(1.0-nu)*(-((A1+C1/2.0)*cos(theta/2.0)+(B1-D1/2.0)*sin(theta/
2.0)-A1*cos(5.0/2.0*theta)/2.0-B1*sin(5.0/2.0*theta)/2.0)/sqrt(r*r*r)*cos(theta)/2.0-(-(A1+
C1/2.0)*sin(theta/2.0)/2.0+(B1-D1/2.0)*cos(theta/2.0)/2.0+5.0/4.0*A1*sin(5.0/2.0*theta)-5.0
/4.0*B1*cos(5.0/2.0*theta))/sqrt(r*r*r)*sin(theta))+(1.0+nu)/(1.0-nu)*(-((A2+C2/2.0)*
cos(theta/2.0)+(B2-D2/2.0)*sin(theta/2.0)-A2*cos(5.0/2.0*theta)/2.0-B2*sin(5.0/2.0*theta)/2.0)/
sqrt(r*r*r)*sin(theta)/2.0+(-(A2+C2/2.0)*sin(theta/2.0)/2.0+(B2-D2/2.0)*cos(theta/2.0)/2.0+
5.0/4.0*A2*sin(5.0/2.0*theta)-5.0/4.0*B2*cos(5.0/2.0*theta))/sqrt(r*r*r)*cos(theta));
      MapleGenVar6 = MapleGenVar7-((B1+D1/2.0)*cos(theta/2.0)+(-A1+C1/2.0)*sin(theta/
2.0)+B1*cos(5.0/2.0*theta)/2.0-A1*sin(5.0/2.0*theta)/2.0)/sqrt(r*r*r)*sin(theta)/2.0+(-(B1+
D1/2.0)*sin(theta/2.0)/2.0+(-A1+C1/2.0)*cos(theta/2.0)/2.0-5.0/4.0*B1*sin(5.0/2.0*theta)
-5.0/4.0*A1*cos(5.0/2.0*theta))/sqrt(r*r*r)*cos(theta)-3.0*(A0/sqrt(r)*sin(theta/2.0)*cos(theta)
/2.0-A0/sqrt(r)*cos(theta/2.0)*sin(theta)/2.0+sqrt(r)*(A1*cos(3.0/2.0*theta)+B1*sin(3.0/2.0*
theta)+C1*cos(theta/2.0)+D1*sin(theta/2.0)))*epsilon*epsilon;
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3*MapleGenVar4;
      M[2*i] = MapleGenVar1*MapleGenVar2;
      MapleGenVar1 = -4.*mu*(mu+lambda) / (2. * mu + lambda);
      MapleGenVar3 = epsilon*epsilon*epsilon;
      MapleGenVar5 = 1/(3.0+3.0*nu);
      MapleGenVar7 = 2.0/(1.0-nu)*(-((B2+D2/2.0)*cos(theta/2.0)+(-A2+C2/2.0)*sin(theta/
2.0)+B2*cos(5.0/2.0*theta)/2.0-A2*sin(5.0/2.0*theta)/2.0)/sqrt(r*r*r)*sin(theta)/2.0+(-(B2+
D2/2.0)*sin(theta/2.0)/2.0+(-A2+C2/2.0)*cos(theta/2.0)/2.0-5.0/4.0*B2*sin(5.0/2.0*theta)
-5.0/4.0*A2*cos(5.0/2.0*theta))/sqrt(r*r*r)*cos(theta))+(1.0+nu)/(1.0-nu)*(-((A1+C1/2.0
)*cos(theta/2.0)+(B1-D1/2.0)*sin(theta/2.0)-A1*cos(5.0/2.0*theta)/2.0-B1*sin(5.0/2.0*theta)/2.0
)/sqrt(r*r*r)*sin(theta)/2.0+(-(A1+C1/2.0)*sin(theta/2.0)/2.0+(B1-D1/2.0)*cos(theta/2.0)/
2.0+5.0/4.0*A1*sin(5.0/2.0*theta)-5.0/4.0*B1*cos(5.0/2.0*theta))/sqrt(r*r*r)*cos(theta));
      MapleGenVar6 = MapleGenVar7-((A2+C2/2.0)*cos(theta/2.0)+(B2-D2/2.0)*sin(theta/2.0
)-A2*cos(5.0/2.0*theta)/2.0-B2*sin(5.0/2.0*theta)/2.0)/sqrt(r*r*r)*cos(theta)/2.0-(-(A2+C2/
2.0)*sin(theta/2.0)/2.0+(B2-D2/2.0)*cos(theta/2.0)/2.0+5.0/4.0*A2*sin(5.0/2.0*theta)-5.0/
4.0*B2*cos(5.0/2.0*theta))/sqrt(r*r*r)*sin(theta)-3.0*(A0/sqrt(r)*sin(theta/2.0)*sin(theta)/2.0+
A0/sqrt(r)*cos(theta/2.0)*cos(theta)/2.0+sqrt(r)*(A2*cos(3.0/2.0*theta)+B2*sin(3.0/2.0*theta)+C2
*cos(theta/2.0)+D2*sin(theta/2.0)))*epsilon*epsilon;
      MapleGenVar4 = MapleGenVar5*MapleGenVar6;
      MapleGenVar2 = MapleGenVar3*MapleGenVar4;
      M[2*i+1 ] = MapleGenVar1*MapleGenVar2;
  } 
  } 
  getfem::mdbrick_plate_source_term<> VOL_F(*ELAS, mf_rhs, F, M);
  
  getfem::mdbrick_plate_clamped_support<> SIMPLE1
    (VOL_F, DIRICHLET_BOUNDARY_NUM, 0,getfem::AUGMENTED_CONSTRAINTS);
       
  SIMPLE = &SIMPLE1 ;

  

  cout << "Setting the Dirichlet condition brick : \n " ;
  
  getfem::mdbrick_Dirichlet<> DIRICHLET_U3(*SIMPLE, DIRICHLET_BOUNDARY_NUM, mf_mult_u3, 1);
  DIRICHLET_U3.rhs().set(exact_sol.mf_u3, exact_sol.U3);
  DIRICHLET_U3.set_constraints_type(getfem::constraints_type(dirichlet_version)); 
  cout << " md_brick DIRICHLET_U3 done     \n" ;
  
  getfem::mdbrick_Dirichlet<> DIRICHLET_THETA(DIRICHLET_U3, DIRICHLET_BOUNDARY_NUM, mf_mult_theta, 2);
  DIRICHLET_THETA.rhs().set(exact_sol.mf_theta, exact_sol.THETA);
  DIRICHLET_THETA.set_constraints_type(getfem::constraints_type(dirichlet_version));  
  cout << " md_brick DIRICHLET_THETA done     \n" ;
  
  getfem::mdbrick_Dirichlet<> DIRICHLET_UT(DIRICHLET_THETA, DIRICHLET_BOUNDARY_NUM, mf_mult_ut, 0);
  DIRICHLET_UT.rhs().set(exact_sol.mf_ut, exact_sol.UT);
  DIRICHLET_UT.set_constraints_type(getfem::constraints_type(dirichlet_version)); 
  cout << " md_brick DIRICHLET_UT done     \n" ;
  
  if (sol_ref == 0 || sol_ref == 1) LAST = &SIMPLE1 ;
  if (sol_ref == 3 )                LAST = &DIRICHLET_UT ;
  
  getfem::mdbrick_plate_closing<> final_model(*LAST, 0, 1);
  
  
  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);
  
  /*affichage de la solution */    
  gmm::resize(U3, mf_u3().nb_dof());
  gmm::copy(ELAS1.get_u3(MS), U3);  

  gmm::resize(UT, mf_ut().nb_dof());
  gmm::copy(ELAS1.get_ut(MS), UT);

  gmm::resize(THETA, mf_theta().nb_dof());
  gmm::copy(ELAS1.get_theta(MS), THETA);  
  cout << "vecteur solution u3 : \n" << U3 << "\n" ;
//   gmm::resize(U, mf_ut().nb_dof() + mf_u3().nb_dof() + mf_theta().nb_dof() );
//   gmm::copy(ELAS1.get_solution(MS), U);
  return (iter.converged());
}

/* compute the error with respect to the exact solution */
void crack_mindlin_problem::compute_error(plain_vector &, plain_vector &U3, plain_vector &THETA ) {
  cout << "Relative error on the vertical displacement u3 :\n" ;
  cout << "L2 ERROR:"
       << getfem::asm_L2_dist(mim, mf_u3(), U3, exact_sol.mf_u3, exact_sol.U3) 
        / getfem::asm_L2_norm(mim, exact_sol.mf_u3, exact_sol.U3)<< "\n";
  cout << "H1 ERROR:"
       << getfem::asm_H1_dist(mim, mf_u3(), U3, exact_sol.mf_u3, exact_sol.U3)
        / getfem::asm_H1_norm(mim, exact_sol.mf_u3, exact_sol.U3) << "\n";
  cout << "Relative error on the section rotations theta = (theta_1, theta_2 ) :\n" ;
  cout << "L2 ERROR:"
       << getfem::asm_L2_dist(mim, mf_theta(), THETA, exact_sol.mf_theta, exact_sol.THETA)
        / getfem::asm_L2_norm(mim, exact_sol.mf_theta, exact_sol.THETA) << "\n";
  cout << "H1 ERROR:"
       << getfem::asm_H1_dist(mim, mf_theta(), THETA, exact_sol.mf_theta, exact_sol.THETA)
        / getfem::asm_H1_norm(mim, exact_sol.mf_theta, exact_sol.THETA) << "\n";
}


/************************************************************
 * subroutine for evaluating Von Mises
 ************************************************************/

template <typename VEC1, typename VEC2>
void calcul_von_mises(const getfem::mesh_fem &mf_u, const VEC1 &U,
		      const getfem::mesh_fem &mf_vm, VEC2 &VM,
		      scalar_type mu=1) {
  /* DU=gf_compute(mfu,U,'gradient',mf_vm);
     
  // from the derivative, we compute the von mises stress
  VM=zeros(1,gf_mesh_fem_get(mf_vm,'nbdof'));
  N=gf_mesh_get(m,'dim');
  for i=1:size(DU,3),
  t=DU(:,:,i);
  E=(t+t')/2;
  VM(i) = sum(E(:).^2) - (1./N)*sum(diag(E))^2;
  end;
  VM = 4*pde.mu{1}^2*VM;
  */
  assert(mf_vm.get_qdim() == 1); 
  unsigned N = mf_u.linked_mesh().dim(); assert(N == mf_u.get_qdim());
  std::vector<scalar_type> DU(mf_vm.nb_dof() * N * N);

  getfem::compute_gradient(mf_u, mf_vm, U, DU);
  
  gmm::resize(VM, mf_vm.nb_dof());
  scalar_type vm_min, vm_max;
  for (size_type i=0; i < mf_vm.nb_dof(); ++i) {
    VM[i] = 0;
    scalar_type sdiag = 0.;
    for (unsigned j=0; j < N; ++j) {
      sdiag += DU[i*N*N + j*N + j];
      for (unsigned k=0; k < N; ++k) {
	scalar_type e = .5*(DU[i*N*N + j*N + k] + DU[i*N*N + k*N + j]);
	VM[i] += e*e;	
      }      
    }
    VM[i] -= 1./N * sdiag * sdiag;
    vm_min = (i == 0 ? VM[0] : std::min(vm_min, VM[i]));
    vm_max = (i == 0 ? VM[0] : std::max(vm_max, VM[i]));
    assert(VM[i] > -1e-6);
  }
  cout << "Von Mises : min=" << 4*mu*mu*vm_min << ", max=" << 4*mu*mu*vm_max << "\n";
  gmm::scale(VM, 4*mu*mu);
}
/************************************************************
 * main program
 ************************************************************/
 
 
int main(int argc, char *argv[]) {

  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.  

  try {
    crack_mindlin_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh") ;
    plain_vector UT, U3, THETA;
    if (!p.solve(UT, U3, THETA))
      GMM_ASSERT1(false, "Solve has failed");
    if (p.sol_ref == 3) p.compute_error(UT, U3, THETA) ;
          
    cout << "post-traitement pour l'affichage :\n" ;
    getfem::mesh mcut;
    p.mls.global_cut_mesh(mcut);

    getfem::stored_mesh_slice sl;
    getfem::mesh mcut_triangles_only;
    sl.build(mcut, 
	     getfem::slicer_build_mesh(mcut_triangles_only), 1);


    getfem::mesh_fem mf(mcut_triangles_only, 1);
    mf.set_classical_discontinuous_finite_element(2, 0.001);
//     mf.set_finite_element
//       	(getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 2, 0.0001)"));
//     getfem::pfem pmf = getfem::fem_descriptor
//         ("FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1, 2, 0.0001), FEM_PK_DISCONTINUOUS(1, 2, 0.0001))") ;
//     mf.set_finite_element(pmf) ; 
    plain_vector V3(mf.nb_dof()), VT(mf.nb_dof()*2), VTHETA(mf.nb_dof()*2);    
    getfem::interpolation(p.mf_ut(), mf, UT, VT);
    getfem::interpolation(p.mf_u3(), mf, U3, V3);
    getfem::interpolation(p.mf_theta(), mf, THETA, VTHETA);
    

    getfem::mesh m3d; getfem::extrude(mcut_triangles_only,m3d,1);
    getfem::base_matrix trans(3,3); 
    trans(0,0) = trans(1,1) = 1; trans(2,2) = p.epsilon;
    m3d.transformation(trans);
    getfem::mesh_fem mf3d(m3d);
    mf3d.set_classical_discontinuous_finite_element(2, 0.001);

    plain_vector V(mf3d.nb_dof()*3);
    bgeot::kdtree tree; tree.reserve(mf.nb_dof());
    for (unsigned i=0; i < mf.nb_dof(); ++i)
      tree.add_point_with_id(mf.point_of_dof(i),i);
    bgeot::kdtree_tab_type pts;
    for (unsigned i=0; i < mf3d.nb_dof(); ++i) {
      base_node P = mf3d.point_of_dof(i);
      base_node P2d0(2), P2d1(2); 
      scalar_type EPS = 1e-6;
      P2d0[0] = P[0]-EPS; P2d0[1] = P[1]-EPS;
      P2d1[0] = P[0]+EPS; P2d1[1] = P[1]+EPS;
      tree.points_in_box(pts, P2d0, P2d1);
      //cout << "P = " << P << ", P2d0=" << P2d0 << " " << P2d1 << ", pts.size=" << pts.size() << "\n";
      assert(pts.size() == 1);
      size_type j=pts[0].i;
      scalar_type x3 = P[2];
      assert(finite(VTHETA[2*j]));
      assert(finite(VT[2*j]));
      V[3*i+0] = VT[2*j+0] + x3 * VTHETA.at(2*j+0);
      V[3*i+1] = VT[2*j+1] + x3 * VTHETA[2*j+1];
      V[3*i+2] = V3[j];
      assert(finite(V[3*i]));
      assert(finite(V[3*i+1]));
      assert(finite(V[3*i+2]));
    }

    
//     getfem::stored_mesh_slice sl;
//     getfem::mesh mcut_refined;
//     sl.build(mcut, 
// 	getfem::slicer_build_mesh(mcut_refined), 4);
//     getfem::mesh_im mim_refined(mcut_refined); 
//     mim_refined.set_integration_method(getfem::int_method_descriptor
// 					("IM_TRIANGLE(6)"));
// 
//     getfem::mesh_fem mf_refined(mcut_refined, p.mf_u3().get_qdim());
//     mf_refined.set_finite_element
//     (getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 1, 0.0001)"));
//     plain_vector W(mf_refined.nb_dof());
//     getfem::interpolation(p.mf_u3(), mf_refined, U, W);

//     plain_vector EXACT(mf_refined.nb_dof());
//     p.exact_sol.mf.set_qdim(2);
//     assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());
//     getfem::interpolation(p.exact_sol.mf, mf_refined, 
// 			p.exact_sol.U, EXACT);

    if (p.PARAM.int_value("VTK_EXPORT"))
      {
	cout << "export to " << p.datafilename + ".vtk" << "..\n";
	getfem::vtk_export exp(p.datafilename + ".vtk",
			       p.PARAM.int_value("VTK_EXPORT")==1);
	exp.exporting(mf3d); 
	exp.write_point_data(mf3d, V, "plate_normal_displacement");
	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename << ".vtk -f ExtractVectorNorm -f "
	  "WarpVector -m BandedSurfaceMap -m Outline\n";
      
// 	getfem::vtk_export exp2("crack_exact.vtk");
// 	exp2.exporting(mf_refined);
// 	exp2.write_point_data(mf_refined, EXACT, "reference solution");
      }
      
    if (p.PARAM.int_value("DX_EXPORT"))
      {
	cout << "export to " << p.datafilename + ".dx" << "..\n";
	getfem::dx_export exp(p.datafilename + ".dx",
			       p.PARAM.int_value("DX_EXPORT")==1);

	/*
	exp.exporting(mf3d); 
	exp.write_point_data(mf3d, V, "plate_normal_displacement");
      	*/
	
	/* opendx ne supporte pas les prismes... */
	getfem::stored_mesh_slice sl_tetra;
	getfem::mesh mtetra;
	sl_tetra.build(m3d, getfem::slicer_build_mesh(mtetra), 1);
	getfem::mesh_fem mftetra(mtetra,3); 
	mftetra.set_classical_discontinuous_finite_element(3, 0.001);
	
	plain_vector Vtetra(mftetra.nb_dof());
	mf3d.set_qdim(3);
	getfem::interpolation(mf3d, mftetra, V, Vtetra);
	
	/* evaluating Von Mises*/
	getfem::mesh_fem mf_vm(mtetra,1) ; 
	mf_vm.set_classical_discontinuous_finite_element(3, 0.001);
	plain_vector VM(mf_vm.nb_dof()) ;
	calcul_von_mises(mftetra, Vtetra, mf_vm, VM, p.mu) ;

	
	exp.exporting(mftetra);
	exp.write_point_data(mftetra, Vtetra, "plate_displacement");
	size_type deb=1;
	if(deb==1){ 
	exp.exporting(mf_vm);
	exp.write_point_data(mf_vm, VM, "von_mises_stress");
	}
	
//     mf.set_finite_element



// 	getfem::dx_export exp2("crack_exact.vtk");
// 	exp2.exporting(mf_refined);
// 	exp2.write_point_data(mf_refined, EXACT, "reference solution");
      }

//       cout << "L2 ERROR:" << getfem::asm_L2_dist(p.mim, p.mf_u(), U, p.exact_sol.mf, p.exact_sol.U) << endl
// 	   << "H1 ERROR:" << getfem::asm_H1_dist(p.mim, p.mf_u(), U, p.exact_sol.mf, p.exact_sol.U)
// 	   << "\n";
//       
//       plain_vector DIFF(EXACT); gmm::add(gmm::scaled(W,-1),DIFF);
//       cout << "OLD ERROR L2:" << getfem::asm_L2_norm(mim_refined,mf_refined,DIFF) 
// 	   << " H1:" << getfem::asm_H1_dist(mim_refined,mf_refined,EXACT,mf_refined,W)
// 	   << "\n";
// 
//       cout << "ex = " << p.exact_sol.U << "\n";
//       cout << "U  = " << gmm::sub_vector(U, gmm::sub_interval(0,8)) << "\n";
    
    cout << "fin du programme atteinte \n" ;
  }
  
    GMM_STANDARD_CATCH_ERROR;

  return 0; 
}



















