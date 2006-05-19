// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2006 Yves Renard, Julien Pommier.
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

#include <getfem_config.h>
#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_fourth_order.h>
#include <getfem_model_solvers.h>
#include <gmm.h>
#include <getfem_superlu.h>
#include <getfem_derivatives.h>
#include <gmm_inoutput.h>
#include <getfem_mesh_im_level_set.h>
#include <getfem_mesh_fem_level_set.h>
#include <getfem_mesh_fem_product.h>
#include <getfem_mesh_fem_global_function.h>
#include <getfem_mesh_fem_sum.h>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;


/******** Exact Solution *******************************/

scalar_type D  = 1.  ;
scalar_type nu = 0.3 ;
scalar_type A1 = 1., A2 = 1. ;
scalar_type b1 = 3. + (A2 / A1) * (24. * nu) / (3. * nu * nu - 6. * nu + 5. ) ;
scalar_type b2 = (3. * nu * nu - 6. * nu - 3. ) / (3. * nu * nu - 6. * nu + 5. ) ;

scalar_type sol_u(const base_node &x){
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 scalar_type theta = 2. * atan( x[1] / ( x[0] + r ) ) ;
 assert(gmm::abs(( x[0] + r )) > 1e-10);
 return pow(r, 1.5) * ( A1 * ( sin( 1.5 * theta ) - b1 * sin(0.5 * theta) )  
                      + A2 * ( cos( 1.5 * theta ) - b2 * cos(0.5 * theta) ); }
 
scalar_type sol_lapl_u(const base_node &x) {
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 assert(gmm::abs(( x[0] + r )) > 1e-10);
 scalar_type theta = 2. * atan( x[1] / ( x[0] + r ) ) ; 
return - 2. * ( A1 * b1 * sin( 0.5 * theta ) + A2 * b2 * cos(0.5 * theta) ) / sqrt(r) ; }

scalar_type sol_f(const base_node &x)
{ return 0. ; }

base_small_vector sol_du(const base_node &x) {
 base_small_vector res(x.size());
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 assert(gmm::abs(( x[0] + r )) > 1e-10);
 scalar_type theta = 2. * atan( x[1] / ( x[0] + r ) ) ; 
 res[0] = sqrt(r) * (1.5 * ( A1 * ( sin( 1.5 * theta ) - 3. * sin( 0.5 * theta) ) + A2 * ( cos(1.5 * theta) - cos(0.5 * theta) )) * cos(theta) - ( A1 * ( 1.5 * cos( 1.5 * theta) - 1.5 * cos( 0.5 * theta))+ A2 * ( - 1.5 * sin( 1.5 * theta) + 0.5 * sin( 0.5 * theta) ) ) * sin ( theta ) ) ;
 res[1] = sqrt(r) * (1.5 * ( A1 * ( sin( 1.5 * theta ) - 3. * sin( 0.5 * theta) ) + A2 * ( cos( 1.5 * theta) - cos( 0.5 * theta) ) ) * sin(theta) + ( A1 * ( 1.5 * cos( 1.5 * theta) - 1.5 * cos( 0.5 * theta) ) + A2 * ( - 1.5 * sin( 1.5 * theta ) + 0.5 * sin( 0.5 * theta ) ) ) * cos( theta ) ) ;
  return res;
}

base_small_vector neumann_val(const base_node &x)
{ base_small_vector res(x.size());
  res[0] = 0. ;
  res[1] = 0. ;
return res ; }

// base_matrix sol_hessian(const base_node &x) {
//   base_matrix m(x.size(), x.size());
//   // remplir la matrice hessienne de la solution exacte 
//   return m;
// }

// base_matrix sol_mtensor(const base_node &x) {
//   // moment de flexion de la solution exacte 
//   base_matrix m = sol_hessian(x), mm(x.size(), x.size());
//   scalar_type l = sol_lapl_u(x);
//   for (size_type i = 0; i < x.size(); ++i) mm(i,i) = l * nu;
//   gmm::scale(m, (1-nu));
//   gmm::add(mm, m);
//   gmm::scale(m, -D);
//   return m;
// }

// base_small_vector sol_bf(const base_node &x)
// { return -D * neumann_val(x); }







/******** Struct for the Bilaplacian problem *******************************/







struct bilaplacian_crack_problem {

  enum { CLAMPED_BOUNDARY_NUM = 0, SIMPLE_SUPPORT_BOUNDARY_NUM = 1,
	 FORCE_BOUNDARY_NUM = 2, MOMENTUM_BOUNDARY_NUM = 3};
  
  getfem::mesh mesh;        /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls;       
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u ; 
  getfem::mesh_fem_level_set mfls_u ; 
  getfem::mesh_fem_global_function mf_sing_u ;
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_u_product ;
  getfem::mesh_fem_sum mf_u_sum ;
  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  
  getfem::mesh_fem mf_mult; /* mesh_fem for the Dirichlet condition.        */
  getfem::mesh_fem mf_rhs;  /* mesh_fem for the right hand side (f(x),..)   */
  
  scalar_type residual;     /* max residual for the iterative solvers       */
  getfem::constraints_type dirichlet_version;
  
  scalar_type cutoff_radius, enr_area_radius;
  int enrichment_option;
  
  std::string datafilename;
  ftool::md_param PARAM;

  bool KL = 1 ;
  
  scalar_type epsilon ;      /* half-plate thickness */
    
  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);

  bilaplacian_crack_problem(void) : mls(mesh), mim(mls), mf_pre_u(mesh), 
			mfls_u(mls, mf_pre_ut), mf_sing_u(mesh)
			mf_partition_of_unity(mesh),
			mf_u_product(mf_partition_of_unity, mf_sing_u),
			mf_u_sum(mesh), mf_rhs(mesh), mf_mult(mesh),
			ls(mesh, 1, true) {} 
};




struct bilaplacian_singular_functions : public global_function, public context_dependencies {
  size_type l;             // singular function number
  const level_set &ls;
  mutable mesher_level_set mls0, mls1;
  mutable size_type cv;
  
  void update_mls(size_type cv_) const { 
    if (cv_ != cv) 
      { cv=cv_; mls0=ls.mls_of_convex(cv, 0); mls1=ls.mls_of_convex(cv, 1); }
  }


  scalar_type sing_function(scalar_type x, scalar_type y) {
    scalar_type r = sqrt(x*x + y*y);
    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));

    scalar_type theta = atan2(y,x); 
 
    switch (l) {
    case 0: {
      return r*sqrt(r)*sin(theta/2);
    } break;
    default: assert(0); 
    }
  }

  /* derivative in the levelset coordinates */
  void sing_function_grad(scalar_type x, scalar_type y, base_small_vector &g) {
    g.resize(2);
    scalar_type r = sqrt(x*x + y*y);
    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));

    scalar_type theta = atan2(y,x);

    switch (l) {
    case 0: {
      g[0] = 3.0/2.0*sqrt(r)*sin(theta/2.0)*cos(theta)-sqrt(r)*cos(theta/2.0)*sin(theta)/2.0;
      g[1] = 3.0/2.0*sqrt(r)*sin(theta/2.0)*sin(theta)+sqrt(r)*cos(theta/2.0)*cos(theta)/2.0;
    } break;
    default: assert(0); 
    }
  }
   
  void sing_function_hess(scalar_type x, scalar_type y, base_matrix &he) {
    he.resize(2,2);
    scalar_type r = sqrt(x*x + y*y);
    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));
    
    scalar_type theta = atan2(y,x);

    switch (l) {
    case 0: {
      he(0,0) = (3.0/4.0/sqrt(r)*sin(theta/2.0)*cos(theta)-1/sqrt(r)*cos(theta/2.0)*
		 sin(theta)/4.0)*cos(theta)-(sqrt(r)*cos(theta/2.0)*cos(theta)/4.0-5.0/4.0*
					     sqrt(r)*sin(theta/2.0)*sin(theta))*sin(theta)/r;
      he(1,0) = (3.0/4.0/sqrt(r)*sin(theta/2.0)*cos(theta)-1/sqrt(r)*cos(theta/2.0)*
		 sin(theta)/4.0)*sin(theta)+(sqrt(r)*cos(theta/2.0)*cos(theta)/4.0-5.0/4.0*
					     sqrt(r)*sin(theta/2.0)*sin(theta))*cos(theta)/r;
      he(0,1) = he(1,0);
      he(1,1) = (3.0/4.0/sqrt(r)*sin(theta/2.0)*sin(theta)+1/sqrt(r)*cos(theta/2.0)*
		 cos(theta)/4.0)*sin(theta)+(sqrt(r)*cos(theta/2.0)*sin(theta)/4.0+5.0/4.0*
					     sqrt(r)*sin(theta/2.0)*cos(theta))*cos(theta)/r;
    } break;
    default: assert(0); 
    }
  }
   
  
  virtual scalar_type val(const fem_interpolation_context& c) const {
    update_mls(c.convex_num());
    scalar_type x = mls1(c.xref()), y = mls0(c.xref());
    scalar_type v = sing_function(x, y);
    return v;
  }

  virtual void grad(const fem_interpolation_context& c,
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
  
  virtual void hess(const fem_interpolation_context& c, base_matrix &he) const {
    update_mls(c.convex_num());
    size_type P = c.xref().size();
    base_small_vector dx(P), dy(P), dfr(2);  
    scalar_type x = mls1.grad(c.xref(), dx), y = mls0.grad(c.xref(), dy);
    base_matrix d2x(P,P), d2y(P,P), hefr(P, P) ;
    mls1.hess(c.xref(), d2x) ;
    mls0.hess(c.xref(), d2y) ;
    if (x*x + y*y < 1e-20) {
      cerr << "Warning, point very close to the singularity. xreal = "
	   << c.xreal() << ", x_crack = " << x << ", y_crack=" << y << "\n";
    } 

    sing_function_grad(x, y, dfr); // gradient in the coord of the levelset
    sing_function_hess(x, y, hefr); // hessian in the coord of the levelset
    
    base_small_vector dref = dfr[0]*dx + dfr[1]*dy; // gradient in the coordinate the ref convex.

    base_vector hh(4);
    // expression of the derivatives in the element reference basis
    hh[0] = hefr(0,0) * gmm::sqr(dx[0]) + hefr(0,1) * (2. * dy[0] * dx[0] )
      + hefr(1,1) * gmm::sqr(dy[0]) + dfr[0] * d2x(0,0) + dfr[1] * d2y(0,0) ;
    hh[3] = hefr(0,0) * gmm::sqr(dx[1]) + hefr(0,1) * (2. * dy[1] * dx[1] )
      + hefr(1,1) * gmm::sqr(dy[1]) + dfr[0] * d2x(1,1) + dfr[1] * d2y(1,1) ;  
    hh[1] = hh[2] = hefr(0, 0) * dx[1] * dx[0] + hefr(0, 1) * ( dy[1] * dx[0] + dx[1] * dy[0])
      + hefr(1, 1) * dy[0] * dy[1] + dfr[0] * d2x(0,1) + dfr[1] * d2y(0, 1) ;

    base_vector &vhe = he;
    // back in the real element basis
    gmm::mult_add(c.B2(), gmm::scaled(dref,-1), hh);
    gmm::mult_add(c.B3(), hh, vhe);
  }
    
  void update_from_context(void) const { cv =  size_type(-1); }

  bilaplacian_singular_functions(size_type l_, const level_set &ls_) : l(l_), ls(ls_) {
    cv = size_type(-1);
    this->add_dependency(ls);
  }

};


/******************* Methods ********************************/





void bilaplacian_crack_problem::init(void) {
  std::string MESH_FILE = PARAM.string_value("MESH_FILE");
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  size_type N = 2 ;

    cout << "MESH_TYPE=" << MESH_TYPE << "\n";
    
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
				      "Enrichment option");
    
    
    /* First step : build the mesh */
    
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(MESH_TYPE);
    N = pgt->dim();
    if (N != 2)
     DAL_THROW(getfem::failure_error, "For a plate problem, N should be 2");
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Number of space steps "));
    getfem::regular_unit_mesh(mesh, nsubdiv, pgt, PARAM.int_value("MESH_NOISED") != 0);

    bgeot::base_matrix M(N,N);
    for (size_type i=0; i < N; ++i) {
      static const char *t[] = {"LX","LY","LZ"};
      M(i,i) = (i<3) ? PARAM.real_value(t[i],t[i]) : 1.0;
    }
    /* scale the unit mesh to [LX,LY,..] and incline it */
    mesh.transformation(M);
    
    base_small_vector tt(N); 
    tt[0] = -0.5 ;
    tt[1] = -0.5;
    mesh.translation(tt); 
  
    
   /* read the parameters   */
  epsilon = PARAM.real_value("EPSILON", "thickness") ;
  int dv = PARAM.int_value("DIRICHLET_VERSION", "Dirichlet version");
  dirichlet_version = getfem::constraints_type(dv);
  datafilename=PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual=PARAM.real_value("RESIDUAL"); if (residual == 0.) residual = 1e-10;
  FT = PARAM.real_value("FT"); if (FT == 0.0) FT = 1.0;
  KL = (PARAM.int_value("KL", "Kirchhoff-Love model or not") != 0);
  D = PARAM.real_value("D", "Flexion modulus");
  if (KL) nu = PARAM.real_value("NU", "Poisson ratio");

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method sppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
    
  mim.set_integration_method(mesh.convex_index(), ppi);
  mf_u.set_finite_element(mesh.convex_index(), pf_u);
  mls.add_level_set(ls);
  mim.set_simplex_im(sppi);
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_partition_of_unity.set_classical_finite_element(1);

  
  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  if (dirichlet_fem_name.size() == 0)
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
  }

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    if (!pf_u->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    mesh.region(SIMPLE_SUPPORT_BOUNDARY_NUM).add(i.cv(), i.f());
    mesh.region(CLAMPED_BOUNDARY_NUM).add(i.cv(), i.f());
    
  }
}

/* compute the error with respect to the exact solution */
void bilaplacian_problem::compute_error(plain_vector &U) {
  std::vector<scalar_type> V(mf_rhs.nb_dof());
  getfem::interpolation(mf_u, mf_rhs, U, V);
  for (size_type i = 0; i < mf_rhs.nb_dof(); ++i)
    V[i] -= sol_u(mf_rhs.point_of_dof(i));
  cout.precision(16);
  cout << "L2 error = " << getfem::asm_L2_norm(mim, mf_rhs, V) << endl
       << "H1 error = " << getfem::asm_H1_norm(mim, mf_rhs, V) << endl
       << "H2 error = " << getfem::asm_H2_norm(mim, mf_rhs, V) << endl
       << "Linfty error = " << gmm::vect_norminf(V) << endl;
}

/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool bilaplacian_crack_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  cout << "Number of dof for u: " << mf_u.nb_dof() << endl;
  
  ls.reinit();  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    ls.values(0)[d] = -0.5 + (ls.get_mesh_fem().point_of_dof(d))[1];
    ls.values(1)[d] = -0.5 + (ls.get_mesh_fem().point_of_dof(d))[0];}
  ls.touch();
  
  mls.adapt();
  mim.adapt();
  mfls_ut.adapt();
  mfls_u3.adapt();
  mfls_theta.adapt();
  
  // setting singularities 
  std::vector<getfem::pglobal_function> utfunc(1) ;                  
  for (size_type i = 0 ; i < 1 ; ++i) {                              
    utfunc[i] = new bilaplacian_singular_functions(i, ls);
  }
  
  mf_sing_u.set_functions(ufunc);
  
 
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
    DAL_WARNING0("There is " << enriched_dofs.card() <<
		 " enriched dofs for the crack tip");
  mf_u_product.set_enrichment(enriched_dofs);
  mf_u_sum.set_mesh_fems(mf_u_product, mfls_u);
  
  
  // Bilaplacian brick.
  getfem::mdbrick_bilaplacian<> BIL(mim, mf_u);
  BIL.D().set(D);
  if (KL) { BIL.set_to_KL(); BIL.nu().set(nu); }
  
  // Defining the normal derivative Dirichlet condition value.
  plain_vector F(nb_dof_rhs*N);
  getfem::interpolation_function(mf_rhs, F, sol_du, CLAMPED_BOUNDARY_NUM);
  
   
  // Normal derivative Dirichlet condition brick.
  
 
  getfem::mdbrick_normal_derivative_Dirichlet<>                   
    NDER_DIRICHLET(BIL, CLAMPED_BOUNDARY_NUM, mf_mult);       
 
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.rhs().set(mf_rhs, F);

  // Defining the Dirichlet condition value.
  gmm::resize(F, nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_u,SIMPLE_SUPPORT_BOUNDARY_NUM);

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<>
    final_model(NDER_DIRICHLET, SIMPLE_SUPPORT_BOUNDARY_NUM, mf_mult);
  final_model.set_constraints_type(dirichlet_version);
  final_model.rhs().set(mf_rhs, F);

  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);

  // Solution extraction
  gmm::copy(BIL.get_solution(MS), U);
  return (iter.converged());

}

int main(int argc, char *argv[]) {

try {
    bilaplacian_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    plain_vector U(p.mf_u.nb_dof());
    if (!p.solve(U)) DAL_THROW(dal::failure_error, "Solve has failed");

    p.compute_error(U);
  
  return 0; 
}
























