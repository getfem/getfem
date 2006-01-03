// -*- c++ -*- (enables emacs c++ mode)
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2005 Yves Renard, Julien Pommier.                    */
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

/**
 * Linear Elastostatic problem with a crack.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_derivatives.h>
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <getfem_mesh_im_level_set.h>
#include <getfem_mesh_fem_level_set.h>
#include <getfem_mesh_fem_product.h>
#include <getfem_mesh_fem_global_function.h>
#include <getfem_spider_fem.h>
#include <getfem_mesh_fem_sum.h>
#include <gmm.h>

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

/**************************************************************************/
/*  Exact solution.                                                       */
/**************************************************************************/

#define VALIDATE_XFEM

#ifdef VALIDATE_XFEM

/* returns sin(theta/2) where theta is the angle
   of 0-(x,y) with the axis Ox */
scalar_type sint2(scalar_type x, scalar_type y) {
  scalar_type r = sqrt(x*x+y*y);
  if (r == 0) return 0;
  else return (y<0 ? -1:1) * sqrt(gmm::abs(r-x)/(2*r));
  // sometimes (gcc3.3.2 -O3), r-x < 0 ....
}
scalar_type cost2(scalar_type x, scalar_type y) {
  scalar_type r = sqrt(x*x+y*y);
  if (r == 0) return 0;
  else return sqrt(gmm::abs(r+x)/(2*r));
}
/* analytical solution for a semi-infinite crack [-inf,a] in an
   infinite plane submitted to +sigma above the crack
   and -sigma under the crack. (The crack is directed along the x axis).
   
   nu and E are the poisson ratio and young modulus
   
   solution taken from "an extended finite elt method with high order
   elts for curved cracks", Stazi, Budyn,Chessa, Belytschko
*/

void elasticite2lame(const scalar_type young_modulus,
		     const scalar_type poisson_ratio, 
		     scalar_type& lambda, scalar_type& mu) {
  mu = young_modulus/(2*(1+poisson_ratio));
  lambda = 2*mu*poisson_ratio/(1-poisson_ratio);
}

void sol_ref_infinite_plane(scalar_type nu, scalar_type E, scalar_type sigma,
			    scalar_type a, scalar_type xx, scalar_type y,
			    base_small_vector& U, int mode,
			    base_matrix *pgrad) {
  scalar_type x  = xx-a; /* the eq are given relatively to the crack tip */
  //scalar_type KI = sigma*sqrt(M_PI*a);
  scalar_type r = std::max(sqrt(x*x+y*y),1e-16);
  scalar_type sqrtr = sqrt(r), sqrtr3 = sqrtr*sqrtr*sqrtr;
  scalar_type cost = x/r, sint = y/r;
  scalar_type theta = atan2(y,x);
  scalar_type s2 = sin(theta/2); //sint2(x,y);
  scalar_type c2 = cos(theta/2); //cost2(x,y);
  // scalar_type c3 = cos(3*theta/2); //4*c2*c2*c2-3*c2; /* cos(3*theta/2) */
  // scalar_type s3 = sin(3*theta/2); //4*s2*c2*c2-s2;  /* sin(3*theta/2) */

  scalar_type lambda, mu;
  elasticite2lame(E,nu,lambda,mu);

  U.resize(2);
  if (pgrad) (*pgrad).resize(2,2);
  scalar_type C= 1./E * (mode == 1 ? 1. : (1+nu));
  if (mode == 1) {
    scalar_type A=2+2*mu/(lambda+2*mu);
    scalar_type B=-2*(lambda+mu)/(lambda+2*mu);
    U[0] = sqrtr/sqrt(2*M_PI) * C * c2 * (A + B*cost);
    U[1] = sqrtr/sqrt(2*M_PI) * C * s2 * (A + B*cost);
    if (pgrad) {
      (*pgrad)(0,0) = C/(2.*sqrt(2*M_PI)*sqrtr)
	* (cost*c2*A-cost*cost*c2*B+sint*s2*A+sint*s2*B*cost+2*c2*B);
      (*pgrad)(1,0) = -C/(2*sqrt(2*M_PI)*sqrtr)
	* (-sint*c2*A+sint*c2*B*cost+cost*s2*A+cost*cost*s2*B);
      (*pgrad)(0,1) = C/(2.*sqrt(2*M_PI)*sqrtr)
	* (cost*s2*A-cost*cost*s2*B-sint*c2*A-sint*c2*B*cost+2*s2*B);
      (*pgrad)(1,1) = C/(2.*sqrt(2*M_PI)*sqrtr)
	* (sint*s2*A-sint*s2*B*cost+cost*c2*A+cost*cost*c2*B);
    }
  } else if (mode == 2) {
    scalar_type C1 = (lambda+3*mu)/(lambda+mu);
    U[0] = sqrtr/sqrt(2*M_PI) * C * s2 * (C1 + 2 + cost);
    U[1] = sqrtr/sqrt(2*M_PI) * C * c2 * (C1 - 2 + cost) * (-1.);
    if (pgrad) {
      (*pgrad)(0,0) = C/(2.*sqrt(2*M_PI)*sqrtr)
	* (cost*s2*C1+2*cost*s2-cost*cost*s2-sint*c2*C1
	   -2*sint*c2-sint*cost*c2+2*s2);
      (*pgrad)(1,0) = C/(2.*sqrt(2*M_PI)*sqrtr)
	* (sint*s2*C1+2*sint*s2-sint*s2*cost+cost*c2*C1
	   +2*cost*c2+cost*cost*c2);
      (*pgrad)(0,1) = -C/(2.*sqrt(2*M_PI)*sqrtr)
	* (cost*c2*C1-2*cost*c2-cost*cost*c2+sint*s2*C1
	   -2*sint*s2+sint*s2*cost+2*c2);
      (*pgrad)(1,1) =  C/(2.*sqrt(2*M_PI)*sqrtr)
	* (-sint*c2*C1+2*sint*c2+sint*cost*c2+cost*s2*C1
	   -2*cost*s2+cost*cost*s2);
    }
  } else if (mode == 100) {
    U[0] = - sqrtr3 * (c2 + 4./3 *(7*mu+3*lambda)/(lambda+mu)*c2*s2*s2
		       -1./3*(7*mu+3*lambda)/(lambda+mu)*c2);
    U[1] = - sqrtr3 * (s2+4./3*(lambda+5*mu)/(lambda+mu)*s2*s2*s2
		       -(lambda+5*mu)/(lambda+mu)*s2);
    if (pgrad) {
      (*pgrad)(0,0) = 2*sqrtr*(-6*cost*c2*mu+7*cost*c2*c2*c2*mu
			       -3*cost*c2*lambda+3*cost*c2*c2*c2*lambda
			       -2*sint*s2*mu
			       +7*sint*s2*c2*c2*mu-sint*s2*lambda
			       +3*sint*s2*c2*c2*lambda)/(lambda+mu);
      (*pgrad)(1,0) = -2*sqrtr*(6*sint*c2*mu-7*sint*c2*c2*c2*mu
				+3*sint*c2*lambda-3*sint*c2*c2*c2*lambda
				-2*cost*s2*mu
				+7*cost*s2*c2*c2*mu-cost*s2*lambda
				+3*cost*s2*c2*c2*lambda)/(lambda+mu);
      (*pgrad)(0,1) = 2*sqrtr*(-2*cost*s2*mu-cost*s2*lambda
			       +cost*s2*c2*c2*lambda+5*cost*s2*c2*c2*mu
			       +4*sint*c2*mu
			       +sint*c2*lambda-sint*c2*c2*c2*lambda
			       -5*sint*c2*c2*c2*mu)/(lambda+mu);
      (*pgrad)(1,1) = 2*sqrtr*(-2*sint*s2*mu-sint*s2*lambda
			       +sint*s2*c2*c2*lambda+5*sint*s2*c2*c2*mu
			       -4*cost*c2*mu
			       -cost*c2*lambda+cost*c2*c2*c2*lambda
			       +5*cost*c2*c2*c2*mu)/(lambda+mu);
    }
  } else if (mode == 101) {
    U[0] = -4*sqrtr3*s2*(-lambda-2*mu+7*lambda*c2*c2
			 +11*mu*c2*c2)/(3*lambda-mu);
    U[1] = -4*sqrtr3*c2*(-3*lambda+3*lambda*c2*c2-mu*c2*c2)/(3*lambda-mu);
    if (pgrad) {
      (*pgrad)(0,0) = -6*sqrtr*(-cost*s2*lambda-2*cost*s2*mu
				+7*cost*s2*lambda*c2*c2
				+11*cost*s2*mu*c2*c2+5*sint*c2*lambda
				+8*sint*c2*mu-7*sint*c2*c2*c2*lambda
				-11*sint*c2*c2*c2*mu)/(3*lambda-mu);
      (*pgrad)(1,0) = -6*sqrtr*(-sint*s2*lambda-2*sint*s2*mu
				+7*sint*s2*lambda*c2*c2
				+11*sint*s2*mu*c2*c2-5*cost*c2*lambda
				-8*cost*c2*mu+7*cost*c2*c2*c2*lambda
				+11*cost*c2*c2*c2*mu)/(3*lambda-mu);
      (*pgrad)(0,1) = -6*sqrtr*(-3*cost*c2*lambda+3*cost*c2*c2*c2*lambda
				-cost*c2*c2*c2*mu-sint*s2*lambda
				+3*sint*s2*lambda*c2*c2
				-sint*s2*mu*c2*c2)/(3*lambda-mu);
      (*pgrad)(1,1) = 6*sqrtr*(3*sint*c2*lambda
			       -3*sint*c2*c2*c2*lambda+sint*c2*c2*c2*mu
			       -cost*s2*lambda+3*cost*s2*lambda*c2*c2
			       -cost*s2*mu*c2*c2)/(3*lambda-mu);
    }

  } else if (mode == 10166666) {

    U[0] = 4*sqrtr3*s2*(-lambda+lambda*c2*c2-3*mu*c2*c2)/(lambda-3*mu);
    U[1] = 4*sqrtr3*c2*(-3*lambda-6*mu+5*lambda*c2*c2+9*mu*c2*c2)/(lambda-3*mu);
    if (pgrad) {
      (*pgrad)(0,0) = 6*sqrtr*(-cost*s2*lambda+cost*s2*lambda*c2*c2-
			       3*cost*s2*mu*c2*c2-2*sint*c2*mu+sint*c2*lambda-
			       sint*c2*c2*c2*lambda
			       +3*sint*c2*c2*c2*mu)/(lambda-3*mu);
      (*pgrad)(1,0) = 6*sqrtr*(-sint*s2*lambda+sint*s2*lambda*c2*c2-
			       3*sint*s2*mu*c2*c2+2*cost*c2*mu-cost*c2*lambda+
			       cost*c2*c2*c2*lambda
			       -3*cost*c2*c2*c2*mu)/(lambda-3*mu);
      (*pgrad)(0,1) = 6*sqrtr*(-3*cost*c2*lambda-6*cost*c2*mu
			       +5*cost*c2*c2*c2*lambda+
			       9*cost*c2*c2*c2*mu-sint*s2*lambda-2*sint*s2*mu+
			       5*sint*s2*lambda*c2*c2
			       +9*sint*s2*mu*c2*c2)/(lambda-3*mu);
      (*pgrad)(1,1) = -6*sqrtr*(3*sint*c2*lambda+6*sint*c2*mu
				-5*sint*c2*c2*c2*lambda-
				9*sint*c2*c2*c2*mu-cost*s2*lambda-2*cost*s2*mu+
				5*cost*s2*lambda*c2*c2
				+9*cost*s2*mu*c2*c2)/(lambda-3*mu);
    }
  } else assert(0);
  if (isnan(U[0]))
    cerr << "raaah not a number ... nu=" << nu << ", E=" << E << ", sig="
	 << sigma << ", a=" << a << ", xx=" << xx << ", y=" << y << ", r="
	 << r << ", sqrtr=" << sqrtr << ", cost=" << cost << ", U=" << U[0]
	 << "," << U[1] << endl;
  assert(!isnan(U[0]));
  assert(!isnan(U[1]));
}

struct exact_solution {
  getfem::mesh_fem_global_function mf;
  getfem::base_vector U;

  exact_solution(getfem::mesh &me) : mf(me) {}
  
  void init(int mode, scalar_type lambda, scalar_type mu,
	    getfem::level_set &ls) {
    std::vector<getfem::pglobal_function> cfun(4);
    for (unsigned j=0; j < 4; ++j)
      cfun[j] = getfem::isotropic_crack_singular_2D(j, ls);
    mf.set_functions(cfun);
    
    mf.set_qdim(1);
    
    U.resize(8); assert(mf.nb_dof() == 4);
    getfem::base_vector::iterator it = U.begin();
    scalar_type coeff=0.;
    switch(mode) {
      case 1: {
	scalar_type A=2+2*mu/(lambda+2*mu), B=-2*(lambda+mu)/(lambda+2*mu);
	/* "colonne" 1: ux, colonne 2: uy */
	*it++ = 0;       *it++ = A-B; /* sin(theta/2) */
	*it++ = A+B;     *it++ = 0;   /* cos(theta/2) */
	*it++ = -B;      *it++ = 0;   /* sin(theta/2)*sin(theta) */ 
	*it++ = 0;       *it++ = B;   /* cos(theta/2)*cos(theta) */
	coeff = 1/sqrt(2*M_PI);
      } break;
      case 2: {
	scalar_type C1 = (lambda+3*mu)/(lambda+mu); 
	*it++ = C1+2-1;   *it++ = 0;
	*it++ = 0;      *it++ = -(C1-2+1);
	*it++ = 0;      *it++ = 1;
	*it++ = 1;      *it++ = 0;
	coeff = 2*(mu+lambda)/(lambda+2*mu)/sqrt(2*M_PI);
      } break;
      default:
	assert(0);
	break;
    }
    U *= coeff;
  }
};

base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N);
  return res;
}

#else

base_small_vector sol_f(const base_node &x) {
  int N = x.size();
  base_small_vector res(N); res[N-1] = x[N-1];
  return res;
}

#endif

/**************************************************************************/
/*  Structure for the crack problem.                                      */
/**************************************************************************/

struct crack_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_level_set mls;       /* the integration methods.              */
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u;
  getfem::mesh_fem mf_mult;
  getfem::mesh_fem_level_set mfls_u; 
  getfem::mesh_fem_global_function mf_sing_u;
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum;
  getfem::spider_fem *spider;
  getfem::mesh_fem mf_us;

  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  // getfem::mesh_fem& mf_u() { return mf_us; }
  

  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure for mixed form     */
#ifdef VALIDATE_XFEM
  exact_solution exact_sol;
#endif
  
  scalar_type lambda, mu;    /* Lamé coefficients.                           */

  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::level_set ls2, ls3;  /* The two level sets defining the additional crack.       */
  
  base_small_vector translation;
  scalar_type theta0;
  scalar_type spider_radius;
  unsigned spider_Nr;
  unsigned spider_Ntheta;
  int spider_K;
  scalar_type residue;       /* max residue for the iterative solvers        */
  bool mixed_pressure, add_crack;
  unsigned dir_with_mult;
  scalar_type cutoff_radius, enr_area_radius;
  int enrichment_option;
  
  std::string datafilename;
  ftool::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  crack_problem(void) : mls(mesh), mim(mls), mf_pre_u(mesh), mf_mult(mesh),
			mfls_u(mls, mf_pre_u), mf_sing_u(mesh),
			mf_partition_of_unity(mesh),
			mf_product(mf_partition_of_unity, mf_sing_u),
			mf_u_sum(mesh), mf_us(mesh), mf_rhs(mesh), mf_p(mesh),
#ifdef VALIDATE_XFEM
			exact_sol(mesh), 
#endif
			ls(mesh, 1, true), ls2(mesh, 1, true),
			ls3(mesh, 1, true) {}
};

/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void crack_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");

  add_crack = (PARAM.int_value("ADDITIONAL_CRACK", "An additional crack ?") != 0);
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
				      "Enrichment option");
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  spider_radius = PARAM.real_value("SPIDER_RADIUS","spider_radius");
  spider_Nr = PARAM.int_value("SPIDER_NR","Spider_Nr ");
  spider_Ntheta = PARAM.int_value("SPIDER_NTHETA","Ntheta ");
  spider_K = PARAM.int_value("SPIDER_K","K ");

  translation.resize(2); 
  translation[0] =0.5;
  translation[1] =0.;
  theta0 =0;
  
  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  base_small_vector tt(N); tt[1] = -0.5;
  mesh.translation(tt); 
  
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residue = PARAM.real_value("RESIDUE"); if (residue == 0.) residue = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  cutoff_radius = PARAM.real_value("CUTOFF", "Cutoff");
  mf_u().set_qdim(N);

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method simp_ppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
		getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
  
  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  if (add_crack) { mls.add_level_set(ls2); mls.add_level_set(ls3); }

  mim.set_simplex_im(simp_ppi, sing_ppi);
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_qdim(N);
  mf_partition_of_unity.set_classical_finite_element(1);
  
//   if (enrichment_option == 3 || enrichment_option == 4) {
//     spider = new getfem::spider_fem(spider_radius, mim, spider_Nr,
// 				    spider_Ntheta, spider_K, translation,
// 				    theta0);
//     mf_us.set_finite_element(mesh.convex_index(),spider->get_pfem());
//     for (dal::bv_visitor_c i(mf_us.convex_index()); !i.finished(); ++i) {
//       if (mf_us.fem_of_element(i)->nb_dof(i) == 0) {
// 	mf_us.set_finite_element(i,0);
//       }
//     }
//   }

  mixed_pressure =
    (PARAM.int_value("MIXED_PRESSURE","Mixed version or not.") != 0);
  dir_with_mult = PARAM.int_value("DIRICHLET_VERSINO");
  if (mixed_pressure) {
    std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
    mf_p.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEM_TYPE_P));
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
#ifdef VALIDATE_XFEM
    mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
#else
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if (un[0] - 1.0 < -1.0E-7) { // new Neumann face
      mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
    } else {
      cout << "normal = " << un << endl;
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
    }
#endif
  }
#ifdef VALIDATE_XFEM
  exact_sol.init(1, lambda, mu, ls);
#endif
}


base_small_vector ls_function(const base_node P, int num = 0) {
  scalar_type x = P[0], y = P[1];
  base_small_vector res(2);
  switch (num) {
    case 0: {
      res[0] = y;
      res[1] = -.5 + x;
    } break;
    case 1: {
      res[0] = gmm::vect_dist2(P, base_node(0.5, 0.)) - .25;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.27;
    } break;
    case 2: {
      res[0] = x - 0.25;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    default: assert(0);
  }
  return res;
}

bool crack_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  ls.reinit();  
  cout << "ls.get_mesh_fem().nb_dof() = " << ls.get_mesh_fem().nb_dof() << "\n";
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_dof(d), 0)[0];
    ls.values(1)[d] = ls_function(ls.get_mesh_fem().point_of_dof(d), 0)[1];
  }
  ls.touch();

  if (add_crack) {
    ls2.reinit();
    for (size_type d = 0; d < ls2.get_mesh_fem().nb_dof(); ++d) {
      ls2.values(0)[d] = ls_function(ls2.get_mesh_fem().point_of_dof(d), 1)[0];
      ls2.values(1)[d] = ls_function(ls2.get_mesh_fem().point_of_dof(d), 1)[1];
    }
    ls2.touch();
    
    ls3.reinit();
    for (size_type d = 0; d < ls3.get_mesh_fem().nb_dof(); ++d) {
      ls3.values(0)[d] = ls_function(ls2.get_mesh_fem().point_of_dof(d), 2)[0];
      ls3.values(1)[d] = ls_function(ls2.get_mesh_fem().point_of_dof(d), 2)[1];
    }
    ls3.touch();
  }

  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
  std::vector<getfem::pglobal_function> vfunc(4);
  for (size_type i = 0; i < 4; ++i)
    vfunc[i] = isotropic_crack_singular_2D(i, ls,
			 (enrichment_option == 2) ? 0.0 : cutoff_radius);
  
  mf_sing_u.set_functions(vfunc);


  if (enrichment_option == 3 || enrichment_option == 4) {
    spider = new getfem::spider_fem(spider_radius, mim, spider_Nr,
				    spider_Ntheta, spider_K, translation,
				    theta0);
    mf_us.set_finite_element(mesh.convex_index(),spider->get_pfem());
    for (dal::bv_visitor_c i(mf_us.convex_index()); !i.finished(); ++i) {
      if (mf_us.fem_of_element(i)->nb_dof(i) == 0) {
	mf_us.set_finite_element(i,0);
      }
    }
    spider->check();
  }



  switch (enrichment_option) {
  case 1 : mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u); break;
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
	DAL_WARNING0("There is " << enriched_dofs.card() <<
		     " enriched dofs for the crack tip");
      mf_product.set_enrichment(enriched_dofs);
      mf_u_sum.set_mesh_fems(mf_product, mfls_u);
    }
    break;
    case 3 : mf_u_sum.set_mesh_fems(mf_us); break;
  case 4 : 
    mf_u_sum.set_mesh_fems(mf_us, mfls_u); 
    break;
  default : mf_u_sum.set_mesh_fems(mfls_u); break;
  }

  U.resize(mf_u().nb_dof());

  if (mixed_pressure) cout << "Number of dof for P: " << mf_p.nb_dof() << endl;
  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<>
    ELAS(mim, mf_u(), mixed_pressure ? 0.0 : lambda, mu);

  getfem::mdbrick_abstract<> *pINCOMP;
  if (mixed_pressure) {
    getfem::mdbrick_linear_incomp<> *incomp
      = new getfem::mdbrick_linear_incomp<>(ELAS, mf_p);
    incomp->penalization_coeff().set(1.0/lambda);
    pINCOMP = incomp;
  } else pINCOMP = &ELAS;

  // Defining the volumic source term.
  plain_vector F(nb_dof_rhs * N);
  for (size_type i = 0; i < nb_dof_rhs; ++i)
      gmm::copy(sol_f(mf_rhs.point_of_dof(i)),
		gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
  
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(*pINCOMP, mf_rhs, F);

  // Defining the Neumann condition right hand side.
  gmm::clear(F);
 
  // Neumann condition brick.
  getfem::mdbrick_source_term<> NEUMANN(VOL_F, mf_rhs, F,NEUMANN_BOUNDARY_NUM);
  
  gmm::clear(F);

  //toto_solution toto(mf_rhs.linked_mesh()); toto.init();
  //assert(toto.mf.nb_dof() == 1);
  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> final_model(NEUMANN, DIRICHLET_BOUNDARY_NUM,
					  mf_mult);
#ifdef VALIDATE_XFEM
  final_model.rhs().set(exact_sol.mf, exact_sol.U);
#endif
  final_model.set_constraints_type(getfem::constraints_type(dir_with_mult));

  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residue, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);

  // Solution extraction
  gmm::copy(ELAS.get_solution(MS), U);

  return (iter.converged());
}
  
/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  DAL_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //getfem::getfem_mesh_level_set_noisy();


  try {
    crack_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");
    plain_vector U(p.mf_u().nb_dof());
    if (!p.solve(U)) DAL_THROW(dal::failure_error,"Solve has failed");

    {
      getfem::mesh mcut;
      p.mls.global_cut_mesh(mcut);
      unsigned Q = p.mf_u().get_qdim();
      getfem::mesh_fem mf(mcut, Q);
      mf.set_classical_discontinuous_finite_element(2, 0.001);
      // mf.set_finite_element
      //	(getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 2, 0.0001)"));
      plain_vector V(mf.nb_dof());

      getfem::interpolation(p.mf_u(), mf, U, V);

      getfem::stored_mesh_slice sl;
      getfem::mesh mcut_refined;

      unsigned NX = p.PARAM.int_value("NX"), nn;
      if (NX < 6) nn = 24;
      else if (NX < 12) nn = 8;
      else if (NX < 30) nn = 3;
      else nn = 1;

      /* choose an adequate slice refinement based on the distance to the crack tip */
      std::vector<bgeot::short_type> nrefine(mcut.convex_index().last_true()+1);
      for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
	scalar_type dmin=0, d;
	base_node Pmin,P;
	for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
	  P = mcut.points_of_convex(cv)[i];
	  d = gmm::vect_norm2(ls_function(P));
	  if (d < dmin || i == 0) { dmin = d; Pmin = P; }
	}

	if (dmin < 1e-5)
	  nrefine[cv] = nn*8;
	else if (dmin < .1) 
	  nrefine[cv] = nn*2;
	else nrefine[cv] = nn;
	if (dmin < .01)
	  cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin << " " << nrefine[cv] << "\n";
      }

      {
	getfem::mesh_slicer slicer(mcut); 
	getfem::slicer_build_mesh bmesh(mcut_refined);
	slicer.push_back_action(bmesh);
	slicer.exec(nrefine, getfem::mesh_region::all_convexes());
      }
      /*
      sl.build(mcut, 
      getfem::slicer_build_mesh(mcut_refined), nrefine);*/

      getfem::mesh_im mim_refined(mcut_refined); 
      mim_refined.set_integration_method(getfem::int_method_descriptor
					 ("IM_TRIANGLE(6)"));

      getfem::mesh_fem mf_refined(mcut_refined, Q);
      mf_refined.set_classical_discontinuous_finite_element(2, 0.0001);
      plain_vector W(mf_refined.nb_dof());

      getfem::interpolation(p.mf_u(), mf_refined, U, W);

#ifdef VALIDATE_XFEM
      p.exact_sol.mf.set_qdim(Q);
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);

      plain_vector DIFF(EXACT); gmm::add(gmm::scaled(W,-1),DIFF);
#endif

      if (p.PARAM.int_value("VTK_EXPORT")) {
	getfem::mesh_fem mf_refined_vm(mcut_refined, 1);
	mf_refined_vm.set_classical_discontinuous_finite_element(1, 0.0001);
	cerr << "mf_refined_vm.nb_dof=" << mf_refined_vm.nb_dof() << "\n";
	plain_vector VM(mf_refined_vm.nb_dof());

	cout << "computing von mises\n";
	getfem::interpolation_von_mises(mf_refined, mf_refined_vm, W, VM);

	plain_vector D(mf_refined_vm.nb_dof() * Q), 
	  DN(mf_refined_vm.nb_dof());
	
#ifdef VALIDATE_XFEM
	getfem::interpolation(mf_refined, mf_refined_vm, DIFF, D);
	for (unsigned i=0; i < DN.size(); ++i) {
	  DN[i] = gmm::vect_norm2(gmm::sub_vector(D, gmm::sub_interval(i*Q, Q)));
	}
#endif

	cout << "export to " << p.datafilename + ".vtk" << "..\n";
	getfem::vtk_export exp(p.datafilename + ".vtk",
			       p.PARAM.int_value("VTK_EXPORT")==1);
	exp.exporting(mf_refined); 
	//exp.write_point_data(mf_refined_vm, DN, "error");
	exp.write_point_data(mf_refined_vm, VM, "von mises stress");
	exp.write_point_data(mf_refined, W, "elastostatic_displacement");
      
#ifdef VALIDATE_XFEM

	plain_vector VM_EXACT(mf_refined_vm.nb_dof());


	/* getfem::mesh_fem_global_function mf(mcut_refined,Q);
	   std::vector<getfem::pglobal_function> cfun(4);
	   for (unsigned j=0; j < 4; ++j)
	   cfun[j] = getfem::isotropic_crack_singular_2D(j, p.ls);
	   mf.set_functions(cfun);
	   getfem::interpolation_von_mises(mf, mf_refined_vm, p.exact_sol.U,
	   VM_EXACT);
	*/


	getfem::interpolation_von_mises(mf_refined, mf_refined_vm, EXACT, VM_EXACT);
	getfem::vtk_export exp2("crack_exact.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined_vm, VM_EXACT, "exact von mises stress");
	exp2.write_point_data(mf_refined, EXACT, "reference solution");
	
#endif

	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename << ".vtk -f "
	  "WarpVector -m BandedSurfaceMap -m Outline\n";
      }


#ifdef VALIDATE_XFEM
      cout << "L2 ERROR:"<< getfem::asm_L2_dist(p.mim, p.mf_u(), U,
						p.exact_sol.mf, p.exact_sol.U)
	   << endl << "H1 ERROR:"
	   << getfem::asm_H1_dist(p.mim, p.mf_u(), U,
				  p.exact_sol.mf, p.exact_sol.U) << "\n";
      
      /* cout << "OLD ERROR L2:" 
	 << getfem::asm_L2_norm(mim_refined,mf_refined,DIFF) 
	 << " H1:" << getfem::asm_H1_dist(mim_refined,mf_refined,
	 EXACT,mf_refined,W)  << "\n";

	 cout << "ex = " << p.exact_sol.U << "\n";
	 cout << "U  = " << gmm::sub_vector(U, gmm::sub_interval(0,8)) << "\n";
      */
#endif
    }

  }
  DAL_STANDARD_CATCH_ERROR;

  return 0; 
}
