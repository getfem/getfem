// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2006-2007 Yves Renard, Julien Pommier.
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

#include "getfem/getfem_config.h"
#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_fourth_order.h"
#include "getfem/getfem_model_solvers.h"
#include "gmm/gmm.h"
#include "getfem/getfem_superlu.h"
#include "getfem/getfem_derivatives.h"
#include "gmm/gmm_inoutput.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_mesh_fem_sum.h"

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_vector; /* dense vector. */
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* dense matrix. */

/* definition of some matrix/vector types. 
 * default types of getfem_model_solvers.h
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;


/******** Exact Solution *******************************/

scalar_type D  = 1.  ;
scalar_type nu = 0.3 ;
scalar_type BB = 0.1 ;
scalar_type AA = BB * (3. * nu + 5.)/ (3. * (nu - 1.))   ;  // (-3.0+nu*nu-2.0*nu)/(nu*nu-2.0*nu+5.0);
scalar_type DD = 0.0 ;
scalar_type CC = DD * (nu + 7.)/ (3. * (nu - 1.))   ;   //  (-8.0*nu+3.0*BB*nu*nu-6.0*nu*BB+15.0*BB)/(nu*nu-2.0*nu+5.0);
 

scalar_type sol_u(const base_node &x){
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 //scalar_type theta = 2. * atan( x[1] / ( x[0] + r ) ) ;
 scalar_type theta = atan2(x[1], x[0]);
 //return sqrt(r*r*r) * (sin(3.0/2.0*theta)+AA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0));
 return sqrt(r*r*r)*(AA*sin(3.0/2.0*theta)+BB*sin(theta/2.0)+CC*cos(3.0/2.0*theta)+DD*cos(theta/2.0));
}
 
scalar_type sol_lapl_u(const base_node &x) {
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 scalar_type theta = atan2(x[1], x[0]);
 return 2.0*(BB*sin(theta/2.0)+DD*cos(theta/2.0))/sqrt(r);
 /* return 9.0/4.0/sqrt(r)*(sin(3.0/2.0*theta)+AA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0))+1/sqrt(r)*(-9.0/4.0*sin(3.0/2.0*theta)-AA*sin(theta/
    2.0)/4.0-9.0/4.0*BB*cos(3.0/2.0*theta)-CC*cos(theta/2.0)/4.0); */ }

scalar_type sol_f(const base_node &)
{ return 0. ; }

base_small_vector sol_du(const base_node &x) {
 base_small_vector res(x.size());
 scalar_type r = sqrt( x[0] * x[0] + x[1] * x[1] ) ;
 scalar_type theta = atan2(x[1], x[0]);
res[0] = 3.0/2.0*sqrt(r)*(AA*sin(3.0/2.0*theta)+BB*sin(theta/2.0)+CC*cos(3.0/2.0*theta)
+DD*cos(theta/2.0))*cos(theta)-sqrt(r)*(3.0/2.0*AA*cos(3.0/2.0*theta)+
BB*cos(theta/2.0)/2.0-3.0/2.0*CC*sin(3.0/2.0*theta)-DD*sin(theta/2.0)/2.0)*sin(theta );

res[1] = 3.0/2.0*sqrt(r)*(AA*sin(3.0/2.0*theta)+BB*sin(theta/2.0)+CC*cos(3.0/2.0*theta)
+DD*cos(theta/2.0))*sin(theta)+sqrt(r)*(3.0/2.0*AA*cos(3.0/2.0*theta)+
BB*cos(theta/2.0)/2.0-3.0/2.0*CC*sin(3.0/2.0*theta)-DD*sin(theta/2.0)/2.0)*cos(theta);

/*
res[0] =  3.0/2.0*sqrt(r)*(sin(3.0/2.0*theta)+AA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0))*cos(theta)-sqrt(r)*(3.0/2.0*cos(3.0/2.0*theta)+AA*
cos(theta/2.0)/2.0-3.0/2.0*BB*sin(3.0/2.0*theta)-CC*sin(theta/2.0)/2.0)*sin(
theta);
 
res[1] = 3.0/2.0*sqrt(r)*(sin(3.0/2.0*theta)+AA*sin(theta/2.0)+BB*cos(3.0/2.0*theta)+CC*cos(theta/2.0))*sin(theta)+sqrt(r)*(3.0/2.0*cos(3.0/2.0*theta)+AA*
cos(theta/2.0)/2.0-3.0/2.0*BB*sin(3.0/2.0*theta)-CC*sin(theta/2.0)/2.0)*cos(
theta);
*/
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



scalar_type eval_fem_gradient_with_finite_differences(getfem::pfem pf, 
					       const base_vector &coeff,
					       size_type cv,
					       bgeot::pgeometric_trans pgt, 
					       bgeot::geotrans_inv_convex &gic,
					       const base_matrix &G, 
					       base_node X0, 
					       scalar_type h, unsigned dg) {
  X0[dg] -= h/2;
  base_node X0ref; gic.invert(X0, X0ref);
  getfem::fem_interpolation_context c(pgt, pf, X0ref, G, cv);

  base_vector val0(1);
  pf->interpolation(c, coeff, val0, 1);

  base_node X1(X0), X1ref; X1[dg] += h;
  gic.invert(X1, X1ref);
  c.set_xref(X1ref);

  base_vector val1(1);
  pf->interpolation(c, coeff, val1, 1);

  return (val1[0] - val0[0])/h;
}

scalar_type eval_fem_hessian_with_finite_differences(getfem::pfem pf, 
					      const base_vector &coeff,
					      size_type cv, 
					      bgeot::pgeometric_trans pgt, 
					      bgeot::geotrans_inv_convex &gic,
					      const base_matrix &G, 
					      base_node X0, 
					      scalar_type h, 
					      unsigned dg, unsigned dh) {
  X0[dh] -= h/2;
  scalar_type Gr0 = 
    eval_fem_gradient_with_finite_differences(pf, coeff, cv, pgt, gic, G, X0, h, dg);
  base_node X1(X0);
  X1[dh] += h;
  scalar_type Gr1 = 
    eval_fem_gradient_with_finite_differences(pf, coeff, cv, pgt, gic, G, X1, h, dg);
  return (Gr1 - Gr0)/h;
}

void validate_fem_derivatives(getfem::pfem pf, unsigned cv, 
			      bgeot::pgeometric_trans pgt, const base_matrix &G) {
  unsigned N = gmm::mat_nrows(G);
  scalar_type h = 1e-5;

  std::vector<base_node> pts(gmm::mat_ncols(G));
  for (unsigned j=0; j < pts.size(); ++j) {
    pts[j].resize(N); gmm::copy(gmm::mat_col(G, j), pts[j]);
  }
  cout << "validate_fem_derivatives: pf = " << &(*pf) << ", nbdof = "<< pf->nb_dof(cv) << ", cv = " << cv << " (~ at " << gmm::mean_value(pts) << ")\n";
  bgeot::geotrans_inv_convex gic(pts, pgt);

  //cout << "pts = " << pts << "\n";
  
  for (unsigned idof = 0; idof < pf->nb_dof(cv); ++idof) {    
    /* choose a random point in the convex */
    base_node X0(N), X0ref;
    base_node w(pgt->nb_points());
    do {
      for (unsigned i=0; i < w.size(); ++i) w[i] = 0.1 + 0.8*gmm::random(); 
      gmm::scale(w, 1/gmm::vect_norm1(w));
      gmm::mult(G, w, X0);

      //cout << "w = " << w << "\n";
      
      gic.invert(X0, X0ref);
      
      // avoid discontinuity lines in the HCT composite element..
      if (gmm::abs(X0ref[0] + X0ref[1] - 1) > 1e-2 &&
	  gmm::abs(X0ref[0] - X0ref[1]) > 1e-2 && 
	  gmm::abs(X0[0]) > 1e-3 && gmm::abs(X0[1])> 1e-3) break;
    } while (1);
    //cout << "testing X0 = " << X0 << " (X0ref=" << X0ref << ")\n";


    base_vector coeff(pf->nb_dof(cv)); coeff[idof] = 1;
    base_matrix grad(1,N), grad_fd(1,N);
    base_matrix hess(1,N*N), hess_fd(1,N*N);

    getfem::fem_interpolation_context c(pgt, pf, X0ref, G, cv);
    pf->interpolation_grad(c, coeff, grad, 1);
    pf->interpolation_hess(c, coeff, hess, 1);

    for (unsigned dg = 0; dg < N; ++dg) {
      grad_fd[dg] = 
	eval_fem_gradient_with_finite_differences(pf, coeff, cv, pgt, gic, G, X0, h, dg);
      for (unsigned dh = 0; dh < N; ++dh) {
	hess_fd(0,dg*N+dh) = 
	  eval_fem_hessian_with_finite_differences(pf, coeff, cv, pgt, gic, G, X0, h, dg, dh);
      }
    }
    
    scalar_type err_grad = 
      gmm::vect_dist2((base_vector&)grad, (base_vector&)grad_fd);
    scalar_type err_hess = 
      gmm::vect_dist2((base_vector&)hess, (base_vector&)hess_fd);
    
    if (err_grad > 1e-4 ||
	err_hess > 1e-4) {
      cout << "validate_fem_derivatives dof=" << idof << "/" << pf->nb_dof(cv) << " -- X0ref = " << X0ref << "\n";

      if (gmm::vect_dist2((base_vector&)grad, (base_vector&)grad_fd) > 1e-4)
	cout << "grad = " << (base_vector&)grad << "\ngrad_fd = " << (base_vector&)grad_fd << "\n";
      cout << "hess = " << (base_vector&)hess << "\nhess_fd = " << (base_vector&)hess_fd << "\n";
      if (err_grad + err_hess > 1.0) { cout << "---------> COMPLETEMENT FAUX!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"; abort(); }
    }
  }
}


void validate_fem_derivatives(const getfem::mesh_fem &mf) {
  bgeot::base_matrix G;
  for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
    //if (mf.nb_dof_of_element(cv) > 12) {
      vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));
      validate_fem_derivatives(mf.fem_of_element(cv), cv, mf.linked_mesh().trans_of_convex(cv), G);
      //}
  }
}








struct bilaplacian_singular_functions : public getfem::global_function, public getfem::context_dependencies {
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
    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    /*scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));
    */

    scalar_type theta = atan2(y,x); 
 
    switch (l) {
    case 0: {
      return r*sqrt(r)*sin(theta/2);
    } break;
    case 1: {
      return r*sqrt(r)*sin(3.0 * theta/2);
    } break;
    case 2: {
      return r*sqrt(r)*cos(3.0 * theta/2);
    } break;
    case 3: {
      return r*sqrt(r)*cos(theta/2);
    } break;
    default: assert(0); 
    }
    return 0;
  }

  /* derivative in the levelset coordinates */
  void sing_function_grad(scalar_type x, scalar_type y, base_small_vector &g) const {
    g.resize(2);
    scalar_type r = sqrt(x*x + y*y);
    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    /*scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));
    */
    scalar_type theta = atan2(y,x);

    switch (l) {
    case 0: {
      g[0] = 3.0/2.0*sqrt(r)*sin(theta/2.0)*cos(theta)-sqrt(r)*cos(theta/2.0)*sin(theta)/2.0;
      g[1] = 3.0/2.0*sqrt(r)*sin(theta/2.0)*sin(theta)+sqrt(r)*cos(theta/2.0)*cos(theta)/2.0;
    } break;
    case 1: {
      g[0] = 3.0/2.0*sqrt(r)* ( sin(3.0/2.0*theta)*cos(theta)-cos(3.0/2.0*theta)*sin(theta) );
      g[1] = 3.0/2.0*sqrt(r)* ( sin(3.0/2.0*theta)*sin(theta)+cos(3.0/2.0*theta)*cos(theta) );
    } break;
    case 2: {
      g[0] = 3.0/2.0*sqrt(r)*sin(3.0/2.0*theta)*sin(theta)+3.0/2.0*sqrt(r)*cos(3.0/2.0*theta)*cos(theta);
      g[1] = 3.0/2.0*sqrt(r)*cos(3.0/2.0*theta)*sin(theta)-3.0/2.0*sqrt(r)*sin(3.0/2.0*theta)*cos(theta);
    } break;
    case 3: {
      g[0] = 3.0/2.0*sqrt(r)*cos(theta/2.0)*cos(theta)+sqrt(r)*sin(theta/2.0)*sin(theta)/2.0;
      g[1] = 3.0/2.0*sqrt(r)*cos(theta/2.0)*sin(theta)-sqrt(r)*sin(theta/2.0)*cos(theta)/2.0;
    } break;
    default: assert(0); 
    }
  }
   
  void sing_function_hess(scalar_type x, scalar_type y, base_matrix &he) const {
    he.resize(2,2);
    scalar_type r = sqrt(x*x + y*y);
    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    /*scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));
    */
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
    case 1: {
      he(0,0) = (3.0/4.0/sqrt(r)*sin(3.0/2.0*theta)*cos(theta)-3.0/4.0/sqrt(r)
                 *cos(3.0/2.0*theta)*sin(theta))*cos(theta)-(3.0/4.0*sqrt(r)*
                  cos(3.0/2.0*theta)*cos(theta)+3.0/4.0*sqrt(r)*sin(3.0/2.0*theta)*sin(theta))*sin(theta)/r;
      he(1,0) = (3.0/4.0/sqrt(r)*sin(3.0/2.0*theta)*cos(theta)-3.0/4.0/sqrt(r)
		 *cos(3.0/2.0*theta)*sin(theta))*sin(theta)+(3.0/4.0*sqrt(r)*
		 cos(3.0/2.0*theta)*cos(theta)+3.0/4.0*sqrt(r)*sin(3.0/2.0*theta)*sin(theta))*cos(theta)/r;
      he(1,1) = (3.0/4.0/sqrt(r)*sin(3.0/2.0*theta)*sin(theta)+3.0/4.0/sqrt(r)
		 *cos(3.0/2.0*theta)*cos(theta))*sin(theta)+(3.0/4.0*sqrt(r)*cos(3.0/2.0*theta)
		 *sin(theta)-3.0/4.0*sqrt(r)*sin(3.0/2.0*theta)*cos(theta))*cos(theta)/r;
      he(0,1) = he(1,0) ;
    } break;
    case 2: {
      he(0,0) = (3.0/4.0/sqrt(r)*sin(3.0/2.0*theta)*sin(theta)+3.0/4.0/sqrt(r)*
	        cos(3.0/2.0*theta)*cos(theta))*cos(theta)-(3.0/4.0*sqrt(r)*
		cos(3.0/2.0*theta)*sin(theta)-3.0/4.0*sqrt(r)*sin(3.0/2.0*theta)*cos(theta))*sin(theta)/r;
      he(1,0) = (3.0/4.0/sqrt(r)*sin(3.0/2.0*theta)*sin(theta)+3.0/4.0/sqrt(r)*
		cos(3.0/2.0*theta)*cos(theta))*sin(theta)+(3.0/4.0*sqrt(r)*cos(3.0/2.0*theta)
	       *sin(theta)-3.0/4.0*sqrt(r)*sin(3.0/2.0*theta)*cos(theta))*cos(theta)/r;
      he(1,1) = (3.0/4.0/sqrt(r)*cos(3.0/2.0*theta)*sin(theta)-3.0/4.0/sqrt(r)*sin(3.0/2.0*theta)
		 *cos(theta))*sin(theta)+(-3.0/4.0*sqrt(r)*sin(3.0/2.0*theta)*sin(theta)
                  -3.0/4.0*sqrt(r)*cos(3.0/2.0*theta)*cos(theta))*cos(theta)/r;
      he(0,1) = he(1,0) ;
    } break;
    case 3: {
      he(0,0) = (3.0/4.0/sqrt(r)*cos(theta/2.0)*cos(theta)+1/sqrt(r)*sin(theta/2.0)*
		 sin(theta)/4.0)*cos(theta)-(-sqrt(r)*sin(theta/2.0)*cos(theta)/4.0-5.0/4.0
					     *sqrt(r)*cos(theta/2.0)*sin(theta))*sin(theta)/r;
      he(1,0) = (3.0/4.0/sqrt(r)*cos(theta/2.0)*cos(theta)+1/sqrt(r)*sin(theta/2.0)*
		 sin(theta)/4.0)*sin(theta)+(-sqrt(r)*sin(theta/2.0)*cos(theta)/4.0-5.0/4.0
					     *sqrt(r)*cos(theta/2.0)*sin(theta))*cos(theta)/r;
      he(1,1) = (3.0/4.0/sqrt(r)*cos(theta/2.0)*sin(theta)-1/sqrt(r)*sin(theta/2.0)*
		 cos(theta)/4.0)*sin(theta)+(-sqrt(r)*sin(theta/2.0)*sin(theta)/4.0+5.0/4.0
					     *sqrt(r)*cos(theta/2.0)*cos(theta))*cos(theta)/r;
      he(0,1) = he(1,0) ;
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
  
  virtual void hess(const getfem::fem_interpolation_context& c, base_matrix &he) const {
    update_mls(c.convex_num());
    size_type P = c.xref().size(); assert(P == 2);
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
    assert(d2x(0,1) == d2x(1,0));
    assert(d2y(0,1) == d2y(1,0));
    // expression of the derivatives in the element reference basis
    for (unsigned i=0; i < 2; ++i) {
      for (unsigned j=0; j < 2; ++j) {
	hh[i*2+j] = (hefr(0, 0) * (dx[i]*dx[j]) + 
		     hefr(0, 1) * (dx[i]*dy[j]) +
		     hefr(1, 0) * (dx[j]*dy[i]) +
		     hefr(1, 1) * (dy[i]*dy[j]) + 
		     dfr[0]*d2x(i,j) + dfr[1]*d2y(i,j));
      }
    }
    base_vector &vhe = he;
    // back in the real element basis
    gmm::mult(c.B3(), hh, vhe);
    gmm::mult_add(c.B32(), gmm::scaled(dref,-1), vhe);
  }
    
  void update_from_context(void) const { cv =  size_type(-1); }

  bilaplacian_singular_functions(size_type l_, const getfem::level_set &ls_) : l(l_), ls(ls_) {
    cv = size_type(-1);
    this->add_dependency(ls);
  }

};

getfem::pglobal_function bilaplacian_crack_singular(size_type i, const getfem::level_set &ls){ 
  return new bilaplacian_singular_functions(i, ls);
}

struct exact_solution {
  getfem::mesh_fem_global_function mf;
  getfem::base_vector U;

  exact_solution(getfem::mesh &me) : mf(me) {}
  
  void init(getfem::level_set &ls) {
    std::vector<getfem::pglobal_function> cfun(4) ;
    for (unsigned j=0; j < 4; ++j)
      cfun[j] = bilaplacian_crack_singular(j, ls) ;
    mf.set_functions(cfun);
    U.resize(4); assert(mf.nb_dof() == 4);
   // scalar_type A1 = 1., nu = 0.3 ;
   // scalar_type b1_ = 3. + (A2 / A1) * (24. * nu) / (3. * nu * nu - 6. * nu + 5. ) ; 
    U[0] = BB ;
    U[1] = AA ;
    U[2] = CC ;
    U[3] = DD ;
  }
  


};

std::string name_of_dof(getfem::pdof_description dof) {
  char s[200];
  sprintf(s, "UnknownDof[%p]", (void*)dof);
  for (unsigned d = 0; d < 4; ++d) {
    if (dof == getfem::lagrange_dof(d)) {
      sprintf(s, "Lagrange[%d]", d); goto found;
    }
    if (dof == getfem::normal_derivative_dof(d)) {
      sprintf(s, "D_n[%d]", d); goto found;
    }
    if (dof == getfem::global_dof(d)) {
      sprintf(s, "GlobalDof[%d]", d);
    }
    if (dof == getfem::mean_value_dof(d)) {
      sprintf(s, "MeanValue[%d]", d);
    }
    if (getfem::dof_xfem_index(dof) != 0) {
      sprintf(s, "Xfem[idx:%d]", int(dof_xfem_index(dof)));
    }
    
    for (unsigned r = 0; r < d; ++r) {
      if (dof == getfem::derivative_dof(d, r)) {
	sprintf(s, "D_%c[%d]", "xyzuvw"[r], d); goto found;
      }
      for (unsigned t = 0; t < d; ++t) {
	if (dof == getfem::second_derivative_dof(d, r, t)) {
	  sprintf(s, "D2%c%c[%d]", "xyzuvw"[r], "xyzuvw"[t], d); 
	  goto found;
	}
      }
    }
  }
 found:
  return s;
}
namespace getfem {
// function for assembling the constraints of the integral matching 
  template<typename MAT, typename VECT1, typename VECT2>
  void asm_normal_derivative_dirichlet_constraints_bis 
 (MAT &H, VECT1 &R, const mesh_im &mim, const mesh_fem &mf_u,
   const mesh_fem &mf_mult, const mesh_fem &mf_r,
   const VECT2 &r_data, const mesh_region &rg, bool R_must_be_derivated, 
   int version) {
    typedef typename gmm::linalg_traits<VECT1>::value_type value_type;
    typedef typename gmm::number_traits<value_type>::magnitude_type magn_type;
    
    rg.from_mesh(mim.linked_mesh()).error_if_not_faces();
    GMM_ASSERT1(mf_r.get_qdim() == 1, 
		"invalid data mesh fem (Qdim=1 required)");
    if (version & ASMDIR_BUILDH) {
      const char *s;
      if (mf_u.get_qdim() == 1 && mf_mult.get_qdim() == 1)
	s = "M(#1,#2)+=comp(Grad(#1).Normal().Grad(#2).Normal())(:,i,i,:,j,j)";
      else
	s = "M(#1,#2)+=comp(vGrad(#1).Normal().vGrad(#2).Normal())(:,j,i,i,:,j,k,k);";
      
      getfem::generic_assembly assem(s);
      assem.push_mi(mim);
      assem.push_mf(mf_mult);
      assem.push_mf(mf_u);
      assem.push_mat(H);
      assem.assembly(rg);
      gmm::clean(H, gmm::default_tol(magn_type())
		 * gmm::mat_maxnorm(H) * magn_type(1000));
    }
    if (version & ASMDIR_BUILDR) {
      if (!R_must_be_derivated) {
	asm_normal_source_term(R, mim, mf_mult, mf_r, r_data, rg);
      } else {
	asm_real_or_complex_1_param
	  (R, mim, mf_mult, mf_r, r_data, rg,
	   "R=data(#2); V(#1)+=comp(Grad(#1).Normal().Grad(#2).Normal())(i,j,k,k).R(j)");
      }
    }
  }
} // end getfem

/******** Struct for the Bilaplacian problem *******************************/


struct bilaplacian_crack_problem {

  enum { CLAMPED_BOUNDARY_NUM = 0, SIMPLE_SUPPORT_BOUNDARY_NUM = 1,
	 FORCE_BOUNDARY_NUM = 2, MOMENTUM_BOUNDARY_NUM = 3,
	 MORTAR_BOUNDARY_IN = 40, MORTAR_BOUNDARY_OUT = 41};
  
  getfem::mesh mesh;        /* the mesh */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::mesh_level_set mls;       
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u, mf_pre_mortar, mf_pre_mortar_deriv ; 
  getfem::mesh_fem_level_set mfls_u, mfls_mortar, mfls_mortar_deriv ; 
  getfem::mesh_fem_global_function mf_sing_u ;
  getfem::mesh_fem mf_partition_of_unity ;
  getfem::mesh_fem_product mf_u_product ;
  getfem::mesh_fem_sum mf_u_sum ;
  getfem::mesh_fem& mf_u() { return mf_u_sum; }
 
  getfem::mesh_fem mf_rhs;  /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_mult; /* mesh_fem for the Dirichlet condition.        */
  
  scalar_type residual;     /* max residual for the iterative solvers       */
  getfem::constraints_type dirichlet_version;
  
  scalar_type cutoff_radius, enr_area_radius;
  int enrichment_option;
  size_type NX;
  
  std::string datafilename;
  bgeot::md_param PARAM;

  bool KL;

  dal::bit_vector pm_convexes; /* convexes inside the enrichment 
				  area when point-wise matching is used.*/
  
  scalar_type epsilon ;      /* half-plate thickness */
  
  exact_solution exact_sol;
    
  bool solve(plain_vector &U);
  void init(void);
  void compute_error(plain_vector &U);

  bilaplacian_crack_problem(void) : ls(mesh, 1, true), 
				    mls(mesh), mim(mls), mf_pre_u(mesh),  
				    mf_pre_mortar(mesh), mf_pre_mortar_deriv(mesh), mfls_u(mls, mf_pre_u), 
				    mfls_mortar(mls, mf_pre_mortar), mfls_mortar_deriv(mls, mf_pre_mortar_deriv),  
				    mf_sing_u(mesh), mf_partition_of_unity(mesh),
				    mf_u_product(mf_partition_of_unity, mf_sing_u), mf_u_sum(mesh),
				    mf_rhs(mesh), mf_mult(mesh), exact_sol(mesh)   
				    { KL = true; } 
};


/*                                                          */
/*****  Methods for class bilaplacian_crack_problem  ********/
/*                                                          */





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
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
				      "Enrichment option");
    enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
    
    /* First step : build the mesh */
    
    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor(MESH_TYPE);
    N = pgt->dim();
    GMM_ASSERT1(N == 2, "For a plate problem, N should be 2");
    std::vector<size_type> nsubdiv(N);
    NX = PARAM.int_value("NX", "Number of space steps ") ;
    std::fill(nsubdiv.begin(),nsubdiv.end(), NX);
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
  KL = (PARAM.int_value("KL", "Kirchhoff-Love model or not") != 0);
  D = PARAM.real_value("D", "Flexion modulus");
  if (KL) nu = PARAM.real_value("NU", "Poisson ratio");

 // Setting the integration methods
  
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method sppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
		getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
    
  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  mim.set_simplex_im(sppi, sing_ppi);
  
  /* Setting the finite element on the mf_u */  
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u); 
  getfem::pfem pf_partition_of_unity = getfem::fem_descriptor(PARAM.string_value("PARTITION_OF_UNITY_FEM_TYPE")) ; 
  mf_partition_of_unity.set_finite_element(mesh.convex_index(), pf_partition_of_unity);     
  mf_pre_mortar.set_finite_element(mesh.convex_index(),
             getfem::fem_descriptor(PARAM.string_value("MORTAR_FEM_TYPE")));
  mf_pre_mortar_deriv.set_finite_element(mesh.convex_index(),
             getfem::fem_descriptor(PARAM.string_value("MORTAR_DERIV_FEM_TYPE")));
  
  // set the mesh_fem of the multipliers (for the dirichlet condition)    
  std::string dirichlet_fem_name = PARAM.string_value("DIRICHLET_FEM_TYPE");
  if (dirichlet_fem_name.size() == 0)
    mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  else {
    cout << "DIRICHLET_FEM_TYPE="  << dirichlet_fem_name << "\n";
    mf_mult.set_finite_element(mesh.convex_index(), 
			       getfem::fem_descriptor(dirichlet_fem_name));
  }

  /* set the finite element on mf_rhs (same as mf_u if DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
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
  
  exact_sol.init(ls);
}


/* compute the error with respect to the exact solution */
void bilaplacian_crack_problem::compute_error(plain_vector &U) {
// plain_vector V(gmm::vect_size(U)) ;
// gmm::clear(V) ;
//   cout << "L2 ERROR : "
//        << getfem::asm_L2_dist(mim, mf_u(), V,
// 			      exact_sol.mf, exact_sol.U) << "\n";
//   cout << "H1 ERROR : "
//        << getfem::asm_H1_dist(mim, mf_u(), V,
// 			      exact_sol.mf, exact_sol.U) << "\n";
//   cout << "H2 ERROR : "
//        << getfem::asm_H2_dist(mim, mf_u(), V, 
//                               exact_sol.mf, exact_sol.U) << "\n"; 

  cout << "L2 ERROR : "
       << getfem::asm_L2_dist(mim, mf_u(), U,
			      exact_sol.mf, exact_sol.U) << "\n";
  cout << "H1 ERROR : "
       << getfem::asm_H1_dist(mim, mf_u(), U,
			      exact_sol.mf, exact_sol.U) << "\n";
  cout << "H2 ERROR : "
       << getfem::asm_H2_dist(mim, mf_u(), U, 
                              exact_sol.mf, exact_sol.U) << "\n"; 

}


/**************************************************************************/
/*  Model.                                                                */
/**************************************************************************/

bool bilaplacian_crack_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  //size_type N = mesh.dim();
  
  // Setting the level-set
  ls.reinit();  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    scalar_type x = ls.get_mesh_fem().point_of_dof(d)[0];
    scalar_type y = ls.get_mesh_fem().point_of_dof(d)[1];
    ls.values(0)[d] = y;
    ls.values(1)[d] = x;
  }
  ls.touch();  
  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
  mfls_mortar.adapt();
  mfls_mortar_deriv.adapt();
  cout << "mfls_u.nb_dof()=" << mfls_u.nb_dof() << "\n";
  
  // setting singularities 
  cout << "setting singularities \n" ;
  std::vector<getfem::pglobal_function> ufunc(4);
  for (size_type i = 0 ; i < ufunc.size() ; ++i) {                              
    ufunc[i] = bilaplacian_crack_singular(i, ls);
  }
  mf_sing_u.set_functions(ufunc);
  
  
  // Setting the enrichment --------------------------------------------/
   
  switch(enrichment_option) {
  case 0 :  // No enrichment
    mf_u_sum.set_mesh_fems(mfls_u);
    break ;
  case 1 : 
    {
      cout << "\npointwise matching\n";
     /* first : selecting the convexes that are completly included in the enrichment area */ 
     for (dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i) {
       pm_convexes.add(i) ;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
       for (unsigned j=0; j < mesh.nb_points_of_convex(i); ++j) {
	 if (gmm::sqr(mesh.points_of_convex(i)[j][0]) + 
	     gmm::sqr(mesh.points_of_convex(i)[j][1]) > 
	     gmm::sqr(enr_area_radius)) {
	   pm_convexes.sup(i); break;
	 }
       }
     }
      
      for (dal::bv_visitor cv(mf_sing_u.convex_index()); !cv.finished(); ++cv) {
	if (!pm_convexes.is_in(cv))
	  mf_sing_u.set_finite_element(cv, 0);
      }
      cout << "mf_sing_u: convex_index() = " << mf_sing_u.convex_index().card() << " convexes\n";

      //mf_u_sum.set_mesh_fems(mfls_u_ext, mf_pre_u); //_ext, mf_sing_u);
      mf_u_sum.set_smart_global_dof_linking(true);
      mf_u_sum.set_mesh_fems(mf_pre_u, mf_sing_u);


      cout << "mf_u_sum.nb_dof = " << mf_u_sum.nb_dof() << "\n";
      cout << "mfls_u.convex_index = " << mfls_u.convex_index() << "\nmf_sing_u: " << mf_sing_u.convex_index() << "\n";
      
    } break ;
  case 2 :  // standard XFEM on a fixed zone
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
          cout << "enriched_dofs: " << enriched_dofs << "\n";
          if (enriched_dofs.card() < 3)
             GMM_WARNING0("There is " << enriched_dofs.card() <<
		   " enriched dofs for the crack tip");
          mf_u_product.set_enrichment(enriched_dofs);
          mf_u_sum.set_mesh_fems(mf_u_product, mfls_u);
          cout << "enrichment done \n" ;}
	  break ;
  case 3 : // Integral matching (mortar)
    {
    cout << "\nIntegral Matching (Mortar)\n" ;
    
    dal::bit_vector cvlist_in_area;
    dal::bit_vector cvlist_out_area;
    for (dal::bv_visitor cv(mesh.convex_index()); 
	   !cv.finished(); ++cv) {
	bool in_area = true;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
	for (unsigned j=0; j < mesh.nb_points_of_convex(cv); ++j) {
	  if (gmm::sqr(mesh.points_of_convex(cv)[j][0] ) + 
	      gmm::sqr(mesh.points_of_convex(cv)[j][1] ) > 
	      gmm::sqr(enr_area_radius)) {
	    in_area = false; break;
	  }
	}

	/* "remove" the global function on convexes outside the enrichment
	   area */
	if (!in_area) {
	  cvlist_out_area.add(cv);
	  mf_sing_u.set_finite_element(cv, 0);
	  mf_u().set_dof_partition(cv, 1);
	} else cvlist_in_area.add(cv);
      }
       
       /* extract the boundary of the enrichment area, from the
	 "inside" point-of-view, and from the "outside"
	 point-of-view */
      getfem::mesh_region r_border, r_enr_out;
      getfem::outer_faces_of_mesh(mesh, r_border);

      getfem::outer_faces_of_mesh(mesh, cvlist_in_area, 
				  mesh.region(MORTAR_BOUNDARY_IN));
      getfem::outer_faces_of_mesh(mesh, cvlist_out_area, 
				  mesh.region(MORTAR_BOUNDARY_OUT));
      for (getfem::mr_visitor v(r_border); !v.finished(); ++v) {
	mesh.region(MORTAR_BOUNDARY_OUT).sup(v.cv(), v.f());
      }
      mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
      //mf_u_sum.set_mesh_fems(mfls_u);

      cout << "cvlist_in_area: " << cvlist_in_area << "\n";
      cout << "mfls_u.nb_dof: " << mfls_u.nb_dof() << "\n";
      cout << "mf_u_sum.nb_dof: " << mf_u_sum.nb_dof() << "\n";
      cout << "MORTAR_BOUNDARY_IN: " << mesh.region(MORTAR_BOUNDARY_IN) << "\n";
      cout << "MORTAR_BOUNDARY_OUT: " << mesh.region(MORTAR_BOUNDARY_OUT) << "\n";
    }
    break ;
  default : 
	GMM_ASSERT1(false, "Enrichment_option parameter is undefined");
	break ;  
	}
  mesh.write_to_file("toto.mesh");
  
  /*
  unsigned Q = mf_u().get_qdim();
  for (unsigned d=0; d < mf_u().nb_dof(); d += Q) {
    printf("dof %4d @ %+6.2f:%+6.2f: ", d, 
    mf_u().point_of_dof(d)[0], mf_u().point_of_dof(d)[1]);
    
    
    const getfem::mesh::ind_cv_ct cvs = mf_u().convex_to_dof(d);
    for (unsigned i=0; i < cvs.size(); ++i) {
      unsigned cv = cvs[i];
      //if (pm_cvlist.is_in(cv)) flag1 = true; else flag2 = true;

      getfem::pfem pf = mf_u().fem_of_element(cv);
      unsigned ld = unsigned(-1);
      for (unsigned dd = 0; dd < mf_u().nb_dof_of_element(cv); dd += Q) {
	if (mf_u().ind_dof_of_element(cv)[dd] == d) {
	  ld = dd/Q; break;
	}
      }
      if (ld == unsigned(-1)) {
	cout << "DOF " << d << "NOT FOUND in " << cv << " BUG BUG\n";
      } else {
	printf(" %3d:%.16s", cv, name_of_dof(pf->dof_types().at(ld)).c_str());
      }
    }
    printf("\n");
  }
  */

  //cout << "validate mf_sing_u():\n"; validate_fem_derivatives(mf_sing_u);

  //cout << "validate mf_u():\n"; validate_fem_derivatives(mf_u());
  
  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  // Bilaplacian brick.
  getfem::mdbrick_bilaplacian<> BIL(mim, mf_u());
  BIL.D().set(D);
  if (KL) { BIL.set_to_KL(); BIL.nu().set(nu); }
  
  // Defining the normal derivative Dirichlet condition value.
  plain_vector F;


  /* WRONG !! 

    F.resize(nb_dof_rhs*N);
  getfem::interpolation_function(mf_rhs, F, sol_du, CLAMPED_BOUNDARY_NUM);  
   
  // Normal derivative Dirichlet condition brick. 
  getfem::mdbrick_normal_derivative_Dirichlet<>                   
    NDER_DIRICHLET(BIL, CLAMPED_BOUNDARY_NUM, mf_mult);       
 
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.rhs().set(mf_rhs, F);
  */

  //gmm::resize(U, mf_u().nb_dof());  return true;


  // Normal derivative Dirichlet condition brick. 
  getfem::mdbrick_normal_derivative_Dirichlet<>                   
    NDER_DIRICHLET(BIL, CLAMPED_BOUNDARY_NUM, mf_mult);    
  NDER_DIRICHLET.set_constraints_type(dirichlet_version);
  NDER_DIRICHLET.R_must_be_derivated(); // hence we give the exact solution , and its gradient will be taken
  NDER_DIRICHLET.rhs().set(exact_sol.mf,exact_sol.U);


  // Defining the Dirichlet condition value.
  gmm::resize(F, nb_dof_rhs);
  getfem::interpolation_function(mf_rhs, F, sol_u,SIMPLE_SUPPORT_BOUNDARY_NUM);

  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<>
    DIRICHLET(NDER_DIRICHLET, SIMPLE_SUPPORT_BOUNDARY_NUM, mf_mult);
  DIRICHLET.rhs().set(exact_sol.mf,exact_sol.U);
  DIRICHLET.set_constraints_type(getfem::constraints_type(dirichlet_version));  
  if (dirichlet_version == getfem::PENALIZED_CONSTRAINTS)
    DIRICHLET.set_penalization_parameter(PARAM.real_value("EPS_DIRICHLET_PENAL")) ;
  getfem::mdbrick_abstract<> *final_model = &DIRICHLET ;
  
  if (enrichment_option == 3) {
     /* add a constraint brick for the mortar junction between
       the enriched area and the rest of the mesh */
    getfem::mesh_fem &mf_mortar = mfls_mortar; 
    getfem::mesh_fem &mf_mortar_deriv = mfls_mortar_deriv;
    getfem::mdbrick_constraint<> &mortar = 
      *(new getfem::mdbrick_constraint<>(DIRICHLET, 0));
    
  //   /* build the list of dof for the "(u-v) lambda" and for the  "\partial_n(u-v) \partial_n lambda" term in the mortar condition */  
//     dal::bit_vector bv_mortar;
//     dal::bit_vector bv_deriv;
//     dal::bit_vector bv_union;
//     sparse_matrix MM(mf_mortar.nb_dof(), mf_mortar.nb_dof());
//     sparse_matrix MD(mf_mortar.nb_dof(), mf_mortar.nb_dof());
//     std::vector<size_type> ind_mortar;
//     std::vector<size_type> ind_deriv;
//     std::vector<size_type> ind_union;
//     getfem::asm_mass_matrix(MM, mim, mf_mortar, MORTAR_BOUNDARY_OUT);

//     getfem::base_vector R(gmm::mat_nrows(MD));
//     getfem::asm_normal_derivative_dirichlet_constraints_bis(MD, R, mim, mf_mortar,
// 	     		mf_mortar, mf_pre_mortar, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;

//     for (size_type i=0; i < mf_mortar.nb_dof(); ++i) {
//       if ( (MM(i,i) > 1e-15) & (MD(i,i) > 1e-15) ) { 
// 	bv_mortar.add(i);
// 	bv_deriv.add(i);
// 	bv_union.add(i);
//       }
//       if ( (MM(i,i) > 1e-15) & (MD(i,i) <= 1e-15) ) { 
// 	bv_mortar.add(i);
// 	bv_union.add(i);
//       }
//       if ( (MM(i,i) <= 1e-15) & (MD(i,i) > 1e-15) ) { 
// 	bv_deriv.add(i);
// 	bv_union.add(i);
//       }
//     }

//     //cout << "matrice des dérivées normales du mesh_fem_mortar : \n" << MD << "\n" ;

//     for (dal::bv_visitor d(bv_mortar); !d.finished(); ++d)
//       ind_mortar.push_back(d);
//     for (dal::bv_visitor d(bv_deriv); !d.finished(); ++d)
//       ind_deriv.push_back(d);
//     for (dal::bv_visitor d(bv_union); !d.finished(); ++d)
//       ind_union.push_back(d);

//   /* This command was not suitable here :
//     dal::bit_vector bv_mortar = 
//       mf_mortar.dof_on_region(MORTAR_BOUNDARY_OUT);  
//       The reason for that is that unfortunately, the method "dof_on_region"
//       sometimes return too much dof when the mesh_fem is enriched. */
  



//     /* building matrices */
//     cout << "Handling mortar junction (" << ind_union.size() << 
//       " dof for the lagrange multiplier)\n";
      
//     sparse_matrix H0(mf_mortar.nb_dof(), mf_u().nb_dof()),
//       H(ind_union.size(), mf_u().nb_dof());

//     cout << "bv_mortar = " << bv_mortar << "\n";
//     cout << "bv_deriv = " << bv_deriv << "\n" ;
//     cout << "bv_union = " << bv_union << "\n" ;
    
//     gmm::sub_index sub_i(ind_union);
//     gmm::sub_index sub_i1(ind_mortar);
//     gmm::sub_index sub_i2(ind_deriv);
//     gmm::sub_interval sub_j(0, mf_u().nb_dof());
//     // build sub_indices of dofs which are either value or derivatives in the matrix of constraints H
//     std::vector<size_type> ind_val_H, ind_deriv_H ;
//     for (unsigned i=0; i< ind_union.size(); ++i) {
//       if ( bv_mortar[ind_union[i]] ) ind_val_H.push_back(i) ;
//       if ( bv_deriv[ind_union[i]] ) ind_deriv_H.push_back(i) ;
//     }   
//     gmm::sub_index sub_val_H(ind_val_H) ;
//     gmm::sub_index sub_deriv_H(ind_deriv_H) ;

//     /* build the mortar constraint matrix -- note that the integration
//        method is conformal to the crack
//      */
//     getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_OUT);
//     gmm::copy(gmm::sub_matrix(H0, sub_i1, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );

//     gmm::clear(H0);
//     getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), MORTAR_BOUNDARY_IN);
//     gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i1, sub_j), -1), gmm::sub_matrix(H, sub_val_H, sub_j));
    
   
//     gmm::clear(H0) ;
//     gmm::resize(R, mat_nrows(H));
//     getfem::asm_normal_derivative_dirichlet_constraints_bis(H0, R, mim, mf_u(),
// 	     		mf_mortar, mf_pre_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
//     gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i2, sub_j), 1. ), gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
    
//     gmm::clear(H0) ;
//     getfem::asm_normal_derivative_dirichlet_constraints_bis(H0, R, mim, mf_u(),
// 			mf_mortar, mf_pre_u, R, MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
//     gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i2, sub_j), -1.), gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
    
//     // -----------------------------------
    
       
    /* New version of the integral matching -----------------------------------------------*/

// selecting nodes indices on the two meth. mult.

// assembling matrices, selceting the lines corresponding to those holding information

   /* build the list of dof for the "(u-v) lambda" condition  */  
    dal::bit_vector bv_mortar;
    dal::bit_vector bv_deriv;
    sparse_matrix MM(mf_pre_mortar.nb_dof(), mf_pre_mortar.nb_dof());
    sparse_matrix MD(mf_pre_mortar_deriv.nb_dof(), mf_pre_mortar_deriv.nb_dof());
    std::vector<size_type> ind_mortar;
    std::vector<size_type> ind_deriv;
    getfem::asm_mass_matrix(MM, mim, mf_pre_mortar, MORTAR_BOUNDARY_OUT);
    getfem::asm_mass_matrix(MD, mim, mf_pre_mortar_deriv, MORTAR_BOUNDARY_OUT);
    //  getfem::base_vector R( mf_pre_mortar_deriv.nb_dof() );
    //  getfem::asm_normal_derivative_dirichlet_constraints_bis(MD, R, mim, mf_pre_mortar_deriv,
    // 	  		mf_pre_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
    for (size_type i=0; i < mf_pre_mortar.nb_dof(); ++i)
      if (MM(i,i) > 1e-15) bv_mortar.add(i);
    for (size_type i=0; i < mf_pre_mortar_deriv.nb_dof(); ++i)
      if (MD(i,i) > 1e-15) bv_deriv.add(i);
    
    for (dal::bv_visitor d(bv_mortar); !d.finished(); ++d)
      ind_mortar.push_back(d);
    for (dal::bv_visitor d(bv_deriv); !d.finished(); ++d)
      ind_deriv.push_back(d);

    // building matrices

    cout << "Handling mortar junction (" << ind_mortar.size() << 
      " dof for the lagrange multiplier of the displacement, " <<
      ind_deriv.size() << " dof for the lagrange multiplier of the derivative)\n";
      
    sparse_matrix H0(mf_pre_mortar.nb_dof(), mf_u().nb_dof()),
      H(ind_mortar.size() + ind_deriv.size(), mf_u().nb_dof()) ; 

    // Defining sub_indexes of the matrices calculated with the 
    // complete set of dofs.
    gmm::sub_index sub_i(ind_mortar);
    gmm::sub_index sub_i1(ind_deriv);
    gmm::sub_interval sub_j(0, mf_u().nb_dof());

    gmm::sub_interval sub_val_H(0, ind_mortar.size()) ;
    gmm::sub_interval sub_deriv_H(ind_mortar.size(), ind_deriv.size()) ;

    cout << "sub_indexes built\n" ;
    /* build the mortar constraint matrix -- note that the integration
       method is conformal to the crack
     */
    getfem::asm_mass_matrix(H0, mim, mf_pre_mortar, mf_u(), MORTAR_BOUNDARY_OUT);
    gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), gmm::sub_matrix(H, sub_val_H, sub_j) );

    gmm::clear(H0);
    getfem::asm_mass_matrix(H0, mim, mf_pre_mortar, mf_u(), MORTAR_BOUNDARY_IN);
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i, sub_j), -1), gmm::sub_matrix(H, sub_val_H, sub_j) );


    cout << "first contraint asm\n" ;
    base_vector R;
    gmm::clear(H0);
    gmm::resize(H0, mf_pre_mortar_deriv.nb_dof(), mf_u().nb_dof() ) ;
    getfem::asm_normal_derivative_dirichlet_constraints
      (H0, R, mim, mf_u(), mf_pre_mortar_deriv, mf_pre_u, base_vector(),
       MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
    gmm::add(gmm::sub_matrix(H0, sub_i1, sub_j),
             gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
   
    cout << "first step \n" ;
    gmm::clear(H0);
    getfem::asm_normal_derivative_dirichlet_constraints
      (H0, R, mim, mf_u(), mf_pre_mortar_deriv, mf_pre_u, R,
       MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
    gmm::add(gmm::sub_matrix(H0, sub_i1, sub_j),
             gmm::sub_matrix(H, sub_deriv_H, sub_j)) ;
    
    
    // ------------------------------------------ end of new version

    /* because of the discontinuous partition of mf_u(), some levelset 
       enriched functions do not contribute any more to the
       mass-matrix (the ones which are null on one side of the
       levelset, when split in two by the mortar partition, may create
       a "null" dof whose base function is all zero.. */
    sparse_matrix M2(mf_u().nb_dof(), mf_u().nb_dof());
    getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());
    //gmm::HarwellBoeing_IO::write("M2.hb", M2);
    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
      if (M2(d,d) < 1e-15) {
	cout << "  removing null mf_u() dof " << d << " @ " << 
	  mf_u().point_of_dof(d) << "\n";	
	unsigned n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }  
    gmm::resize(R, gmm::mat_nrows(H)); 
    mortar.set_constraints(H,R);
    final_model = &mortar;
    gmm::HarwellBoeing_IO::write("H.hb", H);        
    //------------------------------------------------------------------
    // Matching of the normal derivative
//     getfem::mdbrick_constraint<> &mortar_derivative = 
//       *(new getfem::mdbrick_constraint<>(mortar,0));
//       
//     gmm::clear(H0) ;
//     gmm::clear(H) ;
//     gmm::resize(H, ind_mortar.size(), mf_u().nb_dof());
//     getfem::asm_normal_derivative_dirichlet_constraints(H0, R, mim, mf_u(),
// 	     		mf_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_OUT, 0, getfem::ASMDIR_BUILDH) ;
//     gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), H) ;
//     
//     gmm::clear(H0) ;
//     getfem::asm_normal_derivative_dirichlet_constraints(H0, R, mim, mf_u(),
// 			mf_mortar_deriv, mf_pre_u, R, MORTAR_BOUNDARY_IN, 0, getfem::ASMDIR_BUILDH) ;
//     gmm::add(gmm::sub_matrix(H0, sub_i, sub_j), H) ;
//     
//     //cout << "matrice de contraintes : \n" << H << "\n" ;
//     cout << "matching of derivative values done.\n" ;  
//     gmm::resize(R, ind_mortar.size());
//     mortar_derivative.set_constraints(H,R);
//     final_model = &mortar_derivative ;
  }
  

  // Generic solve.
  cout << "Total number of variables : " << final_model->nb_dof() << endl;
  getfem::standard_model_state MS(*final_model);
  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(MS, *final_model, iter);

  // Solution extraction
  gmm::resize(U, mf_u().nb_dof());
  gmm::copy(BIL.get_solution(MS), U);
  return (iter.converged());
}


//function to save a vector for matlab
template<typename VEC> static void vecsave(std::string fname, const VEC& V) {
  std::ofstream f(fname.c_str()); f.precision(16);
  for (size_type i=0; i < V.size(); ++i) f << V[i] << "\n"; 
}

int main(int argc, char *argv[]) {

  try {
    bilaplacian_crack_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    plain_vector U;
    p.mesh.write_to_file("mesh.m") ;
    if (!p.solve(U)) GMM_ASSERT1(false, "Solve has failed");

    p.compute_error(U);

    // visualisation, export au format .vtk
    getfem::mesh mcut;
    p.mls.global_cut_mesh(mcut);
    unsigned Q = p.mf_u().get_qdim();
    assert( Q == 1 ) ;
    getfem::mesh_fem mf(mcut, Q);
    mf.set_classical_discontinuous_finite_element(2, 0.001);
    // mf.set_finite_element
    //	(getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 2, 0.0001)"));
    plain_vector V(mf.nb_dof());

    getfem::interpolation(p.mf_u(), mf, U, V);

    getfem::stored_mesh_slice sl;
    getfem::mesh mcut_refined;

    unsigned NX = p.PARAM.int_value("NX"), nn;
    if (NX < 6) nn = 12;
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
	//d = gmm::vect_norm2(ls_function(P));
	d = gmm::vect_norm2(P) ; 
	if (d < dmin || i == 0) { dmin = d; Pmin = P; }
      }

      if (dmin < 1e-5)
	nrefine[cv] = nn*2;
      else if (dmin < .1) 
	nrefine[cv] = nn*2;
      else nrefine[cv] = nn;
      /*if (dmin < .01)
	cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin << " " << nrefine[cv] << "\n";*/
    }


    getfem::mesh_slicer slicer(mcut); 
    getfem::slicer_build_mesh bmesh(mcut_refined);
    slicer.push_back_action(bmesh);
    slicer.exec(nrefine, getfem::mesh_region::all_convexes());
     
    /*
      sl.build(mcut, 
      getfem::slicer_build_mesh(mcut_refined), nrefine);*/

    getfem::mesh_im mim_refined(mcut_refined); 
    mim_refined.set_integration_method(getfem::int_method_descriptor
				       ("IM_TRIANGLE(6)"));

    getfem::mesh_fem mf_refined(mcut_refined, Q);
    mf_refined.set_classical_discontinuous_finite_element(2, 0.001);
    plain_vector W(mf_refined.nb_dof());

    
    getfem::interpolation(p.mf_u(), mf_refined, U, W);


    int VTK_EXPORT = p.PARAM.int_value("VTK_EXPORT");
    if (VTK_EXPORT ) {
      cout << "exporting solution to " << p.datafilename + ".vtk" << "..\n";
      getfem::vtk_export exp(p.datafilename + ".vtk", false);
      exp.exporting(mf_refined); 
      exp.write_point_data(mf_refined, W, "vertical_displacement");
      cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename  << ".vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
    }
    if ((VTK_EXPORT & (2|4))) {
      p.exact_sol.mf.set_qdim(Q);
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());   
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);
	
      plain_vector DIFF(EXACT); gmm::add(gmm::scaled(W,-1),DIFF);
      if ((VTK_EXPORT & 2)) {
	cout << "exporting exact solution to VTK\n";
	getfem::vtk_export exp2(p.datafilename + "_exact.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined, EXACT, "reference solution");
	cout << "export done, you can view the data file with (for example)\n"
	"mayavi -d " << p.datafilename << "_exact.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
      }
      if ((VTK_EXPORT & 4)) {
	cout << "exporting difference with exact solution to VTK\n";
	getfem::vtk_export exp2(p.datafilename + "_diff.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined, DIFF, "difference solution");
	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename  << "_diff.vtk -f "
	"WarpScalar -m BandedSurfaceMap -m Outline\n";
      }
    }

    int MATLAB_EXPORT = p.PARAM.int_value("MATLAB_EXPORT");
    if (MATLAB_EXPORT) {
      cout << "exporting solution to " << p.datafilename + ".mf" << " and " << p.datafilename << ".U\n";
      mf_refined.write_to_file(p.datafilename + ".mf", true);
      vecsave(p.datafilename + ".U", W);
      p.exact_sol.mf.set_qdim(Q);
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());   // ??
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);
      vecsave(p.datafilename + ".EXACT", EXACT);

      cout << "exporting original mesh to " << p.datafilename + "_base.mesh\n";
      p.mesh.write_to_file(p.datafilename + "_base.mesh");

      cout << "matlab export done, you can view the results with\n";
      cout << "mf=gfMeshFem('load', 'bilaplacian.mf'); U=load('bilaplacian.U')'; "
	"EXACT=load('bilaplacian.EXACT')'; m0=gfMesh('load','bilaplacian_base.mesh');\n";
      cout << "gf_plot(mf, U-EXACT, 'refine',1); hold on; gf_plot_mesh(m0); hold off; colorbar;\n";
    }

    
    //getchar(); 
    
  }
  GMM_STANDARD_CATCH_ERROR;
  return 0; 
}
























