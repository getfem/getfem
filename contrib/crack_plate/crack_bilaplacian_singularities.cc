// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2008 Yves Renard, Julien Pommier.
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
#include "crack_bilaplacian.h"
#include "crack_bilaplacian.h"

void bilaplacian_singular_functions::update_mls(size_type cv_) const { 
  if (cv_ != cv) { 
    cv=cv_; 
    mls0=ls.mls_of_convex(cv, 0); 
    mls1=ls.mls_of_convex(cv, 1); 
  }
}


scalar_type bilaplacian_singular_functions::sing_function(scalar_type x, scalar_type y) const {
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
void bilaplacian_singular_functions::sing_function_grad(scalar_type x, scalar_type y, base_small_vector &g) const {
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
   
void bilaplacian_singular_functions::sing_function_hess(scalar_type x, scalar_type y, base_matrix &he) const {
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
   
  
scalar_type bilaplacian_singular_functions::val(const getfem::fem_interpolation_context& c) const {

  assert(ls.get_mesh_fem().convex_index().is_in(c.convex_num()));
  update_mls(c.convex_num());
  scalar_type x = mls1(c.xref()), y = mls0(c.xref());
  scalar_type v = sing_function(x, y);
  return v;
}

void bilaplacian_singular_functions::grad(const getfem::fem_interpolation_context& c,
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
  
void bilaplacian_singular_functions::hess(const getfem::fem_interpolation_context& c, base_matrix &he) const {
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
    
void bilaplacian_singular_functions::update_from_context(void) const { cv =  size_type(-1); }

bilaplacian_singular_functions::bilaplacian_singular_functions(size_type l_, const getfem::level_set &ls_) : l(l_), ls(ls_) {
  cv = size_type(-1);
  this->add_dependency(ls);
}
