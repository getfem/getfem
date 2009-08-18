// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2007-2009 Yves Renard, Julien Pommier.
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

#include "crack_exact_solution.h"

/* returns sin(theta/2) where theta is the angle
   of 0-(x,y) with the axis Ox */
static scalar_type sint2(scalar_type x, scalar_type y) {
  scalar_type r = sqrt(x*x+y*y);
  if (r == 0) return 0;
  else return (y<0 ? -1:1) * sqrt(gmm::abs(r-x)/(2*r));
  // sometimes (gcc3.3.2 -O3), r-x < 0 ....
}
static scalar_type cost2(scalar_type x, scalar_type y) {
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

/*static void elasticite2lame(const scalar_type young_modulus,
			    const scalar_type poisson_ratio, 
			    scalar_type& lambda, scalar_type& mu) {
  mu = young_modulus/(2*(1+poisson_ratio));
  lambda = 2*mu*poisson_ratio/(1-poisson_ratio);
}
*/

/* solution reguliere pour une fracture en y=0 */
static base_small_vector solution_P2(double lambda, double mu, const base_node& pos, base_matrix* pgrad = 0) {
  double x=pos[0], y = pos[1];
  base_small_vector v(2);
  double c[6] = {0, 0, 0., 0., 0., 0.01};
  double d[6] = {0, .2, .1, 0.05, -.03, 0.06};
  v[0] = c[0]-1/lambda*(lambda+2.0*mu)*d[2]*x-d[1]*y-1/lambda*(lambda+2.0*mu)*d[4]*x*x/2.0-2.0*(lambda+2.0*mu)*d[5]/lambda*x*y+(3.0*lambda+4.0*mu)*d[4]/lambda*y*y/2.0;
  v[1]= d[0]+d[1]*x+d[2]*y+1/lambda*(lambda+2.0*mu)*d[5]*x*x+d[4]*x*y+d[5]*y*y;
  if (pgrad) {
    (*pgrad).resize(2,2);
    (*pgrad)(0,0) = -1/lambda*(lambda+2.0*mu)*d[2]-1/lambda*(lambda+2.0*mu)*d[4]*x-2.0*(lambda+2.0*mu)*d[5]/lambda*y;
    (*pgrad)(1,0) = -d[1]-2.0*(lambda+2.0*mu)*d[5]/lambda*x+(3.0*lambda+4.0*mu)*d[4]/lambda*y;
    (*pgrad)(0,1) = d[1]+d[4]*y+2/lambda*(lambda+2.0*mu)*d[5]*x;
    (*pgrad)(1,1) = d[2]+d[4]*x+2*d[5]*y;
  }
  return v;
}

static void sol_ref_infinite_plane(scalar_type lambda, scalar_type mu, 
				   scalar_type x, scalar_type y,
				   base_small_vector& U, int mode,
				   base_matrix *pgrad) {
  //scalar_type KI = sigma*sqrt(M_PI*a);
  scalar_type r = std::max(sqrt(x*x+y*y),1e-16);
  scalar_type sqrtr = sqrt(r), sqrtr3 = sqrtr*sqrtr*sqrtr;
  scalar_type cost = x/r, sint = y/r;
  scalar_type theta = atan2(y,x);
  scalar_type s2 = sin(theta/2); //sint2(x,y);
  scalar_type c2 = cos(theta/2); //cost2(x,y);
  // scalar_type c3 = cos(3*theta/2); //4*c2*c2*c2-3*c2; /* cos(3*theta/2) */
  // scalar_type s3 = sin(3*theta/2); //4*s2*c2*c2-s2;  /* sin(3*theta/2) */

  U.resize(2);
  if (pgrad) (*pgrad).resize(2,2);
  scalar_type C=1; //1./E * (mode == 1 ? 1. : (1+nu));
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
    cerr << "raaah not a number ...\n";
  assert(!isnan(U[0]));
  assert(!isnan(U[1]));
}


base_small_vector crack_exact_solution_function::eval(const base_node &x, base_matrix *pgrad) const {
  int N = x.size(); base_small_vector res(N); res.fill(0.);
  gmm::clear(res); if (pgrad) gmm::clear(*pgrad);
  if (function_num == 1) {
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], res, 1, pgrad);
  } else if (function_num == 2) {
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], res, 2, pgrad);
  } else if (function_num == 3) {
    //res[0] = (x[1]+2.5)/10;
    res[1] = x[1]/4;
  } else if (function_num == 4) {
    base_small_vector modeI(2), modeII(2);
    base_matrix gradI, gradII;
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeI, 1, (pgrad) ? &gradI : 0);
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeII, 2, (pgrad) ? &gradII : 0);
    res = solution_P2(lambda, mu, x, pgrad);
    res += modeI*4. - modeII*3.;
    if (pgrad) { 
      gmm::add(gmm::scaled(gradI,4.),*pgrad);
      gmm::add(gmm::scaled(gradII,-3.),*pgrad); 
    }
  } else if (function_num == 5) {
    base_small_vector modeI(2), modeII(2), modeZ(2);
    base_matrix gradI, gradII, gradZ;
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeI, 1, (pgrad) ? &gradI : 0);
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeII, 2, (pgrad) ? &gradII : 0);
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeZ, 100, (pgrad) ? &gradZ : 0);
    res = solution_P2(lambda, mu, x, pgrad); 
    res = -1. * res; if (pgrad) { for (size_type k=0; k < (*pgrad).size(); ++k) (*pgrad)[k] = -(*pgrad)[k]; }
    res += modeI*(1.) - modeII*3. - modeZ*.5;
    if (pgrad) { 
      gmm::add(gmm::scaled(gradI,(1.)),*pgrad);
      gmm::add(gmm::scaled(gradII,-3.),*pgrad); 
      gmm::add(gmm::scaled(gradZ,-.5),*pgrad);
    }
  } else if (function_num == 6) {
    base_small_vector modeI(2), modeII(2), modeZI(2), modeZII(2);
    base_matrix gradI, gradII, gradZI, gradZII;
    scalar_type coef[5] = {0.5,1,1.5,1,0};
    res = solution_P2(lambda, mu, x, pgrad); 
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeI, 1, (pgrad) ? &gradI : 0);
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeII, 2, (pgrad) ? &gradII : 0);
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeZI, 100, (pgrad) ? &gradZI : 0);
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], modeZII, 101, (pgrad) ? &gradZII : 0);
    //res = -1. * res; if (pgrad) { for (size_type k=0; k < (*pgrad).size(); ++k) (*pgrad)[k] = -(*pgrad)[k]; }
    res = res*coef[0] + modeI*coef[1] + modeII*coef[2] + modeZI*coef[3] + modeZII*coef[4];
    
    if (pgrad) { 
      if (pgrad) { for (size_type k=0; k < (*pgrad).size(); ++k) (*pgrad)[k] *= coef[0]; }
      gmm::add(gmm::scaled(gradI,coef[1]),*pgrad);
      gmm::add(gmm::scaled(gradII,coef[2]),*pgrad);
      gmm::add(gmm::scaled(gradZI,coef[3]),*pgrad);
      gmm::add(gmm::scaled(gradZII,coef[4]),*pgrad);
    }
  } else if (function_num == 7) {
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], res, 100, pgrad);
    res /= -40; if (pgrad) gmm::scale(*pgrad, -1./40);
  } else if (function_num == 8) {
    sol_ref_infinite_plane(lambda, mu, x[0], x[1], res, 101, pgrad);
    res /= 200; if (pgrad) gmm::scale(*pgrad, 1./200);
  }
  return res;
}

scalar_type crack_exact_solution_function::val(scalar_type x, scalar_type y) const {
  return eval(base_node(x,y), 0)[component_num];
}

base_small_vector crack_exact_solution_function::grad(scalar_type x, scalar_type y) const {
  base_matrix gr(2,2);
  eval(base_node(x,y), &gr);
  base_small_vector v(2);
  for (unsigned i=0; i < 2; ++i) v[i] = gr(i,component_num);
  return v;
}



void crack_exact_solution::init(int function_num, scalar_type lambda, scalar_type mu,
				getfem::level_set &ls) {
  std::vector<getfem::pglobal_function> cfun(2);
  for (size_type i = 0; i < 2; ++i) {
    /* use the singularity */
    getfem::abstract_xy_function *s = 
      new crack_exact_solution_function(function_num, 
					i, lambda, mu);
    cfun[i] = getfem::global_function_on_level_set(ls, *s);
  }
  
  mf.set_functions(cfun);
  
  mf.set_qdim(1);
  
  U.resize(4); assert(mf.nb_dof() == 2);
  U[0] = 1; U[1] = 0;
  U[2] = 0; U[3] = 1;
}
