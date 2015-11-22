/*===========================================================================

 Copyright (C) 2004-2015 Yves Renard

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

#include <getfem/getfem_mesh_fem_global_function.h>
#include <getfem/getfem_level_set.h>

namespace getfem {

  void global_function_fem::init(void) {
    is_pol = is_lag = false; es_degree = 5;
    is_equiv = real_element_defined = true;
    ntarget_dim = 1; // An extension for vectorial elements should be easy
    dim_ = cvr->structure()->dim();
    init_cvs_node();
    base_node P(dim()); P.fill(1./30.);
    for (size_type i = 0; i < functions.size(); ++i) {
      add_node(global_dof(dim()), P);
    }
    std::stringstream nm;
    nm << "GLOBAL_FEM(" << (void*)this << ")";
    debug_name_ = nm.str();
  }

  size_type global_function_fem::nb_dof(size_type) const
  { return functions.size(); }

  size_type global_function_fem::index_of_global_dof
  (size_type /*cv*/, size_type i) const
  { return i; }

  void global_function_fem::base_value(const base_node &, base_tensor &) const
  { GMM_ASSERT1(false, "No base values, real only element."); }
  void global_function_fem::grad_base_value(const base_node &,
                                            base_tensor &) const
  { GMM_ASSERT1(false, "No grad values, real only element."); }
  void global_function_fem::hess_base_value(const base_node &,
                                            base_tensor &) const
  { GMM_ASSERT1(false, "No hess values, real only element."); }

  void global_function_fem::real_base_value(const fem_interpolation_context& c,
                                            base_tensor &t, bool) const {
    mib.resize(2); mib[0] = target_dim();
    mib[1] = short_type(functions.size());
    assert(target_dim() == 1);
    t.adjust_sizes(mib);
    for (size_type i=0; i < functions.size(); ++i) {
      /*cerr << "global_function_fem: real_base_value(" << c.xreal() << ")\n";
      if (c.have_G()) cerr << "G = " << c.G() << "\n";
      else cerr << "no G\n";*/
      t[i] = (*functions[i]).val(c);
    }
  }

  void global_function_fem::real_grad_base_value
  (const fem_interpolation_context& c, base_tensor &t, bool) const {
    mig.resize(3);
    mig[2] = dim(); mig[1] = target_dim();
    mig[0] = short_type(functions.size());
    assert(target_dim() == 1);
    t.adjust_sizes(mig);
    base_small_vector G(dim());
    for (size_type i=0; i < functions.size(); ++i) {
      (*functions[i]).grad(c,G);
      for (unsigned j=0; j < dim(); ++j)
        t[j*functions.size() + i] = G[j];
    }
  }

  void global_function_fem::real_hess_base_value
  (const fem_interpolation_context &c, base_tensor &t, bool) const {
    mih.resize(4);
    mih[3] = mih[2] = dim(); mih[1] = target_dim();
    mih[0] = short_type(functions.size());
    assert(target_dim() == 1);
    t.adjust_sizes(mih);
    base_matrix H(dim(),dim());
    for (size_type i=0; i < functions.size(); ++i) {
      (*functions[i]).hess(c,H);
      for (unsigned k=0; k < dim(); ++k)
        for (unsigned j=0; j < dim(); ++j)
          t.at((k*dim() + j)*functions.size() + i) = H.at((k*dim() + j));
    }
  }

  DAL_SIMPLE_KEY(special_int_globf_fem_key, pfem);

  pfem new_global_function_fem(bgeot::pconvex_ref cvr,
                               const std::vector<pglobal_function> &f) {
    pfem pf = std::make_shared<global_function_fem>(cvr,f);
    dal::pstatic_stored_object_key
      pk = std::make_shared<special_int_globf_fem_key>(pf);
    dal::add_stored_object(pk, pf);
    return pf;
  }

  void mesh_fem_global_function::clear_build_methods() {
    for (std::map<bgeot::pconvex_ref,pfem>::const_iterator
         it = build_methods.begin();
         it != build_methods.end(); ++it)
      del_global_function_fem((*it).second);
    build_methods.clear();
  }
  void mesh_fem_global_function::clear(void) {
    mesh_fem::clear();
    clear_build_methods();
  }

  void mesh_fem_global_function::adapt(void) {
    clear();
    for (dal::bv_visitor cv(linked_mesh().convex_index());
         !cv.finished(); ++cv) {
      bgeot::pconvex_ref cvr =
       bgeot::basic_convex_ref(linked_mesh().trans_of_convex(cv)->convex_ref());

      std::map<bgeot::pconvex_ref,pfem>::iterator it = build_methods.find(cvr);
      pfem pf;
      if (it == build_methods.end()) {
        build_methods[cvr] = pf = new_global_function_fem(cvr,fun);
      } else pf = (*it).second;
      set_finite_element(cv, pf);
    }
    touch();
  }


  parser_xy_function::parser_xy_function(const std::string &sval,
					 const std::string &sgrad,
					 const std::string &shess)
    : f_val(gw_val, sval), f_grad(gw_grad, sgrad), f_hess(gw_hess, shess),
      ptx(1), pty(1), ptr(1), ptt(1) {

    gw_val.add_fixed_size_constant("x", ptx);
    gw_val.add_fixed_size_constant("y", pty);
    gw_val.add_fixed_size_constant("r", ptr);
    gw_val.add_fixed_size_constant("theta", ptt);
    
    gw_grad.add_fixed_size_constant("x", ptx);
    gw_grad.add_fixed_size_constant("y", pty);
    gw_grad.add_fixed_size_constant("r", ptr);
    gw_grad.add_fixed_size_constant("theta", ptt);

    gw_hess.add_fixed_size_constant("x", ptx);
    gw_hess.add_fixed_size_constant("y", pty);
    gw_hess.add_fixed_size_constant("r", ptr);
    gw_hess.add_fixed_size_constant("theta", ptt);

    f_val.compile(); f_grad.compile(); f_hess.compile();
  }
  
  scalar_type
  parser_xy_function::val(scalar_type x, scalar_type y) const {
    ptx[0] = double(x);                   // x
    pty[0] = double(y);                   // y
    ptr[0] = double(sqrt(fabs(x*x+y*y))); // r
    ptt[0] = double(atan2(y,x));          // theta

    const bgeot::base_tensor &t = f_val.eval();
    GMM_ASSERT1(t.size() == 1, "Wrong size of expression result "
                << f_val.expression());
    return t[0];
  }

  base_small_vector
  parser_xy_function::grad(scalar_type x, scalar_type y) const {
    ptx[0] = double(x);                   // x
    pty[0] = double(y);                   // y
    ptr[0] = double(sqrt(fabs(x*x+y*y))); // r
    ptt[0] = double(atan2(y,x));          // theta

    base_small_vector res(2);
    const bgeot::base_tensor &t = f_grad.eval();
    GMM_ASSERT1(t.size() == 2, "Wrong size of expression result "
                << f_grad.expression());
    gmm::copy(t.as_vector(), res);
    return res;
  }

  base_matrix
  parser_xy_function::hess(scalar_type x, scalar_type y) const {
    ptx[0] = double(x);                   // x
    pty[0] = double(y);                   // y
    ptr[0] = double(sqrt(fabs(x*x+y*y))); // r
    ptt[0] = double(atan2(y,x));          // theta

    base_matrix res(2,2);
    const bgeot::base_tensor &t = f_hess.eval();
    GMM_ASSERT1(t.size() == 4, "Wrong size of expression result "
                << f_hess.expression());
    gmm::copy(t.as_vector(), res.as_vector());
    return res;
  }

  /* the basic 4 singular functions for 2D cracks */
  scalar_type
  crack_singular_xy_function::val(scalar_type x, scalar_type y) const {
    scalar_type sgny = (y < 0 ? -1.0 : 1.0);
    scalar_type r = sqrt(x*x + y*y);

    if (r < 1e-10)  return 0;

    /* The absolute value is unfortunately necessary, otherwise, sqrt(-1e-16)
       can be required ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));

    scalar_type res = 0;
    switch(l){
   
  /* First order enrichement displacement field (linear elasticity) */    
      
      case 0 : res = sqrt(r)*sin2; break;
      case 1 : res = sqrt(r)*cos2; break;
      case 2 : res = sin2*y/sqrt(r); break;
      case 3 : res = cos2*y/sqrt(r); break;
  
  /* Second order enrichement of displacement field (linear elasticity) */
      
      case 4 : res = sqrt(r)*r*sin2; break;
      case 5 : res = sqrt(r)*r*cos2; break;
      case 6 : res = sin2*cos2*cos2*r*sqrt(r); break;
      case 7 : res = cos2*cos2*cos2*r*sqrt(r); break;
  
   /* First order enrichement of pressure field (linear elasticity) mixed formulation */

      case 8 : res = -sin2/sqrt(r); break;
      case 9 : res = cos2/sqrt(r); break;

  /* Second order enrichement of pressure field (linear elasticity) mixed formulation */

      case 10 : res = sin2*sqrt(r); break; // cos2*cos2
      case 11 : res = cos2*sqrt(r); break;

  /* First order enrichement of displacement field (nonlinear elasticity)[Rodney Stephenson Journal of elasticity VOL.12 No. 1, January 1982] */

      case 12 : res = r*sin2*sin2; break; 
      case 13 : res = sqrt(r)*sin2; break;

/* First order enrichement of pressure field (nonlinear elasticity)  */

      case 14 : res = sin2/r; break; 
      case 15 : res = cos2/r; break;
 
    default: GMM_ASSERT2(false, "arg");
    }
    return res;
  }


  base_small_vector
  crack_singular_xy_function::grad(scalar_type x, scalar_type y) const {
    scalar_type sgny = (y < 0 ? -1.0 : 1.0);
    scalar_type r = sqrt(x*x + y*y);

    if (r < 1e-10) {
      GMM_WARNING0("Warning, point close to the singularity (r=" << r << ")");
    }

    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));

    base_small_vector res(2);
    switch(l){
 /* First order enrichement displacement field (linear elasticity) */
    case 0 :
      res[0] = -sin2/(2*sqrt(r));
      res[1] = cos2/(2*sqrt(r));
      break;
    case 1 :
      res[0] = cos2/(2*sqrt(r));
      res[1] = sin2/(2*sqrt(r));
      break;
    case 2 :
      res[0] = cos2*((-5*cos2*cos2) + 1. + 4*(cos2*cos2*cos2*cos2))/sqrt(r);
      res[1] = sin2*((-3*cos2*cos2) + 1. + 4*(cos2*cos2*cos2*cos2))/sqrt(r);
      break;
    case 3 :
      res[0] = -cos2*cos2*sin2*((4*cos2*cos2) - 3.)/sqrt(r);
      res[1] = cos2*((4*cos2*cos2*cos2*cos2) + 2. - (5*cos2*cos2))/sqrt(r);
      break;
 /* Second order enrichement of displacement field (linear elasticity) */   
    
    case 4 :
      res[0] = sin2 *((4*cos2*cos2)-3.)*sqrt(r)/2.;
      res[1] = cos2*(5. - (4*cos2*cos2))*sqrt(r)/2. ;
      break;
    case 5 :
      res[0] = cos2*((4*cos2*cos2)-1.)*sqrt(r)/2.;
      res[1] = sin2*((4*cos2*cos2)+1.)*sqrt(r)/2. ;
      break;
   
    case 6 :
      res[0] = sin2*cos2*cos2*sqrt(r)/2.;
      res[1] = cos2*(2. - (cos2*cos2))*sqrt(r)/2.;
      break;
    case 7 :
      res[0] = 3*cos2*cos2*cos2*sqrt(r)/2.;
      res[1] = 3*sin2*cos2*cos2*sqrt(r)/2.;
      break;
    
  /* First order enrichement of pressure field (linear elasticity) mixed formulation */

    case 8 :
      res[0] =sin2*((4*cos2*cos2)-1.)/(2*sqrt(r)*r);
      res[1] =-cos2*((4*cos2*cos2)-3.)/(2*sqrt(r)*r);
      break;
    case 9 :
      res[0] =-cos2*((2*cos2*cos2) - 3.)/(2*sqrt(r)*r);
      res[1] =-sin2*((4*cos2*cos2)-1.)/(2*sqrt(r)*r);
      break;

 /* Second order enrichement of pressure field (linear elasticity) mixed formulation */
    case 10 :
      res[0] = -sin2/(2*sqrt(r));
      res[1] =  cos2/(2*sqrt(r));
      break;
    case 11 :
      res[0] = cos2/(2*sqrt(r));
      res[1] = sin2/(2*sqrt(r));
      break;
    
 /* First order enrichement of displacement field (nonlinear elasticity)[Rodney Stephenson Journal of elasticity VOL.12 No. 1, January 1982] */

    case 12 :
      res[0] = sin2*sin2;
      res[1] = 0.5*cos2*sin2;
      break;
    case 13 :
      res[0] = -sin2/(2*sqrt(r));
      res[1] = cos2/(2*sqrt(r));
      break;

/* First order enrichement of pressure field (****Nonlinear elasticity*****)  */


    case 14 :
      res[0] = -sin2/r;
      res[1] = cos2/(2*r);
      break;
    case 15 :
      res[0] = -cos2/r;
      res[1] = -sin2/(2*r);
      break;


    default: GMM_ASSERT2(false, "oups");
    }
    return res;
  }

  base_matrix crack_singular_xy_function::hess(scalar_type x, scalar_type y) const {
    scalar_type sgny = (y < 0 ? -1.0 : 1.0);
    scalar_type r = sqrt(x*x + y*y);

    if (r < 1e-10) {
      GMM_WARNING0("Warning, point close to the singularity (r=" << r << ")");
    }

    /* ci-dessous: la valeur absolue est malheureusement necessaire,
     * sinon il peut arriver qu'on cherche sqrt(-1e-16) ...
     */
    scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
    scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));

    base_matrix res(2,2);
    switch(l){
    case 0 :
      res(0,0) = (-sin2*x*x + 2.0*cos2*x*y + sin2*y*y) / (4*pow(r, 3.5));
      res(0,1) = (-cos2*x*x - 2.0*sin2*x*y + cos2*y*y) / (4*pow(r, 3.5));
      res(1,0) = res(0, 1);
      res(1,1) = (sin2*x*x - 2.0*cos2*x*y - sin2*y*y) / (4*pow(r, 3.5));
      break;
    case 1 :
      res(0,0) = (-cos2*x*x - 2.0*sin2*x*y + cos2*y*y) / (4*pow(r, 3.5));
      res(0,1) = (sin2*x*x - 2.0*cos2*x*y - sin2*y*y) / (4*pow(r, 3.5));
      res(1,0) = res(0, 1);
      res(1,1) = (cos2*x*x + 2.0*sin2*x*y - cos2*y*y) / (4*pow(r, 3.5));
      break;
    case 2 :
      res(0,0) = 3.0*y*(sin2*x*x + 2.0*cos2*x*y - sin2*y*y) / (4*pow(r, 4.5));
      res(0,1) = (-2.0*sin2*x*x*x - 5.0*cos2*y*x*x + 4.0*sin2*y*y*x + cos2*y*y*y)
        / (4*pow(r, 4.5));
      res(1,0) = res(0, 1);
      res(1,1) = (4.0*cos2*x*x*x - 7.0*sin2*y*x*x - 2.0*cos2*y*y*x - sin2*y*y*y)
        / (4*pow(r, 4.5));
      break;
    case 3 :
      res(0,0) = 3.0*y*(cos2*x*x - 2.0*sin2*x*y - cos2*y*y) / (4*pow(r, 4.5));
      res(0,1) = (-2.0*cos2*x*x*x + 5.0*sin2*y*x*x + 4.0*cos2*y*y*x - sin2*y*y*y)
        / (4*pow(r, 4.5));
      res(1,0) = res(0, 1);
      res(1,1) = (-4.0*sin2*x*x*x - 7.0*cos2*y*x*x + 2.0*sin2*y*y*x - cos2*y*y*y)
        / (4*pow(r, 4.5));
      break;
    default: GMM_ASSERT2(false, "oups");
    }
    return res;
  }

  scalar_type
  cutoff_xy_function::val(scalar_type x, scalar_type y) const {
    scalar_type res = 1;
    switch (fun) {
      case EXPONENTIAL_CUTOFF: {
        res = (a4>0) ? exp(-a4 * gmm::sqr(x*x+y*y)) : 1;
        break;
      }
      case POLYNOMIAL_CUTOFF: {
        assert(r0 > r1);
        scalar_type r = gmm::sqrt(x*x+y*y);

        if (r <= r1) res = 1;
        else if (r >= r0) res = 0;
        else {
          // scalar_type c = 6./(r0*r0*r0 - r1*r1*r1 + 3*r1*r0*(r1-r0));
          // scalar_type k = -(c/6.)*(-pow(r0,3.) + 3*r1*pow(r0,2.));
          // res = (c/3.)*pow(r,3.) - (c*(r0 + r1)/2.)*pow(r,2.) +
          //        c*r0*r1*r + k;
          scalar_type c = 1./pow(r0-r1,3.0);
          res = c*(r*(r*(2.0*r-3.0*(r0+r1))+6.0*r1*r0) + r0*r0*(r0-3.0*r1));
        }
        break;
      }
      case POLYNOMIAL2_CUTOFF: {
        assert(r0 > r1);
        scalar_type r = gmm::sqrt(x*x+y*y);
        if (r <= r1) res = scalar_type(1);
        else if (r >= r0) res = scalar_type(0);
        else {
//           scalar_type c = 1./((-1./30.)*(pow(r1,5) - pow(r0,5))
//                               + (1./6.)*(pow(r1,4)*r0 - r1*pow(r0,4))
//                               - (1./3.)*(pow(r1,3)*pow(r0,2) -
//                                          pow(r1,2)*pow(r0,3)));
//           scalar_type k = 1. - c*((-1./30.)*pow(r1,5) +
//                                   (1./6.)*pow(r1,4)*r0 -
//                                   (1./3.)*pow(r1,3)*pow(r0,2));
//           res = c*( (-1./5.)*pow(r,5) + (1./2.)* (r1+r0)*pow(r,4) -
//                      (1./3.)*(pow(r1,2)+pow(r0,2) + 4.*r0*r1)*pow(r,3) +
//                      r0*r1*(r0+r1)*pow(r,2) - pow(r0,2)*pow(r1,2)*r) + k;
          res = (r*(r*(r*(r*(-6.0*r + 15.0*(r0+r1))
                           - 10.0*(r0*r0 + 4.0*r1*r0 + r1*r1))
                        + 30.0 * r0*r1*(r0+r1)) - 30.0*r1*r1*r0*r0)
                  + r0*r0*r0*(r0*r0-5.0*r1*r0+10*r1*r1)) / pow(r0-r1, 5.0);
        }
        break;
      }
      default : res = 1;
    }
    return res;
  }

  base_small_vector
  cutoff_xy_function::grad(scalar_type x, scalar_type y) const {
    base_small_vector res(2);
    switch (fun) {
    case EXPONENTIAL_CUTOFF: {
      scalar_type r2 = x*x+y*y, ratio = -4.*exp(-a4*r2*r2)*a4*r2;
      res[0] = ratio*x;
      res[1] = ratio*y;
      break;
    }
    case POLYNOMIAL_CUTOFF: {
      scalar_type r = gmm::sqrt(x*x+y*y);
      scalar_type ratio = 0;

      if ( r > r1 && r < r0 ) {
        // scalar_type c = 6./(pow(r0,3.) - pow(r1,3.) + 3*r1*r0*(r1-r0));
        // ratio = c*(r - r0)*(r - r1);
        ratio = 6.*(r - r0)*(r - r1)/pow(r0-r1, 3.);
      }
      res[0] = ratio*x/r;
      res[1] = ratio*y/r;
      break;
    }
    case POLYNOMIAL2_CUTOFF: {
      scalar_type r = gmm::sqrt(x*x+y*y);
      scalar_type ratio = 0;
      if (r > r1 && r < r0) {
//        scalar_type c = 1./((-1./30.)*(pow(r1,5) - pow(r0,5))
//                            + (1./6.)*(pow(r1,4)*r0 - r1*pow(r0,4))
//                            - (1./3.)*(pow(r1,3)*pow(r0,2)
//                            - pow(r1,2)*pow(r0,3)));
//        ratio = - c*gmm::sqr(r-r0)*gmm::sqr(r-r1);
        ratio = -30.0*gmm::sqr(r-r0)*gmm::sqr(r-r1) / pow(r0-r1, 5.0);
      }
      res[0] = ratio*x/r;
      res[1] = ratio*y/r;
      break;
    }
    default :
      res[0] = 0;
      res[1] = 0;
    }
    return res;
  }

  base_matrix
  cutoff_xy_function::hess(scalar_type x, scalar_type y) const {
    base_matrix res(2,2);
    switch (fun) {
    case EXPONENTIAL_CUTOFF: {
      scalar_type r2 = x*x+y*y, r4 = r2*r2;
      res(0,0) = 4.0*a4*(-3.0*x*x - y*y + 4.0*a4*x*x*r4)*exp(-a4*r4);
      res(1,0) = 8.0*a4*x*y*(-1.0 + 2.0*a4*r4)*exp(-a4*r4);
      res(0,1) = res(1,0);
      res(1,1) = 4.0*a4*(-3.0*y*y - x*x + 4.0*a4*y*y*r4)*exp(-a4*r4);
      break;
    }
    case POLYNOMIAL_CUTOFF: {
      scalar_type r2 = x*x+y*y, r = gmm::sqrt(r2), c=6./(pow(r0-r1,3.)*r*r2);
      if ( r > r1 && r < r0 ) {
        res(0,0) = c*(x*x*r2 + r1*r0*y*y - r*r2*(r0+r1) + r2*r2);
        res(1,0) = c*x*y*(r2 - r1*r0);
        res(0,1) = res(1,0);
        res(1,1) = c*(y*y*r2 + r1*r0*x*x - r*r2*(r0+r1) + r2*r2);
      }
      break;
    }
    case POLYNOMIAL2_CUTOFF: {
      scalar_type r2 = x*x+y*y, r = gmm::sqrt(r2), r3 = r*r2;
      if (r > r1 && r < r0) {
        scalar_type dp = -30.0*(r1-r)*(r1-r)*(r0-r)*(r0-r) / pow(r0-r1, 5.0);
        scalar_type ddp = 60.0*(r1-r)*(r0-r)*(r0+r1-2.0*r) / pow(r0-r1, 5.0);
        scalar_type rx= x/r, ry= y/r, rxx= y*y/r3, rxy= -x*y/r3, ryy= x*x/r3;
        res(0,0) = ddp*rx*rx + dp*rxx;
        res(1,0) = ddp*rx*ry + dp*rxy;
        res(0,1) = res(1,0);
        res(1,1) = ddp*ry*ry + dp*ryy;
      }
      break;
    }
    }
    return res;
  }

  cutoff_xy_function::cutoff_xy_function(int fun_num, scalar_type r,
                                         scalar_type r1_, scalar_type r0_)
  {
    fun = fun_num;
    r1 = r1_; r0 = r0_;
    a4 = 0;
    if (r > 0) a4 = pow(2.7/r,4.);
  }


  struct global_function_on_levelset_ :
    public global_function, public context_dependencies {
    const level_set &ls;
    mutable pmesher_signed_distance mls_x, mls_y;
    mutable size_type cv;

    pxy_function fn;

    void update_mls(size_type cv_) const {
      if (cv_ != cv) {
        cv=cv_;
	mls_x = ls.mls_of_convex(cv, 1);
        mls_y = ls.mls_of_convex(cv, 0);
      }
    }

    virtual scalar_type val(const fem_interpolation_context& c) const {
      update_mls(c.convex_num());
      scalar_type x = (*mls_x)(c.xref());
      scalar_type y = (*mls_y)(c.xref());
      if (c.xfem_side() > 0 && y <= 0) y = 1E-13;
      if (c.xfem_side() < 0 && y >= 0) y = -1E-13;
      return fn->val(x,y);
    }
    virtual void grad(const fem_interpolation_context& c,
                      base_small_vector &g) const {
      update_mls(c.convex_num());
      size_type P = c.xref().size();
      base_small_vector dx(P), dy(P), dfr(2);
      scalar_type x = mls_x->grad(c.xref(), dx);
      scalar_type y = mls_y->grad(c.xref(), dy);
      if (c.xfem_side() > 0 && y <= 0) y = 1E-13;
      if (c.xfem_side() < 0 && y >= 0) y = -1E-13;
      base_small_vector gfn = fn->grad(x,y);
      gmm::mult(c.B(), gfn[0]*dx + gfn[1]*dy, g);
    }
    virtual void hess(const fem_interpolation_context&c,
                      base_matrix &h) const {
      update_mls(c.convex_num());
      size_type P = c.xref().size(), N = c.N();

      base_small_vector dx(P), dy(P), dfr(2),  dx_real(N), dy_real(N);
      scalar_type x = mls_x->grad(c.xref(), dx);
      scalar_type y = mls_y->grad(c.xref(), dy);
      if (c.xfem_side() > 0 && y <= 0) y = 1E-13;
      if (c.xfem_side() < 0 && y >= 0) y = -1E-13;
      base_small_vector gfn = fn->grad(x,y);
      base_matrix hfn = fn->hess(x,y);

      base_matrix hx, hy, hx_real(N*N, 1), hy_real(N*N, 1);
      mls_x->hess(c.xref(), hx);
      mls_x->hess(c.xref(), hy);
      gmm::reshape(hx, P*P, 1);
      gmm::reshape(hy, P*P, 1);

      gmm::mult(c.B3(), hx, hx_real);
      gmm::mult(c.B32(), gmm::scaled(dx, -1.0), gmm::mat_col(hx_real, 0));
      gmm::mult(c.B3(), hy, hy_real);
      gmm::mult(c.B32(), gmm::scaled(dy, -1.0), gmm::mat_col(hy_real, 0));
      gmm::mult(c.B(), dx, dx_real);
      gmm::mult(c.B(), dy, dy_real);

      for (size_type i = 0; i < N; ++i)
        for (size_type j = 0; j < N; ++j) {
          h(i, j) = hfn(0,0) * dx_real[i] * dx_real[j]
            + hfn(0,1) * dx_real[i] * dy_real[j]
            + hfn(1,0) * dy_real[i] * dx_real[j]
            + hfn(1,1) * dy_real[i] * dy_real[j]
            + gfn[0] * hx_real(i * N + j, 0) + gfn[1] * hy_real(i*N+j,0);
        }
    }

    void update_from_context(void) const { cv =  size_type(-1); }

    global_function_on_levelset_(const level_set &ls_, const pxy_function &fn_)
      : ls(ls_), fn(fn_) {
      cv = size_type(-1);
      this->add_dependency(ls);
    }

  };

  pglobal_function
  global_function_on_level_set(const level_set &ls,
                               const pxy_function &fn) {
    return std::make_shared<global_function_on_levelset_>(ls, fn);
  }


  struct global_function_on_levelsets_ :
    public global_function, public context_dependencies {
    const std::vector<level_set> &lsets;
    mutable pmesher_signed_distance mls_x, mls_y;
    mutable size_type cv;

    pxy_function fn;

    void update_mls(size_type cv_, size_type n) const {
      if (cv_ != cv) {
        base_node pt(n);
        cv=cv_;
        scalar_type d = scalar_type(-2);
        for (size_type i = 0; i < lsets.size(); ++i) {
          pmesher_signed_distance mls_xx, mls_yy;
          mls_xx = lsets[i].mls_of_convex(cv, 1);
          mls_yy = lsets[i].mls_of_convex(cv, 0);
          scalar_type x = (*mls_xx)(pt), y = (*mls_yy)(pt);
          scalar_type d2 = gmm::sqr(x) + gmm::sqr(y);
          if (d < scalar_type(-1) || d2 < d) {
            d = d2;
            mls_x = mls_xx; mls_y = mls_yy;
          }
        }
      }
    }

    virtual scalar_type val(const fem_interpolation_context& c) const {
      update_mls(c.convex_num(), c.xref().size());
      scalar_type x = (*mls_x)(c.xref());
      scalar_type y = (*mls_y)(c.xref());
      return fn->val(x,y);
    }
    virtual void grad(const fem_interpolation_context& c,
                      base_small_vector &g) const {
      size_type P = c.xref().size();
      update_mls(c.convex_num(), P);
      base_small_vector dx(P), dy(P), dfr(2);
      scalar_type x = mls_x->grad(c.xref(), dx);
      scalar_type y = mls_y->grad(c.xref(), dy);

      base_small_vector gfn = fn->grad(x,y);

      gmm::mult(c.B(), gfn[0]*dx + gfn[1]*dy, g);
    }
    virtual void hess(const fem_interpolation_context&c,
                      base_matrix &h) const {
      size_type P = c.xref().size(), N = c.N();
      update_mls(c.convex_num(), P);

      base_small_vector dx(P), dy(P), dfr(2),  dx_real(N), dy_real(N);
      scalar_type x = mls_x->grad(c.xref(), dx);
      scalar_type y = mls_y->grad(c.xref(), dy);

      base_small_vector gfn = fn->grad(x,y);
      base_matrix hfn = fn->hess(x,y);

      base_matrix hx, hy, hx_real(N*N, 1), hy_real(N*N, 1);
      mls_x->hess(c.xref(), hx);
      mls_x->hess(c.xref(), hy);
      gmm::reshape(hx, P*P, 1);
      gmm::reshape(hy, P*P, 1);

      gmm::mult(c.B3(), hx, hx_real);
      gmm::mult(c.B32(), gmm::scaled(dx, -1.0), gmm::mat_col(hx_real, 0));
      gmm::mult(c.B3(), hy, hy_real);
      gmm::mult(c.B32(), gmm::scaled(dy, -1.0), gmm::mat_col(hy_real, 0));
      gmm::mult(c.B(), dx, dx_real);
      gmm::mult(c.B(), dy, dy_real);

      for (size_type i = 0; i < N; ++i)
        for (size_type j = 0; j < N; ++j) {
          h(i, j) = hfn(0,0) * dx_real[i] * dx_real[j]
            + hfn(0,1) * dx_real[i] * dy_real[j]
            + hfn(1,0) * dy_real[i] * dx_real[j]
            + hfn(1,1) * dy_real[i] * dy_real[j]
            + gfn[0] * hx_real(i * N + j, 0) + gfn[1] * hy_real(i*N+j,0);
        }
    }

    void update_from_context(void) const { cv =  size_type(-1); }

    global_function_on_levelsets_(const std::vector<level_set> &lsets_,
                                  const pxy_function &fn_)
      : lsets(lsets_), fn(fn_) {
      cv = size_type(-1);
      for (size_type i = 0; i < lsets.size(); ++i)
        this->add_dependency(lsets[i]);
    }

  };

  pglobal_function
  global_function_on_level_sets(const std::vector<level_set> &lsets,
                                const pxy_function &fn) {
    return std::make_shared<global_function_on_levelsets_>(lsets, fn);
  }


  void interpolator_on_mesh_fem::init() {
    base_node min, max;
    scalar_type EPS=1E-13;
    cv_stored = size_type(-1);
    boxtree.clear();
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bounding_box(min, max, mf.linked_mesh().points_of_convex(cv),
                   mf.linked_mesh().trans_of_convex(cv));
      for (unsigned k=0; k < min.size(); ++k) { min[k]-=EPS; max[k]+=EPS; }
      boxtree.add_box(min, max, cv);
    }
  }

  bool interpolator_on_mesh_fem::find_a_point(base_node pt, base_node &ptr,
                                              size_type &cv) const {
    bool gt_invertible;
    if (cv_stored != size_type(-1) && gic.invert(pt, ptr, gt_invertible))
      { cv = cv_stored; if (gt_invertible) return true; }

    boxtree.find_boxes_at_point(pt, boxlst);
    bgeot::rtree::pbox_set::const_iterator it = boxlst.begin(),
      ite = boxlst.end();
    for (; it != ite; ++it) {
      gic = bgeot::geotrans_inv_convex
        (mf.linked_mesh().convex((*it)->id),
         mf.linked_mesh().trans_of_convex((*it)->id));
      cv_stored = (*it)->id;
      if (gic.invert(pt, ptr, gt_invertible)) {
        cv = (*it)->id;
        return true;
      }
    }
    return false;
  }

  bool interpolator_on_mesh_fem::eval(const base_node pt, base_vector &val,
                                      base_matrix &grad) const {
    base_node ptref;
    size_type cv;
    base_vector coeff;
    dim_type q = mf.get_qdim(), N = mf.linked_mesh().dim();
    if (find_a_point(pt, ptref, cv)) {
      pfem pf = mf.fem_of_element(cv);
      bgeot::pgeometric_trans pgt =
        mf.linked_mesh().trans_of_convex(cv);
      base_matrix G;
      vectors_to_base_matrix(G, mf.linked_mesh().points_of_convex(cv));
      fem_interpolation_context fic(pgt, pf, ptref, G, cv, short_type(-1));
      slice_vector_on_basic_dof_of_element(mf, U, cv, coeff);
      // coeff.resize(mf.nb_basic_dof_of_element(cv));
      // gmm::copy(gmm::sub_vector
      //          (U,gmm::sub_index(mf.ind_basic_dof_of_element(cv))), coeff);
      val.resize(q);
      pf->interpolation(fic, coeff, val, q);
      grad.resize(q, N);
      pf->interpolation_grad(fic, coeff, grad, q);
      return true;
    } else return false;
  }
}

/* end of namespace getfem  */
