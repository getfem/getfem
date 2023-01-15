/*===========================================================================

 Copyright (C) 2004-2022 Yves Renard
 Copyright (C) 2016-2022 Konstantinos Poulios

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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

#include <getfem/getfem_global_function.h>
#include <getfem/getfem_level_set.h>

namespace getfem {


  // Partial implementation of abstract class global_function_simple

  scalar_type global_function_simple::val
  (const fem_interpolation_context &c) const {
    base_node pt = c.xreal();
    GMM_ASSERT1(pt.size() == dim_, "Point of wrong size (" << pt.size() << ") "
		<< "passed to a global function of dim = "<< dim_ <<".");
    return this->val(pt);
  }

  void global_function_simple::grad
  (const fem_interpolation_context &c, base_small_vector &g) const {
    base_node pt = c.xreal();
    GMM_ASSERT1(pt.size() == dim_, "Point of wrong size (" << pt.size() << ") "
		<< "passed to a global function of dim = "<< dim_ <<".");
    this->grad(pt, g);
  }

  void global_function_simple::hess
  (const fem_interpolation_context &c, base_matrix &h) const {
    base_node pt = c.xreal();
    GMM_ASSERT1(pt.size() == dim_, "Point of wrong size (" << pt.size() << ") "
		<< "passed to a global function of dim = "<< dim_ <<".");
    this->hess(pt, h);
  }

  // Implementation of global_function_parser

  scalar_type global_function_parser::val(const base_node &pt) const {
    const bgeot::base_tensor &t = tensor_val(pt);
    GMM_ASSERT1(t.size() == 1, "Wrong size of expression result "
                << f_val.expression());
    return t[0];
  }

  const base_tensor &global_function_parser::tensor_val(const base_node &pt) const {
    gmm::copy(pt, pt_);
    return f_val.eval();
  }

  void global_function_parser::grad(const base_node &pt, base_small_vector &g) const {
    g.resize(dim_);
    gmm::copy(pt, pt_);
    const bgeot::base_tensor &t = f_grad.eval();
    GMM_ASSERT1(t.size() == dim_, "Wrong size of expression result "
                << f_grad.expression());
    gmm::copy(t.as_vector(), g);

  }

  void global_function_parser::hess(const base_node &pt, base_matrix &h) const {
    h.resize(dim_, dim_);
    gmm::copy(pt, pt_);
    const bgeot::base_tensor &t = f_hess.eval();
    GMM_ASSERT1(t.size() == size_type(dim_*dim_),
                "Wrong size of expression result " << f_hess.expression());
    gmm::copy(t.as_vector(), h.as_vector());
  }

  global_function_parser::global_function_parser(dim_type dim__,
                                                 const std::string &sval,
                                                 const std::string &sgrad,
                                                 const std::string &shess)
    : global_function_simple(dim__),
      f_val(gw, sval), f_grad(gw, sgrad), f_hess(gw, shess) {

    size_type N(dim_);
    pt_.resize(N);
    gmm::fill(pt_, scalar_type(0));
    gw.add_fixed_size_variable("X", gmm::sub_interval(0, N), pt_);
    if (N >= 1) gw.add_macro("x", "X(1)");
    if (N >= 2) gw.add_macro("y", "X(2)");
    if (N >= 3) gw.add_macro("z", "X(3)");
    if (N >= 4) gw.add_macro("w", "X(4)");
    
    f_val.compile();
    f_grad.compile();
    f_hess.compile();
  }

  // Implementation of global_function_sum

  scalar_type global_function_sum::val
  (const fem_interpolation_context &c) const {
    scalar_type res(0);
    for (const auto &f : functions)
      res += f->val(c);
    return res;
  }

  void global_function_sum::grad
  (const fem_interpolation_context &c, base_small_vector &g) const {
    g.resize(dim_);
    gmm::clear(g);
    base_small_vector gg(dim_);
    for (const auto &f : functions) {
      f->grad(c, gg);
      gmm::add(gg, g);
    }
  }

  void global_function_sum::hess
  (const fem_interpolation_context &c, base_matrix &h) const {
    h.resize(dim_, dim_);
    gmm::clear(h);
    base_matrix hh(dim_, dim_);
    for (const auto &f : functions) {
      f->hess(c, hh);
      gmm::add(hh.as_vector(), h.as_vector());
    }
  }

  bool global_function_sum::is_in_support(const base_node &p) const {
    for (const auto &f : functions)
      if (f->is_in_support(p)) return true;
    return false;
  }

  void global_function_sum::bounding_box
  (base_node &bmin_, base_node &bmax_) const {
    if (functions.size() > 0)
      functions[0]->bounding_box(bmin_, bmax_);
    base_node bmin0(dim()), bmax0(dim());
    for (const auto &f : functions) {
      f->bounding_box(bmin0, bmax0);
      for (size_type i=0; i < dim(); ++i) {
        if (bmin0[i] < bmin_[i]) bmin_[i] = bmin0[i];
        if (bmax0[i] > bmax_[i]) bmax_[i] = bmax0[i];
      }
    }
  }

  global_function_sum::global_function_sum(const std::vector<pglobal_function> &funcs)
    : global_function((funcs.size() > 0) ? funcs[0]->dim() : 0), functions(funcs) {
    for (const auto &f : functions)
      GMM_ASSERT1(f->dim() == dim(), "Incompatible dimensions among the provided"
                                     " global functions");
  }

  global_function_sum::global_function_sum(pglobal_function f1, pglobal_function f2)
    : global_function(f1->dim()), functions(2) {
    functions[0] = f1;
    functions[1] = f2;
    GMM_ASSERT1(f1->dim() == dim() && f2->dim() == dim(),
                "Incompatible dimensions between the provided global functions");
  }

  global_function_sum::global_function_sum(pglobal_function f1, pglobal_function f2,
                                           pglobal_function f3)
    : global_function(f1->dim()), functions(3) {
    functions[0] = f1;
    functions[1] = f2;
    functions[2] = f3;
    GMM_ASSERT1(f1->dim() == dim() && f2->dim() == dim() && f3->dim() == dim(),
                "Incompatible dimensions between the provided global functions");
  }

  global_function_sum::global_function_sum(pglobal_function f1, pglobal_function f2,
                                           pglobal_function f3, pglobal_function f4)
    : global_function(f1->dim()), functions(4) {
    functions[0] = f1;
    functions[1] = f2;
    functions[2] = f3;
    functions[3] = f4;
    GMM_ASSERT1(f1->dim() == dim() && f2->dim() == dim() && f3->dim() == dim(),
                "Incompatible dimensions between the provided global functions");
  }


  // Implementation of global_function_product

  scalar_type global_function_product::val
  (const fem_interpolation_context &c) const {
    return f1->val(c) * f2->val(c);
  }

  void global_function_product::grad
  (const fem_interpolation_context &c, base_small_vector &g) const {
    g.resize(dim_);
    base_small_vector gg(dim_);
    f1->grad(c, gg);
    gmm::copy(gmm::scaled(gg, f2->val(c)), g);
    f2->grad(c, gg);
    gmm::add(gmm::scaled(gg, f1->val(c)), g);
  }

  void global_function_product::hess
  (const fem_interpolation_context &c, base_matrix &h) const {
    h.resize(dim_, dim_);
    gmm::clear(h);
    base_matrix hh(dim_, dim_);
    f1->hess(c, hh);
    gmm::copy(gmm::scaled(hh, f2->val(c)), h);
    f2->hess(c, hh);
    gmm::add(gmm::scaled(hh, f1->val(c)), h);
    base_small_vector g1(dim_), g2(dim_);
    f1->grad(c, g1);
    f2->grad(c, g2);
    gmm::rank_one_update(h, g1, g2);
    gmm::rank_one_update(h, g2, g1);
  }

  bool global_function_product::is_in_support(const base_node &p) const {
    return f1->is_in_support(p) && f2->is_in_support(p);
  }

  void global_function_product::bounding_box
  (base_node &bmin_, base_node &bmax_) const {
    base_node bmin0(dim()), bmax0(dim());
    f1->bounding_box(bmin0, bmax0);
    f2->bounding_box(bmin_, bmax_);
    for (size_type i=0; i < dim(); ++i) {
      if (bmin0[i] > bmin_[i]) bmin_[i] = bmin0[i];
      if (bmax0[i] < bmax_[i]) bmax_[i] = bmax0[i];
      if (bmin_[i] > bmax_[i])
        GMM_WARNING1("Global function product with vanishing basis function");
    }
  }

  global_function_product::global_function_product(pglobal_function f1_, pglobal_function f2_)
    : global_function(f1_->dim()), f1(f1_), f2(f2_) {
    GMM_ASSERT1(f2->dim() == dim(), "Incompatible dimensions between the provided"
                                    " global functions");
  }


  // Implementation of global_function_bounded

  bool global_function_bounded::is_in_support(const base_node &pt) const {
    if (has_expr) {
      gmm::copy(pt, pt_);
      const bgeot::base_tensor &t = f_val.eval();
      GMM_ASSERT1(t.size() == 1, "Wrong size of expression result "
                  << f_val.expression());
      return (t[0] > scalar_type(0));
    }
    return true;
  }

  global_function_bounded::global_function_bounded(pglobal_function f_,
                                                   const base_node &bmin_,
                                                   const base_node &bmax_,
                                                   const std::string &is_in_expr)
    : global_function(f_->dim()), f(f_), bmin(bmin_), bmax(bmax_),
      f_val(gw, is_in_expr) {

    has_expr = !is_in_expr.empty();
    if (has_expr) {
      size_type N(dim_);
      pt_.resize(N);
      gmm::fill(pt_, scalar_type(0));
      gw.add_fixed_size_variable("X", gmm::sub_interval(0, N), pt_);
      if (N >= 1) gw.add_macro("x", "X(1)");
      if (N >= 2) gw.add_macro("y", "X(2)");
      if (N >= 3) gw.add_macro("z", "X(3)");
      if (N >= 4) gw.add_macro("w", "X(4)");
      f_val.compile();
    }
  }

  // Implementation of some useful xy functions

  parser_xy_function::parser_xy_function(const std::string &sval,
                                         const std::string &sgrad,
                                         const std::string &shess)
    : f_val(gw, sval), f_grad(gw, sgrad), f_hess(gw, shess),
      ptx(1), pty(1), ptr(1), ptt(1) {

    gw.add_fixed_size_constant("x", ptx);
    gw.add_fixed_size_constant("y", pty);
    gw.add_fixed_size_constant("r", ptr);
    gw.add_fixed_size_constant("theta", ptt);

    f_val.compile();
    f_grad.compile();
    f_hess.compile();
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

  /* the basic singular functions for 2D cracks */
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

    /* The absolute value is unfortunately necessary, otherwise, sqrt(-1e-16)
       can be required ...
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

  base_matrix
  crack_singular_xy_function::hess(scalar_type x, scalar_type y) const {
    scalar_type sgny = (y < 0 ? -1.0 : 1.0);
    scalar_type r = sqrt(x*x + y*y);

    if (r < 1e-10) {
      GMM_WARNING0("Warning, point close to the singularity (r=" << r << ")");
    }

    /* The absolute value is unfortunately necessary, otherwise, sqrt(-1e-16)
       can be required ...
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


  struct global_function_on_levelsets_2D_ :
    public global_function, public context_dependencies {
    const std::vector<level_set> dummy_lsets;
    const std::vector<level_set> &lsets;
    const level_set &ls;
    mutable pmesher_signed_distance mls_x, mls_y;
    mutable size_type cv;

    pxy_function fn;

    void update_mls(size_type cv_, size_type n) const {
      if (cv_ != cv) {
        cv=cv_;
        if (lsets.size() == 0) {
	  mls_x = ls.mls_of_convex(cv, 1);
	  mls_y = ls.mls_of_convex(cv, 0);
        } else {
	  base_node pt(n);
          scalar_type d = scalar_type(-2);
          for (const auto &ls_ : lsets) {
            pmesher_signed_distance mls_xx, mls_yy;
            mls_xx = ls_.mls_of_convex(cv, 1);
            mls_yy = ls_.mls_of_convex(cv, 0);
            scalar_type x = (*mls_xx)(pt), y = (*mls_yy)(pt);
            scalar_type d2 = gmm::sqr(x) + gmm::sqr(y);
            if (d < scalar_type(-1) || d2 < d) {
              d = d2;
              mls_x = mls_xx;
              mls_y = mls_yy;
           }
          }
        }
      }
    }

    virtual scalar_type val(const fem_interpolation_context& c) const {
      update_mls(c.convex_num(), c.xref().size());
      scalar_type x = (*mls_x)(c.xref());
      scalar_type y = (*mls_y)(c.xref());
      if (c.xfem_side() > 0 && y <= 1E-13) y = 1E-13;
      if (c.xfem_side() < 0 && y >= -1E-13) y = -1E-13;
      return fn->val(x,y);
    }
    virtual void grad(const fem_interpolation_context& c,
                      base_small_vector &g) const {
      size_type P = c.xref().size();
      base_small_vector dx(P), dy(P), dfr(2);

      update_mls(c.convex_num(), P);
      scalar_type x = mls_x->grad(c.xref(), dx);
      scalar_type y = mls_y->grad(c.xref(), dy);
      if (c.xfem_side() > 0 && y <= 0) y = 1E-13;
      if (c.xfem_side() < 0 && y >= 0) y = -1E-13;

      base_small_vector gfn = fn->grad(x,y);
      gmm::mult(c.B(), gfn[0]*dx + gfn[1]*dy, g);
    }
    virtual void hess(const fem_interpolation_context&c,
                      base_matrix &h) const {
      size_type P = c.xref().size(), N = c.N();
      base_small_vector dx(P), dy(P), dfr(2),  dx_real(N), dy_real(N);

      update_mls(c.convex_num(), P);
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

    void update_from_context() const { cv =  size_type(-1); }

    global_function_on_levelsets_2D_(const std::vector<level_set> &lsets_,
                                     const pxy_function &fn_)
      : global_function(2), dummy_lsets(0, dummy_level_set()),
        lsets(lsets_), ls(dummy_level_set()), fn(fn_) {
      GMM_ASSERT1(lsets.size() > 0, "The list of level sets should"
                                    " contain at least one level set.");
      cv = size_type(-1);
      for (size_type i = 0; i < lsets.size(); ++i)
        this->add_dependency(lsets[i]);
    }

    global_function_on_levelsets_2D_(const level_set &ls_,
                                     const pxy_function &fn_)
      : global_function(2), dummy_lsets(0, dummy_level_set()),
        lsets(dummy_lsets), ls(ls_), fn(fn_) {
      cv = size_type(-1);
    }

  };

  pglobal_function
  global_function_on_level_sets(const std::vector<level_set> &lsets,
                                const pxy_function &fn) {
    return std::make_shared<global_function_on_levelsets_2D_>(lsets, fn);
  }

  pglobal_function
  global_function_on_level_set(const level_set &ls,
                               const pxy_function &fn) {
    return std::make_shared<global_function_on_levelsets_2D_>(ls, fn);
  }




  // Global function for bspline basis
  const scalar_type eps(1e-12);

  // order k = 3
  scalar_type bsp3_01(scalar_type t) {
    return (t >= -eps && t < 1) ? pow(1.-t,2)
                                : 0;
  }
  scalar_type bsp3_01_der(scalar_type t) {
    return (t >= -eps && t < 1) ? 2.*t-2.
                                : 0;
  }
  scalar_type bsp3_01_der2(scalar_type t) {
    return (t >= eps && t <= 1.-eps) ? 2.
                                     : 0;
  }
  scalar_type bsp3_01_der2_with_hint(scalar_type t, scalar_type t_hint) {
    return (t > -eps && t < 1.+eps && t_hint > 0 && t_hint < 1) ? 2.
                                                                : 0;
  }
  scalar_type bsp3_02(scalar_type t) {
    if (t >= -eps) {
      if (t < 1)
        return (-1.5*t+2.)*t;
      else if (t < 2)
        return 0.5*pow(2.-t,2);
    }
    return 0;
  }
  scalar_type bsp3_02_der(scalar_type t) {
    if (t >= -eps) {
      if (t < 1)
        return -3.*t+2.;
      else if (t < 2)
        return t-2.;
    }
    return 0;
  }
  scalar_type bsp3_02_der2(scalar_type t) {
    if (t >= eps) {
      if (t < 1.-eps)
        return -3.;
      else if (t < 1.+eps)
        return 0;
      else if (t <= 2.-eps)
        return 1.;
    }
    return 0;
  }
  scalar_type bsp3_03(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return 0.5*t*t;
      } else if (t < 2) {
        t -= 1.;
        return t*(1-t)+0.5;
      } else if (t < 3) {
        t = 3.-t;
        return 0.5*t*t;
      }
    }
    return 0;
  }
  scalar_type bsp3_03_der(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return t;
      } else if (t < 2) {
        t -= 1.;
        return 1.-2.*t;
      } else if (t < 3) {
        return t-3.;
      }
    }
    return 0;
  }
  scalar_type bsp3_03_der2(scalar_type t) {
    if (t >= -eps) {
      if (t < eps)
        return 0;
      else if (t <= 1.-eps)
        return 1.;
      else if (t < 1.+eps)
        return 0;
      else if (t <= 2.-eps)
        return -2.;
      else if (t < 2.+eps)
        return 0;
      else if (t <= 3.-eps)
        return 1.;
      else if (t < 3.+eps)
        return 0;
    }
    return 0;
  }

  // order k = 4
  scalar_type bsp4_01(scalar_type t) {
    return (t >= -eps && t < 1) ? pow(1.-t,3)
                                : 0;
  }
  scalar_type bsp4_01_der(scalar_type t) {
    return (t >= -eps && t < 1) ? -3*pow(1.-t,2)
                                : 0;
  }
  scalar_type bsp4_01_der2(scalar_type t) {
    return (t >= -eps && t < 1) ? 6*(1.-t)
                                : 0;
  }
  scalar_type bsp4_02(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return ((7./4.*t-9./2.)*t+3.)*t;
      } else if (t < 2) {
        return 1./4.*pow(2.-t,3);
      }
    }
    return 0;
  }
  scalar_type bsp4_02_der(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (21./4.*t-9.)*t+3.;
      } else if (t < 2) {
        return -3./4.*pow(2.-t,2);
      }
    }
    return 0;
  }
  scalar_type bsp4_02_der2(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return 21./2.*t-9.;
      } else if (t < 2) {
        return 3./2.*(2.-t);
      }
    }
    return 0;
  }
  scalar_type bsp4_03(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (-11./12.*t+3./2.)*t*t;
      } else if (t < 2) {
        t -= 1;
        return ((7./12.*t-5./4.)*t+1./4.)*t+7./12.;
      } else if (t < 3) {
        return 1./6.*pow(3.-t,3);
      }
    }
    return 0;
  }
  scalar_type bsp4_03_der(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (-11./4.*t+3.)*t;
      } else if (t < 2) {
        t -= 1;
        return (7./4.*t-5./2.)*t+1./4.;
      } else if (t < 3) {
        return -1./2.*pow(3.-t,2);
      }
    }
    return 0;
  }
  scalar_type bsp4_03_der2(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return -11./2.*t+3.;
      } else if (t < 2) {
        t -= 1;
        return 7./2.*t-5./2.;
      } else if (t < 3) {
        return 3.-t;
      }
    }
    return 0;
  }
  scalar_type bsp4_04(scalar_type t) {
    if (t > 2) {
      t = 4.-t;
    }
    if (t >= -eps) {
      if (t < 1) {
        return 1./6.*pow(t,3);
      } else if (t <= 2) {
        t = t-1.;
        return ((-1./2.*t+1./2.)*t+1./2.)*t+1./6.;
      }
    }
    return 0;
  }
  scalar_type bsp4_04_der(scalar_type t) {
    scalar_type sgn(1);
    if (t > 2) {
      t = 4.-t;
      sgn = -1.;
    }
    if (t >= -eps) {
      if (t < 1) {
        return 1./2.*pow(t,2)*sgn;
      } else if (t < 2) {
        t = t-1.;
        return ((-3./2.*t+1.)*t+1./2.)*sgn;
      }
    }
    return 0;
  }
  scalar_type bsp4_04_der2(scalar_type t) {
    if (t > 2) {
      t = 4.-t;
    }
    if (t >= -eps) {
      if (t < 1) {
        return t;
      } else if (t < 2) {
        t = t-1.;
        return -3.*t+1.;
      }
    }
    return 0;
  }


  // order k = 5
  scalar_type bsp5_01(scalar_type t) {
    return (t >= -eps && t < 1) ? pow(1.-t,4)
                                : 0;
  }
  scalar_type bsp5_01_der(scalar_type t) {
    return (t >= -eps && t < 1) ? -4.*pow(1.-t,3)
                                : 0;
  }
  scalar_type bsp5_01_der2(scalar_type t) {
    return (t >= -eps && t < 1) ? 12.*pow(1.-t,2)
                                : 0;
  }
  scalar_type bsp5_02(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (((-15./8.*t+7.)*t-9.)*t+4.)*t;
      } else if (t < 2) {
        return 1./8.*pow(2.-t,4);
      }
    }
    return 0;
  }
  scalar_type bsp5_02_der(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return ((-15./2.*t+21.)*t-18.)*t+4.;
      } else if (t < 2) {
        return -1./2.*pow(2.-t,3);
      }
    }
    return 0;
  }
  scalar_type bsp5_02_der2(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (-45./2.*t+42.)*t-18.;
      } else if (t < 2) {
        return 3./2.*pow(2.-t,2);
      }
    }
    return 0;
  }
  scalar_type bsp5_03(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return ((85./72.*t-11./3.)*t+3.)*t*t;
      } else if (t < 2) {
        t = 2.-t;
        return (((-23./72*t+2./9.)*t+1./3.)*t+2./9.)*t+1./18.;
      } else if (t < 3) {
        return 1./18.*pow(3.-t,4);
      }
    }
    return 0;
  }
  scalar_type bsp5_03_der(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return ((85./18.*t-11.)*t+6.)*t;
      } else if (t < 2) {
        t = 2.-t;
        return -(((-23./18*t+2./3.)*t+2./3.)*t+2./9.);
      } else if (t < 3) {
        return -2./9.*pow(3.-t,3);
      }
    }
    return 0;
  }
  scalar_type bsp5_03_der2(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (85./6.*t-22.)*t+6.;
      } else if (t < 2) {
        t = 2.-t;
        return (-23./6*t+4./3.)*t+2./3.;
      } else if (t < 3) {
        return 2./3.*pow(3.-t,2);
      }
    }
    return 0;
  }
  scalar_type bsp5_04(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (-25./72.*t+2./3.)*t*t*t;
      } else if (t < 2) {
        t = 2.-t;
        return (((23./72.*t-5./9.)*t-1./3.)*t+4./9.)*t+4./9.;
      } else if (t < 3) {
        t = t-2.;
        return (((-13./72.*t+5./9.)*t-1./3.)*t-4./9.)*t+4./9.;
      } else if (t < 4) {
        return 1./24.*pow(4.-t,4);
      }
    }
    return 0;
  }
  scalar_type bsp5_04_der(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (-25./18.*t+2.)*t*t;
      } else if (t < 2) {
        t = 2.-t;
        return -(((23./18.*t-5./3.)*t-2./3.)*t+4./9.);
      } else if (t < 3) {
        t = t-2.;
        return ((-13./18.*t+5./3.)*t-2./3.)*t-4./9.;
      } else if (t < 4) {
        return -1./6.*pow(4.-t,3);
      }
    }
    return 0;
  }
  scalar_type bsp5_04_der2(scalar_type t) {
    if (t >= -eps) {
      if (t < 1) {
        return (-25./6.*t+4.)*t;
      } else if (t < 2) {
        t = 2.-t;
        return (23./6.*t-10./3.)*t-2./3.;
      } else if (t < 3) {
        t = t-2.;
        return -((-13./6.*t+10./3.)*t-2./3.);
      } else if (t < 4) {
        return 1./2.*pow(4.-t,2);
      }
    }
    return 0;
  }
  scalar_type bsp5_05(scalar_type t) {
    if (t > scalar_type(2.5)) {
      t = 5.-t;
    }
    if (t >= -eps) {
      if (t < 1) {
        return 1./24.*pow(t,4);
      } else if (t < 2) {
        t = t-1.;
        return (((-1./6.*t+1./6.)*t+1./4.)*t+1./6.)*t+1./24.;
      } else if (t < 3) {
        t = pow(t-2.5,2);
        return (1./4.*t-5./8.)*t+115./192.;
      }
    }
    return 0;
  }
  scalar_type bsp5_05_der(scalar_type t) {
    scalar_type sgn(1);
    if (t > scalar_type(2.5)) {
      t = 5.-t;
      sgn = -1.;
    }
    if (t >= -eps) {
      if (t < 1) {
        return 1./6.*pow(t,3)*sgn;
      } else if (t < 2) {
        t = t-1.;
        return (((-2./3.*t+1./2.)*t+1./2)*t+1./6.)*sgn;
      } else if (t < 3) {
        t = t-2.5;
        return (t*t-5./4.)*t*sgn;
      }
    }
    return 0;
  }
  scalar_type bsp5_05_der2(scalar_type t) {
    if (t > scalar_type(2.5)) {
      t = 5.-t;
    }
    if (t >= -eps) {
      if (t < 1) {
        return 1./2.*pow(t,2);
      } else if (t < 2) {
        t = t-1.;
        return (-2.*t+1.)*t+1./2;
      } else if (t < 3) {
        t = t-2.5;
        return 3*t*t-5./4.;
      }
    }
    return 0;
  }



  class global_function_x_bspline_
    : public global_function_simple, public context_dependencies {
    scalar_type xmin, xmax, xscale;
    std::function<scalar_type(scalar_type)> fx, fx_der, fx_der2;
  public:
    void update_from_context() const {}

    virtual scalar_type val(const base_node &pt) const
    {
      return fx(xscale*(pt[0]-xmin));
    }
    virtual void grad(const base_node &pt, base_small_vector &g) const
    {
      scalar_type dx = xscale*(pt[0]-xmin);
      g.resize(dim_);
      g[0] = xscale * fx_der(dx);
    }
    virtual void hess(const base_node &pt, base_matrix &h) const
    {
      scalar_type dx = xscale*(pt[0]-xmin);
      h.resize(dim_, dim_);
      gmm::clear(h);
      h(0,0) = xscale*xscale * fx_der2(dx);
    }

    virtual bool is_in_support(const base_node &pt) const {
      scalar_type dx = pt[0]-(xmin+xmax)/2;
      return (gmm::abs(dx)+1e-9 < gmm::abs(xmax-xmin)/2);
    }

    virtual void bounding_box(base_node &bmin, base_node &bmax) const {
      GMM_ASSERT1(bmin.size() == dim_ && bmax.size() == dim_,
                  "Wrong dimensions");
      bmin[0] = std::min(xmin,xmax);
      bmax[0] = std::max(xmin,xmax);
    }

    global_function_x_bspline_(const scalar_type &xmin_, const scalar_type &xmax_,
                               const size_type &order, const size_type &xtype)
    : global_function_simple(1), xmin(xmin_), xmax(xmax_),
      xscale(scalar_type(xtype)/(xmax-xmin))
    {
      GMM_ASSERT1(order >= 3 && order <= 5, "Only orders 3 to 5 are supported");
      GMM_ASSERT1(xtype >= 1 && xtype <= order, "Wrong input");
      if (order == 3) {
        if (xtype == 1) {
          fx = bsp3_01;   fx_der = bsp3_01_der;   fx_der2 = bsp3_01_der2;
        } else if (xtype == 2) {
          fx = bsp3_02;   fx_der = bsp3_02_der;   fx_der2 = bsp3_02_der2;
        } else if (xtype == 3) {
          fx = bsp3_03;   fx_der = bsp3_03_der;   fx_der2 = bsp3_03_der2;
        }
      } else if (order == 4) {
        if (xtype == 1) {
          fx = bsp4_01;   fx_der = bsp4_01_der;   fx_der2 = bsp4_01_der2;
        } else if (xtype == 2) {
          fx = bsp4_02;   fx_der = bsp4_02_der;   fx_der2 = bsp4_02_der2;
        } else if (xtype == 3) {
          fx = bsp4_03;   fx_der = bsp4_03_der;   fx_der2 = bsp4_03_der2;
        } else if (xtype == 4) {
          fx = bsp4_04;   fx_der = bsp4_04_der;   fx_der2 = bsp4_04_der2;
        }
      } else if (order == 5) {
        if (xtype == 1) {
          fx = bsp5_01;   fx_der = bsp5_01_der;   fx_der2 = bsp5_01_der2;
        } else if (xtype == 2) {
          fx = bsp5_02;   fx_der = bsp5_02_der;   fx_der2 = bsp5_02_der2;
        } else if (xtype == 3) {
          fx = bsp5_03;   fx_der = bsp5_03_der;   fx_der2 = bsp5_03_der2;
        } else if (xtype == 4) {
          fx = bsp5_04;   fx_der = bsp5_04_der;   fx_der2 = bsp5_04_der2;
        } else if (xtype == 5) {
          fx = bsp5_05;   fx_der = bsp5_05_der;   fx_der2 = bsp5_05_der2;
        }
      }
    }

    virtual ~global_function_x_bspline_()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function x bspline"); }
  };



  class global_function_xy_bspline_
    : public global_function_simple, public context_dependencies {
    scalar_type xmin, ymin, xmax, ymax, xscale, yscale;
    std::function<scalar_type(scalar_type)>
      fx, fy, fx_der, fy_der, fx_der2, fy_der2;
  public:
    void update_from_context() const {}

    virtual scalar_type val(const base_node &pt) const
    {
      return fx(xscale*(pt[0]-xmin))
             * fy(yscale*(pt[1]-ymin));
    }
    virtual void grad(const base_node &pt, base_small_vector &g) const
    {
      scalar_type dx = xscale*(pt[0]-xmin),
                  dy = yscale*(pt[1]-ymin);
      g.resize(dim_);
      g[0] = xscale * fx_der(dx) * fy(dy);
      g[1] = fx(dx) * yscale * fy_der(dy);
    }
    virtual void hess(const base_node &pt, base_matrix &h) const
    {
      scalar_type dx = xscale*(pt[0]-xmin),
                  dy = yscale*(pt[1]-ymin);
      h.resize(dim_, dim_);
      gmm::clear(h);
      h(0,0) = xscale*xscale * fx_der2(dx) * fy(dy);
      h(0,1) = h(1,0) = xscale * fx_der(dx) * yscale * fy_der(dy);
      h(1,1) = fx(dx) * yscale*yscale * fy_der2(dy);
    }

    virtual bool is_in_support(const base_node &pt) const {
      scalar_type dx = pt[0]-(xmin+xmax)/2,
                  dy = pt[1]-(ymin+ymax)/2;
      return (gmm::abs(dx)+1e-9 < gmm::abs(xmax-xmin)/2) &&
             (gmm::abs(dy)+1e-9 < gmm::abs(ymax-ymin)/2);
    }

    virtual void bounding_box(base_node &bmin, base_node &bmax) const {
      GMM_ASSERT1(bmin.size() == dim_ && bmax.size() == dim_,
                  "Wrong dimensions");
      bmin[0] = std::min(xmin,xmax);
      bmin[1] = std::min(ymin,ymax);
      bmax[0] = std::max(xmin,xmax);
      bmax[1] = std::max(ymin,ymax);
    }

    global_function_xy_bspline_(const scalar_type &xmin_, const scalar_type &xmax_,
                                const scalar_type &ymin_, const scalar_type &ymax_,
                                const size_type &order,
                                const size_type &xtype, const size_type &ytype)
    : global_function_simple(2), xmin(xmin_), ymin(ymin_), xmax(xmax_), ymax(ymax_),
      xscale(scalar_type(xtype)/(xmax-xmin)), yscale(scalar_type(ytype)/(ymax-ymin))
    {
      GMM_ASSERT1(order >= 3 && order <= 5, "Wrong input");
      GMM_ASSERT1(xtype >= 1 && xtype <= order &&
                  ytype >= 1 && ytype <= order, "Wrong input");
      if (order == 3) {
        if (xtype == 1) {
          fx = bsp3_01;   fx_der = bsp3_01_der;   fx_der2 = bsp3_01_der2;
        } else if (xtype == 2) {
          fx = bsp3_02;   fx_der = bsp3_02_der;   fx_der2 = bsp3_02_der2;
        } else if (xtype == 3) {
          fx = bsp3_03;   fx_der = bsp3_03_der;   fx_der2 = bsp3_03_der2;
        }

        if (ytype == 1) {
          fy = bsp3_01;   fy_der = bsp3_01_der;   fy_der2 = bsp3_01_der2;
        } else if (ytype == 2) {
          fy = bsp3_02;   fy_der = bsp3_02_der;   fy_der2 = bsp3_02_der2;
        } else if (ytype == 3) {
          fy = bsp3_03;   fy_der = bsp3_03_der;   fy_der2 = bsp3_03_der2;
        }
      } else if (order == 4) {
        if (xtype == 1) {
          fx = bsp4_01;   fx_der = bsp4_01_der;   fx_der2 = bsp4_01_der2;
        } else if (xtype == 2) {
          fx = bsp4_02;   fx_der = bsp4_02_der;   fx_der2 = bsp4_02_der2;
        } else if (xtype == 3) {
          fx = bsp4_03;   fx_der = bsp4_03_der;   fx_der2 = bsp4_03_der2;
        } else if (xtype == 4) {
          fx = bsp4_04;   fx_der = bsp4_04_der;   fx_der2 = bsp4_04_der2;
        }

        if (ytype == 1) {
          fy = bsp4_01;   fy_der = bsp4_01_der;   fy_der2 = bsp4_01_der2;
        } else if (ytype == 2) {
          fy = bsp4_02;   fy_der = bsp4_02_der;   fy_der2 = bsp4_02_der2;
        } else if (ytype == 3) {
          fy = bsp4_03;   fy_der = bsp4_03_der;   fy_der2 = bsp4_03_der2;
        } else if (ytype == 4) {
          fy = bsp4_04;   fy_der = bsp4_04_der;   fy_der2 = bsp4_04_der2;
        }

      } else if (order == 5) {
        if (xtype == 1) {
          fx = bsp5_01;   fx_der = bsp5_01_der;   fx_der2 = bsp5_01_der2;
        } else if (xtype == 2) {
          fx = bsp5_02;   fx_der = bsp5_02_der;   fx_der2 = bsp5_02_der2;
        } else if (xtype == 3) {
          fx = bsp5_03;   fx_der = bsp5_03_der;   fx_der2 = bsp5_03_der2;
        } else if (xtype == 4) {
          fx = bsp5_04;   fx_der = bsp5_04_der;   fx_der2 = bsp5_04_der2;
        } else if (xtype == 5) {
          fx = bsp5_05;   fx_der = bsp5_05_der;   fx_der2 = bsp5_05_der2;
        }

        if (ytype == 1) {
          fy = bsp5_01;   fy_der = bsp5_01_der;   fy_der2 = bsp5_01_der2;
        } else if (ytype == 2) {
          fy = bsp5_02;   fy_der = bsp5_02_der;   fy_der2 = bsp5_02_der2;
        } else if (ytype == 3) {
          fy = bsp5_03;   fy_der = bsp5_03_der;   fy_der2 = bsp5_03_der2;
        } else if (ytype == 4) {
          fy = bsp5_04;   fy_der = bsp5_04_der;   fy_der2 = bsp5_04_der2;
        } else if (ytype == 5) {
          fy = bsp5_05;   fy_der = bsp5_05_der;   fy_der2 = bsp5_05_der2;
        }
      }
    }

    virtual ~global_function_xy_bspline_()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function xy bspline"); }
  };


  pglobal_function
  global_function_bspline(const scalar_type xmin, const scalar_type xmax,
                          const size_type order, const size_type xtype) {
    return std::make_shared<global_function_x_bspline_>
                           (xmin, xmax, order, xtype);
  }

  pglobal_function
  global_function_bspline(const scalar_type xmin, const scalar_type xmax,
                          const scalar_type ymin, const scalar_type ymax,
                          const size_type order,
                          const size_type xtype, const size_type ytype) {
    return std::make_shared<global_function_xy_bspline_>
                           (xmin, xmax, ymin, ymax, order, xtype, ytype);
  }




  // interpolator on mesh_fem

  interpolator_on_mesh_fem::interpolator_on_mesh_fem
  (const mesh_fem &mf_, const std::vector<scalar_type> &U_)
    : mf(mf_), U(U_) {

    if (mf.is_reduced()) {
      gmm::resize(U, mf.nb_basic_dof());
      gmm::mult(mf.extension_matrix(), U_, U);
    }
    base_node bmin, bmax;
    scalar_type EPS=1E-13;
    cv_stored = size_type(-1);
    boxtree.clear();
    for (dal::bv_visitor cv(mf.convex_index()); !cv.finished(); ++cv) {
      bounding_box(bmin, bmax, mf.linked_mesh().points_of_convex(cv),
                   mf.linked_mesh().trans_of_convex(cv));
      for (auto&& val : bmin) val -= EPS;
      for (auto&& val : bmax) val += EPS;
      boxtree.add_box(bmin, bmax, cv);
    }
    boxtree.build_tree();
  }

  bool interpolator_on_mesh_fem::find_a_point(const base_node &pt,
                                              base_node &ptr,
                                              size_type &cv) const {
    bool gt_invertible;
    if (cv_stored != size_type(-1) && gic.invert(pt, ptr, gt_invertible)) {
      cv = cv_stored;
      if (gt_invertible)
        return true;
    }

    boxtree.find_boxes_at_point(pt, boxlst);
    for (const auto &box : boxlst) {
      gic = bgeot::geotrans_inv_convex
        (mf.linked_mesh().convex(box->id),
         mf.linked_mesh().trans_of_convex(box->id));
      cv_stored = box->id;
      if (gic.invert(pt, ptr, gt_invertible)) {
        cv = box->id;
        return true;
      }
    }
    return false;
  }

  bool interpolator_on_mesh_fem::eval(const base_node &pt, base_vector &val,
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
    } else
      return false;
  }
}

/* end of namespace getfem  */
