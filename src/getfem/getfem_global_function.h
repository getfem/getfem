/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2004-2020 Yves Renard
 Copyright (C) 2016-2020 Konstantinos Poulios

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

 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.

===========================================================================*/

/**@file getfem_global_function.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>, J. Pommier
   @date March, 2005.
   @brief Definition of global functions to be used as base or enrichment
          functions in fem.
*/
#ifndef GETFEM_GLOBAL_FUNCTION_H__
#define GETFEM_GLOBAL_FUNCTION_H__

#include "bgeot_rtree.h"
#include "getfem_mesh_fem.h"
#include "getfem_generic_assembly.h"

namespace getfem {

  /// inherit from this class to define new global functions.
  class global_function : virtual public dal::static_stored_object {
  protected:
    const dim_type dim_;
  public:
    dim_type dim() const { return dim_; }

    virtual scalar_type val(const fem_interpolation_context&) const
    { GMM_ASSERT1(false, "this global_function has no value"); }
    virtual void grad(const fem_interpolation_context&, base_small_vector&) const
    { GMM_ASSERT1(false, "this global_function has no gradient"); }
    virtual void hess(const fem_interpolation_context&, base_matrix&) const
    { GMM_ASSERT1(false, "this global_function has no hessian"); }

    virtual bool is_in_support(const base_node & /* pt */ ) const { return true; }
    virtual void bounding_box(base_node &bmin, base_node &bmax) const {
      GMM_ASSERT1(bmin.size() == dim_ && bmax.size() == dim_,
                  "Wrong dimensions");
      for (auto&& xx : bmin) xx = -1e+25;
      for (auto&& xx : bmax) xx = 1e+25;
    }

    global_function(dim_type dim__) : dim_(dim__)
    { DAL_STORED_OBJECT_DEBUG_CREATED(this, "Global function");}
    virtual ~global_function()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function"); }
  };

  typedef std::shared_ptr<const global_function> pglobal_function;


  class global_function_simple : public global_function {
  public:
    // These virtual functions can not be further overriden in derived classes
    virtual scalar_type val(const fem_interpolation_context&) const final;
    virtual void grad(const fem_interpolation_context&, base_small_vector&) const final;
    virtual void hess(const fem_interpolation_context&, base_matrix&) const final;
    // You have to override these instead
    virtual scalar_type val(const base_node &pt) const = 0;
    virtual void grad(const base_node &pt, base_small_vector&) const = 0;
    virtual void hess(const base_node &pt, base_matrix&) const = 0;

    global_function_simple(dim_type dim__) : global_function(dim__)
    { DAL_STORED_OBJECT_DEBUG_CREATED(this, "Global function simple"); }
    virtual ~global_function_simple()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function simple"); }
  };

  class global_function_parser : public global_function_simple {
    ga_workspace gw;
    ga_function f_val, f_grad, f_hess;
    mutable model_real_plain_vector pt_;
  public:
    virtual scalar_type val(const base_node &pt) const;
    virtual const base_tensor &tensor_val(const base_node &pt) const;
    virtual void grad(const base_node &pt, base_small_vector &g) const;
    virtual void hess(const base_node &pt, base_matrix &h) const;

    global_function_parser(dim_type dim_,
                           const std::string &sval,
                           const std::string &sgrad="",
                           const std::string &shess="");
    virtual ~global_function_parser()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function parser"); }
  };


  class global_function_sum : public global_function {
    std::vector<pglobal_function> functions;
  public:
    virtual scalar_type val(const fem_interpolation_context&) const;
    virtual void grad(const fem_interpolation_context&, base_small_vector&) const;
    virtual void hess(const fem_interpolation_context&, base_matrix&) const;
    virtual bool is_in_support(const base_node &p) const;
    virtual void bounding_box(base_node &bmin_, base_node &bmax_) const;
    global_function_sum(const std::vector<pglobal_function> &funcs);
    global_function_sum(pglobal_function f1, pglobal_function f2);
    global_function_sum(pglobal_function f1, pglobal_function f2,
                        pglobal_function f3);
    global_function_sum(pglobal_function f1, pglobal_function f2,
                        pglobal_function f3, pglobal_function f4);
    virtual ~global_function_sum()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function sum"); }
  };

  class global_function_product : public global_function {
    pglobal_function f1, f2;
  public:
    virtual scalar_type val(const fem_interpolation_context&) const;
    virtual void grad(const fem_interpolation_context&, base_small_vector&) const;
    virtual void hess(const fem_interpolation_context&, base_matrix&) const;
    virtual bool is_in_support(const base_node &p) const;
    virtual void bounding_box(base_node &bmin_, base_node &bmax_) const;
    global_function_product(pglobal_function f1_, pglobal_function f2_);
    virtual ~global_function_product()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function product"); }
  };

  class global_function_bounded : public global_function {
    pglobal_function f;
    const base_node bmin, bmax;
    bool has_expr;
    ga_workspace gw;
    ga_function f_val;
    mutable model_real_plain_vector pt_;
  public:
    virtual scalar_type val(const fem_interpolation_context &c) const
    { return f->val(c); }
    virtual void grad(const fem_interpolation_context &c, base_small_vector &g) const
    { f->grad(c, g); }
    virtual void hess(const fem_interpolation_context &c, base_matrix &h) const
    { f->hess(c, h); }

    virtual bool is_in_support(const base_node &) const;
    virtual void bounding_box(base_node &bmin_, base_node &bmax_) const {
      bmin_ = bmin;
      bmax_ = bmax;
    }
    // is_in_expr should evaluate to a positive result inside the support of the
    // function and negative elsewhere
    global_function_bounded(pglobal_function f_,
                            const base_node &bmin_, const base_node &bmax_,
                            const std::string &is_in_expr="");
    virtual ~global_function_bounded()
    { DAL_STORED_OBJECT_DEBUG_DESTROYED(this, "Global function bounded"); }
  };


  /** a general structure for interpolation of a function defined
      by a mesh_fem and a vector U at any point
      (interpolation of value and gradient).
  */

  struct interpolator_on_mesh_fem {
    const mesh_fem &mf;
    std::vector<scalar_type> U;

    mutable bgeot::rtree boxtree;
    mutable size_type cv_stored;
    mutable bgeot::rtree::pbox_set boxlst;
    mutable bgeot::geotrans_inv_convex gic;


    interpolator_on_mesh_fem(const mesh_fem &mf_,
                             const std::vector<scalar_type> &U_);
    bool find_a_point(const base_node &pt, base_node &ptr,
                      size_type &cv) const;
    bool eval(const base_node &pt, base_vector &val, base_matrix &grad) const;
  };

  typedef std::shared_ptr<const interpolator_on_mesh_fem>
    pinterpolator_on_mesh_fem;


  /** below a list of simple functions of (x,y)
      used for building the crack singular functions
  */
  struct abstract_xy_function : virtual public dal::static_stored_object {
    virtual scalar_type val(scalar_type x, scalar_type y) const = 0;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const = 0;
    virtual base_matrix hess(scalar_type x, scalar_type y) const = 0;
    virtual ~abstract_xy_function() {}
  };

  typedef std::shared_ptr<const abstract_xy_function> pxy_function;

  struct parser_xy_function : public abstract_xy_function {
    ga_workspace gw;
    ga_function f_val, f_grad, f_hess;
    mutable model_real_plain_vector ptx, pty, ptr, ptt;

    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    virtual base_matrix hess(scalar_type x, scalar_type y) const;

    parser_xy_function(const std::string &sval,
                       const std::string &sgrad="0;0;",
                       const std::string &shess="0;0;0;0;");
    virtual ~parser_xy_function() {}
  };

  struct crack_singular_xy_function : public abstract_xy_function {
    unsigned l; /* 0 <= l <= 6 */
    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    virtual base_matrix hess(scalar_type x, scalar_type y) const;
    crack_singular_xy_function(unsigned l_) : l(l_) {}
    virtual ~crack_singular_xy_function() {}
  };

  struct cutoff_xy_function : public abstract_xy_function {
    enum { NOCUTOFF = -1,
           EXPONENTIAL_CUTOFF = 0,
           POLYNOMIAL_CUTOFF = 1,
           POLYNOMIAL2_CUTOFF= 2 };
    int fun;
    scalar_type a4, r1, r0;
    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    virtual base_matrix hess(scalar_type x, scalar_type y) const;
    cutoff_xy_function(int fun_num, scalar_type r,
                       scalar_type r1, scalar_type r0);
    virtual ~cutoff_xy_function() {}
  };

  struct interpolated_xy_function : public abstract_xy_function {
    pinterpolator_on_mesh_fem itp;
    size_type component;
    virtual scalar_type val(scalar_type x, scalar_type y) const {
      base_vector v; base_matrix g;
      itp->eval(base_node(x,y), v, g);
      return v[component];
    }
    virtual base_small_vector grad(scalar_type x, scalar_type y) const {
      base_vector v; base_matrix g;
      itp->eval(base_node(x,y), v, g);
      return base_small_vector(g(component,0), g(component,1));
    }
    virtual base_matrix hess(scalar_type, scalar_type) const
    { GMM_ASSERT1(false, "Sorry, to be done ..."); }
    interpolated_xy_function(const pinterpolator_on_mesh_fem &itp_,
                             size_type c) :
      itp(itp_), component(c) {}
    virtual ~interpolated_xy_function() {}
  };

  struct product_of_xy_functions : public abstract_xy_function {
    pxy_function fn1, fn2;
    scalar_type val(scalar_type x, scalar_type y) const {
      return fn1->val(x,y) * fn2->val(x,y);
    }
    base_small_vector grad(scalar_type x, scalar_type y) const {
      return fn1->grad(x,y)*fn2->val(x,y) + fn1->val(x,y)*fn2->grad(x,y);
    }
    virtual base_matrix hess(scalar_type x, scalar_type y) const {
      base_matrix h = fn1->hess(x,y);
      gmm::scale(h, fn2->val(x,y));
      gmm::add(gmm::scaled(fn2->hess(x,y), fn1->val(x,y)), h);
      gmm::rank_two_update(h, fn1->grad(x,y), fn2->grad(x,y));
      return h;
    }
    product_of_xy_functions(pxy_function &fn1_, pxy_function &fn2_)
      : fn1(fn1_), fn2(fn2_) {}
    virtual ~product_of_xy_functions() {}
  };

  struct add_of_xy_functions : public abstract_xy_function {
    pxy_function fn1, fn2;
    scalar_type val(scalar_type x, scalar_type y) const {
      return fn1->val(x,y) + fn2->val(x,y);
    }
    base_small_vector grad(scalar_type x, scalar_type y) const {
      return fn1->grad(x,y) + fn2->grad(x,y);
    }
    virtual base_matrix hess(scalar_type x, scalar_type y) const {
      base_matrix h = fn1->hess(x,y);
      gmm::add(fn2->hess(x,y), h);
      return h;
    }
    add_of_xy_functions(const pxy_function &fn1_, const pxy_function &fn2_)
      : fn1(fn1_), fn2(fn2_) {}
    virtual ~add_of_xy_functions() {}
  };


  /** some useful global functions */
  class level_set;

  pglobal_function
  global_function_on_level_set(const level_set &ls, const pxy_function &fn);

  pglobal_function
  global_function_on_level_sets(const std::vector<level_set> &lsets,
                                const pxy_function &fn);

  pglobal_function
  global_function_bspline(const scalar_type xmin, const scalar_type xmax,
                          const size_type order, const size_type xtype);

  pglobal_function
  global_function_bspline(const scalar_type xmin, const scalar_type xmax,
                          const scalar_type ymin, const scalar_type ymax,
                          const size_type order,
                          const size_type xtype, const size_type ytype);

}  /* end of namespace getfem.                                            */

#endif
