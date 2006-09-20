// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Copyright (C) 2004-2006 Yves Renard
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

/**@file mesh_fem_global_function.h
   @author  Yves Renard <Yves.Renard@insa-toulouse.fr>, J. Pommier
   @date March, 2005.
   @brief Define mesh_fem whose base functions are global function given by the user.
*/
#ifndef GETFEM_GLOBAL_FUNCTION_FEM_H__
#define GETFEM_GLOBAL_FUNCTION_FEM_H__

#include <getfem_fem.h>
#include <getfem_mesh_fem.h>
#include <bgeot_rtree.h>

namespace getfem {
  /// inherit from this class to define new global functions.
  struct global_function : virtual public dal::static_stored_object {
    virtual scalar_type val(const fem_interpolation_context&) const
    { DAL_THROW(dal::failure_error,
		"this global_function has no value"); }
    virtual void grad(const fem_interpolation_context&, base_small_vector &) const
    { DAL_THROW(dal::failure_error,
		"this global_function has no gradient"); }
    virtual void hess(const fem_interpolation_context&, base_matrix &) const
    { DAL_THROW(dal::failure_error,
		"this global_function has no hessian"); }
    virtual ~global_function() {}
  };
  
  typedef boost::intrusive_ptr<const global_function> pglobal_function;

  class global_function_fem : public virtual_fem {
  protected :
    std::vector<pglobal_function> functions;
    mutable bgeot::multi_index mib,mig,mih;
    void init();
  public :

    virtual size_type nb_dof(size_type cv) const;
    virtual size_type index_of_global_dof(size_type cv, size_type i) const;
    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;
    void real_base_value(const fem_interpolation_context& c, 
			 base_tensor &t, bool = true) const;
    void real_grad_base_value(const fem_interpolation_context& c, 
			      base_tensor &t, bool = true) const;
    void real_hess_base_value(const fem_interpolation_context&, 
			      base_tensor &, bool = true) const;

    global_function_fem(bgeot::pconvex_ref cvr_, 
			const std::vector<pglobal_function> &f)
      : functions(f) {
      cvr = cvr_;
      init();
    }
  };
  
  pfem new_global_function_fem(bgeot::pconvex_ref cvr,
			       const std::vector<pglobal_function>& functions);
  
  inline void del_global_function_fem(pfem pf) { dal::del_stored_object(pf); }

  /** mesh_fem whose base functions are global functions (function
      defined on the whole mesh) given by the user. This is much more
      powerful than getfem::external_data_fem.
  */
  class mesh_fem_global_function : public mesh_fem {
  protected :
    mutable std::map<bgeot::pconvex_ref, pfem> build_methods;
    std::vector<pglobal_function> fun;
    void clear_build_methods();
    void init(const std::vector<pglobal_function>& f);
  public :
    void adapt(void);
    void clear(void);
    
    void receipt(const MESH_CLEAR &);
    void receipt(const MESH_DELETE &);
    
    size_type memsize() const { return mesh_fem::memsize(); }
    
    mesh_fem_global_function(mesh &me, dim_type q=1) : mesh_fem(me, q) {}

    void set_functions(pglobal_function f) 
    { fun.resize(1); fun[0]=f; adapt(); }
    void set_functions(pglobal_function f1, pglobal_function f2) 
    { fun.resize(2); fun[0]=f1; fun[1] = f2; adapt(); }
    void set_functions(const std::vector<pglobal_function>& f) 
    { fun = f; adapt(); }
    ~mesh_fem_global_function() { clear_build_methods(); }
  };


  /** a general structure for interpolation of a function defined 
      by a mesh_fem and a vector U at any point 
      (interpolation of value and radient).
  */
  struct interpolator_on_mesh_fem {
    const mesh_fem &mf;
    const std::vector<scalar_type> U;

    mutable bgeot::rtree boxtree;
    mutable size_type cv_stored;
    mutable bgeot::rtree::pbox_set boxlst;
    mutable bgeot::geotrans_inv_convex gic;
    
    
    interpolator_on_mesh_fem(const mesh_fem &mf_, 
			     const std::vector<scalar_type> &U_) :
      mf(mf_), U(U_) { init(); }

    void init();
    bool find_a_point(base_node pt, base_node &ptr,
		      size_type &cv) const;
    bool eval(const base_node pt, base_vector &val, base_matrix &grad) const;
  };


  /* below a list of simple function of (x,y)
     used for building the crack singular functions 
  */
  struct abstract_xy_function {
    virtual scalar_type val(scalar_type x, scalar_type y) const = 0;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const = 0;
    virtual ~abstract_xy_function() {}
  };

  
  struct crack_singular_xy_function : public abstract_xy_function {
    unsigned l; /* 0 <= l <= 3 */
    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    crack_singular_xy_function(unsigned l_) : l(l_) {}
  };

  struct cutoff_xy_function : public abstract_xy_function {
    enum { NOCUTOFF = -1, 
	   EXPONENTIAL_CUTOFF = 0, POLYNOMIAL_CUTOFF = 1, 
	   POLYNOMIAL2_CUTOFF=2 };
    int fun;
    scalar_type a4, r1, r0;
    virtual scalar_type val(scalar_type x, scalar_type y) const;
    virtual base_small_vector grad(scalar_type x, scalar_type y) const;
    cutoff_xy_function(int fun_num, scalar_type r, 
		       scalar_type r1, scalar_type r0);
  };

  struct interpolated_xy_function : public abstract_xy_function {
    interpolator_on_mesh_fem &itp;
    size_type component;
    virtual scalar_type val(scalar_type x, scalar_type y) const {
      base_vector v; base_matrix g;
      itp.eval(base_node(x,y), v, g);
      return v[component];
    }
    virtual base_small_vector grad(scalar_type x, scalar_type y) const {
      base_vector v; base_matrix g;
      itp.eval(base_node(x,y), v, g);
      return base_small_vector(g(0,component), g(1,component));      
    }
    interpolated_xy_function(interpolator_on_mesh_fem &itp_, size_type c) :
      itp(itp_), component(c) {}
  };

  struct product_of_xy_functions : 
    public abstract_xy_function {
    abstract_xy_function &fn1, &fn2;
    scalar_type val(scalar_type x, scalar_type y) const {
      return fn1.val(x,y) * fn2.val(x,y);
    }
    base_small_vector grad(scalar_type x, scalar_type y) const {
      return fn1.grad(x,y)*fn2.val(x,y) + fn1.val(x,y)*fn2.grad(x,y);
    }
    product_of_xy_functions(abstract_xy_function &fn1_,
			    abstract_xy_function &fn2_) 
      : fn1(fn1_), fn2(fn2_) {}
  };



  /*
   * some usefull global functions
   */
  class level_set;

  pglobal_function 
  global_function_on_level_set(const level_set &ls,
			       const abstract_xy_function &fn);

  
}  /* end of namespace getfem.                                            */

#endif
  
