// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_interpolated_fem.h : definition of a finite element
//           method which interpolates a fem on a different mesh.
// Date    : March, 2005.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//           J. Pommier
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

/**@file mesh_fem_global_function.h
   @brief Define mesh_fem whose base functions are global function given by the user.
*/
#ifndef GETFEM_GLOBAL_FUNCTION_FEM_H__
#define GETFEM_GLOBAL_FUNCTION_FEM_H__

#include <getfem_fem.h>
#include <getfem_mesh_fem.h>

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


  /*
   * some usefull global functions
   */
  class level_set;
  pglobal_function isotropic_crack_singular_2D(size_type i,
					       const level_set &ls,
					       scalar_type cutoff_radius = 0,
					       scalar_type cutoff_radius1 = 0,
					       scalar_type cutoff_radius0 = 0,
					       size_type func = 0);



}  /* end of namespace getfem.                                            */

#endif
  
