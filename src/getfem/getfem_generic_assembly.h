/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2013-2014 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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

/** @file getfem_generic_assemblt.h
    @author  Yves Renard <Yves.Renard@insa-lyon.fr>
    @date November 18, 2013.
    @brief A langage for generic assembly of pde boundary value problems.
 */


#ifndef GETFEM_GENERIC_ASSEMBLY_H__ 
#define GETFEM_GENERIC_ASSEMBLY_H__

#include "getfem/getfem_models.h"

namespace getfem {

  #ifdef _WIN32
  #define INFINITY std::numeric_limits<scalar_type>::infinity()
  typedef double (*BoostMathFunction)(double);
  BoostMathFunction const acosh = boost::math::acosh<double>;
  BoostMathFunction const asinh = boost::math::asinh<double>;
  BoostMathFunction const atanh = boost::math::atanh<double>;
  BoostMathFunction const erf = boost::math::erf<double>;
  BoostMathFunction const erfc = boost::math::erfc<double>;
  #endif

  class ga_tree;

  //=========================================================================
  // Structure dealing with predefined operators.
  //=========================================================================

  struct ga_nonlinear_operator {

    typedef std::vector<const base_tensor *> arg_list;

    virtual bool result_size(const arg_list &args,
                             bgeot::multi_index &sizes) const = 0;

    virtual void value(const arg_list &args, base_tensor &result) const = 0;
    
    virtual void derivative(const arg_list &args, size_type i,
                            base_tensor &result) const = 0;

    virtual void second_derivative(const arg_list &args, size_type i,
                                   size_type j, base_tensor &result) const = 0;
  };

  typedef std::map<std::string, ga_nonlinear_operator*> ga_predef_operator_tab;

  //=========================================================================
  // Structure dealing with predefined scalar functions.
  //=========================================================================
  
  typedef scalar_type (*pscalar_func_onearg)(scalar_type);
  typedef scalar_type (*pscalar_func_twoargs)(scalar_type, scalar_type);
  struct ga_interval {
    scalar_type min, max;
    ga_interval(void) { min = -INFINITY; max = +INFINITY; }
    ga_interval(scalar_type a, scalar_type b) { min = a; max = b; }
  };
  
  
  struct ga_predef_function {
    size_type ftype; // 0 : C++ function with a derivative
                     // 1 : C++ function with an expression to be derived
                     // 2 : function defined by an expression
    size_type nbargs;         // One or two arguments
    pscalar_func_onearg f1;   // Function pointer for a one argument function
    pscalar_func_twoargs f2;  // Function pointer for a two arguments function
    std::string expr;
    ga_interval support1, support2;
    ga_interval domain1, domain2;  // Domain of definition of the function
    std::string derivative1, derivative2;
    
    ga_predef_function(void) {}
    ga_predef_function(pscalar_func_onearg f, const ga_interval &s,
                       const ga_interval &dom, const std::string &der)
      : ftype(0), nbargs(1), f1(f), support1(s), domain1(dom),
        derivative1(der) {}
    ga_predef_function(pscalar_func_onearg f, const std::string &e,
                       const ga_interval &s, const ga_interval &dom)
      : ftype(1), nbargs(1), f1(f), expr(e), support1(s), domain1(dom) {}
    ga_predef_function(pscalar_func_twoargs f, const ga_interval &s1,
                       const ga_interval &s2, const ga_interval &dom1,
                       const ga_interval &dom2, const std::string &der1,
                       const std::string &der2)
      : ftype(0), nbargs(2), f2(f), support1(s1), support2(s2), domain1(dom1),
        domain2(dom2), derivative1(der1), derivative2(der2) {}
    ga_predef_function(pscalar_func_twoargs f, const std::string &e,
                       const ga_interval &s1, const ga_interval &s2,
                       const ga_interval &dom1, const ga_interval &dom2)
      : ftype(0), nbargs(2), f2(f), expr(e), support1(s1), support2(s2),
        domain1(dom1), domain2(dom2) {}
  };
  
  typedef std::map<std::string, ga_predef_function> ga_predef_function_tab;


  //=========================================================================
  // Structure dealing with user defined environment : constant, variables,
  // functions, operators.
  //=========================================================================



  class ga_workspace {
    
    const getfem::model *model;

    struct var_description {

      bool is_variable;
      bool is_fem_dofs;
      const mesh_fem *mf;
      gmm::sub_interval I;
      const model_real_plain_vector *V;

      var_description(bool is_var,
                      bool is_fem, 
                      const mesh_fem *mmf,
                      gmm::sub_interval I_,
                      const model_real_plain_vector *v)
        : is_variable(is_var), is_fem_dofs(is_fem), mf(mmf), I(I_), V(v) {}
      var_description() : is_variable(false), is_fem_dofs(false),
                          mf(0), V(0) {}
    };

  public:

    struct tree_description {
      size_type order; // 0: potential, 1: weak form, 2: tangent operator
      std::string name_test1, name_test2;
      const mesh_im *mim;
      size_type region;
      ga_tree *ptree;
      base_vector elem;
      tree_description(void);
    };

  private:

    typedef std::map<std::string, var_description> VAR_SET;

    VAR_SET variables;
    ga_predef_operator_tab user_operators;
    ga_predef_function_tab user_functions;
    std::vector<tree_description> trees;
    std::list<ga_tree *> aux_trees;

    void add_tree(ga_tree &tree, const mesh_im &mim, size_type region,
                  const std::string expr);
    void clear_aux_trees(void);

    model_real_sparse_matrix unreduced_K, *K;
    base_vector unreduced_V, *V;
    scalar_type E;
    bool K_to_delete, V_to_delete;

  public:

    const model_real_sparse_matrix &assembled_matrix(void) const { return *K; }
    model_real_sparse_matrix &assembled_matrix(void) { return *K; }
    scalar_type &assembled_potential(void) { return E; }
    const scalar_type &assembled_potential(void) const { return E; }
    const base_vector &assembled_vector(void) const { return *V; }
    base_vector &assembled_vector(void) { return *V; }
    void set_assembled_matrix(model_real_sparse_matrix &K_)
    { if (K_to_delete) delete K; K = &K_; K_to_delete = false; }
    void set_assembled_vector(base_vector &V_)
    { if (V_to_delete) delete V; V = &V_; V_to_delete = false; }

    model_real_sparse_matrix &unreduced_matrix(void)
    { return unreduced_K; }
    base_vector &unreduced_vector(void) { return unreduced_V; }
    

    void add_expression(const std::string expr, const mesh_im &mim,
                        size_type region = size_type(-1));
    void clear_expressions(void);
    

    void add_aux_tree(ga_tree &tree);
    size_type nb_trees(void);
    tree_description &tree_info(size_type i);

    bool user_operator_exists(const std::string &name) const
    { return user_operators.find(name) != user_operators.end(); }

    bool user_function_exists(const std::string &name) const
    { return user_functions.find(name) != user_functions.end(); }

    const ga_nonlinear_operator &user_operator(const std::string &name) const
    { return *(user_operators.find(name)->second); }

    const ga_predef_function &user_function(const std::string &name) const
    { return user_functions.find(name)->second; }

    // TODO: methods to add a function or an operator
    
    void add_fem_variable(const std::string &name, const mesh_fem &mf,
                          const gmm::sub_interval &I,
                          const model_real_plain_vector &VV) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(true, true, &mf, I, &VV);
    }
    
    void add_fixed_size_variable(const std::string &name,
                                 const gmm::sub_interval &I,
                                 const model_real_plain_vector &VV) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(true, false, 0, I, &VV);
    }

    void add_fem_constant(const std::string &name, const mesh_fem &mf,
                          const model_real_plain_vector &VV) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(false, true, &mf,
                                        gmm::sub_interval(), &VV);
    }
    
    void add_fixed_size_constant(const std::string &name,
                                 const model_real_plain_vector &VV) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(false, false, 0,
                                        gmm::sub_interval(), &VV);
    }

    bool variable_exists(const std::string &name) const {
      if (model)
        return model->variable_exists(name);
      else
        return (variables.find(name) != variables.end());
    }

    bool is_constant(const std::string &name) const {
      if (model)
        return model->is_data(name);
      else {
        VAR_SET::const_iterator it = variables.find(name);
        GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
        return !(it->second.is_variable);
      }
    }

    const gmm::sub_interval &
    interval_of_variable(const std::string &name) const {
      if (model)
        return model->interval_of_variable(name);
      else {
        VAR_SET::const_iterator it = variables.find(name);
        GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
        return it->second.I;
      }
    }

    const mesh_fem *associated_mf(const std::string &name) const {
      if (model)
        return model->pmesh_fem_of_variable(name);
      else {
        VAR_SET::const_iterator it = variables.find(name);
        GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
        return it->second.is_fem_dofs ? it->second.mf : 0;
      }
    }

    size_type qdim(const std::string &name) const {
      const mesh_fem *mf = associated_mf(name);
      size_type n = gmm::vect_size(value(name));
      size_type ndof = mf ? mf->nb_dof() : 0;
      return mf ? associated_mf(name)->get_qdim() * (ndof / n) : n;
    }

    const model_real_plain_vector &value(const std::string &name) const {
      if (model)
        return model->real_variable(name);
      else {
        VAR_SET::const_iterator it = variables.find(name);
        GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
        return *(it->second.V);
      }
    }

    void assembly(size_type order);

    ga_workspace(const getfem::model &md);
    ga_workspace(void);
    ~ga_workspace();

  };



}  /* end of namespace getfem.                                             */


#endif /* GETFEM_GENERIC_ASSEMBLY_H__  */
