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

/** @file getfem_generic_assembly.h
    @author  Yves Renard <Yves.Renard@insa-lyon.fr>
    @date November 18, 2013.
    @brief A langage for generic assembly of pde boundary value problems.
 */


#ifndef GETFEM_GENERIC_ASSEMBLY_H__ 
#define GETFEM_GENERIC_ASSEMBLY_H__

#include <map>
#include "getfem/getfem_models.h"
#include "getfem/getfem_interpolation.h"


#ifdef _WIN32
#include <limits>
#define INFINITY std::numeric_limits<scalar_type>::infinity()
#endif

namespace getfem {

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

    virtual ~ga_nonlinear_operator() {}
  };

  struct ga_predef_operator_tab {
    typedef std::map<std::string, ga_nonlinear_operator*> T;
    std::map<std::string, ga_nonlinear_operator*> tab;
    
   void add_method(const std::string &name, ga_nonlinear_operator *pt)
    { tab[name] = pt; }
    ~ga_predef_operator_tab() {
      for (T::iterator it = tab.begin(); it != tab.end(); ++it)
        delete it->second;
    }
  };

  //=========================================================================
  // For user predefined scalar functions.
  //=========================================================================
  
  typedef scalar_type (*pscalar_func_onearg)(scalar_type);
  typedef scalar_type (*pscalar_func_twoargs)(scalar_type, scalar_type);
  struct ga_interval {
    scalar_type min, max;
    ga_interval(void) { min = -INFINITY; max = +INFINITY; }
    ga_interval(scalar_type a, scalar_type b) { min = a; max = b; }
  };

  void ga_define_function(const std::string name, size_type nb_args,
                          const std::string expr, const std::string der1="",
                          const std::string der2="");
  void ga_define_function(const std::string name, pscalar_func_onearg f,
                          const std::string der1="");
  void ga_define_function(const std::string name, pscalar_func_twoargs f2,
                          const std::string &der1="",
                          const std::string &der2="");
 
  void ga_undefine_function(const std::string name);
  bool ga_function_exists(const std::string name);

  //=========================================================================
  // Structure dealing with user defined environment : constant, variables,
  // functions, operators.
  //=========================================================================

  class ga_workspace {
    
    const model *md;
    const ga_workspace *parent_workspace;

    struct var_description {

      bool is_variable;
      bool is_fem_dofs;
      const mesh_fem *mf;
      gmm::sub_interval I;
      const model_real_plain_vector *V;
      const im_data *imd;

      var_description(bool is_var,
                      bool is_fem, 
                      const mesh_fem *mmf,
                      gmm::sub_interval I_,
                      const model_real_plain_vector *v, const im_data *imd_)
        : is_variable(is_var), is_fem_dofs(is_fem), mf(mmf), I(I_), V(v),
          imd(imd_) {}
      var_description() : is_variable(false), is_fem_dofs(false),
                          mf(0), V(0), imd(0) {}
    };

  public:

    struct tree_description { // CAUTION: Specific copy constructor
      size_type order; // 0: potential, 1: weak form, 2: tangent operator
      std::string name_test1, name_test2;
      std::string interpolate_name_test1, interpolate_name_test2;
      const mesh_im *mim;
      const mesh *m;
      const mesh_region *rg;
      ga_tree *ptree;
      base_vector elem;
      tree_description(void) : ptree(0) {}
      void copy(const tree_description& td);
      tree_description(const tree_description& td) { copy(td); }
      tree_description &operator =(const tree_description& td);
      ~tree_description();
    };

  private:

    std::map<const mesh *, std::list<mesh_region> > registred_mims;
    mesh dummy_mesh;
    mesh dummy_mim;
    mesh_region dummy_region;
    
    mesh_region &register_region(const mesh &m, const mesh_region &region);

    typedef std::map<std::string, var_description> VAR_SET;

    VAR_SET variables;
    std::map<std::string, pinterpolate_transformation> transformations;
    std::vector<tree_description> trees;
    std::list<ga_tree *> aux_trees;

    std::map<std::string, std::vector<std::string> > variable_groups;
    std::map<std::string, std::string> macros;

    struct m_tree {
      ga_tree *ptree;
      size_type meshdim;
      bool ignore_X;
      m_tree(void) : ptree(0) {}
      ~m_tree(void);
    };
    
    mutable std::map<std::string, m_tree> macro_trees;


    void add_tree(ga_tree &tree, const mesh &m, const mesh_im &mim,
                  const mesh_region &rg,
                  const std::string &expr, bool add_derivative = true,
                  bool scalar_expr = true);
    void clear_aux_trees(void);

    struct sparse_matrix_ptr {
      bool todelete;
      model_real_sparse_matrix *ptr;
      model_real_sparse_matrix &operator()(void) { return *ptr; }
      const model_real_sparse_matrix &operator()(void) const { return *ptr; }
      void resize(size_type nb)
      { if (todelete) { gmm::clear(*ptr); gmm::resize(*ptr, nb, nb); } }
      void set_matrix(model_real_sparse_matrix &M)
      { if (todelete) delete ptr; todelete = false; ptr = &M; }
      sparse_matrix_ptr(void):
        todelete(true), ptr(new model_real_sparse_matrix(2,2)) {}
      sparse_matrix_ptr(const sparse_matrix_ptr &smp):
        todelete(smp.todelete), ptr(smp.ptr)
      { if (todelete) ptr = new model_real_sparse_matrix(smp()); }
      sparse_matrix_ptr &operator =(const sparse_matrix_ptr &smp) {
        if (todelete) delete ptr;
        todelete = smp.todelete; ptr = smp.ptr;
        if (todelete) ptr = new model_real_sparse_matrix(smp());
        return *this;
      }
      ~sparse_matrix_ptr() { if (todelete) delete ptr; }
    };

    struct base_vector_ptr {
      bool todelete;
      base_vector *ptr;
      base_vector &operator()(void) { return *ptr; }
      const base_vector &operator()(void) const { return *ptr; }
      void resize(size_type nb)
      { if (todelete) { gmm::clear(*ptr); gmm::resize(*ptr, nb);} }
      void set_vector(base_vector &vector)
      { if (todelete) delete ptr; todelete = false; ptr = &vector; }
      base_vector_ptr(void):
        todelete(true), ptr(new base_vector(2)) {}
      base_vector_ptr(const base_vector_ptr &smp):
        todelete(smp.todelete), ptr(smp.ptr)
      { if (todelete) ptr = new base_vector(smp()); }
      base_vector_ptr &operator =(const base_vector_ptr &smp) {
        if (todelete) delete ptr;
        todelete = smp.todelete; ptr = smp.ptr;
        if (todelete) ptr = new base_vector(smp());
        return *this;
      }
      ~base_vector_ptr() { if (todelete) delete ptr; }
    };

    sparse_matrix_ptr K;
    model_real_sparse_matrix unreduced_K;
    base_vector_ptr V;
    base_vector unreduced_V;
    scalar_type E;
    base_tensor assemb_t;

  public:

    const model_real_sparse_matrix &assembled_matrix(void) const { return K();}
    model_real_sparse_matrix &assembled_matrix(void) { return K(); }
    scalar_type &assembled_potential(void) { return E; }
    const scalar_type &assembled_potential(void) const { return E; }
    const base_vector &assembled_vector(void) const { return V(); }
    base_vector &assembled_vector(void) { return V(); }
    void set_assembled_matrix(model_real_sparse_matrix &K_)
    { K.set_matrix(K_); }
    void set_assembled_vector(base_vector &V_)
    { V.set_vector(V_); }
    base_tensor &assembled_tensor(void) { return assemb_t; }
    const base_tensor &assembled_tensor(void) const { return assemb_t; }

    model_real_sparse_matrix &unreduced_matrix(void)
    { return unreduced_K; }
    base_vector &unreduced_vector(void) { return unreduced_V; }
    
    /** Add an expression, perform the semantic analysis, split into
     *  terms in separated test functions, derive if necessary to obtain
     *  the tangent terms. Return the maximal order found in the expression.
     */
    size_type add_expression(const std::string expr, const mesh_im &mim,
                            const mesh_region &rg=mesh_region::all_convexes());
    /* Internal use */
    void add_function_expression(const std::string expr);
    /* Internal use */
    void add_interpolation_expression(const std::string expr, const mesh &m,
                                      mesh_region rg);

    /** Delete all previously added expressions. */
    void clear_expressions(void);
    

    void add_aux_tree(ga_tree &tree);
    size_type nb_trees(void) const;
    tree_description &tree_info(size_type i);
        
    void add_fem_variable(const std::string &name, const mesh_fem &mf,
                          const gmm::sub_interval &I,
                          const model_real_plain_vector &VV)
    { variables[name] = var_description(true, true, &mf, I, &VV, 0); }
    
    void add_fixed_size_variable(const std::string &name,
                                 const gmm::sub_interval &I,
                                 const model_real_plain_vector &VV)
    { variables[name] = var_description(true, false, 0, I, &VV, 0); }

    void add_fem_constant(const std::string &name, const mesh_fem &mf,
                          const model_real_plain_vector &VV) { 
      variables[name] = var_description(false, true, &mf,
                                        gmm::sub_interval(), &VV, 0);
    }
    
    void add_fixed_size_constant(const std::string &name,
                                 const model_real_plain_vector &VV) {
      variables[name] = var_description(false, false, 0,
                                        gmm::sub_interval(), &VV, 0);
    }

    void add_im_data(const std::string &name, const im_data &imd,
                     const model_real_plain_vector &VV) {
      variables[name] = var_description(false, false, 0,
                                        gmm::sub_interval(), &VV, &imd);
    }

    bool used_variables(model::varnamelist &vl, model::varnamelist &dl,
                        size_type order);

    bool variable_exists(const std::string &name) const {
      return (md && md->variable_exists(name)) ||
        (parent_workspace && parent_workspace->variable_exists(name)) ||
        (variables.find(name) != variables.end());
    }

    const std::string &variable_in_group(const std::string &group_name,
                                         const mesh &m) const;

    void define_variable_group(const std::string &group_name,
                               const std::vector<std::string> &nl);

    bool variable_group_exists(std::string name) const {
      return (variable_groups.find(name) != variable_groups.end()) ||
        (md && md->variable_group_exists(name)) ||
        (parent_workspace && parent_workspace->variable_group_exists(name));
    }

    bool variable_or_group_exists(const std::string &name) const
    { return variable_exists(name) || variable_group_exists(name); }

    const std::vector<std::string>
    &variable_group(const std::string &group_name) const {
      std::map<std::string, std::vector<std::string> >::const_iterator
        it = variable_groups.find(group_name);
      if (it != variable_groups.end())
        return (variable_groups.find(group_name))->second;
      if (md && md->variable_group_exists(group_name))
        return md->variable_group(group_name);
      if (parent_workspace &&
          parent_workspace->variable_group_exists(group_name))
        return parent_workspace->variable_group(group_name);
      GMM_ASSERT1(false, "Undefined variable group " << group_name);
    }

    const std::string &first_variable_of_group(const std::string &name) const {
      const std::vector<std::string> &t = variable_group(name);
      GMM_ASSERT1(t.size(), "Variable group " << name << " has no variable");
      return t[0];
    }
    
    bool macro_exists(const std::string &name) const {
      if (macros.find(name) != macros.end()) return true;
      if (md && md->macro_exists(name)) return true;
      if (parent_workspace &&
          parent_workspace->macro_exists(name)) return true;
      return false;
    }

    void add_macro(const std::string &name, const std::string &expr)
    { macros[name] = expr; }
   
    const std::string& get_macro(const std::string &name) const {
      std::map<std::string, std::string>::const_iterator it=macros.find(name);
      if (it != macros.end()) return it->second;
      if (md && md->macro_exists(name)) return md->get_macro(name);
      if (parent_workspace &&
          parent_workspace->macro_exists(name))
        return parent_workspace->get_macro(name);
      GMM_ASSERT1(false, "Undefined macro");
    }

    ga_tree &macro_tree(const std::string &name, size_type meshdim,
                        bool ignore_X) const;

    void add_interpolate_transformation(const std::string &name,
                                        pinterpolate_transformation ptrans)
    { transformations[name] = ptrans; }

    bool interpolate_transformation_exists(const std::string &name) const {
      return (md && md->interpolate_transformation_exists(name)) ||
        (parent_workspace &&
         parent_workspace->interpolate_transformation_exists(name)) ||
        (transformations.find(name) != transformations.end());
    }

    pinterpolate_transformation
    interpolate_transformation(const std::string &name) const {
      std::map<std::string, pinterpolate_transformation>::const_iterator
        it = transformations.find(name);
      if (it != transformations.end()) return it->second;
      if (md && md->interpolate_transformation_exists(name))
        return md->interpolate_transformation(name);
      if (parent_workspace &&
         parent_workspace->interpolate_transformation_exists(name))
        return parent_workspace->interpolate_transformation(name);
      GMM_ASSERT1(false, "Inexistent transformation " << name);
    }

    bool is_constant(const std::string &name) const {
      VAR_SET::const_iterator it = variables.find(name);
      if (it != variables.end()) return !(it->second.is_variable);
      if (variable_group_exists(name))
        return is_constant(first_variable_of_group(name));
      if (md && md->variable_exists(name)) return md->is_data(name);
      if (parent_workspace && parent_workspace->variable_exists(name))
        return parent_workspace->is_constant(name);
      GMM_ASSERT1(false, "Undefined variable " << name);
    }

    const gmm::sub_interval &
    interval_of_variable(const std::string &name) const {
      VAR_SET::const_iterator it = variables.find(name);
      if (it != variables.end()) return it->second.I;
      if (md && md->variable_exists(name))
        return md->interval_of_variable(name);
      if (parent_workspace && parent_workspace->variable_exists(name))
        return parent_workspace->interval_of_variable(name);
      GMM_ASSERT1(false, "Undefined variable " << name);
    }

    const mesh_fem *associated_mf(const std::string &name) const {
      VAR_SET::const_iterator it = variables.find(name);
      if (it != variables.end())
        return it->second.is_fem_dofs ? it->second.mf : 0;
      if (md && md->variable_exists(name))
        return md->pmesh_fem_of_variable(name);
      if (parent_workspace && parent_workspace->variable_exists(name))
        return parent_workspace->associated_mf(name);
      if (variable_group_exists(name))
        return associated_mf(first_variable_of_group(name));
      GMM_ASSERT1(false, "Undefined variable or group " << name);
    }

    const im_data *associated_im_data(const std::string &name) const {
      VAR_SET::const_iterator it = variables.find(name);
      if (it != variables.end())  return it->second.imd;
      if (md && md->variable_exists(name))
        return md->pim_data_of_variable(name);
      if (parent_workspace && parent_workspace->variable_exists(name))
        return parent_workspace->associated_im_data(name);
      if (variable_group_exists(name)) return 0;
      GMM_ASSERT1(false, "Undefined variable " << name);
    }

    size_type qdim(const std::string &name) const {
      const mesh_fem *mf = associated_mf(name);
      const im_data *imd = associated_im_data(name);
      size_type n = gmm::vect_size(value(name));
      if (mf) {
        size_type ndof = mf->nb_dof();
        GMM_ASSERT1(ndof, "Variable " << name << " with no dof. You probably "
                    "make a wrong initialization of a mesh_fem object");
        return mf->get_qdim() * (n / ndof);
      } else if (imd) {
        size_type q = n / imd->nb_filtered_index();
        GMM_ASSERT1(q % imd->nb_tensor_elem() == 0,
                    "Invalid mesh im data vector");
        return q;
      }
      return n;
    }

    bgeot::multi_index qdims(const std::string &name) const {
      const mesh_fem *mf = associated_mf(name);
      const im_data *imd = associated_im_data(name);
      size_type n = gmm::vect_size(value(name));
      if (mf) {
        size_type ndof = mf->nb_dof();
        GMM_ASSERT1(ndof, "Variable " << name << " with no dof. You probably "
                    "make a wrong initialization of a mesh_fem object");
        bgeot::multi_index mi = mf->get_qdims();
        size_type qmult = n / ndof;
        if (qmult > 1) {
          if (mi.back() == 1) mi.back() *= qmult; else mi.push_back(qmult);
        }
        return mi;
      } else if (imd) {
        bgeot::multi_index mi = imd->tensor_size();
        size_type q = n / imd->nb_filtered_index();
        GMM_ASSERT1(q % imd->nb_tensor_elem() == 0,
                    "Invalid mesh im data vector");
        size_type qmult = q / imd->nb_tensor_elem();
        if (qmult > 1) {
          if (mi.back() == 1) mi.back() *= qmult; else mi.push_back(qmult);
        }
        return mi;
      }
      bgeot::multi_index mi(1); mi[0] = n; return mi;
    }

    const model_real_plain_vector &value(const std::string &name) const {
      VAR_SET::const_iterator it = variables.find(name);
      if (it != variables.end())
        return *(it->second.V);
      if (md && md->variable_exists(name))
        return md->real_variable(name);
      if (parent_workspace && parent_workspace->variable_exists(name))
        return parent_workspace->value(name);
      if (variable_group_exists(name))
        return value(first_variable_of_group(name));
      GMM_ASSERT1(false, "Undefined variable or group " << name);
    }

    void assembly(size_type order);

    ga_workspace(const getfem::model &md_) : md(&md_), parent_workspace(0) {}
    ga_workspace(const ga_workspace &gaw) : md(0), parent_workspace(&gaw) {}
    ga_workspace(void) : md(0), parent_workspace(0) {}
    ~ga_workspace() { clear_expressions(); }

  };

  //=========================================================================
  // Intermediate structure for user function manipulation
  //=========================================================================

  struct ga_instruction_set;

  class ga_function {
    // gerer les arguments et leur modif éventuelle ..
    // Il faut pouvoir deriver les fonctions et obtenir leur gradients
    mutable ga_workspace local_workspace;
    std::string expr;
    mutable ga_instruction_set *gis;
    
  public:
    ga_function(void) : gis(0) {}
    ga_function(const model &md, const std::string &e);
    ga_function(const ga_workspace &workspace_, const std::string &e);
    ga_function(const std::string &e);
    ga_function(const ga_function &gaf);
    ga_function &operator =(const ga_function &gaf);
    ~ga_function();
    const std::string &expression(void) const { return expr; }
    const base_tensor &eval(void) const;
    void derivative(const std::string &variable);
    void compile(void) const;
    ga_workspace &workspace(void) const { return  local_workspace; }

  };

  //=========================================================================
  // Intermediate structure for interpolation functions
  //=========================================================================

  struct ga_interpolation_context {

    virtual const bgeot::stored_point_tab &
    points_for_element(size_type cv, short_type f,
                       std::vector<size_type> &ind) const = 0;
    virtual bool use_pgp(size_type cv) const = 0;
    virtual void store_result(size_type cv, size_type i, base_tensor &t) = 0;
    virtual void finalize(void) = 0;
    virtual const mesh &linked_mesh(void) = 0;
    virtual ~ga_interpolation_context() {}
  };
  

  //=========================================================================
  // Interpolation functions
  //=========================================================================

  
  void ga_interpolation(ga_workspace &workspace,
                        ga_interpolation_context &gic);

  void ga_interpolation_Lagrange_fem
  (ga_workspace &workspace, const mesh_fem &mf, base_vector &result);
  
  void ga_interpolation_Lagrange_fem
  (const getfem::model &md, const std::string &expr, const mesh_fem &mf,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes());

  // Not tested
  void ga_interpolation_mti
  (const getfem::model &md, const std::string &expr, mesh_trans_inv &mti,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes(),
   int extrapolation = 0, const mesh_region &rg_source=mesh_region::all_convexes(),
   size_type nbdof_ = size_type(-1));

  // Not tested
  void ga_interpolation_im_data
  (const getfem::model &md, const std::string &expr, im_data &imd,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes());

  //=========================================================================
  // Interpolation transformation
  //=========================================================================

  /** Add a transformation to the model `md` from mesh `source_mesh` to mesh
      `target_mesh` given by the expression `expr` which corresponds to a
      high-level generic assembly expression which may contains some
      variable of the model. CAUTION: For the moment, the derivative of the
      transformation with respect to the eventual variables used is not
      taken into account in the model solve. 
  */
  void add_interpolate_transformation_from_expression
  (model &md, const std::string &transname, const mesh &source_mesh,
   const mesh &target_mesh, const std::string &expr);

  void add_interpolate_transformation_from_expression
  (ga_workspace &workspace, const std::string &transname,
   const mesh &source_mesh, const mesh &target_mesh, const std::string &expr);
  
}  /* end of namespace getfem.                                             */


#endif /* GETFEM_GENERIC_ASSEMBLY_H__  */
