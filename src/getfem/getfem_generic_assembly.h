/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2013-2016 Yves Renard

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
#if defined(INFINITY)
#undef INFINITY
#endif
#define INFINITY std::numeric_limits<scalar_type>::infinity()
#endif

namespace getfem {

  class ga_tree;

  int ga_check_name_validity(const std::string &name);

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
    typedef std::map<std::string, std::shared_ptr<ga_nonlinear_operator>> T;
    T tab;

    void add_method(const std::string &name,
                    const std::shared_ptr<ga_nonlinear_operator> &pt)
    { tab[name] = pt; }
  };

  //=========================================================================
  // For user predefined scalar functions.
  //=========================================================================

  typedef scalar_type (*pscalar_func_onearg)(scalar_type);
  typedef scalar_type (*pscalar_func_twoargs)(scalar_type, scalar_type);

  void ga_define_function(const std::string name, size_type nb_args,
                          const std::string expr, const std::string der1="",
                          const std::string der2="");
  void ga_define_function(const std::string name, pscalar_func_onearg f,
                          const std::string &der1="");
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
    bool enable_all_md_variables;

    void init();

    struct var_description {

      bool is_variable;
      bool is_fem_dofs;
      const mesh_fem *mf;
      gmm::sub_interval I;
      const model_real_plain_vector *V;
      const im_data *imd;
      bgeot::multi_index qdims;  // For data having a qdim != of the fem
                                 // (dim per dof for dof data)
                                 // and for constant variables.

      size_type qdim() const {
        size_type q = 1;
        for (size_type i = 0; i < qdims.size(); ++i) q *= qdims[i];
        return q;
      }

      var_description(bool is_var, bool is_fem,
                      const mesh_fem *mmf, gmm::sub_interval I_,
                      const model_real_plain_vector *v, const im_data *imd_,
                      size_type Q)
        : is_variable(is_var), is_fem_dofs(is_fem), mf(mmf), I(I_), V(v),
          imd(imd_), qdims(1) {
        GMM_ASSERT1(Q > 0, "Bad dimension");
        qdims[0] = Q;
      }
      var_description() : is_variable(false), is_fem_dofs(false),
                          mf(0), V(0), imd(0), qdims(1) { qdims[0] = 1; }
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
      tree_description()
        : name_test1(""), name_test2(""),
          interpolate_name_test1(""), interpolate_name_test2(""),
          mim(0), m(0), rg(0), ptree(0), elem(0) {}
      void copy(const tree_description& td);
      tree_description(const tree_description& td) { copy(td); }
      tree_description &operator =(const tree_description& td);
      ~tree_description();
    };

    mutable std::set<var_trans_pair> test1, test2;
    var_trans_pair selected_test1, selected_test2;

  private:

    // mesh regions
    std::map<const mesh *, std::list<mesh_region> > registred_mesh_regions;

    const mesh_region &
    register_region(const mesh &m, const mesh_region &region);

    // variables and variable groups
    mutable std::map<std::string, gmm::sub_interval> int_disabled_variables;

    typedef std::map<std::string, var_description> VAR_SET;
    VAR_SET variables;
    std::map<std::string, pinterpolate_transformation> transformations;
    std::map<std::string, pelementary_transformation> elem_transformations;
    std::vector<tree_description> trees;

    std::map<std::string, std::vector<std::string> > variable_groups;
    std::map<std::string, std::string> macros;

    struct m_tree {
      ga_tree *ptree;
      size_type meshdim;
      bool ignore_X;
      m_tree() : ptree(0), meshdim(-1), ignore_X(false) {}
      m_tree(const m_tree& o);
      m_tree &operator =(const m_tree& o);
      ~m_tree();
    };

    mutable std::map<std::string, m_tree> macro_trees;

    void add_tree(ga_tree &tree, const mesh &m, const mesh_im &mim,
                  const mesh_region &rg,
                  const std::string &expr, size_type add_derivative_order = 2,
                  bool scalar_expr = true);


    std::shared_ptr<model_real_sparse_matrix> K;
    model_real_sparse_matrix unreduced_K;
    std::shared_ptr<base_vector> V;
    base_vector unreduced_V;
    scalar_type E = scalar_type(0.);
    base_tensor assemb_t;

  public:

    const model_real_sparse_matrix &assembled_matrix() const { return *K;}
    model_real_sparse_matrix &assembled_matrix() { return *K; }
    scalar_type &assembled_potential() { return E; }
    const scalar_type &assembled_potential() const { return E; }
    const base_vector &assembled_vector() const { return *V; }
    base_vector &assembled_vector() { return *V; }
    void set_assembled_matrix(model_real_sparse_matrix &K_) {
      K = std::shared_ptr<model_real_sparse_matrix>
          (std::shared_ptr<model_real_sparse_matrix>(), &K_);
    }
    void set_assembled_vector(base_vector &V_) {
      V = std::shared_ptr<base_vector>
          (std::shared_ptr<base_vector>(), &V_);
    }
    base_tensor &assembled_tensor() { return assemb_t; }
    const base_tensor &assembled_tensor() const { return assemb_t; }

    model_real_sparse_matrix &unreduced_matrix()
    { return unreduced_K; }
    base_vector &unreduced_vector() { return unreduced_V; }

    /** Add an expression, perform the semantic analysis, split into
     *  terms in separated test functions, derive if necessary to obtain
     *  the tangent terms. Return the maximal order found in the expression.
     */
    size_type add_expression(const std::string expr, const mesh_im &mim,
                             const mesh_region &rg=mesh_region::all_convexes(),
                             size_type add_derivative_order = 2);
    /* Internal use */
    void add_function_expression(const std::string expr);
    /* Internal use */
    void add_interpolation_expression(const std::string expr, const mesh &m,
                                      const mesh_region &rg=mesh_region::all_convexes());
    void add_interpolation_expression(const std::string expr, const mesh_im &mim,
                                      const mesh_region &rg=mesh_region::all_convexes());

    /** Delete all previously added expressions. */
    void clear_expressions();

    /** Print some information about all previously added expressions. */
    void print(std::ostream &str);

    size_type nb_trees() const;
    tree_description &tree_info(size_type i);

    // variables and variable groups
    void add_fem_variable(const std::string &name, const mesh_fem &mf,
                          const gmm::sub_interval &I,
                          const model_real_plain_vector &VV);
    void add_fixed_size_variable(const std::string &name,
                                 const gmm::sub_interval &I,
                                 const model_real_plain_vector &VV);
    void add_fem_constant(const std::string &name, const mesh_fem &mf,
                          const model_real_plain_vector &VV);
    void add_fixed_size_constant(const std::string &name,
                                 const model_real_plain_vector &VV);
    void add_im_data(const std::string &name, const im_data &imd,
                     const model_real_plain_vector &VV);

    bool used_variables(model::varnamelist &vl, model::varnamelist &vl_test1,
                        model::varnamelist &vl_test2, model::varnamelist &dl,
                        size_type order);

    bool variable_exists(const std::string &name) const;

    const std::string &variable_in_group(const std::string &group_name,
                                         const mesh &m) const;

    void define_variable_group(const std::string &group_name,
                               const std::vector<std::string> &nl);

    bool variable_group_exists(std::string name) const;

    bool variable_or_group_exists(const std::string &name) const
    { return variable_exists(name) || variable_group_exists(name); }

    const std::vector<std::string> &
    variable_group(const std::string &group_name) const;

    const std::string& first_variable_of_group(const std::string &name) const;

    bool is_constant(const std::string &name) const;

    bool is_disabled_variable(const std::string &name) const;

    const scalar_type &factor_of_variable(const std::string &name) const;

    const gmm::sub_interval &
    interval_of_disabled_variable(const std::string &name) const;

    const gmm::sub_interval &
    interval_of_variable(const std::string &name) const;

    const mesh_fem *associated_mf(const std::string &name) const;

    const im_data *associated_im_data(const std::string &name) const;

    size_type qdim(const std::string &name) const;

    bgeot::multi_index qdims(const std::string &name) const;

    const model_real_plain_vector &value(const std::string &name) const;


    // macros
    bool macro_exists(const std::string &name) const;

    void add_macro(const std::string &name, const std::string &expr)
    { macros[name] = expr; }

    const std::string& get_macro(const std::string &name) const;

    ga_tree& macro_tree(const std::string &name, size_type meshdim,
                        size_type ref_elt_dim, bool ignore_X) const;


    // interpolate and elementary transformations
    void add_interpolate_transformation(const std::string &name,
                                        pinterpolate_transformation ptrans);

    bool interpolate_transformation_exists(const std::string &name) const;

    pinterpolate_transformation
    interpolate_transformation(const std::string &name) const;

    void add_elementary_transformation(const std::string &name,
                                       pelementary_transformation ptrans)
    { elem_transformations[name] = ptrans; }

    bool elementary_transformation_exists(const std::string &name) const;

    pelementary_transformation
    elementary_transformation(const std::string &name) const;


    // extract terms
    std::string extract_constant_term(const mesh &m);
    std::string extract_order1_term(const std::string &varname);
    std::string extract_order0_term();
    std::string extract_Neumann_term(const std::string &varname);


    void assembly(size_type order);


    ga_workspace(const getfem::model &md_, bool enable_all_variables = false)
      : md(&md_), parent_workspace(0),
        enable_all_md_variables(enable_all_variables)
    { init(); }
    ga_workspace(bool, const ga_workspace &gaw)
      : md(0), parent_workspace(&gaw), enable_all_md_variables(false)
    { init(); }
    ga_workspace()
      : md(0), parent_workspace(0), enable_all_md_variables(false)
    { init(); }
    ~ga_workspace() { clear_expressions(); }

  };

  //=========================================================================
  // Intermediate structure for user function manipulation
  //=========================================================================

  struct ga_instruction_set;

  class ga_function {
    mutable ga_workspace local_workspace;
    std::string expr;
    mutable ga_instruction_set *gis;

  public:
    ga_function() : local_workspace(), expr(""), gis(0) {}
    ga_function(const model &md, const std::string &e);
    ga_function(const ga_workspace &workspace_, const std::string &e);
    ga_function(const std::string &e);
    ga_function(const ga_function &gaf);
    ga_function &operator =(const ga_function &gaf);
    ~ga_function();
    const std::string &expression() const { return expr; }
    const base_tensor &eval() const;
    void derivative(const std::string &variable);
    void compile() const;
    ga_workspace &workspace() const { return  local_workspace; }

  };

  //=========================================================================
  // Intermediate structure for interpolation functions
  //=========================================================================

  struct ga_interpolation_context {

    virtual bgeot::pstored_point_tab
    ppoints_for_element(size_type cv, short_type f,
                        std::vector<size_type> &ind) const = 0;
    inline const bgeot::stored_point_tab &points_for_element
    (size_type cv, short_type f, std::vector<size_type> &ind) const
    { return *ppoints_for_element(cv, f, ind); }
    virtual bool use_pgp(size_type cv) const = 0;
    virtual bool use_mim() const = 0;
    virtual void store_result(size_type cv, size_type i, base_tensor &t) = 0;
    virtual void finalize() = 0;
    virtual const mesh &linked_mesh() = 0;
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

  void ga_interpolation_mti
  (const getfem::model &md, const std::string &expr, mesh_trans_inv &mti,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes(),
   int extrapolation = 0,
   const mesh_region &rg_source=mesh_region::all_convexes(),
   size_type nbdof_ = size_type(-1));

  void ga_interpolation_im_data
  (const getfem::model &md, const std::string &expr, const im_data &imd,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes());

  void ga_interpolation_im_data
  (ga_workspace &workspace, const im_data &imd, base_vector &result,
   const mesh_region &rg=mesh_region::all_convexes());

  //=========================================================================
  // Local projection functions
  //=========================================================================

  /** Make an elementwise L2 projection of an expression with respect
      to the mesh_fem `mf`. This mesh_fem has to be a discontinuous one.
      The expression has to be valid according to the high-level generic
      assembly language possibly including references to the variables
      and data of the model. 
  */
  void ga_local_projection(const getfem::model &md, const mesh_im &mim,
			   const std::string &expr, const mesh_fem &mf,
			   base_vector &result,
			   const mesh_region &rg=mesh_region::all_convexes());

  //=========================================================================
  // Interpolate transformations
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


  /** Create a new instance of a transformation corresponding to the
      interpolation on the neighbour element. Can only be applied to the
      computation on some internal faces of a mesh.
      (mainly for internal use in the constructor of getfem::model)
  */
  pinterpolate_transformation interpolate_transformation_neighbour_instance();

}  /* end of namespace getfem.                                             */


#endif /* GETFEM_GENERIC_ASSEMBLY_H__  */
