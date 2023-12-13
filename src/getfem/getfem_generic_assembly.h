/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2013-2020 Yves Renard

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

/** @file getfem_generic_assembly.h
    @author  Yves Renard <Yves.Renard@insa-lyon.fr>
    @date November 18, 2013.
    @brief A language for generic assembly of pde boundary value problems.
 */


#ifndef GETFEM_GENERIC_ASSEMBLY_H__
#define GETFEM_GENERIC_ASSEMBLY_H__

#include <map>
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_mesh_slice.h"


#ifdef _WIN32
#include <limits>
#if defined(INFINITY)
#undef INFINITY
#endif
#define INFINITY std::numeric_limits<scalar_type>::infinity()
#endif

namespace getfem {

  struct ga_tree;
  class model;
  class ga_workspace;

  typedef gmm::rsvector<scalar_type> model_real_sparse_vector;
  typedef gmm::rsvector<complex_type> model_complex_sparse_vector;
  typedef std::vector<scalar_type> model_real_plain_vector;
  typedef std::vector<complex_type> model_complex_plain_vector;

  typedef gmm::col_matrix<model_real_sparse_vector> model_real_sparse_matrix;
  typedef gmm::col_matrix<model_complex_sparse_vector>
  model_complex_sparse_matrix;

  typedef gmm::row_matrix<model_real_sparse_vector>
  model_real_row_sparse_matrix;
  typedef gmm::row_matrix<model_complex_sparse_vector>
  model_complex_row_sparse_matrix;
  
  // 0 : ok
  // 1 : function or operator name or "X"
  // 2 : reserved prefix Grad, Hess, Div, Test and Test2
  int ga_check_name_validity(const std::string &name);

  //=========================================================================
  //  Virtual interpolate_transformation object.
  //=========================================================================

  struct var_trans_pair {
    std::string varname, transname;
    bool operator <(const var_trans_pair &vt) const {
      return (varname < vt.varname) ||
             (!(varname > vt.varname) && transname < vt.transname);
    }
    var_trans_pair() : varname(), transname() {}
    var_trans_pair(const std::string &v, const std::string &t)
      : varname(v), transname(t) {}
  };

  class APIDECL virtual_interpolate_transformation {

  public:
    virtual void extract_variables
    (const ga_workspace &workspace, std::set<var_trans_pair> &vars,
     bool ignore_data, const mesh &m,
     const std::string &interpolate_name) const = 0;
    virtual void init(const ga_workspace &workspace) const = 0;
    virtual int transform
    (const ga_workspace &workspace, const mesh &m,
     fem_interpolation_context &ctx_x, const base_small_vector &Normal,
     const mesh **m_t, size_type &cv, short_type &face_num,
     base_node &P_ref, base_small_vector &N_y,
     std::map<var_trans_pair, base_tensor> &derivatives,
     bool compute_derivatives) const = 0;
    virtual void finalize() const = 0;
    virtual std::string expression() const { return std::string(); }

    virtual ~virtual_interpolate_transformation() {}
  };

  typedef std::shared_ptr<const virtual_interpolate_transformation>
  pinterpolate_transformation;

  //=========================================================================
  //  Virtual elementary_transformation object.
  //=========================================================================

  class APIDECL virtual_elementary_transformation {

  public:
    
    virtual void give_transformation(const mesh_fem &mf1, const mesh_fem &mf2,
                                     size_type cv, base_matrix &M) const = 0;
    virtual ~virtual_elementary_transformation() {}
  };

  typedef std::shared_ptr<const virtual_elementary_transformation>
  pelementary_transformation;

  //=========================================================================
  //  Virtual secondary_domain object.
  //=========================================================================

  class APIDECL virtual_secondary_domain {
  protected:
    const mesh_im &mim_;
    const mesh_region region;

  public:

    const mesh_im &mim() const { return mim_; }
    virtual const mesh_region &give_region(const mesh &m,
                                     size_type cv, short_type f) const = 0;
    // virtual void init(const ga_workspace &workspace) const = 0;
    // virtual void finalize() const = 0;

    virtual_secondary_domain(const mesh_im &mim__, const mesh_region &region_)
      : mim_(mim__), region(region_) {}
    virtual ~virtual_secondary_domain() {}
  };

  typedef std::shared_ptr<const virtual_secondary_domain> psecondary_domain;

  //=========================================================================
  // Structure dealing with macros.
  //=========================================================================

  class ga_macro {

  protected:
    ga_tree *ptree;
    std::string macro_name_;
    size_type nbp;

  public:
    ga_macro();
    ga_macro(const std::string &name, const ga_tree &t, size_type nbp_);
    ga_macro(const ga_macro &);
    ~ga_macro();
    ga_macro &operator =(const ga_macro &);

    const std::string &name() const { return macro_name_; }
    std::string &name() { return macro_name_; }
    size_type nb_params() const { return nbp; }
    size_type &nb_params() { return nbp; }
    const ga_tree& tree() const { return *ptree; }
    ga_tree& tree() { return *ptree; }
  };


  class ga_macro_dictionary {

  protected:
    const ga_macro_dictionary *parent;
    std::map<std::string, ga_macro> macros;

  public:
    bool macro_exists(const std::string &name) const;
    const ga_macro &get_macro(const std::string &name) const;
    
    void add_macro(const ga_macro &gam);
    void add_macro(const std::string &name, const std::string &expr);
    void del_macro(const std::string &name);
    
    ga_macro_dictionary() : parent(0) {}
    ga_macro_dictionary(bool, const ga_macro_dictionary& gamd)
      : parent(&gamd) {}
    
  };

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
    ga_predef_operator_tab();
  };

  //=========================================================================
  // For user predefined scalar functions.
  //=========================================================================

  typedef scalar_type (*pscalar_func_onearg)(scalar_type);
  typedef scalar_type (*pscalar_func_twoargs)(scalar_type, scalar_type);

  void ga_define_function(const std::string &name, size_type nb_args,
                          const std::string &expr, const std::string &der1="",
                          const std::string &der2="");
  void ga_define_function(const std::string &name, pscalar_func_onearg f,
                          const std::string &der1="");
  void ga_define_function(const std::string &name, pscalar_func_twoargs f2,
                          const std::string &der1="",
                          const std::string &der2="");

  void ga_undefine_function(const std::string &name);
  bool ga_function_exists(const std::string &name);

  //=========================================================================
  // Structure dealing with user defined environment : constant, variables,
  // functions, operators.
  //=========================================================================

  class ga_workspace {

    const model *md;
    const ga_workspace *parent_workspace;
    bool with_parent_variables;
    size_type nb_prim_dof, nb_intern_dof, first_intern_dof, nb_tmp_dof;

    void init();

    struct var_description {

      bool is_variable;
      const mesh_fem *mf;
      const im_data *imd;
      gmm::sub_interval I;
      const model_real_plain_vector *V;
      bgeot::multi_index qdims;  // For data having a qdim different than the
                                 // qdim of the fem or im_data (dim per dof for
                                 // dof data) and for constant variables.
      bool is_internal;

      size_type qdim() const {
        size_type q = 1;
        for (size_type i = 0; i < qdims.size(); ++i) q *= qdims[i];
        return q;
      }

      var_description(bool is_var, const mesh_fem *mf_, const im_data *imd_,
                      gmm::sub_interval I_, const model_real_plain_vector *V_,
                      size_type Q, bool is_intern_=false)
        : is_variable(is_var), mf(mf_), imd(imd_),
          I(I_), V(V_), qdims(1), is_internal(is_intern_)
      {
        GMM_ASSERT1(Q > 0, "Bad dimension");
        qdims[0] = Q;
      }
    };

  public:

    enum operation_type {ASSEMBLY,
                         PRE_ASSIGNMENT,
                         POST_ASSIGNMENT};

    struct tree_description { // CAUTION: Specific copy constructor
      size_type order; //  0 : potential
                       //  1 : residual
                       //  2 : tangent operator
                       // -1 : any
      operation_type operation;
      std::string varname_interpolation; // Where to interpolate
      std::string name_test1, name_test2;
      std::string interpolate_name_test1, interpolate_name_test2;
      std::string secondary_domain;
      const mesh_im *mim;
      const mesh *m;
      const mesh_region *rg;
      ga_tree *ptree;
      tree_description()
        : operation(ASSEMBLY), varname_interpolation(""),
          name_test1(""), name_test2(""),
          interpolate_name_test1(""), interpolate_name_test2(""),
          mim(0), m(0), rg(0), ptree(0) {}
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

    const mesh_region &register_region(const mesh &m, const mesh_region &rg);

    // variables and variable groups
    typedef std::map<std::string, var_description> VAR_SET;
    VAR_SET variables;

    std::map<std::string, gmm::sub_interval> reenabled_var_intervals,
                                             tmp_var_intervals;

    std::map<std::string, pinterpolate_transformation> transformations;
    std::map<std::string, pelementary_transformation> elem_transformations;
    std::map<std::string, psecondary_domain> secondary_domains;
    
    std::vector<tree_description> trees;

    std::map<std::string, std::vector<std::string> > variable_groups;

    ga_macro_dictionary macro_dict;

    // Adds a tree to the workspace. The argument tree is consumed by the
    // function and cannot be reused afterwards.
    void add_tree(ga_tree &tree, const mesh &m, const mesh_im &mim,
                  const mesh_region &rg,
                  const std::string &expr, size_type add_derivative_order,
                  bool scalar_expr, operation_type op_type=ASSEMBLY,
                  const std::string varname_interpolation="");

    std::shared_ptr<model_real_sparse_matrix> K, KQJpr;
    std::shared_ptr<base_vector> V; // reduced residual vector (primary vars + internal vars)
                                    // after condensation it partially holds the condensed residual
                                    // and the internal solution
    model_real_sparse_matrix col_unreduced_K,
                             row_unreduced_K,
                             row_col_unreduced_K;
    base_vector unreduced_V, cached_V;
    base_tensor assemb_t;
    bool include_empty_int_pts = false;

  public:
    // setter functions
    void set_assembled_matrix(model_real_sparse_matrix &K_) {
      K = std::shared_ptr<model_real_sparse_matrix>
          (std::shared_ptr<model_real_sparse_matrix>(), &K_); // alias
    }
    void set_assembled_vector(base_vector &V_) {
      V = std::shared_ptr<base_vector>
          (std::shared_ptr<base_vector>(), &V_); // alias
    }
    // getter functions
    const model_real_sparse_matrix &assembled_matrix() const { return *K; }
    model_real_sparse_matrix &assembled_matrix() { return *K; }
    const base_vector &assembled_vector() const { return *V; }
    base_vector &assembled_vector() { return *V; }
    const base_vector &cached_vector() const { return cached_V; }
    const base_tensor &assembled_tensor() const { return assemb_t; }
    base_tensor &assembled_tensor() { return assemb_t; }
    const scalar_type &assembled_potential() const {
      GMM_ASSERT1(assemb_t.size() == 1, "Bad result size");
      return assemb_t[0];
    }
    scalar_type &assembled_potential() {
      GMM_ASSERT1(assemb_t.size() == 1, "Bad result size");
      return assemb_t[0];
    }
    model_real_sparse_matrix &row_unreduced_matrix()
    { return row_unreduced_K; }
    model_real_sparse_matrix &col_unreduced_matrix()
    { return col_unreduced_K; }
    model_real_sparse_matrix &row_col_unreduced_matrix()
    { return row_col_unreduced_K; }
    base_vector &unreduced_vector() { return unreduced_V; }
    // setter function for condensation matrix
    void set_internal_coupling_matrix(model_real_sparse_matrix &KQJpr_) {
      KQJpr = std::shared_ptr<model_real_sparse_matrix>
              (std::shared_ptr<model_real_sparse_matrix>(), &KQJpr_); // alias
    }
    /** getter function for condensation matrix
        Its size is (nb_primary_dof()+nb_internal_dof()) x nb_primary_dof()
        but its first nb_primary_dof() rows are empty */
    const model_real_sparse_matrix &internal_coupling_matrix() const
    { return *KQJpr; }
    model_real_sparse_matrix &internal_coupling_matrix() { return *KQJpr; }

    /** Add an expression, perform the semantic analysis, split into
     *  terms in separated test functions, derive if necessary to obtain
     *  the tangent terms. Return the maximal order found in the expression.
     */
    size_type add_expression(const std::string &expr, const mesh_im &mim,
                             const mesh_region &rg=mesh_region::all_convexes(),
                             size_type add_derivative_order = 2,
                             const std::string &secondary_dom = "");
    /* Internal use */
    void add_function_expression(const std::string &expr);
    /* Internal use */
    void add_interpolation_expression
    (const std::string &expr, const mesh &m,
     const mesh_region &rg = mesh_region::all_convexes());
    void add_interpolation_expression
    (const std::string &expr, const mesh_im &mim,
     const mesh_region &rg = mesh_region::all_convexes());
    void add_assignment_expression
    (const std::string &dataname, const std::string &expr,
     const mesh_region &rg_ = mesh_region::all_convexes(),
     size_type order = 1, bool before = false);

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
    void add_im_variable(const std::string &name, const im_data &imd,
                         const gmm::sub_interval &I,
                         const model_real_plain_vector &VV);
    void add_internal_im_variable(const std::string &name, const im_data &imd,
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

    bool used_variables(std::vector<std::string> &vl,
                        std::vector<std::string> &vl_test1,
                        std::vector<std::string> &vl_test2,
                        std::vector<std::string> &dl,
                        size_type order);
    bool is_linear(size_type order);

    bool variable_exists(const std::string &name) const;

    bool is_internal_variable(const std::string &name) const;

    const std::string &variable_in_group(const std::string &group_name,
                                         const mesh &m) const;

    void define_variable_group(const std::string &group_name,
                               const std::vector<std::string> &nl);

    bool variable_group_exists(const std::string &name) const;

    bool variable_or_group_exists(const std::string &name) const
    { return variable_exists(name) || variable_group_exists(name); }

    const std::vector<std::string> &
    variable_group(const std::string &group_name) const;

    const std::string& first_variable_of_group(const std::string &name) const;

    bool is_constant(const std::string &name) const;

    bool is_disabled_variable(const std::string &name) const;

    const scalar_type &factor_of_variable(const std::string &name) const;

    const gmm::sub_interval &
    interval_of_variable(const std::string &name) const;

    const mesh_fem *associated_mf(const std::string &name) const;

    const im_data *associated_im_data(const std::string &name) const;

    size_type qdim(const std::string &name) const;

    bgeot::multi_index qdims(const std::string &name) const;

    const model_real_plain_vector &value(const std::string &name) const;
    scalar_type get_time_step() const;

    // macros
    bool macro_exists(const std::string &name) const
    { return macro_dict.macro_exists(name); }

    void add_macro(const std::string &name, const std::string &expr)
    { macro_dict.add_macro(name, expr); }

    void del_macro(const std::string &name) { macro_dict.del_macro(name); }

    const std::string& get_macro(const std::string &name) const;

    const ga_macro_dictionary &macro_dictionary() const { return macro_dict; }


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

    void add_secondary_domain(const std::string &name,
                              psecondary_domain psecdom);

    bool secondary_domain_exists(const std::string &name) const;

    psecondary_domain secondary_domain(const std::string &name) const;

    // extract terms
    std::string extract_constant_term(const mesh &m);
    std::string extract_order1_term(const std::string &varname);
    std::string extract_order0_term();
    std::string extract_Neumann_term(const std::string &varname);

    void assembly(size_type order, bool condensation=false);

    void set_include_empty_int_points(bool include);
    bool include_empty_int_points() const;

    size_type nb_primary_dof() const { return nb_prim_dof; }
    size_type nb_internal_dof() const { return nb_intern_dof; }
    size_type first_internal_dof() const { return first_intern_dof; }
    size_type nb_temporary_dof() const { return nb_tmp_dof; }

    void add_temporary_interval_for_unreduced_variable(const std::string &name);

    void clear_temporary_variable_intervals() {
      tmp_var_intervals.clear();
      nb_tmp_dof = 0;
    }

    const gmm::sub_interval &
    temporary_interval_of_variable(const std::string &name) const {
      static const gmm::sub_interval empty_interval;
      const auto it = tmp_var_intervals.find(name);
      return (it != tmp_var_intervals.end()) ? it->second : empty_interval;
    }

    enum class inherit { NONE, ENABLED, ALL };

    explicit ga_workspace(const getfem::model &md_,
                          const inherit var_inherit=inherit::ENABLED);
    explicit ga_workspace(const ga_workspace &gaw,    // compulsory 2nd arg to avoid
                          const inherit var_inherit); // conflict with copy constructor
    ga_workspace();
    ~ga_workspace();

  };

  // Small tool to make basic substitutions into an assembly string
  std::string ga_substitute(const std::string &expr,
                            const std::map<std::string, std::string> &dict);

  inline std::string ga_substitute(const std::string &expr,
                                  const std::string &o1,const std::string &s1) {
    std::map<std::string, std::string> dict;
    dict[o1] = s1;
    return ga_substitute(expr, dict);
  }

  inline std::string ga_substitute(const std::string &expr,
                                  const std::string &o1,const std::string &s1,
                                  const std::string &o2,const std::string &s2) {
    std::map<std::string, std::string> dict;
    dict[o1] = s1; dict[o2] = s2; 
    return ga_substitute(expr, dict);
  }

  inline std::string ga_substitute(const std::string &expr,
                                  const std::string &o1,const std::string &s1,
                                  const std::string &o2,const std::string &s2,
                                  const std::string &o3,const std::string &s3) {
    std::map<std::string, std::string> dict;
    dict[o1] = s1; dict[o2] = s2; dict[o3] = s3; 
    return ga_substitute(expr, dict);
  }

  inline std::string ga_substitute(const std::string &expr,
                                  const std::string &o1,const std::string &s1,
                                  const std::string &o2,const std::string &s2,
                                  const std::string &o3,const std::string &s3,
                                  const std::string &o4,const std::string &s4) {
    std::map<std::string, std::string> dict;
    dict[o1] = s1; dict[o2] = s2;  dict[o3] = s3; dict[o4] = s4; 
    return ga_substitute(expr, dict);
  }


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
  (ga_workspace &workspace, const im_data &imd, base_vector &result);

  void ga_interpolation_im_data
  (const getfem::model &md, const std::string &expr, const im_data &imd,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes());

  void ga_interpolation_mesh_slice
  (ga_workspace &workspace, const stored_mesh_slice &sl, base_vector &result);

  void ga_interpolation_mesh_slice
  (const getfem::model &md, const std::string &expr, const stored_mesh_slice &sl,
   base_vector &result, const mesh_region &rg=mesh_region::all_convexes());


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

  /** Add a transformation to a workspace `workspace` or a model `md` mapping
      point in mesh `source_mesh` to mesh `target_mesh`, optionally restricted
      to the region `target_region`. The transformation is defined by the
      expression `expr`, which has to be in the high-level generic assembly
      syntax and may contain some variables of the workspace/model.
      CAUTION: For the moment, the derivative of the transformation with
      respect to any of these variables is not taken into account in the model
      solve.
  */
  void add_interpolate_transformation_from_expression
  (ga_workspace &workspace, const std::string &transname,
   const mesh &source_mesh, const mesh &target_mesh, const std::string &expr);
  void add_interpolate_transformation_from_expression
  (ga_workspace &workspace, const std::string &transname,
   const mesh &source_mesh, const mesh &target_mesh,
   size_type target_region, const std::string &expr);
  void add_interpolate_transformation_from_expression
  (model &md, const std::string &transname,
   const mesh &source_mesh, const mesh &target_mesh, const std::string &expr);
  void add_interpolate_transformation_from_expression
  (model &md, const std::string &transname,
   const mesh &source_mesh, const mesh &target_mesh,
   size_type target_region, const std::string &expr);

  /** Add a transformation to the workspace that creates an identity mapping
      between two meshes in deformed state. Conceptually, it can be viewed
      as a transformation from expression Xsource + Usource - Utarget,
      except such an expression cannot be used directly in the transformation
      from expression (function above), as Utarget needs to be interpolated
      though an inversion of the transformation of the target domain.
      Thread safe if added to thread local workspace.
  */
  void add_interpolate_transformation_on_deformed_domains
  (ga_workspace &workspace, const std::string &transname,
   const mesh &source_mesh, const std::string &source_displacements,
   const mesh_region &source_region, const mesh &target_mesh,
   const std::string &target_displacements, const mesh_region &target_region);

  /** The same as above, but adding transformation to the model.
  Note, this version is not thread safe.*/
  void add_interpolate_transformation_on_deformed_domains
  (model &md, const std::string &transname,
   const mesh &source_mesh, const std::string &source_displacements,
   const mesh_region &source_region, const mesh &target_mesh,
   const std::string &target_displacements, const mesh_region &target_region);

  /** Create a new instance of a transformation corresponding to the
      interpolation on the neighbor element. Can only be applied to the
      computation on some internal faces of a mesh.
      (mainly for internal use in the constructor of getfem::model)
  */
  pinterpolate_transformation interpolate_transformation_neighbor_instance();

  /* Add a special interpolation transformation which represents the identity
     transformation but allows to evaluate the expression on another element
     than the current element by polynomial extrapolation. It is used for
     stabilization term in fictitious domain applications. the map elt_cor
     list the element concerned by the transformation and associate them
     to the element on which the extrapolation has to be made. If an element
     is not listed in elt_cor the evaluation is just made on the current
     element.
  */
  void add_element_extrapolation_transformation
  (model &md, const std::string &name, const mesh &sm,
   std::map<size_type, size_type> &elt_corr);

  void add_element_extrapolation_transformation
  (ga_workspace &workspace, const std::string &name, const mesh &sm,
   std::map<size_type, size_type> &elt_corr);
  
  /* Change the correspondence map of an element extrapolation interpolate
     transformation.
  */
  void set_element_extrapolation_correspondence
  (model &md, const std::string &name,
   std::map<size_type, size_type> &elt_corr);
  
  void set_element_extrapolation_correspondence
  (ga_workspace &workspace, const std::string &name,
   std::map<size_type, size_type> &elt_corr);
    
  //=========================================================================
  // Secondary domains
  //=========================================================================

  void add_standard_secondary_domain
  (model &md, const std::string &name, const mesh_im &mim,
   const mesh_region &rg=mesh_region::all_convexes());
  
  void add_standard_secondary_domain
  (ga_workspace &workspace, const std::string &name, const mesh_im &mim,
   const mesh_region &rg=mesh_region::all_convexes());


}  /* end of namespace getfem.                                             */


#endif /* GETFEM_GENERIC_ASSEMBLY_H__  */
