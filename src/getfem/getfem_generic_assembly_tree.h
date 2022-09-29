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

/** @file   getfem_generic_assembly_tree.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @date   November 18, 2013.
    @brief  Definition of the syntax tree and basic operations on it.
            Internal header for the generic assembly language part. 
 */


#ifndef GETFEM_GENERIC_ASSEMBLY_TREE_H__
#define GETFEM_GENERIC_ASSEMBLY_TREE_H__

#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_models.h"
#include "gmm/gmm_blas.h"
#include <iomanip>
#include "getfem/getfem_omp.h"
#include "getfem/dal_singleton.h"
#include "getfem/bgeot_rtree.h"
#include "getfem/bgeot_geotrans_inv.h"
#include "getfem/getfem_copyable_ptr.h"
#ifndef _WIN32
extern "C"{
#include <unistd.h>
}
#endif

#define GA_DEBUG_ASSERT(a, b) GMM_ASSERT1(a, b)
// #define GA_DEBUG_ASSERT(a, b)

#if 1
#  define GA_TIC
#  define GA_TOC(a)
#  define GA_TOCTIC(a)
#else
#  define GA_TIC scalar_type _ga_time_ = gmm::uclock_sec();
#  define GA_TOC(a) { cout <<(a)<<" : "<<gmm::uclock_sec()-_ga_time_<< endl; }
#  define GA_TOCTIC(a) { GA_TOC(a); _ga_time_ = gmm::uclock_sec(); }
#endif

namespace getfem {

  // Basic token types (basic language components)
  enum GA_TOKEN_TYPE {
    GA_INVALID = 0, // invalid token
    GA_END,         // string end
    GA_NAME,        // A variable or user defined nonlinear function name
    GA_SCALAR,      // A real number
    GA_PLUS,        // '+'
    GA_MINUS,       // '-'
    GA_UNARY_MINUS, // '-'
    GA_MULT,        // '*'
    GA_DIV,         // '/'
    GA_COLON,       // ':'
    GA_QUOTE,       // ''' transpose
    GA_COLON_EQ,    // ':=' macro def
    GA_DEF,         // 'Def' macro def
    GA_SYM,         // 'Sym(M)' operator
    GA_SKEW,        // 'Skew(M)' operator
    GA_TRACE,       // 'Trace(M)' operator
    GA_DEVIATOR,    // 'Deviator' operator
    GA_INTERPOLATE, // 'Interpolate' operation
    GA_INTERPOLATE_FILTER, // 'Interpolate_filter' operation
    GA_INTERPOLATE_DERIVATIVE, // 'Interpolate_derivative' operation
    GA_ELEMENTARY,  // 'Elementary' operation (operation at the element level)
    GA_SECONDARY_DOMAIN,  // For the integration on a product of two domains
    GA_XFEM_PLUS,   // Evaluation on the + side of a level-set for fem_level_set
    GA_XFEM_MINUS,  // Evaluation on the - side of a level-set for fem_level_set
    GA_PRINT,       // 'Print' Print the tensor
    GA_DOT,         // '.'
    GA_DOTMULT,     // '.*' componentwise multiplication
    GA_DOTDIV,      // './' componentwise division
    GA_TMULT,       // '@' tensor product
    GA_COMMA,       // ','
    GA_DCOMMA,      // ',,'
    GA_SEMICOLON,   // ';'
    GA_DSEMICOLON,  // ';;'
    GA_LPAR,        // '('
    GA_RPAR,        // ')'
    GA_LBRACKET,    // '['
    GA_RBRACKET,    // ']'
    GA_NB_TOKEN_TYPE
  };

  // Detects Grad_, Hess_ or Div_
  size_type ga_parse_prefix_operator(std::string &name);
  // Detects Test_ and Test2_
  size_type ga_parse_prefix_test(std::string &name);

  // Types of nodes for the syntax tree
  enum GA_NODE_TYPE {
    GA_NODE_VOID = 0,
    GA_NODE_OP,
    GA_NODE_PREDEF_FUNC,
    GA_NODE_SPEC_FUNC,
    GA_NODE_OPERATOR,
    GA_NODE_CONSTANT,
    GA_NODE_NAME,
    GA_NODE_MACRO_PARAM,
    GA_NODE_PARAMS,
    GA_NODE_RESHAPE,
    GA_NODE_CROSS_PRODUCT,
    GA_NODE_SWAP_IND,
    GA_NODE_IND_MOVE_LAST,
    GA_NODE_CONTRACT,
    GA_NODE_ALLINDICES,
    GA_NODE_C_MATRIX,
    GA_NODE_X,
    GA_NODE_ELT_SIZE,
    GA_NODE_ELT_K,
    GA_NODE_ELT_B,
    GA_NODE_NORMAL,
    GA_NODE_VAL,
    GA_NODE_GRAD,
    GA_NODE_HESS,
    GA_NODE_DIVERG,
    GA_NODE_VAL_TEST,
    GA_NODE_GRAD_TEST,
    GA_NODE_HESS_TEST,
    GA_NODE_DIVERG_TEST,
    GA_NODE_INTERPOLATE,
    GA_NODE_INTERPOLATE_FILTER,
    GA_NODE_INTERPOLATE_VAL,
    GA_NODE_INTERPOLATE_GRAD,
    GA_NODE_INTERPOLATE_HESS,
    GA_NODE_INTERPOLATE_DIVERG,
    GA_NODE_INTERPOLATE_VAL_TEST,
    GA_NODE_INTERPOLATE_GRAD_TEST,
    GA_NODE_INTERPOLATE_HESS_TEST,
    GA_NODE_INTERPOLATE_DIVERG_TEST,
    GA_NODE_INTERPOLATE_X,
    GA_NODE_INTERPOLATE_ELT_K,
    GA_NODE_INTERPOLATE_ELT_B,
    GA_NODE_INTERPOLATE_NORMAL,
    GA_NODE_INTERPOLATE_DERIVATIVE,
    GA_NODE_ELEMENTARY,
    GA_NODE_ELEMENTARY_VAL,
    GA_NODE_ELEMENTARY_GRAD,
    GA_NODE_ELEMENTARY_HESS,
    GA_NODE_ELEMENTARY_DIVERG,
    GA_NODE_ELEMENTARY_VAL_TEST,
    GA_NODE_ELEMENTARY_GRAD_TEST,
    GA_NODE_ELEMENTARY_HESS_TEST,
    GA_NODE_ELEMENTARY_DIVERG_TEST,
    GA_NODE_SECONDARY_DOMAIN,
    GA_NODE_SECONDARY_DOMAIN_VAL,
    GA_NODE_SECONDARY_DOMAIN_GRAD,
    GA_NODE_SECONDARY_DOMAIN_HESS,
    GA_NODE_SECONDARY_DOMAIN_DIVERG,
    GA_NODE_SECONDARY_DOMAIN_VAL_TEST,
    GA_NODE_SECONDARY_DOMAIN_GRAD_TEST,
    GA_NODE_SECONDARY_DOMAIN_HESS_TEST,
    GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST,
    GA_NODE_SECONDARY_DOMAIN_X,
    GA_NODE_SECONDARY_DOMAIN_NORMAL,
    GA_NODE_XFEM_PLUS,
    GA_NODE_XFEM_PLUS_VAL,
    GA_NODE_XFEM_PLUS_GRAD,
    GA_NODE_XFEM_PLUS_HESS,
    GA_NODE_XFEM_PLUS_DIVERG,
    GA_NODE_XFEM_PLUS_VAL_TEST,
    GA_NODE_XFEM_PLUS_GRAD_TEST,
    GA_NODE_XFEM_PLUS_HESS_TEST,
    GA_NODE_XFEM_PLUS_DIVERG_TEST,
    GA_NODE_XFEM_MINUS,
    GA_NODE_XFEM_MINUS_VAL,
    GA_NODE_XFEM_MINUS_GRAD,
    GA_NODE_XFEM_MINUS_HESS,
    GA_NODE_XFEM_MINUS_DIVERG,
    GA_NODE_XFEM_MINUS_VAL_TEST,
    GA_NODE_XFEM_MINUS_GRAD_TEST,
    GA_NODE_XFEM_MINUS_HESS_TEST,
    GA_NODE_XFEM_MINUS_DIVERG_TEST,
    GA_NODE_ZERO};

  typedef std::shared_ptr<std::string> pstring;
  // Print error message indicating the position in the assembly string
  void ga_throw_error_msg(pstring expr, size_type pos,
                          const std::string &msg);

# define ga_throw_error(expr, pos, msg)               \
  { std::stringstream ss; ss << msg;                  \
    ga_throw_error_msg(expr, pos, ss.str());          \
    GMM_ASSERT1(false, "Error in assembly string" );  \
  }
  
  // Structure for the tensor associated with a tree node
  struct assembly_tensor {
    bool is_copied;
    int sparsity_; // 0: plain, 1: vectorized base, 2: vectorised grad, ...
    size_type qdim_; // Dimension of the vectorization for sparsity tensors
    base_tensor t;
    assembly_tensor *tensor_copied;

    const base_tensor &org_tensor() const
    { return is_copied ? tensor_copied->org_tensor() : t; }
    base_tensor &org_tensor()
    { return is_copied ? tensor_copied->org_tensor() : t; }

    const base_tensor &tensor() const
    { return (is_copied ? tensor_copied->tensor() : t); }

    base_tensor &tensor()
    { return (is_copied ? tensor_copied->tensor() : t); }

    void set_sparsity(int sp, size_type q)
    { sparsity_ = sp; qdim_ = q; }

    size_type qdim() const { return is_copied ? tensor_copied->qdim() : qdim_; }

    int sparsity() const
    { return is_copied ? tensor_copied->sparsity() : sparsity_; }

    inline void set_to_original() { is_copied = false; }
    inline void set_to_copy(assembly_tensor &t_) {
      is_copied = true; sparsity_ = t_.sparsity_; qdim_ = t_.qdim_;
      t = t_.org_tensor(); tensor_copied = &(t_);
    }

    inline void adjust_sizes(const bgeot::multi_index &ssizes)
    { t.adjust_sizes(ssizes); }

    inline void adjust_sizes()
    { if (t.sizes().size() || t.size() != 1) t.init(); }

    inline void adjust_sizes(size_type i)
    { if (t.sizes().size() != 1 || t.sizes()[0] != i) t.init(i); }

    inline void adjust_sizes(size_type i, size_type j) {
      if (t.sizes().size() != 2 || t.sizes()[0] != i || t.sizes()[1] != j)
        t.init(i, j);
    }

    inline void adjust_sizes(size_type i, size_type j, size_type k) {
      if (t.sizes().size() != 3 || t.sizes()[0] != i || t.sizes()[1] != j
          || t.sizes()[2] != k)
        t.init(i, j, k);
    }
    inline void adjust_sizes(size_type i, size_type j,
                             size_type k, size_type l) {
      if (t.sizes().size() != 3 || t.sizes()[0] != i || t.sizes()[1] != j
          || t.sizes()[2] != k || t.sizes()[3] != l)
       t.init(i, j, k, l);
    }

    void init_scalar_tensor(scalar_type v)
    { set_to_original(); t.adjust_sizes(); t[0] = v; }

    void init_vector_tensor(size_type d)
    { set_to_original(); t.adjust_sizes(d); }

    void init_matrix_tensor(size_type n, size_type m)
    { set_to_original(); t.adjust_sizes(n, m); }

    void init_identity_matrix_tensor(size_type n) {
      init_matrix_tensor(n, n);
      auto itw = t.begin();
      for (size_type i = 0; i < n; ++i)
        for (size_type j = 0; j < n; ++j)
          *itw++ = (i == j) ? scalar_type(1) : scalar_type(0);
    }

    void init_third_order_tensor(size_type n, size_type m,  size_type l)
    { set_to_original(); t.adjust_sizes(n, m, l); }

    void init_fourth_order_tensor(size_type n, size_type m,
                                  size_type l, size_type k)
    { set_to_original(); t.adjust_sizes(n, m, l, k); }

    const bgeot::multi_index &sizes() const { return t.sizes(); }

    assembly_tensor()
    : is_copied(false), sparsity_(0), qdim_(0), tensor_copied(0) {}
  };

  struct ga_tree_node;
  typedef ga_tree_node *pga_tree_node;

  
  struct ga_tree_node {
    GA_NODE_TYPE node_type;
    GA_TOKEN_TYPE op_type;
    assembly_tensor t;
    size_type test_function_type; // -1 = undetermined
                                  // 0 = no test function,
                                  // 1 = first order, 2 = second order,
                                  // 3 = both with always first order in first
    std::string name_test1, name_test2; // variable names corresponding to test
                                  // functions when test_function_type > 0.
    std::string interpolate_name_test1, interpolate_name_test2; // name
                                  // of interpolation transformation if any
    size_type qdim1, qdim2;       // Qdims when test_function_type > 0.
    size_type nbc1, nbc2, nbc3;   // For X (nbc1=coordinate number),
                                  // macros (nbc1=param number, nbc2,nbc3 type))
                                  // and C_MATRIX (nbc1=order).
    size_type pos;                // Position of the first character in string
    pstring expr;                 // Original string, for error messages.
    std::string name;             // variable/constant/function/operator name
    std::string interpolate_name; // For Interpolate : name of transformation
    std::string interpolate_name_der; // For Interpolate derivative:
                                      // name of transformation
    std::string elementary_name;  // For Elementary_transformation :
                                  // name of transformation
    std::string elementary_target;// For Elementary_transformation :
                                  // target variable (for its mesh_fem) 
    size_type der1, der2;         // For functions and nonlinear operators,
                                  // optional derivative or second derivative.
    bool symmetric_op;
    pga_tree_node parent;         // Parent node
    std::vector<pga_tree_node> children; // Children nodes
    scalar_type hash_value;       // Hash value to identify nodes.
    bool marked;                  // For specific use of some algorithms

    inline const base_tensor &tensor() const { return t.tensor(); }
    inline base_tensor &tensor() { return t.tensor(); }
    int sparsity() const { return t.sparsity(); }

    inline size_type nb_test_functions() const {
      if (test_function_type == size_type(-1)) return 0;
      return test_function_type - (test_function_type >= 2 ? 1 : 0);
    }

    inline size_type tensor_order() const
    { return t.sizes().size() - nb_test_functions(); }

    inline size_type tensor_test_size() const {
      size_type st = nb_test_functions();
      return (st >= 1 ? t.sizes()[0] : 1) * (st == 2 ? t.sizes()[1] : 1);
    }

    inline size_type tensor_proper_size() const
    { return t.org_tensor().size() / tensor_test_size(); }

    inline size_type tensor_proper_size(size_type i) const
    { return t.sizes()[nb_test_functions()+i]; }


    void mult_test(const pga_tree_node n0, const pga_tree_node n1);

    bool tensor_is_zero() {
      if (node_type == GA_NODE_ZERO) return true;
      if (node_type != GA_NODE_CONSTANT) return false;
      for (size_type i = 0; i < tensor().size(); ++i)
        if (tensor()[i] != scalar_type(0)) return false;
      return true;
    }

    inline bool is_constant() {
      return (node_type == GA_NODE_CONSTANT ||
              (node_type == GA_NODE_ZERO && test_function_type == 0));
    }
    
    inline void init_scalar_tensor(scalar_type v)
    { t.init_scalar_tensor(v); test_function_type = 0; }

    inline void init_vector_tensor(size_type d)
    { t.init_vector_tensor(d); test_function_type = 0; }

    inline void init_matrix_tensor(size_type n, size_type m)
    { t.init_matrix_tensor(n, m); test_function_type = 0; }

    inline void init_identity_matrix_tensor(size_type n)
    { t.init_identity_matrix_tensor(n); test_function_type = 0; }

    inline void init_third_order_tensor(size_type n, size_type m,  size_type l)
    { t.init_third_order_tensor(n, m, l); test_function_type = 0; }

    inline void init_fourth_order_tensor(size_type n, size_type m,
                                         size_type l, size_type k)
    { t.init_fourth_order_tensor(n, m, l, k); test_function_type = 0; }

    inline void adopt_child(pga_tree_node new_child)
    { children.push_back(new_child); children.back()->parent = this; }

    inline void replace_child(pga_tree_node oldchild,
                              pga_tree_node newchild) {
        bool found = false;
        for (pga_tree_node &child : children)
          if (child == oldchild) { child = newchild; found = true; }
        GMM_ASSERT1(found, "Internal error");
    }

    ga_tree_node()
      : node_type(GA_NODE_VOID), test_function_type(-1), qdim1(0), qdim2(0),
      nbc1(0), nbc2(0), nbc3(0), pos(0), expr(0), der1(0), der2(0),
        symmetric_op(false), hash_value(0) {}
    ga_tree_node(GA_NODE_TYPE ty, size_type p, pstring expr_)
      : node_type(ty), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
      pos(p), expr(expr_), der1(0), der2(0), symmetric_op(false),
      hash_value(0) {}
    ga_tree_node(scalar_type v, size_type p, pstring expr_)
      : node_type(GA_NODE_CONSTANT), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(p), expr(expr_), der1(0), der2(0), symmetric_op(false),
        hash_value(0)
    { init_scalar_tensor(v); }
    ga_tree_node(const char *n, size_type l, size_type p, pstring expr_)
      : node_type(GA_NODE_NAME), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(p), expr(expr_), name(n, l), der1(0), der2(0), symmetric_op(false),
        hash_value(0) {}
    ga_tree_node(GA_TOKEN_TYPE op, size_type p, pstring expr_)
      : node_type(GA_NODE_OP), op_type(op), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(p), expr(expr_), der1(0), der2(0), symmetric_op(false),
        hash_value(0) {}
  };

  struct ga_tree {
    pga_tree_node root, current_node;
    std::string secondary_domain;

    void add_scalar(scalar_type val, size_type pos, pstring expr);
    void add_allindices(size_type pos, pstring expr);
    void add_name(const char *name, size_type length, size_type pos,
                  pstring expr);
    void add_sub_tree(ga_tree &sub_tree);
    void add_params(size_type pos, pstring expr);
    void add_matrix(size_type pos, pstring expr);
    void add_op(GA_TOKEN_TYPE op_type, size_type pos, pstring expr);
    void clear_node_rec(pga_tree_node pnode);
    void clear_node(pga_tree_node pnode);
    void clear() { clear_node_rec(root); root = current_node = nullptr; }
    void clear_children(pga_tree_node pnode);
    void replace_node_by_child(pga_tree_node pnode, size_type i);
    void copy_node(pga_tree_node pnode, pga_tree_node parent,
                   pga_tree_node &child);
    void duplicate_with_operation(pga_tree_node pnode, GA_TOKEN_TYPE op_type);
    void duplicate_with_addition(pga_tree_node pnode)
    { duplicate_with_operation(pnode, GA_PLUS); }
    void duplicate_with_substraction(pga_tree_node pnode)
    { duplicate_with_operation(pnode, GA_MINUS); }
    void insert_node(pga_tree_node pnode, GA_NODE_TYPE node_type);
    void add_child(pga_tree_node pnode, GA_NODE_TYPE node_type = GA_NODE_VOID);
    void swap(ga_tree &tree) {
      std::swap(root, tree.root);
      std::swap(current_node, tree.current_node);
      std::swap(secondary_domain, tree.secondary_domain);
    }

    ga_tree() : root(nullptr), current_node(nullptr), secondary_domain() {}

    ga_tree(const ga_tree &tree) : root(nullptr), current_node(nullptr),
      secondary_domain(tree.secondary_domain)
    { if (tree.root) copy_node(tree.root, nullptr, root); }

    ga_tree &operator =(const ga_tree &tree) {
      clear(); secondary_domain = tree.secondary_domain;
      if (tree.root)
        copy_node(tree.root,nullptr,root);
      return *this;
    }

    ~ga_tree() { clear(); }
  };

  // Test equality or equivalence of two sub trees.
  // version = 0 : strict equality
  //           1 : give the same result
  //           2 : give the same result with transposition of test functions
  bool sub_tree_are_equal
  (const pga_tree_node pnode1, const pga_tree_node pnode2,
   const ga_workspace &workspace, int version);

  // Transform the expression of a node and its sub-nodes in the equivalent
  // assembly string sent to ostream str
  void ga_print_node(const pga_tree_node pnode,
                     std::ostream &str);
  // The same for the whole tree, the result is a std::string
  std::string ga_tree_to_string(const ga_tree &tree);

  // Syntax analysis of an assembly string. Conversion to a tree.
  // No semantic analysis is done. The tree can be inconsistent.
  void ga_read_string(const std::string &expr, ga_tree &tree,
                      const ga_macro_dictionary &macro_dict);
  void ga_read_string_reg(const std::string &expr, ga_tree &tree,
                          ga_macro_dictionary &macro_dict);


} /* end of namespace */


#endif /* GETFEM_GENERIC_ASSEMBLY_TREE_H__  */
