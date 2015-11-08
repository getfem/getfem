/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2013-2015 Yves Renard

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


#include "getfem/getfem_generic_assembly.h"
#include "gmm/gmm_blas.h"
#include <iomanip>
#include "getfem/getfem_omp.h"
#include "getfem/dal_singleton.h"
#include "getfem/bgeot_rtree.h"
#include "getfem/bgeot_geotrans_inv.h"
#include "getfem/getfem_copyable_ptr.h"

/**
   Providing for special Math functions unavailable on Intel or MSVS C++
   compilers
*/

#if defined(_MSC_VER) && _MSC_VER < 1800
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/erf.hpp>
typedef double (*BoostMathFunction)(double);
BoostMathFunction const acosh = boost::math::acosh<double>;
BoostMathFunction const asinh = boost::math::asinh<double>;
BoostMathFunction const atanh = boost::math::atanh<double>;
BoostMathFunction const erf = boost::math::erf<double>;
BoostMathFunction const erfc = boost::math::erfc<double>;
#endif


// #define GA_USES_BLAS // not so interesting, at leat for debian blas

// #define GA_DEBUG_INFO(a) { cout << a << endl; }
#define GA_DEBUG_INFO(a)
#define GA_DEBUG_ASSERT(a, b) GMM_ASSERT1(a, b)
// #define GA_DEBUG_ASSERT(a, b)

#if 1
#  define GA_TIC
#  define GA_TOC(a)
#  define GA_TOCTIC(a)
#else
#  define GA_TIC scalar_type _ga_time_ = gmm::uclock_sec();
#  define GA_TOC(a) { cout << (a) << " : " << gmm::uclock_sec() - _ga_time_ << endl; }
#  define GA_TOCTIC(a) { GA_TOC(a); _ga_time_ = gmm::uclock_sec(); }
#endif

namespace getfem {

  //=========================================================================
  // Lexical analysis for the generic assembly language
  //=========================================================================

  // Basic token types
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
    GA_TRACE,       // 'Trace' operator
    GA_DEVIATOR,    // 'Deviator' operator
    GA_INTERPOLATE, // 'Interpolate' operation
    GA_INTERPOLATE_FILTER, // 'Interpolate_filter' operation
    GA_ELEMENTARY,  // 'Elementary' operation (operation at the element level)
    GA_XFEM_PLUS,   // Évaluation on the + side of a level-set for fem_level_set
    GA_XFEM_MINUS,  // Évaluation on the - side of a level-set for fem_level_set
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

  static GA_TOKEN_TYPE ga_char_type[256];
  static int ga_operator_priorities[GA_NB_TOKEN_TYPE];

  // Initialize ga_char_type and ga_operator_priorities arrays
  static bool init_ga_char_type(void) {
    for (int i = 0; i < 256; ++i) ga_char_type[i] = GA_INVALID;
    ga_char_type['+'] = GA_PLUS;        ga_char_type['-'] = GA_MINUS;
    ga_char_type['*'] = GA_MULT;        ga_char_type['/'] = GA_DIV;
    ga_char_type[':'] = GA_COLON;       ga_char_type['\''] = GA_QUOTE;
    ga_char_type['.'] = GA_DOT;         ga_char_type['@'] = GA_TMULT;
    ga_char_type[','] = GA_COMMA;       ga_char_type[';'] = GA_SEMICOLON;
    ga_char_type['('] = GA_LPAR;        ga_char_type[')'] = GA_RPAR;
    ga_char_type['['] = GA_LBRACKET;    ga_char_type[']'] = GA_RBRACKET;
    ga_char_type['_'] = GA_NAME;
    for (unsigned i = 'a'; i <= 'z'; ++i)  ga_char_type[i] = GA_NAME;
    for (unsigned i = 'A'; i <= 'Z'; ++i)  ga_char_type[i] = GA_NAME;
    for (unsigned i = '0'; i <= '9'; ++i)  ga_char_type[i] = GA_SCALAR;

    for (unsigned i = 0; i < GA_NB_TOKEN_TYPE; ++i)
      ga_operator_priorities[i] = 0;
    ga_operator_priorities[GA_PLUS] = 1;
    ga_operator_priorities[GA_MINUS] = 1;
    ga_operator_priorities[GA_MULT] = 2;
    ga_operator_priorities[GA_DIV] = 2;
    ga_operator_priorities[GA_COLON] = 2;
    ga_operator_priorities[GA_DOT] = 2;
    ga_operator_priorities[GA_DOTMULT] = 2;
    ga_operator_priorities[GA_DOTDIV] = 2;
    ga_operator_priorities[GA_TMULT] = 2;
    ga_operator_priorities[GA_QUOTE] = 3;
    ga_operator_priorities[GA_UNARY_MINUS] = 3;
    ga_operator_priorities[GA_TRACE] = 4;
    ga_operator_priorities[GA_DEVIATOR] = 4;
    ga_operator_priorities[GA_PRINT] = 4;

    return true;
  }

  static bool ga_initialized = init_ga_char_type();

  // Get the next token in the string at position 'pos' end return its type
  static GA_TOKEN_TYPE ga_get_token(const std::string &expr,
                                    size_type &pos,
                                    size_type &token_pos,
                                    size_type &token_length) {
    bool fdot = false, fE = false;
    GMM_ASSERT1(ga_initialized, "Internal error");

    // Ignore white spaces
    while (expr[pos] == ' ' && pos < expr.size()) ++pos;
    token_pos = pos;
    token_length = 0;

    // Detecting end of expression
    if (pos >= expr.size()) return GA_END;

    // Treating the different cases (Operation, name or number)
    GA_TOKEN_TYPE type = ga_char_type[unsigned(expr[pos])];
    ++pos; ++token_length;
    switch (type) {
    case GA_DOT:
      if (pos >= expr.size()) return type;
      if (expr[pos] == '*') { ++pos; ++token_length; return GA_DOTMULT; }
      if (expr[pos] == '/') { ++pos; ++token_length; return GA_DOTDIV; }
      if (ga_char_type[unsigned(expr[pos])] != GA_SCALAR)
        return type;
      fdot = true; type = GA_SCALAR;
    case GA_SCALAR:
      while (pos < expr.size()) {
        GA_TOKEN_TYPE ctype = ga_char_type[unsigned(expr[pos])];
        switch (ctype) {
        case GA_DOT:
          if (fdot) return type;
          fdot = true; ++pos; ++token_length;
          break;
        case GA_NAME:
          if (fE || (expr[pos] != 'E' && expr[pos] != 'e')) return type;
          fE = true; fdot = true; ++pos; ++token_length;
          if (pos < expr.size()) {
            if (expr[pos] == '+' || expr[pos] == '-')
              { ++pos; ++token_length; }
          }
          if (pos >= expr.size()
              || ga_char_type[unsigned(expr[pos])] != GA_SCALAR)
            return GA_INVALID;
          break;
        case GA_SCALAR:
          ++pos; ++token_length; break;
        default:
          return type;
        }
      }
      return type;
    case GA_NAME:
      while (pos < expr.size()) {
        GA_TOKEN_TYPE ctype = ga_char_type[unsigned(expr[pos])];
        if (ctype != GA_SCALAR && ctype != GA_NAME) break;
        ++pos; ++token_length;
      }
      if (expr.compare(token_pos, token_length, "Trace") == 0)
        return GA_TRACE;
      if (expr.compare(token_pos, token_length, "Deviator") == 0)
        return GA_DEVIATOR;
      if (expr.compare(token_pos, token_length, "Interpolate") == 0)
        return GA_INTERPOLATE;
      if (expr.compare(token_pos, token_length, "Interpolate_filter") == 0)
        return GA_INTERPOLATE_FILTER;
      if (expr.compare(token_pos, token_length,
                       "Elementary_transformation") == 0)
        return GA_ELEMENTARY;
      if (expr.compare(token_pos, token_length, "Xfem_plus") == 0)
        return GA_XFEM_PLUS;
      if (expr.compare(token_pos, token_length, "Xfem_minus") == 0)
        return GA_XFEM_MINUS;
      if (expr.compare(token_pos, token_length, "Print") == 0)
        return GA_PRINT;
      return type;
    case GA_COMMA:
      if (pos < expr.size() &&
          ga_char_type[unsigned(expr[pos])] == GA_COMMA) {
        ++pos; return GA_DCOMMA;
      }
      return type;
    case GA_SEMICOLON:
      if (pos < expr.size() &&
          ga_char_type[unsigned(expr[pos])] == GA_SEMICOLON) {
        ++pos; return GA_DSEMICOLON;
      }
      return type;
    default: return type;
    }
  }


  //=========================================================================
  // Tree structure for syntax analysis
  //=========================================================================

  static void ga_throw_error_msg(const std::string &expr, size_type pos,
                                 const std::string &msg) {
    int length_before = 70, length_after = 70;
    if (expr.size()) {
      int first = std::max(0, int(pos)-length_before);
      int last = std::min(int(pos)+length_after, int(expr.size()));
      if (last - first < length_before+length_after)
      first = std::max(0, int(pos)-length_before
                       -(length_before+length_after-last+first));
      if (last - first < length_before+length_after)
        last = std::min(int(pos)+length_after
                        +(length_before+length_after-last+first),
                        int(expr.size()));
      if (first > 0) cerr << "...";
      cerr << expr.substr(first, last-first);
      if (last < int(expr.size())) cerr << "...";
      cerr << endl;
      if (first > 0) cerr << "   ";
      if (int(pos) > first)
        cerr << std::setfill ('-') << std::setw(int(pos)-first) << '-'
             << std::setfill (' ');
      cerr << "^" << endl;
    }
    cerr << msg << endl;
  }

#define ga_throw_error(expr, pos, msg)               \
  { std::stringstream ss; ss << msg;                 \
    ga_throw_error_msg(expr, pos, ss.str());         \
    GMM_ASSERT1(false, "Error in assembly string" ); \
  }

  enum GA_NODE_TYPE {
    GA_NODE_VOID = 0,
    GA_NODE_OP,
    GA_NODE_PREDEF_FUNC,
    GA_NODE_SPEC_FUNC,
    GA_NODE_OPERATOR,
    GA_NODE_CONSTANT,
    GA_NODE_NAME,
    GA_NODE_PARAMS,
    GA_NODE_RESHAPE,
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
    GA_NODE_INTERPOLATE_NORMAL,
    GA_NODE_INTERPOLATE_X,
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

  struct ga_tree_node;
  typedef ga_tree_node *pga_tree_node;

  struct ga_tree_node {
    GA_NODE_TYPE node_type;
    base_tensor t;
    size_type test_function_type; // -1 = undetermined
                                  // 0 = no test function,
                                  // 1 = first order, 2 = second order,
                                  // 3 = both with always first order in first
    std::string name_test1, name_test2; // variable names corresponding to test
                                  // functions when test_function_type > 0.
    std::string interpolate_name_test1, interpolate_name_test2; // name
                                  // of interpolation transformation if any
    size_type qdim1, qdim2;       // Qdims when test_function_type > 0.
    size_type nbc1, nbc2, nbc3;   // For explicit matrices and x.
    size_type pos;                // Position of the first character in string
    std::string name;             // variable/constant/function/operator name
    std::string interpolate_name; // For Interpolate : name of transformation
    std::string interpolate_name_der; // For Interpolate derivative:
                                      // name of transformation
    std::string elementary_name;  // For Elementary_transformation :
                                  // name of transformation
    size_type der1, der2;         // For functions and nonlinear operators,
                                  // optional derivative or second derivative.
    GA_TOKEN_TYPE op_type;
    bool symmetric_op;
    pga_tree_node parent;         // Parent node
    std::vector<pga_tree_node> children; // Children nodes
    scalar_type hash_value;       // Hash value to identify nodes.
    bool marked;                  // For specific use of some algorithms

    inline size_type nb_test_functions(void) const {
      if (test_function_type == size_type(-1)) return 0;
      return test_function_type  - (test_function_type >= 2 ? 1 : 0);
    }

    inline size_type tensor_order(void) const
    { return t.sizes().size() - nb_test_functions(); }

    inline size_type tensor_test_size(void) const {
      size_type st = nb_test_functions();
      return (st >= 1 ? t.sizes()[0] : 1) * (st == 2 ? t.sizes()[1] : 1);
    }

    inline size_type tensor_proper_size(void) const
    { return t.size() / tensor_test_size(); }

    inline size_type tensor_proper_size(size_type i) const
    { return t.sizes()[nb_test_functions()+i]; }


    void mult_test(pga_tree_node n0, pga_tree_node n1,
                   const std::string &expr) {
      size_type test0 = n0->test_function_type, test1 = n1->test_function_type;
      if (test0 && test1 && (test0 == test1 ||
                             test0 >= 3 || test1 >= 3))
        ga_throw_error(expr, pos, "Incompatibility of test functions "
                       "in product.");
      GMM_ASSERT1(test0 != size_type(-1) && test1 != size_type(-1),
                  "internal error");

      test_function_type = test0 + test1;

      size_type st = nb_test_functions();
      bgeot::multi_index mi(st);

      switch (test0) {
      case 1: mi[0] = n0->t.sizes()[0]; break;
      case 2: mi[st-1] = n0->t.sizes()[0]; break;
      case 3: mi[0] = n0->t.sizes()[0]; mi[1] = n0->t.sizes()[1]; break;
      }
      switch (test1) {
      case 1: mi[0] = n1->t.sizes()[0]; break;
      case 2: mi[st-1] = n1->t.sizes()[0]; break;
      case 3: mi[0] = n1->t.sizes()[0]; mi[1] = n1->t.sizes()[1]; break;
      }

      if (n0->name_test1.size()) {
        name_test1 = n0->name_test1; qdim1 = n0->qdim1;
        interpolate_name_test1 = n0->interpolate_name_test1;
      } else {
        name_test1 = n1->name_test1; qdim1 = n1->qdim1;
        interpolate_name_test1 = n1->interpolate_name_test1;
      }

      if (n0->name_test2.size()) {
        name_test2 = n0->name_test2; qdim2 = n0->qdim2;
        interpolate_name_test2 = n0->interpolate_name_test2;
      } else {
        name_test2 = n1->name_test2; qdim2 = n1->qdim2;
        interpolate_name_test2 = n1->interpolate_name_test2;
      }
      t.adjust_sizes(mi);
    }

    bool tensor_is_zero(void) {
      if (node_type == GA_NODE_ZERO) return true;
      if (node_type != GA_NODE_CONSTANT) return false;
      for (size_type i = 0; i < t.size(); ++i)
        if (t[i] != scalar_type(0)) return false;
      return true;
    }

    void init_scalar_tensor(scalar_type v)
    { t.adjust_sizes(); t[0] = v; test_function_type = 0; }

    void init_vector_tensor(size_type d)
    { t.adjust_sizes(d); test_function_type = 0; }

    void init_matrix_tensor(size_type n, size_type m)
    { t.adjust_sizes(n, m); test_function_type = 0; }

    void init_third_order_tensor(size_type n, size_type m,  size_type l) {
      t.adjust_sizes(bgeot::multi_index(n,m,l));
      test_function_type = 0;
    }
    void init_fourth_order_tensor(size_type n, size_type m,
                                  size_type l, size_type k) {
      t.adjust_sizes(bgeot::multi_index(n,m,l,k));
      test_function_type = 0;
    }

    ga_tree_node(void)
      : node_type(GA_NODE_VOID), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(0), der1(0), der2(0), symmetric_op(false), hash_value(0) {}
    ga_tree_node(GA_NODE_TYPE ty, size_type p)
      : node_type(ty), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(p), der1(0), der2(0), symmetric_op(false), hash_value(0) {}
    ga_tree_node(scalar_type v, size_type p)
      : node_type(GA_NODE_CONSTANT), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(p), der1(0), der2(0), symmetric_op(false),
        hash_value(0)
    { init_scalar_tensor(v); }
    ga_tree_node(const char *n, size_type l, size_type p)
      : node_type(GA_NODE_NAME), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(p), name(n, l), der1(0), der2(0), symmetric_op(false),
        hash_value(0) {}
    ga_tree_node(GA_TOKEN_TYPE op, size_type p)
      : node_type(GA_NODE_OP), test_function_type(-1),
        qdim1(0), qdim2(0), nbc1(0), nbc2(0), nbc3(0),
        pos(p), der1(0), der2(0), op_type(op), symmetric_op(false),
        hash_value(0) {}

  };

  struct ga_tree {
    pga_tree_node root, current_node;

    void add_scalar(scalar_type val, size_type pos) {
      while (current_node && current_node->node_type != GA_NODE_OP)
        current_node = current_node->parent;
      if (current_node) {
        pga_tree_node new_node = new ga_tree_node(val, pos);
        current_node->children.push_back(new_node);
        new_node->parent = current_node;
        current_node = new_node;
      }
      else {
        GMM_ASSERT1(root == 0, "Invalid tree operation");
        current_node = root = new ga_tree_node(val, pos);
        root->parent = 0;
      }
    }

    void add_allindices(size_type pos) {
      while (current_node && current_node->node_type != GA_NODE_OP)
        current_node = current_node->parent;
      if (current_node) {
        pga_tree_node new_node = new ga_tree_node(GA_NODE_ALLINDICES, pos);
        current_node->children.push_back(new_node);
        new_node->parent = current_node;
        current_node = new_node;
      }
      else {
        GMM_ASSERT1(root == 0, "Invalid tree operation");
        current_node = root = new ga_tree_node(GA_NODE_ALLINDICES, pos);
        root->parent = 0;
      }
    }

    void add_name(const char *name, size_type length, size_type pos) {
      while (current_node && current_node->node_type != GA_NODE_OP)
        current_node = current_node->parent;
      if (current_node) {
        pga_tree_node new_node = new ga_tree_node(name, length, pos);
        current_node->children.push_back(new_node);
        new_node->parent = current_node;
        current_node = new_node;
      }
      else {
        GMM_ASSERT1(root == 0, "Invalid tree operation");
        current_node = root = new ga_tree_node(name, length, pos);
        root->parent = 0;
      }
    }

    void add_sub_tree(ga_tree &sub_tree) {
      if (current_node && (current_node->node_type == GA_NODE_PARAMS ||
                           current_node->node_type == GA_NODE_INTERPOLATE_FILTER ||
                           current_node->node_type == GA_NODE_C_MATRIX)) {
        GMM_ASSERT1(sub_tree.root, "Invalid tree operation");
        current_node->children.push_back(sub_tree.root);
        sub_tree.root->parent = current_node;
      } else {
        GMM_ASSERT1(sub_tree.root, "Invalid tree operation");
        while (current_node && current_node->node_type != GA_NODE_OP)
          current_node = current_node->parent;
        if (current_node) {
          current_node->children.push_back(sub_tree.root);
          sub_tree.root->parent = current_node;
          current_node = sub_tree.root;
        }
        else {
          GMM_ASSERT1(root == 0, "Invalid tree operation");
          current_node = root = sub_tree.root;
          root->parent = 0;
        }
      }
      sub_tree.root = sub_tree.current_node = 0;
    }

    void add_params(size_type pos) {
      GMM_ASSERT1(current_node, "internal error");
      while (current_node && current_node->parent &&
             current_node->parent->node_type == GA_NODE_OP &&
             ga_operator_priorities[current_node->parent->op_type] >= 4)
        current_node = current_node->parent;
      pga_tree_node new_node = new ga_tree_node(GA_NODE_PARAMS, pos);
      pga_tree_node parent =  current_node->parent;
      if (parent) {
        for (size_type i = 0; i < parent->children.size(); ++i)
          if (parent->children[i] == current_node)
            parent->children[i] = new_node;
      }
      else
        root = new_node;
      new_node->parent = current_node->parent;
      current_node->parent = new_node;
      new_node->children.push_back(current_node);
      current_node = new_node;
    }

    void add_matrix(size_type pos) {
      while (current_node && current_node->node_type != GA_NODE_OP)
        current_node = current_node->parent;
      if (current_node) {
        pga_tree_node new_node = new ga_tree_node(GA_NODE_C_MATRIX, pos);
        current_node->children.push_back(new_node);
        new_node->parent = current_node;
        current_node = new_node;
      }
      else {
        GMM_ASSERT1(root == 0, "Invalid tree operation");
        current_node = root = new ga_tree_node(GA_NODE_C_MATRIX, pos);
        root->parent = 0;
      }
      current_node->nbc1 = current_node->nbc2 = current_node->nbc3 = 0;
    }

    void zip_matrix(const pga_tree_node source_node) {
      GMM_ASSERT1(current_node->node_type == GA_NODE_C_MATRIX &&
                  source_node->node_type == GA_NODE_C_MATRIX,
                  "Internal error");
      size_type target_size = current_node->children.size();
      size_type source_size = source_node->children.size();
      size_type last_dim_size = target_size/source_size;
      GMM_ASSERT1(target_size == source_size*last_dim_size,
                  "Internal error, " << target_size << " != " <<
                  source_size << "*" << last_dim_size);
      std::vector<pga_tree_node> new_children;
      for (size_type i = 0; i < source_size; ++i) {
        for (size_type j = 0; j < last_dim_size; ++j)
          new_children.push_back(current_node->children[i*last_dim_size+j]);
        new_children.push_back(source_node->children[i]);
        source_node->children[i]->parent = current_node;
      }
      source_node->children.resize(0); // so that the destructor of source_node
                                       // will not destruct the children
      current_node->children = new_children;
    }

    void add_op(GA_TOKEN_TYPE op_type, size_type pos) {
      while (current_node && current_node->parent &&
             current_node->parent->node_type == GA_NODE_OP &&
             ga_operator_priorities[current_node->parent->op_type]
             >= ga_operator_priorities[op_type])
        current_node = current_node->parent;
      pga_tree_node new_node = new ga_tree_node(op_type, pos);
      if (current_node) {
        if (op_type == GA_UNARY_MINUS || op_type == GA_TRACE
            || op_type == GA_DEVIATOR || op_type == GA_PRINT) {
          current_node->children.push_back(new_node);
          new_node->parent = current_node;
        } else {
          pga_tree_node parent = current_node->parent;
          if (parent) {
            new_node->parent = parent;
            for (size_type i = 0; i < parent->children.size(); ++i)
              if (parent->children[i] == current_node)
                parent->children[i] = new_node;
          } else {
            root = new_node; new_node->parent = 0;
          }
          new_node->children.push_back(current_node);
          current_node->parent = new_node;
        }
      } else {
        if (root) new_node->children.push_back(root);
        root = new_node; new_node->parent = 0;
      }
      current_node = new_node;
    }

    void clear_node_rec(pga_tree_node pnode) {
      if (pnode) {
        for (size_type i = 0; i < pnode->children.size(); ++i)
          clear_node_rec(pnode->children[i]);
        delete pnode;
        current_node = 0;
      }
    }

    void clear_node(pga_tree_node pnode) {
      if (pnode) {
        pga_tree_node parent = pnode->parent;
        if (parent) {
          for (size_type i = 0, j = 0; i < parent->children.size(); ++i)
            if (parent->children[i] != pnode)
              { parent->children[j] = parent->children[i]; ++j; }
          parent->children.pop_back();
        } else root = 0;
      }
      clear_node_rec(pnode);
    }

    void clear() {
      if (root) clear_node_rec(root); root = current_node = 0;
    }

    void clear_children(pga_tree_node pnode) {
      for (size_type i = 0; i < pnode->children.size(); ++i)
        clear_node_rec(pnode->children[i]);
      pnode->children.resize(0);
    }

    void replace_node_by_child(pga_tree_node pnode, size_type i) {
      GMM_ASSERT1(i < pnode->children.size(), "Internal error");
      pga_tree_node child = pnode->children[i];
      if (pnode->parent) {
        bool found = false;
        for (size_type j = 0; j < pnode->parent->children.size(); ++j)
          if (pnode->parent->children[j] == pnode)
            { pnode->parent->children[j] = child; found = true; }
        GMM_ASSERT1(found, "Internal error");
      } else root = child;
      current_node = 0;
      child->parent = pnode->parent;
      for (size_type j = 0; j < pnode->children.size(); ++j)
        if (j != i) clear_node_rec(pnode->children[j]);
      delete pnode;
    }

    void copy_node(pga_tree_node pnode, pga_tree_node parent,
                   pga_tree_node &child) {
      child = new ga_tree_node();
      *child = *pnode;
      child->parent = parent;
      for (size_type j = 0; j < child->children.size(); ++j)
        child->children[j] = 0;
      for (size_type j = 0; j < child->children.size(); ++j)
        copy_node(pnode->children[j], child, child->children[j]);
    }

    void duplicate_with_operation(pga_tree_node pnode,
                                  GA_TOKEN_TYPE op_type) {
      pga_tree_node newop = new ga_tree_node(op_type, pnode->pos);
      newop->children.resize(2);
      newop->children[0] = pnode;
      newop->parent = pnode->parent;
      if (pnode->parent) {
        for (size_type j = 0; j < pnode->parent->children.size(); ++j)
          if (pnode->parent->children[j] == pnode)
            pnode->parent->children[j] = newop;
      } else root = newop;
      pnode->parent = newop;
      copy_node(pnode, newop, newop->children[1]);
    }

    void duplicate_with_addition(pga_tree_node pnode) {
      duplicate_with_operation(pnode, GA_PLUS);
    }

    void duplicate_with_subtraction(pga_tree_node pnode) {
      duplicate_with_operation(pnode, GA_MINUS);
    }

    void insert_node(pga_tree_node pnode, GA_NODE_TYPE node_type) {
      pga_tree_node newnode = new ga_tree_node();
      newnode->parent = pnode->parent;
      newnode->node_type = node_type;
      if (pnode->parent) {
        for (size_type j = 0; j < pnode->parent->children.size(); ++j)
          if (pnode->parent->children[j] == pnode)
            pnode->parent->children[j] = newnode;
      } else root = newnode;
      newnode->children.push_back(pnode);
      pnode->parent = newnode;
    }

    void add_child(pga_tree_node pnode) {
      pga_tree_node newnode = new ga_tree_node();
      newnode->parent = pnode;
      pnode->children.push_back(newnode);
    }

    void swap(ga_tree &tree)
    { std::swap(root, tree.root); std::swap(current_node, tree.current_node); }

    ga_tree(void) : root(0), current_node(0) {}

    ga_tree(const ga_tree &tree) : root(0), current_node(0)
    { if (tree.root) copy_node(tree.root, 0, root); }

    ga_tree &operator = (const ga_tree &tree)
    { clear(); if (tree.root) copy_node(tree.root, 0, root); return *this; }

    ~ga_tree() { clear(); }
  };


  // Test equality or equivalence of two sub trees.
  // version = 0 : strict equality
  //           1 : give the same result
  //           2 : give the same result with transposition of test functions
  static bool sub_tree_are_equal(pga_tree_node pnode1, pga_tree_node pnode2,
                                 const ga_workspace &workspace,
                                 int version) {
    size_type ntype1 = pnode1->node_type;
    if (ntype1 == GA_NODE_ZERO) ntype1 = GA_NODE_CONSTANT;
    size_type ntype2 = pnode2->node_type;
    if (ntype2 == GA_NODE_ZERO) ntype2 = GA_NODE_CONSTANT;

    if (ntype1 != ntype2) return false;
    if (pnode1->children.size() != pnode2->children.size()) return false;

    switch(ntype1) {
    case GA_NODE_OP:
      if (pnode1->op_type != pnode2->op_type) return false;
      if (pnode1->symmetric_op != pnode2->symmetric_op)  return false;
      break;
    case GA_NODE_OPERATOR:
      if (pnode1->der1 != pnode2->der1 || pnode1->der2 != pnode2->der2)
        return false;
    case GA_NODE_PREDEF_FUNC: case GA_NODE_SPEC_FUNC:
      if (pnode1->name.compare(pnode2->name)) return false;
      break;
    case GA_NODE_CONSTANT: case GA_NODE_ZERO:
      if (pnode1->t.size() != pnode2->t.size()) return false;

      switch(version) {
      case 0: case 1:
        if (pnode1->test_function_type != pnode2->test_function_type)
          return false;
        if ((pnode1->test_function_type & 1) &&
            pnode1->name_test1.compare(pnode2->name_test1) != 0)
          return false;
        if ((pnode1->test_function_type & 2) &&
            pnode1->name_test2.compare(pnode2->name_test2) != 0)
          return false;
        break;
      case 2:
        if ((pnode1->test_function_type == 1 &&
             pnode2->test_function_type == 1) ||
            (pnode1->test_function_type == 2 &&
             pnode2->test_function_type == 2))
          return false;
        if ((pnode1->test_function_type & 1) &&
            pnode1->name_test1.compare(pnode2->name_test2) != 0)
          return false;
        if ((pnode1->test_function_type & 2) &&
            pnode1->name_test2.compare(pnode2->name_test1) != 0)
          return false;
        break;
      }
      if (pnode1->t.size() != 1 &&
          pnode1->t.sizes().size() != pnode2->t.sizes().size()) return false;
      for (size_type i = 0; i < pnode1->t.sizes().size(); ++i)
        if (pnode1->t.sizes()[i] != pnode2->t.sizes()[i]) return false;
      for (size_type i = 0; i < pnode1->t.size(); ++i)
        if (gmm::abs(pnode1->t[i] - pnode2->t[i]) > 1E-25) return false;
      break;
    case GA_NODE_C_MATRIX:
      if (pnode1->nbc1 != pnode2->nbc1 || pnode1->nbc2 != pnode2->nbc2 ||
          pnode1->nbc3 != pnode2->nbc3)
        return false;
      break;
    case GA_NODE_INTERPOLATE_FILTER:
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name) ||
          pnode1->nbc1 != pnode2->nbc1)
        return false;
      break;
    case GA_NODE_INTERPOLATE_NORMAL: case GA_NODE_INTERPOLATE_X:
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name))
        return false;
      break;
    case GA_NODE_INTERPOLATE_DERIVATIVE:
      if (pnode1->interpolate_name_der.compare(pnode2->interpolate_name_der))
        return false;
      // The test continues with what follows
    case GA_NODE_INTERPOLATE_VAL_TEST: case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST: case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_VAL_TEST: case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST: case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST: case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST: case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name) ||
          pnode1->elementary_name.compare(pnode2->elementary_name))
        return false;
      // The test continues with what follows
    case GA_NODE_VAL_TEST: case GA_NODE_GRAD_TEST:
    case GA_NODE_HESS_TEST: case GA_NODE_DIVERG_TEST:
      {
        const mesh_fem *mf1 = workspace.associated_mf(pnode1->name);
        const mesh_fem *mf2 = workspace.associated_mf(pnode2->name);
        switch (version) {
        case 0:
          if (pnode1->name.compare(pnode2->name) ||
              pnode1->test_function_type != pnode2->test_function_type)
            return false;
          break;
        case 1:
          if (mf1 != mf2 ||
              workspace.qdim(pnode1->name) != workspace.qdim(pnode2->name) ||
              pnode1->test_function_type != pnode2->test_function_type)
            return false;
          break;
        case 2:
          if (mf1 != mf2 ||
              workspace.qdim(pnode1->name) != workspace.qdim(pnode2->name) ||
              pnode1->test_function_type == pnode2->test_function_type)
            return false;
          break;
        }
      }
      break;
    case GA_NODE_VAL: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
      if (pnode1->name.compare(pnode2->name)) return false;
      break;
    case GA_NODE_INTERPOLATE_VAL: case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS: case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name) ||
          pnode1->elementary_name.compare(pnode2->elementary_name) ||
          pnode1->name.compare(pnode2->name))
        return false;
      break;
    case GA_NODE_X:
      if (pnode1->nbc1 != pnode2->nbc1) return false;
      break;

    default:break;
    }

    if (version && ntype1 == GA_NODE_OP && pnode1->symmetric_op) {
      if (sub_tree_are_equal(pnode1->children[0], pnode2->children[0],
                             workspace, version) &&
          sub_tree_are_equal(pnode1->children[1], pnode2->children[1],
                             workspace, version))
        return true;
      if (sub_tree_are_equal(pnode1->children[1], pnode2->children[0],
                             workspace, version) &&
          sub_tree_are_equal(pnode1->children[0], pnode2->children[1],
                             workspace, version) )
        return true;
      return false;
    } else {
      for (size_type i = 0; i < pnode1->children.size(); ++i)
        if (!(sub_tree_are_equal(pnode1->children[i], pnode2->children[i],
                                 workspace, version)))
          return false;
    }
    return true;
  }

  static void verify_tree(pga_tree_node pnode, pga_tree_node parent) {
    GMM_ASSERT1(pnode->parent == parent,
                "Invalid tree node " << pnode->node_type);
    for (size_type i = 0; i < pnode->children.size(); ++i)
      verify_tree(pnode->children[i], pnode);
  }


#define ga_valid_operand(expr, pnode)                        \
  {                                                          \
    if (pnode && (pnode->node_type == GA_NODE_PREDEF_FUNC || \
                  pnode->node_type == GA_NODE_SPEC_FUNC ||   \
                  pnode->node_type == GA_NODE_NAME ||        \
                  pnode->node_type == GA_NODE_OPERATOR ||    \
                  pnode->node_type == GA_NODE_ALLINDICES))   \
      ga_throw_error(expr, pnode->pos, "Invalid term");      \
  }

  static void ga_print_constant_tensor(pga_tree_node pnode,
                                       std::ostream &str) {
    size_type nt = pnode->nb_test_functions(); // for printing zero tensors
    switch (pnode->tensor_order()) {
    case 0:
      str << (nt ? scalar_type(0) : pnode->t[0]);
      break;

    case 1:
      str << "[";
      for (size_type i = 0; i < pnode->tensor_proper_size(0); ++i) {
        if (i != 0) str << "; ";
        str << (nt ? scalar_type(0) : pnode->t[i]);
      }
      str << "]";
      break;

    case 2: case 3: case 4:
      {
        size_type ii(0);
        size_type n0 = pnode->tensor_proper_size(0);
        size_type n1 = pnode->tensor_proper_size(1);
        size_type n2 = (pnode->tensor_order() > 2) ? pnode->tensor_proper_size(2) : 1;
        size_type n3 = (pnode->tensor_order() > 3) ? pnode->tensor_proper_size(3) : 1;
        if (n3 > 1) str << "[";
        for (size_type l = 0; l < n3; ++l) {
          if (l != 0) str << ",";
          if (n2 > 1) str << "[";
          for (size_type k = 0; k < n2; ++k) {
            if (k != 0) str << ",";
            if (n1 > 1) str << "[";
            for (size_type j = 0; j < n1; ++j) {
              if (j != 0) str << ",";
              if (n0 > 1) str << "[";
              for (size_type i = 0; i < n0; ++i) {
                if (i != 0) str << ",";
                str << (nt ? scalar_type(0) : pnode->t[ii++]);
              }
              if (n0 > 1) str << "]";
            }
            if (n1 > 1) str << "]";
          }
          if (n2 > 1) str << "]";
        }
        if (n3 > 1) str << "]";
      }
      break;

    case 5: case 6:
      str << "Reshape([";
      for (size_type i = 0; i < pnode->tensor_proper_size(); ++i) {
        if (i != 0) str << "; ";
        str << (nt ? scalar_type(0) : pnode->t[i]);
      }
      str << "]";
      for (size_type i = 0; i < pnode->tensor_order(); ++i) {
        if (i != 0) str << ", ";
        str << pnode->tensor_proper_size(i);
      }
      str << ")";
      break;

    default: GMM_ASSERT1(false, "Invalid tensor dimension");
    }
    GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
  }

  static void ga_print_node(pga_tree_node pnode,
                            std::ostream &str) {
    long prec = str.precision(16);

    bool is_interpolate(false), is_elementary(false);
    bool is_xfem_plus(false), is_xfem_minus(false);
    switch(pnode->node_type) {
    case GA_NODE_INTERPOLATE:
    case GA_NODE_INTERPOLATE_FILTER:
    case GA_NODE_INTERPOLATE_NORMAL:
    case GA_NODE_INTERPOLATE_X:
    case GA_NODE_INTERPOLATE_VAL:
    case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS:
    case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_INTERPOLATE_VAL_TEST:
    case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST:
    case GA_NODE_INTERPOLATE_DIVERG_TEST:
      str << "Interpolate(";
      is_interpolate = true;
      break;
    case GA_NODE_ELEMENTARY:
    case GA_NODE_ELEMENTARY_VAL:
    case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS:
    case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_ELEMENTARY_VAL_TEST:
    case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST:
    case GA_NODE_ELEMENTARY_DIVERG_TEST:
      is_elementary = true;
      str << "Elementary_transformation(";
      break;
    case GA_NODE_XFEM_PLUS:
    case GA_NODE_XFEM_PLUS_VAL:
    case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS:
    case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL_TEST:
    case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST:
    case GA_NODE_XFEM_PLUS_DIVERG_TEST:
      is_xfem_plus = true;
      str << "Xfem_plus(";
      break;
    case GA_NODE_XFEM_MINUS:
    case GA_NODE_XFEM_MINUS_VAL:
    case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS:
    case GA_NODE_XFEM_MINUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL_TEST:
    case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST:
    case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      is_xfem_minus = true;
      str << "Xfem_minus(";
      break;
    default:
      break;
    }

    switch(pnode->node_type) {
    case GA_NODE_GRAD:
    case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_GRAD_TEST:
      str << "Grad_";
      break;
    case GA_NODE_HESS:
    case GA_NODE_INTERPOLATE_HESS:
    case GA_NODE_ELEMENTARY_HESS:
    case GA_NODE_XFEM_PLUS_HESS:
    case GA_NODE_XFEM_MINUS_HESS:
    case GA_NODE_HESS_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST:
      str << "Hess_";
      break;
    case GA_NODE_DIVERG:
    case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_DIVERG:
    case GA_NODE_DIVERG_TEST:
    case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      str << "Div_";
      break;
    default:
      break;
    }

    switch(pnode->node_type) {
    case GA_NODE_OP:
      {
        bool par = false;
        if (pnode->parent) {
          if (pnode->parent->node_type == GA_NODE_OP &&
              (ga_operator_priorities[pnode->op_type] >= 2 ||
               ga_operator_priorities[pnode->op_type]
               < ga_operator_priorities[pnode->parent->op_type]))
            par = true;
          if (pnode->parent->node_type == GA_NODE_PARAMS) par = true;
        }


        if (par) str << "(";
        if (pnode->op_type == GA_UNARY_MINUS) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          str << "-"; ga_print_node(pnode->children[0], str);
        } else if (pnode->op_type == GA_QUOTE) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          ga_print_node(pnode->children[0], str); str << "'";
        } else if (pnode->op_type == GA_TRACE) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          str << "Trace("; ga_print_node(pnode->children[0], str); str << ")";
        } else if (pnode->op_type == GA_DEVIATOR) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree with "
                      << pnode->children.size() << " children instead of 1");
          str << "Deviator("; ga_print_node(pnode->children[0], str); str<<")";
        } else if (pnode->op_type == GA_PRINT) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          str << "Print("; ga_print_node(pnode->children[0], str); str << ")";
        } else {
          // GMM_ASSERT1(pnode->children.size() == 2, "Invalid tree");
          if (!par && pnode->op_type == GA_MULT &&
              (pnode->children.size() == 1 ||
               pnode->test_function_type == size_type(-1) ||
               (pnode->children[0]->tensor_order() == 4 &&
                pnode->children[1]->tensor_order() == 2)))
            { par = true; str << "("; }
          ga_print_node(pnode->children[0], str);
          switch (pnode->op_type) {
          case GA_PLUS: str << "+"; break;
          case GA_MINUS: str << "-"; break;
          case GA_MULT: str << "*"; break;
          case GA_DIV: str << "/"; break;
          case GA_COLON: str << ":"; break;
          case GA_DOT: str << "."; break;
          case GA_DOTMULT: str << ".*"; break;
          case GA_DOTDIV: str << "./"; break;
          case GA_TMULT: str << "@"; break;
          default: GMM_ASSERT1(false, "Invalid or not taken into account "
                               "operation");
          }
          if (pnode->children.size() >= 2)
            ga_print_node(pnode->children[1], str);
          else
            str << "(unknown second argument)";
        }
        if (par) str << ")";
      }
      break;

    case GA_NODE_X:
      if (pnode->nbc1) str << "X(" << pnode->nbc1 << ")"; else str << "X";
      break;
    case GA_NODE_ELT_SIZE: str << "element_size"; break;
    case GA_NODE_ELT_K: str << "element_K"; break;
    case GA_NODE_ELT_B: str << "element_B"; break;
    case GA_NODE_NORMAL: str << "Normal"; break;
    case GA_NODE_INTERPOLATE_FILTER:
      str << "Interpolate_filter(" << pnode->interpolate_name << ",";
      ga_print_node(pnode->children[0], str);
      if (pnode->children.size() == 2)
        {  str << ","; ga_print_node(pnode->children[1], str); }
      else if (pnode->nbc1 != size_type(-1)) str << "," << pnode->nbc1;
      str << ")";
      break;
    case GA_NODE_INTERPOLATE_NORMAL:
      str << "Normal";
      break;
    case GA_NODE_INTERPOLATE_X:
      str << "X";
      break;
    case GA_NODE_INTERPOLATE_DERIVATIVE:
      if (pnode->test_function_type == 1) str << "Test_"; else str << "Test2_";
      str << "Interpolate_derivative(" << pnode->interpolate_name_der << ","
          << pnode->interpolate_name << "," << pnode->name << ")";
      break;
    case GA_NODE_INTERPOLATE:
    case GA_NODE_ELEMENTARY:
    case GA_NODE_XFEM_PLUS:
    case GA_NODE_XFEM_MINUS:
    case GA_NODE_VAL:
    case GA_NODE_INTERPOLATE_VAL:
    case GA_NODE_ELEMENTARY_VAL:
    case GA_NODE_XFEM_PLUS_VAL:
    case GA_NODE_XFEM_MINUS_VAL:
    case GA_NODE_GRAD:
    case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_HESS:
    case GA_NODE_INTERPOLATE_HESS:
    case GA_NODE_ELEMENTARY_HESS:
    case GA_NODE_XFEM_PLUS_HESS:
    case GA_NODE_XFEM_MINUS_HESS:
    case GA_NODE_DIVERG:
    case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_DIVERG:
      str << pnode->name;
      break;
    case GA_NODE_VAL_TEST:
    case GA_NODE_INTERPOLATE_VAL_TEST:
    case GA_NODE_ELEMENTARY_VAL_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST:
    case GA_NODE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_HESS_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST:
    case GA_NODE_DIVERG_TEST:
    case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      if (pnode->test_function_type == 1) str << "Test_"; else str << "Test2_";
      str << pnode->name;
      break;
    case GA_NODE_SPEC_FUNC: str << pnode->name; break;
    case GA_NODE_OPERATOR:
    case GA_NODE_PREDEF_FUNC:
      if (pnode->der1) {
        str << "Derivative_" << pnode->der1 << "_";
        if (pnode->der2) str << pnode->der2 << "_";
      }
      str << pnode->name; break;
    case GA_NODE_ZERO:
      GMM_ASSERT1(pnode->test_function_type != size_type(-1),
                  "Internal error");
      if (pnode->test_function_type) str << "(";
      ga_print_constant_tensor(pnode, str);
      if (pnode->name_test1.size()) {
        GMM_ASSERT1(pnode->qdim1 > 0, "Internal error");
        if (pnode->qdim1 == 1)
          str << "*Test_" << pnode->name_test1;
        else {
          str << "*(Reshape(Test_" << pnode->name_test1 << ","
              << pnode->qdim1<< ")(1))";
        }
      }
      if (pnode->name_test2.size()) {
        GMM_ASSERT1(pnode->qdim2 > 0, "Internal error");
        if (pnode->qdim2 == 1)
          str << "*Test2_" << pnode->name_test2;
        else {
          str << "*(Reshape(Test2_" << pnode->name_test2 << ","
              << pnode->qdim2<< ")(1))";
        }
      }
      if (pnode->test_function_type) str << ")";
      break;

    case GA_NODE_CONSTANT:
      ga_print_constant_tensor(pnode, str);
      break;

    case GA_NODE_ALLINDICES:
      str << ":";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_PARAMS:
      GMM_ASSERT1(pnode->children.size(), "Invalid tree");
      ga_print_node(pnode->children[0], str);
      str << "(";
      for (size_type i = 1; i < pnode->children.size(); ++i)
        { if (i > 1) str << ", "; ga_print_node(pnode->children[i], str); }
      str << ")";
      break;

    case GA_NODE_NAME:
      str << pnode->name;
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_RESHAPE:
      str << "Reshape";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_C_MATRIX:
      {
        GMM_ASSERT1(pnode->children.size(), "Invalid tree");
        size_type nbc1 = pnode->nbc1;
        size_type nbc2 = pnode->nbc2;
        size_type nbc3 = pnode->nbc3;
        size_type nbcl = pnode->children.size()/(nbc1*nbc2*nbc3);
        if (nbc1 > 1) str << "[";
        for (size_type i1 = 0; i1 < nbc1; ++i1) {
          if (i1 != 0) str << ",";
          if (nbc2 > 1) str << "[";
          for (size_type i2 = 0; i2 < nbc2; ++i2) {
            if (i2 != 0) str << ",";
            if (nbc3 > 1) str << "[";
            for (size_type i3 = 0; i3 < nbc3; ++i3) {
              if (i3 != 0) str << ",";
              if (nbcl > 1) str << "[";
              for (size_type i4 = 0; i4 < nbcl; ++i4) {
                if (i4 != 0) str << ",";
                size_type ii = ((i4*nbc3 + i3)*nbc2 + i2)*nbc1 + i1;
                ga_print_node(pnode->children[ii], str);
              }
              if (nbcl > 1) str << "]";
            }
            if (nbc3 > 1) str << "]";
          }
          if (nbc2 > 1) str << "]";
        }
        if (nbc1 > 1) str << "]";
      }
      break;

    default:
      str << "Invalid or not taken into account node type "
           << pnode->node_type;
      break;
    }

    if (is_interpolate)
      str << "," << pnode->interpolate_name << ")";
    else if (is_elementary)
      str << "," << pnode->elementary_name << ")";
    else if (is_xfem_plus || is_xfem_minus)
      str << ")";

    str.precision(prec);
  }

  std::string ga_tree_to_string(const ga_tree &tree) {
    std::stringstream str;
    str.precision(16);
    if (tree.root) verify_tree(tree.root, 0);
    if (tree.root) ga_print_node(tree.root, str); else str << "0";
    return str.str();
  }


  size_type ga_parse_prefix_operator(std::string &name) {
    if (name.size() >= 5 && name.compare(0, 5, "Grad_") == 0)
      { name = name.substr(5); return 1; }
    else if (name.size() >= 5 && name.compare(0, 5, "Hess_") == 0)
      { name = name.substr(5); return 2; }
    else if (name.size() >= 4 && name.compare(0, 4, "Div_") == 0)
      { name = name.substr(4); return 3; }
    return 0;
  }

  size_type ga_parse_prefix_test(std::string &name) {
    if (name.size() >= 5 && name.compare(0, 5, "Test_") == 0)
      { name = name.substr(5); return 1; }
    else if (name.size() >= 6 && name.compare(0, 6, "Test2_") == 0)
      { name = name.substr(6); return 2; }
    return 0;
  }


  //=========================================================================
  // Syntax analysis for the generic assembly langage
  //=========================================================================

  // Read a term with an (implicit) pushdown automaton.
  static GA_TOKEN_TYPE ga_read_term(const std::string &expr, size_type &pos,
                                    ga_tree &tree) {
    size_type token_pos, token_length;
    GA_TOKEN_TYPE t_type;
    int state = 1; // 1 = reading term, 2 = reading after term

    for (;;) {

      t_type = ga_get_token(expr, pos, token_pos, token_length);

      switch (state) {

      case 1:
        switch (t_type) {
        case GA_SCALAR:
          {
            char *endptr; const char *nptr = &(expr[token_pos]);
            scalar_type s_read = ::strtod(nptr, &endptr);
            if (endptr == nptr)
              ga_throw_error(expr, token_pos, "Bad numeric format.");
            tree.add_scalar(s_read, token_pos);
          }
          state = 2; break;

        case GA_COLON:
          tree.add_allindices(token_pos);
          state = 2; break;

        case GA_NAME:
          tree.add_name(&(expr[token_pos]), token_length, token_pos);
          state = 2; break;

        case GA_MINUS: // unary -
          tree.add_op(GA_UNARY_MINUS, token_pos);
        case GA_PLUS:  // unary +
          state = 1; break;

        case GA_TRACE:
          tree.add_op(GA_TRACE, token_pos);
          state = 1; break;

        case GA_DEVIATOR:
          tree.add_op(GA_DEVIATOR, token_pos);
          state = 1; break;

        case GA_INTERPOLATE:
          {
            tree.add_scalar(scalar_type(0), token_pos);
            tree.current_node->node_type = GA_NODE_INTERPOLATE;
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1, "Missing interpolate arguments.");
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "First argument of Interpolate should be a "
                             "variable, test function, X or Normal.");
            tree.current_node->name = std::string(&(expr[token_pos]),
                                                  token_length);

            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_COMMA)
              ga_throw_error(expr, pos, "Bad format for Interpolate "
                             "arguments.");
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "Second argument of Interpolate should be a "
                             "transformation name.");
            tree.current_node->interpolate_name
              = std::string(&(expr[token_pos]), token_length);
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "interpolate arguments.");
            state = 2;
          }
          break;

        case GA_ELEMENTARY:
          {
            tree.add_scalar(scalar_type(0), token_pos);
            tree.current_node->node_type = GA_NODE_ELEMENTARY;
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,
                             "Missing Elementary_transformation arguments.");
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "First argument of Elementary_transformation should be a "
                             "variable or a test function.");
            tree.current_node->name = std::string(&(expr[token_pos]),
                                                  token_length);

            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_COMMA)
              ga_throw_error(expr, pos, "Bad format for Elementary_transformation"
                             "arguments.");
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "Second argument of Elementary_transformation should "
                             "be a transformation name.");
            tree.current_node->elementary_name
              = std::string(&(expr[token_pos]), token_length);
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Elementary_transformation arguments.");
            state = 2;
          }
          break;

        case GA_XFEM_PLUS:
          {
            tree.add_scalar(scalar_type(0), token_pos);
            tree.current_node->node_type = GA_NODE_XFEM_PLUS;
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,
                             "Missing Xfem_plus arguments.");
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "The argument of Xfem_plus should be a "
                             "variable or a test function.");
            tree.current_node->name = std::string(&(expr[token_pos]),
                                                  token_length);
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Xfem_plus argument.");
            state = 2;
          }
          break;

        case GA_XFEM_MINUS:
          {
            tree.add_scalar(scalar_type(0), token_pos);
            tree.current_node->node_type = GA_NODE_XFEM_MINUS;
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,
                             "Missing Xfem_minus arguments.");
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "The argument of Xfem_minus should be a "
                             "variable or a test function.");
            tree.current_node->name = std::string(&(expr[token_pos]),
                                                  token_length);
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Xfem_minus argument.");
            state = 2;
          }
          break;

        case GA_INTERPOLATE_FILTER:
          {
            tree.add_scalar(scalar_type(0), token_pos);
            tree.current_node->node_type = GA_NODE_INTERPOLATE_FILTER;
            tree.current_node->nbc1 = size_type(-1);
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1, "Missing interpolate arguments.");
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos, "First argument of Interpolate_filter "
                             "should be a transformation name.");
            tree.current_node->interpolate_name
              = std::string(&(expr[token_pos]), token_length);
            t_type = ga_get_token(expr, pos, token_pos, token_length);
            if (t_type != GA_COMMA)
              ga_throw_error(expr,pos,
                             "Bad format for Interpolate_filter arguments.");
            ga_tree sub_tree;
            t_type = ga_read_term(expr, pos, sub_tree);
            if (t_type != GA_RPAR && t_type != GA_COMMA)
              ga_throw_error(expr, pos-1,
                             "Bad format for Interpolate_filter arguments.");
            tree.add_sub_tree(sub_tree);
            if (t_type == GA_COMMA) {
               ga_tree sub_tree2;
               t_type = ga_read_term(expr, pos, sub_tree2);
               tree.add_sub_tree(sub_tree2);
            }
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Unbalanced parenthesis.");
            state = 2;
          }
          break;

        case GA_PRINT:
          tree.add_op(GA_PRINT, token_pos);
          state = 1; break;

        case GA_LPAR: // Parenthesed expression
          {
            ga_tree sub_tree;
            GA_TOKEN_TYPE r_type;
            r_type = ga_read_term(expr, pos, sub_tree);
            if (r_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Unbalanced parenthesis.");
            tree.add_sub_tree(sub_tree);
            state = 2;
          }
          break;

        case GA_LBRACKET: // Explicit vector/matrix or tensor
          {
            ga_tree sub_tree;
            GA_TOKEN_TYPE r_type;
            size_type nbc1(0), nbc2(0), nbc3(0), n1(0), n2(0), n3(0);
            size_type tensor_order(1);
            bool foundcomma(false), foundsemi(false), nested_format(false);

            tree.add_matrix(token_pos);
            do {
              r_type = ga_read_term(expr, pos, sub_tree);

              if (sub_tree.root->node_type == GA_NODE_C_MATRIX) {

                // nested format

                if (r_type != GA_COMMA && r_type != GA_RBRACKET)  // in the nested format only "," and "]" are expected
                  ga_throw_error(expr, pos-1, "Bad explicit "
                                 "vector/matrix/tensor format.")
                else if (sub_tree.root->nbc3 != 1)                // the sub-tensor to be merged cannot be of fourth order
                  ga_throw_error(expr, pos-1, "Definition of explicit "
                                 "tensors is limitted to the fourth order. "
                                 "Limit exceeded.")
                else if (foundsemi ||                             // Cannot mix with the non-nested format.
                         (sub_tree.root->children.size() > 1 &&   // The sub-tensor cannot be a column vector [a;b],
                          sub_tree.root->nbc1 == 1))              // the nested format only accepts row vectors [a,b]
                  ga_throw_error(expr, pos-1, "Bad explicit "     // and converts them to column vectors internally
                                 "vector/matrix/tensor format.")  // (see below)

                if (sub_tree.root->children.size() == sub_tree.root->nbc1)
                  sub_tree.root->nbc1 = 1; // convert a row vector [a,b] to a column vector [a;b]

                nested_format = true;

                size_type sub_tensor_order = 3;
                if (sub_tree.root->nbc1 == 1)
                  sub_tensor_order = 1;
                else if (sub_tree.root->nbc2 == 1)
                  sub_tensor_order = 2;

                if (tensor_order == 1) {
                  if (nbc1 != 0 || nbc2 != 0 || nbc3 != 0)
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                  nbc1 = 1;
                  nbc2 = sub_tree.root->nbc1;
                  nbc3 = sub_tree.root->nbc2;
                  tensor_order = sub_tensor_order + 1;
                } else {
                  if ((tensor_order != sub_tensor_order + 1) ||
                      (tensor_order > 2 && nbc2 != sub_tree.root->nbc1) ||
                      (tensor_order > 3 && nbc3 != sub_tree.root->nbc2))
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                  nbc1 += 1;
                }

                tree.zip_matrix(sub_tree.root);
                sub_tree.clear();
                if (r_type == GA_RBRACKET) {
                  tree.current_node->nbc1 = nbc1;
                  tree.current_node->nbc2 = nbc2;
                  tree.current_node->nbc3 = nbc3;
                }

              } else {

                // non-nested format

                tree.add_sub_tree(sub_tree);

                ++n1; ++n2; ++n3;
                if (tensor_order < 2) ++nbc1;
                if (tensor_order < 3) ++nbc2;
                if (tensor_order < 4) ++nbc3;

                if (r_type == GA_COMMA) {
                  if (!foundcomma && tensor_order > 1)
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                  foundcomma = true;
                } else if (r_type == GA_SEMICOLON) {
                  if (n1 != nbc1)
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                  n1 = 0;
                  tensor_order = std::max(tensor_order, size_type(2));
                } else if (r_type == GA_DCOMMA) {
                  if (n1 != nbc1 || n2 != nbc2)
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                  foundsemi = true;
                  n2 = n1 = 0;
                  tensor_order = std::max(tensor_order, size_type(3));
                } else if (r_type == GA_DSEMICOLON) {
                  if (n1 != nbc1 || n2 != nbc2 || n3 != nbc3 ||
                      tensor_order < 3)
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                  n3 = n2 = n1 = 0;
                  tensor_order = std::max(tensor_order, size_type(4));
                } else if (r_type == GA_RBRACKET) {
                  if (n1 != nbc1 || n2 != nbc2 || n3 != nbc3 ||
                      tensor_order == 3)
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                  tree.current_node->nbc1 = nbc1;
                  if (tensor_order == 4) {
                    tree.current_node->nbc2 = nbc2/nbc1;
                    tree.current_node->nbc3 = nbc3/nbc2;
                  } else {
                    tree.current_node->nbc2 = tree.current_node->nbc3 = 1;
                  }
                } else {
                  ga_throw_error(expr, pos-1, "The explicit "
                                 "vector/matrix/tensor components should be "
                                 "separated by ',', ';', ',,' and ';;' and "
                                 "be ended by ']'.");
                }
              }

            } while (r_type != GA_RBRACKET);

            state = 2;
          }
          break;

        default: ga_throw_error(expr, token_pos, "Unexpected token.");
        }
        break;

      case 2:
        switch (t_type) {
        case GA_PLUS: case GA_MINUS: case GA_MULT: case GA_DIV:
        case GA_COLON: case GA_DOT: case GA_DOTMULT: case GA_DOTDIV:
        case GA_TMULT:
          tree.add_op(t_type, token_pos);
          state = 1; break;
        case GA_QUOTE:
          tree.add_op(t_type, token_pos);
          state = 2; break;
        case GA_END: case GA_RPAR: case GA_COMMA: case GA_DCOMMA:
        case GA_RBRACKET: case GA_SEMICOLON: case GA_DSEMICOLON:
          return t_type;
        case GA_LPAR: // Parameter list
          {
            ga_tree sub_tree;
            GA_TOKEN_TYPE r_type;
            tree.add_params(token_pos);
            do {
              r_type = ga_read_term(expr, pos, sub_tree);
              if (r_type != GA_RPAR && r_type != GA_COMMA)
                ga_throw_error(expr, pos-((r_type != GA_END)?1:0),
                               "Parameters should be separated "
                               "by ',' and parameter list ended by ')'.");
              tree.add_sub_tree(sub_tree);
            } while (r_type != GA_RPAR);
            state = 2;
          }
          break;

        default: ga_throw_error(expr, token_pos, "Unexpected token.");
        }
        break;
      }
    }

    return GA_INVALID;
  }

  // Syntax analysis of a string. Conversion to a tree.
  void ga_read_string(const std::string &expr, ga_tree &tree) {
    size_type pos = 0, token_pos, token_length;
    tree.clear();
    GA_TOKEN_TYPE t = ga_get_token(expr, pos, token_pos, token_length);
    if (t == GA_END) return;
    pos = 0;
    t = ga_read_term(expr, pos, tree);
    switch (t) {
    case GA_RPAR: ga_throw_error(expr, pos-1, "Unbalanced parenthesis.");
    case GA_RBRACKET: ga_throw_error(expr, pos-1, "Unbalanced braket.");
    case GA_END: break;
    default: ga_throw_error(expr, pos-1, "Unexpected token.");
    }
  }

  //=========================================================================
  // Structure to gather instructions
  //=========================================================================

  class ga_if_hierarchy : public std::vector<size_type> {

  public:
    void increment(void) { (back())++; }
    void child_of(const ga_if_hierarchy &gih)
    { *this = gih; push_back(0); }
    bool is_compatible(const std::list<ga_if_hierarchy> &gihl) {

      std::list<ga_if_hierarchy>::const_iterator it = gihl.begin();
      for (; it != gihl.end(); ++it) {
        if (it->size() <= size()) {
          bool ok = true;
          for (size_type i = 0; i+1 < it->size(); ++i)
            if ((*it)[i] != (*this)[i]) { ok = false; break; }
          if (it->back() > (*this)[it->size()-1]) { ok = false; break; }
          if (ok) return true;
        }
      }
      return false;
    }

    ga_if_hierarchy(void) : std::vector<size_type>(1) { (*this)[0] = 0; }
  };


  struct ga_instruction {
    virtual int exec(void) = 0;
    virtual ~ga_instruction() {};
  };

  typedef ga_instruction *pga_instruction;
  typedef std::vector<pga_instruction> ga_instruction_list;

  struct ga_instruction_set {

    papprox_integration pai;       // Current approximation method
    fem_interpolation_context ctx; // Current fem interpolation context.
    base_small_vector Normal;      // Outward unit normal vector to the
                                   // boundary in case of boundary integration
    scalar_type elt_size;          // Estimate of the diameter of the element
                                   // if needed.
    bool need_elt_size;
    scalar_type coeff;             // Coefficient for the Gauss point
    size_type nbpt, ipt;           // Number and index of Gauss point
    bgeot::geotrans_precomp_pool gp_pool;
    fem_precomp_pool fp_pool;

    struct region_mim : std::pair<const mesh_im *, const mesh_region *> {
      const mesh_im* mim(void) const { return this->first; }
      const mesh_region* region(void) const { return this->second; }
      region_mim(const mesh_im *mim_, const mesh_region *region_) :
        std::pair<const mesh_im *, const mesh_region *>(mim_, region_) {}
    };

    std::map<std::string, const base_vector *> extended_vars;
    std::map<std::string, base_vector> really_extended_vars;
    std::map<std::string, gmm::sub_interval> var_intervals;
    size_type nb_dof, max_dof;

    struct variable_group_info {
      const mesh_fem *mf;
      gmm::sub_interval Ir, In;
      scalar_type alpha;
      const base_vector *U;
      const std::string *varname;
      variable_group_info(void) : mf(0), U(0), varname(0) {}
    };

    struct interpolate_info {
      size_type pt_type;
      bool has_ctx;
      const mesh *m;
      fem_interpolation_context ctx;
      base_node pt_y;
      base_small_vector Normal;
      std::map<std::string, variable_group_info> groups_info;
      std::map<var_trans_pair, base_tensor> derivatives;
    };

    struct elementary_trans_info {
      base_matrix M;
      const mesh_fem *mf;
      size_type icv;
    };

    std::set<std::string> transformations;

    struct region_mim_instructions {

      const mesh *m;
      ga_if_hierarchy current_hierarchy;
      std::map<std::string, base_vector> local_dofs;
      std::map<std::string, std::list<ga_if_hierarchy> > local_dofs_hierarchy;
      std::map<const mesh_fem *, pfem_precomp> pfps;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy> > pfps_hierarchy;
      std::map<const mesh_fem *, base_tensor> base;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy> > base_hierarchy;
      std::map<const mesh_fem *, base_tensor> grad;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy> > grad_hierarchy;
      std::map<const mesh_fem *, base_tensor> hess;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy> > hess_hierarchy;

      std::map<const mesh_fem *, base_tensor> xfem_plus_base;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy> > xfem_plus_base_hierarchy;
      std::map<const mesh_fem *, base_tensor> xfem_plus_grad;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy> > xfem_plus_grad_hierarchy;
      std::map<const mesh_fem *, base_tensor> xfem_plus_hess;
      std::map<const mesh_fem *, std::list<ga_if_hierarchy> > xfem_plus_hess_hierarchy;
      std::map<const mesh_fem *, base_tensor> xfem_minus_base;
      std::map<const mesh_fem *,std::list<ga_if_hierarchy> > xfem_minus_base_hierarchy;
      std::map<const mesh_fem *, base_tensor> xfem_minus_grad;
      std::map<const mesh_fem *,std::list<ga_if_hierarchy> > xfem_minus_grad_hierarchy;
      std::map<const mesh_fem *, base_tensor> xfem_minus_hess;
      std::map<const mesh_fem *,std::list<ga_if_hierarchy> > xfem_minus_hess_hierarchy;

      std::map<std::string, std::set<std::string> > transformations;
      std::set<std::string> transformations_der;
      std::map<std::string, interpolate_info> interpolate_infos;
      std::map<std::string, elementary_trans_info> elementary_trans_infos;

      ga_instruction_list instructions;
      std::map<scalar_type, std::list<pga_tree_node> > node_list;

      ~region_mim_instructions(void) {
        for (size_type i = 0; i < instructions.size(); ++i) {
          delete instructions[i];
          instructions[i] = 0;
        }
      }
      region_mim_instructions(void): m(0) {}
    };

    std::list<ga_tree> trees; // The trees are stored mainly because they
                              // contain the intermediary tensors.

    typedef std::map<region_mim,  region_mim_instructions> instructions_set;

    instructions_set  whole_instructions;

    ga_instruction_set(void) { max_dof = nb_dof = 0; need_elt_size = false; }
  };


  //=========================================================================
  // Structure dealing with predefined special functions
  // such as mesh_dim, pi, qdim ...
  //=========================================================================

  typedef std::set<std::string> ga_spec_function_tab;
  static ga_spec_function_tab SPEC_FUNCTIONS;

  //=========================================================================
  // Structure dealing with predefined scalar functions.
  //=========================================================================

  static void ga_semantic_analysis(const std::string &, ga_tree &,
                                   const ga_workspace &, size_type, size_type,
                                   bool, bool, int option = 0);
  static void ga_node_analysis(const std::string &, ga_tree &,
                               const ga_workspace &, pga_tree_node, size_type,
                               size_type, bool, bool, int);
  static void ga_derivative(ga_tree &, const ga_workspace &, const mesh &,
                            const std::string &, const std::string &,
                            size_type);
  static void ga_exec(ga_instruction_set &gis, ga_workspace &workspace);
  static void ga_function_exec(ga_instruction_set &gis);
  static void ga_compile(ga_workspace &workspace, ga_instruction_set &gis,
                         size_type order);
  static void ga_compile_function(ga_workspace &workspace,
                                  ga_instruction_set &gis, bool scalar);
  static std::string ga_derivative_scalar_function(const std::string expr,
                                                   const std::string &var);
  static bool ga_is_affine(ga_tree &tree, const ga_workspace &workspace,
                           const std::string &varname,
                           const std::string &interpolatename);

  using instruction_set = omp_distribute<ga_instruction_set>;

  class ga_predef_function {
    size_type ftype_; // 0 : C++ function with C++ derivative(s)
                      // 1 : function defined by an string expression.

    size_type dtype_; // 0 : no derivative(s)
                      // 1 : derivative(s) given by C++ functions
                      // 2 : derivatives(s) given by string expression(s)
                      // 3 : derivatives(s) to be symbolically computed.
    size_type nbargs_;         // One or two arguments
    pscalar_func_onearg f1_;   // Function pointer for a one argument function
    pscalar_func_twoargs f2_;  // Function pointer for a two arguments function
    std::string expr_;
    std::string derivative1_, derivative2_;
    mutable omp_distribute<base_vector> t, u;
    mutable omp_distribute<ga_workspace> workspace;
    copyable_ptr<instruction_set> gis;

    friend void ga_define_function(const std::string name, size_type nbargs,
                                   const std::string expr, const std::string der1,
                                   const std::string der2);
    friend void ga_define_function(const std::string name, pscalar_func_onearg f,
                                   const std::string &der);
    friend void ga_define_function(const std::string name, pscalar_func_twoargs f,
                                   const std::string &der1, const std::string &der2);
  public:
    scalar_type operator()(scalar_type t_, scalar_type u_ = 0.) const {
      switch(ftype_) {
      case 0:
        if (nbargs_ == 2)
          return (*f2_)(t_, u_);
        else
          return (*f1_)(t_);
        break;
      case 1:
        t.thrd_cast()[0] = t_; u.thrd_cast()[0] = u_;
        workspace.thrd_cast().assembled_potential() = scalar_type(0);
        ga_function_exec(*gis);
        return workspace.thrd_cast().assembled_potential();
        break;
      }
      return 0.;
    }

    bool is_affine(const std::string &varname) const {
      if (ftype_ == 1) {
        for (size_type i = 0; i < workspace.thrd_cast().nb_trees(); ++i) {
          ga_workspace::tree_description &td = workspace.thrd_cast().tree_info(i);
          if (!(ga_is_affine(*(td.ptree), workspace, varname, "")))
            return false;
        }
        return true;
      }
      return false;
    }

    size_type ftype() const {return ftype_;}
    size_type dtype() const {return dtype_;}
    size_type nbargs() const {return nbargs_;}
    const std::string &derivative1() const {return derivative1_;}
    const std::string &derivative2() const {return derivative2_;}
    const std::string &expr() const {return expr_;}
    pscalar_func_onearg f1() const {return f1_;}
    pscalar_func_twoargs f2() const {return f2_;}

    ga_predef_function(void) : expr_(""), derivative1_(""), derivative2_(""), gis(nullptr) {}
    ga_predef_function(pscalar_func_onearg f, size_type dtype__ = 0,
                       const std::string &der = "")
      : ftype_(0), dtype_(dtype__), nbargs_(1), f1_(f), expr_(""),
        derivative1_(der), derivative2_("") {}
    ga_predef_function(pscalar_func_twoargs f, size_type dtype__ = 0,
                       const std::string &der1 = "",
                       const std::string &der2 = "")
      : ftype_(0), dtype_(dtype__), nbargs_(2), f2_(f),
        expr_(""), derivative1_(der1), derivative2_(der2), gis(nullptr) {}
    ga_predef_function(const std::string &expr__)
      : ftype_(1), dtype_(3), nbargs_(1), expr_(expr__),
        derivative1_(""), derivative2_(""), t(1, 0.), u(1, 0.), gis(nullptr) {}
  };



  static scalar_type ga_Heaviside(scalar_type t) { return (t >= 0.) ? 1.: 0.; }
  static scalar_type ga_pos_part(scalar_type t) { return (t >= 0.) ? t : 0.; }
  static scalar_type ga_half_sqr_pos_part(scalar_type t)
  { return (t >= 0.) ? 0.5*t*t : 0.; }
  static scalar_type ga_neg_part(scalar_type t) { return (t >= 0.) ? 0. : -t; }
  static scalar_type ga_half_sqr_neg_part(scalar_type t)
  { return (t >= 0.) ? 0. : 0.5*t*t; }
  static scalar_type ga_sinc(scalar_type t) {// cardinal sine function sin(t)/t
    if (gmm::abs(t) < 1E-4) {
      scalar_type t2 = t*t;
      return 1-t2/6.+ t2*t2/120.;
    } else {
      return sin(t)/t;
    }
  }
  static scalar_type ga_sqr(scalar_type t) { return t*t; }
  static scalar_type ga_max(scalar_type t, scalar_type u)
  { return std::max(t,u); }
  static scalar_type ga_min(scalar_type t, scalar_type u)
  { return std::min(t,u); }
  static scalar_type ga_abs(scalar_type t) { return gmm::abs(t); }
  static scalar_type ga_sign(scalar_type t) { return (t >= 0.) ? 1.: -1.; }

  // Derivatives of predefined functions
  static scalar_type ga_der_sinc(scalar_type t) {
    if (gmm::abs(t) < 1E-4) {
      scalar_type t2 = t*t;
      return  -t/3. + t*t2/30. -t*t2*t2/840.;
    } else {
      return (t*cos(t) - sin(t))/(t*t);
    }
  }
  static scalar_type ga_der2_sinc(scalar_type t) {
    if (gmm::abs(t) < 1E-4) {
      scalar_type t2 = t*t;
      return  -1./3. + t2/10. -t2*t2/168.;
    } else {
      return ((2. - t*t)*sin(t) - 2.*t*cos(t))/(t*t*t);
    }
  }
  static scalar_type ga_der_sqrt(scalar_type t) { return 0.5/sqrt(t); }
  // static scalar_type ga_der_sqr(scalar_type t) { return 2*t; }
  static scalar_type ga_der_pow1(scalar_type t, scalar_type u)
  { return u*pow(t,u-1.); }
  static scalar_type ga_der_pow2(scalar_type t, scalar_type u)
  { return pow(t,u)*log(gmm::abs(t)); }
  static scalar_type ga_der_log(scalar_type t) { return 1./t; }
  static scalar_type ga_der_log10(scalar_type t) { return 1./(t*log(10.)); }
  static scalar_type ga_der_tanh(scalar_type t)
  { return 1.-gmm::sqr(tanh(t)); }
  static scalar_type ga_der_asinh(scalar_type t)
  { return 1./(sqrt(t*t+1.)); }
  static scalar_type ga_der_acosh(scalar_type t)
  { return 1./(sqrt(t*t-1.)); }
  static scalar_type ga_der_atanh(scalar_type t)
  { return 1./(1.-t*t); }
  static scalar_type ga_der_cos(scalar_type t)
  { return -sin(t); }
  static scalar_type ga_der_tan(scalar_type t)
  { return 1.+gmm::sqr(tan(t)); }
  static scalar_type ga_der_asin(scalar_type t)
  { return 1./(sqrt(1.-t*t)); }
  static scalar_type ga_der_acos(scalar_type t)
  { return -1./(sqrt(1.-t*t)); }
  static scalar_type ga_der_atan(scalar_type t)
  { return 1./(1.+t*t); }
  static scalar_type ga_der_erf(scalar_type t)
  { return exp(-t*t)*2./sqrt(M_PI); }
  static scalar_type ga_der_erfc(scalar_type t)
  { return -exp(-t*t)*2./sqrt(M_PI); }
  static scalar_type ga_der_neg_part(scalar_type t)
  { return (t >= 0) ? 0. : -1.; }
  static scalar_type ga_der_max1(scalar_type t, scalar_type u)
  { return (t-u >= 0) ? 1. : 0.; }
  static scalar_type ga_der_max2(scalar_type t, scalar_type u)
  { return (u-t >= 0) ? 1. : 0.; }


  typedef std::map<std::string, ga_predef_function> ga_predef_function_tab;
  static ga_predef_function_tab PREDEF_FUNCTIONS;

  //=========================================================================
  // Structure dealing with predefined operators.
  //=========================================================================

  static void ga_init_scalar(bgeot::multi_index &mi) { mi.resize(0); }
  static void ga_init_square_matrix(bgeot::multi_index &mi, size_type N)
  { mi.resize(2); mi[0] = mi[1] = N; }

  // Norm Operator
  struct norm_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() > 2) return false;
      ga_init_scalar(sizes);
      return true;
    }

    void value(const arg_list &args, base_tensor &result) const
    { result[0] = gmm::vect_norm2(args[0]->as_vector()); }

    // Derivative : u/|u|
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const {
      scalar_type no = gmm::vect_norm2(args[0]->as_vector());
      if (no == scalar_type(0))
        gmm::clear(result.as_vector());
      else
        gmm::copy(gmm::scaled(args[0]->as_vector(), scalar_type(1)/no),
                  result.as_vector());
    }

    // Second derivative : (|u|^2 Id - u x u)/|u|^3
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const {
      const base_tensor &t = *args[0];
      size_type N = t.size();
      scalar_type no = gmm::vect_norm2(t.as_vector());
      scalar_type no3 = no*no*no;

      if (no < 1E-25) no = 1E-25; // In order to avoid infinite values

      for (size_type i = 0; i < N; ++i)
        for (size_type j = 0; j < N; ++j) {
          result[j*N+i] = - t[i]*t[j] / no3;
          if (i == j) result[j*N+i] += scalar_type(1)/no;
        }
    }
  };

  // Norm_sqr Operator
  struct norm_sqr_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() > 2) return false;
      ga_init_scalar(sizes);
      return true;
    }

    void value(const arg_list &args, base_tensor &result) const
    { result[0] = gmm::vect_norm2_sqr(args[0]->as_vector()); }

    // Derivative : 2*u
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const {
      gmm::copy(gmm::scaled(args[0]->as_vector(), scalar_type(2)),
                result.as_vector());
    }

    // Second derivative : Id
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const {
      const base_tensor &t = *args[0];
      size_type N = t.size();
      gmm::clear(result.as_vector());
      for (size_type i = 0; i < N; ++i)
        result[i*N+i] = scalar_type(2);
    }
  };

  // Det Operator
  struct det_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() != 2
          || args[0]->sizes()[0] != args[0]->sizes()[1]) return false;
      ga_init_scalar(sizes);
      return true;
    }

    void value(const arg_list &args, base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix M(N, N);
      gmm::copy(args[0]->as_vector(), M.as_vector());
      result[0] = gmm::lu_det(M);

    }

    // Derivative : det(M)M^{-T}
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const { // to be verified
      size_type N = args[0]->sizes()[0];
      base_matrix M(N, N);
      gmm::copy(args[0]->as_vector(), M.as_vector());
      scalar_type det = gmm::lu_inverse(M);
      if (det == scalar_type(0))
        gmm::clear(result.as_vector());
      else {
        base_tensor::iterator it = result.begin();
        for (size_type j = 0; j < N; ++j)
          for (size_type i = 0; i < N; ++i, ++it)
            *it = M(j, i) * det;
        GA_DEBUG_ASSERT(it == result.end(), "Internal error");
      }
    }

    // Second derivative : det(M)(M^{-T}@M^{-T} - M^{-T}_{jk}M^{-T}_{li})
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const { // To be verified
      size_type N = args[0]->sizes()[0];
      base_matrix M(N, N);
      gmm::copy(args[0]->as_vector(), M.as_vector());
      scalar_type det = gmm::lu_inverse(M);
      if (det == scalar_type(0))
        gmm::clear(result.as_vector());
      else {
        base_tensor::iterator it = result.begin();
        for (size_type l = 0; l < N; ++l)
          for (size_type k = 0; k < N; ++k)
            for (size_type j = 0; j < N; ++j)
              for (size_type i = 0; i < N; ++i, ++it)
                *it = (M(j, i) * M(l,k) - M(j,k) * M(l, i)) * det;
        GA_DEBUG_ASSERT(it == result.end(), "Internal error");
      }
    }
  };




  // Inverse Operator (for square matrices)
  struct inverse_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() != 2
          || args[0]->sizes()[0] != args[0]->sizes()[1]) return false;
      ga_init_square_matrix(sizes, args[0]->sizes()[0]);
      return true;
    }

    // Value : M^{-1}
    void value(const arg_list &args, base_tensor &result) const {
      size_type N = args[0]->sizes()[0];
      base_matrix M(N, N);
      gmm::copy(args[0]->as_vector(), M.as_vector());
      gmm::lu_inverse(M);
      gmm::copy(M.as_vector(), result.as_vector());
    }

    // Derivative : -M^{-1}{ik}M^{-1}{lj}  (comes from H -> M^{-1}HM^{-1})
    void derivative(const arg_list &args, size_type,
                    base_tensor &result) const { // to be verified
      size_type N = args[0]->sizes()[0];
      base_matrix M(N, N);
      gmm::copy(args[0]->as_vector(), M.as_vector());
      gmm::lu_inverse(M);
      base_tensor::iterator it = result.begin();
      for (size_type l = 0; l < N; ++l)
        for (size_type k = 0; k < N; ++k)
          for (size_type j = 0; j < N; ++j)
            for (size_type i = 0; i < N; ++i, ++it)
              *it = -M(i,k)*M(l,j);
      GA_DEBUG_ASSERT(it == result.end(), "Internal error");
    }

    // Second derivative :
    // M^{-1}{ik}M^{-1}{lm}M^{-1}{nj} + M^{-1}{im}M^{-1}{mk}M^{-1}{lj}
    // comes from (H,K) -> M^{-1}HM^{-1}KM^{-1} + M^{-1}KM^{-1}HM^{-1}
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const { // to be verified
      size_type N = args[0]->sizes()[0];
      base_matrix M(N, N);
      gmm::copy(args[0]->as_vector(), M.as_vector());
      gmm::lu_inverse(M);
      base_tensor::iterator it = result.begin();
      for (size_type n = 0; n < N; ++n)
        for (size_type m = 0; m < N; ++m)
          for (size_type l = 0; l < N; ++l)
            for (size_type k = 0; k < N; ++k)
              for (size_type j = 0; j < N; ++j)
                for (size_type i = 0; i < N; ++i, ++it)
                  *it = M(i,k)*M(l,m)*M(n,j)+M(i,m)*M(m,k)*M(l,j);
      GA_DEBUG_ASSERT(it == result.end(), "Internal error");
    }
  };




  //=========================================================================
  // Initialization of predefined functions and operators.
  //=========================================================================


  bool init_predef_functions(void) {

    // Predefined functions

    // Power functions and their derivatives
    PREDEF_FUNCTIONS["sqrt"] = ga_predef_function(sqrt, 1, "DER_PDFUNC_SQRT");
    PREDEF_FUNCTIONS["sqr"] = ga_predef_function(ga_sqr, 2, "2*t");
    PREDEF_FUNCTIONS["pow"] = ga_predef_function(pow, 1, "DER_PDFUNC1_POW",
                                                 "DER_PDFUNC2_POW");
    PREDEF_FUNCTIONS["DER_PDFUNC_SQRT"] =
      ga_predef_function(ga_der_sqrt, 2, "-0.25/(t*sqrt(t))");
    PREDEF_FUNCTIONS["DER_PDFUNC1_POW"] =
      ga_predef_function(ga_der_pow1, 2, "u*(u-1)*pow(t,u-2)",
                         "pow(t,u-1)*(u*log(t)+1)");
    PREDEF_FUNCTIONS["DER_PDFUNC2_POW"] =
      ga_predef_function(ga_der_pow2, 2, "pow(t,u-1)*(u*log(t)+1)",
                         "pow(t,u)*sqr(log(t))");

    // Hyperbolic functions
    PREDEF_FUNCTIONS["exp"] = ga_predef_function(exp, 1, "exp");
    PREDEF_FUNCTIONS["log"] = ga_predef_function(log, 1, "DER_PDFUNC_LOG");
    PREDEF_FUNCTIONS["log10"] =
      ga_predef_function(log10, 1, "DER_PDFUNC_LOG10");
    PREDEF_FUNCTIONS["sinh"] = ga_predef_function(sinh, 1, "cosh");
    PREDEF_FUNCTIONS["cosh"] = ga_predef_function(cosh, 1, "sinh");
    PREDEF_FUNCTIONS["tanh"] = ga_predef_function(tanh, 1, "DER_PDFUNC_TANH");
    PREDEF_FUNCTIONS["asinh"] =
      ga_predef_function(asinh, 1, "DER_PDFUNC_ASINH");
    PREDEF_FUNCTIONS["acosh"] =
      ga_predef_function(acosh, 1, "DER_PDFUNC_ACOSH");
    PREDEF_FUNCTIONS["atanh"] =
      ga_predef_function(atanh, 1, "DER_PDFUNC_ATANH");


    PREDEF_FUNCTIONS["DER_PDFUNC_LOG"] =
      ga_predef_function(ga_der_log, 2, "-1/sqr(t)");
    PREDEF_FUNCTIONS["DER_PDFUNC_LOG10"] =
      ga_predef_function(ga_der_log10, 2, "-1/(sqr(t)*log(10))");
    PREDEF_FUNCTIONS["DER_PDFUNC_TANH"] =
      ga_predef_function(ga_der_tanh, 2, "2*tanh(t)*(sqr(tanh(t))-1)");
    PREDEF_FUNCTIONS["DER_PDFUNC_ASINH"] =
      ga_predef_function(ga_der_asinh, 2, "-t/(pow(t*t+1,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ACOSH"] =
      ga_predef_function(ga_der_acosh, 2, "-t/(pow(t*t-1,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ATANH"] =
      ga_predef_function(ga_der_atanh, 2, "-2*t/sqr(1+t*t)");


    // Trigonometric functions
    PREDEF_FUNCTIONS["sin"] = ga_predef_function(sin, 1, "cos");
    PREDEF_FUNCTIONS["cos"] = ga_predef_function(cos, 1, "DER_PDFUNC_COS");
    PREDEF_FUNCTIONS["tan"] = ga_predef_function(tan, 1, "DER_PDFUNC_TAN");
    PREDEF_FUNCTIONS["asin"] = ga_predef_function(asin, 1, "DER_PDFUNC_ASIN");
    PREDEF_FUNCTIONS["acos"] = ga_predef_function(acos, 1, "DER_PDFUNC_ACOS");
    PREDEF_FUNCTIONS["atan"] = ga_predef_function(atan, 1, "DER_PDFUNC_ATAN");
    PREDEF_FUNCTIONS["sinc"] = ga_predef_function(ga_sinc, 1,
                                                  "DER_PDFUNC_SINC");
    PREDEF_FUNCTIONS["DER_PDFUNC_SINC"] = ga_predef_function(ga_der_sinc, 1,
                                                             "DER2_PDFUNC_SINC");
    PREDEF_FUNCTIONS["DER2_PDFUNC_SINC"] = ga_predef_function(ga_der2_sinc);


    PREDEF_FUNCTIONS["DER_PDFUNC_COS"] =
      ga_predef_function(ga_der_cos, 2, "-cos(t)");
    PREDEF_FUNCTIONS["DER_PDFUNC_TAN"] =
      ga_predef_function(ga_der_tan, 2, "2*tan(t)/sqr(cos(t))");
    // PREDEF_FUNCTIONS["DER_PDFUNC_TAN"] =
    //  ga_predef_function(ga_der_tan, 2, "2*tan(t)*(1+sqr(tan(t)))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ASIN"] =
      ga_predef_function(ga_der_asin, 2, "t/(pow(1-t*t,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ACOS"] =
      ga_predef_function(ga_der_acos, 2, "-t/(pow(1-t*t,1.5))");
    PREDEF_FUNCTIONS["DER_PDFUNC_ATAN"] =
      ga_predef_function(ga_der_atan, 2, "2*t/sqr(1-t*t)");


    // Error functions
    PREDEF_FUNCTIONS["erf"]
      = ga_predef_function(erf, 1, "DER_PDFUNC_ERF");
    PREDEF_FUNCTIONS["erfc"]
      = ga_predef_function(erfc, 1, "DER_PDFUNC_ERFC");

    PREDEF_FUNCTIONS["DER_PDFUNC_ERF"] =
      ga_predef_function(ga_der_erf, 2, "exp(-t*t)*2/sqrt(pi)");
    PREDEF_FUNCTIONS["DER_PDFUNC_ERFC"] =
      ga_predef_function(ga_der_erfc, 2, "-exp(-t*t)*2/sqrt(pi)");



    // Miscellaneous functions
    PREDEF_FUNCTIONS["Heaviside"] = ga_predef_function(ga_Heaviside); // ga_predef_function(ga_Heaviside, 2, "(0)");
    PREDEF_FUNCTIONS["sign"] = ga_predef_function(ga_sign);
    PREDEF_FUNCTIONS["abs"] = ga_predef_function(ga_abs, 1, "sign");
    PREDEF_FUNCTIONS["pos_part"]
      = ga_predef_function(ga_pos_part, 1, "Heaviside");
    PREDEF_FUNCTIONS["half_sqr_pos_part"]
      = ga_predef_function(ga_half_sqr_pos_part, 1, "pos_part");
    PREDEF_FUNCTIONS["neg_part"]
      = ga_predef_function(ga_neg_part, 1, "DER_PDFUNC_NEG_PART");
    PREDEF_FUNCTIONS["half_sqr_neg_part"]
      = ga_predef_function(ga_half_sqr_neg_part, 2, "-neg_part(t)");

    PREDEF_FUNCTIONS["max"]
      = ga_predef_function(ga_max, 1, "DER_PDFUNC1_MAX", "DER_PDFUNC2_MAX");
    PREDEF_FUNCTIONS["min"]
      = ga_predef_function(ga_min, 1, "DER_PDFUNC2_MAX", "DER_PDFUNC1_MAX");

    PREDEF_FUNCTIONS["DER_PDFUNC_NEG_PART"]
      = ga_predef_function(ga_der_neg_part, 2, "-Heaviside(-t)");
    PREDEF_FUNCTIONS["DER_PDFUNC1_MAX"] = ga_predef_function(ga_der_max1);
    PREDEF_FUNCTIONS["DER_PDFUNC2_MAX"] = ga_predef_function(ga_der_max2);


    // Predefined special functions

    SPEC_FUNCTIONS.insert("pi");
    SPEC_FUNCTIONS.insert("meshdim");
    SPEC_FUNCTIONS.insert("qdim");
    SPEC_FUNCTIONS.insert("qdims");
    SPEC_FUNCTIONS.insert("Id");

    // Predefined operators
    ga_predef_operator_tab &PREDEF_OPERATORS
      = dal::singleton<ga_predef_operator_tab>::instance();

    PREDEF_OPERATORS.add_method("Norm", new norm_operator());
    PREDEF_OPERATORS.add_method("Norm_sqr", new norm_sqr_operator());
    PREDEF_OPERATORS.add_method("Det", new det_operator());
    PREDEF_OPERATORS.add_method("Inv", new inverse_operator());
    return true;
  }

  static bool predef_functions_initialized = init_predef_functions();

  bool ga_function_exists(const std::string name) {
    return PREDEF_FUNCTIONS.find(name) != PREDEF_FUNCTIONS.end();
  }


  void ga_define_function(const std::string name, size_type nbargs,
                          const std::string expr, const std::string der1,
                          const std::string der2) {
    GMM_ASSERT1(nbargs >= 1 && nbargs <= 2, "Generic assembly only allows "
                "the definition of scalar function with one or two arguments");
    { // Only for syntax analysis
      base_vector t(1);
      ga_workspace workspace;
      workspace.add_fixed_size_variable("t", gmm::sub_interval(0,1), t);
      if (nbargs == 2)
        workspace.add_fixed_size_variable("u", gmm::sub_interval(0,1), t);
      workspace.add_function_expression(expr);
    }

    GMM_ASSERT1(PREDEF_FUNCTIONS.find(name) == PREDEF_FUNCTIONS.end(),
                "Already defined function " << name);
    PREDEF_FUNCTIONS[name] = ga_predef_function(expr);
    ga_predef_function &F = PREDEF_FUNCTIONS[name];
    GMM_ASSERT1(!me_is_multithreaded_now(),
                "functions should not be defined in multi-threaded code");
    F.gis = std::unique_ptr<instruction_set>(new instruction_set());
    for (size_type thread = 0; thread < num_threads(); ++thread)
    {
      F.workspace(thread).add_fixed_size_variable("t", gmm::sub_interval(0,1), F.t(thread));
      if (nbargs == 2)
        F.workspace(thread).add_fixed_size_variable("u", gmm::sub_interval(0,1), F.u(thread));
      F.workspace(thread).add_function_expression(expr);
      ga_compile_function(F.workspace(thread), (*F.gis)(thread), true);
    }
    F.nbargs_ = nbargs;
    if (nbargs == 1) {
      if (der1.size()) { F.derivative1_ = der1; F.dtype_ = 2; }
    } else {
      if (der1.size() && der2.size()) {
        F.derivative1_ = der1;  F.derivative2_ = der2; F.dtype_ = 2;
      }
    }
  }

  void ga_define_function(const std::string name, pscalar_func_onearg f,
                          const std::string &der) {
    PREDEF_FUNCTIONS[name] = ga_predef_function(f, 1, der);
    ga_predef_function &F = PREDEF_FUNCTIONS[name];
    if (der.size() == 0) F.dtype_ = 0;
    else if (!(ga_function_exists(der))) F.dtype_ = 2;
  }

  void ga_define_function(const std::string name, pscalar_func_twoargs f,
                          const std::string &der1, const std::string &der2) {
    PREDEF_FUNCTIONS[name] = ga_predef_function(f, 1, der1, der2);
    ga_predef_function &F = PREDEF_FUNCTIONS[name];
    if (der1.size() == 0 || der2.size() == 0)
      F.dtype_ = 0;
    else if (!(ga_function_exists(der1)) || !(ga_function_exists(der2)))
      F.dtype_ = 2;
  }

  void ga_undefine_function(const std::string name) {
    ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
    if (it != PREDEF_FUNCTIONS.end()) {
      PREDEF_FUNCTIONS.erase(name);
      std::string name0 = "DER_PDFUNC_" + name;
      ga_undefine_function(name0);
      std::string name1 = "DER_PDFUNC1_" + name;
      ga_undefine_function(name1);
      std::string name2 = "DER_PDFUNC2_" + name;
      ga_undefine_function(name2);
    }
  }

  //=========================================================================
  // Instructions for compilation: basic optimized operations on tensors
  //=========================================================================



  struct ga_instruction_extract_local_im_data : public ga_instruction {
    base_tensor &t;
    const im_data &imd;
    papprox_integration &pai;
    const base_vector &U;
    fem_interpolation_context &ctx;
    size_type qdim, cv_old;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: extract local im data");
      size_type cv = ctx.convex_num();
      if (cv != cv_old) {
        cv_old = cv;
        GMM_ASSERT1(imd.linked_mesh_im().int_method_of_element(cv)
                    ->approx_method() == pai, "Im data have to be used only "
                    "on their original integration method.");
      }
      size_type ipt = imd.filtered_index_of_point(cv, ctx.ii());
      GMM_ASSERT1(ipt != size_type(-1),
                  "Im data with no data on the current integration point.");
      gmm::copy(gmm::sub_vector(U, gmm::sub_interval(ipt*qdim, qdim)),
                t.as_vector());
      return 0;
    }
    ga_instruction_extract_local_im_data
    (base_tensor &t_, const im_data &imd_, const base_vector &U_,
     papprox_integration &pai_, fem_interpolation_context &ctx_,
     size_type qdim_)
      : t(t_), imd(imd_), pai(pai_), U(U_), ctx(ctx_), qdim(qdim_),
        cv_old(-1) {}
  };

  struct ga_instruction_slice_local_dofs : public ga_instruction {
    const mesh_fem &mf;
    const base_vector &U;
    fem_interpolation_context &ctx;
    base_vector &coeff;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Slice local dofs");
      slice_vector_on_basic_dof_of_element(mf, U, ctx.convex_num(), coeff);
      return 0;
    }
    ga_instruction_slice_local_dofs(const mesh_fem &mf_, const base_vector &U_,
                                    fem_interpolation_context &ctx_,
                                    base_vector &coeff_)
      : mf(mf_), U(U_), ctx(ctx_), coeff(coeff_) {}
  };

  struct ga_instruction_update_pfp : public ga_instruction {
    const mesh_fem &mf;
    fem_interpolation_context &ctx;
    fem_precomp_pool &fp_pool;
    pfem_precomp &pfp;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Pfp update");
      if (ctx.have_pgp()) {
        pfem pf = mf.fem_of_element(ctx.convex_num());
        if (!pfp || pf != pfp->get_pfem() ||
            &(ctx.pgp()->get_point_tab()) != &(pfp->get_point_tab())) {
          if (pf->is_on_real_element())
            pfp = 0;
          else {
            pfp = fp_pool(pf, &(ctx.pgp()->get_point_tab()));
          }
        }
      } else {
        pfp = 0;
      }
      return 0;
    }

    ga_instruction_update_pfp(const mesh_fem &mf_, pfem_precomp &pfp_,
                              fem_interpolation_context &ctx_,
                              fem_precomp_pool &fp_pool_)
      : mf(mf_), ctx(ctx_) , fp_pool(fp_pool_), pfp(pfp_) {}
  };

  struct ga_instruction_first_ind_tensor : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    size_type qdim;
    const mesh_fem *mfn, **mfg;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: adapt first index of tensor");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      // cout << "mfg = " << mfg << endl;
      GA_DEBUG_ASSERT(&mf, "Internal error");
      // cout << "ctx.is_convex_num_valid() = " << int(ctx.is_convex_num_valid()) << endl;
      size_type cv_1 = ctx.is_convex_num_valid()
        ? ctx.convex_num() : mf.convex_index().first_true();
      if (!mf.convex_index().is_in(cv_1))
        cv_1 = mf.convex_index().first_true();


      // cout << "cv_1 " << cv_1 << endl;
      pfem pf = mf.fem_of_element(cv_1);
      GMM_ASSERT1(pf, "An element without finite element methode defined");
      // cout << "pf = " << pf << endl;
      size_type Qmult = qdim / pf->target_dim();
      size_type s = pf->nb_dof(cv_1) * Qmult;
      if (t.sizes()[0] != s)
        { bgeot::multi_index mi = t.sizes(); mi[0] = s; t.adjust_sizes(mi); }
      // cout << "here " << endl;
      return 0;
    }

    ga_instruction_first_ind_tensor(base_tensor &t_,
                                    fem_interpolation_context &ctx_,
                                    size_type qdim_, const mesh_fem *mfn_,
                                    const mesh_fem **mfg_)
      : t(t_),  ctx(ctx_), qdim(qdim_), mfn(mfn_), mfg(mfg_) {}
  };

  struct ga_instruction_second_ind_tensor : public ga_instruction_first_ind_tensor {

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: adapt second index of tensor");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      size_type cv_1 = ctx.is_convex_num_valid()
        ? ctx.convex_num() : mf.convex_index().first_true();
      pfem pf = mf.fem_of_element(cv_1);
      GMM_ASSERT1(pf, "An element without finite element methode defined");
      size_type Qmult = qdim / pf->target_dim();
      size_type s = pf->nb_dof(cv_1) * Qmult;
      if (t.sizes()[1] != s)
        { bgeot::multi_index mi = t.sizes(); mi[1] = s; t.adjust_sizes(mi); }
      return 0;
    }

    ga_instruction_second_ind_tensor(base_tensor &t_, fem_interpolation_context &ctx_,
                                    size_type qdim_, const mesh_fem *mfn_,
                                    const mesh_fem **mfg_)
   : ga_instruction_first_ind_tensor(t_, ctx_, qdim_, mfn_, mfg_)
   {}
;

  };

  struct ga_instruction_two_first_ind_tensor : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx1;
    fem_interpolation_context &ctx2;
    size_type qdim1;
    const mesh_fem *mfn1, **mfg1;
    size_type qdim2;
    const mesh_fem *mfn2, **mfg2;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: adapt two first indices of tensor");
      const mesh_fem &mf1 = *(mfg1 ? *mfg1 : mfn1);
      const mesh_fem &mf2 = *(mfg2 ? *mfg2 : mfn2);
      size_type cv_1 = ctx1.is_convex_num_valid()
        ? ctx1.convex_num() : mf1.convex_index().first_true();
      size_type cv_2 = ctx2.is_convex_num_valid()
        ? ctx2.convex_num() : mf2.convex_index().first_true();
      pfem pf1 = mf1.fem_of_element(cv_1);
      GMM_ASSERT1(pf1, "An element without finite element method defined");
      pfem pf2 = mf2.fem_of_element(cv_2);
      GMM_ASSERT1(pf2, "An element without finite element method defined");
      size_type Qmult1 = qdim1 / pf1->target_dim();
      size_type s1 = pf1->nb_dof(cv_1) * Qmult1;
      size_type Qmult2 = qdim2 / pf2->target_dim();
      size_type s2 = pf2->nb_dof(cv_2) * Qmult2;
      if (t.sizes()[0] != s1 || t.sizes()[1] != s2) {
        bgeot::multi_index mi = t.sizes();
        mi[0] = s1; mi[1] = s2;
        t.adjust_sizes(mi);
      }
      return 0;
    }

    ga_instruction_two_first_ind_tensor
    (base_tensor &t_, fem_interpolation_context &ctx1_,
     fem_interpolation_context &ctx2_,
     size_type qdim1_, const mesh_fem *mfn1_, const mesh_fem **mfg1_,
     size_type qdim2_, const mesh_fem *mfn2_, const mesh_fem **mfg2_)
      : t(t_),  ctx1(ctx1_), ctx2(ctx2_), qdim1(qdim1_), mfn1(mfn1_),
        mfg1(mfg1_), qdim2(qdim2_), mfn2(mfn2_), mfg2(mfg2_) {}
  };


  struct ga_instruction_X_component : public ga_instruction {
    scalar_type &t;
    fem_interpolation_context &ctx;
    size_type n;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: X component");
      t = ctx.xreal()[n];
      return 0;
    }

    ga_instruction_X_component(scalar_type &t_,
                               fem_interpolation_context &ctx_, size_type n_)
      : t(t_),  ctx(ctx_), n(n_) {}
  };

  struct ga_instruction_X : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: X");
      GA_DEBUG_ASSERT(t.size() == ctx.xreal().size(), "dimensions mismatch");
      gmm::copy(ctx.xreal(), t.as_vector());
      return 0;
    }

    ga_instruction_X(base_tensor &t_, fem_interpolation_context &ctx_)
      : t(t_),  ctx(ctx_) {}
  };

  struct ga_instruction_Normal : public ga_instruction {
    base_tensor &t;
    base_small_vector &Normal;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Normal");
      GMM_ASSERT1(t.size() == Normal.size(), "Invalid outward unit normal "
                  "vector. Possible reasons: not on boundary or "
                  "transformation failed.");
      gmm::copy(Normal, t.as_vector());
      return 0;
    }
    ga_instruction_Normal(base_tensor &t_, base_small_vector &Normal_)
      : t(t_), Normal(Normal_)  {}
  };

  struct ga_instruction_element_size : public ga_instruction {
    base_tensor &t;
    scalar_type &es;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: element_size");
      GMM_ASSERT1(t.size() == 1, "Invalid element size.");
      t[0] = es;
      return 0;
    }
    ga_instruction_element_size(base_tensor &t_, scalar_type &es_)
      : t(t_), es(es_)  {}
  };

  struct ga_instruction_element_K : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: element_K");
      GMM_ASSERT1(t.size() == (ctx.K()).size(), "Invalid tensor size.");
      gmm::copy(ctx.K().as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_element_K(base_tensor &t_, fem_interpolation_context &ct)
      : t(t_), ctx(ct)  {}
  };

  struct ga_instruction_element_B : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: element_B");
      GMM_ASSERT1(t.size() == (ctx.B()).size(), "Invalid tensor size.");
      gmm::copy(ctx.B().as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_element_B(base_tensor &t_, fem_interpolation_context &ct)
      : t(t_), ctx(ct)  {}
  };

  struct ga_instruction_val_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;

    virtual int exec(void) { // --> t(ndof,target_dim)
      GA_DEBUG_INFO("Instruction: compute value of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      ctx.pf()->real_base_value(ctx, t);
      return 0;
    }

    ga_instruction_val_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_xfem_plus_val_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;

    virtual int exec(void) { // --> t(ndof,target_dim)
      GA_DEBUG_INFO("Instruction: compute value of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(1);
      ctx.pf()->real_base_value(ctx, t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_plus_val_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_xfem_minus_val_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;

    virtual int exec(void) { // --> t(ndof,target_dim)
      GA_DEBUG_INFO("Instruction: compute value of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(-1);
      ctx.pf()->real_base_value(ctx, t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_minus_val_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_grad_base : public ga_instruction_val_base {

    virtual int exec(void) { // --> t(ndof,target_dim,N)
      GA_DEBUG_INFO("Instruction: compute gradient of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      ctx.pf()->real_grad_base_value(ctx, t);
      return 0;
    }

    ga_instruction_grad_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_,pfp_)
    {}
  };

  struct ga_instruction_xfem_plus_grad_base : public ga_instruction_val_base {

    virtual int exec(void) { // --> t(ndof,target_dim,N)
      GA_DEBUG_INFO("Instruction: compute gradient of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(1);
      ctx.pf()->real_grad_base_value(ctx, t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_plus_grad_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_,pfp_)
    {}
  };

  struct ga_instruction_xfem_minus_grad_base : public ga_instruction_val_base {

    virtual int exec(void) { // --> t(ndof,target_dim,N)
      GA_DEBUG_INFO("Instruction: compute gradient of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(-1);
      ctx.pf()->real_grad_base_value(ctx, t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_minus_grad_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_,pfp_)
    {}
  };


  struct ga_instruction_hess_base : public ga_instruction_val_base {

    virtual int exec(void) { // --> t(ndof,target_dim,N,N)
      GA_DEBUG_INFO("Instruction: compute Hessian of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      ctx.pf()->real_hess_base_value(ctx, t);
      return 0;
    }

    ga_instruction_hess_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_xfem_plus_hess_base : public ga_instruction_val_base {

    virtual int exec(void) { // --> t(ndof,target_dim,N,N)
      GA_DEBUG_INFO("Instruction: compute Hessian of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(1);
      ctx.pf()->real_hess_base_value(ctx, t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_plus_hess_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_xfem_minus_hess_base : public ga_instruction_val_base {

    virtual int exec(void) { // --> t(ndof,target_dim,N,N)
      GA_DEBUG_INFO("Instruction: compute Hessian of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      int old_xfem_side = ctx.xfem_side();
      ctx.set_xfem_side(-1);
      ctx.pf()->real_hess_base_value(ctx, t);
      ctx.set_xfem_side(old_xfem_side);
      return 0;
    }

    ga_instruction_xfem_minus_hess_base(base_tensor &tt, fem_interpolation_context &ct,
                            const mesh_fem &mf_, pfem_precomp &pfp_)
    : ga_instruction_val_base(tt, ct, mf_, pfp_)
    {}
  };

  struct ga_instruction_val : public ga_instruction {
    base_tensor &t;
    const base_tensor &Z;
    const base_vector &coeff;
    size_type qdim;
    // Z(ndof,target_dim), coeff(Qmult,ndof) --> t(target_dim*Qmult)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: variable value");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(t.size() == qdim, "dimensions mismatch");
      GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*Qmult,
                      "Wrong size for coeff vector");

      gmm::clear(t.as_vector());
      for (size_type j = 0; j < ndof; ++j) {
        for (size_type q = 0; q < Qmult; ++q) {
          scalar_type co = coeff[j*Qmult+q];
          for (size_type r = 0; r < target_dim; ++r)
            t[r + q*target_dim] += co * Z[j + r*ndof];
        }
      }
      GA_DEBUG_INFO("Instruction: end of variable value");
      return 0;
    }

    ga_instruction_val(base_tensor &tt, const base_tensor &Z_,
                       const base_vector &co, size_type q)
      : t(tt), Z(Z_), coeff(co), qdim(q) {}
  };

  struct ga_instruction_grad : public ga_instruction_val {
    // Z(ndof,target_dim,N), coeff(Qmult,ndof) --> t(target_dim*Qmult,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: gradient");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT((qdim == 1 && t.sizes()[0] == N) ||
                      (t.sizes()[1] == N && t.sizes()[0] == qdim) ||
                      (N == 1 && t.sizes()[0] == qdim),
                      "dimensions mismatch");
      GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*Qmult,
                      "Wrong size for coeff vector");
      gmm::clear(t.as_vector());
      for (size_type q = 0; q < Qmult; ++q) {
        base_tensor::const_iterator it = Z.begin();
        for (size_type k = 0; k < N; ++k)
          for (size_type r = 0; r < target_dim; ++r)
            for (size_type j = 0; j < ndof; ++j, ++it)
              t[r + q*target_dim + k*qdim] += coeff[j*Qmult+q] * (*it);
      }
      return 0;
    }

    ga_instruction_grad(base_tensor &tt, const base_tensor &Z_,
                       const base_vector &co, size_type q)
    : ga_instruction_val(tt, Z_, co, q)
    {}

  };

  struct ga_instruction_hess : public ga_instruction_val {
    // Z(ndof,target_dim,N*N), coeff(Qmult,ndof) --> t(target_dim*Qmult,N,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Hessian");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = size_type(round(sqrt(scalar_type(Z.sizes()[2]))));
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT((qdim == 1 && t.sizes()[0] == N && t.sizes()[1] == N) ||
                      (t.sizes()[1] == N && t.sizes()[2] == N
                       && t.sizes()[0] == qdim), "dimensions mismatch");
      GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*Qmult,
                      "Wrong size for coeff vector");
      gmm::clear(t.as_vector());
      for (size_type q = 0; q < Qmult; ++q) {
        base_tensor::const_iterator it = Z.begin();
        for (size_type k = 0; k < N; ++k)
          for (size_type l = 0; l < N; ++l)
            for (size_type r = 0; r < target_dim; ++r)
              for (size_type j = 0; j < ndof; ++j, ++it)
                t[r + q*target_dim + k*qdim + l*qdim*N]
                  += coeff[j*Qmult+q] * (*it);
      }
      return 0;
    }

    ga_instruction_hess(base_tensor &tt, const base_tensor &Z_,
                       const base_vector &co, size_type q)
    : ga_instruction_val(tt, Z_, co, q)
    {}
  };

  struct ga_instruction_diverg : public ga_instruction_val {
    // Z(ndof,target_dim,N), coeff(Qmult,ndof) --> t(1)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: divergence");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(Qmult*target_dim == N && (Qmult == 1 || target_dim == 1),
                      "Dimensions mismatch for divergence operator");
      GA_DEBUG_ASSERT(gmm::vect_size(coeff) == ndof*Qmult,
                      "Wrong size for coeff vector");

      t[0] = scalar_type(0);
      base_tensor::const_iterator it = Z.begin();
      if (Qmult == 1)
        for (size_type k = 0; k < N; ++k) {
          if (k) it += (N*ndof + 1);
          for (size_type j = 0; j < ndof; ++j) {
            if (j) ++it;
            t[0] += coeff[j] * (*it);
          }
        }
      else // if (target_dim() == 1)
        for (size_type k = 0; k < N; ++k) {
          if (k) ++it;
          for (size_type j = 0; j < ndof; ++j) {
            if (j) ++it;
            t[0] += coeff[j*N+k] * (*it);
          }
        }
      return 0;
    }

    ga_instruction_diverg(base_tensor &tt, const base_tensor &Z_,
                          const base_vector &co, size_type q)
    : ga_instruction_val(tt, Z_, co, q)
    {}
  };


  struct ga_instruction_copy_val_base : public ga_instruction {
    base_tensor &t;
    const base_tensor &Z;
    size_type qdim;
    // Z(ndof,target_dim) --> t(Qmult*ndof,Qmult*target_dim)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: value of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                      "Wrong size for base vector");
      if (Qmult == 1) {
        gmm::copy(Z.as_vector(), t.as_vector());
      } else {
        gmm::clear(t.as_vector());
        base_tensor::const_iterator itZ = Z.begin();
        size_type s = t.sizes()[0], ss = s * Qmult, sss = s+1;

        // Performs t(i*Qmult+j, k*Qmult + j) = Z(i,k);
        for (size_type k = 0; k < target_dim; ++k) {
          base_tensor::iterator it = t.begin() + (ss * k);
          for (size_type i = 0; i < ndof; ++i, ++itZ) {
            if (i) it += Qmult;
            base_tensor::iterator it2 = it;
            *it2 = *itZ;
            for (size_type j = 1; j < Qmult; ++j) { it2 += sss; *it2 = *itZ; }
          }
        }
      }
      return 0;
    }

    ga_instruction_copy_val_base(base_tensor &tt, const base_tensor &Z_, size_type q)
      : t(tt), Z(Z_), qdim(q) {}
  };

  struct ga_instruction_copy_grad_base : public ga_instruction_copy_val_base {
    // Z(ndof,target_dim,N) --> t(Qmult*ndof,Qmult*target_dim,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: gradient of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                      "Wrong size for gradient vector");
      if (Qmult == 1) {
        gmm::copy(Z.as_vector(), t.as_vector());
      } else {
        gmm::clear(t.as_vector());
        base_tensor::const_iterator itZ = Z.begin();
        size_type s = t.sizes()[0], ss = s * Qmult, sss = s+1;
        size_type ssss=ss*target_dim;

        // Performs t(i*Qmult+j, k*Qmult + j, l) = Z(i,k,l);
        for (size_type l = 0; l < N; ++l)
          for (size_type k = 0; k < target_dim; ++k) {
            base_tensor::iterator it = t.begin() + (ss * k + ssss*l);
            for (size_type i = 0; i < ndof; ++i, ++itZ) {
              if (i) it += Qmult;
              base_tensor::iterator it2 = it;
              *it2 = *itZ;
              for (size_type j = 1; j < Qmult; ++j) { it2 += sss; *it2 = *itZ; }
            }
          }
      }
      return 0;
    }

     ga_instruction_copy_grad_base(base_tensor &tt, const base_tensor &Z_, size_type q)
     : ga_instruction_copy_val_base(tt,Z_,q)
     {}
  };

  struct ga_instruction_copy_hess_base : public ga_instruction_copy_val_base {
    // Z(ndof,target_dim,N*N) --> t(Qmult*ndof,Qmult*target_dim,N,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Hessian of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N2 = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(t.size() == Z.size() * Qmult * Qmult,
                      "Wrong size for Hessian vector");
      if (Qmult == 1) {
        gmm::copy(Z.as_vector(), t.as_vector());
      } else {
        gmm::clear(t.as_vector());

        base_tensor::const_iterator itZ = Z.begin();
        size_type s = t.sizes()[0], ss = s * Qmult, sss = s+1;
        size_type ssss=ss*target_dim;

        // Performs t(i*Qmult+j, k*Qmult + j, l, m) = Z(i,k,l*N+m)
        for (size_type l = 0; l < N2; ++l)
          for (size_type k = 0; k < target_dim; ++k) {
            base_tensor::iterator it = t.begin() + (ss * k + ssss*l);
            for (size_type i = 0; i < ndof; ++i, ++itZ) {
              if (i) it += Qmult;
              base_tensor::iterator it2 = it;
              *it2 = *itZ;
              for (size_type j = 1; j < Qmult; ++j) { it2 += sss; *it2 = *itZ; }
            }
          }
      }
      return 0;
    }

    ga_instruction_copy_hess_base(base_tensor &tt, const base_tensor &Z_, size_type q)
    : ga_instruction_copy_val_base(tt, Z_, q)
    {}
  };

  struct ga_instruction_copy_diverg_base : public ga_instruction_copy_val_base {
    // Z(ndof,target_dim,N) --> t(Qmult*ndof)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: divergence of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GA_DEBUG_ASSERT(Qmult*target_dim == N && (Qmult == 1 || target_dim == 1),
                      "Dimensions mismatch for divergence operator");
      GA_DEBUG_ASSERT(t.size() == ndof * Qmult,
                      "Wrong size for divergence vector");
      gmm::clear(t.as_vector());
      base_tensor::const_iterator itZ = Z.begin();
      if (Qmult == 1) { // target_dim == N
        // Performs t(i) = Trace(Z(i,:,:))
        for (size_type l = 0; l < N; ++l) {
          base_tensor::iterator it = t.begin();
          if (l) itZ += target_dim*ndof+1;
          for (size_type i = 0; i < ndof; ++i) {
            if (i) { ++it; ++itZ; }
            *it += *itZ;
          }
        }
      } else { // Qmult == N
        // Performs t(i*Qmult+j) = Z(i,1,j)
        for (size_type j = 0; j < N; ++j) {
          base_tensor::iterator it = t.begin() + j;
          if (j) ++itZ;
          for (size_type i = 0; i < ndof; ++i) {
            if (i) { it += Qmult; ++itZ; }
            *it += *itZ;
          }
        }
      }
      return 0;
    }

    ga_instruction_copy_diverg_base(base_tensor &tt, const base_tensor &Z_, size_type q)
    : ga_instruction_copy_val_base(tt, Z_, q)
    {}
  };


  struct ga_instruction_elementary_transformation {
    const base_vector &coeff_in;
    base_vector coeff_out;
    pelementary_transformation elemtrans;
    const mesh_fem &mf;
    fem_interpolation_context &ctx;
    base_matrix &M;
    const mesh_fem **mf_M;
    size_type &icv;

    void do_transformation(void) {
      size_type nn = gmm::vect_size(coeff_in);
      if (M.size() == 0 || icv != ctx.convex_num() || &mf != *mf_M) {
        gmm::resize(M, nn, nn);
        *mf_M = &mf; icv = ctx.convex_num();
        elemtrans->give_transformation(mf, icv, M);
      }
      coeff_out.resize(nn);
      gmm::mult(M, coeff_in, coeff_out); // remember: coeff == coeff_out
    }

    ga_instruction_elementary_transformation
    (const base_vector &co, pelementary_transformation e,
     const mesh_fem &mf_, fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : coeff_in(co), elemtrans(e), mf(mf_), ctx(ctx_),
        M(M_), mf_M(mf_M_), icv(icv_) {}
    ~ga_instruction_elementary_transformation() {};
  };

  struct ga_instruction_elementary_transformation_val
    : public ga_instruction_val, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim), coeff_in(Qmult,ndof) --> coeff_out --> t(target_dim*Qmult)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: variable value with elementary "
                    "transformation");
      do_transformation();
      return ga_instruction_val::exec();
    }

    ga_instruction_elementary_transformation_val
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_val(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_grad
    : public ga_instruction_grad, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim,N), coeff_in(Qmult,ndof) --> coeff_out --> t(target_dim*Qmult,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: gradient with elementary transformation");
      do_transformation();
      return ga_instruction_grad::exec();
    }

    ga_instruction_elementary_transformation_grad
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_grad(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_hess
    : public ga_instruction_hess, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim,N*N), coeff_in(Qmult,ndof) --> coeff_out --> t(target_dim*Qmult,N,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Hessian with elementary transformation");
      do_transformation();
      return ga_instruction_hess::exec();
    }

    ga_instruction_elementary_transformation_hess
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_hess(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_diverg
    : public ga_instruction_diverg, ga_instruction_elementary_transformation {
    // Z(ndof,target_dim,N), coeff_in(Qmult,ndof) --> coeff_out --> t(1)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: divergence with elementary transformation");
      do_transformation();
      return ga_instruction_diverg::exec();
    }

    ga_instruction_elementary_transformation_diverg
    (base_tensor &tt, const base_tensor &Z_, const base_vector &co, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_diverg(tt, Z_, coeff_out, q),
        ga_instruction_elementary_transformation(co, e, mf_, ctx_, M_,
                                                 mf_M_, icv_) {}
  };

  struct ga_instruction_update_group_info : public ga_instruction {
    ga_workspace &workspace;
    ga_instruction_set &gis;
    ga_instruction_set::interpolate_info &inin;
    const std::string gname;
    ga_instruction_set::variable_group_info &vgi;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Update group info for "+gname);
      if (vgi.varname &&
          &(workspace.associated_mf(*(vgi.varname))->linked_mesh())==inin.m)
        return 0;
      const std::string &varname
        = inin.m ? workspace.variable_in_group(gname, *(inin.m))
                 : workspace.first_variable_of_group(gname);
      vgi.mf = workspace.associated_mf(varname);
      vgi.Ir = gis.var_intervals[varname];
      vgi.In = workspace.interval_of_variable(varname);
      vgi.alpha = workspace.factor_of_variable(varname);
      vgi.U = gis.extended_vars[varname];
      vgi.varname = &varname;
      return 0;
    }

    ga_instruction_update_group_info
    (ga_workspace &workspace_, ga_instruction_set &gis_,
     ga_instruction_set::interpolate_info &inin_, const std::string &gname_,
     ga_instruction_set::variable_group_info &vgi_) :
      workspace(workspace_), gis(gis_), inin(inin_), gname(gname_),
      vgi(vgi_) {}
  };

  struct ga_instruction_interpolate_filter : public ga_instruction {
    base_tensor &t;
    ga_instruction_set::interpolate_info &inin;
    size_type pt_type;
    int nb;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated filter");
      if ((pt_type == size_type(-1) && inin.pt_type) ||
          (pt_type != size_type(-1) && inin.pt_type == pt_type)) {
        GA_DEBUG_INFO("Instruction: interpolated filter: pass");
        return 0;
      }
      else {
        GA_DEBUG_INFO("Instruction: interpolated filter: filtered");
        gmm::clear(t.as_vector());
        return nb;
      }
      return 0;
    }

    ga_instruction_interpolate_filter
    (base_tensor &t_, ga_instruction_set::interpolate_info &inin_,
     size_type ind_, int nb_)
      : t(t_), inin(inin_), pt_type(ind_), nb(nb_) {}
  };


  struct ga_instruction_interpolate : public ga_instruction {
    base_tensor &t;
    const mesh **m;
    const mesh_fem *mfn, **mfg;
    const base_vector *Un, **Ug;
    fem_interpolation_context &ctx;
    base_vector coeff;
    size_type qdim;

    virtual int exec(void) {
      GMM_ASSERT1(ctx.is_convex_num_valid(), "No valid element for the "
                  "transformation. Probably transformation failed");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      const base_vector &U = *(Ug ? *Ug : Un);
      GMM_ASSERT1(&(mf.linked_mesh()) == *m, "Interpolation of a variable "
        "on another mesh than the one it is defined on");
      slice_vector_on_basic_dof_of_element(mf, U, ctx.convex_num(), coeff);
      ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      return 0;
    }

    ga_instruction_interpolate
    (base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_,
     const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q)
      : t(tt), m(m_), mfn(mfn_), mfg(mfg_), Un(Un_), Ug(Ug_),
        ctx(ctx_), qdim(q) {}
  };

  struct ga_instruction_interpolate_val : public ga_instruction_interpolate {
    // --> t(target_dim*Qmult)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated variable value");
      ga_instruction_interpolate::exec();
      ctx.pf()->interpolation(ctx, coeff, t.as_vector(), dim_type(qdim));
      // cout << "interpolate " << &U << " result : " << t.as_vector() << endl;
      return 0;
    }

    ga_instruction_interpolate_val(base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_,
     const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_,ctx_, q)
    {}
  };

  struct ga_instruction_interpolate_grad : public ga_instruction_interpolate {
    // --> t(target_dim*Qmult,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated variable grad");
      ga_instruction_interpolate::exec();
      base_matrix v(qdim, ctx.N());
      ctx.pf()->interpolation_grad(ctx, coeff, v, dim_type(qdim));
      gmm::copy(v.as_vector(), t.as_vector());
      return 0;
    }

    ga_instruction_interpolate_grad(base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_,
     const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_, ctx_, q)
    {}
  };

  struct ga_instruction_interpolate_hess : public ga_instruction_interpolate {
    // --> t(target_dim*Qmult,N,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated variable hessian");
      ga_instruction_interpolate::exec();
      base_matrix v(qdim, ctx.N()*ctx.N());
      ctx.pf()->interpolation_hess(ctx, coeff, v, dim_type(qdim));
      gmm::copy(v.as_vector(), t.as_vector());
      return 0;
    }

    ga_instruction_interpolate_hess(base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_,
     const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_, ctx_, q)
    {}
  };

  struct ga_instruction_interpolate_diverg : public ga_instruction_interpolate {
    // --> t(1)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated variable divergence");
      ga_instruction_interpolate::exec();
      ctx.pf()->interpolation_diverg(ctx, coeff, t[0]);
      return 0;
    }

    ga_instruction_interpolate_diverg(base_tensor &tt, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_,
     const base_vector *Un_, const base_vector **Ug_,
     fem_interpolation_context &ctx_, size_type q)
    : ga_instruction_interpolate(tt, m_, mfn_, mfg_, Un_, Ug_, ctx_, q)
    {}
  };

  struct ga_instruction_interpolate_base {
    base_tensor ZZ;
    const mesh **m;
    const mesh_fem *mfn, **mfg;
    fem_interpolation_context &ctx;

    virtual int exec(void) {
      GMM_ASSERT1(ctx.is_convex_num_valid(), "No valid element for the "
                  "transformation. Probably transformation failed");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      GMM_ASSERT1(&(mf.linked_mesh()) == *m, "Interpolation of a variable "
        "on another mesh than the one it is defined on");
      ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      GMM_ASSERT1(ctx.pf(), "Undefined finite element method");
      return 0;
    }

    ga_instruction_interpolate_base
    (const mesh **m_, const mesh_fem *mfn_, const mesh_fem **mfg_,
     fem_interpolation_context &ctx_)
      : m(m_), mfn(mfn_), mfg(mfg_), ctx(ctx_) {}
  };

  struct ga_instruction_interpolate_val_base
    : public ga_instruction_copy_val_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim) --> t(Qmult*ndof,Qmult*target_dim)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated base value");
      ga_instruction_interpolate_base::exec();
      ctx.pf()->real_base_value(ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_val_base::exec();
    }

    ga_instruction_interpolate_val_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, fem_interpolation_context &ctx_, size_type q)
      : ga_instruction_copy_val_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ctx_) {}
  };

  struct ga_instruction_interpolate_grad_base
    : public ga_instruction_copy_grad_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim,N) --> t(Qmult*ndof,Qmult*target_dim,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated base grad");
      ga_instruction_interpolate_base::exec();
      ctx.pf()->real_grad_base_value(ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_grad_base::exec();
    }

    ga_instruction_interpolate_grad_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, fem_interpolation_context &ctx_, size_type q)
      : ga_instruction_copy_grad_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ctx_) {}
  };

  struct ga_instruction_interpolate_hess_base
    : public ga_instruction_copy_hess_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim,N*N) --> t(Qmult*ndof,Qmult*target_dim,N,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated base hessian");
      ga_instruction_interpolate_base::exec();
      ctx.pf()->real_hess_base_value(ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_hess_base::exec();
    }

    ga_instruction_interpolate_hess_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, fem_interpolation_context &ctx_, size_type q)
      : ga_instruction_copy_hess_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ctx_) {}
  };

  struct ga_instruction_interpolate_diverg_base
    : public ga_instruction_copy_diverg_base, ga_instruction_interpolate_base {
    // ctx --> Z(ndof,target_dim,N*N) --> t(Qmult*ndof)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: interpolated base divergence");
      ga_instruction_interpolate_base::exec();
      ctx.pf()->real_grad_base_value(ctx, ZZ); // remember Z == ZZ
      return ga_instruction_copy_diverg_base::exec();
    }

    ga_instruction_interpolate_diverg_base
    (base_tensor &t_, const mesh **m_, const mesh_fem *mfn_,
     const mesh_fem **mfg_, fem_interpolation_context &ctx_, size_type q)
      : ga_instruction_copy_diverg_base(t_, ZZ, q),
        ga_instruction_interpolate_base(m_, mfn_, mfg_, ctx_) {}
  };


  struct ga_instruction_elementary_transformation_base {
    base_tensor t_in;
    base_tensor &t_out;
    pelementary_transformation elemtrans;
    const mesh_fem &mf;
    fem_interpolation_context &ctx;
    base_matrix &M;
    const mesh_fem **mf_M;
    size_type &icv;

    void do_transformation(size_type n) {
      if (M.size() == 0 || icv != ctx.convex_num() || &mf != *mf_M) {
        gmm::resize(M, n, n);
        *mf_M = &mf; icv = ctx.convex_num();
        elemtrans->give_transformation(mf, icv, M);
      }
      t_out.mat_reduction(t_in, M, 0);
    }

    ga_instruction_elementary_transformation_base
    (base_tensor &t_, pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : t_out(t_), elemtrans(e), mf(mf_), ctx(ctx_),
        M(M_), mf_M(mf_M_), icv(icv_) {}
  };

  struct ga_instruction_elementary_transformation_val_base
    : public ga_instruction_copy_val_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim) --> t_in --> t_out(Qmult*ndof,Qmult*target_dim)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: value of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_val_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_val_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_val_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_grad_base
    : public ga_instruction_copy_grad_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim,N) --> t_in --> t_out(Qmult*ndof,Qmult*target_dim,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: gradient of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_grad_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_grad_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_grad_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_hess_base
    : public ga_instruction_copy_hess_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim,N*N) --> t_in --> t_out(Qmult*ndof,Qmult*target_dim,N,N)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Hessian of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_hess_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_hess_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_hess_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
  };

  struct ga_instruction_elementary_transformation_diverg_base
    : public ga_instruction_copy_diverg_base,
             ga_instruction_elementary_transformation_base {
    // Z(ndof,target_dim,N) --> t_in --> t_out(Qmult*ndof)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: divergence of test functions with elementary "
                    "transformation");
      size_type ndof = Z.sizes()[0];
      size_type Qmult = qdim / Z.sizes()[1];
      t_in.adjust_sizes(t_out.sizes());
      ga_instruction_copy_diverg_base::exec();
      do_transformation(ndof*Qmult);
      return 0;
    }

    ga_instruction_elementary_transformation_diverg_base
    (base_tensor &t_, const base_tensor &Z_, size_type q,
     pelementary_transformation e, const mesh_fem &mf_,
     fem_interpolation_context &ctx_, base_matrix &M_,
     const mesh_fem **mf_M_, size_type &icv_)
      : ga_instruction_copy_diverg_base(t_in, Z_, q),
        ga_instruction_elementary_transformation_base(t_, e, mf_, ctx_, M_,
                                                      mf_M_, icv_) {}
  };


  struct ga_instruction_add : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: addition");
      GA_DEBUG_ASSERT(t.size() == tc1.size(),
                      "internal error " << t.size() << " != " << tc1.size());
      GA_DEBUG_ASSERT(t.size() == tc2.size(),
                      "internal error " << t.size() << " != " << tc2.size());
      gmm::add(tc1.as_vector(), tc2.as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_add(base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_add_to : public ga_instruction {
    base_tensor &t, &tc1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: addition");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "internal error " << t.size()
                      << " incompatible with " << tc1.size());
      gmm::add(tc1.as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_add_to(base_tensor &t_, base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_sub : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: subtraction");
      GA_DEBUG_ASSERT(t.size() == tc1.size() && t.size() == tc2.size(),
                      "internal error");
      gmm::add(tc1.as_vector(), gmm::scaled(tc2.as_vector(), scalar_type(-1)),
               t.as_vector());
      return 0;
    }
    ga_instruction_sub(base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_opposite : public ga_instruction {
    base_tensor &t;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: multiplication with -1");
      gmm::scale(t.as_vector(), scalar_type(-1));
      return 0;
    }
    ga_instruction_opposite(base_tensor &t_) : t(t_) {}
  };

  struct ga_instruction_print_tensor : public ga_instruction {
    base_tensor &t; pga_tree_node pnode;
    fem_interpolation_context &ctx;
    size_type &nbpt, &ipt;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: tensor print");
      cout << "Print term "; ga_print_node(pnode, cout);
      cout << " on Gauss point " << ipt << "/" << nbpt << " of element "
           << ctx.convex_num() << ": " << t << endl;
      return 0;
    }
    ga_instruction_print_tensor(base_tensor &t_, pga_tree_node pnode_,
                                fem_interpolation_context &ctx_,
                                size_type &nbpt_, size_type &ipt_)
      : t(t_), pnode(pnode_), ctx(ctx_), nbpt(nbpt_), ipt(ipt_) {}
  };

  struct ga_instruction_copy_tensor : public ga_instruction {
    base_tensor &t, &tc1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: tensor copy");
      gmm::copy(tc1.as_vector(), t.as_vector());
      return 0;
    }
    ga_instruction_copy_tensor(base_tensor &t_, base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_copy_tensor_possibly_void : public ga_instruction {
    base_tensor &t, &tc1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: tensor copy possibly void");
      if (tc1.size())
        gmm::copy(tc1.as_vector(), t.as_vector());
      else
        gmm::clear(t.as_vector());
      return 0;
    }
    ga_instruction_copy_tensor_possibly_void(base_tensor &t_,
                                             base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_copy_scalar : public ga_instruction {
    scalar_type &t; const scalar_type &t1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: scalar copy");
      t = t1;
      return 0;
    }
    ga_instruction_copy_scalar(scalar_type &t_, const scalar_type &t1_)
      : t(t_), t1(t1_) {}
  };

  struct ga_instruction_copy_vect : public ga_instruction {
    base_vector &t; const base_vector &t1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: fixed size tensor copy");
      gmm::copy(t1, t);
      return 0;
    }
    ga_instruction_copy_vect(base_vector &t_, const base_vector &t1_)
      : t(t_), t1(t1_) {}
  };

  struct ga_instruction_trace : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    size_type n;
    // tc1(:,:,...,n,n) --> t(:,:,...)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Trace");
      GA_DEBUG_ASSERT(t.size()*n*n == tc1.size(), "Wrong sizes");
      size_type s = t.size() * (n+1);
      base_tensor::iterator it = t.begin();
      base_tensor::const_iterator it1 = tc1.begin();
      for (; it != t.end(); ++it, ++it1) {
        *it = scalar_type(0);
        base_tensor::const_iterator it2 = it1;
        *it += *it2;
        for (size_type i = 1; i < n; ++i) { it2 += s; *it += *it2; }
      }
      return 0;
    }

    ga_instruction_trace(base_tensor &t_, base_tensor &tc1_, size_type n_)
      : t(t_), tc1(tc1_), n(n_) {}
  };

  struct ga_instruction_deviator : public ga_instruction {
    base_tensor &t;
    const base_tensor &tc1;
    size_type n;
    // tc1(:,:,...,n,n) --> t(:,:,...,n,n)
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: Deviator");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");

      gmm::copy(tc1.as_vector(), t.as_vector());

      size_type nb = t.size()/(n*n);
      size_type s = nb * (n+1), j = 0;
      base_tensor::iterator it = t.begin();
      base_tensor::const_iterator it1 = tc1.begin();
      for (; j < nb; ++it, ++it1, ++j) {
        scalar_type tr(0);
        base_tensor::const_iterator it2 = it1;
        tr += *it2;
        for (size_type i = 1; i < n; ++i) { it2 += s; tr += *it2; }
        tr /= scalar_type(n);

        base_tensor::iterator it3 = it;
        *it3 -= tr;
        for (size_type i = 1; i < n; ++i) { it3 += s; *it3 -= tr; }
      }
      return 0;
    }

    ga_instruction_deviator(base_tensor &t_, base_tensor &tc1_, size_type n_)
      : t(t_), tc1(tc1_), n(n_) {}
  };

  struct ga_instruction_transpose : public ga_instruction {
    base_tensor &t, &tc1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: transpose");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      size_type order = t.sizes().size();
      size_type s1 = t.sizes()[order-2], s2 = t.sizes()[order-1];
      size_type s = t.size() / (s1*s2);
      for (size_type i = 0; i < s1;  ++i)
        for (size_type j = 0; j < s2;  ++j) {
          base_tensor::iterator it = t.begin() + s*(i + s1*j);
          base_tensor::iterator it1 = tc1.begin() + s*(j + s2*i);
          for (size_type k = 0; k < s; ++k) *it++ = *it1++;
        }
      return 0;
    }
    ga_instruction_transpose(base_tensor &t_, base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_transpose_test : public ga_instruction {
    base_tensor &t, &tc1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: copy tensor and transpose test functions");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      GA_DEBUG_ASSERT(t.sizes().size() >= 2, "Wrong sizes");

      size_type s1 = t.sizes()[0], s2 = t.sizes()[1], s3 = s1*s2;
      size_type s = t.size() / s3;
      base_tensor::iterator it = t.begin();
      for (size_type k = 0; k < s; ++k)
        for (size_type j = 0; j < s2;  ++j)
          for (size_type i = 0; i < s1; ++i, ++it)
            *it = tc1[j+s2*i+k*s3];
      return 0;
    }
    ga_instruction_transpose_test(base_tensor &t_, base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };


  struct ga_instruction_scalar_add : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: scalar addition");
      t = c + d;
      return 0;
    }
    ga_instruction_scalar_add(scalar_type &t_, const scalar_type &c_,
                              const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_sub : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: scalar subtraction");
      t = c - d;
      return 0;
    }
    ga_instruction_scalar_sub(scalar_type &t_, const scalar_type &c_,
                              const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_scalar_mult : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: scalar multiplication");
      t = c * d;
      return 0;
    }
    ga_instruction_scalar_scalar_mult(scalar_type &t_, const scalar_type &c_,
                                      const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_scalar_div : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: scalar division");
      t = c / d;
      return 0;
    }
    ga_instruction_scalar_scalar_div(scalar_type &t_, const scalar_type &c_,
                                     const  scalar_type &d_)
      : t(t_), c(c_), d(d_) {}
  };

  struct ga_instruction_scalar_mult : public ga_instruction {
    base_tensor &t, &tc1;
    const scalar_type &c;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: multiplication of a tensor by a scalar");
      gmm::copy(gmm::scaled(tc1.as_vector(), c), t.as_vector());
      return 0;
    }
    ga_instruction_scalar_mult(base_tensor &t_, base_tensor &tc1_,
                               const scalar_type &c_)
      : t(t_), tc1(tc1_), c(c_) {}
  };

  struct ga_instruction_scalar_div : public ga_instruction {
    base_tensor &t, &tc1;
    const scalar_type &c;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: division of a tensor by a scalar");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");

      base_tensor::iterator it = t.begin(), it1 = tc1.begin();
      for (; it != t.end(); ++it, ++it1) *it = *it1/c;
      return 0;
    }
    ga_instruction_scalar_div(base_tensor &t_, base_tensor &tc1_,
                               const scalar_type &c_)
      : t(t_), tc1(tc1_), c(c_) {}
  };

  struct ga_instruction_dotmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: componentwise multiplication");
      size_type s2 = tc2.size(), s1_1 = tc1.size() / s2;
      GA_DEBUG_ASSERT(t.size() == s1_1*s2, "Wrong sizes");

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s2; ++i)
        for (size_type m = 0; m < s1_1; ++m, ++it)
          *it = tc1[m+s1_1*i] * tc2[i];
      return 0;
    }
    ga_instruction_dotmult(base_tensor &t_, base_tensor &tc1_,
                           base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_dotdiv : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: componentwise division");
      size_type s2 = tc2.size(), s1_1 = tc1.size() / s2;
      GA_DEBUG_ASSERT(t.size() == s1_1*s2, "Wrong sizes");

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s2; ++i)
        for (size_type m = 0; m < s1_1; ++m, ++it)
          *it = tc1[m+s1_1*i] / tc2[i];
      return 0;
    }
    ga_instruction_dotdiv(base_tensor &t_, base_tensor &tc1_,
                          base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ami Bni -> Cmni
  struct ga_instruction_dotmult_spec : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: specific componentwise "
                           "multiplication");
      size_type s2_1 = tc2.sizes()[0], s2_2 = tc2.size() / s2_1;
      size_type s1_1 = tc1.size() / s2_2;

      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s2_2; ++i)
        for (size_type n = 0; n < s2_1; ++n)
          for (size_type m = 0; m < s1_1; ++m, ++it)
            *it = tc1[m+s1_1*i] * tc2[n+s2_1*i];
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_dotmult_spec(base_tensor &t_, base_tensor &tc1_,
                           base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Amij Bjk -> Cmik
  struct ga_instruction_matrix_mult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: matrix multiplication");
      size_type order = tc2.sizes().size();
      size_type s2_1 = tc2.sizes()[order-2];
      size_type s2_2 = tc2.sizes()[order-1];
      size_type s1 = tc1.size() / s2_1;
      size_type s2 = tc2.size() / (s2_1*s2_2);

      base_tensor::iterator it = t.begin();
      for (size_type k = 0; k < s2_2; ++k)
        for (size_type i = 0; i < s1; ++i)
          for (size_type m = 0; m < s2; ++m, ++it) {
            *it = scalar_type(0);
            for (size_type j = 0; j < s2_1; ++j)
              *it += tc1[i+j*s1] * tc2[m+j*s2+k*s2_1*s2];
          }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_matrix_mult(base_tensor &t_, base_tensor &tc1_,
                               base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Amij Bnjk -> Cmnik
  struct ga_instruction_matrix_mult_spec : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: specific matrix multiplication");
      size_type s1_1 = tc1.sizes()[0];
      size_type s1_2 = tc1.sizes()[1];
      // size_type s1_3 = tc1.sizes()[2];
      size_type s2_1 = tc2.sizes()[0];
      size_type s2_2 = tc2.sizes()[1];
      size_type s2_3 = tc2.sizes()[2];

      base_tensor::iterator it = t.begin();
      for (size_type k = 0; k < s2_3; ++k)
        for (size_type i = 0; i < s1_2; ++i)
          for (size_type n = 0; n < s2_1; ++n)
            for (size_type m = 0; m < s1_1; ++m, ++it) {
              *it = scalar_type(0);
              for (size_type j = 0; j < s2_2; ++j)
                *it += tc1[m+i*s1_1+j*s1_1*s1_2] * tc2[n+j*s2_1+k*s2_1*s2_2];
            }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_matrix_mult_spec(base_tensor &t_, base_tensor &tc1_,
                               base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };


  // Performs Ani Bmi -> Cmn for i index of size 2  Unroll loop test.
  struct ga_instruction_reduction_2 : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: reduction operation of size 2 optimized ");
      size_type s1 = tc1.size()/2, s2 = tc2.size()/2;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      base_tensor::iterator it1=tc1.begin(), it2=tc2.begin(), it2end=it2 + s2;
      for (base_tensor::iterator it = t.begin(); it != t.end(); ++it) {
        *it = (*it1)*(*it2) + it1[s1] * it2[s2];
        ++it2; if (it2 == it2end) { it2 = tc2.begin(), ++it1; }
      }
      return 0;
    }
    ga_instruction_reduction_2(base_tensor &t_, base_tensor &tc1_,
                               base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ani Bmi -> Cmn
  struct ga_instruction_reduction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: reduction operation of size " << nn);
      #ifdef GA_USES_BLAS
      int m = int(tc1.size()/nn), k = int(nn), n = int(tc2.size()/nn);
      int lda = m, ldb = n, ldc = m;
      char T = 'T', N = 'N';
      scalar_type alpha(1), beta(0);
      gmm::dgemm_(&N, &T, &m, &n, &k, &alpha, &(tc1[0]), &lda, &(tc2[0]), &ldb,
                  &beta, &(t[0]), &ldc);
      #else
      size_type s1 = tc1.size()/nn, s2 = tc2.size()/nn;
      GA_DEBUG_ASSERT(t.size() == s1*s2, "Internal error");

      base_tensor::iterator it1=tc1.begin(), it2=tc2.begin(), it2end=it2 + s2;
      for (base_tensor::iterator it = t.begin(); it != t.end(); ++it) {
        base_tensor::iterator it11 = it1, it22 = it2;
        scalar_type a = (*it11) * (*it22);
        for (size_type i = 1; i < nn; ++i)
          { it11 += s1; it22 += s2; a += (*it11) * (*it22); }
        *it = a;
        ++it2; if (it2 == it2end) { it2 = tc2.begin(), ++it1; }
      }
      #endif
      return 0;
    }
    ga_instruction_reduction(base_tensor &t_, base_tensor &tc1_,
                             base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Amij Bnj -> Cmni
  struct ga_instruction_spec_reduction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: specific reduction operation of "
                    "size " << nn);
      size_type s1 = tc1.sizes()[0], s11 = tc1.size() / (s1*nn), s111 = s1*s11;
      size_type s2 = tc2.sizes()[0];
      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < s11; ++i)
        for (size_type n = 0; n < s2; ++n)
          for (size_type m = 0; m < s1; ++m, ++it) {
            *it = scalar_type(0);
            for (size_type j = 0; j < nn; ++j)
              *it += tc1[m+i*s1+j*s111] * tc2[n+j*s2];
          }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec_reduction(base_tensor &t_, base_tensor &tc1_,
                                  base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Amik Bnjk -> Cmnij
  struct ga_instruction_spec2_reduction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type nn;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: second specific reduction operation of "
                    "size " << nn);
      size_type s1 = tc1.sizes()[0], s11 = tc1.size() / (s1*nn), s111 = s1*s11;
      size_type s2 = tc2.sizes()[0], s22 = tc2.size() / (s2*nn), s222 = s2*s22;
      base_tensor::iterator it = t.begin();
      for (size_type j = 0; j < s22; ++j)
        for (size_type i = 0; i < s11; ++i)
          for (size_type m = 0; m < s1; ++m)
            for (size_type n = 0; n < s2; ++n, ++it) {
              *it = scalar_type(0);
              for (size_type k = 0; k < nn; ++k)
                *it += tc1[m+i*s1+k*s111] * tc2[n+j*s2+k*s222];
            }
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec2_reduction(base_tensor &t_, base_tensor &tc1_,
                                   base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), nn(n_) {}
  };

  // Performs Aij Bjk -> Cijkl
  struct ga_instruction_simple_tmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: simple tensor product");
      size_type s1 = tc1.size();
      GA_DEBUG_ASSERT(t.size() == s1 * tc2.size(), "Wrong sizes");
      base_tensor::iterator it2=tc2.begin(), it1=tc1.begin(), it1end=it1 + s1;
      for (base_tensor::iterator it = t.begin(); it != t.end(); ++it) {
        *it = *(it2) * (*it1);
        ++it1; if (it1 == it1end) { it1 = tc1.begin(), ++it2; }
      }
      return 0;
    }
    ga_instruction_simple_tmult(base_tensor &t_, base_tensor &tc1_,
                                base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  // Performs Ami Bnj -> Cmnij
  struct ga_instruction_spec_tmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type s1_2, s2_2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: specific tensor product");
      GA_DEBUG_ASSERT(t.size() == tc1.size() * tc2.size(), "Wrong sizes");
      size_type s1_1 = tc1.size() / s1_2;
      size_type s2_1 = tc2.size() / s2_2;

      base_tensor::iterator it = t.begin();
      for (size_type j = 0; j < s2_2; ++j)
        for (size_type i = 0; i < s1_2; ++i)
          for (size_type n = 0; n < s2_1; ++n)
            for (size_type m = 0; m < s1_1; ++m, ++it)
              *it = tc1[m+i*s1_1] * tc2[n+j*s2_1];
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec_tmult(base_tensor &t_, base_tensor &tc1_,
                              base_tensor &tc2_, size_type s1_2_,
                              size_type s2_2_)
      : t(t_), tc1(tc1_), tc2(tc2_), s1_2(s1_2_), s2_2(s2_2_) {}
  };

  // Performs Ai Bmj -> Cmij
  struct ga_instruction_spec2_tmult : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: second specific tensor product");
      GA_DEBUG_ASSERT(t.size() == tc1.size() * tc2.size(), "Wrong sizes");
      size_type s1 = tc1.size();
      size_type s2_1 = tc2.sizes()[0], s2_2 = tc2.size() / s2_1;

      base_tensor::iterator it = t.begin();
      for (size_type j = 0; j < s2_2; ++j)
        for (size_type i = 0; i < s1; ++i)
          for (size_type m = 0; m < s2_1; ++m, ++it)
            *it = tc1[i] * tc2[m+j*s2_1];
      GA_DEBUG_ASSERT(it == t.end(), "Wrong sizes");
      return 0;
    }
    ga_instruction_spec2_tmult(base_tensor &t_, base_tensor &tc1_,
                              base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };



  struct ga_instruction_simple_c_matrix : public ga_instruction {
    base_tensor &t;
    std::vector<scalar_type *> components;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: gathering components for explicit "
                    "matrix");
      GA_DEBUG_ASSERT(t.size() == components.size(), "Wrong sizes");
      for (size_type i = 0; i < components.size(); ++i)
        t[i] = *(components[i]);
      return 0;
    }
    ga_instruction_simple_c_matrix(base_tensor &t_,
                                   std::vector<scalar_type *> &components_)
      : t(t_), components(components_) {}
  };

  struct ga_instruction_c_matrix_with_tests : public ga_instruction {
    base_tensor &t;
    std::vector<base_tensor *> components;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: gathering components for explicit "
                    "matrix with tests functions");
      size_type s = t.size() / components.size();
      GA_DEBUG_ASSERT(s, "Wrong sizes");
      base_tensor::iterator it = t.begin();
      for (size_type i = 0; i < components.size(); ++i) {
        base_tensor &t1 = *(components[i]);
        if (t1.size() > 1) {
          GA_DEBUG_ASSERT(t1.size() == s, "Wrong sizes, " << t1.size() << " != " << s);
          for (size_type j = 0; j < s; ++j) *it++ = t1[j];
        } else {
          for (size_type j = 0; j < s; ++j) *it++ = t1[0];
        }
      }
      return 0;
    }
    ga_instruction_c_matrix_with_tests(base_tensor &t_,
                                       std::vector<base_tensor *>  &components_)
      : t(t_), components(components_) {}
  };

  struct ga_instruction_eval_func_1arg_1res : public ga_instruction {
    scalar_type &t;
    const scalar_type &c;
    pscalar_func_onearg f1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on a scalar");
      t = (*f1)(c);
      return 0;
    }
    ga_instruction_eval_func_1arg_1res(scalar_type &t_, const scalar_type &c_,
                                       pscalar_func_onearg f1_)
      : t(t_), c(c_), f1(f1_) {}
  };

  struct ga_instruction_eval_func_1arg_1res_expr : public ga_instruction {
    scalar_type &t;
    const scalar_type &c;
    const ga_predef_function &F;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on a scalar");
      t = F(c);
      return 0;
    }
    ga_instruction_eval_func_1arg_1res_expr(scalar_type &t_,
                                            const scalar_type &c_,
                                            const ga_predef_function &F_)
      : t(t_), c(c_), F(F_) {}
  };

  struct ga_instruction_eval_func_1arg : public ga_instruction {
    base_tensor &t, &tc1;
    pscalar_func_onearg f1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on tensor");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f1)(tc1[i]);
      return 0;
    }
    ga_instruction_eval_func_1arg(base_tensor &t_, base_tensor &c_,
                                  pscalar_func_onearg f1_)
      : t(t_), tc1(c_), f1(f1_) {}
  };

  struct ga_instruction_eval_func_1arg_expr : public ga_instruction {
    base_tensor &t, &tc1;
    const ga_predef_function &F;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a one argument "
                    "predefined function on tensor");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[i]);
      return 0;
    }
    ga_instruction_eval_func_1arg_expr(base_tensor &t_, base_tensor &c_,
                                       const ga_predef_function &F_)
      : t(t_), tc1(c_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg_1res : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    pscalar_func_twoargs f2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on two scalar");
      t = (*f2)(c, d);
      return 0;
    }
    ga_instruction_eval_func_2arg_1res(scalar_type &t_, const scalar_type &c_,
                                       const scalar_type &d_,
                                       pscalar_func_twoargs f2_)
      : t(t_), c(c_), d(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_1res_expr : public ga_instruction {
    scalar_type &t;
    const scalar_type &c, &d;
    const ga_predef_function &F;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                           "predefined function on two scalar");
      t = F(c, d);
      return 0;
    }
    ga_instruction_eval_func_2arg_1res_expr(scalar_type &t_,
                                            const scalar_type &c_,
                                            const scalar_type &d_,
                                            const ga_predef_function &F_)
      : t(t_), c(c_), d(d_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg_first_scalar : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    pscalar_func_twoargs f2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one scalar and one tensor");
      GA_DEBUG_ASSERT(t.size() == tc2.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f2)(tc1[0], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg_first_scalar
    (base_tensor &t_, base_tensor &c_, base_tensor &d_,
     pscalar_func_twoargs f2_)
      : t(t_), tc1(c_), tc2(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_first_scalar_expr
    : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    const ga_predef_function &F;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one scalar and one tensor");
      GA_DEBUG_ASSERT(t.size() == tc2.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[0], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg_first_scalar_expr
    (base_tensor &t_, base_tensor &c_, base_tensor &d_,
     const ga_predef_function &F_)
      : t(t_), tc1(c_), tc2(d_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg_second_scalar : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    pscalar_func_twoargs f2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one tensor and one scalar");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f2)(tc1[i], tc2[0]);
      return 0;
    }
    ga_instruction_eval_func_2arg_second_scalar(base_tensor &t_,
                                                base_tensor &c_,
                                                base_tensor &d_,
                                                pscalar_func_twoargs f2_)
      : t(t_), tc1(c_), tc2(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_second_scalar_expr
    : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    const ga_predef_function &F;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on one tensor and one scalar");
      GA_DEBUG_ASSERT(t.size() == tc1.size(), "Wrong sizes");
      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[i], tc2[0]);
      return 0;
    }
    ga_instruction_eval_func_2arg_second_scalar_expr
    (base_tensor &t_, base_tensor &c_, base_tensor &d_,
     const ga_predef_function &F_)
      : t(t_), tc1(c_), tc2(d_), F(F_) {}
  };

  struct ga_instruction_eval_func_2arg : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    pscalar_func_twoargs f2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on two tensors");
      GA_DEBUG_ASSERT(t.size() == tc1.size() && t.size() == tc2.size(),
                      "Wrong sizes");

      for (size_type i = 0; i < t.size(); ++i)  t[i] = (*f2)(tc1[i], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg(base_tensor &t_, base_tensor &c_,
                                  base_tensor &d_, pscalar_func_twoargs f2_)
      : t(t_), tc1(c_), tc2(d_), f2(f2_) {}
  };

  struct ga_instruction_eval_func_2arg_expr : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    const ga_predef_function &F;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: evaluation of a two arguments "
                    "predefined function on two tensors");
      GA_DEBUG_ASSERT(t.size() == tc1.size() && t.size() == tc2.size(),
                      "Wrong sizes");

      for (size_type i = 0; i < t.size(); ++i)  t[i] = F(tc1[i], tc2[i]);
      return 0;
    }
    ga_instruction_eval_func_2arg_expr(base_tensor &t_, base_tensor &c_,
                                       base_tensor &d_,
                                       const ga_predef_function &F_)
      : t(t_), tc1(c_), tc2(d_), F(F_) {}
  };

  struct ga_instruction_eval_OP : public ga_instruction {
    base_tensor &t;
    const ga_nonlinear_operator &OP;
    ga_nonlinear_operator::arg_list args;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: operator evaluation");
      OP.value(args, t);
      return 0;
    }
    ga_instruction_eval_OP(base_tensor &t_, const ga_nonlinear_operator &OP_,
                           ga_nonlinear_operator::arg_list &args_)
      : t(t_), OP(OP_), args(args_) {}
  };

  struct ga_instruction_eval_derivative_OP : public ga_instruction {
    base_tensor &t;
    const ga_nonlinear_operator &OP;
    ga_nonlinear_operator::arg_list args;
    size_type der1;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: operator derivative evaluation");
      OP.derivative(args, der1, t);
      return 0;
    }
    ga_instruction_eval_derivative_OP(base_tensor &t_,
                                      const ga_nonlinear_operator &OP_,
                                      ga_nonlinear_operator::arg_list &args_,
                                      size_type der1_)
      : t(t_), OP(OP_), args(args_), der1(der1_) {}
  };

  struct ga_instruction_eval_second_derivative_OP : public ga_instruction {
    base_tensor &t;
    const ga_nonlinear_operator &OP;
    ga_nonlinear_operator::arg_list args;
    size_type der1, der2;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: operator second derivative evaluation");
      OP.second_derivative(args, der1, der2, t);
      return 0;
    }
    ga_instruction_eval_second_derivative_OP
    (base_tensor &t_, const ga_nonlinear_operator &OP_,
     ga_nonlinear_operator::arg_list &args_, size_type der1_, size_type der2_)
      : t(t_), OP(OP_), args(args_), der1(der1_), der2(der2_) {}
  };

  struct ga_instruction_tensor_slice : public ga_instruction {
    base_tensor &t, &tc1;
    bgeot::multi_index mi, indices;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: tensor slice");
      size_type order = t.sizes().size();
      for (bgeot::multi_index mi3(order); !mi3.finished(t.sizes());
           mi3.incrementation(t.sizes())) {
        for (size_type j = 0; j < order; ++j)
          mi[indices[j]] = mi3[j];
        t(mi3) = tc1(mi);
      }
      return 0;
    }
    ga_instruction_tensor_slice(base_tensor &t_, base_tensor &tc1_,
                                bgeot::multi_index &mi_,
                                bgeot::multi_index &indices_)
      : t(t_), tc1(tc1_), mi(mi_), indices(indices_)  {}
  };

  struct ga_instruction_transformation_call : public ga_instruction {
    ga_workspace &workspace;
    ga_instruction_set::interpolate_info &inin;
    pinterpolate_transformation trans;
    fem_interpolation_context &ctx;
    base_small_vector &Normal;
    const mesh &m;
    base_matrix G;
    bool compute_der;

    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: call interpolate transformation");
      base_node P_ref;
      size_type cv;
      short_type face_num;
      gmm::clear(inin.Normal);
      inin.pt_type = trans->transform(workspace, m, ctx, Normal, &(inin.m), cv,
                                      face_num, P_ref, inin.Normal,
                                      inin.derivatives, compute_der);
      if (inin.pt_type) {
        if (cv != size_type(-1)) {
          bgeot::vectors_to_base_matrix(G, (inin.m)->points_of_convex(cv));
          inin.ctx = fem_interpolation_context((inin.m)->trans_of_convex(cv),
                                               0, P_ref, G, cv, face_num);
          inin.has_ctx = true;
          if (face_num != short_type(-1)) {
            inin.Normal = bgeot::compute_normal(inin.ctx, face_num);
            gmm::scale(inin.Normal, 1.0/gmm::vect_norm2(inin.Normal));
          } else
            inin.Normal.resize(0);
          inin.pt_y = inin.ctx.xreal();
        } else {
          inin.ctx = fem_interpolation_context();
          inin.pt_y = P_ref;
          inin.has_ctx = false;
        }
      } else {
        inin.ctx = fem_interpolation_context();
        inin.Normal.resize(0);
        inin.pt_y.resize(0);
        inin.has_ctx = false;
      }
      GA_DEBUG_INFO("Instruction: end of call interpolate transformation");
      return 0;
    }
    ga_instruction_transformation_call
    (ga_workspace &w, ga_instruction_set::interpolate_info &i,
     pinterpolate_transformation t, fem_interpolation_context &ctxx,
     base_small_vector &No, const mesh &mm, bool compute_der_)
      : workspace(w), inin(i), trans(t), ctx(ctxx), Normal(No), m(mm),
        compute_der(compute_der_) {}
  };


  struct ga_instruction_scalar_assembly : public ga_instruction {
    base_tensor &t;
    scalar_type &E, &coeff;
     virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: scalar term assembly\n");
      E += t[0] * coeff;
      return 0;
     }
    ga_instruction_scalar_assembly(base_tensor &t_, scalar_type &E_,
                                   scalar_type &coeff_)
      : t(t_), E(E_), coeff(coeff_) {}
  };

  struct ga_instruction_fem_vector_assembly : public ga_instruction {
    base_tensor &t;
    base_vector &Vr, &Vn;
    fem_interpolation_context &ctx;
    const gmm::sub_interval &Ir, &In;
    const mesh_fem *mfn, **mfg;
    scalar_type &coeff;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: vector term assembly for fem variable\n");
      const mesh_fem &mf = *(mfg ? *mfg : mfn);
      GMM_ASSERT1(&mf, "Internal error");
      const gmm::sub_interval &I = mf.is_reduced() ? Ir : In;
      base_vector &V = mf.is_reduced() ? Vr : Vn;
      size_type cv_1 = ctx.is_convex_num_valid()
        ? ctx.convex_num() : mf.convex_index().first_true();
      GA_DEBUG_ASSERT(V.size() >= I.first() + mf.nb_basic_dof(),
                      "Bad assembly vector size");
      mesh_fem::ind_dof_ct ct = mf.ind_basic_dof_of_element(cv_1);
      for (size_type i = 0; i < ct.size(); ++i)
        V[I.first()+ct[i]] += t[i] * coeff;
      return 0;
    }
    ga_instruction_fem_vector_assembly(base_tensor &t_, base_vector &Vr_,
                                       base_vector &Vn_,
                                       fem_interpolation_context &ctx_,
                                       const gmm::sub_interval &Ir_,
                                       const gmm::sub_interval &In_,
                                       const mesh_fem *mfn_,
                                       const mesh_fem **mfg_,
                                       scalar_type &coeff_)
      : t(t_), Vr(Vr_), Vn(Vn_), ctx(ctx_), Ir(Ir_), In(In_), mfn(mfn_),
        mfg(mfg_), coeff(coeff_) {}
  };

  struct ga_instruction_vector_assembly : public ga_instruction {
    base_tensor &t;
    base_vector &V;
    const gmm::sub_interval &I;
    scalar_type &coeff;
     virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: vector term assembly for "
                    "fixed size variable\n");
      gmm::add(gmm::scaled(t.as_vector(), coeff), gmm::sub_vector(V, I));
      return 0;
     }
    ga_instruction_vector_assembly(base_tensor &t_, base_vector &V_,
                                   const gmm::sub_interval &I_,
                                   scalar_type &coeff_)
      : t(t_), V(V_), I(I_), coeff(coeff_) {}
  };

  template <class MAT>
  struct ga_instruction_matrix_assembly : public ga_instruction {
    base_tensor &t;
    MAT &Kr, &Kn;
    fem_interpolation_context &ctx1, &ctx2;
    const gmm::sub_interval &Ir1, &Ir2;
    const gmm::sub_interval &In1, &In2;
    const mesh_fem *mfn1, *mfn2;
    const mesh_fem **mfg1, **mfg2;
    const scalar_type &coeff, &alpha1, &alpha2;
    size_type &nbpt, &ipt;
    base_vector &elem;
    bool interpolate;
    virtual int exec(void) {
      GA_DEBUG_INFO("Instruction: matrix term assembly\n");
      const mesh_fem &mf1 = *(mfg1 ? *mfg1 : mfn1);
      const mesh_fem &mf2 = *(mfg2 ? *mfg2 : mfn2);
      bool reduced = (&mf1 && mf1.is_reduced()) || (&mf2 && mf2.is_reduced());
      const gmm::sub_interval &I1 = reduced ? Ir1 : In1;
      const gmm::sub_interval &I2 = reduced ? Ir2 : In2;
      MAT &K = reduced ? Kr : Kn;
      if (ipt == 0 || interpolate) {
        gmm::resize(elem, t.size());
        gmm::copy(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      } else {
        gmm::add(gmm::scaled(t.as_vector(), coeff*alpha1*alpha2), elem);
      }
      if (ipt == nbpt-1 || interpolate) {
        size_type s1 = t.sizes()[0], s2 = t.sizes()[1];
        mesh_fem::ind_dof_ct ct1;
        if (&mf1) {
          if (!ctx1.is_convex_num_valid()) return 0;
          ct1 = mf1.ind_basic_dof_of_element(ctx1.convex_num());
          GA_DEBUG_ASSERT(ct1.size() == s1, "Internal error");
        }
        mesh_fem::ind_dof_ct ct2;
        if (&mf2) {
          if (!ctx2.is_convex_num_valid()) return 0;
          ct2 = mf2.ind_basic_dof_of_element(ctx2.convex_num());
          GA_DEBUG_ASSERT(ct2.size() == s2,
                          "Internal error, " << ct2.size() << " != " << s2);
        }

        GA_DEBUG_ASSERT(I1.size() && I2.size(), "Internal error");

        scalar_type ninf = gmm::vect_norminf(elem);
        if (ninf == scalar_type(0)) return 0;
        scalar_type threshold = ninf * 1E-14;
        for (size_type i1 = 0; i1 < s1; ++i1)
          for (size_type i2 = 0; i2 < s2; ++i2) {
            scalar_type e = elem[i2*s1+i1];
            if (gmm::abs(e) > threshold) {
              size_type j1 = I1.first()+((&mf1) ? ct1[i1] : i1);
              size_type j2 = I2.first()+((&mf2) ? ct2[i2] : i2);
              K(j1, j2) += e;
            }
          }
      }
      return 0;
    }
    ga_instruction_matrix_assembly(base_tensor &t_,
                                   MAT &Kr_, MAT &Kn_,
                                   fem_interpolation_context &ctx1_,
                                   fem_interpolation_context &ctx2_,
                                   const gmm::sub_interval &Ir1_,
                                   const gmm::sub_interval &In1_,
                                   const gmm::sub_interval &Ir2_,
                                   const gmm::sub_interval &In2_,
                                   const mesh_fem *mfn1_,
                                   const mesh_fem **mfg1_,
                                   const mesh_fem *mfn2_,
                                   const mesh_fem **mfg2_,
                                   const scalar_type &coeff_,
                                   const scalar_type &alpha2_,
                                   const scalar_type &alpha1_,
                                   size_type &nbpt_,
                                   size_type &ipt_,  base_vector &elem_,
                                   bool interpolate_)
      : t(t_), Kr(Kr_), Kn(Kn_), ctx1(ctx1_), ctx2(ctx2_),
        Ir1(Ir1_), Ir2(Ir2_), In1(In1_), In2(In2_),
        mfn1(mfn1_), mfn2(mfn2_), mfg1(mfg1_), mfg2(mfg2_),
        coeff(coeff_), alpha1(alpha1_), alpha2(alpha2_),
        nbpt(nbpt_), ipt(ipt_), elem(elem_), interpolate(interpolate_) {}
  };


  const mesh ga_workspace::dummy_mesh = mesh();
  const mesh_im ga_workspace::dummy_mim = mesh_im();
  const mesh_region ga_workspace::dummy_region = mesh_region();

  //=========================================================================
  // Structure dealing with user defined environment : constant, variables,
  // functions, operators.
  //=========================================================================

  void ga_workspace::add_fem_variable
  (const std::string &name, const mesh_fem &mf,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    variables[name] = var_description(true, true, &mf, I, &VV, 0, 1);
  }

  void ga_workspace::add_fixed_size_variable
  (const std::string &name,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    variables[name] = var_description(true, false, 0, I, &VV, 0,
                                      dim_type(gmm::vect_size(VV)));
  }

  void ga_workspace::add_fem_constant
  (const std::string &name, const mesh_fem &mf,
   const model_real_plain_vector &VV) {
    GMM_ASSERT1(mf.nb_dof(), "The provided mesh_fem of variable" << name
                             << "has zero degrees of freedom.");
    size_type Q = gmm::vect_size(VV)/mf.nb_dof();
    if (Q == 0) Q = size_type(1);
    variables[name] = var_description(false, true, &mf,
                                      gmm::sub_interval(), &VV, 0, Q);
  }

  void ga_workspace::add_fixed_size_constant
  (const std::string &name, const model_real_plain_vector &VV) {
    variables[name] = var_description(false, false, 0,
                                      gmm::sub_interval(), &VV, 0,
                                      gmm::vect_size(VV));
  }

  void ga_workspace::add_im_data(const std::string &name, const im_data &imd,
                                 const model_real_plain_vector &VV) {
    variables[name] = var_description
      (false, false, 0, gmm::sub_interval(), &VV, &imd,
       gmm::vect_size(VV)/(imd.nb_filtered_index() * imd.nb_tensor_elem()));
  }

  const mesh_region &ga_workspace::register_region(const mesh &m,
                                                   const mesh_region &region) {
    if (&m == &dummy_mesh) return dummy_region;

    std::list<mesh_region> &lmr = registred_mims[&m];
    for (std::list<mesh_region>::iterator it = lmr.begin();
         it != lmr.end(); ++it) {
      if (it->compare(m, region, m)) return *it;
    }
    lmr.push_back(region);
    return lmr.back();
  }

  bgeot::multi_index ga_workspace::qdims(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) {
      const mesh_fem *mf =  it->second.is_fem_dofs ? it->second.mf : 0;
      const im_data *imd = it->second.imd;
      size_type n = it->second.qdim();
      if (mf) {
        bgeot::multi_index mi = mf->get_qdims();
        if (n > 1 || it->second.qdims.size() > 1) {
          size_type i = 0;
          if (mi.back() == 1) { mi.back() *= it->second.qdims[0]; ++i; }
          for (; i < it->second.qdims.size(); ++i)
            mi.push_back(it->second.qdims[i]);
        }
        return mi;
      } else if (imd) {
        bgeot::multi_index mi = imd->tensor_size();
        size_type q = n / imd->nb_filtered_index();
        GMM_ASSERT1(q % imd->nb_tensor_elem() == 0,
                    "Invalid mesh im data vector");
        if (n > 1 || it->second.qdims.size() > 1) {
          size_type i = 0;
          if (mi.back() == 1) { mi.back() *= it->second.qdims[0]; ++i; }
          for (; i < it->second.qdims.size(); ++i)
            mi.push_back(it->second.qdims[i]);
        }
        return mi;
      }
      return it->second.qdims;
    }
    if (md && md->variable_exists(name))
      return md->qdims_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->qdims(name);
    if (variable_group_exists(name))
      return qdims(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  size_type ga_workspace::qdim(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) {
      const mesh_fem *mf =  it->second.is_fem_dofs ? it->second.mf : 0;
      const im_data *imd = it->second.imd;
      size_type n = it->second.qdim();
      if (mf) {
        return n * mf->get_qdim();
      } else if (imd) {
        return n * imd->tensor_size().total_size();
      }
      return n;
    }
    if (md && md->variable_exists(name))
      return md->qdim_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->qdim(name);
    if (variable_group_exists(name))
      return qdim(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  static bool ga_extract_variables(pga_tree_node pnode,
                                   ga_workspace &workspace,
                                   const mesh &m,
                                   std::set<var_trans_pair> &vars,
                                   bool ignore_data) {
    bool expand_groups = !ignore_data;
    bool found_var = false;
    if (pnode->node_type == GA_NODE_VAL ||
        pnode->node_type == GA_NODE_GRAD ||
        pnode->node_type == GA_NODE_HESS ||
        pnode->node_type == GA_NODE_DIVERG ||
        pnode->node_type == GA_NODE_INTERPOLATE_VAL ||
        pnode->node_type == GA_NODE_INTERPOLATE_GRAD ||
        pnode->node_type == GA_NODE_INTERPOLATE_HESS ||
        pnode->node_type == GA_NODE_INTERPOLATE_DIVERG ||
        pnode->node_type == GA_NODE_ELEMENTARY_VAL ||
        pnode->node_type == GA_NODE_ELEMENTARY_GRAD ||
        pnode->node_type == GA_NODE_ELEMENTARY_HESS ||
        pnode->node_type == GA_NODE_ELEMENTARY_DIVERG ||
        pnode->node_type == GA_NODE_XFEM_PLUS_VAL ||
        pnode->node_type == GA_NODE_XFEM_PLUS_GRAD ||
        pnode->node_type == GA_NODE_XFEM_PLUS_HESS ||
        pnode->node_type == GA_NODE_XFEM_PLUS_DIVERG ||
        pnode->node_type == GA_NODE_XFEM_MINUS_VAL ||
        pnode->node_type == GA_NODE_XFEM_MINUS_GRAD ||
        pnode->node_type == GA_NODE_XFEM_MINUS_HESS ||
        pnode->node_type == GA_NODE_XFEM_MINUS_DIVERG) {
      bool group = workspace.variable_group_exists(pnode->name);
      bool iscte = (!group) && workspace.is_constant(pnode->name);
      if (!iscte) found_var = true;
      if (!ignore_data || !iscte) {
        if (group && expand_groups) {
          const std::vector<std::string> &t
            = workspace.variable_group(pnode->name);
          for (size_type i = 0; i < t.size(); ++i)
            vars.insert(var_trans_pair(t[i], pnode->interpolate_name));

        } else
          vars.insert(var_trans_pair(pnode->name, pnode->interpolate_name));
      }
    }
    if (pnode->node_type == GA_NODE_INTERPOLATE_VAL ||
        pnode->node_type == GA_NODE_INTERPOLATE_GRAD ||
        pnode->node_type == GA_NODE_INTERPOLATE_HESS ||
        pnode->node_type == GA_NODE_INTERPOLATE_DIVERG ||
        pnode->node_type == GA_NODE_INTERPOLATE_VAL_TEST ||
        pnode->node_type == GA_NODE_INTERPOLATE_GRAD_TEST ||
        pnode->node_type == GA_NODE_INTERPOLATE_HESS_TEST ||
        pnode->node_type == GA_NODE_INTERPOLATE_DIVERG_TEST ||
        pnode->node_type == GA_NODE_INTERPOLATE_NORMAL ||
        pnode->node_type == GA_NODE_INTERPOLATE_X) {
      workspace.interpolate_transformation(pnode->interpolate_name)
        ->extract_variables(workspace, vars, ignore_data, m,
                            pnode->interpolate_name);
    }
    for (size_type i = 0; i < pnode->children.size(); ++i)
      found_var = ga_extract_variables(pnode->children[i], workspace, m,
                                       vars, ignore_data)
        || found_var;
    return found_var;
  }

  static size_type ref_elt_dim_of_mesh(const mesh &m) {
    if (m.convex_index().card())
      return m.trans_of_convex(m.convex_index().first())->dim();
    else
      return size_type(0);
  }

  void ga_workspace::add_tree(ga_tree &tree, const mesh &m,
                              const mesh_im &mim, const mesh_region &rg,
                              const std::string &expr,
                              size_type add_derivative_order,
                              bool function_expr) {
    if (tree.root) {

      // Eliminate the term if it corresponds to disabled variables
      if ((tree.root->test_function_type >= 1 &&
           is_disabled_variable(tree.root->name_test1)) ||
          (tree.root->test_function_type >= 2 &&
           is_disabled_variable(tree.root->name_test2))) {
        // cout << "disabling term ";  ga_print_node(tree.root, cout); cout << endl;
        return;
      }
      // cout << "add tree with tests functions of " <<  tree.root->name_test1
      //      << " and " << tree.root->name_test2 << endl;
      //      ga_print_node(tree.root, cout); cout << endl;
      bool remain = true;
      size_type order = 0, ind_tree = 0;

      switch(tree.root->test_function_type) {
      case 0: order = 0; break;
      case 1: order = 1; break;
      case 3: order = 2; break;
      default: GMM_ASSERT1(false, "Inconsistent term "
                           << tree.root->test_function_type);
      }

      bool found = false;
      for (size_type i = 0; i < trees.size(); ++i) {
        if (trees[i].mim == &mim && trees[i].m == &m &&
            trees[i].order == order &&
            trees[i].name_test1.compare(tree.root->name_test1) == 0 &&
            trees[i].interpolate_name_test1.compare
            (tree.root->interpolate_name_test1) == 0 &&
            trees[i].name_test2.compare(tree.root->name_test2) == 0 &&
            trees[i].interpolate_name_test2.compare
            (tree.root->interpolate_name_test2) == 0 &&
            trees[i].rg == &rg) {
          ga_tree &ftree = *(trees[i].ptree);

          ftree.insert_node(ftree.root, GA_NODE_OP);
          ftree.root->op_type = GA_PLUS;
          ftree.add_child(ftree.root);
          ftree.copy_node(tree.root, ftree.root, ftree.root->children[1]);
          ga_semantic_analysis("", ftree, *this, m.dim(),
                               ref_elt_dim_of_mesh(m), false, function_expr);
          found = true;
          break;
        }
      }

      if (!found) {
        ind_tree = trees.size(); remain = false;
        trees.push_back(tree_description());
        trees.back().mim = &mim; trees.back().m = &m;
        trees.back().rg = &rg;
        trees.back().ptree = new ga_tree;
        trees.back().ptree->swap(tree);
        pga_tree_node root = trees.back().ptree->root;
        trees.back().name_test1 = root->name_test1;
        trees.back().name_test2 = root->name_test2;
        trees.back().interpolate_name_test1 = root->interpolate_name_test1;
        trees.back().interpolate_name_test2 = root->interpolate_name_test2;
        trees.back().order = order;
      }

      if (order < add_derivative_order) {
        std::set<var_trans_pair> expr_variables;
        ga_extract_variables((remain ? tree : *(trees[ind_tree].ptree)).root,
                             *this, m, expr_variables, true);

        for (std::set<var_trans_pair>::iterator it = expr_variables.begin();
             it != expr_variables.end(); ++it) {
          if (!(is_constant(it->first))) {
            ga_tree dtree = (remain ? tree : *(trees[ind_tree].ptree));
            // cout << "Derivation with respect to " << it->first << " : "
            //     << it->second << " of " << ga_tree_to_string(dtree) << endl;
            GA_TIC;
            ga_derivative(dtree, *this, m, it->first, it->second, 1+order);
            // cout << "Result : " << ga_tree_to_string(dtree) << endl;
            GA_TOCTIC("Derivative time");
            ga_semantic_analysis(expr, dtree, *this, m.dim(),
                                 ref_elt_dim_of_mesh(m), false, function_expr);
            GA_TOCTIC("Analysis after Derivative time");
            // cout << "after analysis "  << ga_tree_to_string(dtree) << endl;
            add_tree(dtree, m, mim, rg, expr, add_derivative_order,
                     function_expr);
          }
        }
      }
    }
  }

  ga_workspace::m_tree::~m_tree(void) { if (ptree) delete ptree; }
  ga_workspace::m_tree::m_tree(const m_tree& o)
    : ptree(o.ptree), meshdim(o.meshdim), ignore_X(o.ignore_X)
  { if (o.ptree) ptree = new ga_tree(*(o.ptree)); }
  ga_workspace::m_tree &ga_workspace::m_tree::operator =(const m_tree& o) {
    ptree = o.ptree; meshdim = o.meshdim; ignore_X = o.ignore_X;
    if (o.ptree) ptree = new ga_tree(*(o.ptree));
    return *this;
  }

  ga_tree &ga_workspace::macro_tree(const std::string &name,
                                    size_type meshdim, size_type ref_elt_dim,
                                    bool ignore_X) const {
    GMM_ASSERT1(macro_exists(name), "Undefined macro");
    auto it = macro_trees.find(name);
    bool to_be_analyzed = false;
    m_tree *mt = 0;

    if (it == macro_trees.end()) {
      mt = &(macro_trees[name]);
      to_be_analyzed = true;
    } else {
      mt = &(it->second);
      GMM_ASSERT1(mt->ptree, "Recursive definition of macro " << name);
      if (mt->meshdim != meshdim || mt->ignore_X != ignore_X) {
        to_be_analyzed = true;
        delete mt->ptree; mt->ptree = 0;
      }
    }
    if (to_be_analyzed) {
      ga_tree tree;
      ga_read_string(get_macro(name), tree);
      ga_semantic_analysis(get_macro(name), tree, *this, meshdim, ref_elt_dim,
                           false, ignore_X);
      GMM_ASSERT1(tree.root, "Invalid macro");
      mt->ptree = new ga_tree(tree);
      mt->meshdim = meshdim;
      mt->ignore_X = ignore_X;
    }
    return *(mt->ptree);
  }

  size_type ga_workspace::add_expression(const std::string expr,
                                         const mesh_im &mim,
                                         const mesh_region &rg_,
                                         size_type add_derivative_order) {
    const mesh_region &rg = register_region(mim.linked_mesh(), rg_);
    // cout << "adding expression " << expr << endl;
    GA_TIC;
    size_type max_order = 0;
    std::vector<ga_tree> ltrees(1);
    ga_read_string(expr, ltrees[0]);
    // cout << "read : " << ga_tree_to_string(ltrees[0])  << endl;
    ga_semantic_analysis(expr, ltrees[0], *this, mim.linked_mesh().dim(),
                         ref_elt_dim_of_mesh(mim.linked_mesh()),
                         false, false, 1);
    // cout << "analysed : " << ga_tree_to_string(ltrees[0]) << endl;
    GA_TOC("First analysis time");
    if (ltrees[0].root) {
      if (test1.size() > 1 || test2.size() > 1) {
        size_type ntest2 = std::max(size_type(1), test2.size());
        size_type nb_ltrees = test1.size()*ntest2;
        ltrees.resize(nb_ltrees);
        for (size_type i = 1; i < nb_ltrees; ++i) ltrees[i] = ltrees[0];
        std::set<var_trans_pair>::iterator it1 = test1.begin();
        for (size_type i = 0; i < test1.size(); ++i, ++it1) {
          std::set<var_trans_pair>::iterator it2 = test2.begin();
          for (size_type j = 0; j < ntest2; ++j) {
            selected_test1 = *it1;
            if (test2.size()) selected_test2 = *it2++;
            // cout << "analysis with " << selected_test1.first << endl;
            ga_semantic_analysis(expr, ltrees[i*ntest2+j], *this,
                                 mim.linked_mesh().dim(),
                                 ref_elt_dim_of_mesh(mim.linked_mesh()),
                                 false, false, 2);
            // cout <<"split: "<< ga_tree_to_string(ltrees[i*ntest2+j]) << endl;
          }
        }
      }

      for (size_type i = 0; i < ltrees.size(); ++i) {
        if (ltrees[i].root) {
          // cout << "adding tree " << ga_tree_to_string(ltrees[i]) << endl;
          max_order = std::max(ltrees[i].root->nb_test_functions(), max_order);
          add_tree(ltrees[i], mim.linked_mesh(), mim, rg, expr,
                   add_derivative_order);
        }
      }
    }
    GA_TOC("Time for add expression");
    return max_order;
  }

  void ga_workspace::add_function_expression(const std::string expr) {
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, 1, 1, false, true);
    if (tree.root) {
      // GMM_ASSERT1(tree.root->nb_test_functions() == 0,
      //            "Invalid function expression");
      add_tree(tree, dummy_mesh, dummy_mim, dummy_region, expr, 0);
    }
  }

  void ga_workspace::add_interpolation_expression(const std::string expr,
                                                  const mesh &m,
                                                  const mesh_region &rg_) {
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, m.dim(), ref_elt_dim_of_mesh(m),
                         false, false);
    if (tree.root) {
      GMM_ASSERT1(tree.root->nb_test_functions() == 0,
                  "Invalid expression containing test functions");
      add_tree(tree, m, dummy_mim, rg, expr, 0, false);
    }
  }

  void ga_workspace::add_interpolation_expression(const std::string expr,
                                                  const mesh_im &mim,
                                                  const mesh_region &rg_) {
    const mesh &m = mim.linked_mesh();
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, m.dim(), ref_elt_dim_of_mesh(m),
                         false, false);
    if (tree.root) {
      GMM_ASSERT1(tree.root->nb_test_functions() == 0,
                  "Invalid expression containing test functions");
      add_tree(tree, m, mim, rg, expr, 0, false);
    }
  }

  size_type ga_workspace::nb_trees(void) const { return trees.size(); }

  ga_workspace::tree_description &ga_workspace::tree_info(size_type i)
  { return trees[i]; }

  bool ga_workspace::used_variables(model::varnamelist &vl,
                                    model::varnamelist &vl_test1,
                                    model::varnamelist &vl_test2,
                                    model::varnamelist &dl,
                                    size_type order) {
    bool islin = true;
    std::set<var_trans_pair> vll, dll;
    for (size_type i = 0; i < vl.size(); ++i)
      vll.insert(var_trans_pair(vl[i], ""));
    for (size_type i = 0; i < dl.size(); ++i)
      dll.insert(var_trans_pair(dl[i], ""));

    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      std::set<var_trans_pair> dllaux;
      bool fv = ga_extract_variables(td.ptree->root, *this, *(td.m),
                                     dllaux, false);

      if (td.order == order) {
        for (std::set<var_trans_pair>::iterator it = dllaux.begin();
             it!=dllaux.end(); ++it)
          dll.insert(*it);
      }
      switch (td.order) {
      case 0:  break;
      case 1:
        if (td.order == order) {
          if (variable_group_exists(td.name_test1)) {
            const std::vector<std::string> &t= variable_group(td.name_test1);
            for (size_type j = 0; j < t.size(); ++j) {
              vll.insert(var_trans_pair(t[j], td.interpolate_name_test1));
            }
          } else {
            vll.insert(var_trans_pair(td.name_test1,
                                      td.interpolate_name_test1));
          }
          bool found = false;
          for (size_type j = 0; j < vl_test1.size(); ++j)
            if (td.name_test1.compare(vl_test1[j]) == 0)
              found = true;
          if (!found)
            vl_test1.push_back(td.name_test1);
        }
        break;
      case 2:
        if (td.order == order) {
          if (variable_group_exists(td.name_test1)) {
            const std::vector<std::string> &t= variable_group(td.name_test1);
            for (size_type j = 0; j < t.size(); ++j) {
              vll.insert(var_trans_pair(t[j], td.interpolate_name_test1));
              }
          } else {
            vll.insert(var_trans_pair(td.name_test1,
                                      td.interpolate_name_test1));
          }
          if (variable_group_exists(td.name_test2)) {
            const std::vector<std::string> &t= variable_group(td.name_test2);
            for (size_type j = 0; j < t.size(); ++j) {
              vll.insert(var_trans_pair(t[j], td.interpolate_name_test2));
            }
          } else {
            vll.insert(var_trans_pair(td.name_test2,
                                      td.interpolate_name_test2));
          }
          bool found = false;
          for (size_type j = 0; j < vl_test1.size(); ++j)
            if ((td.name_test1.compare(vl_test1[j]) == 0) &&
                (td.name_test2.compare(vl_test2[j]) == 0))
              found = true;
          if (!found) {
            vl_test1.push_back(td.name_test1);
            vl_test2.push_back(td.name_test2);
          }
        }
        if (fv) islin = false;
        break;
      }
    }
    vl.clear();
    for (auto it = vll.begin(); it!=vll.end(); ++it)
      if (vl.size() == 0 || it->first.compare(vl.back()))
        vl.push_back(it->first);
    dl.clear();
    for (auto it = dll.begin(); it!=dll.end(); ++it)
      if (dl.size() == 0 || it->first.compare(dl.back()))
        dl.push_back(it->first);

    return islin;
  }

  void ga_workspace::define_variable_group(const std::string &group_name,
                                         const std::vector<std::string> &nl) {
    GMM_ASSERT1(!(variable_exists(group_name)), "The name of a group of "
                "variables cannot be the same as a variable name");

    std::set<const mesh *> ms;
    bool is_data_ = false;
    for (size_type i = 0; i < nl.size(); ++i) {
      if (i == 0)
        is_data_ = is_constant(nl[i]);
      else {
        GMM_ASSERT1(is_data_ == is_constant(nl[i]),
                    "It is not possible to mix variables and data in a group");
      }
      GMM_ASSERT1(variable_exists(nl[i]),
                  "All variables in a group have to exist in the model");
      const mesh_fem *mf = associated_mf(nl[i]);
      GMM_ASSERT1(mf, "Variables in a group should be fem variables");
      GMM_ASSERT1(ms.find(&(mf->linked_mesh())) == ms.end(),
                  "Two variables in a group cannot share the same mesh");
      ms.insert(&(mf->linked_mesh()));
    }
    variable_groups[group_name] = nl;
  }


  const std::string &ga_workspace::variable_in_group
  (const std::string &group_name, const mesh &m) const {
    if (variable_group_exists(group_name)) {
      const std::vector<std::string> &t=variable_group(group_name);
      for (size_type j = 0; j < t.size(); ++j)
        if (&((associated_mf(t[j]))->linked_mesh()) == &m) return t[j];
      GMM_ASSERT1(false, "No variable in this group for the given mesh");
    } else return group_name;
  }


  void ga_workspace::assembly(size_type order) {

    GA_TIC;
    ga_instruction_set gis;
    ga_compile(*this, gis, order);
    size_type ndof = gis.nb_dof, max_dof =  gis.max_dof;
    GA_TOCTIC("Compile time");

    if (order == 2) {
      K.resize(max_dof);
      gmm::clear(unreduced_K); gmm::resize(unreduced_K, ndof, ndof);
    }
    if (order == 1) {
      V.resize(max_dof);
      gmm::clear(unreduced_V); gmm::resize(unreduced_V, ndof);
    }
    E = 0;
    GA_TOCTIC("Init time");

    ga_exec(gis, *this);
    GA_TOCTIC("Exec time");

    if (order == 1) {
      MPI_SUM_VECTOR(assembled_vector());
      MPI_SUM_VECTOR(unreduced_V);
    }

    // Deal with reduced fems.
    if (order) {
      std::list<ga_tree>::iterator it = gis.trees.begin();
      std::set<std::string> vars_vec_done;
      std::set<std::pair<std::string, std::string> > vars_mat_done;
      for (; it != gis.trees.end(); ++it) {
        if (it->root) {
          if (order == 1) {
            const std::string &name = it->root->name_test1;
            bool is_group = variable_group_exists(name);
            const std::vector<std::string> *t
              = is_group ? &(variable_group(name)) : 0;
            const std::string *vname = &(is_group ? (*t)[0] : name);

            for (size_type i = 0; (i==0) || (is_group && i < t->size()); ++i) {
              if (i > 0) vname = &((*t)[i]);
              const mesh_fem *mf = associated_mf(*vname);
              if (mf && mf->is_reduced() &&
                  vars_vec_done.find(*vname) == vars_vec_done.end()) {
                gmm::mult_add(gmm::transposed(mf->extension_matrix()),
                              gmm::sub_vector(unreduced_V,
                                              gis.var_intervals[*vname]),
                              gmm::sub_vector(V(),
                                              interval_of_variable(*vname)));
                vars_vec_done.insert(*vname);
              }
            }
          } else {
            std::string &name1 = it->root->name_test1;
            std::string &name2 = it->root->name_test2;
            bool is_group1 = variable_group_exists(name1);
            const std::vector<std::string> *t1
              = is_group1 ? &(variable_group(name1)) : 0;
            const std::string *vname1 = &(is_group1 ? (*t1)[0] : name1);
            bool is_group2 = variable_group_exists(name2);
            const std::vector<std::string> *t2
              = is_group2 ? &(variable_group(name2)) : 0;
            const std::string *vname2 = &(is_group2 ? (*t2)[0] : name2);
            for (size_type i1 = 0; (i1==0) || (is_group1 && i1 < t1->size());
                 ++i1) {
              if (i1 > 0) vname1 = &((*t1)[i1]);
              for (size_type i2 = 0; (i2==0) || (is_group2 && i2 < t2->size());
                   ++i2) {
                if (i2 > 0) vname2 = &((*t2)[i2]);
                const mesh_fem *mf1 = associated_mf(*vname1);
                const mesh_fem *mf2 = associated_mf(*vname2);
                if (((mf1 && mf1->is_reduced())
                     || (mf2 && mf2->is_reduced()))) {
                  std::pair<std::string, std::string> p(*vname1, *vname2);
                  if (vars_mat_done.find(p) == vars_mat_done.end()) {
                    gmm::sub_interval uI1 = gis.var_intervals[*vname1];
                    gmm::sub_interval uI2 = gis.var_intervals[*vname2];
                    gmm::sub_interval I1 = interval_of_variable(*vname1);
                    gmm::sub_interval I2 = interval_of_variable(*vname2);
                    if ((mf1 && mf1->is_reduced()) &&
                        (mf2 && mf2->is_reduced())) {
                      model_real_sparse_matrix aux(I1.size(), uI2.size());
                      model_real_row_sparse_matrix M(I1.size(), I2.size());
                      gmm::mult(gmm::transposed(mf1->extension_matrix()),
                                gmm::sub_matrix(unreduced_K, uI1, uI2), aux);
                      gmm::mult(aux, mf2->extension_matrix(), M);
                      gmm::add(M, gmm::sub_matrix(K(), I1, I2));
                    } else if (mf1 && mf1->is_reduced()) {
                      model_real_sparse_matrix M(I1.size(), I2.size());
                      gmm::mult(gmm::transposed(mf1->extension_matrix()),
                                gmm::sub_matrix(unreduced_K, uI1, uI2), M);
                      gmm::add(M, gmm::sub_matrix(K(), I1, I2));
                    } else {
                      model_real_row_sparse_matrix M(I1.size(), I2.size());
                      gmm::mult(gmm::sub_matrix(unreduced_K, uI1, uI2),
                                mf2->extension_matrix(), M);
                      gmm::add(M, gmm::sub_matrix(K(), I1, I2));
                    }
                    vars_mat_done.insert(p);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  const gmm::sub_interval &
  ga_workspace::interval_of_disabled_variable(const std::string &name) const {
    std::map<std::string, gmm::sub_interval>::const_iterator
      it1 = int_disabled_variables.find(name);
    if (it1 != int_disabled_variables.end()) return it1->second;
    if (md->is_affine_dependent_variable(name))
      return interval_of_disabled_variable(md->org_variable(name));

    size_type first = md->nb_dof();
    for (std::map<std::string, gmm::sub_interval>::const_iterator it
           = int_disabled_variables.begin();
         it !=  int_disabled_variables.end(); ++it) {
      first = std::max(first, it->second.last());
    }
    int_disabled_variables[name]
      = gmm::sub_interval(first, gmm::vect_size(value(name)));
    return int_disabled_variables[name];
  }

  std::map<std::string, gmm::sub_interval> int_disabled_variables;

  void ga_workspace::clear_expressions(void) {
    trees.clear();
    macro_trees.clear();
  }

  void ga_workspace::print(std::ostream &str) {
    for (size_type i = 0; i < trees.size(); ++i)
      if (trees[i].ptree->root) {
        cout << "Expression tree " << i << " of order " <<
                trees[i].ptree->root->nb_test_functions() << " :" << endl;
        ga_print_node(trees[i].ptree->root, str);
        cout << endl;
      }
  }

  void ga_workspace::tree_description::copy(const tree_description& td) {
    order = td.order;
    name_test1 = td.name_test1;
    name_test2 = td.name_test2;
    interpolate_name_test1 = td.interpolate_name_test1;
    interpolate_name_test2 = td.interpolate_name_test2;
    mim = td.mim;
    m = td.m;
    rg = td.rg;
    ptree = 0;
    elem = td.elem;
    if (td.ptree) ptree = new ga_tree(*(td.ptree));
  }

  ga_workspace::tree_description &ga_workspace::tree_description::operator =
  (const ga_workspace::tree_description& td)
  { if (ptree) delete ptree; ptree = 0; copy(td); return *this; }
  ga_workspace::tree_description::~tree_description()
  { if (ptree) delete ptree; ptree = 0; }

  //=========================================================================
  // Some hash code functions for node identification
  //=========================================================================

  static scalar_type ga_hash_code(const std::string &s) {
    scalar_type c(0);
    for (size_type i = 0; i < s.size(); ++i)
      c += sin(M_E+scalar_type(s[i])+M_PI*scalar_type(i+1));
    return c;
  }

  static scalar_type ga_hash_code(const base_tensor &t) {
    scalar_type c(0);
    for (size_type i = 0; i < t.size(); ++i)
      c += sin(M_E+t[i]+M_E*M_E*scalar_type(i+1));
    return c;
  }

  static scalar_type ga_hash_code(GA_NODE_TYPE e) {
    return cos(M_E + scalar_type((e == GA_NODE_ZERO) ? GA_NODE_CONSTANT : e));
  }

  static scalar_type ga_hash_code(pga_tree_node pnode) {
    scalar_type c = ga_hash_code(pnode->node_type);

    switch (pnode->node_type) {
    case GA_NODE_CONSTANT: case GA_NODE_ZERO:
      c += ga_hash_code(pnode->t);
      if (pnode->test_function_type & 1)
        c += 34.731 * ga_hash_code(pnode->name_test1);
      if (pnode->test_function_type & 2)
        c += 34.731 * ga_hash_code(pnode->name_test2);
      break;

    case GA_NODE_OP: c += scalar_type(pnode->op_type)*M_E*M_PI*M_PI; break;
    case GA_NODE_X: c += scalar_type(pnode->nbc1) + M_E*M_PI; break;
    case GA_NODE_VAL: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
    case GA_NODE_VAL_TEST: case GA_NODE_GRAD_TEST:
    case GA_NODE_HESS_TEST: case GA_NODE_DIVERG_TEST:
      c += ga_hash_code(pnode->name); break;

    case GA_NODE_INTERPOLATE_FILTER:
      c += 1.73*ga_hash_code(pnode->interpolate_name)
        + 0.84*M_PI*scalar_type(pnode->nbc1);
      break;
    case GA_NODE_INTERPOLATE_DERIVATIVE:
      c += 2.321*ga_hash_code(pnode->interpolate_name_der);
      // No break. The hash code is completed with the next item
    case GA_NODE_INTERPOLATE_VAL: case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS: case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_INTERPOLATE_VAL_TEST: case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST: case GA_NODE_INTERPOLATE_DIVERG_TEST:
      c += 1.33*(1.22+ga_hash_code(pnode->name))
        + 1.66*ga_hash_code(pnode->interpolate_name);
      break;
    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_ELEMENTARY_VAL_TEST: case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
      c += 1.33*(1.22+ga_hash_code(pnode->name))
        + 2.63*ga_hash_code(pnode->elementary_name);
      break;
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL_TEST: case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL_TEST: case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST: case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      c += 1.33*(1.22+ga_hash_code(pnode->name));
      break;
    case GA_NODE_INTERPOLATE_NORMAL: case GA_NODE_INTERPOLATE_X:
      c += M_PI*1.33*ga_hash_code(pnode->interpolate_name);
      break;
    case GA_NODE_PREDEF_FUNC: case GA_NODE_SPEC_FUNC: case GA_NODE_OPERATOR:
      c += ga_hash_code(pnode->name)
        + tanh(scalar_type(pnode->der1)/M_PI + scalar_type(pnode->der2)*M_PI);
      break;
    default: break;
    }
    return c;
  }




  // 0 : ok
  // 1 : function or operator name or "X"
  // 2 : reserved prefix Grad, Hess, Div, Test and Test2
  // 3 : reserved prefix Dot and Previous
  int ga_check_name_validity(const std::string &name) {

    if (!(name.compare("X")) ||
        !(name.compare("Normal")) ||
        !(name.compare("Reshape")))
      return 1;


    if (name.compare(0, 11, "Derivative_") == 0)
      return 2;

    ga_predef_operator_tab &PREDEF_OPERATORS
      = dal::singleton<ga_predef_operator_tab>::instance(0);
    ga_predef_function_tab::const_iterator it=PREDEF_FUNCTIONS.find(name);
    if (it != PREDEF_FUNCTIONS.end())
      return 1;

    if (SPEC_FUNCTIONS.find(name) != SPEC_FUNCTIONS.end())
      return 1;

    if (PREDEF_OPERATORS.tab.find(name) != PREDEF_OPERATORS.tab.end())
      return 1;

    if (name.size() >= 5 && name.compare(0, 5, "Grad_") == 0)
      return 2;

    if (name.size() >= 5 && name.compare(0, 5, "Hess_") == 0)
      return 2;

    if (name.size() >= 4 && name.compare(0, 4, "Div_") == 0)
      return 2;

    if (name.size() >= 6 && name.compare(0, 6, "Test2_") == 0)
      return 2;

    if (name.size() >= 5 && name.compare(0, 5, "Test_") == 0)
      return 2;

//     if (name.size() >= 4 && name.compare(0, 4, "Dot_") == 0)
//       return 3;
//     if (name.size() >= 5 && name.compare(0, 5, "Dot2_") == 0)
//       return 3;

//     if (name.size() >= 9 && name.compare(0, 9, "Previous_") == 0)
//       return 3;
//     if (name.size() >= 10 && name.compare(0, 10, "Previous2_") == 0)
//       return 3;
//     if (name.size() >= 12 && name.compare(0, 12, "Previous1_2_") == 0)
//       return 3;


    return 0;
  }


  //=========================================================================
  // Semantic analysis, tree simplification and tree enrichment
  //    - Control tensor sizes for operations, operator or function call
  //    - Compute all constant operations (i.e. non element dependent)
  //    - Build a ready to use tree for derivation/compilation
  //=========================================================================
  // option = 0 : strict analysis,
  //          1 : do not complain about incompatible test functions but
  //              store them,
  //          2 : cut incompatible test function branches with respect to the
  //              one in workspace.selected_pair
  static void ga_semantic_analysis(const std::string &expr, ga_tree &tree,
                                   const ga_workspace &workspace,
                                   size_type meshdim,
                                   size_type ref_elt_dim,
                                   bool eval_fixed_size,
                                   bool ignore_X, int option) {
    GMM_ASSERT1(predef_functions_initialized, "Internal error");
    if (!(tree.root)) return;
    if (option == 1) { workspace.test1.clear(); workspace.test2.clear(); }
    // cout << "semantic analysis of " << ga_tree_to_string(tree) << endl;
    ga_node_analysis(expr, tree, workspace, tree.root, meshdim, ref_elt_dim,
                     eval_fixed_size, ignore_X, option);
    if (tree.root && option == 2) {
      if (((tree.root->test_function_type & 1) &&
           (tree.root->name_test1.compare(workspace.selected_test1.first)
            || tree.root->interpolate_name_test1.compare
            (workspace.selected_test1.second)))
          ||
          ((tree.root->test_function_type & 2) &&
           (tree.root->name_test2.compare(workspace.selected_test2.first)
            || tree.root->interpolate_name_test2.compare
            (workspace.selected_test2.second))))
        tree.clear();
    }
    // cout << "semantic analysis done " << endl;
    ga_valid_operand(expr, tree.root);
  }

  static void ga_node_analysis(const std::string &expr, ga_tree &tree,
                               const ga_workspace &workspace,
                               pga_tree_node pnode, size_type meshdim,
                               size_type ref_elt_dim, bool eval_fixed_size,
                               bool ignore_X, int option) {
    bool all_cte = true, all_sc = true;
    pnode->symmetric_op = false;

    for (size_type i = 0; i < pnode->children.size(); ++i) {
      ga_node_analysis(expr, tree, workspace, pnode->children[i], meshdim,
                       ref_elt_dim, eval_fixed_size, ignore_X, option);
      all_cte = all_cte && (pnode->children[i]->node_type == GA_NODE_CONSTANT);
      all_sc = all_sc && (pnode->children[i]->tensor_proper_size() == 1);
      GMM_ASSERT1(pnode->children[i]->test_function_type != size_type(-1),
                  "internal error on child " << i);
      if (pnode->node_type != GA_NODE_PARAMS)
        ga_valid_operand(expr, pnode->children[i]);
    }

    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bgeot::multi_index mi;
    const bgeot::multi_index &size0 = child0 ? child0->t.sizes() : mi;
    const bgeot::multi_index &size1 = child1 ? child1->t.sizes() : mi;
    size_type dim0 = child0 ? child0->tensor_order() : 0;
    size_type dim1 = child1 ? child1->tensor_order() : 0;

    // cout << "child1 = " << child1 << endl;
    // cout << "child0 = " << child0 << endl;
    // cout << "nbch = " << nbch << endl;
    // cout<<"begin analysis of node "; ga_print_node(pnode, cout); cout<<endl;


    switch (pnode->node_type) {
    case GA_NODE_PREDEF_FUNC: case GA_NODE_OPERATOR: case GA_NODE_SPEC_FUNC :
    case GA_NODE_CONSTANT: case GA_NODE_X: case GA_NODE_ELT_SIZE:
    case GA_NODE_ELT_K:  case GA_NODE_ELT_B:
    case GA_NODE_NORMAL: case GA_NODE_RESHAPE:
    case GA_NODE_INTERPOLATE_X: case GA_NODE_INTERPOLATE_NORMAL:
      pnode->test_function_type = 0; break;

    case GA_NODE_ALLINDICES: pnode->test_function_type = 0; break;
    case GA_NODE_VAL:
      if (eval_fixed_size && !(workspace.associated_mf(pnode->name))
          && !(workspace.associated_im_data(pnode->name))) {
        gmm::copy(workspace.value(pnode->name), pnode->t.as_vector());
        pnode->node_type = GA_NODE_CONSTANT;
      }
      break;

    case GA_NODE_ZERO: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
    case GA_NODE_INTERPOLATE_VAL:  case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS: case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
      break;

    case GA_NODE_VAL_TEST: case GA_NODE_GRAD_TEST:
    case GA_NODE_HESS_TEST: case GA_NODE_DIVERG_TEST:
    case GA_NODE_INTERPOLATE_VAL_TEST: case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST: case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_INTERPOLATE_DERIVATIVE:
    case GA_NODE_ELEMENTARY_VAL_TEST: case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST: case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST: case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST: case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      {
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        size_type t_type = pnode->test_function_type;
        if (t_type == 1) {
          pnode->name_test1 = pnode->name;
          pnode->interpolate_name_test1 = pnode->interpolate_name;
          pnode->interpolate_name_test2 = pnode->name_test2 = "";
          pnode->qdim1 = (mf ? workspace.qdim(pnode->name)
                          : gmm::vect_size(workspace.value(pnode->name)));
          if (option == 1)
            workspace.test1.insert
              (var_trans_pair(pnode->name_test1,
                              pnode->interpolate_name_test1));
          if (!(pnode->qdim1))
            ga_throw_error(expr, pnode->pos, "Invalid null size of variable");
        } else {
          pnode->interpolate_name_test1 = pnode->name_test1 = "";
          pnode->name_test2 = pnode->name;
          pnode->interpolate_name_test2 = pnode->interpolate_name;
          pnode->qdim2 = (mf ? workspace.qdim(pnode->name)
                          : gmm::vect_size(workspace.value(pnode->name)));
          if (option == 1)
            workspace.test2.insert
              (var_trans_pair(pnode->name_test2,
                              pnode->interpolate_name_test2));
          if (!(pnode->qdim2))
            ga_throw_error(expr, pnode->pos, "Invalid null size of variable");
        }
        if (!mf) {
          size_type n = workspace.qdim(pnode->name);
          if (!n)
            ga_throw_error(expr, pnode->pos, "Invalid null size of variable");
          if (n == 1) {
            pnode->init_vector_tensor(1);
            pnode->t[0] = scalar_type(1);
            pnode->test_function_type = t_type;
          } else {
            pnode->init_matrix_tensor(n,n);
            pnode->test_function_type = t_type;
            for (size_type i = 0; i < n; ++i)
              for (size_type j = 0; j < n; ++j)
                pnode->t(i,j) = (i == j) ? scalar_type(1) : scalar_type(0);
          }
        }
      }
      break;

    case GA_NODE_INTERPOLATE:
      if (!(pnode->name.compare("Normal"))) {
        pnode->node_type = GA_NODE_INTERPOLATE_NORMAL;
        pnode->init_vector_tensor(meshdim);
        break;
      }
      if (!(pnode->name.compare("X"))) {
        pnode->node_type = GA_NODE_INTERPOLATE_X;
        pnode->init_vector_tensor(meshdim);
        break;
      }
      // else continue with what follows
    case GA_NODE_ELEMENTARY: // and ... case GA_NODE_INTERPOLATE:
    case GA_NODE_XFEM_PLUS:
    case GA_NODE_XFEM_MINUS:
      {
        int ndt = (pnode->node_type == GA_NODE_INTERPOLATE) ? 1 : 0
          + (pnode->node_type == GA_NODE_ELEMENTARY) ? 2 : 0
          + (pnode->node_type == GA_NODE_XFEM_PLUS) ? 3 : 0
          + (pnode->node_type == GA_NODE_XFEM_MINUS) ? 4 : 0;
        std::string op__name =
          (pnode->node_type == GA_NODE_INTERPOLATE) ? "Interpolation" : ""
          + (pnode->node_type == GA_NODE_ELEMENTARY) ? "Elementary transformation" : ""
          + (pnode->node_type == GA_NODE_XFEM_PLUS) ? "Xfem_plus" : ""
          + (pnode->node_type == GA_NODE_XFEM_MINUS) ? "Xfem_minus" : "";

        std::string name = pnode->name;
        size_type prefix_id = ga_parse_prefix_operator(name);
        size_type test = ga_parse_prefix_test(name);
        pnode->name = name;

        // Group must be tested and it should be a fem variable
        if (!(workspace.variable_or_group_exists(name)))
            ga_throw_error(expr, pnode->pos,
                           "Unknown variable or group of variables");

        const mesh_fem *mf = workspace.associated_mf(name);
        if (!mf)
          ga_throw_error(expr, pnode->pos, op__name
                        << " can only apply to finite element variables/data");

        size_type q = workspace.qdim(name), n = mf->linked_mesh().dim();
        if (!q) ga_throw_error(expr, pnode->pos,
                               "Invalid null size of variable");

        bgeot::multi_index mii = workspace.qdims(name);
        if (mii.size() > 6)
          ga_throw_error(expr, pnode->pos,
                         "Tensor with too many dimensions. Limited to 6");

        if (test == 1) {
          pnode->name_test1 = name;
          pnode->interpolate_name_test1 = pnode->interpolate_name;
          if (option == 1)
            workspace.test1.insert
              (var_trans_pair(pnode->name_test1,
                              pnode->interpolate_name_test1));
          pnode->qdim1 = workspace.qdim(name);
          if (!(pnode->qdim1))
            ga_throw_error(expr, pnode->pos,
                           "Invalid null size of variable");
        } else if (test == 2) {
          pnode->name_test2 = name;
          pnode->interpolate_name_test2 = pnode->interpolate_name;
          if (option == 1)
            workspace.test2.insert
              (var_trans_pair(pnode->name_test2,
                              pnode->interpolate_name_test2));
          pnode->qdim2 = workspace.qdim(name);
          if (!(pnode->qdim2))
            ga_throw_error(expr, pnode->pos,
                           "Invalid null size of variable");
        }

        switch (prefix_id) {
        case 0: // value
          if (!test) {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_VAL; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_VAL; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_VAL; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_VAL; break;
            default: GMM_ASSERT1(false, "internal error");
            }
          } else {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_VAL_TEST; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_VAL_TEST; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_VAL_TEST; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_VAL_TEST; break;
            default: GMM_ASSERT1(false, "internal error");
            }
            if (q == 1 && mii.size() <= 1) {
              mii.resize(1);
              mii[0] = 2;
            } else
              mii.insert(mii.begin(), 2);
          }
          break;
        case 1: // grad
          if (!test) {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_GRAD; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_GRAD; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_GRAD; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_GRAD; break;
            default: GMM_ASSERT1(false, "internal error");
            }
            if (n > 1) {
              if (q == 1 && mii.size() == 1) mii[0] = n;
              else mii.push_back(n);
            }
          } else {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_GRAD_TEST; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_GRAD_TEST; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_GRAD_TEST; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_GRAD_TEST; break;
            default: GMM_ASSERT1(false, "internal error");
            }
            if (q == 1 && mii.size() <= 1) {
              mii.resize(1);
              mii[0] = 2;
            } else
              mii.insert(mii.begin(), 2);
            if (n > 1) mii.push_back(n);
          }
          break;
        case 2: // Hessian
          if (!test) {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_HESS; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_HESS; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_HESS; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_HESS; break;
            default: GMM_ASSERT1(false, "internal error");
            }
            if (n > 1) {
              if (q == 1 && mii.size() == 1) { mii[0] = n;  mii.push_back(n); }
              else { mii.push_back(n); mii.push_back(n); }
            }
          } else {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_HESS_TEST; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_HESS_TEST; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_HESS_TEST; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_HESS_TEST; break;
            default: GMM_ASSERT1(false, "internal error");
            }
            if (q == 1 && mii.size() <= 1) {
              mii.resize(1);
              mii[0] = 2;
            } else
              mii.insert(mii.begin(), 2);
            if (n > 1) { mii.push_back(n); mii.push_back(n); }
          }
          break;
        case 3: // divergence
          if (q != n)
            ga_throw_error(expr, pnode->pos,
                           "Divergence operator requires fem qdim ("
                           << q << ") to be equal to dim (" << n << ")");
          if (!test) {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_DIVERG; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_DIVERG; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_DIVERG; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_DIVERG; break;
            default: GMM_ASSERT1(false, "internal error");
            }
            mii.resize(1);
            mii[0] = 1;
          } else {
            switch (ndt) {
            case 1: pnode->node_type = GA_NODE_INTERPOLATE_DIVERG_TEST; break;
            case 2: pnode->node_type = GA_NODE_ELEMENTARY_DIVERG_TEST; break;
            case 3: pnode->node_type = GA_NODE_XFEM_PLUS_DIVERG_TEST; break;
            case 4: pnode->node_type = GA_NODE_XFEM_MINUS_DIVERG_TEST; break;
            default: GMM_ASSERT1(false, "internal error");
            }
            mii.resize(1);
            mii[0] = 2;
          }
          break;
        }
        pnode->t.adjust_sizes(mii);
        pnode->test_function_type = test;

        if (ndt == 1) {
          if (!(workspace.interpolate_transformation_exists
                (pnode->interpolate_name)))  {
            ga_throw_error(expr, pnode->pos,
                           "Unknown interpolate transformation");
          }
        } else if (ndt == 2) {
          if (!(workspace.elementary_transformation_exists
                (pnode->elementary_name))) {
            ga_throw_error(expr, pnode->pos,
                           "Unknown elementary transformation");
          }
        }
      }
      break;

    case GA_NODE_INTERPOLATE_FILTER:
      {
        if (pnode->children.size() == 2) {
          bool valid = (child1->node_type == GA_NODE_CONSTANT);
          int n = valid ? int(round(child1->t[0])) : -1;
          if (n < 0 || n > 100 || child1->tensor_order() > 0)
            ga_throw_error(expr, pnode->pos, "The third argument of "
                           "Interpolate_filter should be a (small) "
                           "non-negative integer.");
          pnode->nbc1 = size_type(n);
          tree.clear_node(child1);
        }
        if (!(workspace.interpolate_transformation_exists
              (pnode->interpolate_name)))
          ga_throw_error(expr, pnode->pos,
                         "Unknown interpolate transformation");
        pnode->t = child0->t;
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
        pnode->interpolate_name_test1 = child0->interpolate_name_test1;
        pnode->interpolate_name_test2 = child0->interpolate_name_test2;
        pnode->qdim1 = child0->qdim1;
        pnode->qdim2 = child0->qdim2;
      }
      break;


    case GA_NODE_OP:
      switch(pnode->op_type) {

      case GA_PLUS: case GA_MINUS:
        {
          if (pnode->op_type == GA_PLUS) pnode->symmetric_op = true;
          size_type c_size = std::min(size0.size(), size1.size());
          bool compatible = true;

          size_type f_ind = 0;
          if (child0->test_function_type &&
              child1->test_function_type == child0->test_function_type)
            f_ind = (child0->test_function_type == 3) ? 2:1;

          for (size_type i = f_ind; i < c_size; ++i)
            if (size0[i] != size1[i]) compatible = false;
          for (size_type i = c_size; i < size0.size(); ++i)
            if (size0[i] != 1) compatible = false;
          for (size_type i = c_size; i < size1.size(); ++i)
            if (size1[i] != 1) compatible = false;

          if (!compatible)
            ga_throw_error(expr, pnode->pos, "Addition or subtraction of "
                           "expressions of different sizes: "
                           << size0 << " != " << size1);
          if (child0->test_function_type || child1->test_function_type) {
            switch (option) {
            case 0: case 2:
              if (child0->name_test1.compare(child1->name_test1) ||
                  child0->name_test2.compare(child1->name_test2) ||
                  child0->interpolate_name_test1.compare
                  (child1->interpolate_name_test1) ||
                  child0->interpolate_name_test2.compare
                  (child1->interpolate_name_test2))
                compatible = false;
              break;
            case 1: break;
            default: GMM_ASSERT1(false, "Unknown option");
            }
          }

          if (child0->test_function_type != child1->test_function_type ||
              (!compatible && option != 2))
            ga_throw_error(expr, pnode->pos, "Addition or subtraction of "
                           "incompatible test functions");
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->test_function_type = 0;
            pnode->t = pnode->children[0]->t;
            if (pnode->op_type == GA_MINUS)
              pnode->t -= pnode->children[1]->t;
            else
              pnode->t += pnode->children[1]->t;
            tree.clear_children(pnode);
          } else {
            pnode->t = child0->t;
            pnode->test_function_type = child0->test_function_type;
            pnode->name_test1 = child0->name_test1;
            pnode->name_test2 = child0->name_test2;
            pnode->interpolate_name_test1 = child0->interpolate_name_test1;
            pnode->interpolate_name_test2 = child0->interpolate_name_test2;
            pnode->qdim1 = child0->qdim1;
            pnode->qdim2 = child0->qdim2;

            if (option == 2 && !compatible) {
              bool child0_compatible = true, child1_compatible = true;
              if (pnode->test_function_type & 1) {
                if (child0->name_test1.compare(workspace.selected_test1.first)
                    || child0->interpolate_name_test1.compare
                    (workspace.selected_test1.second))
                  child0_compatible = false;
                if (child1->name_test1.compare(workspace.selected_test1.first)
                    || child1->interpolate_name_test1.compare
                    (workspace.selected_test1.second))
                  child1_compatible = false;
              }
              if (pnode->test_function_type & 2) {
                if (child0->name_test2.compare(workspace.selected_test2.first)
                    || child0->interpolate_name_test2.compare
                    (workspace.selected_test2.second))
                  child0_compatible = false;
                if (child1->name_test2.compare(workspace.selected_test2.first)
                    || child1->interpolate_name_test1.compare
                    (workspace.selected_test2.second))
                  child1_compatible = false;
              }
              if (child0_compatible) {
                tree.replace_node_by_child(pnode, 0);
                pnode = child0;
              } else if (child1_compatible) {
                if (pnode->op_type == GA_MINUS) {
                  pnode->op_type = GA_UNARY_MINUS;
                  pnode->t = child1->t;
                  pnode->test_function_type = child1->test_function_type;
                  pnode->name_test1 = child1->name_test1;
                  pnode->name_test2 = child1->name_test2;
                  pnode->interpolate_name_test1=child1->interpolate_name_test1;
                  pnode->interpolate_name_test2=child1->interpolate_name_test2;
                  pnode->qdim1 = child1->qdim1;
                  pnode->qdim2 = child1->qdim2;
                  tree.clear_node(child0);
                } else {
                  tree.replace_node_by_child(pnode, 1);
                  pnode = child1;
                }
              }
            }

            // simplification if one of the two operands is constant and zero
            if (child0->tensor_is_zero()) {
              if (pnode->op_type == GA_MINUS) {
                pnode->op_type = GA_UNARY_MINUS;
                tree.clear_node(child0);
              } else {
                tree.replace_node_by_child(pnode, 1);
                pnode = child1;
              }
            } else if (child1->tensor_is_zero()) {
              tree.replace_node_by_child(pnode, 0);
              pnode = child0;
            }
          }
        }
        break;

      case GA_DOTMULT: case GA_DOTDIV:
        {
          if (pnode->op_type == GA_DOTMULT) pnode->symmetric_op = true;
          bool compatible = true;
          if (child0->tensor_proper_size() != child1->tensor_proper_size())
            compatible = false;

          if (child0->tensor_proper_size() != 1) {
            if (child0->tensor_order() != child1->tensor_order())
              compatible = false;

            for (size_type i = 0; i < child0->tensor_order(); ++i)
              if (child0->tensor_proper_size(i)!=child1->tensor_proper_size(i))
                compatible = false;
          }

          if (!compatible)
            ga_throw_error(expr, pnode->pos,
                           "Arguments of different sizes for .* or ./");

          if (pnode->op_type == GA_DOTDIV && child1->test_function_type)
            ga_throw_error(expr, pnode->pos,
                           "Division by test functions is not allowed");

          pnode->mult_test(child0, child1, expr);
          mi = pnode->t.sizes();
          for (size_type i = 0; i < child0->tensor_order(); ++i)
            mi.push_back(child0->tensor_proper_size(i));
          pnode->t.adjust_sizes(mi);

          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->test_function_type = 0;
            if (pnode->op_type == GA_DOTMULT) {
              for (size_type i = 0; i < child0->t.size(); ++i)
                pnode->t[i] = child0->t[i] * child1->t[i];
            } else {
              for (size_type i = 0; i < child0->t.size(); ++i) {
                if (child1->t[i] == scalar_type(0))
                  ga_throw_error(expr, pnode->pos, "Division by zero.");
                pnode->t[i] = child0->t[i] / child1->t[i];
              }
            }
            tree.clear_children(pnode);
          } else {
            if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
              gmm::clear(pnode->t.as_vector());
              pnode->node_type = GA_NODE_ZERO;
              tree.clear_children(pnode);
            }
            if (child1->tensor_is_zero() && pnode->op_type == GA_DOTDIV)
              ga_throw_error(expr, pnode->pos, "Division by zero.");
          }
        }
        break;

      case GA_UNARY_MINUS:
        pnode->t = child0->t;
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
        pnode->interpolate_name_test1 = child0->interpolate_name_test1;
        pnode->interpolate_name_test2 = child0->interpolate_name_test2;
        pnode->qdim1 = child0->qdim1;
        pnode->qdim2 = child0->qdim2;
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->test_function_type = 0;
          gmm::scale(pnode->t.as_vector(), scalar_type(-1));
          tree.clear_children(pnode);
        } else if (child0->node_type == GA_NODE_ZERO) {
          tree.replace_node_by_child(pnode, 0);
          pnode = child0;
        }
        break;

      case GA_QUOTE:
        if (dim0 > 2)
          ga_throw_error(expr, pnode->pos, "Transpose operator is for "
                         "vectors or matrices only.");
        mi = size0;
        if (child0->tensor_proper_size() == 1)
          { tree.replace_node_by_child(pnode, 0); pnode = child0; break; }
        else if (dim0 == 2) std::swap(mi.back(), mi[size0.size()-2]);
        else { size_type N = mi.back(); mi.back() = 1; mi.push_back(N); }

        pnode->t.adjust_sizes(mi);
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
        pnode->interpolate_name_test1 = child0->interpolate_name_test1;
        pnode->interpolate_name_test2 = child0->interpolate_name_test2;
        pnode->qdim1 = child0->qdim1;
        pnode->qdim2 = child0->qdim2;
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->test_function_type = 0;
          if (dim0 == 2) {
            for (size_type i = 0; i < mi.back(); ++i)
              for (size_type j = 0; j < mi[size0.size()-2]; ++j)
                pnode->t(j, i) = child0->t(i,j);
          } else if (dim0 == 1) {
            for (size_type i = 0; i < mi.back(); ++i)
              pnode->t(0, i) = child0->t[i];
          }
          tree.clear_children(pnode);
        } else if (child0->node_type == GA_NODE_ZERO) {
          pnode->node_type = GA_NODE_ZERO;
          gmm::clear(pnode->t.as_vector());
          tree.clear_children(pnode);
        }
        break;

      case GA_TRACE:
        {
          mi = size0;
          size_type N = (child0->tensor_proper_size() == 1) ? 1 : mi.back();

          if ((dim0 != 2 && child0->tensor_proper_size() != 1) ||
              (dim0 == 2 && mi[mi.size()-2] != N))
            ga_throw_error(expr, pnode->pos,
                           "Trace operator is for square matrices only.");

          if (dim0 == 2) { mi.pop_back(); mi.pop_back(); }
          pnode->t.adjust_sizes(mi);
          pnode->test_function_type = child0->test_function_type;
          pnode->name_test1 = child0->name_test1;
          pnode->name_test2 = child0->name_test2;
          pnode->interpolate_name_test1 = child0->interpolate_name_test1;
          pnode->interpolate_name_test2 = child0->interpolate_name_test2;
          pnode->qdim1 = child0->qdim1;
          pnode->qdim2 = child0->qdim2;
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->test_function_type = 0;
            if (dim0 == 2) {
              pnode->t[0] = scalar_type(0);
              for (size_type i = 0; i < N; ++i)
                pnode->t[0] += child0->t(i,i);
            } else {
              pnode->t[0] += child0->t[0];
            }
            tree.clear_children(pnode);
          } else if (child0->node_type == GA_NODE_ZERO) {
            pnode->node_type = GA_NODE_ZERO;
            gmm::clear(pnode->t.as_vector());
            tree.clear_children(pnode);
          }
        }
        break;

      case GA_DEVIATOR:
        {
          mi = size0;
          size_type N = (child0->tensor_proper_size() == 1) ? 1 : mi.back();

          if ((dim0 != 2 && child0->tensor_proper_size() != 1) ||
              (dim0 == 2 && mi[mi.size()-2] != N))
            ga_throw_error(expr, pnode->pos,
                           "Deviator operator is for square matrices only.");

          pnode->t.adjust_sizes(mi);
          pnode->test_function_type = child0->test_function_type;
          pnode->name_test1 = child0->name_test1;
          pnode->name_test2 = child0->name_test2;
          pnode->interpolate_name_test1 = child0->interpolate_name_test1;
          pnode->interpolate_name_test2 = child0->interpolate_name_test2;
          pnode->qdim1 = child0->qdim1;
          pnode->qdim2 = child0->qdim2;
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->test_function_type = 0;
            if (dim0 == 2) {
              scalar_type tr(0);
              gmm::copy(child0->t.as_vector(), pnode->t.as_vector());
              for (size_type i = 0; i < N; ++i)
                tr += child0->t(i,i);
              for (size_type i = 0; i < N; ++i)
                pnode->t(i,i) -= tr / scalar_type(N);
            } else {
              pnode->t[0] = scalar_type(0);
            }
            tree.clear_children(pnode);
          } else if (child0->node_type == GA_NODE_ZERO) {
            pnode->node_type = GA_NODE_ZERO;
            gmm::clear(pnode->t.as_vector());
            tree.clear_children(pnode);
          }
        }
        break;

      case GA_PRINT:
        {
          pnode->t = child0->t;
          pnode->test_function_type = child0->test_function_type;
          pnode->name_test1 = child0->name_test1;
          pnode->name_test2 = child0->name_test2;
          pnode->interpolate_name_test1 = child0->interpolate_name_test1;
          pnode->interpolate_name_test2 = child0->interpolate_name_test2;
          pnode->qdim1 = child0->qdim1;
          pnode->qdim2 = child0->qdim2;
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            cout << "Print constant term "; ga_print_node(child0, cout);
            cout << ": " << pnode->t << endl;
            tree.clear_children(pnode);
          } else if (child0->node_type == GA_NODE_ZERO) {
            pnode->node_type = GA_NODE_ZERO;
            gmm::clear(pnode->t.as_vector());
            cout << "Print zero term "; ga_print_node(child0, cout);
            cout << ": " << pnode->t << endl;
            tree.clear_children(pnode);
          }
        }
        break;

      case GA_DOT:
        if (dim1 > 1)
          ga_throw_error(expr, pnode->pos, "The second argument of the dot "
                         "product has to be a vector.")
        else {
          size_type s0 = dim0 == 0 ? 1 : size0.back();
          size_type s1 = dim1 == 0 ? 1 : size1.back();
          if (s0 != s1) ga_throw_error(expr, pnode->pos, "Dot product "
                                       "of expressions of different sizes ("
                                       << s0 << " != " << s1 << ").");
          if (child0->tensor_order() <= 1) pnode->symmetric_op = true;
          pnode->mult_test(child0, child1, expr);
          if (dim0 > 1) {
            mi = pnode->t.sizes();
            for (size_type i = 1; i < dim0; ++i)
              mi.push_back(child0->tensor_proper_size(i-1));
            pnode->t.adjust_sizes(mi);
          }

          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            gmm::clear(pnode->t.as_vector());
            size_type k = 0;
            for (size_type i = 0, j = 0; i < child0->t.size(); ++i) {
             pnode->t[j] += child0->t[i] * child1->t[k];
             ++j; if (j == pnode->t.size()) { j = 0; ++k; }
            }
            GMM_ASSERT1(k == child1->t.size(), "Internal error");
            tree.clear_children(pnode);
          } else {
            if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
              gmm::clear(pnode->t.as_vector());
              pnode->node_type = GA_NODE_ZERO;
              tree.clear_children(pnode);
            }
          }
        }
        break;

      case GA_COLON:
        if (dim1 > 2)
          ga_throw_error(expr, pnode->pos,
                          "Frobenius product acts only on matrices.")
        else {
          size_type s00 = (dim0 == 0) ? 1
            : (dim0 == 1 ? size0.back() : size0[size0.size()-2]);
          size_type s01 = (dim0 >= 2) ? size0.back() : 1;
          size_type s10 = (dim1 == 0) ? 1
            : (dim1 == 1 ? size1.back() : size1[size1.size()-2]);
          size_type s11 = (dim1 >= 2) ? size1.back() : 1;
          if (s00 != s10 || s01 != s11)
            ga_throw_error(expr, pnode->pos, "Frobenius product "
                           "of expressions of different sizes ("
                           << s00 << "," << s01 << " != " << s10 << ","
                           << s11 << ").");
          if (child0->tensor_order() <= 2) pnode->symmetric_op = true;
          pnode->mult_test(child0, child1, expr);
          if (dim0 > 2) {
            mi = pnode->t.sizes();
            for (size_type i = 2; i < dim0; ++i)
              mi.push_back(child0->tensor_proper_size(i-2));
            pnode->t.adjust_sizes(mi);
          }
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            gmm::clear(pnode->t.as_vector());
            size_type k = 0;
            for (size_type i = 0, j = 0; i < child0->t.size(); ++i) {
             pnode->t[j] += child0->t[i] * child1->t[k];
             ++j; if (j == pnode->t.size()) { j = 0; ++k; }
            }
            GMM_ASSERT1(k == child1->t.size(), "Internal error");
            tree.clear_children(pnode);
          } else {
            if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
              gmm::clear(pnode->t.as_vector());
              pnode->node_type = GA_NODE_ZERO;
              tree.clear_children(pnode);
            }
          }
        }
        break;

      case GA_TMULT:
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->test_function_type = 0;
          if (child0->t.size() == 1 && child1->t.size() == 1) {
            pnode->init_scalar_tensor
              (child0->t[0] * child1->t[0]);
          } else if (child0->t.size() == 1) {
            pnode->t = child1->t;
            gmm::scale(pnode->t.as_vector(), scalar_type(child0->t[0]));
          } else if (child1->t.size() == 1) {
            pnode->t = child0->t;
            gmm::scale(pnode->t.as_vector(), scalar_type(child1->t[0]));
          } else {
            if (dim0+dim1 > 6)
              ga_throw_error(expr, pnode->pos, "Unauthorized "
                              "tensor multiplication.");
            for (size_type i = 0; i < dim0; ++i)
              mi.push_back(child0->t.size(i));
            for (size_type i = 0; i < dim1; ++i)
              mi.push_back(child1->t.size(i));
            pnode->t.adjust_sizes(mi);
            size_type n0 = child0->t.size();
            size_type n1 = child1->t.size();
            for (size_type i = 0; i < n0; ++i)
              for (size_type j = 0; j < n1; ++j)
                pnode->t[i+j*n0] = child0->t[i] * child1->t[j];
          }
          tree.clear_children(pnode);
        } else {
          pnode->mult_test(child0, child1, expr);
          mi = pnode->t.sizes();
          if (child0->tensor_proper_size() != 1
              || child1->tensor_proper_size() != 1) {
            if (child0->tensor_proper_size() == 1) {
              for (size_type i = 0; i < dim1; ++i)
                mi.push_back(child1->tensor_proper_size(i));
            } else if (child1->t.size() == 1) {
              for (size_type i = 0; i < dim0; ++i)
                mi.push_back(child0->tensor_proper_size(i));
            } else {
              if (dim0+dim1 > 6)
                ga_throw_error(expr, pnode->pos, "Unauthorized "
                                "tensor multiplication.");
              for (size_type i = 0; i < dim0; ++i)
                mi.push_back(child0->tensor_proper_size(i));
              for (size_type i = 0; i < dim1; ++i)
                mi.push_back(child1->tensor_proper_size(i));
            }
            pnode->t.adjust_sizes(mi);
          }
          if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
            gmm::clear(pnode->t.as_vector());
            pnode->node_type = GA_NODE_ZERO;
            tree.clear_children(pnode);
          }
        }
        break;

      case GA_MULT:
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->test_function_type = 0;
          if (child0->tensor_proper_size() == 1 &&
              child1->tensor_proper_size() == 1) {
            pnode->init_scalar_tensor(child0->t[0]*child1->t[0]);
          } else if (child0->tensor_proper_size() == 1) {
            pnode->t = child1->t;
            gmm::scale(pnode->t.as_vector(), child0->t[0]);
          } else if (child1->tensor_proper_size() == 1) {
            pnode->t = child0->t;
            gmm::scale(pnode->t.as_vector(), child1->t[0]);
          } else if (dim0 == 2 && dim1 == 1) {
            size_type m = child0->t.size(0), n = child0->t.size(1);
            if (n != child1->t.size(0))
              ga_throw_error(expr, pnode->pos,
                             "Incompatible sizes in matrix-vector "
                             "multiplication (" << n << " != "
                             << child1->t.size(0) << ").");
            pnode->init_vector_tensor(m);
            gmm::clear(pnode->t.as_vector());
            for (size_type i = 0; i < m; ++i)
              for (size_type j = 0; j < n; ++j)
                pnode->t[i] += child0->t(i,j) * child1->t[j];
          } else if (dim0 == 2 && dim1 == 2) {
            size_type m = child0->t.size(0);
            size_type n = child0->t.size(1);
            size_type p = child1->t.size(1);
            if (n != child1->t.size(0))
              ga_throw_error(expr, pnode->pos,
                             "Incompatible sizes in matrix-matrix "
                             "multiplication (" << n << " != "
                             << child1->t.size(0) << ").");
            pnode->init_matrix_tensor(m,p);
            gmm::clear(pnode->t.as_vector());
            for (size_type i = 0; i < m; ++i)
              for (size_type j = 0; j < n; ++j)
                for (size_type k = 0; k < p; ++k)
                  pnode->t(i,k) += child0->t(i,j) * child1->t(j,k);
          }
          else if (dim0 == 4 && dim1 == 2) {
            size_type m = child0->t.size(0), n = child0->t.size(1);
            size_type o = child0->t.size(2), p = child0->t.size(3);
            if (o != child1->t.size(0) || p != child1->t.size(1))
              ga_throw_error(expr, pnode->pos,
                             "Incompatible sizes in tensor-matrix "
                             "multiplication (" << o << "," << p << " != "
                             << child1->t.size(0) << "," << child1->t.size(1)
                             << ").");
            pnode->init_matrix_tensor(m,n);
            gmm::clear(pnode->t.as_vector());
            for (size_type i = 0; i < m; ++i)
              for (size_type j = 0; j < n; ++j)
                for (size_type k = 0; k < o; ++k)
                  for (size_type l = 0; l < p; ++l)
                    pnode->t(i,j) += child0->t(i,j,k,l) * child1->t(k,l);
          } else ga_throw_error(expr, pnode->pos,
                                 "Unauthorized multiplication.");
          tree.clear_children(pnode);
        } else {
          pnode->mult_test(child0, child1, expr);
          mi = pnode->t.sizes();

          if (child0->tensor_proper_size() == 1 &&
              child1->tensor_proper_size() == 1) {
            pnode->symmetric_op = true;
          } else if (child0->tensor_proper_size() == 1) {
            pnode->symmetric_op = true;
            for (size_type i = 0; i < dim1; ++i)
                mi.push_back(child1->tensor_proper_size(i));
          } else if (child1->tensor_proper_size() == 1) {
            pnode->symmetric_op = true;
            for (size_type i = 0; i < dim0; ++i)
                mi.push_back(child0->tensor_proper_size(i));
          } else if (child0->tensor_order() == 2 &&
                     child1->tensor_order() == 1) {
            size_type m = child0->tensor_proper_size(0);
            size_type n = child0->tensor_proper_size(1);
            mi.push_back(m);
            if (n != child1->tensor_proper_size(0))
              ga_throw_error(expr, pnode->pos,
                             "Incompatible sizes in matrix-vector "
                             "multiplication (" << n << " != "
                             << child1->tensor_proper_size(0) << ").");
          } else if (child0->tensor_order() == 2 &&
                     child1->tensor_order() == 2) {
            size_type m = child0->tensor_proper_size(0);
            size_type n = child0->tensor_proper_size(1);
            size_type p = child1->tensor_proper_size(1);
            mi.push_back(m); mi.push_back(p);
            if (n != child1->tensor_proper_size(0))
              ga_throw_error(expr, pnode->pos,
                             "Incompatible sizes in matrix-matrix "
                             "multiplication (" << n << " != "
                             << child1->tensor_proper_size(0) << ").");
          }
          else if (pnode->children[0]->tensor_order() == 4 &&
                   pnode->children[1]->tensor_order() == 2) {
            size_type m = child0->tensor_proper_size(0);
            size_type n = child0->tensor_proper_size(1);
            size_type o = child0->tensor_proper_size(2);
            size_type p = child0->tensor_proper_size(3);
            mi.push_back(m); mi.push_back(n);
            if (o != child1->tensor_proper_size(0) ||
                p != child1->tensor_proper_size(1))
              ga_throw_error(expr, pnode->pos,
                             "Incompatible sizes in tensor-matrix "
                             "multiplication (" << o << "," << p << " != "
                             << child1->tensor_proper_size(0) << ","
                             << child1->tensor_proper_size(1) << ").");
          } else ga_throw_error(expr, pnode->pos,
                                "Unauthorized multiplication.");
          pnode->t.adjust_sizes(mi);
          // Simplifications
          if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
            gmm::clear(pnode->t.as_vector());
            pnode->node_type = GA_NODE_ZERO;
            tree.clear_children(pnode);
          } else if (child0->node_type == GA_NODE_CONSTANT &&
                     child0->t.size() == 1 && child0->t[0] == scalar_type(1)) {
            tree.replace_node_by_child(pnode, 1);
            pnode = child1;
          } else if (child1->node_type == GA_NODE_CONSTANT &&
                     child1->t.size() == 1 && child1->t[0] == scalar_type(1)) {
            tree.replace_node_by_child(pnode, 0);
            pnode = child0;
          }
        }
        break;

      case GA_DIV:
        if (child1->tensor_order() > 0)
          ga_throw_error(expr, pnode->pos,
                         "Only the division by a scalar is allowed.");
        if (child1->test_function_type)
          ga_throw_error(expr, pnode->pos,
                         "Division by test functions is not allowed.");
        if (child1->node_type == GA_NODE_CONSTANT &&
            child1->t[0] == scalar_type(0))
          ga_throw_error(expr, pnode->children[1]->pos, "Division by zero");

        pnode->t = child0->t;
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
        pnode->interpolate_name_test1 = child0->interpolate_name_test1;
        pnode->interpolate_name_test2 = child0->interpolate_name_test2;
        pnode->qdim1 = child0->qdim1;
        pnode->qdim2 = child0->qdim2;

        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->t = pnode->children[0]->t;
          pnode->test_function_type = 0;
          gmm::scale(pnode->t.as_vector(),
                     scalar_type(1) / pnode->children[1]->t[0]);
          tree.clear_children(pnode);
        } else if (child0->tensor_is_zero()) {
          gmm::clear(pnode->t.as_vector());
          pnode->node_type = GA_NODE_ZERO;
          tree.clear_children(pnode);
        } else if (child1->node_type == GA_NODE_CONSTANT &&
                   child1->t.size() == 1 && child1->t[0] == scalar_type(1)) {
          tree.replace_node_by_child(pnode, 0);
          pnode = child0;
        }
        break;

      default:GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX:
      {
        if (!all_sc) {
          ga_throw_error(expr, pnode->pos, "Constant vector/matrix/tensor "
                         "components should be scalar valued.");
        }

        size_type nbc1 = pnode->nbc1, nbc2 = pnode->nbc2, nbc3 = pnode->nbc3;
        size_type nbl = pnode->children.size() / (nbc1*nbc2*nbc3);
        if (all_cte) pnode->node_type = GA_NODE_CONSTANT;
        pnode->test_function_type = 0;
        for (size_type i = 0; i < pnode->children.size(); ++i) {
          if (pnode->children[i]->test_function_type) {
            if (pnode->test_function_type == 0) {
              pnode->test_function_type=pnode->children[i]->test_function_type;
              pnode->name_test1 = pnode->children[i]->name_test1;
              pnode->name_test2 = pnode->children[i]->name_test2;
              pnode->interpolate_name_test1
                = pnode->children[i]->interpolate_name_test1;
              pnode->interpolate_name_test2
                = pnode->children[i]->interpolate_name_test2;
              pnode->qdim1 = pnode->children[i]->qdim1;
              pnode->qdim2 = pnode->children[i]->qdim2;
            } else {
              if (pnode->test_function_type !=
                  pnode->children[i]->test_function_type ||
                  pnode->name_test1.compare(pnode->children[i]->name_test1) ||
                  pnode->name_test2.compare(pnode->children[i]->name_test2) ||
                  pnode->interpolate_name_test1.compare
                  (pnode->children[i]->interpolate_name_test1) ||
                  pnode->interpolate_name_test2.compare
                  (pnode->children[i]->interpolate_name_test2))
                ga_throw_error(expr, pnode->pos, "Inconsistent use of test "
                               "function in constant matrix.");
            }
          }
        }
        mi.resize(0);
        if (pnode->test_function_type) mi.push_back(2);
        if (pnode->test_function_type >= 3) mi.push_back(2);
        if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1 && nbl == 1) {
          pnode->t.adjust_sizes(mi);
          if (all_cte) pnode->t[0] = child0->t[0];
        } else {
          mi.push_back(nbl);
          if (nbc3 != 1) mi.push_back(nbc3);
          if (nbc2 != 1) mi.push_back(nbc2);
          if (nbc1 != 1) mi.push_back(nbc1);
          pnode->t.adjust_sizes(mi);
          if (all_cte) {
            size_type n = 0;
            if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1)
              for (size_type i = 0; i < nbl; ++i)
                pnode->t[i] = pnode->children[i]->t[0];
            else if (nbc2 == 1 && nbc3 == 1) // TODO: verify order
              for (size_type i = 0; i < nbl; ++i)
                for (size_type j = 0; j < nbc1; ++j)
                  pnode->t(i,j) = pnode->children[n++]->t[0];
            else if (nbc3 == 1) // TODO: verify order
              for (size_type i = 0; i < nbl; ++i)
                for (size_type j = 0; j < nbc2; ++j)
                  for (size_type k = 0; k < nbc1; ++k)
                    pnode->t(i,j,k) = pnode->children[n++]->t[0];
            else // TODO: verify order
              for (size_type i = 0; i < nbl; ++i)
                for (size_type j = 0; j < nbc3; ++j)
                  for (size_type k = 0; k < nbc2; ++k)
                    for (size_type l = 0; l < nbc1; ++l)
                      pnode->t(i,j,k,l) = pnode->children[n++]->t[0];
          }
        }
        if (all_cte) tree.clear_children(pnode);
      }
      break;


    case GA_NODE_NAME:
      {
        std::string name = pnode->name;

        if (!ignore_X && !(name.compare("X"))) {
          pnode->node_type = GA_NODE_X;
          pnode->nbc1 = 0;
          pnode->init_vector_tensor(meshdim);
          break;
        }
        if (!(name.compare("element_size"))) {
          pnode->node_type = GA_NODE_ELT_SIZE;
          pnode->init_scalar_tensor(0);
          break;
        }
        if (!(name.compare("element_K"))) {
          pnode->node_type = GA_NODE_ELT_K;
          pnode->init_matrix_tensor(meshdim, ref_elt_dim);
          break;
        }
        if (!(name.compare("element_B"))) {
          pnode->node_type = GA_NODE_ELT_B;
          pnode->init_matrix_tensor(ref_elt_dim, meshdim);
          break;
        }
        if (!(name.compare("Normal"))) {
          pnode->node_type = GA_NODE_NORMAL;
          pnode->init_vector_tensor(meshdim);
          break;
        }
        if (!(name.compare("Reshape"))) {
          pnode->node_type = GA_NODE_RESHAPE;
          pnode->init_vector_tensor(meshdim);
          break;
        }

        if (name.compare(0, 11, "Derivative_") == 0) {
          name = name.substr(11);
          bool valid = true;
          pnode->der1 = 1; pnode->der2 = 0;
          char *p;
          size_type d = strtol(name.c_str(), &p, 10);
          size_type s = p - name.c_str();
          if (s > 0) {
            pnode->der1 = d;
            if (name[s] != '_') valid = false; else
              name = name.substr(s+1);
          }
          d = strtol(name.c_str(), &p, 10);
          s = p - name.c_str();
          if (s > 0) {
            pnode->der2 = d;
            if (name[s] != '_') valid = false; else
              name = name.substr(s+1);
          }
          if (!valid || pnode->der1 == 0)
            ga_throw_error(expr, pnode->pos, "Invalid derivative format");
        }

        ga_predef_operator_tab &PREDEF_OPERATORS
          = dal::singleton<ga_predef_operator_tab>::instance(0);
        ga_predef_function_tab::const_iterator it=PREDEF_FUNCTIONS.find(name);
        if (it != PREDEF_FUNCTIONS.end()) {
          // Predefined function found
          pnode->node_type = GA_NODE_PREDEF_FUNC;
          pnode->name = name;
          pnode->test_function_type = 0;
          if (pnode->der1) {
            if (pnode->der1 > it->second.nbargs()
                || pnode->der2 > it->second.nbargs())
              ga_throw_error(expr, pnode->pos, "Invalid derivative.");
            const ga_predef_function &F = it->second;
            if (F.ftype() == 0 && !(pnode->der2)) {
              pnode->name = ((pnode->der1 == 1) ?
                             F.derivative1() : F.derivative2());
              pnode->der1 = pnode->der2 = 0;
            }
          }
        } else if (SPEC_FUNCTIONS.find(name) != SPEC_FUNCTIONS.end()) {
          // Special function found
          if (pnode->der1)
            ga_throw_error(expr, pnode->pos, "Special functions do not "
                           "support derivatives.");
          pnode->node_type = GA_NODE_SPEC_FUNC;
          pnode->name = name;
          pnode->test_function_type = 0;
          if (!name.compare("pi")) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->init_scalar_tensor(M_PI);
          } else if (!name.compare("meshdim")) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->init_scalar_tensor(scalar_type(meshdim));
          }
        } else if (PREDEF_OPERATORS.tab.find(name)
                   != PREDEF_OPERATORS.tab.end()) {
          // Nonlinear operator found
          pnode->node_type = GA_NODE_OPERATOR;
          pnode->name = name;
          pnode->test_function_type = 0;
        } else if (workspace.macro_exists(name)) {
          GMM_ASSERT1(pnode->der1 == 0 && pnode->der2 == 0,
                      "Derivativation of a macro is not allowed");
          pga_tree_node parent = pnode->parent;
          size_type ind_in_parent = size_type(-1);
          if (parent) {
            for (size_type i = 0; i < parent->children.size(); ++i)
              if (parent->children[i] == pnode)
                ind_in_parent = i;
            GMM_ASSERT1(ind_in_parent != size_type(-1), "Internal error");
          }
          ga_tree &ma_tree
            = workspace.macro_tree(name, meshdim, ref_elt_dim, ignore_X);
          pga_tree_node &newnode = (ind_in_parent == size_type(-1))
            ? tree.root : pnode->parent->children[ind_in_parent];
          tree.copy_node(ma_tree.root, pnode->parent, newnode);
          delete pnode;
          pnode = newnode;
        } else {
          // Search for a variable name with optional gradient, Hessian,
          // divergence or test functions

          size_type prefix_id = ga_parse_prefix_operator(name);
          size_type test = ga_parse_prefix_test(name);

          if (!(workspace.variable_exists(name)))
            ga_throw_error(expr, pnode->pos, "Unknown variable, function, "
                           "operator or data " + name);

          if (pnode->der1)
            ga_throw_error(expr, pnode->pos, "Derivative is for functions or "
                           "operators, not for variables. Use Grad instead.");
          pnode->name = name;

          const mesh_fem *mf = workspace.associated_mf(name);
          const im_data *imd = workspace.associated_im_data(name);

          if (test && workspace.is_constant(name) &&
              !(workspace.is_disabled_variable(name)))
            ga_throw_error(expr, pnode->pos, "Test functions of constants "
                           "are not allowed.");
          if (test == 1) {
            pnode->name_test1 = name;
            pnode->interpolate_name_test1 = "";
            if (option == 1)
              workspace.test1.insert
                (var_trans_pair(pnode->name_test1,
                                pnode->interpolate_name_test1));
            pnode->qdim1 = mf ? workspace.qdim(name)
                              : gmm::vect_size(workspace.value(name));
            if (!(pnode->qdim1))
              ga_throw_error(expr, pnode->pos,
                             "Invalid null size of variable");
          } else if (test == 2) {
            pnode->name_test2 = name;
            pnode->interpolate_name_test2 = "";
            if (option == 1)
              workspace.test2.insert
                (var_trans_pair(pnode->name_test2,
                                pnode->interpolate_name_test2));
            pnode->qdim2 = mf ? workspace.qdim(name)
                              : gmm::vect_size(workspace.value(name));
            if (!(pnode->qdim2))
              ga_throw_error(expr, pnode->pos,
                             "Invalid null size of variable");
          }

          if (!mf && (test || !imd)) {
            if (prefix_id)
              ga_throw_error(expr, pnode->pos, "Gradient, Hessian or Divergence"
                             " cannot be evaluated for fixed size data.");
            if (test)
              pnode->node_type = GA_NODE_VAL_TEST;
            else if (eval_fixed_size)
              pnode->node_type = GA_NODE_CONSTANT;
            else
              pnode->node_type = GA_NODE_VAL;

            size_type n = gmm::vect_size(workspace.value(name));
            if (n == 1) {
              if (test) {
                pnode->init_vector_tensor(1);
                pnode->t[0]=scalar_type(1);
              }
              else pnode->init_scalar_tensor(workspace.value(name)[0]);
            } else {
              if (test) {
                pnode->init_matrix_tensor(n,n);
                for (size_type i = 0; i < n; ++i)
                  for (size_type j = 0; j < n; ++j)
                    pnode->t(i,j) = (i == j) ? scalar_type(1) : scalar_type(0);
              } else {
                pnode->init_vector_tensor(n);
                gmm::copy(workspace.value(name), pnode->t.as_vector());
              }
            }
          } else if (!test && imd) {
            if (prefix_id)
              ga_throw_error(expr, pnode->pos, "Gradient, Hessian or Divergence"
                              " cannot be evaluated for im data.");
            pnode->node_type = GA_NODE_VAL;
            pnode->t.adjust_sizes(workspace.qdims(name));
          } else {
            size_type q = workspace.qdim(name);
            size_type n = mf->linked_mesh().dim();
            bgeot::multi_index mii = workspace.qdims(name);

            if (!q) ga_throw_error(expr, pnode->pos,
                                   "Invalid null size of variable " << name);
            if (mii.size() > 6)
              ga_throw_error(expr, pnode->pos,
                            "Tensor with too much dimensions. Limited to 6");

            switch (prefix_id) {
            case 0: // value
              pnode->node_type = test ? GA_NODE_VAL_TEST : GA_NODE_VAL;
              // For Test nodes a first dimension of size equal to 2 has to be
              // prepended by convention (to be adapted later)
              if (test && q == 1 && mii.size() <= 1) {
                mii.resize(1);
                mii[0] = 2;
              } else if (test) {
                mii.insert(mii.begin(), 2);
                pnode->t.adjust_sizes(mii);
              }
              break;
            case 1: // grad
              pnode->node_type = test ? GA_NODE_GRAD_TEST : GA_NODE_GRAD;
              if (test) {
                if (q == 1 && mii.size() <= 1) {
                  mii.resize(1);
                  mii[0] = 2;
                } else
                  mii.insert(mii.begin(), 2);
              }
              if (n > 1) {
                if (mii.size() == 1 && mii[0] == 1) mii[0] = n;
                else mii.push_back(n);
              }
              break;
            case 2: // Hessian
              pnode->node_type = test ? GA_NODE_HESS_TEST : GA_NODE_HESS;
              if (test) {
                if (q == 1 && mii.size() <= 1) {
                  mii.resize(1);
                  mii[0] = 2;
                } else
                  mii.insert(mii.begin(), 2);
              }
              if (n > 1) {
                if (mii.size() == 1 && mii[0] == 1) mii[0] = n;
                else mii.push_back(n);
                mii.push_back(n);
              }
              break;
            case 3: // divergence
              pnode->node_type = test ? GA_NODE_DIVERG_TEST : GA_NODE_DIVERG;
              if (q != n)
                ga_throw_error(expr, pnode->pos,
                               "Divergence operator can only be applied to"
                               "Fields with qdim (" << q << ") equal to dim ("
                               << n << ")");
              mii.resize(1);
              mii[0] = test ? 2 : 1;
              break;
            }
            pnode->t.adjust_sizes(mii);
          }
          pnode->test_function_type = test;
        }
      }
      break;

    case GA_NODE_PARAMS:
      if (child0->node_type == GA_NODE_X) {
        child0->init_scalar_tensor(0);
        if (pnode->children.size() != 2)
          ga_throw_error(expr, child1->pos, "X stands for the coordinates on "
                         "the real elements. It accepts only one index.");
        if (!(child1->node_type == GA_NODE_CONSTANT) || child1->t.size() != 1)
          ga_throw_error(expr, child1->pos, "Index for X has to be constant "
                         "and of size 1.");
        child0->nbc1 = size_type(round(child1->t[0]));
        if (child0->nbc1 == 0 || child0->nbc1 > meshdim)
          ga_throw_error(expr, child1->pos, "Index for X not convenient. "
                         "Found " << child0->nbc1 << " with meshdim = "
                         << meshdim);
        tree.replace_node_by_child(pnode, 0);
        pnode = child0;

      } else if (child0->node_type == GA_NODE_RESHAPE) {
        if (pnode->children.size() < 3)
          ga_throw_error(expr, child1->pos,
                         "Not enough parameters for Reshape");
        if (pnode->children.size() > 8)
          ga_throw_error(expr, child1->pos,
                         "Too many parameters for Reshape");
        pnode->t = child1->t;
        pnode->test_function_type = child1->test_function_type;
        pnode->name_test1 = child1->name_test1;
        pnode->name_test2 = child1->name_test2;
        pnode->interpolate_name_test1 = child1->interpolate_name_test1;
        pnode->interpolate_name_test2 = child1->interpolate_name_test2;
        pnode->qdim1 = child1->qdim1;
        pnode->qdim2 = child1->qdim2;
        mi.resize(0);
        for (size_type i = 0; i < pnode->nb_test_functions(); ++i)
          mi.push_back(size1[i]);

        for (size_type i = 2; i < pnode->children.size(); ++i) {
          if (pnode->children[i]->node_type != GA_NODE_CONSTANT)
            ga_throw_error(expr, pnode->children[i]->pos, "Reshape sizes "
                           "should be constant positive integers.");
          mi.push_back(size_type(round(pnode->children[i]->t[0])));
          if (mi.back() == 0)
            ga_throw_error(expr, pnode->children[i]->pos, "Wrong zero size "
                           "for Reshape.");
        }
        size_type total_size(1);
        for (size_type i = 0; i < mi.size(); ++i)
          total_size *= mi[i];
        if (total_size != pnode->t.size())
           ga_throw_error(expr, pnode->pos, "Invalid sizes for reshape.");
        pnode->t.adjust_sizes(mi);

        if (child1->node_type == GA_NODE_CONSTANT) {
          pnode->node_type = GA_NODE_CONSTANT;
          tree.clear_children(pnode);
        } else if (child1->node_type == GA_NODE_ZERO) {
          pnode->node_type = GA_NODE_ZERO;
          tree.clear_children(pnode);
        }
      } else if (child0->node_type == GA_NODE_PREDEF_FUNC) {

        // Evaluation of a predefined function

        for (size_type i = 1; i < pnode->children.size(); ++i)
          ga_valid_operand(expr, pnode->children[i]);
        std::string name = child0->name;
        ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
        const ga_predef_function &F = it->second;
        size_type nbargs = F.nbargs();
        if (nbargs+1 != pnode->children.size()) {
            ga_throw_error(expr, pnode->pos, "Bad number of arguments for "
                "predefined function " << name << ". Found "
                 << pnode->children.size()-1 << ", should be "<<nbargs << ".");
        }
        pnode->test_function_type = 0;
        pga_tree_node child2 = (nbargs == 2) ? pnode->children[2] : child1;
        all_cte = child1->node_type == GA_NODE_CONSTANT;
        if (nbargs == 2)
          all_cte = all_cte && (child2->node_type == GA_NODE_CONSTANT);
        if (child1->test_function_type || child2->test_function_type)
          ga_throw_error(expr, pnode->pos, "Test functions cannot be passed "
                         "as argument of a predefined function.");
        if (child1->tensor_order() > 2 || child2->tensor_order() > 2)
          ga_throw_error(expr, pnode->pos, "Sorry, function can be applied "
                         "to scalar, vector and matrices only.");
        size_type s1 = child1->t.size();
        size_type s2 = (nbargs == 2) ? child2->t.size() : s1;
        if (s1 != s2 && (s1 != 1 || s2 != 1))
          ga_throw_error(expr, pnode->pos,
                         "Invalid argument size for a scalar function.");

        if (nbargs == 1) {
          pnode->t = child1->t;
        } else {
          if (s1 == s2) {
            pnode->t = child1->t;
          } else if (s1 == 1) {
            pnode->t = child2->t;
          } else {
            pnode->t = child1->t;
          }
        }

        if (all_cte) {
          if (pnode->der1)
            GMM_ASSERT1(false, "Sorry, to be done");
          pnode->node_type = GA_NODE_CONSTANT;
          if (nbargs == 1) {
            for (size_type i = 0; i < s1; ++i)
              pnode->t[i] = F(child1->t[i]);
          } else {
            if (s1 == s2) {
              for (size_type i = 0; i < s1; ++i)
                pnode->t[i] = F(child1->t[i], child2->t[i]);
            } else if (s1 == 1) {
              for (size_type i = 0; i < s2; ++i)
                pnode->t[i] = F(child1->t[0], child2->t[i]);
            } else {
              for (size_type i = 0; i < s1; ++i)
                pnode->t[i] = F(child1->t[i], child2->t[0]);
            }
          }
          tree.clear_children(pnode);
        }
      } else if (child0->node_type == GA_NODE_SPEC_FUNC) {

        // Special constant functions: meshdim, qdim(u) ...

        for (size_type i = 1; i < pnode->children.size(); ++i)
          ga_valid_operand(expr, pnode->children[i]);
        if (pnode->children.size() != 2)
          ga_throw_error(expr, pnode->pos,
                         "One and only one argument is allowed for function "
                         +child0->name+".");

        if (!(child0->name.compare("qdim"))) {
          if (child1->node_type != GA_NODE_VAL)
            ga_throw_error(expr, pnode->pos, "The argument of qdim "
                           "function can only be a variable name.");
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->init_scalar_tensor(scalar_type(workspace.qdim(child1->name)));
          if (pnode->t[0] <= 0)
            ga_throw_error(expr, pnode->pos,
                           "Invalid null size of variable");
        } else if (!(child0->name.compare("qdims"))) {
          if (child1->node_type != GA_NODE_VAL)
            ga_throw_error(expr, pnode->pos, "The argument of qdim "
                           "function can only be a variable name.");
          pnode->node_type = GA_NODE_CONSTANT;
          bgeot::multi_index mii = workspace.qdims(child1->name);
          if (mii.size() > 6)
            ga_throw_error(expr, pnode->pos,
                           "Tensor with too much dimensions. Limited to 6");
          if (mii.size() == 0 || scalar_type(mii[0]) <= 0)
            ga_throw_error(expr, pnode->pos,
                           "Invalid null size of variable");
          if (mii.size() == 1)
            pnode->init_scalar_tensor(scalar_type(mii[0]));
          if (mii.size() >= 1) {
            pnode->init_vector_tensor(mii.size());
            for (size_type i = 0; i < mii.size(); ++i)
              pnode->t[i] = scalar_type(mii[i]);
          }
        } else if (!(child0->name.compare("Id"))) {
          bool valid = (child1->node_type == GA_NODE_CONSTANT);
          int n = valid ? int(round(child1->t[0])) : -1;
          if (n <= 0 || n > 100 || child1->tensor_order() > 0)
            ga_throw_error(expr, pnode->pos, "The argument of Id "
                           "should be a (small) positive integer.");
          pnode->node_type = GA_NODE_CONSTANT;
          if (n == 1)
            pnode->init_scalar_tensor(scalar_type(1));
          else {
            pnode->init_matrix_tensor(n,n);
            for (int i = 0; i < n; ++i) pnode->t(i,i) = scalar_type(1);
          }
        } else ga_throw_error(expr, pnode->children[0]->pos,
                              "Unknown special function.");
        tree.clear_children(pnode);
      } else if (child0->node_type == GA_NODE_OPERATOR) {

        // Call to a nonlinear operator

        for (size_type i = 1; i < pnode->children.size(); ++i)
          ga_valid_operand(expr, pnode->children[i]);
        all_cte = true;
        ga_nonlinear_operator::arg_list args;
        for (size_type i = 1; i < pnode->children.size(); ++i) {
          all_cte = all_cte
            && (pnode->children[i]->node_type == GA_NODE_CONSTANT);
          args.push_back(&(pnode->children[i]->t));
          if (pnode->children[i]->node_type == GA_NODE_ALLINDICES)
            ga_throw_error(expr, pnode->children[i]->pos,
                           "Colon operator is not allowed in nonlinear "
                           "operator call.");
          if (pnode->children[i]->test_function_type)
            ga_throw_error(expr, pnode->pos, "Test functions cannot be passed "
                           "as argument of a nonlinear operator.");
          if (pnode->children[i]->tensor_order() > 2)
            ga_throw_error(expr, pnode->pos, "Sorry, arguments to nonlinear "
                        "operators should only be scalar, vector or matrices");
        }
        ga_predef_operator_tab &PREDEF_OPERATORS
          = dal::singleton<ga_predef_operator_tab>::instance(0);
        ga_predef_operator_tab::T::iterator it
          = PREDEF_OPERATORS.tab.find(child0->name);
        const ga_nonlinear_operator &OP = *(it->second);
        mi.resize(0);
        if (!(OP.result_size(args, mi)))
          ga_throw_error(expr, pnode->pos,
                         "Wrong number or wrong size of arguments for the "
                         "call of nonlinear operator " + child0->name);

        pnode->test_function_type = 0;

        if (child0->der1 > args.size() || child0->der2 > args.size())
           ga_throw_error(expr, child0->pos,
                         "Invalid derivative number for nonlinear operator "
                          + child0->name);

        if (child0->der1 && child0->der2 == 0) {
          for (size_type i = 0; i < args[child0->der1-1]->sizes().size(); ++i)
            mi.push_back(args[child0->der1-1]->sizes()[i]);
          pnode->t.adjust_sizes(mi);
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            OP.derivative(args, child0->der1, pnode->t);
            tree.clear_children(pnode);
          }
        } else if (child0->der1 && child0->der2) {
          for (size_type i = 0; i < args[child0->der1-1]->sizes().size(); ++i)
            mi.push_back(args[child0->der1-1]->sizes()[i]);
          for (size_type i = 0; i < args[child0->der2-1]->sizes().size(); ++i)
            mi.push_back(args[child0->der2-1]->sizes()[i]);
          pnode->t.adjust_sizes(mi);
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            OP.second_derivative(args, child0->der1, child0->der2,
                                          pnode->t);
            tree.clear_children(pnode);
          }
        } else {
          pnode->t.adjust_sizes(mi);
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            OP.value(args, pnode->t);
            tree.clear_children(pnode);
          }
        }
      } else {

        // Access to components of a tensor
        all_cte = (child0->node_type == GA_NODE_CONSTANT);
        // cout << "child0->tensor_order() = " << child0->tensor_order();
        // cout << endl << "child0->t.sizes() = " << child0->t.sizes() << endl;
        if (pnode->children.size() != child0->tensor_order() + 1)
          ga_throw_error(expr, pnode->pos, "Bad number of indices.");
        for (size_type i = 1; i < pnode->children.size(); ++i)
          if (pnode->children[i]->node_type != GA_NODE_ALLINDICES &&
              (pnode->children[i]->node_type != GA_NODE_CONSTANT ||
               pnode->children[i]->t.size() != 1))
            ga_throw_error(expr, pnode->children[i]->pos,
                            "Indices should be constant integers or colon.");

        bgeot::multi_index mi1(size0.size()), mi2, indices;
        for (size_type i = 0; i < child0->tensor_order(); ++i) {
          if (pnode->children[i+1]->node_type == GA_NODE_ALLINDICES) {
            mi2.push_back(child0->tensor_proper_size(i));
            indices.push_back(i);
            mi1[i] = 0;
          } else {
            mi1[i] = size_type(round(pnode->children[i+1]->t[0])-1);
            if (mi1[i] >= child0->tensor_proper_size(i))
              ga_throw_error(expr, pnode->children[i+1]->pos,
                             "Index out of range, " << mi1[i]+1 << " > "
                             << child0->tensor_proper_size(i) << " .");
          }
        }
        mi.resize(0);
        for (size_type i = 0; i < child0->nb_test_functions(); ++i)
          mi.push_back(child0->t.sizes()[i]);
        for (size_type i = 0; i < mi2.size(); ++i) mi.push_back(mi2[i]);
        pnode->t.adjust_sizes(mi);
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
        pnode->interpolate_name_test1 = child0->interpolate_name_test1;
        pnode->interpolate_name_test2 = child0->interpolate_name_test2;
        pnode->qdim1 = child0->qdim1;
        pnode->qdim2 = child0->qdim2;

        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          for (bgeot::multi_index mi3(mi2.size()); !mi3.finished(mi2);
               mi3.incrementation(mi2)) {
            for (size_type j = 0; j < mi2.size(); ++j) {
              mi1[indices[j]] = mi3[j];
            }
            pnode->t(mi3) = pnode->children[0]->t(mi1);
          }
          tree.clear_children(pnode);
        } else {
          if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
            gmm::clear(pnode->t.as_vector());
            pnode->node_type = GA_NODE_ZERO;
            tree.clear_children(pnode);
          }
        }
      }
      break;

    default:GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                        << " in semantic analysis. Internal error.");
    }
    // cout << " begin hash code = " << pnode->hash_value << endl;
    pnode->hash_value = ga_hash_code(pnode);
    // cout << "node_type = " << pnode->node_type << " op_type = "
    //      << pnode->op_type << " proper hash code = " << pnode->hash_value;
    for (size_type i = 0; i < pnode->children.size(); ++i) {
      pnode->hash_value += (pnode->children[i]->hash_value)
        * 1.0101 * (pnode->symmetric_op ? scalar_type(1) : scalar_type(i+1));
    }
    // cout << " final hash code = " << pnode->hash_value << endl;
  }




  //=========================================================================
  // Extract a sub tree which consists of the corresponding node and of
  // the terms multiplying this term, but not the term in addition.
  // The aim is to expand an expression is a sum of elementary factors.
  // Complains if a nonlinear term is encountered.
  //=========================================================================
  static void ga_extract_factor(ga_tree &result_tree, pga_tree_node pnode,
                                pga_tree_node &new_pnode) {

    result_tree.clear();
    result_tree.copy_node(pnode, 0,  result_tree.root);
    new_pnode = result_tree.root;

    bool minus_sign = false;

    pga_tree_node pnode_child = pnode;
    pnode = pnode->parent;

    while (pnode) {

      size_type nbch = pnode->children.size();
      pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
      pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;


      switch (pnode->node_type) {


      case GA_NODE_OP:
        switch(pnode->op_type) {
        case GA_PLUS:
          // Nothing to do
          break;
        case GA_MINUS:
          if (child1 == pnode_child) minus_sign = !(minus_sign);
          // A remaining minus sign is added at the end if necessary.
          break;
        case GA_UNARY_MINUS: case GA_QUOTE: case GA_TRACE: case GA_DEVIATOR:
        case GA_PRINT:
          // Copy of the term
          result_tree.insert_node(result_tree.root, pnode->node_type);
          result_tree.root->op_type = pnode->op_type;
          result_tree.root->pos = pnode->pos;
          break;
        case GA_DOT: case GA_MULT: case GA_COLON: case GA_TMULT:
        case GA_DOTMULT: case GA_DIV: case GA_DOTDIV:
          // Copy of the term and other child
          result_tree.insert_node(result_tree.root, pnode->node_type);
          result_tree.root->op_type = pnode->op_type;
          result_tree.root->pos = pnode->pos;
          if (child0 == pnode_child) {
            result_tree.root->children.resize(2);
            result_tree.copy_node(child1, result_tree.root,
                                  result_tree.root->children[1]);
          } else if (child1 == pnode_child) {
            result_tree.root->children.resize(2);
            result_tree.root->children[1] = result_tree.root->children[0];
            result_tree.copy_node(child0, result_tree.root,
                                  result_tree.root->children[0]);
          } else GMM_ASSERT1(false, "Corrupted tree");
          break;
        default: GMM_ASSERT1(false, "Unexpected operation. Internal error.");
        }
        break;

      case GA_NODE_PARAMS:
        GMM_ASSERT1(child0->node_type == GA_NODE_RESHAPE, "Cannot extract a "
                    "factor which is a parameter of a nonlinear "
                    "operator/function");
        GMM_ASSERT1(child1 == pnode_child, "Cannot extract a factor of a "
                    "Reshape size parameter");
        // Copy of the term and other children
        result_tree.insert_node(result_tree.root, pnode->node_type);
        result_tree.root->pos = pnode->pos;
        result_tree.root->children.resize(pnode->children.size());
        result_tree.root->children[1] = result_tree.root->children[0];
        for (size_type i = 0; i < pnode->children.size(); ++i) {
          if (i != 1) {
            result_tree.copy_node(pnode->children[i], result_tree.root,
                                  result_tree.root->children[i]);
          }
        }
        break;

      case GA_NODE_C_MATRIX:
        result_tree.insert_node(result_tree.root, pnode->node_type);
        result_tree.root->pos = pnode->pos;
        result_tree.root->children.resize(pnode->children.size());
        for (size_type i = 0; i < pnode->children.size(); ++i)
          if (pnode_child == pnode->children[i]) {
            result_tree.root->children[i] = result_tree.root->children[0];
            result_tree.root->children[0] = 0;
          }

        for (size_type i = 0; i < pnode->children.size(); ++i) {
          if (pnode_child == pnode->children[i]) {
            pnode->children[i] = new ga_tree_node(GA_NODE_ZERO, pnode->pos);
            pnode->children[i]->init_scalar_tensor(scalar_type(0));
            pnode->children[i]->parent = pnode;
          }
        }
        break;

      default: GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                           << " in extract constant term. Internal error.");
      }

      pnode_child = pnode;
      pnode = pnode->parent;
    }

    if (minus_sign) {
      result_tree.insert_node(result_tree.root, GA_NODE_OP);
      result_tree.root->op_type = GA_UNARY_MINUS;
      result_tree.root->pos = pnode->children[0]->pos;
    }
  }

  //=========================================================================
  // Extract the constant term of degree 1 expressions
  //=========================================================================

  static bool ga_node_extract_constant_term
  (ga_tree &tree, pga_tree_node pnode, const ga_workspace &workspace,
   const mesh &m) {
    bool is_constant = true;
    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    // pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bool child_0_is_constant = (nbch <= 0) ? true :
      ga_node_extract_constant_term(tree, pnode->children[0], workspace, m);
    bool child_1_is_constant = (nbch <= 1) ? true :
      ga_node_extract_constant_term(tree, pnode->children[1], workspace, m);

    switch (pnode->node_type) {
    case GA_NODE_ZERO: is_constant = false; break;

    case GA_NODE_ELEMENTARY_VAL_TEST: case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST: case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST: case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST: case GA_NODE_XFEM_MINUS_DIVERG_TEST:
    case GA_NODE_VAL_TEST: case GA_NODE_GRAD_TEST: case GA_NODE_PREDEF_FUNC:
    case GA_NODE_HESS_TEST: case GA_NODE_DIVERG_TEST: case GA_NODE_RESHAPE:
    case GA_NODE_ELT_SIZE: case GA_NODE_ELT_K: case GA_NODE_ELT_B:
    case GA_NODE_CONSTANT: case GA_NODE_X:
    case GA_NODE_NORMAL: case GA_NODE_OPERATOR:
      is_constant = true; break;

    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
    case GA_NODE_VAL: case GA_NODE_GRAD: case GA_NODE_HESS:
    case GA_NODE_DIVERG:
      is_constant = workspace.is_constant(pnode->name);
      break;

    case GA_NODE_INTERPOLATE_VAL:
    case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS:
    case GA_NODE_INTERPOLATE_DIVERG:
      {
        if (!(workspace.is_constant(pnode->name))) {
          is_constant = false; break;
        }
        std::set<var_trans_pair> vars;
        workspace.interpolate_transformation(pnode->interpolate_name)
          ->extract_variables(workspace, vars, true, m,
                              pnode->interpolate_name);
        for (std::set<var_trans_pair>::iterator it=vars.begin();
             it != vars.end(); ++it) {
          if (!(workspace.is_constant(it->first)))
            { is_constant = false; break; }
        }
      }
      break;

    case GA_NODE_INTERPOLATE_FILTER:
      if (!child_0_is_constant) { is_constant = false; break; }
      // No break intentionally
    case GA_NODE_INTERPOLATE_VAL_TEST:
    case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST:
    case GA_NODE_INTERPOLATE_X:
    case GA_NODE_INTERPOLATE_NORMAL:
    case GA_NODE_INTERPOLATE_DERIVATIVE:
      {
        std::set<var_trans_pair> vars;
        workspace.interpolate_transformation(pnode->interpolate_name)
          ->extract_variables(workspace, vars, true, m,
                              pnode->interpolate_name);
        for (std::set<var_trans_pair>::iterator it=vars.begin();
             it != vars.end(); ++it) {
          if (!(workspace.is_constant(it->first)))
            { is_constant = false; break; }
        }
      }
      break;

    case GA_NODE_OP:
      switch(pnode->op_type) {
        case GA_PLUS: case GA_MINUS:
          if (!child_0_is_constant && !child_1_is_constant)
            { is_constant = false; break; }
          break;

      case GA_UNARY_MINUS: case GA_QUOTE: case GA_TRACE: case GA_DEVIATOR:
      case GA_PRINT:
        is_constant = child_0_is_constant;
        break;

      case GA_DOT: case GA_MULT: case GA_COLON: case GA_TMULT:
      case GA_DOTMULT: case GA_DIV: case GA_DOTDIV:
        is_constant = (child_0_is_constant && child_1_is_constant);
        break;

      default: GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX:
      for (size_type i = 0; i < pnode->children.size(); ++i) {
        if (!(ga_node_extract_constant_term(tree, pnode->children[i],
                                            workspace, m)))
          { is_constant = false; break; }
      }
      break;

    case GA_NODE_PARAMS:
      if (child0->node_type == GA_NODE_RESHAPE) {
        is_constant = child_1_is_constant;
      } else if (child0->node_type == GA_NODE_PREDEF_FUNC) {
        for (size_type i = 1; i < pnode->children.size(); ++i) {
          if (!(ga_node_extract_constant_term(tree, pnode->children[i],
                                              workspace, m)))
            { is_constant = false; break; }
        }
      } else if (child0->node_type == GA_NODE_SPEC_FUNC) {
        GMM_ASSERT1(false, "internal error");
      } else if (child0->node_type == GA_NODE_OPERATOR) {
        for (size_type i = 1; i < pnode->children.size(); ++i) {
          if (!(ga_node_extract_constant_term(tree, pnode->children[i],
                                              workspace, m)))
            { is_constant = false; break; }
        }
      } else {
        is_constant = child_0_is_constant;
      }
      break;

    default: GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                         << " in extract constant term. Internal error.");
    }

    if (!is_constant) {
      pnode->node_type = GA_NODE_ZERO;
      tree.clear_children(pnode);
    }
    return is_constant;
  }

  std::string ga_workspace::extract_constant_term(const mesh &m) {
    std::string constant_term;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];

      if (td.order == 1) {
        ga_tree local_tree = *(td.ptree);
        if (local_tree.root)
          ga_node_extract_constant_term(local_tree, local_tree.root, *this, m);
        if (local_tree.root)
          ga_semantic_analysis("", local_tree, *this, m.dim(),
                               ref_elt_dim_of_mesh(m), false, false);
        if (local_tree.root && local_tree.root->node_type != GA_NODE_ZERO) {
          constant_term += "-("+ga_tree_to_string(local_tree)+")";
        }
      }
    }
    return constant_term;
  }

  //=========================================================================
  // Extract the order zero term
  //=========================================================================

  std::string ga_workspace::extract_order0_term(void) {
    std::string term;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      if (td.order == 0) {
        ga_tree &local_tree = *(td.ptree);
        if (term.size())
          term += "+("+ga_tree_to_string(local_tree)+")";
        else
          term = "("+ga_tree_to_string(local_tree)+")";
      }
    }
    return term;
  }


  //=========================================================================
  // Extract the order one term corresponding to a certain test function
  //=========================================================================

  std::string ga_workspace::extract_order1_term(const std::string &varname) {
    std::string term;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      if (td.order == 1 && td.name_test1.compare(varname) == 0) {
        ga_tree &local_tree = *(td.ptree);
        if (term.size())
          term += "+("+ga_tree_to_string(local_tree)+")";
        else
          term = "("+ga_tree_to_string(local_tree)+")";
      }
    }
    return term;
  }


  //=========================================================================
  // Extract Neumann terms
  //=========================================================================

  static std::string ga_extract_one_Neumann_term
  (const std::string &varname,
    ga_workspace &workspace, pga_tree_node pnode) {

    size_type N = workspace.qdim(varname);
    const mesh_fem *mf = workspace.associated_mf(varname);
    GMM_ASSERT1(mf, "Works only with fem variables.");
    size_type meshdim = mf->linked_mesh().dim();
    ga_tree factor;
    pga_tree_node new_pnode;
    ga_extract_factor(factor, pnode, new_pnode);
    std::string result;
    pga_tree_node nnew_pnode = new_pnode;

    int cas = new_pnode->node_type == GA_NODE_GRAD_TEST ? 0 : 1;
    // Allow to detect Trace(Grad_Test_u)
    if (cas == 0 && new_pnode->parent &&
        new_pnode->parent->node_type == GA_NODE_OP &&
        new_pnode->parent->op_type == GA_TRACE) {
      cas = 2; nnew_pnode = new_pnode->parent;
    }
    bool ok = true;
    pga_tree_node colon_pnode = 0;
    bool quote_before_colon = false;

    // A:Grad_Test_u --> A*Normal if A is a matrix
    // Grad_Test_u:A --> A*Normal if A is a matrix
    // A*Div_Test_u  --> A*Normal if A is a scalar
    // Div_Test_u*A  --> Normal*A if A is a scalar
    // A*(Grad_Test_u)' --> (A)'*Normal if A is a matrix
    // intercaled scalar multplications and divisions are taken into account
    while (nnew_pnode->parent) {
      pga_tree_node previous_node = nnew_pnode;
      nnew_pnode = nnew_pnode->parent;

      if (nnew_pnode->node_type == GA_NODE_OP &&
          nnew_pnode->op_type == GA_MULT &&
          nnew_pnode->children[0] == previous_node &&
          nnew_pnode->children[1]->tensor_proper_size() == 1) {
      } else if (nnew_pnode->node_type == GA_NODE_OP &&
                 nnew_pnode->op_type == GA_MULT &&
                 nnew_pnode->children[1] == previous_node &&
                 nnew_pnode->children[0]->tensor_proper_size() == 1) {

      } else if (nnew_pnode->node_type == GA_NODE_OP &&
                 nnew_pnode->op_type == GA_DIV &&
                 nnew_pnode->children[0] == previous_node &&
                 nnew_pnode->children[1]->tensor_proper_size() == 1) {

      } else if (nnew_pnode->node_type == GA_NODE_OP &&
                 nnew_pnode->op_type == GA_COLON &&
                 nnew_pnode->children[0] == previous_node &&
                 nnew_pnode->children[1]->tensor_order() == 2 &&
                 colon_pnode == 0 && cas == 0) {
        std::swap(nnew_pnode->children[0],  nnew_pnode->children[1]);
        colon_pnode = nnew_pnode;
      } else if (nnew_pnode->node_type == GA_NODE_OP &&
                 nnew_pnode->op_type == GA_COLON &&
                 nnew_pnode->children[1] == previous_node &&
                 nnew_pnode->children[0]->tensor_order() == 2 &&
                 colon_pnode == 0 && cas == 0) {
        colon_pnode = nnew_pnode;
      } else if (nnew_pnode->node_type == GA_NODE_OP &&
                 nnew_pnode->op_type == GA_QUOTE &&
                 colon_pnode == 0 && cas == 0 && !quote_before_colon) {
        quote_before_colon = true;
      } else ok = false;
    }

    // cout << "analyzing factor : " <<  ga_tree_to_string(factor) << endl;
    // cout << "ok = " << int(ok) << endl;
    // cout << "colon_pnode = " << colon_pnode << endl;

    if (ok && cas == 0 && !colon_pnode) ok = false;

    if (N == 1) {
      new_pnode->node_type = GA_NODE_NORMAL;
      result = "(" + ga_tree_to_string(factor) + ")";
    } else if (ok) {
      switch (cas) {
      case 0:
        new_pnode->node_type = GA_NODE_NORMAL;
        colon_pnode->op_type = GA_MULT;
        if (quote_before_colon) {
          factor.insert_node(colon_pnode->children[0], GA_NODE_OP);
          colon_pnode->children[0]->op_type = GA_QUOTE;
          nnew_pnode = new_pnode->parent;
          while(nnew_pnode->node_type != GA_NODE_OP ||
                nnew_pnode->op_type != GA_QUOTE)
            nnew_pnode = nnew_pnode->parent;
          factor.replace_node_by_child(nnew_pnode,0);
        }
        break;
      case 1:
        new_pnode->node_type = GA_NODE_NORMAL;
        break;
      case 2:
        new_pnode->parent->node_type = GA_NODE_NORMAL;
        factor.clear_children(new_pnode->parent);
        break;
      }
      result = "(" + ga_tree_to_string(factor) + ")";

    } else {
      // General case

      result = "[";
      bgeot::multi_index mi(2); mi[0] = N; mi[1] = meshdim;

      for (size_type i = 0; i < N; ++i) {
        factor.clear_children(new_pnode);
        new_pnode->node_type = GA_NODE_C_MATRIX;
        new_pnode->nbc1 = meshdim;
        new_pnode->nbc2 = new_pnode->nbc3 = 1;
        new_pnode->t.adjust_sizes(mi);
        new_pnode->children.resize(N*meshdim);
        for (size_type j = 0; j < N; ++j) {
          for (size_type k = 0; k < meshdim; ++k) {
            if (j == i) {
              pga_tree_node param_node = new_pnode->children[k*N+j]
                = new ga_tree_node(GA_NODE_PARAMS, pnode->pos);
              new_pnode->children[k*N+j]->parent = new_pnode;
              param_node->children.resize(2);
              param_node->children[0] = new ga_tree_node(GA_NODE_NORMAL,
                                                         pnode->pos);
              param_node->children[0]->parent = param_node;
              param_node->children[1] = new ga_tree_node(GA_NODE_CONSTANT,
                                                         pnode->pos);
              param_node->children[1]->parent = param_node;
              param_node->children[1]->init_scalar_tensor(scalar_type(k));

            } else {
              new_pnode->children[k*N+j] = new ga_tree_node(GA_NODE_ZERO,
                                                            pnode->pos);
              new_pnode->children[k*N+j]->init_scalar_tensor(scalar_type(0));
              new_pnode->children[k*N+j]->parent = new_pnode;
            }
          }
        }
        result += "(" + ga_tree_to_string(factor) + ")";
        if (i < N-1) result += ";";
      }
      result += "]";
      GMM_TRACE2("Warning, generic Neumann term used: " << result);
    }

    return result;
  }


  static void ga_extract_Neumann_term_rec
  (ga_tree &tree, const std::string &varname,
   ga_workspace &workspace, pga_tree_node pnode, std::string &result) {

    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_extract_Neumann_term_rec(tree, varname, workspace,
                                  pnode->children[i], result);

    switch (pnode->node_type) {
    case GA_NODE_DIVERG_TEST: case GA_NODE_GRAD_TEST:
    case GA_NODE_ELEMENTARY_GRAD_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
      if (pnode->name.compare(varname) == 0) {
        if (result.size()) result += " + ";
        result += ga_extract_one_Neumann_term(varname, workspace, pnode);
      }
      break;
    case GA_NODE_INTERPOLATE_GRAD_TEST: case GA_NODE_INTERPOLATE_DIVERG_TEST:
          if (pnode->name.compare(varname) == 0)
            GMM_ASSERT1(false, "Do not know how to extract a "
                        "Neumann term with an interpolate transformation");
       break;
    default: break;
    }
  }

  std::string ga_workspace::extract_Neumann_term(const std::string &varname) {
    std::string result;
    for (size_type i = 0; i < trees.size(); ++i) {
      ga_workspace::tree_description &td =  trees[i];
      if (td.order == 1 && td.name_test1.compare(varname) == 0) {
        ga_tree &local_tree = *(td.ptree);
        if (local_tree.root)
          ga_extract_Neumann_term_rec(local_tree, varname, *this,
                                      local_tree.root, result);
      }
    }
    return result;
  }

  //=========================================================================
  // Derivation algorithm: derivation of a tree with respect to a variable
  //   The result tree is not ready to use. It has to be passed again in
  //   ga_semantic_analysis for enrichment.
  //=========================================================================

  static bool ga_node_mark_tree_for_variable
  (pga_tree_node pnode, const ga_workspace &workspace, const mesh &m,
   const std::string &varname,
   const std::string &interpolatename) {
    bool marked = false;
    for (size_type i = 0; i < pnode->children.size(); ++i)
      if (ga_node_mark_tree_for_variable(pnode->children[i], workspace, m,
                                         varname, interpolatename))
        marked = true;

    bool plain_node(pnode->node_type == GA_NODE_VAL ||
                    pnode->node_type == GA_NODE_GRAD ||
                    pnode->node_type == GA_NODE_HESS ||
                    pnode->node_type == GA_NODE_DIVERG);
    bool interpolate_node(pnode->node_type == GA_NODE_INTERPOLATE_VAL ||
                          pnode->node_type == GA_NODE_INTERPOLATE_GRAD ||
                          pnode->node_type == GA_NODE_INTERPOLATE_HESS ||
                          pnode->node_type == GA_NODE_INTERPOLATE_DIVERG);
    bool elementary_node(pnode->node_type == GA_NODE_ELEMENTARY_VAL ||
                         pnode->node_type == GA_NODE_ELEMENTARY_GRAD ||
                         pnode->node_type == GA_NODE_ELEMENTARY_HESS ||
                         pnode->node_type == GA_NODE_ELEMENTARY_DIVERG);
    bool xfem_node(pnode->node_type == GA_NODE_XFEM_PLUS_VAL ||
                   pnode->node_type == GA_NODE_XFEM_PLUS_GRAD ||
                   pnode->node_type == GA_NODE_XFEM_PLUS_HESS ||
                   pnode->node_type == GA_NODE_XFEM_PLUS_DIVERG ||
                   pnode->node_type == GA_NODE_XFEM_MINUS_VAL ||
                   pnode->node_type == GA_NODE_XFEM_MINUS_GRAD ||
                   pnode->node_type == GA_NODE_XFEM_MINUS_HESS ||
                   pnode->node_type == GA_NODE_XFEM_MINUS_DIVERG);
    bool interpolate_test_node
      (pnode->node_type == GA_NODE_INTERPOLATE_VAL_TEST ||
       pnode->node_type == GA_NODE_INTERPOLATE_GRAD_TEST ||
       pnode->node_type == GA_NODE_INTERPOLATE_HESS_TEST ||
       pnode->node_type == GA_NODE_INTERPOLATE_DIVERG_TEST);

    if ((plain_node || interpolate_node || elementary_node || xfem_node) &&
        (pnode->name.compare(varname) == 0 &&
         pnode->interpolate_name.compare(interpolatename) == 0)) marked = true;

    if (interpolate_node || interpolate_test_node ||
        pnode->node_type == GA_NODE_INTERPOLATE_NORMAL ||
        pnode->node_type == GA_NODE_INTERPOLATE_X) {
      std::set<var_trans_pair> vars;
      workspace.interpolate_transformation(pnode->interpolate_name)
        ->extract_variables(workspace, vars, true,
                            m, pnode->interpolate_name);
      for (std::set<var_trans_pair>::iterator it=vars.begin();
           it != vars.end(); ++it) {
        if (it->first.compare(varname) == 0 &&
            it->second.compare(interpolatename) == 0) marked = true;
      }
    }
    pnode->marked = marked;
    return marked;
  }

  static void ga_node_derivation(ga_tree &tree, const ga_workspace &workspace,
                                 const mesh &m,
                                 pga_tree_node pnode,
                                 const std::string &varname,
                                 const std::string &interpolatename,
                                 size_type order) {

    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bool mark0 = ((nbch > 0) ? child0->marked : false);
    bool mark1 = ((nbch > 1) ? child1->marked : false);
    bgeot::multi_index mi;

    switch (pnode->node_type) {
    case GA_NODE_VAL: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
      mi.resize(1); mi[0] = 2;
      for (size_type i = 0; i < pnode->tensor_order(); ++i)
        mi.push_back(pnode->tensor_proper_size(i));
      pnode->t.adjust_sizes(mi);
      if (pnode->node_type == GA_NODE_VAL)
        pnode->node_type = GA_NODE_VAL_TEST;
      else if (pnode->node_type == GA_NODE_GRAD)
        pnode->node_type = GA_NODE_GRAD_TEST;
      else if (pnode->node_type == GA_NODE_HESS)
        pnode->node_type = GA_NODE_HESS_TEST;
      else if (pnode->node_type == GA_NODE_DIVERG)
        pnode->node_type = GA_NODE_DIVERG_TEST;
      pnode->test_function_type = order;
      break;

    case GA_NODE_INTERPOLATE_VAL:
    case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS:
    case GA_NODE_INTERPOLATE_DIVERG:
      {
        bool is_val(pnode->node_type == GA_NODE_INTERPOLATE_VAL);
        bool is_grad(pnode->node_type == GA_NODE_INTERPOLATE_GRAD);
        bool is_hess(pnode->node_type == GA_NODE_INTERPOLATE_HESS);
        bool is_diverg(pnode->node_type == GA_NODE_INTERPOLATE_DIVERG);

        bool ivar = (pnode->name.compare(varname) == 0 &&
                     pnode->interpolate_name.compare(interpolatename) == 0);
        bool itrans = !ivar;
        if (!itrans) {
          std::set<var_trans_pair> vars;
          workspace.interpolate_transformation(pnode->interpolate_name)
            ->extract_variables(workspace, vars, true, m,
                                pnode->interpolate_name);
          for (std::set<var_trans_pair>::iterator it=vars.begin();
               it != vars.end(); ++it) {
            if (it->first.compare(varname) == 0 &&
                it->second.compare(interpolatename) == 0)
              itrans = true;
          }
        }

        pga_tree_node pnode_trans = pnode;
        if (is_hess) {
          GMM_ASSERT1(!itrans, "Sorry, cannot derive a hessian once more");
        } else if (itrans && ivar) {
          tree.duplicate_with_addition(pnode);
          pnode_trans = pnode->parent->children[1];
        }

        if (ivar) {
          mi.resize(1); mi[0] = 2;
          for (size_type i = 0; i < pnode->tensor_order(); ++i)
            mi.push_back(pnode->tensor_proper_size(i));
          pnode->t.adjust_sizes(mi);
          if (is_val) // --> t(Qmult*ndof,Qmult*target_dim)
            pnode->node_type = GA_NODE_INTERPOLATE_VAL_TEST;
          else if (is_grad) // --> t(Qmult*ndof,Qmult*target_dim,N)
            pnode->node_type = GA_NODE_INTERPOLATE_GRAD_TEST;
          else if (is_hess) // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pnode->node_type = GA_NODE_INTERPOLATE_HESS_TEST;
          else if (is_diverg) // --> t(Qmult*ndof)
            pnode->node_type = GA_NODE_INTERPOLATE_DIVERG_TEST;
          pnode->test_function_type = order;
        }

        if (itrans) {
          const mesh_fem *mf = workspace.associated_mf(pnode_trans->name);
          size_type q = workspace.qdim(pnode_trans->name);
          size_type n = mf->linked_mesh().dim();
          bgeot::multi_index mii = workspace.qdims(pnode_trans->name);

          if (is_val)  // --> t(target_dim*Qmult,N)
            pnode_trans->node_type = GA_NODE_INTERPOLATE_GRAD;
          else if (is_grad || is_diverg)  // --> t(target_dim*Qmult,N,N)
            pnode_trans->node_type = GA_NODE_INTERPOLATE_HESS;

          if (n > 1) {
            if (q == 1 && mii.size() <= 1) { mii.resize(1); mii[0] = n; }
            else mii.push_back(n);

            if (is_grad || is_diverg) mii.push_back(n);
          }
          pnode_trans->t.adjust_sizes(mii);
          tree.duplicate_with_operation(pnode_trans,
                                        (n > 1) ? GA_DOT : GA_MULT);
          pga_tree_node pnode_der = pnode_trans->parent->children[1];
          pnode_der->node_type = GA_NODE_INTERPOLATE_DERIVATIVE;
          if (n == 1)
            pnode_der->init_vector_tensor(2);
          else
            pnode_der->init_matrix_tensor(2, n);
          pnode_der->test_function_type = order;
          pnode_der->name = varname;
          pnode_der->interpolate_name_der = pnode_der->interpolate_name;
          pnode_der->interpolate_name = interpolatename;

          if (is_diverg) { // --> t(Qmult*ndof)
            tree.insert_node(pnode_trans->parent, GA_NODE_OP);
            pga_tree_node pnode_tr = pnode_trans->parent->parent;
            pnode_tr->op_type = GA_TRACE;
            pnode_tr->init_vector_tensor(2);
//            pnode_tr->test_function_type = order;
//            pnode_tr->name_test1 = pnode_trans->name_test1;
//            pnode_tr->name_test2 = pnode_trans->name_test2;
          }
        }
      }
      break;

    case GA_NODE_INTERPOLATE_VAL_TEST:
    case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_DIVERG_TEST:
      {
        bool is_val(pnode->node_type == GA_NODE_INTERPOLATE_VAL_TEST);
        bool is_grad(pnode->node_type == GA_NODE_INTERPOLATE_GRAD_TEST);
        bool is_diverg(pnode->node_type == GA_NODE_INTERPOLATE_DIVERG_TEST);

        pga_tree_node pnode_trans = pnode;
        const mesh_fem *mf = workspace.associated_mf(pnode_trans->name);
        size_type q = workspace.qdim(pnode_trans->name);
        size_type n = mf->linked_mesh().dim();
        bgeot::multi_index mii = workspace.qdims(pnode_trans->name);
        if (is_val) // --> t(Qmult*ndof,Qmult*target_dim,N)
          pnode_trans->node_type = GA_NODE_INTERPOLATE_GRAD_TEST;
        else if (is_grad || is_diverg) // --> t(Qmult*ndof,Qmult*target_dim,N,N)
          pnode_trans->node_type = GA_NODE_INTERPOLATE_HESS_TEST;

        if (q == 1 && mii.size() <= 1) { mii.resize(1); mii[0] = 2; }
        else mii.insert(mii.begin(), 2);

        if (n > 1) {
          mii.push_back(n);
          if (is_grad || is_diverg) mii.push_back(n);
        }
        pnode_trans->t.adjust_sizes(mii);
        // pnode_trans->test_function_type = order;
        tree.duplicate_with_operation(pnode_trans,
                                      (n > 1 ? GA_DOT : GA_MULT));
        pga_tree_node pnode_der = pnode_trans->parent->children[1];
        pnode_der->node_type = GA_NODE_INTERPOLATE_DERIVATIVE;
        if (n == 1)
          pnode_der->init_vector_tensor(2);
        else
          pnode_der->init_matrix_tensor(2, n);
        pnode_der->test_function_type = order;
        pnode_der->name = varname;
        pnode_der->interpolate_name_der = pnode_der->interpolate_name;
        pnode_der->interpolate_name = interpolatename;

        if (is_diverg) { // --> t(Qmult*ndof)
          tree.insert_node(pnode_trans->parent, GA_NODE_OP);
          pga_tree_node pnode_tr = pnode_trans->parent->parent;
          pnode_tr->op_type = GA_TRACE;
          pnode_tr->init_vector_tensor(2);
//          pnode_tr->test_function_type = order;
//          pnode_tr->name_test1 = pnode_trans->name_test1;
//          pnode_tr->name_test2 = pnode_trans->name_test2;
        }
      }
      break;

    case GA_NODE_INTERPOLATE_HESS_TEST:
      GMM_ASSERT1(false, "Sorry, cannot derive a hessian once more");
      break;

    case GA_NODE_INTERPOLATE_X:
      {
        size_type n = m.dim();
        pga_tree_node pnode_der = pnode;
        pnode_der->node_type = GA_NODE_INTERPOLATE_DERIVATIVE;
        if (n == 1)
          pnode_der->init_vector_tensor(2);
        else
          pnode_der->init_matrix_tensor(2, n);
        pnode_der->test_function_type = order;
        pnode_der->name = varname;
        pnode_der->interpolate_name_der  = pnode_der->interpolate_name;
        pnode_der->interpolate_name = interpolatename;
      }
      break;

    case GA_NODE_INTERPOLATE_NORMAL:
      GMM_ASSERT1(false, "Sorry, cannot derive the interpolated Normal");
      break;

    case GA_NODE_INTERPOLATE_DERIVATIVE:
      GMM_ASSERT1(false, "Sorry, second order transformation derivative "
                  "not taken into account");
      break;

    case GA_NODE_INTERPOLATE_FILTER:
      ga_node_derivation(tree, workspace, m, child0, varname,
                         interpolatename, order);
      break;

    case GA_NODE_ELEMENTARY_VAL:
    case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS:
    case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL:
    case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS:
    case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL:
    case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS:
    case GA_NODE_XFEM_MINUS_DIVERG:
      mi.resize(1); mi[0] = 2;
      for (size_type i = 0; i < pnode->tensor_order(); ++i)
        mi.push_back(pnode->tensor_proper_size(i));
      pnode->t.adjust_sizes(mi);
      switch(pnode->node_type) {
      case GA_NODE_ELEMENTARY_VAL:
        pnode->node_type = GA_NODE_ELEMENTARY_VAL_TEST; break;
      case GA_NODE_ELEMENTARY_GRAD:
        pnode->node_type = GA_NODE_ELEMENTARY_GRAD_TEST; break;
      case GA_NODE_ELEMENTARY_HESS:
        pnode->node_type = GA_NODE_ELEMENTARY_HESS_TEST; break;
      case GA_NODE_ELEMENTARY_DIVERG:
        pnode->node_type = GA_NODE_ELEMENTARY_DIVERG_TEST; break;
      case GA_NODE_XFEM_PLUS_VAL:
        pnode->node_type = GA_NODE_XFEM_PLUS_VAL_TEST; break;
      case GA_NODE_XFEM_PLUS_GRAD:
        pnode->node_type = GA_NODE_XFEM_PLUS_GRAD_TEST; break;
      case GA_NODE_XFEM_PLUS_HESS:
        pnode->node_type = GA_NODE_XFEM_PLUS_HESS_TEST; break;
      case GA_NODE_XFEM_PLUS_DIVERG:
        pnode->node_type = GA_NODE_XFEM_PLUS_DIVERG_TEST; break;
      case GA_NODE_XFEM_MINUS_VAL:
        pnode->node_type = GA_NODE_XFEM_MINUS_VAL_TEST; break;
      case GA_NODE_XFEM_MINUS_GRAD:
        pnode->node_type = GA_NODE_XFEM_MINUS_GRAD_TEST; break;
      case GA_NODE_XFEM_MINUS_HESS:
        pnode->node_type = GA_NODE_XFEM_MINUS_HESS_TEST; break;
      case GA_NODE_XFEM_MINUS_DIVERG:
        pnode->node_type = GA_NODE_XFEM_MINUS_DIVERG_TEST; break;
      default : GMM_ASSERT1(false, "internal error");
      }
      pnode->test_function_type = order;
      break;

    case GA_NODE_OP:
      switch(pnode->op_type) {
        case GA_PLUS: case GA_MINUS:
          if (mark0 && mark1) {
            ga_node_derivation(tree, workspace, m, child0, varname,
                               interpolatename, order);
            ga_node_derivation(tree, workspace, m, child1, varname,
                               interpolatename, order);
          } else if (mark0) {
            ga_node_derivation(tree, workspace, m, child0, varname,
                               interpolatename, order);
            tree.replace_node_by_child(pnode, 0);
          } else {
            ga_node_derivation(tree, workspace, m, child1, varname,
                               interpolatename, order);
            if (pnode->op_type == GA_MINUS) {
              pnode->op_type = GA_UNARY_MINUS;
              tree.clear_node(child0);
            }
            else
              tree.replace_node_by_child(pnode, 1);
          }
          break;

      case GA_UNARY_MINUS: case GA_QUOTE: case GA_TRACE: case GA_DEVIATOR:
      case GA_PRINT:
        ga_node_derivation(tree, workspace, m, child0, varname,
                           interpolatename, order);
        break;

      case GA_DOT: case GA_MULT: case GA_COLON: case GA_TMULT:
      case GA_DOTMULT:
        if (mark0 && mark1) {
          if (sub_tree_are_equal(child0, child1, workspace, 0) &&
              (pnode->op_type != GA_MULT || child0->tensor_order() < 2)) {
            ga_node_derivation(tree, workspace, m, child1, varname,
                               interpolatename, order);
            tree.insert_node(pnode, GA_NODE_OP);
            pnode->parent->op_type = GA_MULT;
            tree.add_child(pnode->parent);
            pga_tree_node pnode_cte = pnode->parent->children[1];
            pnode_cte->node_type = GA_NODE_CONSTANT;
            pnode_cte->init_scalar_tensor(scalar_type(2));
          } else {
            tree.duplicate_with_addition(pnode);
            if ((pnode->op_type == GA_COLON && child0->tensor_order() == 2) ||
                (pnode->op_type == GA_DOT && child0->tensor_order() == 1) ||
                pnode->op_type == GA_DOTMULT ||
                (child0->tensor_proper_size()== 1 &&
                 child1->tensor_proper_size()== 1))
              std::swap(pnode->children[0], pnode->children[1]);
            ga_node_derivation(tree, workspace, m, child0, varname,
                               interpolatename, order);
            ga_node_derivation(tree, workspace, m,
                               pnode->parent->children[1]->children[1],
                               varname, interpolatename, order);
          }
        } else if (mark0) {
          ga_node_derivation(tree, workspace, m, child0, varname,
                             interpolatename, order);
        } else
          ga_node_derivation(tree, workspace, m, child1, varname,
                             interpolatename, order);
        break;

      case GA_DIV: case GA_DOTDIV:
        if (mark1) {
          if (pnode->children[0]->node_type == GA_NODE_CONSTANT)
            gmm::scale(pnode->children[0]->t.as_vector(), scalar_type(-1));
          else {
            if (mark0) {
              tree.duplicate_with_subtraction(pnode);
              ga_node_derivation(tree, workspace, m, child0, varname,
                                 interpolatename, order);
              pnode = pnode->parent->children[1];
            } else {
              tree.insert_node(pnode, GA_NODE_OP);
              pnode->parent->op_type = GA_UNARY_MINUS;
            }
          }
          tree.insert_node(pnode->children[1], GA_NODE_PARAMS);
          pga_tree_node pnode_param = pnode->children[1];
          tree.add_child(pnode_param);
          std::swap(pnode_param->children[0], pnode_param->children[1]);
          pnode_param->children[0]->node_type = GA_NODE_PREDEF_FUNC;
          pnode_param->children[0]->name = "sqr";
          tree.insert_node(pnode, GA_NODE_OP);
          pga_tree_node pnode_mult = pnode->parent;
          if (pnode->op_type == GA_DOTDIV)
            pnode_mult->op_type = GA_DOTMULT;
          else
            pnode_mult->op_type = GA_MULT;
          pnode_mult->children.push_back(0);
          tree.copy_node(pnode_param->children[1],
                         pnode_mult, pnode_mult->children[1]);
          ga_node_derivation(tree, workspace, m, pnode_mult->children[1],
                             varname, interpolatename, order);
        } else {
          ga_node_derivation(tree, workspace, m, child0, varname,
                             interpolatename, order);
        }
        break;

      default: GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX:
      for (size_type i = 0; i < pnode->children.size(); ++i) {
        if (pnode->children[i]->marked)
          ga_node_derivation(tree, workspace, m, pnode->children[i],
                             varname, interpolatename, order);
        else {
          pnode->children[i]->init_scalar_tensor(scalar_type(0));
          pnode->children[i]->node_type = GA_NODE_ZERO;
          tree.clear_children(pnode->children[i]);
        }
      }
      break;

    case GA_NODE_PARAMS:
      if (child0->node_type == GA_NODE_RESHAPE) {
        ga_node_derivation(tree, workspace, m, pnode->children[1],
                           varname, interpolatename, order);
      } else if (child0->node_type == GA_NODE_PREDEF_FUNC) {
        std::string name = child0->name;
        ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
        const ga_predef_function &F = it->second;

        if (F.nbargs() == 1) {
          switch (F.dtype()) {
          case 0:
            GMM_ASSERT1(false, "Cannot derive function " << child0->name
                        << ". No derivative provided or not derivable function.");
          case 1:
            child0->name = F.derivative1();
            break;
          case 2: case 3:
            {
              child0->name = "DER_PDFUNC_" + child0->name;
              if (!(ga_function_exists(child0->name))) {
                if (F.dtype() == 2)
                  ga_define_function(child0->name, 1, F.derivative1());
                else {
                  std::string expr = ga_derivative_scalar_function(F.expr(),"t");
                  ga_define_function(child0->name, 1, expr);
                }
              }
              // Inline extension if the derivative is affine (for instance
              // for sqr)
              const ga_predef_function &Fp = PREDEF_FUNCTIONS[child0->name];
              if (Fp.is_affine("t")) {
                scalar_type b = Fp(scalar_type(0));
                scalar_type a = Fp(scalar_type(1)) - b;
                pnode->node_type = GA_NODE_OP;
                pnode->op_type = GA_MULT;
                child0->init_scalar_tensor(a);
                child0->node_type = (a == scalar_type(0)) ? GA_NODE_ZERO : GA_NODE_CONSTANT;
                if (b != scalar_type(0)) {
                  tree.insert_node(pnode, GA_NODE_OP);
                  pnode->parent->op_type = (b > 0) ? GA_PLUS : GA_MINUS;
                  tree.add_child(pnode->parent);
                  pga_tree_node pnode_cte = pnode->parent->children[1];
                  pnode_cte->node_type = GA_NODE_CONSTANT;
                  pnode_cte->t = pnode->t;
                  std::fill(pnode_cte->t.begin(), pnode_cte->t.end(), gmm::abs(b));
                  pnode = pnode->parent;
                }
              }
            }
            break;
          }
          if (pnode->children.size() >= 2) {
            tree.insert_node(pnode, GA_NODE_OP);
            pga_tree_node pnode_op = pnode->parent;
            if (child1->tensor_order() == 0)
              pnode_op->op_type = GA_MULT;
            else
              pnode_op->op_type = GA_DOTMULT;
            pnode_op->children.push_back(0);
            tree.copy_node(child1, pnode_op, pnode_op->children[1]);
            ga_node_derivation(tree, workspace, m, pnode_op->children[1],
                               varname, interpolatename, order);
          }
        } else {
          pga_tree_node child2 = pnode->children[2];

          if (child1->marked && child2->marked)
            tree.duplicate_with_addition(pnode);

          if (child1->marked) {
            switch (F.dtype()) {
            case 0:
              GMM_ASSERT1(false, "Cannot derive function " << child0->name
                          << ". No derivative provided");
            case 1:
              child0->name = F.derivative1();
              break;
            case 2:
              child0->name = "DER_PDFUNC1_" + child0->name;
              if (!(ga_function_exists(child0->name)))
                ga_define_function(child0->name, 2, F.derivative1());
              break;
            case 3:
              child0->name = "DER_PDFUNC1_" + child0->name;
              if (!(ga_function_exists(child0->name))) {
                std::string expr = ga_derivative_scalar_function(F.expr(), "t");
                ga_define_function(child0->name, 2, expr);
              }
              break;
            }
            tree.insert_node(pnode, GA_NODE_OP);
            pga_tree_node pnode_op = pnode->parent;
            if (child1->tensor_order() == 0)
              pnode_op->op_type = GA_MULT;
            else
              pnode_op->op_type = GA_DOTMULT;
            pnode_op->children.push_back(0);
            tree.copy_node(child1, pnode_op, pnode_op->children[1]);
            ga_node_derivation(tree, workspace, m, pnode_op->children[1],
                               varname, interpolatename, order);
          }
          if (child2->marked) {
            if (child1->marked && child2->marked)
              pnode = pnode->parent->parent->children[1];
            child0 = pnode->children[0]; child1 = pnode->children[1];
            child2 = pnode->children[2];

            switch (F.dtype()) {
            case 0:
              GMM_ASSERT1(false, "Cannot derive function " << child0->name
                          << ". No derivative provided");
            case 1:
              child0->name = F.derivative2();
              break;
            case 2:
              child0->name = "DER_PDFUNC2_" + child0->name;
              if (!(ga_function_exists(child0->name)))
                ga_define_function(child0->name, 2, F.derivative2());
              break;
            case 3:
              child0->name = "DER_PDFUNC2_" + child0->name;
              if (!(ga_function_exists(child0->name))) {
                std::string expr = ga_derivative_scalar_function(F.expr(), "u");
                ga_define_function(child0->name, 2, expr);
              }
              break;
            }
            tree.insert_node(pnode, GA_NODE_OP);
            pga_tree_node pnode_op = pnode->parent;
            if (child2->tensor_order() == 0)
              pnode_op->op_type = GA_MULT;
            else
              pnode_op->op_type = GA_DOTMULT;
            pnode_op->children.push_back(0);
            tree.copy_node(child2, pnode_op, pnode_op->children[1]);
            ga_node_derivation(tree, workspace, m, pnode_op->children[1],
                               varname, interpolatename, order);
          }
        }
      } else if (child0->node_type == GA_NODE_SPEC_FUNC) {
        GMM_ASSERT1(false, "internal error");
      } else if (child0->node_type == GA_NODE_OPERATOR) {
        if (child0->der2)
          GMM_ASSERT1(false, "Error in derivation of the assembly string. "
                      "Cannot derive again operator " <<  child0->name);

        size_type nbargs_der = 0;
        for (size_type i = 1; i < pnode->children.size(); ++i)
          if (pnode->children[i]->marked) ++nbargs_der;
        pga_tree_node pnode2 = 0;

        size_type j = 0;
        for (size_type i = 1; i < pnode->children.size(); ++i) {
          if (pnode->children[i]->marked) {
            ++j;
            if (j != nbargs_der) {
              tree.insert_node(pnode, GA_NODE_OP);
              pga_tree_node pnode_op = pnode->parent;
              pnode_op->node_type = GA_NODE_OP;
              pnode_op->op_type = GA_PLUS;
              pnode_op->children.push_back(0);
              tree.copy_node(pnode, pnode_op , pnode_op->children[1]);
              pnode2 = pnode_op->children[1];
            }
            else pnode2 = pnode;

            if (child0->der1)
              pnode2->children[0]->der2 = i;
            else
              pnode2->children[0]->der1 = i;
            tree.insert_node(pnode2, GA_NODE_OP);
            pga_tree_node pnode_op = pnode2->parent;
            // calcul de l'ordre de reduction
            size_type red = pnode->children[i]->tensor_order();
            switch (red) {
            case 0 : pnode_op->op_type = GA_MULT; break;
            case 1 : pnode_op->op_type = GA_DOT; break;
            case 2 : pnode_op->op_type = GA_COLON; break;
            default: GMM_ASSERT1(false, "Error in derivation of the assembly "
                                 "string. Bad reduction order.")
            }
            pnode_op->children.push_back(0);
            tree.copy_node(pnode->children[i], pnode_op,
                           pnode_op->children[1]);
            ga_node_derivation(tree, workspace, m, pnode_op->children[1],
                               varname, interpolatename, order);

            if (pnode2->children[0]->name.compare("Norm_sqr") == 0
                && pnode2->children[0]->der1 == 1) {
              pnode2->node_type = GA_NODE_OP;
              pnode2->op_type = GA_MULT;
              pnode2->children[0]->node_type = GA_NODE_CONSTANT;
              pnode2->children[0]->init_scalar_tensor(scalar_type(2));
            }


          }
        }

      } else {
        ga_node_derivation(tree, workspace, m, child0, varname,
                           interpolatename, order);
      }
      break;

    default: GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                         << " in derivation. Internal error.");
    }
  }

  // The tree is modified. Should be copied first and passed to
  // ga_semantic_analysis after for enrichment
  static void ga_derivative(ga_tree &tree, const ga_workspace &workspace,
                            const mesh &m, const std::string &varname,
                            const std::string &interpolatename,
                            size_type order) {
    // cout << "compute derivative of " << ga_tree_to_string(tree)
    //  << " with respect to " << varname << " : " << interpolatename << endl;
    if (!(tree.root)) return;
    if (ga_node_mark_tree_for_variable(tree.root, workspace, m, varname,
                                       interpolatename))
      ga_node_derivation(tree, workspace, m, tree.root, varname,
                         interpolatename, order);
    else
      tree.clear();
    // cout << "derived tree : " << ga_tree_to_string(tree) << endl;
  }

  static void ga_replace_test_by_cte(pga_tree_node pnode,  bool full_replace) {
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_replace_test_by_cte(pnode->children[i], full_replace);
    GMM_ASSERT1(pnode->node_type != GA_NODE_GRAD_TEST, "Invalid tree");
    GMM_ASSERT1(pnode->node_type != GA_NODE_HESS_TEST, "Invalid tree");
    GMM_ASSERT1(pnode->node_type != GA_NODE_DIVERG_TEST, "Invalid tree");
    if (pnode->node_type == GA_NODE_VAL_TEST) {
      pnode->node_type = GA_NODE_CONSTANT;
      if (full_replace) pnode->init_scalar_tensor(scalar_type(1));
    }
  }

  static std::string ga_derivative_scalar_function(const std::string expr,
                                                   const std::string &var) {
    base_vector t(1), u(1);
    ga_workspace workspace;
    workspace.add_fixed_size_variable("t", gmm::sub_interval(0,1), t);
    workspace.add_fixed_size_variable("u", gmm::sub_interval(0,1), u);
    workspace.add_function_expression(expr);
    GMM_ASSERT1(workspace.nb_trees() <= 1, "Internal error");
    if (workspace.nb_trees()) {
      ga_tree tree = *(workspace.tree_info(0).ptree);
      ga_derivative(tree, workspace, *((const mesh *)(0)), var, "", 1);
      if (tree.root) {
        ga_replace_test_by_cte(tree.root, true);
        ga_semantic_analysis(expr, tree, workspace, 1, 1, false, true);
      }
      return ga_tree_to_string(tree);
    } else return "0";
  }

  static bool ga_node_is_affine(pga_tree_node pnode) {

    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bool mark0 = ((nbch > 0) ? child0->marked : false);
    bool mark1 = ((nbch > 1) ? child1->marked : false);

    switch (pnode->node_type) {
    case GA_NODE_VAL: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
    case GA_NODE_INTERPOLATE_VAL: case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS: case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_INTERPOLATE_DERIVATIVE:
    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
      return true;
    case GA_NODE_OP:
      switch(pnode->op_type) {
      case GA_PLUS: case GA_MINUS:
        if (mark0 && mark1)
          return ga_node_is_affine(child0) &&
            ga_node_is_affine(child1);
        if (mark0) return ga_node_is_affine(child0);
        return ga_node_is_affine(child1);

      case GA_UNARY_MINUS: case GA_QUOTE: case GA_TRACE: case GA_DEVIATOR:
      case GA_PRINT: case GA_NODE_INTERPOLATE_FILTER:
        return ga_node_is_affine(child0);

      case GA_DOT: case GA_MULT: case GA_COLON: case GA_TMULT:
      case GA_DOTMULT:
        if (mark0 && mark1) return false;
        if (mark0) return ga_node_is_affine(child0);
        return ga_node_is_affine(child1);

      case GA_DIV: case GA_DOTDIV:
        if (mark1) return false;
        if (mark0) return ga_node_is_affine(child0);

      default: GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX:
        for (size_type i = 0; i < pnode->children.size(); ++i)
          if (pnode->children[i]->marked &&
              !(ga_node_is_affine(pnode->children[i])))
            return false;
        return true;

    case GA_NODE_PARAMS:
      if (child0->node_type == GA_NODE_RESHAPE)
        return ga_node_is_affine(child1);
      if (child0->node_type == GA_NODE_PREDEF_FUNC)
        return false;
      if (child0->node_type == GA_NODE_OPERATOR)
        return false;
      return ga_node_is_affine(child0);

    default: GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                         << " in derivation. Internal error.");
    }
  }

  static bool ga_is_affine(ga_tree &tree, const ga_workspace &workspace,
                           const std::string &varname,
                           const std::string &interpolatename) {
    const mesh &m = *((const mesh *)(0));
    if (tree.root && ga_node_mark_tree_for_variable(tree.root, workspace, m,
                                                    varname, interpolatename))
      return ga_node_is_affine(tree.root);
    return true;
  }


  //=========================================================================
  // Compilation of assembly trees into a list of basic instructions
  //=========================================================================

  static void add_interval_to_gis(ga_workspace &workspace, std::string varname,
                                  ga_instruction_set &gis) {
    if (workspace.variable_group_exists(varname)) {
      const std::vector<std::string> &vg = workspace.variable_group(varname);
      for (size_type i = 0; i < vg.size(); ++i)
        add_interval_to_gis(workspace, vg[i], gis);
    } else if (gis.var_intervals.find(varname) == gis.var_intervals.end()) {
      const mesh_fem *mf = workspace.associated_mf(varname);
      size_type nd = mf ? mf->nb_basic_dof() :
        gmm::vect_size(workspace.value(varname));
      gis.var_intervals[varname]=gmm::sub_interval(gis.nb_dof, nd);
      gis.nb_dof += nd;
    }
  }

  static void extend_variable_in_gis(ga_workspace &workspace,
                                     std::string varname,
                                     ga_instruction_set &gis) {
    if (workspace.variable_group_exists(varname)) {
      const std::vector<std::string> &vg = workspace.variable_group(varname);
      for (size_type i = 0; i < vg.size(); ++i)
        extend_variable_in_gis(workspace, vg[i], gis);
    } else if (gis.extended_vars.find(varname)==gis.extended_vars.end()) {
      const mesh_fem *mf = workspace.associated_mf(varname);
      if (mf->is_reduced()) {
        base_vector U(mf->nb_basic_dof());
        mf->extend_vector(workspace.value(varname), U);
        gis.really_extended_vars[varname] = U;
        gis.extended_vars[varname] = &(gis.really_extended_vars[varname]);
      } else {
        gis.extended_vars[varname] = &(workspace.value(varname));
      }
    }
  }

  static void ga_clear_node_list
  (pga_tree_node pnode, std::map<scalar_type,
   std::list<pga_tree_node> > &node_list) {
    std::list<pga_tree_node> &loc_node_list = node_list[pnode->hash_value];
    for (std::list<pga_tree_node>::iterator it = loc_node_list.begin();
         it != loc_node_list.end(); ) {
      if (*it == pnode) it = loc_node_list.erase(it); else ++it;
    }
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_clear_node_list(pnode->children[i], node_list);
  }

  static void ga_compile_node(pga_tree_node pnode,
                              ga_workspace &workspace,
                              ga_instruction_set &gis,
                              ga_instruction_set::region_mim_instructions &rmi,
                              const mesh &m, bool function_case,
                              ga_if_hierarchy &if_hierarchy) {

    if (pnode->node_type == GA_NODE_PREDEF_FUNC ||
        pnode->node_type == GA_NODE_OPERATOR ||
        pnode->node_type == GA_NODE_SPEC_FUNC ||
        pnode->node_type == GA_NODE_CONSTANT ||
        pnode->node_type == GA_NODE_ALLINDICES ||
        // pnode->node_type == GA_NODE_ZERO ||   // zero nodes can still have test functions
        pnode->node_type == GA_NODE_RESHAPE) return;

    // cout << "compiling "; ga_print_node(pnode, cout); cout << endl;

    pga_instruction pgai = 0;
    ga_if_hierarchy *pif_hierarchy = &if_hierarchy;
    ga_if_hierarchy new_if_hierarchy;

    const mesh_fem *mf1 = 0, *mf2 = 0;
    const mesh_fem **mfg1 = 0, **mfg2 = 0;
    fem_interpolation_context *pctx1 = 0, *pctx2 = 0;

    if (pnode->test_function_type) {
      if (pnode->name_test1.size())
        mf1 = workspace.associated_mf(pnode->name_test1);
      if (mf1) {
        pctx1 = &(gis.ctx);
        const std::string &intn1 = pnode->interpolate_name_test1;
        if (intn1.size()) pctx1 = &(rmi.interpolate_infos[intn1].ctx);
        if (intn1.size()
            && workspace.variable_group_exists(pnode->name_test1)) {
          ga_instruction_set::variable_group_info &vgi =
            rmi.interpolate_infos[intn1].groups_info[pnode->name_test1];
          mfg1 = &(vgi.mf); mf1 = 0;
        }
      }
      if (pnode->name_test2.size())
        mf2 = workspace.associated_mf(pnode->name_test2);
      if (mf2) {
        pctx2 = &(gis.ctx);
        const std::string &intn2 = pnode->interpolate_name_test2;
        if (intn2.size()) pctx2 = &(rmi.interpolate_infos[intn2].ctx);
        if (intn2.size()
            && workspace.variable_group_exists(pnode->name_test2)) {
          ga_instruction_set::variable_group_info &vgi =
            rmi.interpolate_infos[intn2].groups_info[pnode->name_test2];
          mfg2 = &(vgi.mf); mf2 = 0;
        }
      }
    }

    if (pnode->test_function_type == 1) {
      if (mf1 || mfg1)
        pgai = new ga_instruction_first_ind_tensor
          (pnode->t, *pctx1, pnode->qdim1, mf1, mfg1);
    } else if (pnode->test_function_type == 2) {
      if (mf2 || mfg2)
        pgai = new ga_instruction_first_ind_tensor
          (pnode->t, *pctx2, pnode->qdim2, mf2, mfg2);
    } else if (pnode->test_function_type == 3) {
      if ((mf1 || mfg1) && (mf2 || mfg2))
        pgai = new ga_instruction_two_first_ind_tensor
          (pnode->t, *pctx1, *pctx2, pnode->qdim1, mf1, mfg1,
           pnode->qdim2, mf2, mfg2);
      else if (mf1 || mfg1)
        pgai = new ga_instruction_first_ind_tensor
          (pnode->t, *pctx1, pnode->qdim1, mf1, mfg1);
      else if (mf2 || mfg2)
        pgai = new ga_instruction_second_ind_tensor
          (pnode->t, *pctx2, pnode->qdim2, mf2, mfg2);
    }
    if (pgai) rmi.instructions.push_back(pgai);

    // Optimization: detect if an equivalent node has already been compiled
    if (rmi.node_list.find(pnode->hash_value) != rmi.node_list.end()) {
      std::list<pga_tree_node> &node_list = rmi.node_list[pnode->hash_value];
      for (std::list<pga_tree_node>::iterator it = node_list.begin();
           it != node_list.end(); ++it) {
        // cout << "found potential equivalent nodes ";
        // ga_print_node(pnode, cout);
        // cout << " and "; ga_print_node(*it, cout); cout << endl;
        if (sub_tree_are_equal(pnode, *it, workspace, 1)) {
          // cout << "confirmed no transpose" << endl;
          if (pnode->t.size() == 1) {
            pgai = new ga_instruction_copy_scalar(pnode->t[0], (*it)->t[0]);
          } else {
            pgai = new ga_instruction_copy_tensor(pnode->t, (*it)->t);
          }
          rmi.instructions.push_back(pgai);
          return;
        }
        if (sub_tree_are_equal(pnode, *it, workspace, 2)) {
          // cout << "confirmed with transpose" << endl;
          if (pnode->nb_test_functions() == 2) {
            pgai = new ga_instruction_transpose_test(pnode->t, (*it)->t);
          } else {
            pgai = new ga_instruction_copy_tensor(pnode->t, (*it)->t);
          }
          rmi.instructions.push_back(pgai);
          return;
        }
        cerr << "Detected wrong equivalent nodes: ";
        ga_print_node(pnode, cerr);
        cerr << " and "; ga_print_node(*it, cout);
        cerr << " (no problem, but hash code would be adapted) " << endl;
      }
    }

    size_type interpolate_filter_inst = rmi.instructions.size();
    if (pnode->node_type == GA_NODE_INTERPOLATE_FILTER) {
      pgai = 0;
      rmi.instructions.push_back(pgai);
      if_hierarchy.increment();
      new_if_hierarchy.child_of(if_hierarchy);
      pif_hierarchy = &new_if_hierarchy;
    }

    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_compile_node(pnode->children[i], workspace, gis, rmi, m,
                      function_case, *pif_hierarchy);

    if (pnode->node_type == GA_NODE_INTERPOLATE_FILTER) {
      const std::string &intn = pnode->interpolate_name;
      ga_instruction_set::interpolate_info &inin = rmi.interpolate_infos[intn];
      pgai = new ga_instruction_interpolate_filter
        (pnode->t, inin, pnode->nbc1,
         int(rmi.instructions.size() - interpolate_filter_inst));
      rmi.instructions[interpolate_filter_inst] = pgai;
      pgai = new ga_instruction_copy_tensor(pnode->t, pnode->children[0]->t);
      rmi.instructions.push_back(pgai);
      ga_clear_node_list(pnode->children[0], rmi.node_list);
    }

    static scalar_type minus = -scalar_type(1);
    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bgeot::multi_index mi;
    const bgeot::multi_index &size0 = child0 ? child0->t.sizes() : mi;
    // const bgeot::multi_index &size1 = child1 ? child1->t.sizes() : mi;
    size_type dim0 = child0 ? child0->tensor_order() : 0;
    size_type dim1 = child1 ? child1->tensor_order() : 0;

    switch (pnode->node_type) {

    case GA_NODE_PREDEF_FUNC: case GA_NODE_OPERATOR: case GA_NODE_SPEC_FUNC:
    case GA_NODE_CONSTANT: case GA_NODE_ALLINDICES: case GA_NODE_ZERO:
    case GA_NODE_RESHAPE: case GA_NODE_INTERPOLATE_FILTER:
      break;

    case GA_NODE_X:
      GMM_ASSERT1(!function_case,
                  "No use of X is allowed in scalar functions");
      if (pnode->nbc1) {
        GA_DEBUG_ASSERT(pnode->t.size() == 1, "dimensions mismatch");
        GMM_ASSERT1(pnode->nbc1 <= m.dim(),
                    "Bad index for X in expression");
        pgai = new ga_instruction_X_component
            (pnode->t[0], gis.ctx, pnode->nbc1-1);
      } else {
        if (pnode->t.size() != m.dim())
          pnode->init_vector_tensor(m.dim());
        pgai = new ga_instruction_X(pnode->t, gis.ctx);
      }
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_ELT_SIZE:
      GMM_ASSERT1(!function_case,
                  "No use of element_size is allowed in functions");
      if (pnode->t.size() != 1) pnode->init_scalar_tensor(0);
      pgai = new ga_instruction_element_size(pnode->t, gis.elt_size);
      gis.need_elt_size = true;
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_ELT_K:
      GMM_ASSERT1(!function_case,
                  "No use of element_K is allowed in functions");
      pgai = new ga_instruction_element_K(pnode->t, gis.ctx);
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_ELT_B:
      GMM_ASSERT1(!function_case,
                  "No use of element_B is allowed in functions");
      pgai = new ga_instruction_element_B(pnode->t, gis.ctx);
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_NORMAL:
      GMM_ASSERT1(!function_case,
                  "No use of Normal is allowed in functions");
      if (pnode->t.size() != m.dim())
        pnode->init_vector_tensor(m.dim());
      pgai = new ga_instruction_Normal(pnode->t, gis.Normal);
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_INTERPOLATE_NORMAL:
      GMM_ASSERT1(!function_case,
                  "No use of Interpolate is allowed in functions");
      if (pnode->t.size() != m.dim())
        pnode->init_vector_tensor(m.dim());
      pgai = new ga_instruction_Normal(pnode->t,
                    rmi.interpolate_infos[pnode->interpolate_name].Normal);
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_INTERPOLATE_X:
      GMM_ASSERT1(!function_case,
                  "No use of Interpolate is allowed in functions");
      if (pnode->t.size() != m.dim())
        pnode->init_vector_tensor(m.dim());
      pgai = new ga_instruction_Normal(pnode->t,
                    rmi.interpolate_infos[pnode->interpolate_name].pt_y);
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_VAL: case GA_NODE_GRAD:
    case GA_NODE_HESS: case GA_NODE_DIVERG:
    case GA_NODE_ELEMENTARY_VAL: case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_ELEMENTARY_HESS: case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_VAL: case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_PLUS_HESS: case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_VAL: case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_XFEM_MINUS_HESS: case GA_NODE_XFEM_MINUS_DIVERG:
      if (function_case) {
        GMM_ASSERT1(pnode->node_type != GA_NODE_ELEMENTARY_VAL &&
                    pnode->node_type != GA_NODE_ELEMENTARY_GRAD &&
                    pnode->node_type != GA_NODE_ELEMENTARY_HESS &&
                    pnode->node_type != GA_NODE_ELEMENTARY_DIVERG,
                    "No elementary transformation is allowed in functions");
        GMM_ASSERT1(pnode->node_type != GA_NODE_XFEM_PLUS_VAL &&
                    pnode->node_type != GA_NODE_XFEM_PLUS_GRAD &&
                    pnode->node_type != GA_NODE_XFEM_PLUS_HESS &&
                    pnode->node_type != GA_NODE_XFEM_PLUS_DIVERG,
                    "Xfem_plus not allowed in functions");
        GMM_ASSERT1(pnode->node_type != GA_NODE_XFEM_MINUS_VAL &&
                    pnode->node_type != GA_NODE_XFEM_MINUS_GRAD &&
                    pnode->node_type != GA_NODE_XFEM_MINUS_HESS &&
                    pnode->node_type != GA_NODE_XFEM_MINUS_DIVERG,
                    "Xfem_plus not allowed in functions");
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        const im_data *imd = workspace.associated_im_data(pnode->name);
        GMM_ASSERT1(!mf,"No fem expression is allowed in function expression");
        GMM_ASSERT1(!imd, "No integration method data is allowed in "
                    "function expression");
        if (gmm::vect_size(workspace.value(pnode->name)) == 1)
          pgai = new ga_instruction_copy_scalar
            (pnode->t[0], (workspace.value(pnode->name))[0]);
        else
          pgai = new ga_instruction_copy_vect
            (pnode->t.as_vector(), workspace.value(pnode->name));
        rmi.instructions.push_back(pgai);
      } else {
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        const im_data *imd = workspace.associated_im_data(pnode->name);

        if (imd) {
          pgai = new ga_instruction_extract_local_im_data
            (pnode->t, *imd, workspace.value(pnode->name), gis.pai, gis.ctx,
             workspace.qdim(pnode->name));
          rmi.instructions.push_back(pgai);
        } else {

          GMM_ASSERT1(mf, "Internal error");

          GMM_ASSERT1(&(mf->linked_mesh()) == &(m),
                      "The finite element of variable " << pnode->name <<
                      " has to be defined on the same mesh than the "
                      "integration method used");

          // An instruction for extracting local dofs of the variable.
          if (rmi.local_dofs.find(pnode->name) == rmi.local_dofs.end() ||
              !(if_hierarchy.is_compatible
                (rmi.local_dofs_hierarchy[pnode->name]))) {
            rmi.local_dofs[pnode->name] = base_vector(1);
            rmi.local_dofs_hierarchy[pnode->name].push_back(if_hierarchy);
            extend_variable_in_gis(workspace, pnode->name, gis);
            // cout << "local dof of " << pnode->name << endl;
            pgai = new ga_instruction_slice_local_dofs
              (*mf, *(gis.extended_vars[pnode->name]), gis.ctx,
               rmi.local_dofs[pnode->name]);
            rmi.instructions.push_back(pgai);
          }

          // An instruction for pfp update
          if (rmi.pfps.find(mf) == rmi.pfps.end() ||
              !(if_hierarchy.is_compatible(rmi.pfps_hierarchy[mf]))) {
            rmi.pfps[mf] = 0;
            rmi.pfps_hierarchy[mf].push_back(if_hierarchy);
            pgai = new ga_instruction_update_pfp
              (*mf,  rmi.pfps[mf], gis.ctx, gis.fp_pool);
            rmi.instructions.push_back(pgai);
          }

          // At this level Xfem needs some particular storage ...

          // An instruction for the base value
          pgai = 0;
          switch (pnode->node_type) {
          case GA_NODE_VAL: case GA_NODE_ELEMENTARY_VAL:
            if (rmi.base.find(mf) == rmi.base.end() ||
               !(if_hierarchy.is_compatible(rmi.base_hierarchy[mf]))) {
              rmi.base_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_val_base
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_VAL:
            if (rmi.xfem_plus_base.find(mf) == rmi.xfem_plus_base.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_base_hierarchy[mf]))) {
              rmi.xfem_plus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_plus_val_base
                (rmi.xfem_plus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_VAL:
            if (rmi.xfem_minus_base.find(mf) == rmi.xfem_minus_base.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_base_hierarchy[mf]))) {
              rmi.xfem_minus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_minus_val_base
                (rmi.xfem_minus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_GRAD: case GA_NODE_DIVERG:
          case GA_NODE_ELEMENTARY_GRAD: case GA_NODE_ELEMENTARY_DIVERG:
            if (rmi.grad.find(mf) == rmi.grad.end() ||
                !(if_hierarchy.is_compatible(rmi.grad_hierarchy[mf]))) {
              rmi.grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_grad_base
                (rmi.grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_GRAD: case GA_NODE_XFEM_PLUS_DIVERG:
            if (rmi.xfem_plus_grad.find(mf) == rmi.xfem_plus_grad.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_grad_hierarchy[mf]))) {
              rmi.xfem_plus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_plus_grad_base
                (rmi.xfem_plus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_GRAD: case GA_NODE_XFEM_MINUS_DIVERG:
            if (rmi.xfem_minus_grad.find(mf) == rmi.xfem_minus_grad.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_grad_hierarchy[mf]))) {
              rmi.xfem_minus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_minus_grad_base
                (rmi.xfem_minus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_HESS: case GA_NODE_ELEMENTARY_HESS:
            if (rmi.hess.find(mf) == rmi.hess.end() ||
                !(if_hierarchy.is_compatible(rmi.hess_hierarchy[mf]))) {
              rmi.hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_hess_base
                (rmi.hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_HESS:
            if (rmi.xfem_plus_hess.find(mf) == rmi.xfem_plus_hess.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_hess_hierarchy[mf]))) {
              rmi.xfem_plus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_plus_hess_base
                (rmi.xfem_plus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_HESS:
            if (rmi.xfem_minus_hess.find(mf) == rmi.xfem_minus_hess.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_hess_hierarchy[mf]))) {
              rmi.xfem_minus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_minus_hess_base
                (rmi.xfem_minus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;

          default : GMM_ASSERT1(false, "Internal error");
          }
          if (pgai) rmi.instructions.push_back(pgai);

          // The eval instruction
          switch (pnode->node_type) {
          case GA_NODE_VAL: // --> t(target_dim*Qmult)
            pgai = new ga_instruction_val
              (pnode->t, rmi.base[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_GRAD: // --> t(target_dim*Qmult,N)
            pgai = new ga_instruction_grad
              (pnode->t, rmi.grad[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_HESS: // --> t(target_dim*Qmult,N,N)
            pgai = new ga_instruction_hess
              (pnode->t, rmi.hess[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_DIVERG: // --> t(1)
            pgai = new ga_instruction_diverg
              (pnode->t, rmi.grad[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_PLUS_VAL: // --> t(target_dim*Qmult)
            pgai = new ga_instruction_val
              (pnode->t, rmi.xfem_plus_base[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_PLUS_GRAD: // --> t(target_dim*Qmult,N)
            pgai = new ga_instruction_grad
              (pnode->t, rmi.xfem_plus_grad[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_PLUS_HESS: // --> t(target_dim*Qmult,N,N)
            pgai = new ga_instruction_hess
              (pnode->t, rmi.xfem_plus_hess[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_PLUS_DIVERG: // --> t(1)
            pgai = new ga_instruction_diverg
              (pnode->t, rmi.xfem_plus_grad[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_MINUS_VAL: // --> t(target_dim*Qmult)
            pgai = new ga_instruction_val
              (pnode->t, rmi.xfem_minus_base[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_MINUS_GRAD: // --> t(target_dim*Qmult,N)
            pgai = new ga_instruction_grad
              (pnode->t, rmi.xfem_minus_grad[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_MINUS_HESS: // --> t(target_dim*Qmult,N,N)
            pgai = new ga_instruction_hess
              (pnode->t, rmi.xfem_minus_hess[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_XFEM_MINUS_DIVERG: // --> t(1)
            pgai = new ga_instruction_diverg
              (pnode->t, rmi.xfem_minus_grad[mf],
               rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
            break;
          case GA_NODE_ELEMENTARY_VAL:
            { // --> t(target_dim*Qmult)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_val
                (pnode->t, rmi.base[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_GRAD:
            { // --> t(target_dim*Qmult,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_grad
                (pnode->t, rmi.grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_HESS:
            { // --> t(target_dim*Qmult,N,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_hess
                (pnode->t, rmi.hess[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_DIVERG:
            { // --> t(1)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_diverg
                (pnode->t, rmi.grad[mf],
                 rmi.local_dofs[pnode->name], workspace.qdim(pnode->name),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          default: break;
          }
          rmi.instructions.push_back(pgai);
        }
      }
      break;

    case GA_NODE_INTERPOLATE_VAL: case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_INTERPOLATE_HESS: case GA_NODE_INTERPOLATE_DIVERG:
      {
        extend_variable_in_gis(workspace, pnode->name, gis);

        const mesh_fem *mfn = workspace.associated_mf(pnode->name), **mfg = 0;
        const std::string &intn = pnode->interpolate_name;
        const base_vector *Un = gis.extended_vars[pnode->name], **Ug = 0;
        fem_interpolation_context *pctx = &(rmi.interpolate_infos[intn].ctx);
        const mesh **m2 = &(rmi.interpolate_infos[intn].m);
        if (workspace.variable_group_exists(pnode->name)) {
          ga_instruction_set::variable_group_info &vgi =
            rmi.interpolate_infos[intn].groups_info[pnode->name];
          mfg = &(vgi.mf); mfn = 0; Ug = &(vgi.U); Un = 0;
        }

        if (pnode->node_type == GA_NODE_INTERPOLATE_VAL) {
          pgai = new ga_instruction_interpolate_val // --> t(target_dim*Qmult)
            (pnode->t, m2, mfn, mfg, Un, Ug, *pctx, workspace.qdim(pnode->name));
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_GRAD) {
          pgai = new ga_instruction_interpolate_grad // --> t(target_dim*Qmult,N)
            (pnode->t, m2, mfn, mfg, Un, Ug, *pctx, workspace.qdim(pnode->name));
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_HESS) {
          pgai = new ga_instruction_interpolate_hess // --> t(target_dim*Qmult,N,N)
            (pnode->t, m2, mfn, mfg, Un, Ug, *pctx, workspace.qdim(pnode->name));
        } else { // --> t(1)
          pgai = new ga_instruction_interpolate_diverg
            (pnode->t, m2, mfn, mfg, Un, Ug, *pctx, workspace.qdim(pnode->name));
        }
        rmi.instructions.push_back(pgai);
      }
      break;

    case GA_NODE_INTERPOLATE_DERIVATIVE:
      GMM_ASSERT1(!function_case,
                  "No use of Interpolate is allowed in functions");
      pgai = new ga_instruction_copy_tensor_possibly_void
        (pnode->t,
         rmi.interpolate_infos[pnode->interpolate_name_der]
         .derivatives[var_trans_pair(pnode->name, pnode->interpolate_name)]);
      rmi.instructions.push_back(pgai);
      break;

    case GA_NODE_VAL_TEST: case GA_NODE_GRAD_TEST:
    case GA_NODE_HESS_TEST: case GA_NODE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_VAL_TEST: case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST: case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST: case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST: case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      // GMM_ASSERT1(!function_case,
      //            "Test functions not allowed in functions");
      {
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        if (mf) {
          GMM_ASSERT1(&(mf->linked_mesh()) == &(m),
                      "The finite element of variable " << pnode->name <<
                      " and the applied integration method have to be"
                      " defined on the same mesh");

          // An instruction for pfp update
          if (rmi.pfps.find(mf) == rmi.pfps.end() ||
              !(if_hierarchy.is_compatible(rmi.pfps_hierarchy[mf]))) {
            rmi.pfps[mf] = 0;
            rmi.pfps_hierarchy[mf].push_back(if_hierarchy);
            pgai = new ga_instruction_update_pfp
              (*mf,  rmi.pfps[mf], gis.ctx, gis.fp_pool);
            rmi.instructions.push_back(pgai);
          }

          // An instruction for the base value
          pgai = 0;
          switch (pnode->node_type) {
          case GA_NODE_VAL_TEST: case GA_NODE_ELEMENTARY_VAL_TEST:
             if (rmi.base.find(mf) == rmi.base.end() ||
                !(if_hierarchy.is_compatible(rmi.base_hierarchy[mf]))) {
              rmi.base_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_val_base
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
             }
             break;
          case GA_NODE_XFEM_PLUS_VAL_TEST:
            if (rmi.xfem_plus_base.find(mf) == rmi.xfem_plus_base.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_base_hierarchy[mf]))) {
              rmi.xfem_plus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_plus_val_base
                (rmi.xfem_plus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_VAL_TEST:
            if (rmi.xfem_minus_base.find(mf) == rmi.xfem_minus_base.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_base_hierarchy[mf]))) {
              rmi.xfem_minus_base_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_minus_val_base
                (rmi.xfem_minus_base[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_GRAD_TEST: case GA_NODE_DIVERG_TEST:
          case GA_NODE_ELEMENTARY_GRAD_TEST:
          case GA_NODE_ELEMENTARY_DIVERG_TEST:
            if (rmi.grad.find(mf) == rmi.grad.end() ||
                !(if_hierarchy.is_compatible(rmi.grad_hierarchy[mf]))) {
              rmi.grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_grad_base
                (rmi.grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_GRAD_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
            if (rmi.xfem_plus_grad.find(mf) == rmi.xfem_plus_grad.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_grad_hierarchy[mf]))) {
              rmi.xfem_plus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_plus_grad_base
                (rmi.xfem_plus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_GRAD_TEST:
          case GA_NODE_XFEM_MINUS_DIVERG_TEST:
            if (rmi.xfem_minus_grad.find(mf) == rmi.xfem_minus_grad.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_grad_hierarchy[mf]))) {
              rmi.xfem_minus_grad_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_minus_grad_base
                (rmi.xfem_minus_grad[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_HESS_TEST: case GA_NODE_ELEMENTARY_HESS_TEST:
            if (rmi.hess.find(mf) == rmi.hess.end() ||
                !(if_hierarchy.is_compatible(rmi.hess_hierarchy[mf]))) {
              rmi.hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_hess_base
                (rmi.hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_PLUS_HESS_TEST:
            if (rmi.xfem_plus_hess.find(mf) == rmi.xfem_plus_hess.end() ||
             !(if_hierarchy.is_compatible(rmi.xfem_plus_hess_hierarchy[mf]))) {
              rmi.xfem_plus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_plus_hess_base
                (rmi.xfem_plus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;
          case GA_NODE_XFEM_MINUS_HESS_TEST:
            if (rmi.xfem_minus_hess.find(mf) == rmi.xfem_minus_hess.end() ||
            !(if_hierarchy.is_compatible(rmi.xfem_minus_hess_hierarchy[mf]))) {
              rmi.xfem_minus_hess_hierarchy[mf].push_back(if_hierarchy);
              pgai = new ga_instruction_xfem_minus_hess_base
                (rmi.xfem_minus_hess[mf], gis.ctx, *mf, rmi.pfps[mf]);
            }
            break;

          default : GMM_ASSERT1(false, "Internal error");
          }
          if (pgai) rmi.instructions.push_back(pgai);

          // The copy of the real_base_value
          switch(pnode->node_type) {
          case GA_NODE_VAL_TEST: // --> t(Qmult*ndof,Qmult*target_dim)
            pgai = new ga_instruction_copy_val_base
              (pnode->t, rmi.base[mf], mf->get_qdim());
            break;
          case GA_NODE_GRAD_TEST: // --> t(Qmult*ndof,Qmult*target_dim,N)
            pgai = new ga_instruction_copy_grad_base
              (pnode->t, rmi.grad[mf], mf->get_qdim());
            break;
          case GA_NODE_HESS_TEST: // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pgai = new ga_instruction_copy_hess_base
              (pnode->t, rmi.hess[mf], mf->get_qdim());
            break;
          case GA_NODE_DIVERG_TEST: // --> t(Qmult*ndof)
            pgai = new ga_instruction_copy_diverg_base
              (pnode->t, rmi.grad[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_VAL_TEST: // -->t(Qmult*ndof,Qmult*target_dim)
            pgai = new ga_instruction_copy_val_base
              (pnode->t, rmi.xfem_plus_base[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_GRAD_TEST: // --> t(Qmult*ndof,Qmult*target_dim,N)
            pgai = new ga_instruction_copy_grad_base
              (pnode->t, rmi.xfem_plus_grad[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_HESS_TEST: // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pgai = new ga_instruction_copy_hess_base
              (pnode->t, rmi.xfem_plus_hess[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_PLUS_DIVERG_TEST: // --> t(Qmult*ndof)
            pgai = new ga_instruction_copy_diverg_base
              (pnode->t, rmi.xfem_plus_grad[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_VAL_TEST: // -->t(Qmult*ndof,Qmult*target_dim)
            pgai = new ga_instruction_copy_val_base
              (pnode->t, rmi.xfem_minus_base[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_GRAD_TEST: // --> t(Qmult*ndof,Qmult*target_dim,N)
            pgai = new ga_instruction_copy_grad_base
              (pnode->t, rmi.xfem_minus_grad[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_HESS_TEST: // --> t(Qmult*ndof,Qmult*target_dim,N,N)
            pgai = new ga_instruction_copy_hess_base
              (pnode->t, rmi.xfem_minus_hess[mf], mf->get_qdim());
            break;
          case GA_NODE_XFEM_MINUS_DIVERG_TEST: // --> t(Qmult*ndof)
            pgai = new ga_instruction_copy_diverg_base
              (pnode->t, rmi.xfem_minus_grad[mf], mf->get_qdim());
            break;
          case GA_NODE_ELEMENTARY_VAL_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_val_base
                (pnode->t, rmi.base[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_GRAD_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_grad_base
                (pnode->t, rmi.grad[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_HESS_TEST:
            { // --> t(Qmult*ndof,Qmult*target_dim,N,N)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_hess_base
                (pnode->t, rmi.hess[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          case GA_NODE_ELEMENTARY_DIVERG_TEST:
            { // --> t(Qmult*ndof)
              ga_instruction_set::elementary_trans_info &eti
                = rmi.elementary_trans_infos[pnode->elementary_name];
              pgai = new ga_instruction_elementary_transformation_diverg_base
                (pnode->t, rmi.grad[mf], mf->get_qdim(),
                 workspace.elementary_transformation(pnode->elementary_name),
                 *mf, gis.ctx, eti.M, &(eti.mf), eti.icv);
            }
            break;
          default: break;
          }
          rmi.instructions.push_back(pgai);
        }
        add_interval_to_gis(workspace, pnode->name, gis);
        gis.max_dof = std::max
          (gis.max_dof, workspace.interval_of_variable(pnode->name).last());
      }
      break;

    case GA_NODE_INTERPOLATE_VAL_TEST: case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST: case GA_NODE_INTERPOLATE_DIVERG_TEST:
      {
        const mesh_fem *mfn = workspace.associated_mf(pnode->name), **mfg = 0;
        const std::string &intn = pnode->interpolate_name;
        fem_interpolation_context *pctx = &(rmi.interpolate_infos[intn].ctx);
        const mesh **m2 = &(rmi.interpolate_infos[intn].m);
        if (workspace.variable_group_exists(pnode->name)) {
          ga_instruction_set::variable_group_info &vgi =
            rmi.interpolate_infos[intn].groups_info[pnode->name];
          mfg = &(vgi.mf); mfn = 0;
        }

        if (pnode->node_type == GA_NODE_INTERPOLATE_VAL_TEST) {
          // --> t(Qmult*ndof,Qmult*target_dim)
          pgai = new ga_instruction_interpolate_val_base
            (pnode->t, m2, mfn, mfg, *pctx, workspace.qdim(pnode->name));
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_GRAD_TEST) {
           // --> t(Qmult*ndof,Qmult*target_dim,N)
          pgai = new ga_instruction_interpolate_grad_base
            (pnode->t, m2, mfn, mfg, *pctx, workspace.qdim(pnode->name));
        } else if (pnode->node_type == GA_NODE_INTERPOLATE_HESS_TEST) {
           // --> t(Qmult*ndof,Qmult*target_dim,N,N)
          pgai = new ga_instruction_interpolate_hess_base
            (pnode->t, m2, mfn, mfg, *pctx, workspace.qdim(pnode->name));
        } else { // if (pnode->node_type == GA_NODE_INTERPOLATE_DIVERG_TEST) {
           // --> t(Qmult*ndof)
          pgai = new ga_instruction_interpolate_diverg_base
            (pnode->t, m2, mfn, mfg, *pctx, workspace.qdim(pnode->name));
        }
        rmi.instructions.push_back(pgai);
      }
      break;

     case GA_NODE_OP:
       switch(pnode->op_type) {

       case GA_PLUS:
         if (pnode->t.size() == 1) {
           GA_DEBUG_ASSERT(child0->t.size() == 1,
                           "Internal error: child0 not scalar");
           GA_DEBUG_ASSERT(child1->t.size() == 1,
                           "Internal error: child1 not scalar");
           pgai = new ga_instruction_scalar_add
             (pnode->t[0], child0->t[0], child1->t[0]);
         } else {
           pgai = new ga_instruction_add(pnode->t, child0->t, child1->t);
         }
         rmi.instructions.push_back(pgai);
         break;

       case GA_MINUS:
         if (pnode->t.size() == 1) {
           // GA_DEBUG_ASSERT(pnode->nb_test_functions() == 0,
           //              "Internal error: non zero number of test functions");
           GA_DEBUG_ASSERT(child0->t.size() == 1,
                           "Internal error: child0 not scalar");
           GA_DEBUG_ASSERT(child1->t.size() == 1,
                           "Internal error: child1 not scalar");
           pgai = new ga_instruction_scalar_sub
             (pnode->t[0], child0->t[0], child1->t[0]);
         } else {
           pgai = new ga_instruction_sub(pnode->t, child0->t, child1->t);
         }
         rmi.instructions.push_back(pgai);
         break;

       case GA_UNARY_MINUS:
         if (pnode->t.size() == 1) {
           GA_DEBUG_ASSERT(child0->t.size() == 1, "Internal error");
           pgai = new ga_instruction_scalar_scalar_mult
             (pnode->t[0], child0->t[0], minus);
         } else {
           pgai = new ga_instruction_scalar_mult(pnode->t, child0->t, minus);
         }
         rmi.instructions.push_back(pgai);
         break;


       case GA_DOT: case GA_COLON: case GA_MULT:
         {
           size_type s1 = (child0->t.size()*child1->t.size())/pnode->t.size();
           size_type s2 = size_type(round(sqrt(scalar_type(s1))));

           pgai = 0;
           if (pnode->op_type == GA_DOT || pnode->op_type == GA_COLON ||
               (pnode->op_type == GA_MULT && dim0 == 4) ||
               (pnode->op_type == GA_MULT && dim1 <= 1) ||
               child0->t.size() == 1 || child1->t.size() == 1) {

             if (child0->t.size() == 1 && child1->t.size() == 1) {
               pgai = new ga_instruction_scalar_scalar_mult
                 (pnode->t[0], child0->t[0], child1->t[0]);
             }
             else if (child0->t.size() == 1)
               pgai = new ga_instruction_scalar_mult
                 (pnode->t, child1->t, child0->t[0]);
             else if (child1->t.size() == 1)
               pgai = new ga_instruction_scalar_mult
                 (pnode->t, child0->t, child1->t[0]);
             else if (pnode->test_function_type < 3) {
               if (child0->tensor_proper_size() == 1)
                 pgai = new ga_instruction_simple_tmult
                   (pnode->t, child0->t, child1->t);
               else {
                 if (s2 == 2) // Unroll loop test ... to be extended
                   pgai = new ga_instruction_reduction_2
                     (pnode->t, child0->t, child1->t);
                 else
                   pgai = new ga_instruction_reduction
                     (pnode->t, child0->t, child1->t, s2);
               }
             } else {
               if (child1->test_function_type == 1 ||
                   child1->test_function_type == 3) {
                 if (child1->test_function_type == 3 ||
                     child1->tensor_proper_size() <= s2) {
                   if (s2 == 2) // Unroll loop test ... to be extended
                     pgai = new ga_instruction_reduction_2
                       (pnode->t, child0->t, child1->t);
                   else
                     pgai = new ga_instruction_reduction
                       (pnode->t, child0->t, child1->t, s2);
                 } else {
                   pgai = new ga_instruction_spec_reduction
                     (pnode->t, child1->t, child0->t, s2);
                 }
               } else if (child1->test_function_type == 0 ||
                          (child0->tensor_proper_size() == s2 &&
                           child1->tensor_proper_size() == s2)) {
                 if (s2 == 2) // Unroll loop test ... to be extended
                   pgai = new ga_instruction_reduction_2
                     (pnode->t, child1->t, child0->t);
                 else
                   pgai = new ga_instruction_reduction
                     (pnode->t, child1->t, child0->t, s2);
               } else {
                 if (child0->tensor_proper_size() == s2)
                   pgai = new ga_instruction_reduction
                     (pnode->t, child1->t, child0->t, s2);
                 else if (child1->tensor_proper_size() == s2)
                   pgai = new ga_instruction_spec_reduction
                     (pnode->t, child0->t, child1->t, s2);
                 else
                   pgai = new ga_instruction_spec2_reduction
                     (pnode->t, child0->t, child1->t, s2);
               }
             }


           } else { // GA_MULT

             if (pnode->test_function_type < 3) {
               if (child1->tensor_proper_size() == 1)
                 pgai = new ga_instruction_simple_tmult
                   (pnode->t, child1->t, child0->t);
               else if (child0->tensor_proper_size() == 1)
                 pgai = new ga_instruction_simple_tmult
                   (pnode->t, child0->t, child1->t);
               else {
                 if (dim0 == 2)
                   pgai = new ga_instruction_matrix_mult
                     (pnode->t, child0->t, child1->t);
               }
             } else {
               if (child1->tensor_proper_size() == 1) {
                 if (child1->test_function_type == 0 ||
                     child1->test_function_type == 1)
                   pgai = new ga_instruction_simple_tmult
                     (pnode->t, child1->t, child0->t);
                 else
                   pgai = new ga_instruction_spec_tmult
                     (pnode->t, child0->t, child1->t,
                      child0->tensor_proper_size(),
                      child1->tensor_proper_size());
               } else if (child0->tensor_proper_size() == 1) {
                 if (child0->test_function_type == 0 ||
                     child0->test_function_type == 1)
                   pgai = new ga_instruction_simple_tmult
                     (pnode->t, child0->t, child1->t);
                 else
                   pgai = new ga_instruction_spec_tmult
                     (pnode->t, child1->t, child0->t,
                      child1->tensor_proper_size(),
                      child0->tensor_proper_size());
               } else if (dim0 == 2) {
                 if (child1->test_function_type != 2)
                   pgai = new ga_instruction_matrix_mult
                     (pnode->t, child0->t, child1->t);
                 else
                   pgai = new ga_instruction_matrix_mult_spec
                     (pnode->t, child0->t, child1->t);
               }
             }

           }
           GMM_ASSERT1(pgai, "Internal error");
           rmi.instructions.push_back(pgai);
         }
         break;

       case GA_DIV:
         if (child0->t.size() == 1 && child1->t.size() == 1) {
           pgai = new ga_instruction_scalar_scalar_div
             (pnode->t[0], child0->t[0], child1->t[0]);
         } else if (child1->t.size() == 1) {
           pgai = new ga_instruction_scalar_div
             (pnode->t, child0->t, child1->t[0]);
         } else GMM_ASSERT1(false, "Internal error");
         rmi.instructions.push_back(pgai);
         break;

       case GA_PRINT:
         pgai = new ga_instruction_copy_tensor(pnode->t, child0->t);
         rmi.instructions.push_back(pgai);
         pgai = new ga_instruction_print_tensor(pnode->t, child0, gis.ctx,
                                                gis.nbpt, gis.ipt);
         rmi.instructions.push_back(pgai);
         break;

       case GA_TRACE:
         {
           size_type N = (child0->tensor_proper_size() == 1) ? 1:size0.back();
           pgai = new ga_instruction_trace(pnode->t, child0->t, N);
           rmi.instructions.push_back(pgai);
         }
         break;

       case GA_DEVIATOR:
         {
           size_type N = (child0->tensor_proper_size() == 1) ? 1:size0.back();
           pgai = new ga_instruction_deviator(pnode->t, child0->t, N);
           rmi.instructions.push_back(pgai);
         }
         break;

       case GA_QUOTE:
         if (pnode->tensor_proper_size() != 1) {
           pgai = new ga_instruction_transpose(pnode->t, child0->t);
           rmi.instructions.push_back(pgai);
         } else {
           pgai = new ga_instruction_copy_tensor(pnode->t, child0->t);
           rmi.instructions.push_back(pgai);
         }
         break;

       case GA_DOTMULT:

         if (child0->t.size() == 1 && child1->t.size() == 1) {
           pgai = new ga_instruction_scalar_scalar_mult
             (pnode->t[0], child0->t[0], child1->t[0]);
         } else if (child0->t.size() == 1)
           pgai = new ga_instruction_scalar_mult
             (pnode->t, child1->t, child0->t[0]);
         else if (child1->t.size() == 1)
           pgai = new ga_instruction_scalar_mult
             (pnode->t, child0->t, child1->t[0]);
         else if (child1->test_function_type == 0)
           pgai = new ga_instruction_dotmult
             (pnode->t, child0->t, child1->t);
         else if (child0->test_function_type == 0)
           pgai = new ga_instruction_dotmult
             (pnode->t, child1->t, child0->t);
         else if (child0->test_function_type == 1)
           pgai = new ga_instruction_dotmult_spec
             (pnode->t, child0->t, child1->t);
         else
           pgai = new ga_instruction_dotmult_spec
             (pnode->t, child1->t, child0->t);

         rmi.instructions.push_back(pgai);
         break;


       case GA_DOTDIV:
         if (child0->t.size() == 1 && child1->t.size() == 1) {
           pgai = new ga_instruction_scalar_scalar_div
             (pnode->t[0], child0->t[0], child1->t[0]);
         } else if (child1->t.size() == 1) {
           pgai = new ga_instruction_scalar_div
             (pnode->t, child0->t, child1->t[0]);
         } else if (child1->test_function_type == 0) {
           pgai = new ga_instruction_dotdiv
             (pnode->t, child0->t, child1->t);
         } else GMM_ASSERT1(false, "Internal error");
         rmi.instructions.push_back(pgai);
         break;


       case GA_TMULT:
         if (child0->t.size() == 1 && child1->t.size() == 1) {
           pgai = new ga_instruction_scalar_scalar_mult
             (pnode->t[0], child0->t[0], child1->t[0]);
         } else if (child0->t.size() == 1)
           pgai = new ga_instruction_scalar_mult
             (pnode->t, child1->t, child0->t[0]);
         else if (child1->t.size() == 1)
           pgai = new ga_instruction_scalar_mult
             (pnode->t, child0->t, child1->t[0]);
         else if (child1->test_function_type == 0)
           pgai = new ga_instruction_simple_tmult
             (pnode->t, child0->t, child1->t);
         else if (child1->tensor_proper_size() == 1)
           pgai = new ga_instruction_spec2_tmult
             (pnode->t, child0->t, child1->t);
         else
           pgai = new ga_instruction_spec_tmult
             (pnode->t, child0->t, child1->t,
              child0->tensor_proper_size(),
              child1->tensor_proper_size());

         rmi.instructions.push_back(pgai);
         break;

       default:GMM_ASSERT1(false, "Unexpected operation. Internal error.");
       }
       break;

    case GA_NODE_C_MATRIX:
      {
        size_type nbc1 = pnode->nbc1, nbc2 = pnode->nbc2, nbc3 = pnode->nbc3;
        size_type nbl = pnode->children.size() / (nbc1*nbc2*nbc3);
        if (pnode->test_function_type) {
          std::vector<base_tensor *> components(pnode->children.size());

          if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < pnode->children.size(); ++i)
              components[i]  = &(pnode->children[i]->t);
          } else if (nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc1; ++j)
                components[i+j*nbl] = &(pnode->children[i*nbc1+j]->t);
          } else {
            size_type n = 0;
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc3; ++j)
                for (size_type k = 0; k < nbc2; ++k)
                  for (size_type l = 0; l < nbc1; ++l)
                    components[i+j*nbl+k*nbl*nbc3+l*nbc2*nbc3*nbl]
                      = &(pnode->children[n++]->t);
          }
          pgai = new ga_instruction_c_matrix_with_tests(pnode->t, components);
        } else {
          std::vector<scalar_type *> components(pnode->children.size());
          if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < pnode->children.size(); ++i)
              components[i]  = &(pnode->children[i]->t[0]);
          } else if (nbc2 == 1 && nbc3 == 1) {
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc1; ++j)
                components[i+j*nbl] = &(pnode->children[i*nbc1+j]->t[0]);
          } else {
            size_type n = 0;
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc3; ++j)
                for (size_type k = 0; k < nbc2; ++k)
                  for (size_type l = 0; l < nbc1; ++l)
                    components[i+j*nbl+k*nbl*nbc3+l*nbc2*nbc3*nbl]
                      = &(pnode->children[n++]->t[0]);
          }
          pgai = new ga_instruction_simple_c_matrix(pnode->t, components);
        }
        rmi.instructions.push_back(pgai);
      }
      break;

    case GA_NODE_PARAMS:
      if (child0->node_type == GA_NODE_RESHAPE) {
        pgai = new ga_instruction_copy_tensor(pnode->t, child1->t);
        rmi.instructions.push_back(pgai);
      } else if (child0->node_type == GA_NODE_PREDEF_FUNC) {

        std::string name = child0->name;
        ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
        const ga_predef_function &F = it->second;
        size_type nbargs = F.nbargs();
        pga_tree_node child2 = (nbargs == 2) ? pnode->children[2] : child1;

        if (nbargs == 1) {
          if (child1->t.size() == 1) {
            if (F.ftype() == 0)
              pgai = new ga_instruction_eval_func_1arg_1res
                (pnode->t[0], child1->t[0], F.f1());
            else
              pgai = new ga_instruction_eval_func_1arg_1res_expr
                (pnode->t[0], child1->t[0], F);
          } else {
            if (F.ftype() == 0)
              pgai = new ga_instruction_eval_func_1arg
                (pnode->t, child1->t, F.f1());
            else
              pgai = new ga_instruction_eval_func_1arg_expr
                (pnode->t, child1->t, F);
          }
        } else {
          if (child1->t.size() == 1 && child2->t.size() == 1) {
            if (F.ftype() == 0)
              pgai = new ga_instruction_eval_func_2arg_1res
                (pnode->t[0], child1->t[0], child2->t[0], F.f2());
            else
              pgai = new ga_instruction_eval_func_2arg_1res_expr
                (pnode->t[0], child1->t[0], child2->t[0], F);
          } else if (child1->t.size() == 1) {
            if (F.ftype() == 0)
              pgai = new ga_instruction_eval_func_2arg_first_scalar
                (pnode->t, child1->t, child2->t, F.f2());
            else
              pgai = new ga_instruction_eval_func_2arg_first_scalar_expr
                (pnode->t, child1->t, child2->t, F);
          } else if (child2->t.size() == 1) {
            if (F.ftype() == 0)
              pgai = new ga_instruction_eval_func_2arg_second_scalar
                (pnode->t, child1->t, child2->t, F.f2());
            else
              pgai = new ga_instruction_eval_func_2arg_second_scalar_expr
                (pnode->t, child1->t, child2->t, F);
          } else {
            if (F.ftype() == 0)
              pgai = new ga_instruction_eval_func_2arg
                (pnode->t, child1->t, child2->t, F.f2());
            else
              pgai = new ga_instruction_eval_func_2arg_expr
                (pnode->t, child1->t, child2->t, F);
          }
        }
        rmi.instructions.push_back(pgai);

      } else if (child0->node_type == GA_NODE_SPEC_FUNC) {

        GMM_ASSERT1(false, "Internal error");

      } else if (child0->node_type == GA_NODE_OPERATOR) {

        ga_predef_operator_tab &PREDEF_OPERATORS
          = dal::singleton<ga_predef_operator_tab>::instance(0);
        ga_predef_operator_tab::T::iterator it
          = PREDEF_OPERATORS.tab.find(child0->name);
        const ga_nonlinear_operator &OP = *(it->second);
        ga_nonlinear_operator::arg_list args;
        for (size_type i = 1; i < pnode->children.size(); ++i)
          args.push_back(&(pnode->children[i]->t));

        if (child0->der1 && child0->der2 == 0) {
           pgai = new ga_instruction_eval_derivative_OP
             (pnode->t, OP, args, child0->der1);
        } else if (child0->der1 && child0->der2) {
          pgai = new ga_instruction_eval_second_derivative_OP
             (pnode->t, OP, args, child0->der1, child0->der2);
        } else {
          pgai = new ga_instruction_eval_OP(pnode->t, OP, args);
        }
        rmi.instructions.push_back(pgai);

      } else { // Access to a component of the tensor
        bgeot::multi_index mi1(size0.size()), indices;
        if (pnode->t.size() == 1) {
          for (size_type i = 0; i < child0->tensor_order(); ++i)
            mi1[i] = size_type(round(pnode->children[i+1]->t[0])-1);
          pgai = new ga_instruction_copy_scalar(pnode->t[0], child0->t(mi1));
        } else {
          size_type nb_test = pnode->nb_test_functions();
          for (size_type i = 0; i < nb_test; ++i) indices.push_back(i);
          for (size_type i = 0; i < child0->tensor_order(); ++i) {
            if (pnode->children[i+1]->node_type != GA_NODE_ALLINDICES)
              mi1[i+nb_test] = size_type(round(pnode->children[i+1]->t[0])-1);
            else
              indices.push_back(i+nb_test);
          }
          pgai = new ga_instruction_tensor_slice(pnode->t, child0->t,
                                                 mi1, indices);
        }
        rmi.instructions.push_back(pgai);
      }

      break;

    default:GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                        << " in compilation. Internal error.");
    }
    rmi.node_list[pnode->hash_value].push_back(pnode);
  }

  static void ga_compile_function(ga_workspace &workspace,
                                  ga_instruction_set &gis, bool scalar) {
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      ga_workspace::tree_description &td = workspace.tree_info(i);

      gis.trees.push_back(*(td.ptree));
      pga_tree_node root = gis.trees.back().root;
      if (root) {
        GMM_ASSERT1(!scalar || (root->t.size() == 1),
                    "The result of the given expression is not a scalar");
        ga_instruction_set::region_mim rm(td.mim, td.rg);
        gis.whole_instructions[rm].m = td.m;
        ga_if_hierarchy if_hierarchy;
        ga_compile_node(root, workspace, gis,
                        gis.whole_instructions[rm],*(td.m),true,if_hierarchy);

        gis.coeff = scalar_type(1);
        pga_instruction pgai = 0;
        if (scalar) {
          pgai = new ga_instruction_scalar_assembly
            (root->t, workspace.assembled_potential(), gis.coeff);

        } else {
          workspace.assembled_tensor() = root->t;
          pgai = new ga_instruction_add_to
            (workspace.assembled_tensor(), root->t);
        }
        gis.whole_instructions[rm].instructions.push_back(pgai);
      }
    }
  }

  static bool ga_node_used_interpolates
  (pga_tree_node pnode, ga_workspace &workspace,
   std::map<std::string, std::set<std::string> > &interpolates,
   std::set<std::string> &interpolates_der) {
    bool found = false;
    bool intrpl(pnode->node_type == GA_NODE_INTERPOLATE_VAL ||
                pnode->node_type == GA_NODE_INTERPOLATE_GRAD ||
                pnode->node_type == GA_NODE_INTERPOLATE_HESS ||
                pnode->node_type == GA_NODE_INTERPOLATE_DIVERG);
    bool intrpl_test(pnode->node_type == GA_NODE_INTERPOLATE_VAL_TEST ||
                     pnode->node_type == GA_NODE_INTERPOLATE_GRAD_TEST ||
                     pnode->node_type == GA_NODE_INTERPOLATE_HESS_TEST ||
                     pnode->node_type == GA_NODE_INTERPOLATE_DIVERG_TEST);

    if (intrpl || intrpl_test ||
        pnode->node_type == GA_NODE_INTERPOLATE_FILTER ||
        pnode->node_type == GA_NODE_INTERPOLATE_NORMAL ||
        pnode->node_type == GA_NODE_INTERPOLATE_X) {
      interpolates[pnode->interpolate_name].size();
      if (intrpl || intrpl_test) {
        if (workspace.variable_group_exists(pnode->name))
          interpolates[pnode->interpolate_name].insert(pnode->name);
      }
      found = true;
    }
    if (pnode->node_type == GA_NODE_INTERPOLATE_DERIVATIVE) {
      interpolates_der.insert(pnode->interpolate_name_der);
      interpolates[pnode->interpolate_name_der].size();
      if (workspace.variable_group_exists(pnode->name))
          interpolates[pnode->interpolate_name_der].insert(pnode->name);
    }
    for (size_type i = 0; i < pnode->children.size(); ++i)
      found = ga_node_used_interpolates(pnode->children[i], workspace,
                                        interpolates, interpolates_der)
        || found;
    return found;
  }


  static void ga_compile_interpolate_trans
  (pga_tree_node pnode, ga_workspace &workspace, ga_instruction_set &gis,
   ga_instruction_set::region_mim_instructions &rmi, const mesh &m) {

    std::set<std::string> interpolates_der;
    std::map<std::string, std::set<std::string> > transformations;
    ga_node_used_interpolates(pnode, workspace, transformations,
                              interpolates_der);

    for (std::map<std::string, std::set<std::string> >::iterator
         it = transformations.begin(); it != transformations.end(); ++it) {
      bool compute_der = (interpolates_der.find(it->first)
                          != interpolates_der.end());
      if (rmi.transformations.find(it->first) == rmi.transformations.end() ||
          (compute_der && rmi.transformations_der.find(it->first)
           == rmi.transformations_der.end())) {
        rmi.transformations[it->first].size();
        gis.transformations.insert(it->first);
        if (compute_der) rmi.transformations_der.insert(it->first);
        pga_instruction pgai = new ga_instruction_transformation_call
          (workspace, rmi.interpolate_infos[it->first],
           workspace.interpolate_transformation(it->first), gis.ctx,
           gis.Normal, m, compute_der);
        rmi.instructions.push_back(pgai);
      }

      for (std::set<std::string>::iterator itt = it->second.begin();
           itt != it->second.end(); ++itt) {
        if (rmi.transformations[it->first].find(*itt)
            == rmi.transformations[it->first].end()) {
          pga_instruction pgai = new ga_instruction_update_group_info
            (workspace, gis, rmi.interpolate_infos[it->first],
             *itt, rmi.interpolate_infos[it->first].groups_info[*itt]);
          rmi.instructions.push_back(pgai);
          rmi.transformations[it->first].insert(*itt);
        }
      }
    }
  }

  static void ga_compile_interpolation(ga_workspace &workspace,
                                       ga_instruction_set &gis) {
    gis.transformations.clear();
    gis.whole_instructions.clear();
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      ga_workspace::tree_description &td = workspace.tree_info(i);
      if (td.order == 0) {
        gis.trees.push_back(*(td.ptree));

        // Semantic analysis mainly to evaluate fixed size variables and data
        const mesh *m = td.m;
        GMM_ASSERT1(m, "Internal error");
        ga_semantic_analysis("", gis.trees.back(), workspace, m->dim(),
                             ref_elt_dim_of_mesh(*m), true, false);
        pga_tree_node root = gis.trees.back().root;
        if (root) {
          // Compile tree
          ga_instruction_set::region_mim rm(td.mim, td.rg);
          ga_instruction_set::region_mim_instructions &rmi
            = gis.whole_instructions[rm];
          rmi.m = td.m;
          // rmi.interpolate_infos.clear();
          ga_compile_interpolate_trans(root, workspace, gis, rmi, *(td.m));
          ga_compile_node(root, workspace, gis,rmi, *(td.m), false,
                          rmi.current_hierarchy);

          // After compile tree
          workspace.assembled_tensor() = root->t;
          pga_instruction pgai = new ga_instruction_add_to
            (workspace.assembled_tensor(), root->t);
          rmi.instructions.push_back(pgai);
        }
      }
    }
  }

  static void ga_compile(ga_workspace &workspace, ga_instruction_set &gis,
                         size_type order) {
    gis.transformations.clear();
    gis.whole_instructions.clear();
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      ga_workspace::tree_description &td = workspace.tree_info(i);
      if (td.order == order) {
        gis.trees.push_back(*(td.ptree));

        // Semantic analysis mainly to evaluate fixed size variables and data
        ga_semantic_analysis("", gis.trees.back(), workspace,
                             td.mim->linked_mesh().dim(),
                             ref_elt_dim_of_mesh(td.mim->linked_mesh()),
                             true, false);
        pga_tree_node root = gis.trees.back().root;
        if (root) {
          // Compiling tree
          // cout << "Will compile "; ga_print_node(root, cout); cout << endl;

          ga_instruction_set::region_mim rm(td.mim, td.rg);
          ga_instruction_set::region_mim_instructions &rmi
            = gis.whole_instructions[rm];
          rmi.m = td.m;
          // rmi.interpolate_infos.clear();
          ga_compile_interpolate_trans(root, workspace, gis, rmi, *(td.m));
          ga_compile_node(root, workspace, gis, rmi, *(td.m), false,
                          rmi.current_hierarchy);
          // cout << "compilation finished "; ga_print_node(root, cout);
          // cout << endl;

          // Addition of an assembly instruction
          pga_instruction pgai = 0;
          switch(order) {
          case 0:
            pgai = new ga_instruction_scalar_assembly
              (root->t, workspace.assembled_potential(), gis.coeff);
            break;
          case 1:
            {
              const mesh_fem *mf = workspace.associated_mf(root->name_test1);
              const mesh_fem **mfg = 0;
              add_interval_to_gis(workspace, root->name_test1, gis);

              if (mf) {
                fem_interpolation_context *pctx = &(gis.ctx);
                const std::string &intn1 = root->interpolate_name_test1;
                if (intn1.size()) pctx = &(rmi.interpolate_infos[intn1].ctx);

                const gmm::sub_interval *Ir = 0, *In = 0;
                if (intn1.size() &&
                    workspace.variable_group_exists(root->name_test1)) {
                  ga_instruction_set::variable_group_info &vgi =
                    rmi.interpolate_infos[intn1].groups_info[root->name_test1];
                  Ir = &(vgi.Ir); In = &(vgi.In);
                  mfg = &(vgi.mf); mf = 0;
                } else {
                  Ir = &(gis.var_intervals[root->name_test1]);
                  In = &(workspace.interval_of_variable(root->name_test1));
                }
                pgai = new ga_instruction_fem_vector_assembly
                  (root->t, workspace.unreduced_vector(),
                   workspace.assembled_vector(), *pctx, *Ir, *In, mf, mfg,
                   gis.coeff);
              } else {
                pgai = new ga_instruction_vector_assembly
                    (root->t, workspace.assembled_vector(),
                     workspace.interval_of_variable(root->name_test1),
                     gis.coeff);
              }
            }
            break;
          case 2:
            {
              const mesh_fem *mf1 = workspace.associated_mf(root->name_test1);
              const mesh_fem *mf2 = workspace.associated_mf(root->name_test2);
              const mesh_fem **mfg1 = 0, **mfg2 = 0;
              const std::string &intn1 = root->interpolate_name_test1;
              const std::string &intn2 = root->interpolate_name_test2;
              fem_interpolation_context *pctx1 = &(gis.ctx);
              bool interpolate = false;
              if (intn1.size()) {
                pctx1
                  = &(rmi.interpolate_infos[root->interpolate_name_test1].ctx);
                interpolate = true;
              }
              fem_interpolation_context *pctx2 = &(gis.ctx);
              if (intn2.size()) {
                pctx2
                  = &(rmi.interpolate_infos[root->interpolate_name_test2].ctx);
                interpolate = true;
              }

              add_interval_to_gis(workspace, root->name_test1, gis);
              add_interval_to_gis(workspace, root->name_test2, gis);

              const gmm::sub_interval *Ir1 = 0, *In1 = 0, *Ir2 = 0, *In2 = 0;
              const scalar_type *alpha1 = 0, *alpha2 = 0;

              if (intn1.size() &&
                  workspace.variable_group_exists(root->name_test1)) {
                ga_instruction_set::variable_group_info &vgi =
                  rmi.interpolate_infos[intn1].groups_info[root->name_test1];
                Ir1 = &(vgi.Ir); In1 = &(vgi.In);
                mfg1 = &(vgi.mf); mf1 = 0;
                alpha1 = &(vgi.alpha);
              } else {
                alpha1 = &(workspace.factor_of_variable(root->name_test1));
                Ir1 = &(gis.var_intervals[root->name_test1]);
                In1 = &(workspace.interval_of_variable(root->name_test1));
              }

              if (intn2.size() &&
                  workspace.variable_group_exists(root->name_test2)) {
                ga_instruction_set::variable_group_info &vgi =
                  rmi.interpolate_infos[intn2].groups_info[root->name_test2];
                Ir2 = &(vgi.Ir); In2 = &(vgi.In);
                mfg2 = &(vgi.mf); mf2 = 0;
                alpha2 = &(vgi.alpha);
              } else {
                alpha2 = &(workspace.factor_of_variable(root->name_test2));
                Ir2 = &(gis.var_intervals[root->name_test2]);
                In2 = &(workspace.interval_of_variable(root->name_test2));
              }


              pgai = new ga_instruction_matrix_assembly
                <model_real_sparse_matrix>
                (root->t, workspace.unreduced_matrix(),
                 workspace.assembled_matrix(), *pctx1, *pctx2,
                 *Ir1, *In1, *Ir2, *In2, mf1, mfg1, mf2, mfg2,
                 gis.coeff, *alpha1, *alpha2, gis.nbpt, gis.ipt,
                 td.elem, interpolate);
              break;
            }
          }
          if (pgai) gis.whole_instructions[rm].instructions.push_back(pgai);
        }
      }
    }
  }


  //=========================================================================
  // Execution of a compiled set of assembly terms
  //=========================================================================


  static void ga_function_exec(ga_instruction_set &gis) {

    ga_instruction_set::instructions_set::iterator it
      = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {
      ga_instruction_list &gil = it->second.instructions;
      for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
    }
  }

  static void ga_interpolation_exec(ga_instruction_set &gis,
                                    ga_workspace &workspace,
                                    ga_interpolation_context &gic) {
    base_matrix G;
    base_small_vector un, up;

    std::set<std::string>::iterator iti = gis.transformations.begin();
    for (; iti != gis.transformations.end(); ++iti)
      (workspace.interpolate_transformation(*iti))->init(workspace);

    ga_instruction_set::instructions_set::iterator it
      = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {

      const getfem::mesh_im &mim = *(it->first.mim());
      const mesh_region &region = *(it->first.region());
      const getfem::mesh &m = *(it->second.m);
      GMM_ASSERT1(&m == &(gic.linked_mesh()),
                  "Incompatibility of meshes in interpolation");
      size_type P = m.dim();
      ga_instruction_list &gil = it->second.instructions;

      // iteration on elements (or faces of elements)
      std::vector<size_type> ind;
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);
      for (getfem::mr_visitor v(rg, m); !v.finished(); ++v) {
        if (gic.use_mim()) {
          if (!mim.convex_index().is_in(v.cv())) continue;
          gis.pai = mim.int_method_of_element(v.cv())->approx_method();
        }

        ind.resize(0);
        const bgeot::stored_point_tab &spt
          = gic.points_for_element(v.cv(), v.f(), ind);

        if (&spt && ind.size() && spt.size()) {
          bgeot::vectors_to_base_matrix(G, m.points_of_convex(v.cv()));
          size_type N = G.nrows();
          bgeot::pgeometric_trans pgt = m.trans_of_convex(v.cv());

          if (gis.ctx.have_pgp() && gis.ctx.pgt() == pgt) {
            gis.ctx = fem_interpolation_context(gis.ctx.pgp(), 0, 0, G,
                                                v.cv(), v.f());
          } else {
            if (!(gic.use_pgp(v.cv()))) {
              gis.ctx = fem_interpolation_context(pgt, 0, spt[0], G,
                                                  v.cv(), v.f());
            } else {
              bgeot::pgeotrans_precomp pgp = gis.gp_pool(pgt, &spt);
              gis.ctx = fem_interpolation_context(pgp, 0, 0, G,
                                                  v.cv(), v.f());
            }
          }

          if (gis.need_elt_size)
            gis.elt_size = m.convex_radius_estimate(v.cv()) * scalar_type(2);

          // iterations on interpolation points
          gis.nbpt = spt.size();
          for (size_type ii = 0; ii < ind.size(); ++ii) {
            gis.ipt = ind[ii];
            if (gis.ctx.have_pgp()) gis.ctx.set_ii(gis.ipt);
            else gis.ctx.set_xref(spt[gis.ipt]);

            if (ii == 0 || !(pgt->is_linear())) {
              // Computation of unit normal vector in case of a boundary
              if (v.f() != short_type(-1)) {
                const base_matrix& B = gis.ctx.B();
                up.resize(N); un.resize(P);
                gmm::copy(pgt->normals()[v.f()], un);
                gmm::mult(B, un, up);
                scalar_type nup = gmm::vect_norm2(up);
                gmm::scale(up,1.0/nup);
                gis.Normal = up;
              } else gis.Normal.resize(0);
            }
            gmm::clear(workspace.assembled_tensor().as_vector());
            for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
            gic.store_result(v.cv(), gis.ipt, workspace.assembled_tensor());
          }
        }
      }
    }
    iti = gis.transformations.begin();
    for (; iti != gis.transformations.end(); ++iti)
      (workspace.interpolate_transformation(*iti))->finalize();

    gic.finalize();
  }

  static void ga_interpolation_single_point_exec
  (ga_instruction_set &gis, ga_workspace &workspace,
   const fem_interpolation_context &ctx_x, const base_small_vector &Normal,
   const mesh &interp_mesh) {
    gis.ctx = ctx_x;
    gis.Normal = Normal;
    gmm::clear(workspace.assembled_tensor().as_vector());
    gis.nbpt = 1;
    gis.ipt = 0;

    ga_instruction_set::instructions_set::iterator
      it = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {
      const getfem::mesh &m = *(it->second.m);
      GMM_ASSERT1(&m == &interp_mesh,
                  "Incompatibility of meshes in interpolation");
      ga_instruction_list &gil = it->second.instructions;
      for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
    }
  }

  static void ga_exec(ga_instruction_set &gis, ga_workspace &workspace) {
    base_matrix G;
    base_small_vector un, up;
    scalar_type J(0);

    std::set<std::string>::iterator iti = gis.transformations.begin();
    for (; iti != gis.transformations.end(); ++iti)
      (workspace.interpolate_transformation(*iti))->init(workspace);

    ga_instruction_set::instructions_set::iterator it
      = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {

      const getfem::mesh_im &mim = *(it->first.mim());
      const getfem::mesh &m = *(it->second.m);

      GMM_ASSERT1(&m == &(mim.linked_mesh()), "Incompatibility of meshes");
      size_type P = m.dim();
      ga_instruction_list &gil = it->second.instructions;
      const mesh_region &region = *(it->first.region());

      // iteration on elements (or faces of elements)
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);
      for (getfem::mr_visitor v(rg, m); !v.finished(); ++v) {
        if (mim.convex_index().is_in(v.cv())) {
          // cout << "proceed with element " << v.cv() << endl;
          bgeot::vectors_to_base_matrix(G, m.points_of_convex(v.cv()));
          size_type N = G.nrows();
          bgeot::pgeometric_trans pgt = m.trans_of_convex(v.cv());
          pintegration_method pim = mim.int_method_of_element(v.cv());
          if (pim->type() == IM_NONE) continue;
          GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods cannot "
                      "be used in high level generic assembly");
          const bgeot::stored_point_tab &spt
            = pim->approx_method()->integration_points();
          if (spt.size()) {
            if (gis.ctx.have_pgp() && gis.pai == pim->approx_method() &&
                gis.ctx.pgt() == pgt) {
              gis.ctx = fem_interpolation_context(gis.ctx.pgp(), 0, 0, G,
                                                  v.cv(), v.f());
            } else {
              if (pim->approx_method()->is_built_on_the_fly()) {
                gis.ctx = fem_interpolation_context(pgt, 0, spt[0], G,
                                                    v.cv(), v.f());
              } else {
                bgeot::pgeotrans_precomp pgp = gis.gp_pool(pgt, &spt);
                gis.ctx = fem_interpolation_context(pgp, 0, 0, G,
                                                    v.cv(), v.f());
              }
            }
            gis.pai = pim->approx_method();
            if (gis.need_elt_size)
              gis.elt_size = m.convex_radius_estimate(v.cv()) * scalar_type(2);
            // iterations on Gauss points
            gis.nbpt = gis.pai->nb_points_on_convex();
            size_type first_ind = 0;
            if (v.f() != short_type(-1)) {
              gis.nbpt = gis.pai->nb_points_on_face(v.f());
              first_ind = gis.pai->ind_first_point_on_face(v.f());
            }
            for (gis.ipt = 0; gis.ipt < gis.nbpt; ++(gis.ipt)) {
              if (gis.ctx.have_pgp()) gis.ctx.set_ii(first_ind+gis.ipt);
              else gis.ctx.set_xref(spt[first_ind+gis.ipt]);
              if (gis.ipt == 0 || !(pgt->is_linear())) {
                J = gis.ctx.J();
                // Computation of unit normal vector in case of a boundary
                if (v.f() != short_type(-1)) {
                  up.resize(N); un.resize(P);
                  gmm::copy(pgt->normals()[v.f()], un);
                  gmm::mult(gis.ctx.B(), un, up);
                  scalar_type nup = gmm::vect_norm2(up);
                  J *= nup;
                  gmm::scale(up,1.0/nup);
                  gis.Normal = up;
                } else gis.Normal.resize(0);
              }
              gis.coeff = J * gis.pai->coeff(first_ind+gis.ipt);
              for (size_type j = 0; j < gil.size(); ++j) j += gil[j]->exec();
            }
          }
        }
      }
      GA_DEBUG_INFO("-----------------------------");
    }
    iti = gis.transformations.begin();
    for (; iti != gis.transformations.end(); ++iti)
      (workspace.interpolate_transformation(*iti))->finalize();
  }

  //=========================================================================
  // User defined functions
  //=========================================================================

  ga_function::ga_function(const ga_workspace &workspace_,
                           const std::string &e)
    : local_workspace(true, workspace_), expr(e), gis(0) {}

  ga_function::ga_function(const model &md, const std::string &e)
    : local_workspace(md), expr(e), gis(0) {}

  ga_function::ga_function(const std::string &e)
    : local_workspace(), expr(e), gis(0) {}

  ga_function::ga_function(const ga_function &gaf)
    : local_workspace(gaf.local_workspace), expr(gaf.expr), gis(0)
  { if (gaf.gis) compile(); }

  void ga_function::compile(void) const {
    if (gis) delete gis;
    gis = new ga_instruction_set;
    local_workspace.clear_expressions();
    local_workspace.add_function_expression(expr);
    ga_compile_function(local_workspace, *gis, false);
  }

  ga_function &ga_function::operator =(const ga_function &gaf) {
    if (gis) delete gis; gis = 0;
    local_workspace = gaf.local_workspace; expr = gaf.expr;
    if (gaf.gis) compile();
    return *this;
  }

  ga_function::~ga_function() { if (gis) delete gis; gis = 0; }

  const base_tensor &ga_function::eval(void) const {
    GMM_ASSERT1(gis, "Uncompiled function");
    gmm::clear(local_workspace.assembled_tensor().as_vector());
    ga_function_exec(*gis);
    return local_workspace.assembled_tensor();
  }

  void ga_function::derivative(const std::string &var) {
    GMM_ASSERT1(gis, "Uncompiled function");
    if (local_workspace.nb_trees()) {
      ga_tree tree = *(local_workspace.tree_info(0).ptree);
      ga_derivative(tree, local_workspace, *((const mesh *)(0)), var, "", 1);
      if (tree.root) {
        ga_semantic_analysis(expr, tree, local_workspace, 1, 1, false, true);
        // To be improved to suppress test functions in the expression ...
        // ga_replace_test_by_cte do not work in all operations like
        // vector components x(1)
        // ga_replace_test_by_cte(tree.root, false);
        // ga_semantic_analysis(expr, tree, local_workspace, 1, 1,
        //                      false, true);
      }
      expr = ga_tree_to_string(tree);
    }
    if (gis) delete gis; gis = 0;
    compile();
  }

  //=========================================================================
  // Interpolation functions
  //=========================================================================

  // general Interpolation
  void ga_interpolation(ga_workspace &workspace,
                        ga_interpolation_context &gic) {
    ga_instruction_set gis;
    ga_compile_interpolation(workspace, gis);
    ga_interpolation_exec(gis, workspace, gic);
  }

  // Interpolation on a Lagrange fem on the same mesh
  struct ga_interpolation_context_fem_same_mesh
    : public ga_interpolation_context {
    base_vector &result;
    std::vector<int> dof_count;
    const mesh_fem &mf;
    bool initialized;
    size_type s;

    virtual const bgeot::stored_point_tab &
    points_for_element(size_type cv, short_type f,
                       std::vector<size_type> &ind) const {
      pfem pf = mf.fem_of_element(cv);
      GMM_ASSERT1(pf->is_lagrange(),
                  "Only Lagrange fems can be used in interpolation");

      if (f != short_type(-1)) {

        for (size_type i = 0;
             i < pf->node_convex(cv).structure()->nb_points_of_face(f); ++i)
          ind.push_back
            (pf->node_convex(cv).structure()->ind_points_of_face(f)[i]);
      } else {
        for (size_type i = 0; i < pf->node_convex(cv).nb_points(); ++i)
          ind.push_back(i);
      }

      return *(pf->node_tab(cv));
    }

    virtual bool use_pgp(size_type) const { return true; }
    virtual bool use_mim(void) const { return false; }

    void init_(size_type si, size_type q, size_type qmult) {
      s = si;
      gmm::resize(result, qmult * mf.nb_basic_dof());
      gmm::clear(result);
      gmm::resize(dof_count, mf.nb_basic_dof()/q);
      gmm::clear(dof_count);
      initialized = true;
    }

    virtual void store_result(size_type cv, size_type i, base_tensor &t) {
      size_type si = t.size();
      size_type q = mf.get_qdim();
      size_type qmult = si / q;
      GMM_ASSERT1( (si % q) == 0, "Incompatibility between the mesh_fem and "
                   "the size of the expression to be interpolated");
      if (!initialized) { init_(si, q, qmult); }
      GMM_ASSERT1(s == si, "Internal error");
      size_type idof = mf.ind_basic_dof_of_element(cv)[i*q];
      gmm::add(t.as_vector(),
               gmm::sub_vector(result, gmm::sub_interval(qmult*idof, s)));
      (dof_count[idof/q])++;
    }

    virtual void finalize(void) {
      std::vector<size_type> data(3);
      data[0] = initialized ? result.size() : 0;
      data[1] = initialized ? dof_count.size() : 0;
      data[2] = initialized ? s : 0;
      MPI_MAX_VECTOR(data);
      if (!initialized) {
        if (data[0]) {
          gmm::resize(result, data[0]);
          gmm::clear(result);
          gmm::resize(dof_count, data[1]);
          gmm::clear(dof_count);
          s = data[2];
        } else {
          gmm::clear(result);
        }
      }
      if (initialized)
        GMM_ASSERT1(gmm::vect_size(result) == data[0] &&  s == data[2] &&
                  gmm::vect_size(dof_count) == data[1], "Incompatible sizes");
      MPI_SUM_VECTOR(result);
      MPI_SUM_VECTOR(dof_count);
      for (size_type i = 0; i < dof_count.size(); ++i)
        if (dof_count[i])
          gmm::scale(gmm::sub_vector(result, gmm::sub_interval(s*i, s)),
                     scalar_type(1) / scalar_type(dof_count[i]));
    }

    virtual const mesh &linked_mesh(void) { return mf.linked_mesh(); }

    ga_interpolation_context_fem_same_mesh(const mesh_fem &mf_, base_vector &r)
      : result(r), mf(mf_), initialized(false) {
      GMM_ASSERT1(!(mf.is_reduced()),
                  "Interpolation on reduced fem is not allowed");
    }
  };

  void ga_interpolation_Lagrange_fem
  (ga_workspace &workspace, const mesh_fem &mf, base_vector &result) {
    ga_interpolation_context_fem_same_mesh gic(mf, result);
    ga_interpolation(workspace, gic);
  }

  void ga_interpolation_Lagrange_fem
  (const getfem::model &md, const std::string &expr, const mesh_fem &mf,
   base_vector &result, const mesh_region &rg) {
    ga_workspace workspace(md);
    workspace.add_interpolation_expression(expr, mf.linked_mesh(), rg);
    ga_interpolation_Lagrange_fem(workspace, mf, result);
  }

  // Interpolation on a cloud of points
  struct ga_interpolation_context_mti
    : public ga_interpolation_context {
    base_vector &result;
    const mesh_trans_inv &mti;
    bool initialized;
    size_type s, nbdof;


    virtual const bgeot::stored_point_tab &
    points_for_element(size_type cv, short_type,
                       std::vector<size_type> &ind) const {
      std::vector<size_type> itab;
      mti.points_on_convex(cv, itab);
      std::vector<base_node> pt_tab(itab.size());
      for (size_type i = 0; i < itab.size(); ++i) {
        pt_tab[i] = mti.reference_coords()[itab[i]];
        ind.push_back(i);
      }
      return *(store_point_tab(pt_tab));
    }

    virtual bool use_pgp(size_type) const { return false; }
    virtual bool use_mim(void) const { return false; }

    virtual void store_result(size_type cv, size_type i, base_tensor &t) {
      size_type si = t.size();
      if (!initialized) {
        s = si;
        gmm::resize(result, s * nbdof);
        gmm::clear(result);
        initialized = true;
      }
      GMM_ASSERT1(s == si, "Internal error");
      size_type ipt = mti.point_on_convex(cv, i);
      size_type dof_t = mti.id_of_point(ipt);
      gmm::add(t.as_vector(),
               gmm::sub_vector(result, gmm::sub_interval(s*dof_t, s)));
    }

    virtual void finalize(void) {
      std::vector<size_type> data(2);
      data[0] = initialized ? result.size() : 0;
      data[1] = initialized ? s : 0;
      MPI_MAX_VECTOR(data);
      if (!initialized) {
        if (data[0]) {
          gmm::resize(result, data[0]);
          gmm::clear(result);
          s = data[1];
        } else {
           gmm::clear(result);
        }
      }
      if (initialized)
        GMM_ASSERT1(gmm::vect_size(result) == data[0] &&  s == data[1],
                    "Incompatible sizes");
      MPI_SUM_VECTOR(result);
    }

    virtual const mesh &linked_mesh(void) { return mti.linked_mesh(); }

    ga_interpolation_context_mti(const mesh_trans_inv &mti_, base_vector &r,
                                 size_type nbdof_ = size_type(-1))
      : result(r), mti(mti_), initialized(false), nbdof(nbdof_) {
      if (nbdof == size_type(-1)) nbdof = mti.nb_points();
    }
  };

  // Distribute to be parallelized
  void ga_interpolation_mti
  (const getfem::model &md, const std::string &expr, mesh_trans_inv &mti,
   base_vector &result, const mesh_region &rg, int extrapolation,
   const mesh_region &rg_source, size_type nbdof) {

    ga_workspace workspace(md);
    workspace.add_interpolation_expression(expr, mti.linked_mesh(), rg);

    mti.distribute(extrapolation, rg_source);
    ga_interpolation_context_mti gic(mti, result, nbdof);
    ga_interpolation(workspace, gic);
  }


  // Interpolation on a im_data
  struct ga_interpolation_context_im_data
    : public ga_interpolation_context {
    base_vector &result;
    const im_data &imd;
    bool initialized;
    size_type s;

    virtual const bgeot::stored_point_tab &
    points_for_element(size_type cv, short_type f,
                       std::vector<size_type> &ind) const {
      pintegration_method pim =imd.linked_mesh_im().int_method_of_element(cv);
      if (pim->type() == IM_NONE) return *(bgeot::pstored_point_tab(0));
      GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods cannot "
                  "be used in high level generic assembly");
      size_type i_start(0), i_end(0);
      if (f == short_type(-1))
        i_end = pim->approx_method()->nb_points_on_convex();
      else {
        i_start = pim->approx_method()->ind_first_point_on_face(f);
        i_end = i_start + pim->approx_method()->nb_points_on_face(f);
      }
      for (size_type i = i_start; i < i_end; ++i) ind.push_back(i);
      return pim->approx_method()->integration_points();
    }

    virtual bool use_pgp(size_type cv) const {
      pintegration_method pim =imd.linked_mesh_im().int_method_of_element(cv);
      if (pim->type() == IM_NONE) return false;
      GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods cannot "
                  "be used in high level generic assembly");
      return !(pim->approx_method()->is_built_on_the_fly());
    }
    virtual bool use_mim(void) const { return true; }

    virtual void store_result(size_type cv, size_type i, base_tensor &t) {
      size_type si = t.size();
      if (!initialized) {
        s = si;
        GMM_ASSERT1(imd.tensor_size() == t.sizes() ||
                    (imd.tensor_size().size() == size_type(1) &&
                     imd.tensor_size()[0] == size_type(1) &&
                     si == size_type(1)),
                    "Im_data tensor size " << imd.tensor_size() <<
                    " does not match the size of the interpolated "
                    "expression " << t.sizes() << ".");
        gmm::resize(result, s * imd.nb_filtered_index());
        gmm::clear(result);
        initialized = true;
      }
      GMM_ASSERT1(s == si, "Internal error");
      size_type ipt = imd.filtered_index_of_point(cv, i);
      GMM_ASSERT1(ipt != size_type(-1),
                  "Im data with no data on the current integration point.");
      gmm::add(t.as_vector(),
               gmm::sub_vector(result, gmm::sub_interval(s*ipt, s)));
    }

    virtual void finalize(void) {
      std::vector<size_type> data(2);
      data[0] = initialized ? result.size() : 0;
      data[1] = initialized ? s : 0;
      MPI_MAX_VECTOR(data);
      if (!initialized) {
        if (data[0]) {
          gmm::resize(result, data[0]);
          gmm::clear(result);
          s = data[1];
        } else {
           gmm::clear(result);
        }
      }
      if (initialized)
        GMM_ASSERT1(gmm::vect_size(result) == data[0] &&  s == data[1],
                    "Incompatible sizes");
      MPI_SUM_VECTOR(result);
    }

    virtual const mesh &linked_mesh(void)
    { return imd.linked_mesh_im().linked_mesh(); }

    ga_interpolation_context_im_data(const im_data &imd_, base_vector &r)
      : result(r), imd(imd_), initialized(false) { }
  };

  void ga_interpolation_im_data
  (const getfem::model &md, const std::string &expr, const im_data &imd,
   base_vector &result, const mesh_region &rg) {

    ga_workspace workspace(md);
    workspace.add_interpolation_expression
      (expr, imd.linked_mesh_im(), rg);

    ga_interpolation_context_im_data gic(imd, result);
    ga_interpolation(workspace, gic);
  }


  //=========================================================================
  // Interpolate transformation with an expression
  //=========================================================================

  class  interpolate_transformation_expression
    : public virtual_interpolate_transformation, public context_dependencies {

    const mesh &source_mesh;
    const mesh &target_mesh;
    std::string expr;
    mutable bgeot::rtree element_boxes;
    mutable bool recompute_elt_boxes;
    mutable ga_workspace local_workspace;
    mutable ga_instruction_set local_gis;
    mutable bgeot::geotrans_inv_convex gic;
    mutable base_node P;
    mutable std::set<var_trans_pair> used_vars;
    mutable std::set<var_trans_pair> used_data;
    mutable std::map<var_trans_pair,
                     std::pair<ga_workspace,
                               ga_instruction_set> > compiled_derivatives;
    mutable bool extract_variable_done;
    mutable bool extract_data_done;

  public:
    void update_from_context(void) const {
      recompute_elt_boxes = true;
    }

    void extract_variables(const ga_workspace &workspace,
                           std::set<var_trans_pair> &vars,
                           bool ignore_data, const mesh &/* m */,
                           const std::string &/* interpolate_name */) const {
      if ((ignore_data && !extract_variable_done) ||
          (!ignore_data && !extract_data_done)) {
        used_vars.clear();
        ga_workspace aux_workspace;
        aux_workspace = ga_workspace(true, workspace);
        aux_workspace.clear_expressions();
        aux_workspace.add_interpolation_expression(expr, source_mesh);
        for (size_type i = 0; i < aux_workspace.nb_trees(); ++i)
          ga_extract_variables(aux_workspace.tree_info(i).ptree
                               ->root, aux_workspace, source_mesh,
                               ignore_data ? used_vars : used_data,
                               ignore_data);
        if (ignore_data)
          extract_variable_done = true;
        else
          extract_data_done = true;
      }
      if (ignore_data)
        vars.insert(used_vars.begin(), used_vars.end());
      else
        vars.insert(used_data.begin(), used_data.end());
    }

    void init(const ga_workspace &workspace) const {
      size_type N = target_mesh.dim();

      // Expression compilation
      local_workspace = ga_workspace(true, workspace);
      local_workspace.clear_expressions();

      local_workspace.add_interpolation_expression(expr, source_mesh);
      local_gis = ga_instruction_set();
      ga_compile_interpolation(local_workspace, local_gis);

      // In fact, transformations are not allowed  ... for future compatibility
      std::set<std::string>::iterator iti = local_gis.transformations.begin();
      for (; iti != local_gis.transformations.end(); ++iti)
        (local_workspace.interpolate_transformation(*iti))
          ->init(local_workspace);

      if (!extract_variable_done) {
        std::set<var_trans_pair> vars;
        extract_variables(workspace, vars, true, source_mesh, "");
      }

      for (std::set<var_trans_pair>::iterator it = used_vars.begin();
           it != used_vars.end(); ++it) {
        std::pair<ga_workspace, ga_instruction_set>
          &pwi = compiled_derivatives[*it];
        pwi.first = local_workspace;
        pwi.second = ga_instruction_set();
        if (pwi.first.nb_trees()) {
          ga_tree tree = *(pwi.first.tree_info(0).ptree);
          ga_derivative(tree, pwi.first, source_mesh,
                        it->first, it->second, 1);
          if (tree.root)
            ga_semantic_analysis(expr, tree, local_workspace, 1, 1,
                                 false, true);
          ga_compile_interpolation(pwi.first, pwi.second);
        }
      }

      // Element_boxes update (if necessary)
      if (recompute_elt_boxes) {

        element_boxes.clear();
        base_node bmin(N), bmax(N);
        for (dal::bv_visitor cv(target_mesh.convex_index());
             !cv.finished(); ++cv) {

          bgeot::pgeometric_trans pgt = target_mesh.trans_of_convex(cv);

          size_type nbd_t = pgt->nb_points();
          if (nbd_t) {
            gmm::copy(target_mesh.points_of_convex(cv)[0], bmin);
            gmm::copy(bmin, bmax);
          } else {
            gmm::clear(bmin);
            gmm::clear(bmax);
          }
          for (short_type ip = 1; ip < nbd_t; ++ip) {
            // size_type ind = target_mesh.ind_points_of_convex(cv)[ip];
            const base_node &pt = target_mesh.points_of_convex(cv)[ip];

            for (size_type k = 0; k < N; ++k) {
              bmin[k] = std::min(bmin[k], pt[k]);
              bmax[k] = std::max(bmax[k], pt[k]);
            }
          }

          scalar_type h = bmax[0] - bmin[0];
          for (size_type k = 1; k < N; ++k) h = std::max(h, bmax[k]-bmin[k]);
          if (pgt->is_linear()) h *= 1E-8;
          for (size_type k = 0; k < N; ++k)
            { bmin[k] -= h*0.2; bmax[k] += h*0.2; }

          element_boxes.add_box(bmin, bmax, cv);
        }
        recompute_elt_boxes = false;
      }
    }

    void finalize(void) const {
      std::set<std::string>::iterator iti = local_gis.transformations.begin();
      for (; iti != local_gis.transformations.end(); ++iti)
        (local_workspace.interpolate_transformation(*iti))->finalize();
      local_gis = ga_instruction_set();
    }


    int transform(const ga_workspace &/*workspace*/, const mesh &m,
                  fem_interpolation_context &ctx_x,
                  const base_small_vector &Normal,
                  const mesh **m_t,
                  size_type &cv, short_type &face_num, base_node &P_ref,
                  base_small_vector &/*N_y*/,
                  std::map<var_trans_pair, base_tensor> &derivatives,
                  bool compute_derivatives) const {
      int ret_type = 0;

      ga_interpolation_single_point_exec(local_gis, local_workspace, ctx_x,
                                         Normal, m);

      GMM_ASSERT1(local_workspace.assembled_tensor().size() == m.dim(),
                  "Wrong dimension of the tranformation expression");
      P.resize(m.dim());
      gmm::copy(local_workspace.assembled_tensor().as_vector(), P);

      bgeot::rtree::pbox_set bset;
      element_boxes.find_boxes_at_point(P, bset);
      *m_t = &target_mesh;

      while (bset.size()) {
        bgeot::rtree::pbox_set::iterator it = bset.begin(), itmax = it;

        if (bset.size() > 1) {
          // Searching the box for which the point is the most in the interior
          scalar_type rate_max = scalar_type(-1);
          for (; it != bset.end(); ++it) {

            scalar_type rate_box = scalar_type(1);
            for (size_type i = 0; i < m.dim(); ++i) {
              scalar_type h = (*it)->max[i] - (*it)->min[i];
              if (h > scalar_type(0)) {
                scalar_type rate
                  = std::min((*it)->max[i] - P[i], P[i] - (*it)->min[i]) / h;
                rate_box = std::min(rate, rate_box);
              }
            }
            if (rate_box > rate_max) {
              itmax = it;
              rate_max = rate_box;
            }
          }
        }

        cv = (*itmax)->id;
        gic.init(target_mesh.points_of_convex(cv),
                 target_mesh.trans_of_convex(cv));

        bool converged = true;
        bool is_in = gic.invert(P, P_ref, converged, 1E-4);
        // cout << "cv = " << cv << " P = " << P << " P_ref = " << P_ref << endl;
        // cout << " is_in = " << int(is_in) << endl;
        // for (size_type iii = 0;
        //     iii < target_mesh.points_of_convex(cv).size(); ++iii)
        //  cout << target_mesh.points_of_convex(cv)[iii] << endl;

        if (is_in && converged) {
          face_num = short_type(-1); // Should detect potential faces ?
          ret_type = 1;
          break;
        }

        if (bset.size() == 1) break;
        bset.erase(itmax);
      }

      // Note on derivatives of the transformation : for efficiency and
      // simplicity reasons, the derivative should be computed with
      // the value of corresponding test functions. This means that
      // for a transformation F(u) the computed derivative is F'(u).Test_u
      // including the Test_u.
      if (compute_derivatives) { // To be tested both with the computation
                                 // of derivative. Could be optimized ?
        for (std::map<var_trans_pair, base_tensor>::iterator
             itd = derivatives.begin(); itd != derivatives.end(); ++itd) {
          std::pair<ga_workspace, ga_instruction_set>
            &pwi = compiled_derivatives[itd->first];

          gmm::clear(pwi.first.assembled_tensor().as_vector());
          ga_function_exec(pwi.second);
          itd->second = pwi.first.assembled_tensor();
        }
      }
      return ret_type;
    }

    interpolate_transformation_expression(const mesh &sm, const mesh &tm,
                                          const std::string &expr_)
      : source_mesh(sm), target_mesh(tm), expr(expr_),
        recompute_elt_boxes(true), extract_variable_done(false),
        extract_data_done(false)
    { this->add_dependency(tm); }

  };


  void add_interpolate_transformation_from_expression
  (model &md, const std::string &name, const mesh &sm, const mesh &tm,
   const std::string &expr) {
    pinterpolate_transformation p
      = new interpolate_transformation_expression(sm, tm, expr);

    md.add_interpolate_transformation(name, p);
  }

  void add_interpolate_transformation_from_expression
  (ga_workspace &workspace, const std::string &name, const mesh &sm,
   const mesh &tm, const std::string &expr) {
    pinterpolate_transformation p
      = new interpolate_transformation_expression(sm, tm, expr);

    workspace.add_interpolate_transformation(name, p);
  }

  //=========================================================================
  // Interpolate transformation on neighbour element (for internal faces)
  //=========================================================================


  struct gauss_pt_corresp {
    bgeot::pgeometric_trans pgt1, pgt2;
    papprox_integration pai;
    std::vector<size_type> nodes;
  };

  bool operator <(const gauss_pt_corresp &gpc1,
		  const gauss_pt_corresp &gpc2) {
    if (gpc1.pai != gpc2.pai)
      { if (gpc1.pai  <  gpc2.pai ) return true; else return false; }
    if (gpc1.nodes.size() !=  gpc2.nodes.size()) {
      if (gpc1.nodes.size() < gpc2.nodes.size())
	return true; else return false;
    }
    for (size_type i = 0; i < gpc1.nodes.size(); ++i) {
      if (gpc1.nodes[i] != gpc2.nodes[i])
	{ if (gpc1.nodes[i] < gpc2.nodes[i]) return true; else return false; }
    }
    if (gpc1.pgt1 != gpc2.pgt1)
      { if (gpc1.pgt1 <  gpc2.pgt1) return true; else return false; }
    if (gpc1.pgt2 !=  gpc2.pgt2)
      { if (gpc1.pgt2 <  gpc2.pgt2) return true; else return false; }
    return false;
  }

  std::map<gauss_pt_corresp, bgeot::pstored_point_tab> stored_corresp; // a mettre dans une structure intermédiaire ?
  

  static bgeot::pstored_point_tab furnish_pspt(const mesh_im &mim,
						size_type cv,
						short_type f) {
    const mesh& m = mim.linked_mesh();
    GMM_ASSERT1(f != short_type(-1) && cv != size_type (-1)
		&& m.convex_index().is_in(cv), "Invalid convex or face id");
    auto adj_face = m.adjacent_face(cv, f);
    GMM_ASSERT1(adj_face.cv != size_type(-1), "No adjacent face");
    pintegration_method pim = mim.int_method_of_element(cv);
    GMM_ASSERT1(pim->type() == IM_APPROX, "Unvalid integration method");
    papprox_integration pai = pim->approx_method();
    GMM_ASSERT1(!(pai->is_built_on_the_fly()), "Do not call "
		"this function for build on the fly integration methods");


    // Fill a gauss_pt_corresp structure.
    gauss_pt_corresp gpc;
    gpc.pgt1 = m.trans_of_convex(cv);
    gpc.pgt2 = m.trans_of_convex(adj_face.cv);
    gpc.pai = pai;
    auto inds_pt1 = m.ind_points_of_face_of_convex(cv, f);
    auto inds_pt2 = m.ind_points_of_face_of_convex(adj_face.cv, adj_face.f);
    auto str1 = gpc.pgt1->structure();
    auto str2 = gpc.pgt2->structure();
    size_type nbptf1 = str1->nb_points_of_face(f);
    size_type nbptf2 = str2->nb_points_of_face(adj_face.f);
    gpc.nodes.resize(nbptf1*2);
    for (size_type i = 0; i < nbptf1; ++i)  {
      gpc.nodes[2*i] = str1->ind_points_of_face(f)[i];
      bool found = false;
      for (size_type j = 0; j < nbptf2; ++j) {
	if (inds_pt2[j] == inds_pt1[i]) {
	  gpc.nodes[2*i+1] = str2->ind_points_of_face(adj_face.f)[j];
	  found = true;
	  break;
	}
      }
      GMM_ASSERT1(found, "Internal error");
    }
    
    auto itm = stored_corresp.find(gpc);
    if (itm != stored_corresp.end()) return itm->second;
    else {
      size_type nbpt = pai->nb_points_on_face(f);
      base_matrix G;
      bgeot::geotrans_inv_convex gic;
      gic.init(m.points_of_convex(adj_face.cv), gpc.pgt2);
      size_type first_ind = pai->ind_first_point_on_face(f);
      const bgeot::stored_point_tab &spt = pai->integration_points();
      bgeot::vectors_to_base_matrix(G, m.points_of_convex(cv));
      fem_interpolation_context ctx(gpc.pgt1, 0, spt[0], G, cv, f);
      std::vector<base_node> P_ref(nbpt);
      
      for (size_type i = 0; i < nbpt; ++i) {
	ctx.set_xref(spt[first_ind+i]);
	bool converged = true;
	bool is_in = gic.invert(ctx.xreal(), P_ref[i], converged, 1E-4);
	GMM_ASSERT1(is_in && converged, "Geometric transformation inversion "
		    "has failed in neighbour transformation");
      }
      bgeot::pstored_point_tab pspt = store_point_tab(P_ref);
      stored_corresp[gpc] = pspt;
      return pspt;
    }
  }


  class  interpolate_transformation_neighbour
    : public virtual_interpolate_transformation, public context_dependencies {

  public:
    void update_from_context(void) const {}
    void extract_variables(const ga_workspace &/* workspace */,
                           std::set<var_trans_pair> &/* vars */,
                           bool /* ignore_data */, const mesh &/* m */,
                           const std::string &/* interpolate_name */) const {}
    void init(const ga_workspace &/* workspace */) const {}
    void finalize(void) const {}

    int transform(const ga_workspace &/*workspace*/, const mesh &m_x,
                  fem_interpolation_context &ctx_x,
                  const base_small_vector &/*Normal*/, const mesh **m_t,
                  size_type &cv, short_type &face_num, base_node &P_ref,
                  base_small_vector &/*N_y*/,
                  std::map<var_trans_pair, base_tensor> &/*derivatives*/,
                  bool compute_derivatives) const {

      int ret_type = 0;
      *m_t = &m_x;
      size_type cv_x = ctx_x.convex_num();
      short_type face_x = ctx_x.face_num();
      GMM_ASSERT1(face_x != short_type(-1), "Neighbour transformation can "
                  "only be applied to internal faces");

      auto adj_face = m_x.adjacent_face(cv_x, face_x);

      if (adj_face.cv != size_type(-1)) {
        bgeot::geotrans_inv_convex gic;
        gic.init(m_x.points_of_convex(adj_face.cv),
                 m_x.trans_of_convex(adj_face.cv));
        bool converged = true;
        bool is_in = gic.invert(ctx_x.xreal(), P_ref, converged, 1E-4);
        GMM_ASSERT1(is_in && converged, "Geometric transformation inversion "
                    "has failed in neighbour transformation");
        face_num = adj_face.f;
        cv = adj_face.cv;
        ret_type = 1;
      }
      GMM_ASSERT1(!compute_derivatives,
                  "No derivative for this transformation");
      return ret_type;
    }

    interpolate_transformation_neighbour(void) { }

  };

  // Should be added by default

  void add_interpolate_transformation_neighbour(model &md) {
    pinterpolate_transformation p = new interpolate_transformation_neighbour();
    md.add_interpolate_transformation("neighbour_elt", p);
  }

  void add_interpolate_transformation_neighbour(ga_workspace &workspace) {
    pinterpolate_transformation p = new interpolate_transformation_neighbour();
    workspace.add_interpolate_transformation("neighbour_elt", p);
  }


} /* end of namespace */
