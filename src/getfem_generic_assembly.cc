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
 
===========================================================================*/

/** 
  Provising for special Math functions unavailable on Intel or MSVS C++
  compilers
*/
#ifdef _WIN32
#include <limits>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/erf.hpp>
#endif

#include "getfem/getfem_generic_assembly.h"
#include "gmm/gmm_blas.h"
#include <iomanip>
#include "getfem/getfem_omp.h"



extern "C" void daxpy_(const int *n, const double *alpha, const double *x,
		       const int *incx, double *y, const int *incy);


namespace getfem {
  
  //=========================================================================
  // Lexical analysis for the generic assembly langage
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
    GA_TRACE,       // 'Trace' Trace operator
    GA_PRINT,       // 'Print' Print the tensor
    GA_DOT,         // '.'
    GA_DOTMULT,     // '.*' componentwize multiplication
    GA_DOTDIV,      // './' componentwize division
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
    ga_operator_priorities[GA_TRACE] = 3;
    ga_operator_priorities[GA_PRINT] = 3;
    ga_operator_priorities[GA_UNARY_MINUS] = 3;

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
    if (expr.size()) {
      int first = std::max(0, int(pos)-40);
      int last = std::min(int(pos)+20, int(expr.size()));
      if (last - first < 60)
      first = std::max(0, int(pos)-40-(60-last+first));
      if (last - first < 60)
        last = std::min(int(pos)+20+(60-last+first),int(expr.size()));
      
      if (first > 0) cerr << "...";
      cerr << expr.substr(first, last-first);
      if (last < int(expr.size())) cerr << "...";
      cerr << endl;
      if (first > 0) cerr << "   ";
      if (int(pos) > first) cerr << std::setw(int(pos)-first) << ' ';
      cerr << '|' << endl;
    }
    cerr << msg << endl;
  }

#define ga_throw_error(expr, pos, msg)               \
  { ga_throw_error_msg(expr, pos, msg);              \
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
    GA_NODE_ALLINDICES,
    GA_NODE_C_MATRIX,
    GA_NODE_VAL,
    GA_NODE_GRAD,
    GA_NODE_HESS,
    GA_NODE_TEST,
    GA_NODE_GRAD_TEST,
    GA_NODE_HESS_TEST,
    GA_NODE_ZERO};

  struct ga_tree_node;
  typedef ga_tree_node *pga_tree_node;

  struct ga_tree_node {
    GA_NODE_TYPE node_type;
    base_tensor t;
    size_type test_function_type; // -1 = undetermined
                                  // 0 = no test function,
                                  // 1 = first order, 2 = second order,
                                  // 3 = both with first order in first
                                  // 4 = both with first order in second
    std::string name_test1, name_test2; // variable names corresponding to test
                                  // functions when test_function_type > 0.
    size_type qdim1, qdim2;       // Qdims when test_function_type > 0.
    size_type nbc1, nbc2, nbc3;   // For explicit matrices.
    size_type pos;                // Position of the first character in string
    std::string name;             // variable/constant/function/operator name 
    size_type der1, der2;         // For functions and nonlinear operators,
                                  // optional derivative or second derivative.
    GA_TOKEN_TYPE op_type;
    pga_tree_node parent;         // Parent node
    std::vector<pga_tree_node> children; // Children nodes
    scalar_type hash_value;       // Hash value to identify nodes.
    bool marked;                  // For specific use of some algorithms
    
    inline size_type nb_test_functions(void) const {
      if (test_function_type == size_type(-1)) return 0;
      return (test_function_type ? 1 : 0) + (test_function_type >= 3 ? 1 : 0);
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
                   const std::string &expr, const ga_workspace &workspace) {
      size_type test0 = n0->test_function_type, test1 = n1->test_function_type;
      if (test0 && test1 && (test0 == test1 ||
                             test0 >= 3 || test1 >= 3))
        ga_throw_error(expr, pos, "Incompatibility of test functions "
                        " in product.");
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
      
      if (n0->name_test1.size())
        { name_test1 = n0->name_test1; qdim1 = n0->qdim1; }
      else
        { name_test1 = n1->name_test1; qdim1 = n1->qdim1; }
      
      if (n0->name_test2.size())
        { name_test2 = n0->name_test2; qdim2 = n0->qdim2; }
      else
        { name_test2 = n1->name_test2; qdim2 = n1->qdim2; }

      if (test1 == 2 && test0 &&
          workspace.associated_mf(name_test1)
          != workspace.associated_mf(name_test2)) {
        test_function_type = 4;
        std::swap(mi[0], mi[1]);
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
    
    void init_scalar_tensor(scalar_type v) {
      t.adjust_sizes(bgeot::multi_index());
      t[0] = v;
      test_function_type = 0;
    }
    void init_vector_tensor(size_type d) {
      bgeot::multi_index mi(1);
      mi[0]=d; t.adjust_sizes(mi);
      test_function_type = 0;
    }
    void init_matrix_tensor(size_type n, size_type m) {
      t.adjust_sizes(bgeot::multi_index(n,m));
      test_function_type = 0;
    }
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
      : node_type(GA_NODE_VOID), test_function_type(-1), qdim1(0), qdim2(0),
        der1(0), der2(0), hash_value(0) {}
    ga_tree_node(GA_NODE_TYPE ty, size_type p)
      : node_type(ty), test_function_type(-1), qdim1(0), qdim2(0),
        pos(p), der1(0), der2(0), hash_value(0) {}
    ga_tree_node(scalar_type v, size_type p)
      : node_type(GA_NODE_CONSTANT), test_function_type(-1), qdim1(0),
        qdim2(0), pos(p), der1(0), der2(0), hash_value(0)
    { init_scalar_tensor(v); }
    ga_tree_node(const char *n, size_type l, size_type p)
      : node_type(GA_NODE_NAME), test_function_type(-1), qdim1(0), qdim2(0),
        pos(p), name(n, l), der1(0), der2(0), hash_value(0) {}
    ga_tree_node(GA_TOKEN_TYPE op, size_type p)
      : node_type(GA_NODE_OP), test_function_type(-1), qdim1(0), qdim2(0),
        pos(p), der1(0), der2(0), op_type(op), hash_value(0) {}
    
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
                           current_node->node_type == GA_NODE_C_MATRIX)) {
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
    
    void add_op(GA_TOKEN_TYPE op_type, size_type pos) {
      while (current_node && current_node->parent &&
             current_node->parent->node_type == GA_NODE_OP &&
             ga_operator_priorities[current_node->parent->op_type]
             >= ga_operator_priorities[op_type])
        current_node = current_node->parent;
      pga_tree_node new_node = new ga_tree_node(op_type, pos);
      if (current_node) {
        if (op_type == GA_UNARY_MINUS || op_type == GA_TRACE
            || op_type == GA_PRINT) {
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

    void clear(void)
    { if (root) clear_node_rec(root); root = current_node = 0; }
    
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

    void duplicate_with_addition(pga_tree_node pnode) {
      pga_tree_node plus = new ga_tree_node(GA_PLUS, pnode->pos);
      plus->children.resize(2);
      plus->children[0] = pnode;
      plus->parent = pnode->parent;
      if (pnode->parent) {
        for (size_type j = 0; j < pnode->parent->children.size(); ++j)
          if (pnode->parent->children[j] == pnode)
            pnode->parent->children[j] = plus;
      } else root = plus; 
      pnode->parent = plus;
      copy_node(pnode, plus, plus->children[1]);
    }

    void insert_node(pga_tree_node pnode) {
      pga_tree_node newnode = new ga_tree_node;
      newnode->parent = pnode->parent;
      if (pnode->parent) {
        for (size_type j = 0; j < pnode->parent->children.size(); ++j)
          if (pnode->parent->children[j] == pnode)
            pnode->parent->children[j] = newnode;
      } else root = newnode;
      newnode->children.push_back(pnode);
      pnode->parent = newnode;
    }

    void add_children(pga_tree_node pnode) {
      pga_tree_node newnode = new ga_tree_node;
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


  static bool sub_tree_are_equal(pga_tree_node pnode1, pga_tree_node pnode2) {
    if (pnode1->node_type != pnode2->node_type) return false; // TODO: Cas des constant nulles
    if (pnode1->children.size() != pnode2->children.size()) return false;

    switch(pnode1->node_type) {
    case GA_NODE_OP:
      if (pnode1->op_type != pnode2->op_type) return false;
      break;
    case GA_NODE_OPERATOR:
      if (pnode1->der1 != pnode2->der1 || pnode1->der2 != pnode2->der2)
        return false;
    case GA_NODE_PREDEF_FUNC: case GA_NODE_SPEC_FUNC:
      if (pnode1->name.compare(pnode2->name)) return false;
      break;
    case GA_NODE_CONSTANT: case GA_NODE_ZERO:
      if (pnode1->t.size() != pnode2->t.size()) return false;
      if (pnode1->t.size() != 1 &&
          pnode1->t.sizes().size() != pnode2->t.sizes().size()) return false;
      for (size_type i = 0; i < pnode1->t.sizes().size(); ++i)
        if (pnode1->t.sizes()[i] != pnode2->t.sizes()[i]) return false;
      for (size_type i = 0; i < pnode1->t.size(); ++i)
        if (gmm::abs(pnode1->t[i] - pnode2->t[i]) > 1E-25) return false;
      break;
    case GA_NODE_C_MATRIX:
      if (pnode1->nbc1 != pnode2->nbc1 || pnode1->nbc2 != pnode2->nbc2
          ||   pnode1->nbc3 != pnode2->nbc3)  return false;
      break;
    case GA_NODE_TEST: case GA_NODE_GRAD_TEST: case GA_NODE_HESS_TEST:
      if (pnode2->test_function_type != pnode2->test_function_type)
        return false;
    case GA_NODE_VAL: case GA_NODE_GRAD: case GA_NODE_HESS:
      if (pnode1->name.compare(pnode2->name)) return false;
      break;

    default:break;
    }
    for (size_type i = 0; i < pnode1->children.size(); ++i)
      if (!(sub_tree_are_equal(pnode1->children[i], pnode2->children[i])))
          return false;
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

  static void ga_print_constant_tensor(pga_tree_node pnode) {
    size_type nt = pnode->nb_test_functions(); // for printing zero tensors
    switch (pnode->tensor_order()) {
    case 0:
      cout << (nt ? scalar_type(0) : pnode->t[0]);
      break;
    case 1:
      cout << "[";
      for (size_type i = 0; i < pnode->tensor_proper_size(0); ++i) {
        if (i != 0) cout << "; ";
        cout << (nt ? scalar_type(0) : pnode->t[i]);
      }
      cout << "]";
      break;
    case 2:
      cout << "[";
      for (size_type i = 0; i < pnode->tensor_proper_size(0); ++i) {
        if (i != 0) cout << "; ";
        for (size_type j = 0; j < pnode->tensor_proper_size(1); ++j) {
          if (j != 0) cout << ", ";
          cout << (nt ? scalar_type(0) : pnode->t(i,j));
        }
      }
      cout << "]";
      break;
    case 3:
      cout << "[";
      for (size_type i = 0; i < pnode->tensor_proper_size(0); ++i) {
        if (i != 0) cout << ",, ";
        for (size_type j = 0; j < pnode->tensor_proper_size(1); ++j) {
          if (j != 0) cout << "; "; 
          for (size_type k = 0; k < pnode->tensor_proper_size(2); ++k) {
            if (k != 0) cout << ", "; 
            cout << (nt ? scalar_type(0) : pnode->t(i,j,k));
          }
        }
      }
      cout << "]";
      break;
    case 4:
      cout << "[";
      for (size_type i = 0; i < pnode->tensor_proper_size(0); ++i) {
        if (i != 0) cout << ";; ";
        for (size_type j = 0; j < pnode->tensor_proper_size(1); ++j) {
          if (j != 0) cout << ",, "; 
          for (size_type k = 0; k < pnode->tensor_proper_size(2); ++k) {
            if (k != 0) cout << "; "; 
            for (size_type l = 0; l < pnode->tensor_proper_size(3); ++l) {
              if (l != 0) cout << ", ";
              cout << (nt ? scalar_type(0) : pnode->t(i,j,k,l));
            }
          }
        }
      }
      cout << "]";
      break;
    default: GMM_ASSERT1(false, "Invalid tensor dimension");
    }
    GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
  }

  static void ga_print_node(pga_tree_node pnode) {
    
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

        
        if (par) cout << "(";
        if (pnode->op_type == GA_UNARY_MINUS) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          cout << "-"; ga_print_node(pnode->children[0]); 
        } else if (pnode->op_type == GA_QUOTE) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          ga_print_node(pnode->children[0]); cout << "'";
        } else if (pnode->op_type == GA_TRACE) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          cout << "Trace("; ga_print_node(pnode->children[0]); cout << ")";
        } else if (pnode->op_type == GA_PRINT) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          cout << "Print("; ga_print_node(pnode->children[0]); cout << ")";
        } else {
          GMM_ASSERT1(pnode->children.size() == 2, "Invalid tree");
          if (pnode->op_type == GA_MULT &&
              (pnode->test_function_type == size_type(-1) ||
               (pnode->children[0]->tensor_order() == 4 &&
                pnode->children[1]->tensor_order() == 2)))
            { par = true; cout << "("; }
          ga_print_node(pnode->children[0]);
          switch (pnode->op_type) {
          case GA_PLUS: cout << "+"; break;
          case GA_MINUS: cout << "-"; break;
          case GA_MULT: cout << "*"; break;
          case GA_DIV: cout << "/"; break;
          case GA_COLON: cout << ":"; break;
          case GA_DOT: cout << "."; break;
          case GA_DOTMULT: cout << ".*"; break;
          case GA_DOTDIV: cout << "./"; break;
          case GA_TMULT: cout << "@"; break;
          default: GMM_ASSERT1(false, "Invalid or not taken into account "
                               "operation");
          }
          ga_print_node(pnode->children[1]);
        }
        if (par) cout << ")";
      }
      break;
      
    case GA_NODE_VAL: cout << pnode->name; break;
    case GA_NODE_GRAD: cout << "Grad_" << pnode->name; break;
    case GA_NODE_HESS: cout << "Hess_" << pnode->name; break;
    case GA_NODE_TEST: 
      if (pnode->test_function_type == 1) cout << "Test_" << pnode->name;
      else cout << "Test2_" << pnode->name;
      break;
    case GA_NODE_GRAD_TEST:
      if (pnode->test_function_type == 1) cout << "Grad_Test_" << pnode->name;
      else cout << "Grad_Test2_" << pnode->name;
      break;
    case GA_NODE_HESS_TEST:
      if (pnode->test_function_type == 1) cout << "Hess_Test_" << pnode->name;
      else cout << "Hess_Test2_" << pnode->name;
      break;
    case GA_NODE_SPEC_FUNC: cout << pnode->name; break;
    case GA_NODE_OPERATOR:
    case GA_NODE_PREDEF_FUNC:
      if (pnode->der1) {
        cout << "Derivative_" << pnode->der1 << "_";
        if (pnode->der2) cout << pnode->der2 << "_";
      }
      cout << pnode->name; break;
    case GA_NODE_ZERO:
      GMM_ASSERT1(pnode->test_function_type != size_type(-1),
                  "Internal error");
      if (pnode->test_function_type) cout << "(";
      ga_print_constant_tensor(pnode);
      if (pnode->name_test1.size()) {
        GMM_ASSERT1(pnode->qdim1 > 0, "Internal error");
        if (pnode->qdim1 == 1)
          cout << "*Test_" << pnode->name_test1;
        else {
          cout << "*([ 0";
          for (size_type i = 1; i < pnode->qdim1; ++i) cout << ", 0";
          cout << "].Test_" << pnode->name_test1 << ")";
        }
      }
      if (pnode->name_test2.size()) {
        GMM_ASSERT1(pnode->qdim2 > 0, "Internal error");
        if (pnode->qdim2 == 1)
          cout << "*Test2_" << pnode->name_test2;
        else {
          cout << "*([ 0";
          for (size_type i = 1; i < pnode->qdim2; ++i) cout << ", 0";
          cout << "].Test2_" << pnode->name_test2 << ")";
        }
      }
      if (pnode->test_function_type) cout << ")";
      break;

    case GA_NODE_CONSTANT: ga_print_constant_tensor(pnode); break;

    case GA_NODE_ALLINDICES:
      cout << ":";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_PARAMS:
      GMM_ASSERT1(pnode->children.size(), "Invalid tree");
      ga_print_node(pnode->children[0]);
      cout << "(";
      for (size_type i = 1; i < pnode->children.size(); ++i)
        { if (i > 1) cout << ", "; ga_print_node(pnode->children[i]); }
      cout << ")";
      break;

    case GA_NODE_NAME:
      cout << pnode->name;
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_C_MATRIX:
      GMM_ASSERT1(pnode->children.size(), "Invalid tree");
      cout << "[";
      for (size_type i = 0; i < pnode->children.size(); ++i) {
        if (i > 0) {
          if (i%pnode->nbc1 != 0) cout << ", ";
          else {
            if (pnode->nbc2 > 1 || pnode->nbc3 > 1) {
              if (i%(pnode->nbc1*pnode->nbc2) != 0) cout << "; ";
              else if (i%(pnode->nbc1*pnode->nbc2*pnode->nbc3) != 0)
                cout << ",, ";
              else cout << ";; ";
            } else cout << "; ";
          }
        }
        ga_print_node(pnode->children[i]);
      }
      cout << "]";
      break;

    default:
      cout << "Invalid or not taken into account node type "
           << pnode->node_type;
      break;
    }
  }
 

  static void ga_print_tree(const ga_tree &tree) {
    size_type pr = cout.precision(16); // Not exception safe
    if (tree.root) verify_tree(tree.root, 0);
    if (tree.root) ga_print_node(tree.root); else cout << "Empty tree";
    cout << endl;
    cout.precision(pr);
  }

  static void ga_print_hash_value_node(pga_tree_node pnode) {
    cout << "node_type = " << pnode->node_type << " hash value = "
         << pnode->hash_value << endl;
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_print_hash_value_node(pnode->children[i]);
  }

  static void ga_print_hash_values(const ga_tree &tree) {
    size_type pr = cout.precision(16); // Not exception safe
    if (tree.root) ga_print_hash_value_node(tree.root);
    else cout << "Empty tree";
    cout << endl;
    cout.precision(pr);
  }


  //=========================================================================
  // Syntax analysis for the generic assembly langage
  //=========================================================================
             
  // Read a term with an (implicit) pushdown automaton.
  static GA_TOKEN_TYPE ga_read_term(const std::string &expr,
                                    size_type &pos, ga_tree &tree) {
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
            size_type nbc1 = 0, nbc2 = 0, nbc3 = 0, n1 = 0, n2 = 0, n3 = 0;
            bool foundsemi = false, founddcomma = false, founddsemi = false;
            tree.add_matrix(token_pos);
            for(;;) {
              r_type = ga_read_term(expr, pos, sub_tree);
              ++n1; ++n2; ++n3;
              if (!foundsemi) ++nbc1;
              if (!founddcomma) ++nbc2;
              if (!founddsemi) ++nbc3;

              switch(r_type) {
              case GA_COMMA: break;
              case GA_SEMICOLON: foundsemi = true; n1 = 0; break;
              case GA_DCOMMA: founddcomma = true; n2 = 0; n1 = 0; break;
              case GA_DSEMICOLON:
                founddsemi = true; n3 = 0; n2 = 0; n1 = 0; break;
              case GA_RBRACKET:
                if (n1 != nbc1 || n2 != nbc2 || n3 != nbc3 ||
                    (founddcomma && !founddsemi) ||
                    (!founddcomma && founddsemi))
                  ga_throw_error(expr, pos-1, "Bad explicit "
                                  "vector/matrix/tensor format. ");
                tree.current_node->nbc1 = nbc1;
                if (founddcomma) {
                  tree.current_node->nbc2 = nbc2/nbc1;
                  tree.current_node->nbc3 = nbc3/nbc2;
                } else {
                  tree.current_node->nbc2 = tree.current_node->nbc3 = 1;
                }
                break;
              default:
                ga_throw_error(expr, pos-1, "The explicit "
                                "vector/matrix/tensor components should be "
                                "separated by ',', ';', ',,' and ';;' and "
                                "be ended by ']'.");
                break;
              }
         
              tree.add_sub_tree(sub_tree);
              if (r_type == GA_RBRACKET) break;
            }
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
            for(;;) {
              r_type = ga_read_term(expr, pos, sub_tree);
              if (r_type != GA_RPAR && r_type != GA_COMMA)
                ga_throw_error(expr, pos-((r_type != GA_END)?1:0),
                               "Parameters should be separated "
                               "by ',' and parameter list ended by ')'.");
              tree.add_sub_tree(sub_tree);
              if (r_type == GA_RPAR) break;
            }
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
  // Structure dealing with predefined special functions
  // such as mesh_dim, pi, qdim ...
  //=========================================================================

  typedef std::set<std::string> ga_spec_function_tab;
  static ga_spec_function_tab SPEC_FUNCTIONS;

  //=========================================================================
  // Structure dealing with predefined scalar functions.
  //=========================================================================

  static scalar_type ga_Heaveside(scalar_type t) { return (t >= 0.) ? 1.: 0.; }
  static scalar_type ga_pos_part(scalar_type t) { return (t >= 0.) ? t : 0.; }
  static scalar_type ga_neg_part(scalar_type t) { return (t >= 0.) ? 0. : -t; }
  static scalar_type ga_sqr(scalar_type t) { return t*t; }
  static scalar_type ga_max(scalar_type t, scalar_type u)
  { return std::max(t,u); }
  static scalar_type ga_min(scalar_type t, scalar_type u)
  { return std::min(t,u); }
  static scalar_type ga_abs(scalar_type t) { return gmm::abs(t); }
  static scalar_type ga_sign(scalar_type t) { return (t >= 0.) ? 1.: -1.; }
  
  // Derivatives of predefined functions
  static scalar_type ga_der_sqrt(scalar_type t) { return -0.5/sqrt(t); }
  static scalar_type ga_der_sqr(scalar_type t) { return 2*t; }
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
  static scalar_type ga_der_heaviside(scalar_type t)
  { return (t == 0) ? INFINITY : 0.; }
  static scalar_type ga_der_neg_part(scalar_type t)
  { return (t >= 0) ? 0. : -1.; }
  static scalar_type ga_der_max1(scalar_type t, scalar_type u)
  { return (t-u >= 0) ? 1. : 0.; }
  static scalar_type ga_der_max2(scalar_type t, scalar_type u)
  { return (u-t >= 0) ? 1. : 0.; }


  static ga_predef_function_tab PREDEF_FUNCTIONS;

  //=========================================================================
  // Structure dealing with predefined operators.
  //=========================================================================

  static ga_predef_operator_tab PREDEF_OPERATORS;

  void ga_init_scalar(bgeot::multi_index &mi) { mi.resize(0); }

  // Norm Operator
  struct norm_operator : public ga_nonlinear_operator {
    bool result_size(const arg_list &args, bgeot::multi_index &sizes) const {
      if (args.size() != 1 || args[0]->sizes().size() > 1) return false;
      ga_init_scalar(sizes);
      return true;
    }
    
    void value(const arg_list &args, base_tensor &result) const
    { result[0] = gmm::vect_norm2(args[0]->as_vector()); }

    // Derivative : u/|u|
    void derivative(const arg_list &args, size_type, base_tensor &result) const {
      scalar_type no = gmm::vect_norm2(args[0]->as_vector());
      if (no == scalar_type(0))
        gmm::clear(result.as_vector());
      else
        gmm::copy(gmm::scaled(args[0]->as_vector(), scalar_type(1)/no),
                  result.as_vector());
    }

    // Second derivative : (|u|^2 Id - u x u)/|u|^3
    void second_derivative(const arg_list &args, size_type, size_type,
                           base_tensor &result) const { // To be verified
      const base_tensor &t = *args[0];
      size_type N = t.size();
      scalar_type no = gmm::vect_norm2(t.as_vector());
      scalar_type no2 = no*no, no3 = no*no2;

      if (no < 1E-25) no = 1E-25; // To avoid infinite values
      
      for (size_type i = 0; i < N; ++i)
        for (size_type j = 0; j < N; ++j) {
          result[j*N+i] = - t[i]*t[j] / no3;
          if (i == j) result[j*N+i] += scalar_type(1)/no;
        }
    }
    

  };


  //=========================================================================
  // Initialization of predefined functions and operators.
  //=========================================================================

  
  bool init_predef_functions(void) {

    // Predefined functions

    ga_interval R;
    // Power functions and their derivatives
    PREDEF_FUNCTIONS["sqrt"] =
      ga_predef_function(sqrt, R, ga_interval(0, INFINITY),
                         "DER_PDFUNC_SQRT_");
    PREDEF_FUNCTIONS["sqr"] =
      ga_predef_function(ga_sqr, R, R, "DER_PDFUNC_SQR_");
    PREDEF_FUNCTIONS["pow"] =
      ga_predef_function(pow, R, R, R, R, "DER_PDFUNC_POW1_",
                         "DER_PDFUNC_POW2_");

    PREDEF_FUNCTIONS["DER_PDFUNC_SQRT_"] =
      ga_predef_function(ga_der_sqrt, "-0.5/sqrt(t)", R,
                         ga_interval(0, INFINITY));
    PREDEF_FUNCTIONS["DER_PDFUNC_SQR_"] =
      ga_predef_function(ga_der_sqr, "2*t", R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_POW1_"] =
      ga_predef_function(ga_der_pow1, "u*pow(t,u-1)", R, R, R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_POW2_"] =
      ga_predef_function(ga_der_pow2, "pow(t,u)*log(t)", R, R,
                         ga_interval(0, INFINITY), R);


    // Hyperbolic functions
    PREDEF_FUNCTIONS["exp"] = ga_predef_function(exp, R, R, "exp");
    PREDEF_FUNCTIONS["log"] =
      ga_predef_function(log, R, ga_interval(0, INFINITY), "1/t");
    PREDEF_FUNCTIONS["log10"] =
      ga_predef_function(log10, R, ga_interval(0, INFINITY),
                         "1/(t*log(10))");
    PREDEF_FUNCTIONS["sinh"] = ga_predef_function(sinh, R, R, "cosh");
    PREDEF_FUNCTIONS["cosh"] = ga_predef_function(cosh, R, R, "sinh");
    PREDEF_FUNCTIONS["tanh"] =
      ga_predef_function(tanh, R, R, "DER_PDFUNC_TANH_");
    PREDEF_FUNCTIONS["asinh"]
      = ga_predef_function(asinh, R, R, "DER_PDFUNC_ASINH_");
    PREDEF_FUNCTIONS["acosh"] =
      ga_predef_function(acosh, R, ga_interval(1, INFINITY),
                         "DER_PDFUNC_ACOSH_");
    PREDEF_FUNCTIONS["atanh"]
      = ga_predef_function(atanh, R, R,"DER_PDFUNC_ATANH_");


    PREDEF_FUNCTIONS["DER_PDFUNC_LOG_"] =
      ga_predef_function(ga_der_log, "1/t", R, ga_interval(0, INFINITY));
    PREDEF_FUNCTIONS["DER_PDFUNC_LOG10_"] =
      ga_predef_function(ga_der_log10, "1/(t*log(10))", R,
                         ga_interval(0, INFINITY));
    PREDEF_FUNCTIONS["DER_PDFUNC_TANH_"] =
      ga_predef_function(ga_der_tanh, "1-sqr(tanh(t))", R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_ASINH_"] =
      ga_predef_function(ga_der_asinh, "1/(sqrt(t*t+1))", R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_ACOSH_"] =
      ga_predef_function(ga_der_acosh, "1/(sqrt(t*t-1))", R,
                         ga_interval(1, INFINITY));
    PREDEF_FUNCTIONS["DER_PDFUNC_ATANH_"] =
      ga_predef_function(ga_der_atanh, "1/(1-t*t)", R, R);


    // Trigonometric functions
    PREDEF_FUNCTIONS["sin"] = ga_predef_function(sin, R, R, "cos");
    PREDEF_FUNCTIONS["cos"] = ga_predef_function(cos, R, R, "DER_PDFUNC_COS_");
    PREDEF_FUNCTIONS["tan"] = ga_predef_function(tan, R, R, "DER_PDFUNC_TAN_");
    PREDEF_FUNCTIONS["asin"]
      = ga_predef_function(asin, R, ga_interval(-1., 1.), "DER_PDFUNC_ASIN_");
    PREDEF_FUNCTIONS["acos"]
      = ga_predef_function(acos, R, ga_interval(-1., 1.), "DER_PDFUNC_ACOS_");
    PREDEF_FUNCTIONS["atan"]
      = ga_predef_function(atan, R, R, "DER_PDFUNC_ATAN_");
    
    PREDEF_FUNCTIONS["DER_PDFUNC_COS_"] =
      ga_predef_function(ga_der_cos, "-sin(t)", R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_TAN_"] =
      ga_predef_function(ga_der_tan, "1+sqr(tan(t))", R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_ASIN_"] =
      ga_predef_function(ga_der_asin, "1/(sqrt(1-t*t))", R,
                         ga_interval(-1., 1.));
    PREDEF_FUNCTIONS["DER_PDFUNC_ACOS_"] =
      ga_predef_function(ga_der_acos, "-1/(sqrt(1-t*t))", R,
                         ga_interval(-1., 1.));
    PREDEF_FUNCTIONS["DER_PDFUNC_ATAN_"] =
      ga_predef_function(ga_der_atan, "1/(1+t*t)", R, R);


    // Error functions
    PREDEF_FUNCTIONS["erf"]
      = ga_predef_function(erf, R, R, "DER_PDFUNC_ERF_");
    PREDEF_FUNCTIONS["erfc"]
      = ga_predef_function(erfc, R, R, "DER_PDFUNC_ERFC_");

    PREDEF_FUNCTIONS["DER_PDFUNC_ERF_"] =
      ga_predef_function(ga_der_erf, "exp(-t*t)*2/sqrt(pi)", R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_ERFC_"] =
      ga_predef_function(ga_der_erfc, "-exp(-t*t)*2/sqrt(pi)", R, R);

    

    // Miscellaneous functions
    PREDEF_FUNCTIONS["Heaveside"]
      = ga_predef_function(ga_Heaveside, ga_interval(0., INFINITY), R,
                           "DER_PDFUNC_HEAVISIDE_");
    PREDEF_FUNCTIONS["sign"]
      = ga_predef_function(ga_sign, R, R, "DER_PDFUNC_SIGN_");
    PREDEF_FUNCTIONS["abs"] = ga_predef_function(ga_abs, R, R, "sign");    
    PREDEF_FUNCTIONS["pos_part"]
      = ga_predef_function(ga_pos_part, ga_interval(0., INFINITY), R,
                           "Heaveside");
    PREDEF_FUNCTIONS["neg_part"]
      = ga_predef_function(ga_neg_part, ga_interval(-INFINITY, 0.), R,
                           "DER_PDFUNC_NEG_PART_");
    PREDEF_FUNCTIONS["max"] // TODO: the intervals could be precised
      = ga_predef_function(ga_max, R, R, R, R,
                           "DER_PDFUNC_MAX1_", "DER_PDFUNC_MAX2_");
    PREDEF_FUNCTIONS["min"] // TODO: the intervals could be precised
      = ga_predef_function(ga_min, R, R, R, R,
                           "DER_PDFUNC_MAX2_", "DER_PDFUNC_MAX1_");

    PREDEF_FUNCTIONS["DER_PDFUNC_HEAVISIDE_"] =
      ga_predef_function(ga_der_heaviside, "Dirac(0)", ga_interval(0., 0.), R);
    PREDEF_FUNCTIONS["DER_PDFUNC_SIGN_"] =
      ga_predef_function(ga_der_heaviside, "2*Dirac(0)",
                         ga_interval(0., 0.), R);
    PREDEF_FUNCTIONS["DER_PDFUNC_NEG_PART_"] =
      ga_predef_function(ga_der_neg_part, "-Heaveside(-t)",
                         ga_interval(-INFINITY, 0.), R);
    PREDEF_FUNCTIONS["DER_PDFUNC_MAX1_"] =
      ga_predef_function(ga_der_max1, "Heaveside(t-u)", R, R, R, R);
    PREDEF_FUNCTIONS["DER_PDFUNC_MAX2_"] =
      ga_predef_function(ga_der_max2, "Heaveside(u-t)", R, R, R, R);


    // Predefined special functions

    SPEC_FUNCTIONS.insert("pi");
    SPEC_FUNCTIONS.insert("meshdim");
    SPEC_FUNCTIONS.insert("qdim");
    SPEC_FUNCTIONS.insert("Id");

    // Predefined operators

    PREDEF_OPERATORS["Norm"] = new norm_operator();

    return true;
  }

  static bool predef_functions_initialized = init_predef_functions();


  //=========================================================================
  // Instructions for compilation: basic optimized operations on tensors
  //=========================================================================

  // #define GA_DEBUG_INSTRUCTION(a) { cout << a << endl; }
  #define GA_DEBUG_INSTRUCTION(a)

  struct ga_instruction {
    virtual void exec(void) = 0;
  };


  struct ga_instruction_slice_local_dofs : public ga_instruction {
    const mesh_fem &mf;
    const base_vector &U;
    fem_interpolation_context &ctx;
    base_vector &coeff;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: Slice");
      slice_vector_on_basic_dof_of_element(mf, U, ctx.convex_num(), coeff);
    }
    ga_instruction_slice_local_dofs(const mesh_fem &mf_, const base_vector &U_,
                                    fem_interpolation_context &ctx_,
                                    base_vector &coeff_)
      : mf(mf_), U(U_), ctx(ctx_), coeff(coeff_) {}
  };

  struct ga_instruction_update_pfp : public ga_instruction {
    const mesh_fem &mf;
    fem_interpolation_context &ctx;
    papprox_integration &pai;
    fem_precomp_pool &fp_pool;
    pfem_precomp &pfp;

    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: Pfp update");
      if (ctx.have_pgp()) {
        pfem pf = mf.fem_of_element(ctx.convex_num());
        if (pf != ctx.pf()) {
          if (pf->is_on_real_element()) 
            pfp = 0;
          else {
            
            pfp = fp_pool(pf, &(pai->integration_points()));
          }
        }
      } else {
        pfp = 0;
      }
    }

    ga_instruction_update_pfp(const mesh_fem &mf_, pfem_precomp &pfp_,
                              fem_interpolation_context &ctx_,
                              papprox_integration &pai_,
                              fem_precomp_pool &fp_pool_)
      : mf(mf_), ctx(ctx_), pai(pai_), fp_pool(fp_pool_), pfp(pfp_) {}
  };

  struct ga_instruction_first_ind_tensor : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    size_type qdim;
    const mesh_fem &mf;
    
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: Adapt first tensor index");
      pfem pf = mf.fem_of_element(ctx.convex_num());
      size_type Qmult = qdim / pf->target_dim();
      size_type s = pf->nb_dof(ctx.convex_num()) * Qmult;
      if (t.sizes()[0] != s)
        { bgeot::multi_index mi = t.sizes(); mi[0] = s; t.adjust_sizes(mi); }
    }

    ga_instruction_first_ind_tensor(base_tensor &t_,
                                    fem_interpolation_context &ctx_,
                                    size_type qdim_, const mesh_fem &mf_)
      : t(t_),  ctx(ctx_), qdim(qdim_), mf(mf_) {}
  };

  struct ga_instruction_two_first_ind_tensor : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    size_type qdim1;
    const mesh_fem &mf1;
    size_type qdim2;
    const mesh_fem &mf2;

    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: Adapt two first tensor indices");
      pfem pf1 = mf1.fem_of_element(ctx.convex_num());
      pfem pf2 = mf2.fem_of_element(ctx.convex_num());
      size_type Qmult1 = qdim1 / pf1->target_dim();
      size_type s1 = pf1->nb_dof(ctx.convex_num()) * Qmult1;
      size_type Qmult2 = qdim2 / pf2->target_dim();
      size_type s2 = pf2->nb_dof(ctx.convex_num()) * Qmult2;
      if (t.sizes()[0] != s1 || t.sizes()[1] != s2) {
        bgeot::multi_index mi = t.sizes();
        mi[0] = s1; mi[1] = s2;
        t.adjust_sizes(mi);
      }
    }

    ga_instruction_two_first_ind_tensor(base_tensor &t_,
                                        fem_interpolation_context &ctx_,
                                        size_type qdim1_, const mesh_fem &mf1_,
                                        size_type qdim2_, const mesh_fem &mf2_)
      : t(t_),  ctx(ctx_), qdim1(qdim1_), mf1(mf1_),
        qdim2(qdim2_), mf2(mf2_) {}
  };


  struct ga_instruction_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: compute value of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      ctx.pf()->real_base_value(ctx, t);
    }
    ga_instruction_base(base_tensor &tt, fem_interpolation_context &ct,
                        const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_grad_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: compute gradient of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      ctx.pf()->real_grad_base_value(ctx, t);
    }
    ga_instruction_grad_base(base_tensor &tt,
                             fem_interpolation_context &ct,
                             const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_hess_base : public ga_instruction {
    base_tensor &t;
    fem_interpolation_context &ctx;
    const mesh_fem &mf;
    pfem_precomp &pfp;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: compute Hessian of base functions");
      if (pfp) ctx.set_pfp(pfp);
      else ctx.set_pf(mf.fem_of_element(ctx.convex_num()));
      ctx.pf()->real_hess_base_value(ctx, t);
    }
    ga_instruction_hess_base(base_tensor &tt,
                             fem_interpolation_context &ct,
                             const mesh_fem &mf_, pfem_precomp &pfp_)
      : t(tt), ctx(ct), mf(mf_), pfp(pfp_) {}
  };

  struct ga_instruction_val : public ga_instruction {
    base_tensor &t;
    base_tensor &Z;
    const base_vector &coeff;
    size_type qdim;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: variable value");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type Qmult = qdim / target_dim;
      GMM_ASSERT1(t.size() == qdim, "dimensions mismatch"); // à supprimer ?
      GMM_ASSERT1(gmm::vect_size(coeff) == ndof*Qmult,      // à supprimer ?
                  "Wrong size for coeff vector");
      gmm::clear(t.as_vector());
      for (size_type j = 0; j < ndof; ++j) {
        for (size_type q = 0; q < Qmult; ++q) {
          scalar_type co = coeff[j*Qmult+q];
          for (size_type r = 0; r < target_dim; ++r) // à optimiser !
            t[r + q*target_dim] += co * Z[j + r*ndof];
        }
      }
    }
    ga_instruction_val(base_tensor &tt, base_tensor &Z_,
                       const base_vector &co, size_type q)
      : t(tt), Z(Z_), coeff(co), qdim(q) {}
  };


  struct ga_instruction_grad : public ga_instruction {
    base_tensor &t;
    base_tensor &Z;
    const base_vector &coeff;
    size_type qdim;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: gradient");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      GMM_ASSERT1(t.sizes()[1] == N && t.sizes()[0] == qdim, // à supprimer ?
                  "dimensions mismatch");
      gmm::clear(t.as_vector());
      for (size_type q = 0; q < Qmult; ++q) {
        base_tensor::const_iterator it = Z.begin();
        for (size_type k = 0; k < N; ++k)
          for (size_type r = 0; r < target_dim; ++r) // à optimiser ?
            for (size_type j = 0; j < ndof; ++j, ++it)
              t(r + q*target_dim, k) += coeff[j*Qmult+q] * (*it);
      }
    }

    ga_instruction_grad(base_tensor &tt, base_tensor &Z_,
                        const base_vector &co, size_type q)
      : t(tt), Z(Z_), coeff(co), qdim(q) {}
  };

  struct ga_instruction_hess : public ga_instruction {
    base_tensor &t;
    base_tensor &Z;
    const base_vector &coeff;
    size_type qdim;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: Hessian");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = size_type(round(sqrt(scalar_type(Z.sizes()[2]))));
      size_type Qmult = qdim / target_dim;
      GMM_ASSERT1(t.sizes()[1] == N && t.sizes()[2] == N // à supprimer ?
                  && t.sizes()[0] == qdim, "dimensions mismatch");
      gmm::clear(t.as_vector());
      for (size_type q = 0; q < Qmult; ++q) {
        base_tensor::const_iterator it = Z.begin();
        for (size_type k = 0; k < size_type(N); ++k)
          for (size_type l = 0; l < size_type(N); ++l)
            for (size_type r = 0; r < target_dim; ++r) // à optimiser ?
              for (size_type j = 0; j < ndof; ++j, ++it)
                t(r + q*target_dim, k, l) += coeff[j*Qmult+q] * (*it);
      }
    }

    ga_instruction_hess(base_tensor &tt, base_tensor &Z_,
                        const base_vector &co, size_type q)
      : t(tt), Z(Z_), coeff(co), qdim(q) {}
  };

  struct ga_instruction_copy_base : public ga_instruction {
    base_tensor &t;
    base_tensor &Z;
    size_type qdim;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: value of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type Qmult = qdim / target_dim;
      if (Qmult == 1) {
        gmm::copy(Z.as_vector(), t.as_vector());
      } else {
        gmm::clear(t.as_vector());
        base_tensor::const_iterator it = Z.begin();
        for (size_type i = 0; i < ndof; ++i, ++it) // TODO: à optimiser
          for (size_type j = 0; j < Qmult; ++j)
            t(i*Qmult+j, j) = *it;
      }
    }
    ga_instruction_copy_base(base_tensor &tt, base_tensor &Z_, size_type q)
      : t(tt), Z(Z_), qdim(q) {}
  };

  struct ga_instruction_copy_grad : public ga_instruction {
    base_tensor &t;
    base_tensor &Z;
    size_type qdim;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: gradient of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = Z.sizes()[2];
      size_type Qmult = qdim / target_dim;
      if (Qmult == 1) {
        gmm::copy(Z.as_vector(), t.as_vector());
      } else {
        gmm::clear(t.as_vector());
        base_tensor::const_iterator it = Z.begin();
        for (size_type k = 0; k < N; ++k)
          for (size_type i = 0; i < ndof; ++i, ++it)
            for (size_type j = 0; j < Qmult; ++j) // TODO: à optimiser
              t(i*Qmult+j, j, k) = *it;
      }
    }
    ga_instruction_copy_grad(base_tensor &tt, base_tensor &Z_, size_type q)
      : t(tt), Z(Z_), qdim(q) {}
  };

  struct ga_instruction_copy_hess : public ga_instruction {
    base_tensor &t;
    base_tensor &Z;
    size_type qdim;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: Hessian of test functions");
      size_type ndof = Z.sizes()[0];
      size_type target_dim = Z.sizes()[1];
      size_type N = size_type(round(sqrt(scalar_type(Z.sizes()[2]))));
      size_type Qmult = qdim / target_dim;
      if (Qmult == 1) {
        gmm::copy(Z.as_vector(), t.as_vector());
      } else {
        gmm::clear(t.as_vector());
        base_tensor::const_iterator it = Z.begin();
        for (size_type l = 0; l < N; ++l)
          for (size_type k = 0; k < N; ++k)
            for (size_type i = 0; i < ndof; ++i, ++it)
              for (size_type j = 0; j < Qmult; ++j) // TODO: à optimiser
                t(i*Qmult+j, j, k, l) = *it;
      }
    }
    ga_instruction_copy_hess(base_tensor &tt, base_tensor &Z_, size_type q)
      : t(tt), Z(Z_), qdim(q) {}
  };

  struct ga_instruction_add : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual void exec(void) { // TODO: use a daxpy
      GA_DEBUG_INSTRUCTION("Instruction: addition");
      GMM_ASSERT1(t.size() == tc1.size() && t.size() == tc2.size(),
                  "internal error");
      gmm::add(tc1.as_vector(), tc2.as_vector(), t.as_vector());
    }
    ga_instruction_add(base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_add_transp_test : public ga_instruction {
    base_tensor &t, &tc;
    virtual void exec(void) { // TODO: use a daxpy
      GA_DEBUG_INSTRUCTION("Instruction: add with a transpose in test "
                           "functions");
      GMM_ASSERT1(t.size() == tc.size(), "internal error");
      size_type s1 = tc.sizes()[0], s2 = tc.sizes()[1], s = s1*s2;

      for (base_tensor::iterator it = t.begin(), itc = tc.begin();
           it != t.end(); itc += s) {
        for (size_type i = 0; i < s1; ++i)
          for (size_type j = 0; j < s2; ++j)
            *it++ += itc[i+j*s1];
      }
    }
    ga_instruction_add_transp_test(base_tensor &t_, base_tensor &tc_)
      : t(t_), tc(tc_) {}
  };

  struct ga_instruction_sub : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    virtual void exec(void) { // TODO: use a daxpy
      GA_DEBUG_INSTRUCTION("Instruction: substraction");
      GMM_ASSERT1(t.size() == tc1.size() && t.size() == tc2.size(),
                  "internal error");
      gmm::add(tc1.as_vector(), gmm::scaled(tc2.as_vector(), scalar_type(-1)),
               t.as_vector());
    }
    ga_instruction_sub(base_tensor &t_, base_tensor &tc1_, base_tensor &tc2_)
      : t(t_), tc1(tc1_), tc2(tc2_) {}
  };

  struct ga_instruction_minus : public ga_instruction {
    base_tensor &t, &tc1;
    virtual void exec(void) { // TODO: use a daxpy
      GA_DEBUG_INSTRUCTION("Instruction: unary minus");
      GMM_ASSERT1(t.size() == tc1.size(), "internal error");
      gmm::copy(gmm::scaled(tc1.as_vector(), scalar_type(-1)), t.as_vector());
    }
    ga_instruction_minus(base_tensor &t_, base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_opposite : public ga_instruction {
    base_tensor &t;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: multiplication with -1");
      gmm::scale(t.as_vector(), scalar_type(-1));
    }
    ga_instruction_opposite(base_tensor &t_) : t(t_) {}
  };

  struct ga_instruction_print_tensor : public ga_instruction {
    base_tensor &t; pga_tree_node pnode;
    fem_interpolation_context &ctx;
    size_type &nbpt, &ipt;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: tensor print");
      cout << "Print term "; ga_print_node(pnode);
      cout << " on Gauss point " << ipt << "/" << nbpt << " of element "
           << ctx.convex_num() << ": " << t << endl;
    }
    ga_instruction_print_tensor(base_tensor &t_, pga_tree_node pnode_,
                                fem_interpolation_context &ctx_,
                                size_type &nbpt_, size_type &ipt_)
      : t(t_), pnode(pnode_), ctx(ctx_), nbpt(nbpt_), ipt(ipt_) {}
  };

  struct ga_instruction_copy_tensor : public ga_instruction {
    base_tensor &t, &tc1;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: tensor copy");
      gmm::copy(tc1.as_vector(), t.as_vector());
    }
    ga_instruction_copy_tensor(base_tensor &t_, base_tensor &tc1_)
      : t(t_), tc1(tc1_) {}
  };

  struct ga_instruction_reduction : public ga_instruction {
    base_tensor &t, &tc1, &tc2;
    size_type n;
    virtual void exec(void) { // TODO: use a dgemm
      GA_DEBUG_INSTRUCTION("Instruction: reduction operation of size " << n);
      // cout << "tc1 = " << tc1 << endl;
      // cout << "tc2 = " << tc2 << endl;
      size_type s1 = tc1.size()/n, s2 = tc2.size()/n, s = t.size();
      GMM_ASSERT1(s == s1*s2, "internal error");
      
      base_tensor::iterator it1=tc1.begin(), it2=tc2.begin(), it2end=it2 + s2;
      for (base_tensor::iterator it = t.begin(); it != t.end(); ++it) {
        scalar_type a(0);
        base_tensor::iterator it11 = it1, it22 = it2;
        for (size_type i = 0; i < n; ++i)
          { a+= (*it11) * (*it22); it11 += s1; it22 += s2; }
        *it = a;  
        ++it2; if (it2 == it2end) { it2 = tc2.begin(), ++it1; }
      }
      // cout << "result t = " << t << endl;
    }
    ga_instruction_reduction(base_tensor &t_, base_tensor &tc1_,
                             base_tensor &tc2_, size_type n_)
      : t(t_), tc1(tc1_), tc2(tc2_), n(n_) {}
  };

  struct ga_instruction_scalar_assembly : public ga_instruction {
    base_tensor &t;
    scalar_type &E, &coeff;
     virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: scalar term assembly ");
      E += t[0] * coeff;
     }
    ga_instruction_scalar_assembly(base_tensor &t_, scalar_type &E_,
                                   scalar_type &coeff_)
      : t(t_), E(E_), coeff(coeff_) {}
  };

  struct ga_instruction_vector_assembly : public ga_instruction {
    base_tensor &t;
    base_vector &V;
    fem_interpolation_context &ctx;
    gmm::sub_interval &I;
    const mesh_fem &mf;
    scalar_type &coeff;
     virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: vector term assembly ");
      mesh_fem::ind_dof_ct ct = mf.ind_basic_dof_of_element(ctx.convex_num());
      for (size_type i = 0; i < ct.size(); ++i)
        V[I.first()+ct[i]] += t[i] * coeff;
        
     }
    ga_instruction_vector_assembly(base_tensor &t_, base_vector &V_,
                                   fem_interpolation_context &ctx_,
                                   gmm::sub_interval &I_,
                                   const mesh_fem &mf_,
                                   scalar_type &coeff_)
      : t(t_), V(V_), ctx(ctx_), I(I_), mf(mf_), coeff(coeff_) {}
  };

  struct ga_instruction_matrix_assembly : public ga_instruction {
    base_tensor &t;
    model_real_sparse_matrix &K;
    fem_interpolation_context &ctx;
    gmm::sub_interval &I1, &I2;
    const mesh_fem &mf1, &mf2;
    scalar_type &coeff;
    size_type &nbpt, &ipt;
    base_vector &elem;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: matrix term assembly ");
      if (nbpt == 1) {
        mesh_fem::ind_dof_ct ct1
          = mf1.ind_basic_dof_of_element(ctx.convex_num());
        mesh_fem::ind_dof_ct ct2
          = mf2.ind_basic_dof_of_element(ctx.convex_num());
        size_type s1 = ct1.size(), s2 = ct2.size(); 
        for (size_type i1 = 0; i1 < s1; ++i1)
          for (size_type i2 = 0; i2 < s2; ++i2)
            K(I1.first()+ct1[i1], I2.first()+ct2[i2]) += t[i2*s1+i1]*coeff;
      } else {
        if (ipt == 0) {
          gmm::resize(elem, t.size());
          gmm::copy(gmm::scaled(t.as_vector(), coeff), elem);
        } else {
          gmm::add(gmm::scaled(t.as_vector(), coeff), elem);
        }
        if (ipt == nbpt-1) {
          mesh_fem::ind_dof_ct ct1
            = mf1.ind_basic_dof_of_element(ctx.convex_num());
          mesh_fem::ind_dof_ct ct2
            = mf2.ind_basic_dof_of_element(ctx.convex_num());
          scalar_type threshold = gmm::vect_norminf(elem) * 1E-13;
          size_type s1 = ct1.size(), s2 = ct2.size();
          for (size_type i1 = 0; i1 < s1; ++i1)
            for (size_type i2 = 0; i2 < s2; ++i2) {
              scalar_type e = elem[i2*s1+i1];
              if (gmm::abs(e) > threshold)
                K(I1.first()+ct1[i1], I2.first()+ct2[i2]) += e;
            }
        }
      }
    }
    ga_instruction_matrix_assembly(base_tensor &t_,
                                   model_real_sparse_matrix &K_,
                                   fem_interpolation_context &ctx_,
                                   gmm::sub_interval &I1_,
                                   gmm::sub_interval &I2_,
                                   const mesh_fem &mf1_,
                                   const mesh_fem &mf2_,
                                   scalar_type &coeff_, size_type &nbpt_,
                                   size_type &ipt_,  base_vector &elem_)
      : t(t_), K(K_), ctx(ctx_), I1(I1_), I2(I2_), mf1(mf1_), mf2(mf2_),
        coeff(coeff_), nbpt(nbpt_), ipt(ipt_), elem(elem_) {}
  };

  struct ga_instruction_transp_matrix_assembly : public ga_instruction {
    base_tensor &t;
    model_real_sparse_matrix &K;
    fem_interpolation_context &ctx;
    gmm::sub_interval &I1, &I2;
    const mesh_fem &mf1, &mf2;
    scalar_type &coeff;
    size_type &nbpt, &ipt;
    base_vector &elem;
    virtual void exec(void) {
      GA_DEBUG_INSTRUCTION("Instruction: (transposed) matrix term assembly ");
      if (nbpt == 1) {
        mesh_fem::ind_dof_ct ct1
          = mf1.ind_basic_dof_of_element(ctx.convex_num());
        mesh_fem::ind_dof_ct ct2
          = mf2.ind_basic_dof_of_element(ctx.convex_num());
        size_type s1 = ct1.size(), s2 = ct2.size(); 
        for (size_type i1 = 0; i1 < s1; ++i1)
          for (size_type i2 = 0; i2 < s2; ++i2)
            K(I1.first()+ct1[i1], I2.first()+ct2[i2]) += t[i1*s2+i2]*coeff;
      } else {
        if (ipt == 0) {
          gmm::resize(elem, t.size());
          gmm::copy(gmm::scaled(t.as_vector(), coeff), elem);
        } else {
          gmm::add(gmm::scaled(t.as_vector(), coeff), elem);
        }
        if (ipt == nbpt-1) {
          mesh_fem::ind_dof_ct ct1
            = mf1.ind_basic_dof_of_element(ctx.convex_num());
          mesh_fem::ind_dof_ct ct2
            = mf2.ind_basic_dof_of_element(ctx.convex_num());
          size_type s1 = ct1.size(), s2 = ct2.size();
          scalar_type threshold = gmm::vect_norminf(elem) * 1E-13;
          for (size_type i1 = 0; i1 < s1; ++i1)
            for (size_type i2 = 0; i2 < s2; ++i2) {
              scalar_type e = elem[i1*s2+i2];
              if (gmm::abs(e) > threshold)
                K(I1.first()+ct1[i1], I2.first()+ct2[i2]) += e;
            }
        }
      }
    }
    ga_instruction_transp_matrix_assembly(base_tensor &t_,
                                          model_real_sparse_matrix &K_,
                                          fem_interpolation_context &ctx_,
                                          gmm::sub_interval &I1_,
                                          gmm::sub_interval &I2_,
                                          const mesh_fem &mf1_,
                                          const mesh_fem &mf2_,
                                          scalar_type &coeff_,
                                          size_type &nbpt_,
                                          size_type &ipt_,  base_vector &elem_)
      :  t(t_), K(K_), ctx(ctx_), I1(I1_), I2(I2_), mf1(mf1_), mf2(mf2_),
         coeff(coeff_), nbpt(nbpt_), ipt(ipt_), elem(elem_)  {}
  };


  //=========================================================================
  // Structure to gather instructions
  //=========================================================================

  typedef ga_instruction *pga_instruction;
  typedef std::vector<pga_instruction> ga_instruction_list;
  
  struct ga_instruction_set {

    papprox_integration pai;       // Current approximation method
    fem_interpolation_context ctx; // Current fem interpolation context.
    scalar_type coeff;             // Coefficient for the Gauss point
    size_type nbpt, ipt;           // Number and index of Gauss point
    bgeot::geotrans_precomp_pool gp_pool;
    fem_precomp_pool fp_pool;

    typedef std::pair<const mesh_im *, size_type> region_mim;

    std::map<std::string, const base_vector *> extended_vars;
    std::map<std::string, base_vector> really_extended_vars;
    std::map<std::string, gmm::sub_interval> var_intervals;
    size_type nb_dof;
   

    struct region_mim_instructions {

      std::map<std::string, base_vector> local_dofs;
      std::map<const mesh_fem *, pfem_precomp> pfps;
      std::map<const mesh_fem *, base_tensor> base;
      std::map<const mesh_fem *, base_tensor> grad;
      std::map<const mesh_fem *, base_tensor> hess;
      ga_instruction_list instructions;

      ~region_mim_instructions(void) {
        for (size_type i = 0; i < instructions.size(); ++i)
          delete instructions[i];
      }
    };

    std::list<ga_tree> trees; // The trees are stored mainly because they
                              // contain the intermediary tensors.

    typedef std::map<region_mim,  region_mim_instructions> instructions_set;
    
    instructions_set  whole_instructions;

    ga_instruction_set(void) { nb_dof = 0; }
  };

  


  //=========================================================================
  // Structure dealing with user defined environment : constant, variables,
  // functions, operators.
  //=========================================================================

  class ga_instruction_set;

  static void ga_semantic_analysis(const std::string &, ga_tree &,
                                   const ga_workspace &, bool);
  static void ga_split_tree(const std::string &, ga_tree &,
                            ga_workspace &, pga_tree_node);
  static void ga_derivation(ga_tree &, const ga_workspace &,
                            const std::string &, size_type);
  static bool ga_node_mark_tree_for_variable(pga_tree_node,
                                             const std::string &);
  static void ga_exec(ga_instruction_set &gis);
  static void ga_compile(ga_workspace &workspace, ga_instruction_set &gis,
                         size_type order);

  static void ga_extract_variables(pga_tree_node pnode,
                                   std::set<std::string> &vars) {
    if (pnode->node_type == GA_NODE_VAL ||
        pnode->node_type == GA_NODE_GRAD ||
        pnode->node_type == GA_NODE_HESS)
      vars.insert(pnode->name);
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_extract_variables(pnode->children[i], vars);
  }

  ga_workspace::tree_description::tree_description(void)
  { ptree = 0; }

  void ga_workspace::add_tree(ga_tree &tree, const mesh_im &mim,
                              size_type region, const std::string expr) {
    if (tree.root) {
      
      cout << "adding tree: "; ga_print_tree(tree);
      cout << "With tests functions of " <<  tree.root->name_test1
           << " and " << tree.root->name_test2 << endl;
      bool remain = true;
      size_type order = 0, ind_tree = 0;
      switch(tree.root->test_function_type) {
      case 0: order = 0; break;
      case 1: order = 1; break;
      case 3: order = 2; break;
      case 4: order = 2; break;
      default: GMM_ASSERT1(false, "Inconsistent term");
      }
      
      if (tree.root->tensor_order() > 0)
        ga_throw_error(expr, tree.root->pos, "Incorrect term. Each term "
                       "should be reduced to a scalar in order to perform "
                       "the assembly.");

      bool found = false;
      for (size_type i = 0; i < trees.size(); ++i) {
        if (trees[i].mim == &mim && trees[i].order == order &&
            trees[i].name_test1.compare(tree.root->name_test1) == 0 &&
            trees[i].name_test2.compare(tree.root->name_test2) == 0) {
          ga_tree &ftree = *(trees[i].ptree);
            
          ftree.insert_node(ftree.root);
          ftree.root->node_type = GA_NODE_OP;
          ftree.root->op_type = GA_PLUS;
          ftree.add_children(ftree.root);
          ftree.copy_node(tree.root, ftree.root, ftree.root->children[1]);
          ga_semantic_analysis("", ftree, *this, false);
          found = true;
          break;
        }
      }

      if (!found) {
        ind_tree = trees.size(); remain = false;
        trees.push_back(tree_description());
        trees.back().mim = &mim; trees.back().region = region;
        trees.back().ptree = new ga_tree;
        trees.back().ptree->swap(tree);
        pga_tree_node root = trees.back().ptree->root;
        trees.back().name_test1 = root->name_test1;
        trees.back().name_test2 = root->name_test2;
        trees.back().order = order;
      }
      
      if (order <= 1) {
        std::set<std::string> expr_variables;
        ga_extract_variables((remain ? tree : *(trees[ind_tree].ptree)).root,
                             expr_variables);
        
        for (std::set<std::string>::iterator it = expr_variables.begin();
             it != expr_variables.end(); ++it) {
          if (!(is_constant(*it))) {
            ga_tree dtree = (remain ? tree : *(trees[ind_tree].ptree));
            cout << (order == 0 ? "First " : "Second ");
            cout << "derivative of tree: "; ga_print_tree(dtree);
            cout << "with respect to " << *it << endl;
            ga_derivation(dtree, *this, *it, 1+order);
            cout << "result: "; ga_print_tree(dtree);
            ga_semantic_analysis(expr, dtree, *this, false);
            add_tree(dtree, mim, region, expr);
          }
        }
      }
    }
  }

  void ga_workspace::clear_aux_trees(void) {
    for (std::list<ga_tree *>::iterator it = aux_trees.begin();
         it != aux_trees.end(); ++it)
      delete(*it);
    aux_trees.clear();
  }

  void ga_workspace::add_expression(const std::string expr,
                                    const mesh_im &mim, size_type region) {
    ga_tree tree;
    ga_read_string(expr, tree);
    ga_semantic_analysis(expr, tree, *this, false);
    if (tree.root) {
      ga_split_tree(expr, tree, *this, tree.root);
      
      for (std::list<ga_tree *>::iterator it = aux_trees.begin();
           it != aux_trees.end(); ++it)
        add_tree(*(*it), mim, region, expr);
      
      add_tree(tree, mim, region, expr);
      clear_aux_trees();
    }
  }

  void ga_workspace::clear_expressions(void) {
    clear_aux_trees();
    for (size_type i = 0; i < trees.size(); ++i)
      if (trees[i].ptree) delete trees[i].ptree;
    trees.clear();
  }
  
  void ga_workspace::add_aux_tree(ga_tree &tree) {
    ga_tree *a_tree = new ga_tree(tree); aux_trees.push_back(a_tree);
  }
  
  size_type ga_workspace::nb_trees(void) { return trees.size(); }
  
  ga_workspace::tree_description &ga_workspace::tree_info(size_type i)
  { return trees[i]; }
  
  // TODO: methods to add a function or an operator
  
  ga_workspace::ga_workspace(const getfem::model &md) : model(&md) {}
  ga_workspace::ga_workspace(void) : model(0) {}
  ga_workspace::~ga_workspace() { clear_expressions(); }

  void ga_workspace::assembly(size_type order) {
    ga_instruction_set gis;
    ga_compile(*this, gis, order);
    size_type ndof = gis.nb_dof;
    gmm::resize(K, ndof, ndof); gmm::clear(K);
    gmm::resize(V, ndof); gmm::clear(V);
    E = 0;
    ga_exec(gis);
  }



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
    return (e != GA_NODE_ZERO) ? cos(M_E + scalar_type(e))
      : ga_hash_code(GA_NODE_CONSTANT);
  }

  static scalar_type ga_hash_code(pga_tree_node pnode) {
    scalar_type c = ga_hash_code(pnode->node_type);
    switch (pnode->node_type) {
    case GA_NODE_CONSTANT: case GA_NODE_ZERO:
      c += ga_hash_code(pnode->t); break;
      
    case GA_NODE_VAL: case GA_NODE_GRAD: case GA_NODE_HESS:
    case GA_NODE_TEST:
    case GA_NODE_GRAD_TEST: case GA_NODE_HESS_TEST:
      c += ga_hash_code(pnode->name); break;

    case GA_NODE_PREDEF_FUNC:
    case GA_NODE_SPEC_FUNC:
    case GA_NODE_OPERATOR:
      c += ga_hash_code(pnode->name)
        + tanh(scalar_type(pnode->der1)/M_PI + scalar_type(pnode->der2)*M_PI);
      break;
    default: break;
    }
    return c;
  }
  

  //=========================================================================
  // Semantic analysis, tree simplification and tree enrichment
  //    - Control tensor sizes for operations, operator or function call
  //    - Compute all constant operations (i.e. non element dependent)
  //    - Build a ready to use tree for derivation/compilation
  //=========================================================================

  static void ga_node_analysis(const std::string &expr, ga_tree &tree,
                               const ga_workspace &workspace,
                               pga_tree_node pnode, bool eval_fixed_size) {
    
    bool all_cte = true, all_sc = true; // all_primal = true;
    for (size_type i = 0; i < pnode->children.size(); ++i) {
      ga_node_analysis(expr, tree, workspace, pnode->children[i],
                       eval_fixed_size);
      all_cte = all_cte && (pnode->children[i]->node_type == GA_NODE_CONSTANT);
      all_sc = all_sc && pnode->children[i]->t.size() == 1;
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

    if (pnode->node_type != GA_NODE_OP ||
        (pnode->op_type != GA_PLUS && pnode->op_type != GA_MINUS))
      pnode->marked = false;
    else {
      if (pnode->parent) pnode->marked = pnode->parent->marked;
      else pnode->marked = true;
    }

    switch (pnode->node_type) {
    case GA_NODE_PREDEF_FUNC: case GA_NODE_OPERATOR: case GA_NODE_SPEC_FUNC :
    case GA_NODE_CONSTANT:
      pnode->test_function_type = 0; break;

    case GA_NODE_ALLINDICES: pnode->test_function_type = 0; break;
    case GA_NODE_VAL:
      if (eval_fixed_size && !(workspace.associated_mf(pnode->name))) {
        gmm::copy(workspace.value(pnode->name), pnode->t.as_vector());
        pnode->node_type = GA_NODE_CONSTANT;
      }
      break;

    case GA_NODE_ZERO: case GA_NODE_GRAD: case GA_NODE_HESS: break;

    case GA_NODE_TEST: case GA_NODE_GRAD_TEST: case GA_NODE_HESS_TEST: {
      const mesh_fem *mf = workspace.associated_mf(pnode->name);
      if (pnode->test_function_type == 1) {
        pnode->name_test1 = pnode->name;
        pnode->name_test2 = "";
        pnode->qdim1 = (mf ? workspace.qdim(pnode->name)
                        : gmm::vect_size(workspace.value(pnode->name)));
      } else {
        pnode->name_test1 = "";
        pnode->name_test2 = pnode->name;
        pnode->qdim2 = (mf ? workspace.qdim(pnode->name)
                        : gmm::vect_size(workspace.value(pnode->name)));
      }
      if (!mf) {
        size_type n = workspace.qdim(pnode->name);
        if (n == 1) {
          pnode->init_vector_tensor(1);
          pnode->t[0] = scalar_type(1);
        } else {
          pnode->init_matrix_tensor(n,n);
          for (size_type i = 0; i < n; ++i)
            for (size_type j = 0; j < n; ++j)
              pnode->t(i,j) = (i == j) ? scalar_type(1) : scalar_type(0);
        }
      }
    }  break;

    case GA_NODE_OP:
      switch(pnode->op_type) {

      case GA_PLUS: case GA_MINUS:
        {
          size_type c_size = std::min(size0.size(), size1.size());
          bool compatible = true;
          for (size_type i = 0; i < c_size; ++i)
            if (size0[i] != size1[i]) compatible = false; 
          for (size_type i = c_size; i < size0.size(); ++i)
            if (size0[i] != 1) compatible = false;
          for (size_type i = c_size; i < size1.size(); ++i)
            if (size1[i] != 1) compatible = false;

          if (child0->test_function_type || child1->test_function_type) {
            if (std::max(child0->test_function_type, size_type(3))
                != std::max(child1->test_function_type, size_type(3)) ||
                (!(pnode->marked) && 
                 (child0->name_test1.compare(child1->name_test1) ||
                  child0->name_test2.compare(child1->name_test2))))
              compatible = false;
          }
          if (!compatible)
            ga_throw_error(expr, pnode->pos, "Addition or substraction of "
                            "incompatible expressions or of different sizes");
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
            pnode->test_function_type
              = std::min(child0->test_function_type,
                         child1->test_function_type);
            pnode->name_test1 = child0->name_test1;
            pnode->name_test2 = child0->name_test2;
            pnode->qdim1 = child0->qdim1;
            pnode->qdim2 = child0->qdim2;

            // simplification if one of the two operands is constant and zero
            if (child0->tensor_is_zero())
              tree.replace_node_by_child(pnode, 1);
            else if (child1->tensor_is_zero())
              tree.replace_node_by_child(pnode, 0);
          }
        }
        break;

      case GA_DOTMULT: case GA_DOTDIV:
        {
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
            ga_throw_error(expr, pnode->pos, "Arguments of different sizes"
                           "for .* or ./");
          
          pnode->mult_test(child0, child1, expr, workspace);
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
        pnode->qdim1 = child0->qdim1;
        pnode->qdim2 = child0->qdim2;
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->test_function_type = 0;
          gmm::scale(pnode->t.as_vector(), scalar_type(-1));
          tree.clear_children(pnode);
        } else if (child0->node_type == GA_NODE_ZERO) {
          tree.replace_node_by_child(pnode, 0);
        }
        break;

      case GA_QUOTE:
        if (dim0 > 2)
          ga_throw_error(expr, pnode->pos, "Transpose operator is for "
                         "vectors or matrices only.");
        mi = size0;
        if (child0->tensor_proper_size() == 1)
          { tree.replace_node_by_child(pnode, 0); break; }
        else if (dim0 == 2) std::swap(mi.back(), mi[size0.size()-2]);
        else { size_type N = mi.back(); mi.back() = 1; mi.push_back(N); }

        pnode->t.adjust_sizes(mi);
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
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

      case GA_PRINT:
        {
          pnode->t = child0->t; 
          pnode->test_function_type = child0->test_function_type;
          pnode->name_test1 = child0->name_test1;
          pnode->name_test2 = child0->name_test2;
          pnode->qdim1 = child0->qdim1;
          pnode->qdim2 = child0->qdim2;
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            cout << "Print constant term "; ga_print_node(child0);
            cout << ": " << pnode->t << endl;
            ga_throw_error_msg(expr, child0->pos,
                               "Generic assembly Print command");
            tree.clear_children(pnode);
          } else if (child0->node_type == GA_NODE_ZERO) {
            pnode->node_type = GA_NODE_ZERO;
            gmm::clear(pnode->t.as_vector());
            cout << "Print zero term "; ga_print_node(child0);
            cout << ": " << pnode->t << endl;
            tree.clear_children(pnode);
          }
        }
        break;

      case GA_DOT:
        if (dim1 > 1)
          ga_throw_error(expr, pnode->pos, "The second argument of the dot "
                         "product have to be a vector.")
        else {
          size_type s0 = dim0 == 0 ? 1 : size0.back();
          size_type s1 = dim1 == 0 ? 1 : size1.back();
          if (s0 != s1) ga_throw_error(expr, pnode->pos, "Dot product "
                                       "of expressions of different sizes");
          pnode->mult_test(child0, child1, expr, workspace);
          if (pnode->test_function_type == 4 && dim0 <= 1 &&
              child1->test_function_type == 2 &&
              child0->test_function_type == 1)
            pnode->test_function_type = 3;
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
                            "of expressions of different sizes");
          pnode->mult_test(child0, child1, expr, workspace);
          if (pnode->test_function_type == 4 && dim0 <= 2 &&
              child1->test_function_type == 2 &&
              child0->test_function_type == 1)
            pnode->test_function_type = 3;
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
            if (dim0+dim1 > 4)
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
          pnode->mult_test(child0, child1, expr, workspace);
          if (pnode->test_function_type == 4 && dim0 == 0 &&
              child1->test_function_type == 2 &&
              child0->test_function_type == 1)
            pnode->test_function_type = 3;
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
              if (dim0+dim1 > 4)
                ga_throw_error(expr, pnode->pos, "Unauthorized "
                                "tensor multiplication.");
              for (size_type i = 0; i < dim0; ++i)
                mi.push_back(child0->t.size(i));
              for (size_type i = 0; i < dim1; ++i)
                mi.push_back(child1->t.size(i));
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
              ga_throw_error(expr, pnode->pos, "Incompatible sizes in "
                              "matrix-vector multiplication.");
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
              ga_throw_error(expr, pnode->pos, "Incompatible sizes in "
                              "matrix-matrix multiplication.");
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
              ga_throw_error(expr, pnode->pos, "Incompatible sizes in "
                              "tensor-matrix multiplication.");
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
          pnode->mult_test(child0, child1, expr, workspace);
          mi = pnode->t.sizes();

          if (child0->tensor_proper_size() == 1 &&
              child1->tensor_proper_size() == 1) {
          } else if (child0->tensor_proper_size() == 1) {
            for (size_type i = 0; i < dim1; ++i)
                mi.push_back(child1->tensor_proper_size(i));
          } else if (child1->tensor_proper_size() == 1) {
            for (size_type i = 0; i < dim0; ++i)
                mi.push_back(child0->tensor_proper_size(i));
          } else if (child0->tensor_order() == 2 &&
                     child1->tensor_order() == 1) {
            size_type m = child0->tensor_proper_size(0);
            size_type n = child0->tensor_proper_size(1);
            mi.push_back(m);
            if (n != child1->tensor_proper_size(0))
              ga_throw_error(expr, pnode->pos, "Incompatible sizes in "
                              "matrix-vector multiplication.");
          } else if (child0->tensor_order() == 2 &&
                     child1->tensor_order() == 2) {
            size_type m = child0->tensor_proper_size(0);
            size_type n = child0->tensor_proper_size(1);
            size_type p = child1->tensor_proper_size(1);
            mi.push_back(m); mi.push_back(p);
            if (n != child1->tensor_proper_size(0))
              ga_throw_error(expr, pnode->pos, "Incompatible sizes in "
                              "matrix-matrix multiplication.");
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
              ga_throw_error(expr, pnode->pos, "Incompatible sizes in "
                              "tensor-matrix multiplication.");
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
          } else if (child1->node_type == GA_NODE_CONSTANT &&
                     child1->t.size() == 1 && child1->t[0] == scalar_type(1)) {
            tree.replace_node_by_child(pnode, 0);
          }
        }
        break;

      case GA_DIV:
        if (child1->tensor_order() > 0)
          ga_throw_error(expr, pnode->pos, "Only the division by a scalar "
                          "is allowed.");
        if (child1->test_function_type)
          ga_throw_error(expr, pnode->pos, "Division by test functions "
                          "is not allowed.");
        if (child1->node_type == GA_NODE_CONSTANT &&
            child1->t[0] == scalar_type(0))
          ga_throw_error(expr, pnode->children[1]->pos, "Division by zero");
        
        pnode->t = child0->t;
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
        pnode->qdim1 = child0->qdim1;
        pnode->qdim2 = child0->qdim2;
        
        if (all_cte) {          
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->t = pnode->children[0]->t;
          pnode->test_function_type = 0;
          gmm::scale(pnode->t.as_vector(),
                     scalar_type(1) / pnode->children[1]->t[0]);
          tree.clear_children(pnode);
        }
        if (child0->tensor_is_zero()) {
          gmm::clear(pnode->t.as_vector());
          pnode->node_type = GA_NODE_ZERO;
          tree.clear_children(pnode);
        } else if (child1->node_type == GA_NODE_CONSTANT &&
                   child1->t.size() == 1 && child1->t[0] == scalar_type(1)) {
          tree.replace_node_by_child(pnode, 0);
        }
        break;
        
      default:GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX:
      {
        if (!all_sc) 
          ga_throw_error(expr, pnode->pos, "Constant vector/matrix/tensor "
                          "components should be scalar valued.");
//         if (!all_primal)
//           ga_throw_error(expr, pnode->pos, "Test functions are not allowed "
//                           "in constant vector/matrix/tensor components.");
        
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
              pnode->qdim1 = pnode->children[i]->qdim1;
              pnode->qdim2 = pnode->children[i]->qdim2;
            } else {
              if (std::max(pnode->test_function_type, size_type(3))
                  != std::max(pnode->children[i]->test_function_type,
                              size_type(3)) ||
                  pnode->name_test1.compare(pnode->children[i]->name_test1)
                  != 0 ||
                  pnode->name_test2.compare(pnode->children[i]->name_test2)
                  != 0)
                ga_throw_error(expr, pnode->pos, "Inconsistent use of test "
                               "function in constant matrix.");
            }
          }
        }
        mi.resize(0);
        if (pnode->test_function_type) mi.push_back(1);
        if (pnode->test_function_type >= 3) mi.push_back(1);
        if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1 && nbl == 1) {
          pnode->t.adjust_sizes(mi);
          if (all_cte) pnode->t[0] = child0->t[0];
        } else if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
          mi.push_back(nbl);
          pnode->t.adjust_sizes(mi);
          if (all_cte)
            for (size_type i = 0; i < nbl; ++i)
              pnode->t[i] = pnode->children[i]->t[0];
        } else if (nbc2 == 1 && nbc3 == 1) {
          mi.push_back(nbl); mi.push_back(nbc1); 
          pnode->t.adjust_sizes(mi);
          if (all_cte)
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc1; ++j)
                pnode->t(i,j) = pnode->children[i*nbc1+j]->t[0];
        } else {
          mi.push_back(nbl); mi.push_back(nbc3);
          mi.push_back(nbc2); mi.push_back(nbc1);
          pnode->t.adjust_sizes(mi);
          size_type n = 0;
          if (all_cte)
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc3; ++j)
                for (size_type k = 0; k < nbc2; ++k)
                  for (size_type l = 0; l < nbc1; ++l)
                    pnode->t(i,j,k,l) = pnode->children[n++]->t[0];
        }
        if (all_cte) tree.clear_children(pnode);
      }
      break;


    case GA_NODE_NAME:
      {
        std::string name = pnode->name;

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

        ga_predef_function_tab::const_iterator it=PREDEF_FUNCTIONS.find(name);
        if (it != PREDEF_FUNCTIONS.end() ||
            workspace.user_function_exists(name)) {
          // Predefined function found
          pnode->node_type = GA_NODE_PREDEF_FUNC;
          pnode->name = name;
          pnode->test_function_type = 0;
          if (pnode->der1) {
            if (pnode->der1 > it->second.nbargs
                || pnode->der2 > it->second.nbargs)
              ga_throw_error(expr, pnode->pos, "Special functions do not "
                             "support derivatives.");
            const ga_predef_function &F = (it != PREDEF_FUNCTIONS.end()) ?
              it->second : workspace.user_function(name);
            if (F.ftype == 0 && !(pnode->der2)) {
              pnode->name = ((pnode->der1 == 1) ?
                             F.derivative1 : F.derivative2);
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
          }
        } else if (PREDEF_OPERATORS.find(name) != PREDEF_OPERATORS.end()
                   || workspace.user_operator_exists(name)) {
          // Nonlinear operator found
          pnode->node_type = GA_NODE_OPERATOR;
          pnode->name = name;
          pnode->test_function_type = 0;
        } else {
          // Search for a variable name with optional gradient, Hessian 
          // or test functions

          int val_grad_or_hess = 0;
          if (name.size() >= 5 && name.compare(0, 5, "Grad_") == 0)
            { val_grad_or_hess = 1; name = name.substr(5); }
          else if (name.size() >= 5 && name.compare(0, 5, "Hess_") == 0)
            { val_grad_or_hess = 2; name = name.substr(5); }
          int test = 0;
          if (name.size() >= 6 && name.compare(0, 6, "Test2_") == 0)
            { test = 2; name = name.substr(6); }
          else if (name.size() >= 5 && name.compare(0, 5, "Test_") == 0)
            { test = 1; name = name.substr(5); }
          
          if (!(workspace.variable_exists(name)))
            ga_throw_error(expr, pnode->pos,
                           "Unknown variable, function, operator or data");

          if (pnode->der1)
            ga_throw_error(expr, pnode->pos, "Derivative is for functions or "
                           "operators, not for variables. Use Grad instead.");
          pnode->name = name;
          
          const mesh_fem *mf = workspace.associated_mf(name);
          if (!test) {
            if (!mf) {
              if (val_grad_or_hess)
                ga_throw_error(expr, pnode->pos, "Gradient or Hessian cannot "
                                "be evaluated for fixed size data.");
              if (eval_fixed_size)
                pnode->node_type = GA_NODE_CONSTANT;
              else
                pnode->node_type = GA_NODE_VAL;
              size_type n = gmm::vect_size(workspace.value(name));
              if (n == 1) {
                pnode->init_scalar_tensor(workspace.value(name)[0]);
              } else {
                pnode->init_vector_tensor(n);
                gmm::copy(workspace.value(name), pnode->t.as_vector());
              }
            } else {
              size_type q = workspace.qdim(name), n = mf->linked_mesh().dim();
              switch (val_grad_or_hess) {
              case 0: // value
                pnode->node_type = GA_NODE_VAL;
                if (q == 1)
                  pnode->init_scalar_tensor(scalar_type(0));
                else
                  pnode->init_vector_tensor(q);
                break;
              case 1: // grad
                pnode->node_type = GA_NODE_GRAD;
                if (q == 1 && n == 1)
                  pnode->init_scalar_tensor(scalar_type(0));
                else if (q == 1)
                  pnode->init_vector_tensor(n);
                else
                  pnode->init_matrix_tensor(q, n);
                break;
              case 2: // Hessian
                pnode->node_type = GA_NODE_HESS;
                if (q == 1 && n == 1)
                  pnode->init_scalar_tensor(scalar_type(0));
                else if (q == 1)
                  pnode->init_matrix_tensor(n,n);
                else
                  pnode->init_third_order_tensor(q, n, n);
                break;
              }
            }
          } else {
            if (workspace.is_constant(name))
              ga_throw_error(expr, pnode->pos, "Test functions of constant "
                              "are not allowed.");
            if (test == 1) {
              pnode->name_test1 = name;
              pnode->qdim1
                = (mf ? workspace.qdim(name)
                   : gmm::vect_size(workspace.value(name)));
            } else {
              pnode->name_test2 = name;
              pnode->qdim2
                = (mf ? workspace.qdim(name)
                   : gmm::vect_size(workspace.value(name)));
            }
          
            if (!mf) {
              if (val_grad_or_hess)
                ga_throw_error(expr, pnode->pos, "Gradient or Hessian cannot "
                                "be evaluated for fixed size variables.");
              pnode->node_type = GA_NODE_TEST;
              size_type n = gmm::vect_size(workspace.value(name));
              if (n == 1) {
                pnode->init_vector_tensor(1);
                pnode->t[0] = scalar_type(1);
                pnode->test_function_type = test;
              } else {
                pnode->init_matrix_tensor(n,n);
                pnode->test_function_type = test;
                for (size_type i = 0; i < n; ++i)
                  for (size_type j = 0; j < n; ++j)
                    pnode->t(i,j) = (i == j) ? scalar_type(1) : scalar_type(0);
              }
            } else {
              size_type q = workspace.qdim(name), n = mf->linked_mesh().dim();
              switch (val_grad_or_hess) {
              case 0: // value
                pnode->node_type = GA_NODE_TEST;
                if (q == 1)
                  pnode->init_vector_tensor(1);
                else
                  pnode->init_matrix_tensor(1,q);
                pnode->test_function_type = test;
                break;
              case 1: // grad
                pnode->node_type = GA_NODE_GRAD_TEST;
                if (q == 1 && n == 1)
                  pnode->init_vector_tensor(1);
                else if (q == 1)
                  pnode->init_matrix_tensor(1,n);
                else
                  pnode->init_third_order_tensor(1,q,n);
                pnode->test_function_type = test;
                break;
              case 2: // hessian
                pnode->node_type = GA_NODE_HESS_TEST;
                if (q == 1 && n == 1)
                  pnode->init_vector_tensor(1);
                else if (q == 1)
                  pnode->init_third_order_tensor(1,n,n);
                else
                  pnode->init_fourth_order_tensor(1,q,n,n);
                pnode->test_function_type = test;
                break;
              }
            }
          }
        }
      }
      break;

    case GA_NODE_PARAMS:

      if (child0->node_type == GA_NODE_PREDEF_FUNC) {
        
        // Évaluation of predefined function
        
        for (size_type i = 1; i < pnode->children.size(); ++i)
          ga_valid_operand(expr, pnode->children[i]);
        std::string name = child0->name;
        ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
        const ga_predef_function &F = (it != PREDEF_FUNCTIONS.end()) ?
          it->second : workspace.user_function(name);
        size_type nbargs = F.nbargs;
        GMM_ASSERT1(F.ftype != 2, "Not already taken into account");
        if (nbargs+1 != pnode->children.size()) {
            std::stringstream msg;
            msg << "Bad number of arguments for predefined function "
                << name << ". Found " << pnode->children.size()-1
                << " should be " << nbargs << ".";
            ga_throw_error(expr, pnode->pos, msg.str());
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
              pnode->t[i] = (*(F.f1))(child1->t[i]);
          } else {
            if (s1 == s2) {
              for (size_type i = 0; i < s1; ++i)
                pnode->t[i] = (*(F.f2))(child1->t[i], child2->t[i]);
            } else if (s1 == 1) {
              for (size_type i = 0; i < s2; ++i)
                pnode->t[i] = (*(F.f2))(child1->t[0], child2->t[i]);
            } else {
              for (size_type i = 0; i < s1; ++i)
                pnode->t[i] = (*(F.f2))(child1->t[i], child2->t[0]);
            }
          }
          tree.clear_children(pnode);
        }
      } else if (child0->node_type == GA_NODE_SPEC_FUNC) {

        // Special constant functions: meshdim(u), qdim(u) ...
        
        for (size_type i = 1; i < pnode->children.size(); ++i)
          ga_valid_operand(expr, pnode->children[i]);
        if (pnode->children.size() != 2)
          ga_throw_error(expr, pnode->pos,
                         "One and only one argument is allowed for function "
                         +child0->name+".");

        if (!(child0->name.compare("meshdim"))) {
          bool valid = (child1->node_type == GA_NODE_VAL);
          const mesh_fem *mf = valid ? workspace.associated_mf(child1->name):0;
          if (!mf)
            ga_throw_error(expr, pnode->pos, "The argument of mesh_dim "
                           "function can only be a fem variable name.");
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->init_scalar_tensor(scalar_type(mf->linked_mesh().dim()));
        } else if (!(child0->name.compare("qdim"))) {
          if (child1->node_type != GA_NODE_VAL)
            ga_throw_error(expr, pnode->pos, "The argument of qdim "
                           "function can only be a variable name.");
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->init_scalar_tensor(scalar_type(workspace.qdim(child1->name)));
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
        ga_predef_operator_tab::iterator it
          = PREDEF_OPERATORS.find(child0->name);
        const ga_nonlinear_operator &OP = (it != PREDEF_OPERATORS.end()) ?
          *(it->second) : workspace.user_operator(child0->name);
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
                              "Index out of range.");
          }
        }
        mi = size0; mi.resize(child0->nb_test_functions());
        for (size_type i = 0; i < mi2.size(); ++i) mi.push_back(mi2[i]);
        pnode->t.adjust_sizes(mi);
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test1 = child0->name_test1;
        pnode->name_test2 = child0->name_test2;
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
    pnode->hash_value = ga_hash_code(pnode);
    for (size_type i = 0; i < pnode->children.size(); ++i) {
      pnode->hash_value += (pnode->children[i]->hash_value)
        * M_PI * scalar_type(i+1); 
    }
  }

  static void ga_semantic_analysis(const std::string &expr, ga_tree &tree,
                                   const ga_workspace &workspace,
                                   bool eval_fixed_size) {
    GMM_ASSERT1(predef_functions_initialized, "Internal error");
    if (!(tree.root)) return;
    ga_node_analysis(expr, tree, workspace, tree.root, eval_fixed_size);
    ga_valid_operand(expr, tree.root);
  }

  //=========================================================================
  // Splitting of the terms which depend on different test functions.
  // To be performed just after a semantic analysis.
  //=========================================================================

  static void ga_split_tree(const std::string &expr, ga_tree &tree,
                            ga_workspace &workspace,
                            pga_tree_node pnode) {
    
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_split_tree(expr, tree, workspace, pnode->children[i]);
    
    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;

    switch (pnode->node_type) {

    case GA_NODE_OP:
      switch(pnode->op_type) {

      case GA_PLUS: case GA_MINUS:
        {
          bool compatible = true;
          if (child0->test_function_type || child1->test_function_type) {
            if (std::max(child0->test_function_type, size_type(3))
                != std::max(child1->test_function_type, size_type(3)) ||
                (child0->name_test1.compare(child1->name_test1) ||
                 child0->name_test2.compare(child1->name_test2)))
              compatible = false;
          }
          
          if (!compatible) {
            ga_tree aux_tree;
            aux_tree.root = child1;
            child1->parent = 0;
            
            pnode->children.pop_back();
            tree.replace_node_by_child(pnode, 0);
            
            workspace.add_aux_tree(aux_tree);
          }
        }
        break;
      default: break;
      }
      break;
    default: break;
    }
  }


  //=========================================================================
  // Derivation algorithm: derivation of a tree with respect to a variable
  //   The result tree is not ready to use. It has to be passed again in
  //   ga_semantic_analysis for enrichment.
  //=========================================================================

  static bool ga_node_mark_tree_for_variable(pga_tree_node pnode,
                                             const std::string &varname) {
    bool marked = false;
    for (size_type i = 0; i < pnode->children.size(); ++i)
      if (ga_node_mark_tree_for_variable(pnode->children[i], varname))
        marked = true;
    
    if ((pnode->node_type == GA_NODE_VAL ||
        pnode->node_type == GA_NODE_GRAD ||
        pnode->node_type == GA_NODE_HESS) &&
        (pnode->name.compare(varname) == 0)) marked = true;

    pnode->marked = marked;
    return marked;
  }
  

  static void ga_node_derivation(ga_tree &tree, const ga_workspace &workspace,
                                 pga_tree_node pnode,
                                 const std::string &varname, size_type order) {

    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bool mark0 = ((nbch > 0) ? child0->marked : false);
    bool mark1 = ((nbch > 1) ? child1->marked : false);
    bgeot::multi_index mi;

    switch (pnode->node_type) {
    case GA_NODE_VAL:
      mi.resize(1); mi[0] = 1;
      for (size_type i = 0; i < pnode->tensor_order(); ++i)
        mi.push_back(pnode->tensor_proper_size(i));
      pnode->t.adjust_sizes(mi);
      pnode->node_type = GA_NODE_TEST;
      pnode->test_function_type = order;
      break;
    case GA_NODE_GRAD:
      mi.resize(1); mi[0] = 1;
      for (size_type i = 0; i < pnode->tensor_order(); ++i)
        mi.push_back(pnode->tensor_proper_size(i));
      pnode->t.adjust_sizes(mi);   
      pnode->node_type = GA_NODE_GRAD_TEST;
      pnode->test_function_type = order;
      break;
    case GA_NODE_HESS:
      mi.resize(1); mi[0] = 1;
      for (size_type i = 0; i < pnode->tensor_order(); ++i)
        mi.push_back(pnode->tensor_proper_size(i));
      pnode->t.adjust_sizes(mi);   
      pnode->node_type = GA_NODE_HESS_TEST;
      pnode->test_function_type = order;
      break;
    case GA_NODE_OP:
      switch(pnode->op_type) {
        case GA_PLUS: case GA_MINUS:
          if (mark0 && mark1) {
            ga_node_derivation(tree, workspace, child0, varname, order);
            ga_node_derivation(tree, workspace, child1, varname, order);
          } else if (mark0) {
            ga_node_derivation(tree, workspace, child0, varname, order);
            tree.replace_node_by_child(pnode, 0);
          } else {
            ga_node_derivation(tree, workspace, child1, varname, order);
            if (pnode->op_type == GA_MINUS) {
              pnode->op_type = GA_UNARY_MINUS;
              tree.clear_node(child0);
            }
            else
              tree.replace_node_by_child(pnode, 1);
          }
          break;

      case GA_UNARY_MINUS: case GA_QUOTE: case GA_TRACE: case GA_PRINT:
        ga_node_derivation(tree, workspace, child0, varname, order);
        break;

      case GA_DOT: case GA_MULT: case GA_COLON: case GA_TMULT:
      case GA_DOTMULT:
        if (mark0 && mark1) {
          if (sub_tree_are_equal(child0, child1) &&
              (pnode->op_type != GA_MULT || child0->tensor_order() < 2)) {
            ga_node_derivation(tree, workspace, child1, varname, order);
            tree.insert_node(pnode);
            pnode->parent->node_type = GA_NODE_OP;
            pnode->parent->op_type = GA_MULT;
            tree.add_children(pnode->parent);
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
            ga_node_derivation(tree, workspace, child0, varname, order);
            ga_node_derivation(tree, workspace,
                               pnode->parent->children[1]->children[1],
                               varname, order);
          }
        } else if (mark0) {
          ga_node_derivation(tree, workspace, child0, varname, order);
        } else
          ga_node_derivation(tree, workspace, child1, varname, order);
        break;
          
      case GA_DIV: case GA_DOTDIV:
        if (mark1) {
          if (mark0) {
            tree.duplicate_with_addition(pnode);
            ga_node_derivation(tree, workspace, child0, varname, order);
            pnode->parent->op_type = GA_MINUS;
            pnode = pnode->parent->children[1];
          }
          tree.insert_node(pnode->children[1]);
          pga_tree_node pnode_param = pnode->children[1];
          pnode_param->node_type = GA_NODE_PARAMS;
          tree.add_children(pnode_param);
          std::swap(pnode_param->children[0], pnode_param->children[1]);
          pnode_param->children[0]->node_type = GA_NODE_PREDEF_FUNC;
          pnode_param->children[0]->name = "sqr";
          tree.insert_node(pnode);
          pga_tree_node pnode_mult = pnode->parent;
          pnode_mult->node_type = GA_NODE_OP;
          if (pnode->op_type == GA_DOTDIV)
            pnode_mult->op_type = GA_DOTMULT;
          else
            pnode_mult->op_type = GA_MULT;
          pnode_mult->children.push_back(0);
          tree.copy_node(pnode_param->children[1],
                         pnode_mult, pnode_mult->children[1]);
          ga_node_derivation(tree, workspace,pnode_mult->children[1],
                             varname, order);
        } else {
          ga_node_derivation(tree, workspace, child0, varname, order);
        }
        break;

      default: GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX:
      for (size_type i = 0; i < pnode->children.size(); ++i) {
        if (pnode->children[i]->marked)
          ga_node_derivation(tree, workspace, pnode->children[i],
                             varname, order);
        else {
          pnode->children[i]->init_scalar_tensor(scalar_type(0));
          pnode->children[i]->node_type = GA_NODE_ZERO;
          tree.clear_children(pnode->children[i]);
        }
      }
      break;

    case GA_NODE_PARAMS:
      if (child0->node_type == GA_NODE_PREDEF_FUNC) {
        std::string name = child0->name;
        ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
        const ga_predef_function &F = (it != PREDEF_FUNCTIONS.end()) ?
          it->second : workspace.user_function(name);

        GMM_ASSERT1(F.ftype == 0, "Not already taken into account");

        if (F.nbargs == 1) {
          child0->name = F.derivative1;
          tree.insert_node(pnode);
          pga_tree_node pnode_op = pnode->parent;
          pnode_op->node_type = GA_NODE_OP;
          if (child1->tensor_order() == 0)
            pnode_op->op_type = GA_MULT;
          else
            pnode_op->op_type = GA_DOTMULT;
          pnode_op->children.push_back(0);
          tree.copy_node(child1, pnode_op, pnode_op->children[1]);
          ga_node_derivation(tree, workspace, pnode_op->children[1],
                             varname, order);
        } else {
          pga_tree_node child2 = pnode->children[2];
          
          if (child1->marked && child2->marked)
            tree.duplicate_with_addition(pnode);
          
          if (child1->marked) {
            child0->name = F.derivative1;
            tree.insert_node(pnode);
            pga_tree_node pnode_op = pnode->parent;
            pnode_op->node_type = GA_NODE_OP;
            if (child1->tensor_order() == 0)
              pnode_op->op_type = GA_MULT;
            else
              pnode_op->op_type = GA_DOTMULT;
            pnode_op->children.push_back(0);
            tree.copy_node(child1, pnode_op, pnode_op->children[1]);
            ga_node_derivation(tree, workspace, pnode_op->children[1],
                               varname, order);
          }
          if (child2->marked) {
            if (child1->marked && child2->marked)
              pnode = pnode->parent->parent->children[1];
            child0 = pnode->children[0]; child1 = pnode->children[1];
            child2 = pnode->children[2];
            child0->name = F.derivative2;
            tree.insert_node(pnode);
            pga_tree_node pnode_op = pnode->parent;
            pnode_op->node_type = GA_NODE_OP;
            if (child2->tensor_order() == 0)
              pnode_op->op_type = GA_MULT;
            else
              pnode_op->op_type = GA_DOTMULT;
            pnode_op->children.push_back(0);
            tree.copy_node(child2, pnode_op, pnode_op->children[1]);
            ga_node_derivation(tree, workspace, pnode_op->children[1],
                               varname, order);
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
              tree.insert_node(pnode);
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
            tree.insert_node(pnode2);
            pga_tree_node pnode_op = pnode2->parent;
            pnode_op->node_type = GA_NODE_OP;
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
            tree.copy_node(pnode->children[i], pnode_op ,
                           pnode_op->children[1]);
            ga_node_derivation(tree, workspace, pnode_op->children[1], 
                               varname, order);
          }
        }

      } else {
        ga_node_derivation(tree, workspace, child0, varname, order);
      }
      break;

    default: GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                         << " in derivation. Internal error.");
    }
  }

  // The tree is modified. Should be copied first and passed to
  // ga_semantic_analysis after for enrichment
  static void ga_derivation(ga_tree &tree, const ga_workspace &workspace,
                     const std::string &varname, size_type order) {
    if (!(tree.root)) return;
    if (ga_node_mark_tree_for_variable(tree.root, varname))
      ga_node_derivation(tree, workspace, tree.root, varname, order);
    else
      tree.clear();
  }



  //=========================================================================
  // Compilation of assembly trees into a list of basic instructions
  //=========================================================================

  static void ga_compile_node(pga_tree_node pnode,
                              ga_workspace &workspace,
                              ga_instruction_set &gis,
                              ga_instruction_set::region_mim_instructions &rmi,
                              const mesh_im *mim) {
    
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_compile_node(pnode->children[i], workspace, gis, rmi, mim);

    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bgeot::multi_index mi;
    // const bgeot::multi_index &size0 = child0 ? child0->t.sizes() : mi;
    // const bgeot::multi_index &size1 = child1 ? child1->t.sizes() : mi;
    size_type dim0 = child0 ? child0->tensor_order() : 0;
    // size_type dim1 = child1 ? child1->tensor_order() : 0;

    pga_instruction pgai = 0;
    switch (pnode->test_function_type) {
    case 1:
      {
        const mesh_fem *mf = workspace.associated_mf(pnode->name_test1);
        pgai = new ga_instruction_first_ind_tensor
          (pnode->t, gis.ctx, pnode->qdim1, *mf);
        rmi.instructions.push_back(pgai);
      }
      break;
    case 2:
      {
        const mesh_fem *mf = workspace.associated_mf(pnode->name_test2);
        pgai = new ga_instruction_first_ind_tensor
          (pnode->t, gis.ctx, pnode->qdim2, *mf);
        rmi.instructions.push_back(pgai);
      }
      break;
    case 3:
      {
        const mesh_fem *mf1 = workspace.associated_mf(pnode->name_test1);
        const mesh_fem *mf2 = workspace.associated_mf(pnode->name_test2);
        pgai = new ga_instruction_two_first_ind_tensor
          (pnode->t, gis.ctx, pnode->qdim1, *mf1, pnode->qdim2, *mf2);
        rmi.instructions.push_back(pgai);
      }
      break;
    case 4:
      {
        const mesh_fem *mf1 = workspace.associated_mf(pnode->name_test1);
        const mesh_fem *mf2 = workspace.associated_mf(pnode->name_test2);
        pgai = new ga_instruction_two_first_ind_tensor
          (pnode->t, gis.ctx, pnode->qdim2, *mf2, pnode->qdim1, *mf1);
        rmi.instructions.push_back(pgai);
      }
    }
    
    switch (pnode->node_type) {

    case GA_NODE_PREDEF_FUNC: case GA_NODE_OPERATOR: case GA_NODE_SPEC_FUNC:
    case GA_NODE_CONSTANT: case GA_NODE_ALLINDICES: case GA_NODE_ZERO:
      break;

    case GA_NODE_VAL: case GA_NODE_GRAD: case GA_NODE_HESS:
      {
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        GMM_ASSERT1(mf, "Internal error");

        // An instruction for extracting local dofs of the variable.
        if (rmi.local_dofs.find(pnode->name) == rmi.local_dofs.end()) {
          rmi.local_dofs[pnode->name] = base_vector(1);
          
          if (gis.extended_vars.find(pnode->name) == gis.extended_vars.end()) {
            if (mf->is_reduced()) {
              base_vector U(mf->nb_basic_dof());
              mf->extend_vector(workspace.value(pnode->name), U);
              gis.really_extended_vars[pnode->name] = U;
              gis.extended_vars[pnode->name]
                = &(gis.really_extended_vars[pnode->name]);
            } else {
              gis.extended_vars[pnode->name]
                = &(workspace.value(pnode->name));
            }
          }
          pgai = new ga_instruction_slice_local_dofs
            (*mf, *(gis.extended_vars[pnode->name]), gis.ctx,
             rmi.local_dofs[pnode->name]);
          rmi.instructions.push_back(pgai);
        }

        // An instruction for pfp update
        if (rmi.pfps.find(mf) == rmi.pfps.end()) {
          rmi.pfps[mf] = 0;
          pgai = new ga_instruction_update_pfp
            (*mf,  rmi.pfps[mf], gis.ctx, gis.pai, gis.fp_pool);
          rmi.instructions.push_back(pgai);
        }

        // An intruction for the base value
        pgai = 0;
        if (pnode->node_type == GA_NODE_VAL) {
          if (rmi.base.find(mf) == rmi.base.end())
            pgai = new ga_instruction_base
              (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
        } else if (pnode->node_type == GA_NODE_GRAD) {
          if (rmi.grad.find(mf) == rmi.grad.end())
            pgai = new ga_instruction_grad_base
              (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
        } else {
          if (rmi.hess.find(mf) == rmi.hess.end())
            pgai = new ga_instruction_hess_base
              (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
        }
        if (pgai) rmi.instructions.push_back(pgai);
                          
        // The eval instruction
        if (pnode->node_type == GA_NODE_VAL)
          pgai = new ga_instruction_val
            (pnode->t, rmi.base[mf],
             rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
        else if (pnode->node_type == GA_NODE_GRAD)
          pgai = new ga_instruction_grad
            (pnode->t, rmi.base[mf],
             rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
        else
          pgai = new ga_instruction_hess
            (pnode->t, rmi.base[mf],
             rmi.local_dofs[pnode->name], workspace.qdim(pnode->name));
        rmi.instructions.push_back(pgai);
      }
      break;

    case GA_NODE_TEST: case GA_NODE_GRAD_TEST: case GA_NODE_HESS_TEST:
      {
        const mesh_fem *mf = workspace.associated_mf(pnode->name);
        if (mf) {

          // An instruction for pfp update
          if (rmi.pfps.find(mf) == rmi.pfps.end()) {
            rmi.pfps[mf] = 0;
            pgai = new ga_instruction_update_pfp
              (*mf,  rmi.pfps[mf], gis.ctx, gis.pai, gis.fp_pool);
            rmi.instructions.push_back(pgai);
          }

          // An intruction for the base value
          pgai = 0;
          if (pnode->node_type == GA_NODE_TEST) {
            if (rmi.base.find(mf) == rmi.base.end())
              pgai = new ga_instruction_base
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
          } else if (pnode->node_type == GA_NODE_GRAD_TEST) {
            if (rmi.grad.find(mf) == rmi.grad.end())
              pgai = new ga_instruction_grad_base
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
          } else {
            if (rmi.hess.find(mf) == rmi.hess.end())
              pgai = new ga_instruction_hess_base
                (rmi.base[mf], gis.ctx, *mf, rmi.pfps[mf]);
          }  
          if (pgai) rmi.instructions.push_back(pgai);
          
          // The copy of the real_base_value
          if (pnode->node_type == GA_NODE_TEST)
            pgai = new ga_instruction_copy_base
              (pnode->t, rmi.base[mf], mf->get_qdim());
          else if (pnode->node_type == GA_NODE_GRAD_TEST)
            pgai = new ga_instruction_copy_grad
              (pnode->t, rmi.base[mf], mf->get_qdim());
          else
            pgai = new ga_instruction_copy_hess
              (pnode->t, rmi.base[mf], mf->get_qdim());
          rmi.instructions.push_back(pgai);

          if (gis.var_intervals.find(pnode->name) == gis.var_intervals.end()) {
            size_type nd = mf->nb_basic_dof();
            gis.var_intervals[pnode->name] = gmm::sub_interval(gis.nb_dof, nd);
            gis.nb_dof += nd;
          }
        } else {
          if (gis.var_intervals.find(pnode->name) == gis.var_intervals.end()) {
            size_type nd = gmm::vect_size(workspace.value(pnode->name));
            gis.var_intervals[pnode->name] = gmm::sub_interval(gis.nb_dof, nd);
            gis.nb_dof += nd;
          }
        }
      }
      break;
      
     case GA_NODE_OP:
       switch(pnode->op_type) {

       case GA_PLUS:
         if (child0->test_function_type == child1->test_function_type) {
           pgai = new ga_instruction_add(pnode->t, child0->t, child1->t);
         } else if (child0->test_function_type > child1->test_function_type) {
           pgai = new ga_instruction_copy_tensor(pnode->t, child1->t);
           rmi.instructions.push_back(pgai);
           pgai = new ga_instruction_add_transp_test(pnode->t, child0->t);
         } else {
           pgai = new ga_instruction_copy_tensor(pnode->t, child0->t);
           rmi.instructions.push_back(pgai);
           pgai = new ga_instruction_add_transp_test(pnode->t, child1->t);
         }
         rmi.instructions.push_back(pgai);
         break;

       case GA_MINUS:
         if (child0->test_function_type == child1->test_function_type) {
           pgai = new ga_instruction_sub(pnode->t, child0->t, child1->t);
         } else if (child0->test_function_type > child1->test_function_type) {
           pgai = new ga_instruction_minus(pnode->t, child1->t);
           rmi.instructions.push_back(pgai);
           pgai = new ga_instruction_add_transp_test(pnode->t, child0->t);
         } else {
           pgai = new ga_instruction_minus(pnode->t, child0->t);
           rmi.instructions.push_back(pgai);
           pgai = new ga_instruction_add_transp_test(pnode->t, child1->t);
           rmi.instructions.push_back(pgai);
           pgai = new ga_instruction_opposite(pnode->t);
         }
         rmi.instructions.push_back(pgai);
         break;

       case GA_UNARY_MINUS:
         pgai = new ga_instruction_minus(pnode->t, child0->t);
         rmi.instructions.push_back(pgai);
         break;

       case GA_COLON: 
         {
           size_type s1 = (child0->t.size()*child1->t.size())/pnode->t.size();
           size_type s2 = size_type(round(sqrt(scalar_type(s1))));
           GMM_ASSERT1(s1 == s2*s2, "Very strange !");
           if (pnode->test_function_type == 3 && dim0 <= 2 &&
               child1->test_function_type == 2 &&
               child0->test_function_type == 1)
             pgai = new ga_instruction_reduction
               (pnode->t, child1->t, child0->t, s2);
           else
             pgai = new ga_instruction_reduction
               (pnode->t, child0->t, child1->t, s2);
           rmi.instructions.push_back(pgai);
         }
         break;

       case GA_DOT:
         {
           size_type s1 = (child0->t.size()*child1->t.size())/pnode->t.size();
           size_type s2 = size_type(round(sqrt(scalar_type(s1))));
           GMM_ASSERT1(s1 == s2*s2, "Very strange !");
           if (pnode->test_function_type == 3 && dim0 <= 1 &&
              child1->test_function_type == 2 &&
              child0->test_function_type == 1)
             pgai = new ga_instruction_reduction
               (pnode->t, child1->t, child0->t, s2);
           else
             pgai = new ga_instruction_reduction
               (pnode->t, child0->t, child1->t, s2);
           rmi.instructions.push_back(pgai);
         }
         break;

       case GA_MULT:
         {
           size_type s1 = (child0->t.size()*child1->t.size())/pnode->t.size();
           size_type s2 = size_type(round(sqrt(scalar_type(s1))));
           GMM_ASSERT1(s1 == s2*s2, "Very strange !");
           if (pnode->test_function_type == 3 && dim0 == 0 &&
               child1->test_function_type == 2 &&
              child0->test_function_type == 1)
             pgai = new ga_instruction_reduction
               (pnode->t, child1->t, child0->t, s2);
           else
             pgai = new ga_instruction_reduction
               (pnode->t, child0->t, child1->t, s2);
           rmi.instructions.push_back(pgai);
         }
         break;

        case GA_PRINT:
         pgai = new ga_instruction_copy_tensor(pnode->t, child0->t);
         rmi.instructions.push_back(pgai);
         pgai = new ga_instruction_print_tensor(pnode->t, child0, gis.ctx,
                                                gis.nbpt, gis.ipt);
         rmi.instructions.push_back(pgai);
         break;


//       case GA_DOTMULT: case GA_DOTDIV:
//         {
//           bool compatible = true;
//           if (child0->tensor_proper_size() != child1->tensor_proper_size())
//             compatible = false;

//           if (child0->tensor_proper_size() != 1) {
//             if (child0->tensor_order() != child1->tensor_order())
//               compatible = false;
            
//             for (size_type i = 0; i < child0->tensor_order(); ++i)
//               if (child0->tensor_proper_size(i)!=child1->tensor_proper_size(i))
//                 compatible = false;
//           }

//           if (!compatible)
//             ga_throw_error(expr, pnode->pos, "Arguments of different sizes"
//                            "for .* or ./");
          
//           pnode->mult_test(child0, child1, expr, workspace);
//           mi = pnode->t.sizes();
//           for (size_type i = 0; i < child0->tensor_order(); ++i)
//             mi.push_back(child0->tensor_proper_size(i));
//           pnode->t.adjust_sizes(mi);
          
//           if (all_cte) {
//             pnode->node_type = GA_NODE_CONSTANT;
//             pnode->test_function_type = 0;
//             if (pnode->op_type == GA_DOTMULT) {
//               for (size_type i = 0; i < child0->t.size(); ++i)
//                 pnode->t[i] = child0->t[i] * child1->t[i];
//             } else {
//               for (size_type i = 0; i < child0->t.size(); ++i) {
//                 if (child1->t[i] == scalar_type(0))
//                   ga_throw_error(expr, pnode->pos, "Division by zero.");
//                 pnode->t[i] = child0->t[i] / child1->t[i];
//               }
//             }
//             tree.clear_children(pnode);
//           } else {
//             if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
//               gmm::clear(pnode->t.as_vector());
//               pnode->node_type = GA_NODE_ZERO;
//               tree.clear_children(pnode);
//             }
//             if (child1->tensor_is_zero() && pnode->op_type == GA_DOTDIV)
//               ga_throw_error(expr, pnode->pos, "Division by zero.");
//           }
//         }
//         break;

//       case GA_QUOTE:
//         if (dim0 > 2)
//           ga_throw_error(expr, pnode->pos, "Transpose operator is for "
//                          "vectors or matrices only.");
//         mi = size0;
//         if (child0->tensor_proper_size() == 1)
//           { tree.replace_node_by_child(pnode, 0); break; }
//         else if (dim0 == 2) std::swap(mi.back(), mi[size0.size()-2]);
//         else { size_type N = mi.back(); mi.back() = 1; mi.push_back(N); }

//         pnode->t.adjust_sizes(mi);
//         pnode->test_function_type = child0->test_function_type;
//         pnode->name_test1 = child0->name_test1;
//         pnode->name_test2 = child0->name_test2;
//         pnode->qdim1 = child0->qdim1;
//         pnode->qdim2 = child0->qdim2;
//         if (all_cte) {
//           pnode->node_type = GA_NODE_CONSTANT;
//           pnode->test_function_type = 0;
//           if (dim0 == 2) {
//             for (size_type i = 0; i < mi.back(); ++i)
//               for (size_type j = 0; j < mi[size0.size()-2]; ++j)
//                 pnode->t(j, i) = child0->t(i,j);
//           } else if (dim0 == 1) {
//             for (size_type i = 0; i < mi.back(); ++i)
//               pnode->t(0, i) = child0->t[i];
//           }
//           tree.clear_children(pnode);
//         } else if (child0->node_type == GA_NODE_ZERO) {
//           pnode->node_type = GA_NODE_ZERO;
//           gmm::clear(pnode->t.as_vector());
//           tree.clear_children(pnode);
//         }
//         break;

//       case GA_TRACE:
//         {
//           mi = size0;
//           size_type N = (child0->tensor_proper_size() == 1) ? 1 : mi.back();
          
//           if ((dim0 != 2 && child0->tensor_proper_size() != 1) ||
//               (dim0 == 2 && mi[mi.size()-2] != N))
//             ga_throw_error(expr, pnode->pos,
//                            "Trace operator is for square matrices only.");
          
//           if (dim0 == 2) { mi.pop_back(); mi.pop_back(); }
//           pnode->t.adjust_sizes(mi);
//           pnode->test_function_type = child0->test_function_type;
//           pnode->name_test1 = child0->name_test1;
//           pnode->name_test2 = child0->name_test2;
//           pnode->qdim1 = child0->qdim1;
//           pnode->qdim2 = child0->qdim2;
//           if (all_cte) {
//             pnode->node_type = GA_NODE_CONSTANT;
//             pnode->test_function_type = 0;
//             if (dim0 == 2) {
//               pnode->t[0] = scalar_type(0);
//               for (size_type i = 0; i < N; ++i)
//                 pnode->t[0] += child0->t(i,i);
//             } else {
//               pnode->t[0] += child0->t[0];
//             }
//             tree.clear_children(pnode);
//           } else if (child0->node_type == GA_NODE_ZERO) {
//             pnode->node_type = GA_NODE_ZERO;
//             gmm::clear(pnode->t.as_vector());
//             tree.clear_children(pnode);
//           }
//         }
//         break;

//       case GA_TMULT:
//         if (all_cte) {
//           pnode->node_type = GA_NODE_CONSTANT;
//           pnode->test_function_type = 0;
//           if (child0->t.size() == 1 && child1->t.size() == 1) {
//             pnode->init_scalar_tensor
//               (child0->t[0] * child1->t[0]);
//           } else if (child0->t.size() == 1) {
//             pnode->t = child1->t;
//             gmm::scale(pnode->t.as_vector(), scalar_type(child0->t[0]));
//           } else if (child1->t.size() == 1) {
//             pnode->t = child0->t;
//             gmm::scale(pnode->t.as_vector(), scalar_type(child1->t[0]));
//           } else {
//             if (dim0+dim1 > 4)
//               ga_throw_error(expr, pnode->pos, "Unauthorized "
//                               "tensor multiplication.");
//             for (size_type i = 0; i < dim0; ++i)
//               mi.push_back(child0->t.size(i));
//             for (size_type i = 0; i < dim1; ++i)
//               mi.push_back(child1->t.size(i));
//             pnode->t.adjust_sizes(mi);
//             size_type n0 = child0->t.size();
//             size_type n1 = child1->t.size();
//             for (size_type i = 0; i < n0; ++i)
//               for (size_type j = 0; j < n1; ++j)
//                 pnode->t[i+j*n0] = child0->t[i] * child1->t[j];
//           }
//           tree.clear_children(pnode);
//         } else {
//           pnode->mult_test(child0, child1, expr, workspace);
//           mi = pnode->t.sizes();
//           if (child0->tensor_proper_size() != 1
//               || child1->tensor_proper_size() != 1) {
//             if (child0->tensor_proper_size() == 1) {
//               for (size_type i = 0; i < dim1; ++i)
//                 mi.push_back(child1->tensor_proper_size(i));
//             } else if (child1->t.size() == 1) {
//               for (size_type i = 0; i < dim0; ++i)
//                 mi.push_back(child0->tensor_proper_size(i));
//             } else {
//               if (dim0+dim1 > 4)
//                 ga_throw_error(expr, pnode->pos, "Unauthorized "
//                                 "tensor multiplication.");
//               for (size_type i = 0; i < dim0; ++i)
//                 mi.push_back(child0->t.size(i));
//               for (size_type i = 0; i < dim1; ++i)
//                 mi.push_back(child1->t.size(i));
//             }
//             pnode->t.adjust_sizes(mi);
//           }
//           if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
//             gmm::clear(pnode->t.as_vector());
//             pnode->node_type = GA_NODE_ZERO;
//             tree.clear_children(pnode);
//           }
//         }
//         break;
        


//       case GA_DIV:
//         if (child1->tensor_order() > 0)
//           ga_throw_error(expr, pnode->pos, "Only the division by a scalar "
//                           "is allowed.");
//         if (child1->test_function_type)
//           ga_throw_error(expr, pnode->pos, "Division by test functions "
//                           "is not allowed.");
//         if (child1->node_type == GA_NODE_CONSTANT &&
//             child1->t[0] == scalar_type(0))
//           ga_throw_error(expr, pnode->children[1]->pos, "Division by zero");
        
//         pnode->t = child0->t;
//         pnode->test_function_type = child0->test_function_type;
//         pnode->name_test1 = child0->name_test1;
//         pnode->name_test2 = child0->name_test2;
//         pnode->qdim1 = child0->qdim1;
//         pnode->qdim2 = child0->qdim2;
        
//         if (all_cte) {          
//           pnode->node_type = GA_NODE_CONSTANT;
//           pnode->t = pnode->children[0]->t;
//           pnode->test_function_type = 0;
//           gmm::scale(pnode->t.as_vector(),
//                      scalar_type(1) / pnode->children[1]->t[0]);
//           tree.clear_children(pnode);
//         }
//         if (child0->tensor_is_zero()) {
//           gmm::clear(pnode->t.as_vector());
//           pnode->node_type = GA_NODE_ZERO;
//           tree.clear_children(pnode);
//         } else if (child1->node_type == GA_NODE_CONSTANT &&
//                    child1->t.size() == 1 && child1->t[0] == scalar_type(1)) {
//           tree.replace_node_by_child(pnode, 0);
//         }
//         break;
        
       default:GMM_ASSERT1(false, "Unexpected operation. Internal error.");
       }
       break;

    case GA_NODE_C_MATRIX:
      // TODO: bien penser que dans le cas de fonctions test doubles, certaines
      //       composantes peuvent être inversées
      // TOBEDONE !!
      break;
//       {
//         if (!all_sc) 
//           ga_throw_error(expr, pnode->pos, "Constant vector/matrix/tensor "
//                           "components should be scalar valued.");
// //         if (!all_primal)
// //           ga_throw_error(expr, pnode->pos, "Test functions are not allowed "
// //                           "in constant vector/matrix/tensor components.");
        
//         size_type nbc1 = pnode->nbc1, nbc2 = pnode->nbc2, nbc3 = pnode->nbc3;
//         size_type nbl = pnode->children.size() / (nbc1*nbc2*nbc3);
//         if (all_cte) pnode->node_type = GA_NODE_CONSTANT;
//         pnode->test_function_type = 0;
//         for (size_type i = 0; i < pnode->children.size(); ++i) {
//           if (pnode->children[i]->test_function_type) {
//             if (pnode->test_function_type == 0) {
//               pnode->test_function_type=pnode->children[i]->test_function_type;
//               pnode->name_test1 = pnode->children[i]->name_test1;
//               pnode->name_test2 = pnode->children[i]->name_test2;
//               pnode->qdim1 = pnode->children[i]->qdim1;
//               pnode->qdim2 = pnode->children[i]->qdim2;
//             } else {
//               if (pnode->test_function_type
//                   != pnode->children[i]->test_function_type ||
//                   pnode->name_test1.compare(pnode->children[i]->name_test1)
//                   != 0 ||
//                   pnode->name_test2.compare(pnode->children[i]->name_test2)
//                   != 0)
//                 ga_throw_error(expr, pnode->pos, "Inconsistent use of test "
//                                "function in constant matrix.");
//             }
//           }
//         }
//         mi.resize(0);
//         if (pnode->test_function_type) mi.push_back(1);
//         if (pnode->test_function_type == 3) mi.push_back(1);
//         if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1 && nbl == 1) {
//           pnode->t.adjust_sizes(mi);
//           if (all_cte) pnode->t[0] = child0->t[0];
//         } else if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
//           mi.push_back(nbl);
//           pnode->t.adjust_sizes(mi);
//           if (all_cte)
//             for (size_type i = 0; i < nbl; ++i)
//               pnode->t[i] = pnode->children[i]->t[0];
//         } else if (nbc2 == 1 && nbc3 == 1) {
//           mi.push_back(nbl); mi.push_back(nbc1); 
//           pnode->t.adjust_sizes(mi);
//           if (all_cte)
//             for (size_type i = 0; i < nbl; ++i)
//               for (size_type j = 0; j < nbc1; ++j)
//                 pnode->t(i,j) = pnode->children[i*nbc1+j]->t[0];
//         } else {
//           mi.push_back(nbl); mi.push_back(nbc3);
//           mi.push_back(nbc2); mi.push_back(nbc1);
//           pnode->t.adjust_sizes(mi);
//           size_type n = 0;
//           if (all_cte)
//             for (size_type i = 0; i < nbl; ++i)
//               for (size_type j = 0; j < nbc3; ++j)
//                 for (size_type k = 0; k < nbc2; ++k)
//                   for (size_type l = 0; l < nbc1; ++l)
//                     pnode->t(i,j,k,l) = pnode->children[n++]->t[0];
//         }
//         if (all_cte) tree.clear_children(pnode);
//       }
//       break;


//     case GA_NODE_PARAMS:

//       if (child0->node_type == GA_NODE_PREDEF_FUNC) {
        
//         // Évaluation of predefined function
        
//         for (size_type i = 1; i < pnode->children.size(); ++i)
//           ga_valid_operand(expr, pnode->children[i]);
//         std::string name = child0->name;
//         ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
//         const ga_predef_function &F = (it != PREDEF_FUNCTIONS.end()) ?
//           it->second : workspace.user_function(name);
//         size_type nbargs = F.nbargs;
//         GMM_ASSERT1(F.ftype != 2, "Not already taken into account");
//         if (nbargs+1 != pnode->children.size()) {
//             std::stringstream msg;
//             msg << "Bad number of arguments for predefined function "
//                 << name << ". Found " << pnode->children.size()-1
//                 << " should be " << nbargs << ".";
//             ga_throw_error(expr, pnode->pos, msg.str());
//         }
//         pnode->test_function_type = 0;
//         pga_tree_node child2 = (nbargs == 2) ? pnode->children[2] : child1;
//         all_cte = child1->node_type == GA_NODE_CONSTANT;
//         if (nbargs == 2)
//           all_cte = all_cte && (child2->node_type == GA_NODE_CONSTANT);
//         if (child1->test_function_type || child2->test_function_type)
//           ga_throw_error(expr, pnode->pos, "Test functions cannot be passed "
//                          "as argument of a predefined function.");
//         if (child1->tensor_order() > 2 || child2->tensor_order() > 2)
//           ga_throw_error(expr, pnode->pos, "Sorry, function can be applied "
//                          "to scalar, vector and matrices only.");
//         size_type s1 = child1->t.size();
//         size_type s2 = (nbargs == 2) ? child2->t.size() : s1;
//         if (s1 != s2 && (s1 != 1 || s2 != 1))
//           ga_throw_error(expr, pnode->pos,
//                          "Invalid argument size for a scalar function.");

//         if (nbargs == 1) {
//           pnode->t = child1->t;
//         } else {
//           if (s1 == s2) {
//             pnode->t = child1->t;
//           } else if (s1 == 1) {
//             pnode->t = child2->t;
//           } else {
//             pnode->t = child1->t;
//           }
//         }

//         if (all_cte) {
//           if (pnode->der1)
//             GMM_ASSERT1(false, "Sorry, to be done");
//           pnode->node_type = GA_NODE_CONSTANT;
//           if (nbargs == 1) {
//             for (size_type i = 0; i < s1; ++i)
//               pnode->t[i] = (*(F.f1))(child1->t[i]);
//           } else {
//             if (s1 == s2) {
//               for (size_type i = 0; i < s1; ++i)
//                 pnode->t[i] = (*(F.f2))(child1->t[i], child2->t[i]);
//             } else if (s1 == 1) {
//               for (size_type i = 0; i < s2; ++i)
//                 pnode->t[i] = (*(F.f2))(child1->t[0], child2->t[i]);
//             } else {
//               for (size_type i = 0; i < s1; ++i)
//                 pnode->t[i] = (*(F.f2))(child1->t[i], child2->t[0]);
//             }
//           }
//           tree.clear_children(pnode);
//         }
//       } else if (child0->node_type == GA_NODE_SPEC_FUNC) {

//         // Special constant functions: meshdim(u), qdim(u) ...
        
//         for (size_type i = 1; i < pnode->children.size(); ++i)
//           ga_valid_operand(expr, pnode->children[i]);
//         if (pnode->children.size() != 2)
//           ga_throw_error(expr, pnode->pos,
//                          "One and only one argument is allowed for function "
//                          +child0->name+".");

//         if (!(child0->name.compare("meshdim"))) {
//           bool valid = (child1->node_type == GA_NODE_VAL);
//           const mesh_fem *mf = valid ? workspace.associated_mf(child1->name):0;
//           if (!mf)
//             ga_throw_error(expr, pnode->pos, "The argument of mesh_dim "
//                            "function can only be a fem variable name.");
//           pnode->node_type = GA_NODE_CONSTANT;
//           pnode->init_scalar_tensor(scalar_type(mf->linked_mesh().dim()));
//         } else if (!(child0->name.compare("qdim"))) {
//           if (child1->node_type != GA_NODE_VAL)
//             ga_throw_error(expr, pnode->pos, "The argument of qdim "
//                            "function can only be a variable name.");
//           pnode->node_type = GA_NODE_CONSTANT;
//           pnode->init_scalar_tensor(scalar_type(workspace.qdim(child1->name)));
//         } else if (!(child0->name.compare("Id"))) {
//           bool valid = (child1->node_type == GA_NODE_CONSTANT);
//           int n = valid ? int(round(child1->t[0])) : -1;
//           if (n <= 0 || n > 100 || child1->tensor_order() > 0)
//             ga_throw_error(expr, pnode->pos, "The argument of Id "
//                            "should be a (small) positive integer.");
//           pnode->node_type = GA_NODE_CONSTANT;
//           if (n == 1)
//             pnode->init_scalar_tensor(scalar_type(1));
//           else {
//             pnode->init_matrix_tensor(n,n);
//             for (int i = 0; i < n; ++i) pnode->t(i,i) = scalar_type(1);
//           }
//         } else ga_throw_error(expr, pnode->children[0]->pos,
//                               "Unknown special function.");
//         tree.clear_children(pnode);
//       } else if (child0->node_type == GA_NODE_OPERATOR) {

//         // Call to a nonlinear operator

//         for (size_type i = 1; i < pnode->children.size(); ++i)
//           ga_valid_operand(expr, pnode->children[i]);
//         all_cte = true;
//         ga_nonlinear_operator::arg_list args;
//         for (size_type i = 1; i < pnode->children.size(); ++i) {
//           all_cte = all_cte
//             && (pnode->children[i]->node_type == GA_NODE_CONSTANT);
//           args.push_back(&(pnode->children[i]->t));
//           if (pnode->children[i]->node_type == GA_NODE_ALLINDICES)
//             ga_throw_error(expr, pnode->children[i]->pos,
//                            "Colon operator is not allowed in nonlinear "
//                            "operator call.");
//           if (pnode->children[i]->test_function_type)
//             ga_throw_error(expr, pnode->pos, "Test functions cannot be passed "
//                            "as argument of a nonlinear operator.");
//           if (pnode->children[i]->tensor_order() > 2)
//             ga_throw_error(expr, pnode->pos, "Sorry, arguments to nonlinear "
//                         "operators should only be scalar, vector or matrices");
//         }
//         ga_predef_operator_tab::iterator it
//           = PREDEF_OPERATORS.find(child0->name);
//         const ga_nonlinear_operator &OP = (it != PREDEF_OPERATORS.end()) ?
//           *(it->second) : workspace.user_operator(child0->name);
//         mi.resize(0);
//         if (!(OP.result_size(args, mi)))
//           ga_throw_error(expr, pnode->pos,
//                          "Wrong number or wrong size of arguments for the "
//                          "call of nonlinear operator " + child0->name);
        
//         pnode->test_function_type = 0;

//         if (child0->der1 > args.size() || child0->der2 > args.size())
//            ga_throw_error(expr, child0->pos,
//                          "Invalid derivative number for nonlinear operator "
//                           + child0->name);

//         if (child0->der1 && child0->der2 == 0) {
//           for (size_type i = 0; i < args[child0->der1-1]->sizes().size(); ++i)
//             mi.push_back(args[child0->der1-1]->sizes()[i]);
//           pnode->t.adjust_sizes(mi);
//           if (all_cte) {
//             pnode->node_type = GA_NODE_CONSTANT;
//             OP.derivative(args, child0->der1, pnode->t);
//             tree.clear_children(pnode);
//           }
//         } else if (child0->der1 && child0->der2) {
//           for (size_type i = 0; i < args[child0->der1-1]->sizes().size(); ++i)
//             mi.push_back(args[child0->der1-1]->sizes()[i]);
//           for (size_type i = 0; i < args[child0->der2-1]->sizes().size(); ++i)
//             mi.push_back(args[child0->der2-1]->sizes()[i]);
//           pnode->t.adjust_sizes(mi);
//           if (all_cte) {
//             pnode->node_type = GA_NODE_CONSTANT;
//             OP.second_derivative(args, child0->der1, child0->der2,
//                                           pnode->t);
//             tree.clear_children(pnode);
//           }
//         } else {
//           pnode->t.adjust_sizes(mi);
//           if (all_cte) {
//             pnode->node_type = GA_NODE_CONSTANT;
//             OP.value(args, pnode->t);
//             tree.clear_children(pnode);
//           }
//         }
//       } else {

//         // Access to components of a tensor
        
//         all_cte = (child0->node_type == GA_NODE_CONSTANT);
//         if (pnode->children.size() != child0->tensor_order() + 1)
//           ga_throw_error(expr, pnode->pos, "Bad number of indices.");
//         for (size_type i = 1; i < pnode->children.size(); ++i)
//           if (pnode->children[i]->node_type != GA_NODE_ALLINDICES &&
//               (pnode->children[i]->node_type != GA_NODE_CONSTANT ||
//                pnode->children[i]->t.size() != 1))
//             ga_throw_error(expr, pnode->children[i]->pos,
//                             "Indices should be constant integers or colon.");
        
//         bgeot::multi_index mi1(size0.size()), mi2, indices;
//         for (size_type i = 0; i < child0->tensor_order(); ++i) {
//           if (pnode->children[i+1]->node_type == GA_NODE_ALLINDICES) {
//             mi2.push_back(child0->tensor_proper_size(i));
//             indices.push_back(i);
//             mi1[i] = 0;
//           } else {
//             mi1[i] = size_type(round(pnode->children[i+1]->t[0])-1);
//             if (mi1[i] >= child0->tensor_proper_size(i))
//               ga_throw_error(expr, pnode->children[i+1]->pos,
//                               "Index out of range.");
//           }
//         }
//         mi = size0; mi.resize(child0->nb_test_functions());
//         for (size_type i = 0; i < mi2.size(); ++i) mi.push_back(mi2[i]);
//         pnode->t.adjust_sizes(mi);
//         pnode->test_function_type = child0->test_function_type;
//         pnode->name_test1 = child0->name_test1;
//         pnode->name_test2 = child0->name_test2;
//         pnode->qdim1 = child0->qdim1;
//         pnode->qdim2 = child0->qdim2;
        
//         if (all_cte) {
//           pnode->node_type = GA_NODE_CONSTANT;
//           for (bgeot::multi_index mi3(mi2.size()); !mi3.finished(mi2);
//                mi3.incrementation(mi2)) {
//             for (size_type j = 0; j < mi2.size(); ++j) {
//               mi1[indices[j]] = mi3[j];
//             }
//             pnode->t(mi3) = pnode->children[0]->t(mi1);
//           }
//           tree.clear_children(pnode);
//         } else {
//           if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
//             gmm::clear(pnode->t.as_vector());
//             pnode->node_type = GA_NODE_ZERO;
//             tree.clear_children(pnode);
//           }
//         }
//       }
//       break;

    default:GMM_ASSERT1(false, "Unexpected node type " << pnode->node_type
                        << " in compilation. Internal error.");
    }
  }
  

  static void ga_compile(ga_workspace &workspace, ga_instruction_set &gis,
                         size_type order) {
    for (size_type i = 0; i < workspace.nb_trees(); ++i) {
      ga_workspace::tree_description &td = workspace.tree_info(i);
      if (td.order == order) {
        gis.trees.push_back(*(td.ptree));
        
        // Semantic analysis mainly to evaluate fixed size variables and data
        ga_semantic_analysis("", gis.trees.back(), workspace, true);
        pga_tree_node root = gis.trees.back().root;
        if (root) {
          cout << "compiling tree: "; ga_print_tree(gis.trees.back());
          // Compiling tree
          ga_instruction_set::region_mim rm(td.mim, td.region);
          ga_compile_node(root, workspace, gis,
                          gis.whole_instructions[rm], td.mim);
         
          // Addition of an assembly instruction
          pga_instruction pgai = 0;
          switch(order) {
          case 0:
            pgai = new ga_instruction_scalar_assembly(root->t, workspace.E,
                                                      gis.coeff);
            break;
          case 1:
            GMM_ASSERT1(gis.var_intervals.find(root->name_test1)
                        != gis.var_intervals.end(), "Internal error");
            pgai = new ga_instruction_vector_assembly
              (root->t, workspace.V, gis.ctx,
               gis.var_intervals[root->name_test1],
               *(workspace.associated_mf(root->name_test1)), gis.coeff);
            break;
          case 2:
            if (gis.trees.back().root->test_function_type == 4)
              pgai = new ga_instruction_transp_matrix_assembly
                (root->t, workspace.K, gis.ctx,
                 gis.var_intervals[root->name_test1],
                 gis.var_intervals[root->name_test2],
                 *(workspace.associated_mf(root->name_test1)),
                 *(workspace.associated_mf(root->name_test2)),
                 gis.coeff, gis.nbpt, gis.ipt, td.elem);
            else
              pgai = new ga_instruction_matrix_assembly
                (root->t, workspace.K, gis.ctx,
                 gis.var_intervals[root->name_test1],
                 gis.var_intervals[root->name_test2],
                 *(workspace.associated_mf(root->name_test1)),
                 *(workspace.associated_mf(root->name_test2)),
                 gis.coeff, gis.nbpt, gis.ipt, td.elem);
            break;
          }
          if (pgai) gis.whole_instructions[rm].instructions.push_back(pgai);
        }
      }
    }
  }
  

  //=========================================================================
  // Execution of a compiled set of assembly terms
  //=========================================================================

  
  static void ga_exec(ga_instruction_set &gis) {
    base_matrix G;
    base_small_vector un, up;

    ga_instruction_set::instructions_set::iterator it
      = gis.whole_instructions.begin();
    for (; it != gis.whole_instructions.end(); ++it) {

      const getfem::mesh_im &mim = *(it->first.first);
      const getfem::mesh &mesh = mim.linked_mesh();
      size_type P = mesh.dim();
      ga_instruction_list &gil = it->second.instructions;
      mesh_region region(it->first.second);
      
      // iteration on elements (or faces of elements)
      for (getfem::mr_visitor v(region, mesh); !v.finished(); ++v) {
        if (mim.convex_index().is_in(v.cv())) {
          // cout << "proceed with element " << v.cv() << endl;
          bgeot::vectors_to_base_matrix(G, mesh.points_of_convex(v.cv()));
          size_type N = G.nrows();
          bgeot::pgeometric_trans pgt = mesh.trans_of_convex(v.cv());
          pintegration_method pim = mim.int_method_of_element(v.cv());
          GMM_ASSERT1(pim->type() == IM_APPROX, "Sorry, exact methods cannot "
                      "be used in high level generic assembly");
          gis.pai = pim->approx_method();
          const bgeot::stored_point_tab &spt = gis.pai->integration_points();
          if (spt.size()) {
            if (gis.ctx.have_pgp() && gis.ctx.pgt() == pgt) {
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
          
            // Computation of unit normal vector in case of a boundary
            const base_matrix& B = gis.ctx.B();
            scalar_type J = gis.ctx.J();
            if (v.f() != short_type(-1)) {
              up.resize(N); un.resize(P);
              gmm::copy(pgt->normals()[v.f()], un);
              gmm::mult(B, un, up);
              scalar_type nup = gmm::vect_norm2(up);
              J *= nup;
              gmm::scale(up,1.0/nup); // up sera utile plus tard ...
            }

            // iterations on Gauss points
            gis.nbpt = gis.pai->nb_points_on_convex();
            size_type first_ind = 0;
            if (v.f() != short_type(-1)) {
              gis.nbpt = gis.pai->nb_points_on_face(v.f());
              first_ind = gis.pai->ind_first_point_on_face(v.f());
            }
            for (gis.ipt = 0; gis.ipt < gis.nbpt; ++(gis.ipt)) {
              gis.coeff = J * gis.pai->coeff(first_ind+gis.ipt);

              // cout << "proceed with Gauss point " << i <<"/"<< nbpt << endl;
              if (gis.ctx.have_pgp()) gis.ctx.set_ii(first_ind+gis.ipt);
              else gis.ctx.set_xref(spt[first_ind+gis.ipt]);
              for (size_type j = 0; j < gil.size(); ++j) gil[j]->exec();
            }
          }
        }
      }
    }
  }



} /* end of namespace */
