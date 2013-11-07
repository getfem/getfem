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

#include "getfem/getfem_models.h"
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
    GA_DOT,         // '.'
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
    ga_char_type[':'] = GA_COLON;       ga_char_type['.'] = GA_DOT;
    ga_char_type['@'] = GA_TMULT;       ga_char_type[','] = GA_COMMA;
    ga_char_type[';'] = GA_SEMICOLON;   ga_char_type['('] = GA_LPAR;
    ga_char_type[')'] = GA_RPAR;        ga_char_type['['] = GA_LBRACKET;
    ga_char_type[']'] = GA_RBRACKET;    ga_char_type['_'] = GA_NAME;
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
    ga_operator_priorities[GA_TMULT] = 2;
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

    // Ignore white spaces
    while (expr[pos] == ' ' && pos < expr.size()) ++pos;
    token_pos = pos;
    token_length = 0;

    // Detecting end of expression
    if (pos >= expr.size()) return GA_END;

    // Treating the different cases (Operation, name or number)
    GA_TOKEN_TYPE type = ga_char_type[unsigned(expr[pos++])];
    ++token_length;
    switch (type) {
    case GA_DOT:
      if (pos >= expr.size() || ga_char_type[unsigned(expr[pos])] != GA_SCALAR)
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

  static void ga_throw_error(const std::string &expr, size_type pos,
                              const std::string &msg) {
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
    cerr << '|' << endl << msg << endl;
    GMM_ASSERT1(false, "Error in assembly string" );
  }

  enum GA_NODE_TYPE {
    GA_NODE_VOID = 0,
    GA_NODE_OP,
    GA_NODE_PREDEF_FUNC,
    GA_NODE_CONSTANT,
    GA_NODE_NAME,
    GA_NODE_PARAMS,
    GA_NODE_ALLINDICES,
    GA_NODE_C_MATRIX,
    GA_NODE_VAL,
    GA_NODE_GRAD,
    GA_NODE_HESS,
    GA_NODE_FIXED_SIZE_TEST,
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
                                  // 0 = no test function, 1 = first order
                                  // 2 = second order, 3 = both
    std::string name_test0, name_test1; // variable names corresponding to
                              // test functions when test_function_type > 0.
    size_type qdim0, qdim1;   // Qdims when test_function_type > 0.
    size_type nbc1, nbc2, nbc3, pos;
    std::string name;             // variable/constant/function/operator name 
    GA_TOKEN_TYPE op_type;
    pga_tree_node parent; // only one for the moment ...
    std::vector<pga_tree_node> children;
    
    inline size_type nb_test_functions(void) const {
      if (test_function_type == size_type(-1)) return 0;
      return test_function_type - (test_function_type >= 2 ? 1 : 0);
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
                             test0 == 3 || test1 == 3))
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
      t.adjust_sizes(mi);
      
      if (n0->name_test0.size())
        { name_test0 = n0->name_test0; qdim0 = n0->qdim0; }
      else
        { name_test0 = n1->name_test0; qdim0 = n1->qdim0; }
      
      if (n0->name_test1.size())
        { name_test1 = n0->name_test1; qdim1 = n0->qdim1; }
      else
        { name_test1 = n1->name_test1; qdim1 = n1->qdim1; }
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

    ga_tree_node(void): node_type(GA_NODE_VOID), test_function_type(-1),
                        qdim0(0), qdim1(0)  {}
    ga_tree_node(GA_NODE_TYPE ty, size_type p)
      : node_type(ty), test_function_type(-1), qdim0(0), qdim1(0), pos(p) {}
    ga_tree_node(scalar_type v, size_type p)
      : node_type(GA_NODE_CONSTANT), test_function_type(-1), qdim0(0),
        qdim1(0), pos(p)
    { init_scalar_tensor(v); }
    ga_tree_node(const char *n, size_type l, size_type p)
      : node_type(GA_NODE_NAME), test_function_type(-1), qdim0(0), qdim1(0),
        pos(p), name(n, l) {}
    ga_tree_node(GA_TOKEN_TYPE op, size_type p)
      : node_type(GA_NODE_OP), test_function_type(-1), qdim0(0), qdim1(0),
        pos(p), op_type(op) {}
    
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
        if (op_type == GA_UNARY_MINUS) {
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

    void clear_node(pga_tree_node pnode) {
      for (size_type i = 0; i < pnode->children.size(); ++i)
        clear_node(pnode->children[i]);
      delete pnode;
      current_node = 0;
    }

    void clear(void)
    { if (root) clear_node(root); root = current_node = 0; }
    
    void clear_children(pga_tree_node pnode) {
      for (size_type i = 0; i < pnode->children.size(); ++i)
        clear_node(pnode->children[i]);
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
        if (j != i) clear_node(pnode->children[j]);
      delete pnode;
    }

    ga_tree(void) : root(0), current_node(0) {}
    ~ga_tree() { clear(); }
  };

  static void verify_tree(pga_tree_node pnode, pga_tree_node parent) {
    GMM_ASSERT1(pnode->parent == parent,
                "Invalid tree node " << pnode->node_type);
    for (size_type i = 1; i < pnode->children.size(); ++i)
      verify_tree(pnode->children[i], pnode);
  }

  static void ga_print_constant_tensor(pga_tree_node pnode) {
    size_type nt = pnode->nb_test_functions(); // for printing zero tensors
    switch (pnode->tensor_order()) {
    case 0:
      cout << pnode->t[0];
      break;
    case 1:
      cout << "[";
      for (size_type i = 0; i < pnode->t.size(0); ++i) {
        if (i != 0) cout << "; ";
        cout << (nt ? scalar_type(0) : pnode->t[i]);
      }
      cout << "]";
      break;
    case 2:
      cout << "[";
      for (size_type i = 0; i < pnode->t.size(0); ++i) {
        if (i != 0) cout << "; ";
        for (size_type j = 0; j < pnode->t.size(1); ++j) {
          if (j != 0) cout << ", ";
          cout << (nt ? scalar_type(0) : pnode->t(i,j));
        }
      }
      cout << "]";
      break;
    case 3:
      cout << "[";
      for (size_type i = 0; i < pnode->t.size(0); ++i) {
        if (i != 0) cout << ",, ";
        for (size_type j = 0; j < pnode->t.size(1); ++j) {
          if (j != 0) cout << "; "; 
          for (size_type k = 0; k < pnode->t.size(2); ++k) {
            if (k != 0) cout << ", "; 
            cout << (nt ? scalar_type(0) : pnode->t(i,j,k));
          }
        }
      }
      cout << "]";
      break;
    case 4:
      cout << "[";
      for (size_type i = 0; i < pnode->t.size(0); ++i) {
        if (i != 0) cout << ";; ";
        for (size_type j = 0; j < pnode->t.size(1); ++j) {
          if (j != 0) cout << ",, "; 
          for (size_type k = 0; k < pnode->t.size(2); ++k) {
            if (k != 0) cout << "; "; 
            for (size_type l = 0; l < pnode->t.size(3); ++l) {
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
        bool par = (pnode->parent &&
                    pnode->parent->node_type == GA_NODE_OP &&
                    (ga_operator_priorities[pnode->op_type] >= 2 ||
                     ga_operator_priorities[pnode->op_type]
                     < ga_operator_priorities[pnode->parent->op_type]));
        if (par) cout << "(";
        if (pnode->op_type == GA_UNARY_MINUS) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          cout << "-"; ga_print_node(pnode->children[0]);
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
    case GA_NODE_FIXED_SIZE_TEST: cout << "Test_" << pnode->name; break;
    case GA_NODE_TEST: cout << "Test_" << pnode->name; break;
    case GA_NODE_GRAD_TEST: cout << "Grad_Test_" << pnode->name; break;
    case GA_NODE_HESS_TEST: cout << "Hess_Test_" << pnode->name; break;
    case GA_NODE_PREDEF_FUNC: cout << pnode->name; break;
    case GA_NODE_ZERO:
      GMM_ASSERT1(pnode->test_function_type != size_type(-1),
                  "Internal error");
      if (pnode->test_function_type) cout << "(";
      ga_print_constant_tensor(pnode);
      if (pnode->test_function_type & 1) {
        GMM_ASSERT1(pnode->qdim0 > 0, "Internal error");
        if (pnode->qdim0 == 1)
          cout << "*Test_" << pnode->name_test0;
        else {
          cout << "*([ 0";
          for (size_type i = 1; i < pnode->qdim0; ++i) cout << ", 0";
          cout << "].Test_" << pnode->name_test0 << ")";
        }
      }
      if (pnode->test_function_type & 2) {
        GMM_ASSERT1(pnode->qdim1 > 0, "Internal error");
        if (pnode->qdim1 == 1)
          cout << "*Test2_" << pnode->name_test1;
        else {
          cout << "*([ 0";
          for (size_type i = 1; i < pnode->qdim1; ++i) cout << ", 0";
          cout << "].Test2_" << pnode->name_test1 << ")";
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
    if (tree.root) verify_tree(tree.root, 0);
    if (tree.root) ga_print_node(tree.root); else cout << "Empty tree";
    cout << endl;
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

        case GA_LBRACKET: // Constant vector/matrix or tensor
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
                  ga_throw_error(expr, pos-1, "Bad constant "
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
                ga_throw_error(expr, pos-1, "The constant "
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
        case GA_COLON: case GA_DOT: case GA_TMULT:
          tree.add_op(t_type, token_pos);
          state = 1; break;
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
    size_type pos = 0;
    tree.clear();
    GA_TOKEN_TYPE t = ga_read_term(expr, pos, tree);
    switch (t) {
    case GA_RPAR: ga_throw_error(expr, pos-1, "Unbalanced parenthesis.");
    case GA_RBRACKET: ga_throw_error(expr, pos-1, "Unbalanced braket.");
    case GA_END: break;
    default: ga_throw_error(expr, pos-1, "Unexpected token.");
    }
  }

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

  

  typedef std::map<std::string, ga_predef_function> ga_predef_function_tab;

  static ga_predef_function_tab PREDEF_FUNCTIONS;

  bool init_predef_functions(void) {
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

    return true;
  }


  static bool predef_functions_initialized = init_predef_functions();


  //=========================================================================
  // Structure dealing with variables and constants.
  //=========================================================================
  
  class ga_variables {
    
    const getfem::model *model;


    struct var_description {

      bool is_variable;
      bool is_fem_dofs;
      const mesh_fem *mf;
      const model_real_plain_vector *V;

      var_description(bool is_var,
                      bool is_fem, 
                      const mesh_fem *mmf,
                      const model_real_plain_vector *v)
        : is_variable(is_var), is_fem_dofs(is_fem), mf(mmf), V(v) {}
      var_description() : is_variable(false), is_fem_dofs(false),
                          mf(0), V(0) {}

    };

    typedef std::map<std::string, var_description> VAR_SET;

    VAR_SET variables;

  public:

    void add_fem_variable(const std::string &name, const mesh_fem &mf,
                          const model_real_plain_vector &V) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(true, true, &mf, &V);
    }
    
    void add_fixed_size_variable(const std::string &name,
                                 const model_real_plain_vector &V) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(true, false, 0, &V);
    }

    void add_fem_constant(const std::string &name, const mesh_fem &mf,
                          const model_real_plain_vector &V) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(false, true, &mf, &V);
    }
    
    void add_fixed_size_constant(const std::string &name,
                                 const model_real_plain_vector &V) {
      GMM_ASSERT1(!model, "Invalid use");
      variables[name] = var_description(false, false, 0, &V);
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

    ga_variables(const getfem::model &md) : model(&md) {}
    ga_variables(void) : model(0) {}

  };


  //=========================================================================
  // Semantic analysis, tree simplification and tree enrichment
  //=========================================================================

  static void ga_node_analysis(const std::string &expr, ga_tree &tree,
                               const ga_variables &vars,
                               pga_tree_node pnode) {
    
    bool all_cte = true, all_sc = true, all_primal = true;
    for (size_type i = 0; i < pnode->children.size(); ++i) {
      ga_node_analysis(expr, tree, vars, pnode->children[i]);
      all_cte = all_cte && (pnode->children[i]->node_type == GA_NODE_CONSTANT);
      all_sc = all_sc && pnode->children[i]->t.size() == 1;
      GMM_ASSERT1(pnode->children[i]->test_function_type != size_type(-1),
                  "internal error on child " << i);
      all_primal = all_primal && (pnode->children[i]->test_function_type == 0);
    }
    
    size_type nbch = pnode->children.size();
    pga_tree_node child0 = (nbch > 0) ? pnode->children[0] : 0;
    pga_tree_node child1 = (nbch > 1) ? pnode->children[1] : 0;
    bgeot::multi_index mi;
    const bgeot::multi_index &size0 = child0 ? child0->t.sizes() : mi;
    const bgeot::multi_index &size1 = child1 ? child1->t.sizes() : mi;
    size_type dim0 = child0 ? child0->tensor_order() : 0;
    size_type dim1 = child1 ? child1->tensor_order() : 0;

    switch (pnode->node_type) {
    case GA_NODE_CONSTANT: case GA_NODE_ALLINDICES: case GA_NODE_VAL:
      pnode->test_function_type = 0; break;

    case GA_NODE_GRAD: case GA_NODE_HESS: case GA_NODE_FIXED_SIZE_TEST:
    case GA_NODE_TEST: case GA_NODE_GRAD_TEST: case GA_NODE_HESS_TEST: break;

    case GA_NODE_NAME: // TODO: gerer les noms de fonctions et d'opérateur
      {
        std::string name = pnode->name;

        // Search for a predefined function
        ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
        if (it != PREDEF_FUNCTIONS.end()) {
          pnode->node_type = GA_NODE_PREDEF_FUNC;
          pnode->name = name;
          pnode->test_function_type = 0;
        } else if (!name.compare("pi")) {
          // Predefined constant pi
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->init_scalar_tensor(M_PI);
        } else {
          // Search for a variable name with optional gradient, Hessian 
          // or test functions
          int val_grad_or_hess = 0;
          if (name.size() >= 5 && name.compare(0, 5, "Grad_") == 0)
            { val_grad_or_hess = 1; name = name.substr(5); }
          else if (name.size() >= 5 && name.compare(0, 5, "Hess_") == 0)
            { val_grad_or_hess = 2; name = name.substr(5); }
          int test = 0;
          if (name.size() >= 6 && name.compare(0, 5, "Test1_") == 0)
            { test = 2; name = name.substr(6); }
          else if (name.size() >= 5 && name.compare(0, 5, "Test_") == 0)
            { test = 1; name = name.substr(5); }
          
          if (!(vars.variable_exists(name)))
            ga_throw_error(expr, pnode->pos, "Unknown variable or constant");
          
          const mesh_fem *mf = vars.associated_mf(name);
          if (!test) {
            if (!mf) {
              if (val_grad_or_hess)
                ga_throw_error(expr, pnode->pos, "Gradient or Hessian cannot "
                                "be evaluated for fixed size data.");
              pnode->node_type = GA_NODE_CONSTANT;
              size_type n = gmm::vect_size(vars.value(name));
              if (n == 1) {
                pnode->init_scalar_tensor(vars.value(name)[0]);
              } else {
                pnode->init_vector_tensor(n);
                gmm::copy(vars.value(name), pnode->t.as_vector());
              }
            } else {
              size_type q = vars.qdim(name), n = mf->linked_mesh().dim();
              pnode->name = name;
              switch (val_grad_or_hess) {
              case 0: // value
                // TODO: sortir une instruction d'évaluation de variable
                pnode->node_type = GA_NODE_VAL;
                if (q == 1)
                  pnode->init_scalar_tensor(scalar_type(0));
                else
                  pnode->init_vector_tensor(q);
                break;
              case 1: // grad
                // TODO: sortir une instruction d'évaluation de grad
                pnode->node_type = GA_NODE_GRAD;
                if (q == 1 && n == 1)
                  pnode->init_scalar_tensor(scalar_type(0));
                else if (q == 1)
                  pnode->init_vector_tensor(n);
                else
                  pnode->init_matrix_tensor(q, n);
                break;
              case 2: // Hessian
                // TODO: sortir une instruction d'évaluation de Hessien
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
            if (vars.is_constant(name))
              ga_throw_error(expr, pnode->pos, "Test functions of constant "
                              "are not allowed.");
            pnode->name = name;
            if (test == 1) {
              pnode->name_test0 = name;
              pnode->qdim0
                = (mf ? vars.qdim(name) : gmm::vect_size(vars.value(name)));
            } else {
              pnode->name_test1 = name;
              pnode->qdim1
                = (mf ? vars.qdim(name) : gmm::vect_size(vars.value(name)));
            }
          
            if (!mf) {
              if (val_grad_or_hess)
                ga_throw_error(expr, pnode->pos, "Gradient or Hessian cannot "
                                "be evaluated for fixed size variables.");
              pnode->node_type = GA_NODE_FIXED_SIZE_TEST;
              // TODO: Instruction : mise à l'identité du tenseur (init).
              size_type n = gmm::vect_size(vars.value(name));
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
              size_type q = vars.qdim(name), n = mf->linked_mesh().dim();
              switch (val_grad_or_hess) {
              case 0: // value
                // TODO: sortir une instruction type Base_value
                // et un tenseur intermédiaire
                // TODO gérer la taille locale du fem (inconnu) (mise à 1)
                pnode->node_type = GA_NODE_TEST;
                if (q == 1)
                  pnode->init_vector_tensor(1);
                else
                  pnode->init_matrix_tensor(1,q);
                pnode->test_function_type = test;
                break;
              case 1: // grad
                // TODO: sortir une instruction type Grad_Base_value
                // et un tenseur intermédiaire
                // TODO gérer la taille locale du fem (inconnu) (mise à 1)
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
                // TODO: sortir une instruction type Hess_Base_value
                // et un tenseur intermédiaire
                // TODO gérer la taille locale du fem (inconnu) (mise à 1)
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
            if (child0->test_function_type != child1->test_function_type ||
                child0->name_test0.compare(child1->name_test0) ||
                child0->name_test1.compare(child1->name_test1))
              compatible = false;
          }
          // TODO: prendre en compte l'addition de différentes parties
          // de formulation faible avec différentes fonctions test.
          // Dans ce cas, c'est l'assemblage qui gère l'addition. 
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
            pnode->t = child1->t;
            pnode->test_function_type = child1->test_function_type;
            pnode->name_test0 = child1->name_test0;
            pnode->name_test1 = child1->name_test1;
            pnode->qdim0 = child1->qdim0;
            pnode->qdim1 = child1->qdim1;

            // simplification if one of the two operand is constant and zero
            if (child0->tensor_is_zero())
              tree.replace_node_by_child(pnode, 1);
            else if (child1->tensor_is_zero())
              tree.replace_node_by_child(pnode, 0);
          }
        }
        break;

      case GA_UNARY_MINUS:
        pnode->t = child0->t;
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test0 = child0->name_test0;
        pnode->name_test1 = child0->name_test1;
        pnode->qdim0 = child0->qdim0;
        pnode->qdim1 = child0->qdim1;
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->test_function_type = 0;
          gmm::scale(pnode->t.as_vector(), scalar_type(-1));
          tree.clear_children(pnode);
        } else if (child0->node_type == GA_NODE_ZERO) {
          tree.replace_node_by_child(pnode, 0);
        }
        break;

      case GA_DOT:
        if (dim0 > 1 || dim1 > 1)
          ga_throw_error(expr, pnode->pos,
                          "Dot product acts only on vectors.");
        else {
          size_type s0 = dim0 == 0 ? 1 : size0.back();
          size_type s1 = dim1 == 0 ? 1 : size1.back();
          if (s0 != s1)
            ga_throw_error(expr, pnode->pos, "Dot product "
                            "of expressions of different sizes");
          
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->init_scalar_tensor
              (gmm::vect_sp(child0->t.as_vector(), child1->t.as_vector()));
            tree.clear_children(pnode);
          } else {
            pnode->mult_test(child0, child1, expr);
            if (child0->tensor_is_zero() || child1->tensor_is_zero()) {
              gmm::clear(pnode->t.as_vector());
              pnode->node_type = GA_NODE_ZERO;
              tree.clear_children(pnode);
            }
          }
        }
        break;

      case GA_COLON:
        if (dim0 > 2 || dim1 > 2)
          ga_throw_error(expr, pnode->pos,
                          "Frobenius product acts only on matrices.");
        else {
          size_type s00 = (dim0 == 0) ? 1
            : (dim0 == 1 ? size0.back() : size0[size0.size()-2]);
          size_type s01 = (dim0 == 2) ? size0.back() : 1;
          size_type s10 = (dim1 == 0) ? 1
            : (dim1 == 1 ? size1.back() : size1[size1.size()-2]);
          size_type s11 = (dim1 == 2) ? size1.back() : 1;
          if (s00 != s10 || s01 != s11)
            ga_throw_error(expr, pnode->pos, "Frobenius product "
                            "of expressions of different sizes");

          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->init_scalar_tensor
              (gmm::vect_sp(pnode->children[0]->t.as_vector(),
                            pnode->children[1]->t.as_vector()));
            tree.clear_children(pnode);
          } else {
            pnode->mult_test(child0, child1, expr);
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
          pnode->mult_test(child0, child1, expr);
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
        pnode->name_test0 = child0->name_test0;
        pnode->name_test1 = child0->name_test1;
        pnode->qdim0 = child0->qdim0;
        pnode->qdim1 = child0->qdim1;
        
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
        if (!all_primal)
          ga_throw_error(expr, pnode->pos, "Test functions are not allowed "
                          "in constant vector/matrix/tensor components.");
        
        size_type nbc1 = pnode->nbc1, nbc2 = pnode->nbc2, nbc3 = pnode->nbc3;
        size_type nbl = pnode->children.size() / (nbc1*nbc2*nbc3);
        if (all_cte) pnode->node_type = GA_NODE_CONSTANT;
        pnode->test_function_type = 0;
        if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1 && nbl == 1) {
          pnode->init_scalar_tensor(child0->t[0]);
        } else if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
          pnode->init_vector_tensor(nbl);
          if (all_cte)
            for (size_type i = 0; i < nbl; ++i)
              pnode->t[i] = pnode->children[i]->t[0];
        } else if (nbc2 == 1 && nbc3 == 1) {
          pnode->init_matrix_tensor(nbl, nbc1);
          if (all_cte)
            for (size_type i = 0; i < nbl; ++i)
              for (size_type j = 0; j < nbc1; ++j)
                pnode->t(i,j) = pnode->children[i*nbc1+j]->t[0];
        } else {
          pnode->init_fourth_order_tensor(nbl, nbc3, nbc2, nbc1);
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

    case GA_NODE_PARAMS:

      if (child0->node_type == GA_NODE_PREDEF_FUNC) {
        // Évaluation of predefined function
        std::string name = child0->name;
        ga_predef_function_tab::iterator it = PREDEF_FUNCTIONS.find(name);
        GMM_ASSERT1(it != PREDEF_FUNCTIONS.end(), "internal error");
        size_type nbargs = it->second.nbargs;
        GMM_ASSERT1(it->second.ftype != 2, "Not already taken into account");
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
        size_type s1 = child1->t.size();
        size_type s2 = (nbargs == 2) ? child2->t.size() : s1;
        if (s1 != s2 && (s1 != 1 || s2 != 1))
          ga_throw_error(expr, pnode->pos,
                         "Invalid argument size for a scalar function.");

        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          if (nbargs == 1) {
            pnode->t = child1->t;
            for (size_type i = 0; i < s1; ++i)
              pnode->t[i] = (*(it->second.f1))(child1->t[i]);
          } else {
            if (s1 == s2) {
              pnode->t = child1->t;
              for (size_type i = 0; i < s1; ++i)
                pnode->t[i] = (*(it->second.f2))(child1->t[i], child2->t[i]);
            } else if (s1 == 1) {
              pnode->t = child2->t;
              for (size_type i = 0; i < s2; ++i)
                pnode->t[i] = (*(it->second.f2))(child1->t[0], child2->t[i]);
            } else {
              pnode->t = child1->t;
              for (size_type i = 0; i < s1; ++i)
                pnode->t[i] = (*(it->second.f2))(child1->t[i], child2->t[0]);
            }
          }
          tree.clear_children(pnode);
        }
      } else { // Access to components of a tensor
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
            mi1[i] = size_type(::round(pnode->children[i+1]->t[0])-1);
            if (mi1[i] >= child0->tensor_proper_size(i))
              ga_throw_error(expr, pnode->children[i+1]->pos,
                              "Index out of range.");
          }
        }
        mi = size0; mi.resize(child0->nb_test_functions());
        for (size_type i = 0; i < mi2.size(); ++i) mi.push_back(mi2[i]);
        pnode->t.adjust_sizes(mi);
        pnode->test_function_type = child0->test_function_type;
        pnode->name_test0 = child0->name_test0;
        pnode->name_test1 = child0->name_test1;
        pnode->qdim0 = child0->qdim0;
        pnode->qdim1 = child0->qdim1;
        
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

    default:GMM_ASSERT1(false, "Unexpected node type. Internal error.");
    }
    
  }

} // end of namespace getfem


//=========================================================================
// Small test for debug
//=========================================================================

#include "getfem/getfem_regular_meshes.h"

namespace getfem {

  void lex_analysis(void) {

    // std::string expr="([1,2;3,4]@[1,2;1,2])(:,2,1,1)(1)+ [1,2;3,4](1,:)(2)"; // should give 4
    // std::string expr="[1,2;3,4]@[1,2;1,2]*[2,3;2,1]/4 + [1,2;3,1]*[1;1](1)"; // should give [4, 8; 12, 13]
    // std::string expr="[1,2;3,a](2,:) + b(:)"; // should give [6, 9]
    // std::string expr="[1,1;1,2,,1,1;1,2;;1,1;1,2,,1,1;1,3](:,:,:,2)";
    std::string expr="sin([pi;2*pi])";
    // std::string expr = "p*Trace(Grad_Test_u) + Test_u(1,2)+1.0E-1";
    // std::string expr = "(3*(1*Grad_u)).Grad_Test_u*2 + 0*[1;2].Grad_Test_u + c*Grad_Test_u(1) + [u;1](1)*Test_u";
    // std::string expr = "-(4+(2*3)+2*(1+2))/-(-3+5)"; // should give 8
    // std::string expr="[1,2;3,4]@[1,2;1,2]*(Grad_u@Grad_u)/4 + [1,2;3,1]*[1;1](1)";
    
    ga_variables vars; 

    model_real_plain_vector a(1); a[0] = 3.0;
    vars.add_fixed_size_constant("a", a);
    model_real_plain_vector b(2); b[0] = 3.0; b[1] = 6.0;
    vars.add_fixed_size_constant("b", b);
    model_real_plain_vector c(1); c[0] = 1.0;
    vars.add_fixed_size_constant("c", c);
    
    getfem::mesh m;

    bgeot::pgeometric_trans pgt = 
      bgeot::geometric_trans_descriptor("GT_PK(2,1)");
    size_type N = pgt->dim();
    std::vector<size_type> nsubdiv(N);
    std::fill(nsubdiv.begin(),nsubdiv.end(), 10);
    getfem::regular_unit_mesh(m, nsubdiv, pgt);

    getfem::mesh_fem mf(m);
    getfem::pfem pf_u = getfem::fem_descriptor("FEM_PK(2,1)");
    mf.set_finite_element(m.convex_index(), pf_u);

    std::vector<scalar_type> U(mf.nb_dof());

    vars.add_fem_variable("u", mf, U);


    ga_tree tree;
    ga_read_string(expr, tree);
    ga_print_tree(tree);
    ga_node_analysis(expr, tree, vars, tree.root);
    ga_print_tree(tree);
  }
  

} /* end of namespace */
