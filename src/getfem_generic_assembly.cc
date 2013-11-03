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

  static bool initialized = init_ga_char_type();

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

  enum GA_NODE_TYPE {
    GA_NODE_VOID = 0,
    GA_NODE_OP,
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
    GA_NODE_HESS_TEST};

  struct ga_tree_node {
    GA_NODE_TYPE node_type;
    base_tensor t;
    size_type test_function_type; // 0 = no_type function, 1 = first order
                                  // 2 = second order, 3 = both
    size_type nbc1, nbc2, nbc3, pos;
    std::string name;
    GA_TOKEN_TYPE op_type;
    bool valid;
    ga_tree_node *parent; // only one for the moment ...
    std::vector<ga_tree_node *> children;
    
    inline size_type tensor_order(void) const {
      return t.sizes().size() - test_function_type
        + (test_function_type >= 2 ? 1 : 0);
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
    void init_order_three_tensor(size_type n, size_type m,  size_type l) {
      t.adjust_sizes(bgeot::multi_index(n,m,l));
      test_function_type = 0;
    }
    void init_order_four_tensor(size_type n, size_type m,
                                size_type l, size_type k) {
      t.adjust_sizes(bgeot::multi_index(n,m,l,k));
      test_function_type = 0;
    }

    ga_tree_node(void): node_type(GA_NODE_VOID), test_function_type(0),
                        valid(true)  {}
    ga_tree_node(GA_NODE_TYPE ty, size_type p)
      : node_type(ty), test_function_type(0), pos(p), valid(true) {}
    ga_tree_node(scalar_type v, size_type p)
      : node_type(GA_NODE_CONSTANT), test_function_type(0), pos(p), valid(true)
    { init_scalar_tensor(v); }
    ga_tree_node(const char *n, size_type l, size_type p)
      : node_type(GA_NODE_NAME), test_function_type(0), pos(p),
        name(n, l), valid(true) {}
    ga_tree_node(GA_TOKEN_TYPE op, size_type p)
      : node_type(GA_NODE_OP), test_function_type(0), pos(p),
        op_type(op), valid(true) {}
    
  };

  typedef ga_tree_node *pga_tree_node;

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
      if (pnode->valid) {
        pnode->valid = false;
        for (size_type i = 0; i < pnode->children.size(); ++i)
          clear_node(pnode->children[i]);
        delete pnode;
        current_node = 0;
      }
    }

    void clear(void)
    { if (root) clear_node(root); root = current_node = 0; }
    
    void clear_children(pga_tree_node pnode) {
      for (size_type i = 0; i < pnode->children.size(); ++i)
        clear_node(pnode->children[i]);
      pnode->children.resize(0);
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

  static void ga_print_node(pga_tree_node pnode) {
    
    switch(pnode->node_type) {
    case GA_NODE_OP: 
      {
        bool par = (pnode->parent &&
                    pnode->parent->node_type == GA_NODE_OP &&
                    ga_operator_priorities[pnode->op_type]
                    < ga_operator_priorities[pnode->parent->op_type]);
        if (par) cout << "(";
        if (pnode->op_type == GA_UNARY_MINUS) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          cout << "-"; ga_print_node(pnode->children[0]);
        } else {
          GMM_ASSERT1(pnode->children.size() == 2, "Invalid tree");
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

    case GA_NODE_CONSTANT:
      switch (pnode->tensor_order()) {
      case 0:
        cout << pnode->t[0];
        break;
      case 1:
        cout << "[[";
        for (size_type i = 0; i < pnode->t.size(0); ++i)
          { if (i != 0) cout << "; "; cout << pnode->t[i]; }
        cout << "]]";
        break;
      case 2:
        cout << "[[";
        for (size_type i = 0; i < pnode->t.size(0); ++i) {
          if (i != 0) cout << "; ";
          for (size_type j = 0; j < pnode->t.size(1); ++j)
            { if (j != 0) cout << ", "; cout << pnode->t(i,j); }
        }
        cout << "]]";
        break;
      case 4:
        cout << "[[";
        for (size_type i = 0; i < pnode->t.size(0); ++i) {
          if (i != 0) cout << ";; ";
          for (size_type j = 0; j < pnode->t.size(1); ++j) {
            if (j != 0) cout << ",, "; 
            for (size_type k = 0; k < pnode->t.size(2); ++k) {
              if (k != 0) cout << "; "; 
              for (size_type l = 0; l < pnode->t.size(3); ++l)
                { if (l != 0) cout << ", "; cout << pnode->t(i,j,k,l); }
            }
          }
        }
        cout << "]]";
        break;
      default: GMM_ASSERT1(false, "Invalid tensor dimension");
      }
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

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

    default: cout << "Invalid or not taken into account node type"; break;
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
  
  
  static void ga_syntax_error(const std::string &expr, size_type pos,
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
             
  // Read a term with a pushdown automaton.
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
              ga_syntax_error(expr, token_pos, "Bad numeric format.");
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
              ga_syntax_error(expr, pos-1, "Unbalanced parenthesis.");
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
                  ga_syntax_error(expr, pos-1, "Bad constant "
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
                ga_syntax_error(expr, pos-1, "The constant "
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

        default: ga_syntax_error(expr, token_pos, "Unexpected token.");
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
                ga_syntax_error(expr, pos-((r_type != GA_END)?1:0),
                               "Parameters should be separated "
                               "by ',' and parameter list ended by ')'.");
              tree.add_sub_tree(sub_tree);
              if (r_type == GA_RPAR) break;
            }
            state = 2;
          }
          break;
          
        default: ga_syntax_error(expr, token_pos, "Unexpected token.");
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
    case GA_RPAR: ga_syntax_error(expr, pos-1, "Unbalanced parenthesis.");
    case GA_RBRACKET: ga_syntax_error(expr, pos-1, "Unbalanced braket.");
    case GA_END: break;
    default: ga_syntax_error(expr, pos-1, "Unexpected token.");
    }
  }

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

  // TODO: should also detect and eliminate void operations such as
  //       0+expr, 1*expr, 1@expr, I*expr, 0*expr (all elinitate in that case)
    static void ga_node_analysis(const std::string &expr, ga_tree &tree,
                                 const ga_variables &vars,
                                 pga_tree_node pnode) {

    bool all_cte = true, all_sc = true;
    for (size_type i = 0; i < pnode->children.size(); ++i) {
      ga_node_analysis(expr, tree, vars, pnode->children[i]);
      all_cte = all_cte && (pnode->children[i]->node_type == GA_NODE_CONSTANT);
      all_sc = all_sc && pnode->children[i]->t.size() == 1;
    }
    
    switch (pnode->node_type) {
    case GA_NODE_CONSTANT: case GA_NODE_ALLINDICES: case GA_NODE_VAL:
    case GA_NODE_GRAD: case GA_NODE_HESS: case GA_NODE_FIXED_SIZE_TEST:
    case GA_NODE_TEST: case GA_NODE_GRAD_TEST: case GA_NODE_HESS_TEST: break;

    case GA_NODE_NAME: // TODO: gerer les noms de fonctions
      {
        std::string name = pnode->name;
        int val_grad_or_hess = 0;
        if (name.size() >= 5 && name.compare(0, 5, "Grad_") == 0)
          { val_grad_or_hess = 1; name = name.substr(5); }
        else if (name.size() >= 5 && name.compare(0, 5, "Hess_") == 0)
          { val_grad_or_hess = 2; name = name.substr(5); }
        int test = 0;
        if (name.size() >= 6 && name.compare(0, 5, "Test2_") == 0)
          { test = 2; name = name.substr(6); }
        else if (name.size() >= 5 && name.compare(0, 5, "Test_") == 0)
          { test = 1; name = name.substr(5); }
        
        if (!(vars.variable_exists(name)))
          ga_syntax_error(expr, pnode->pos, "Unknown variable or constant");

        const mesh_fem *mf = vars.associated_mf(name);
        if (!test) {
          if (!mf) {
            if (val_grad_or_hess)
              ga_syntax_error(expr, pnode->pos, "Gradient or Hessian cannot "
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
            case 2: // hessian
              // TODO: sortir une instruction d'évaluation de Hessien
              pnode->node_type = GA_NODE_HESS;
              if (q == 1 && n == 1)
                pnode->init_scalar_tensor(scalar_type(0));
              else if (q == 1)
                pnode->init_vector_tensor(n*n);
              else
                pnode->init_matrix_tensor(q, n*n);
              break;
            }
          }
        } else {
          if (vars.is_constant(name))
            ga_syntax_error(expr, pnode->pos, "Test functions of constant are "
                            "not allowed.");
          pnode->name = name;
          if (!mf) {
            if (val_grad_or_hess)
              ga_syntax_error(expr, pnode->pos, "Gradient or hessian cannot "
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
                pnode->init_order_three_tensor(1,q,n);
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
                pnode->init_matrix_tensor(1,n*n);
              else
                pnode->init_order_three_tensor(1,q,n*n);
              pnode->test_function_type = test;
              break;
            }
          }
        }
      }
      break;

    case GA_NODE_OP:
      switch(pnode->op_type) {
      case GA_PLUS: case GA_MINUS:
        {
          const bgeot::multi_index &size0 = pnode->children[0]->t.sizes();
          const bgeot::multi_index &size1 = pnode->children[1]->t.sizes();
          size_type c_size = std::min(size0.size(), size1.size());
          bool compatible = true;
          for (size_type i = 0; i < c_size; ++i)
            if (size0[i] != size1[i]) compatible = false; 
          for (size_type i = c_size; i < size0.size(); ++i)
            if (size0[i] != 1) compatible = false;
          for (size_type i = c_size; i < size1.size(); ++i)
            if (size1[i] != 1) compatible = false;
          if (!compatible)
            ga_syntax_error(expr, pnode->pos, "Addition or substraction "
                            "of expressions of different sizes");
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->t = pnode->children[0]->t;
            if (pnode->op_type == GA_MINUS)
              pnode->t -= pnode->children[1]->t;
            else
              pnode->t += pnode->children[1]->t;
            tree.clear_children(pnode);
          }
        }
        break;

      case GA_UNARY_MINUS:
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->t = pnode->children[0]->t;
          gmm::scale(pnode->t.as_vector(), scalar_type(-1));
          tree.clear_children(pnode);
        }
        break;

      case GA_DOT:
        {
          size_type dim0 = pnode->children[0]->tensor_order();
          size_type dim1 = pnode->children[1]->tensor_order();
          if (dim0 > 1 || dim1 > 1)
            ga_syntax_error(expr, pnode->pos,
                            "Dot product acts only on vectors.");
          const bgeot::multi_index &size0 = pnode->children[0]->t.sizes();
          const bgeot::multi_index &size1 = pnode->children[1]->t.sizes();
          size_type c_size = std::min(dim0, dim1);
          bool compatible = true;
          for (size_type i = 0; i < c_size; ++i)
            if (size0[i] != size1[i]) compatible = false; 
          for (size_type i = c_size; i < dim0; ++i)
            if (size0[i] != 1) compatible = false;
          for (size_type i = c_size; i < dim1; ++i)
            if (size1[i] != 1) compatible = false;
          if (!compatible)
            ga_syntax_error(expr, pnode->pos, "Dot product "
                            "of expressions of different sizes");

          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->init_scalar_tensor
              (gmm::vect_sp(pnode->children[0]->t.as_vector(),
                            pnode->children[1]->t.as_vector()));
            tree.clear_children(pnode);
          }
        }
        break;

      case GA_COLON:
        {
          size_type dim0 = pnode->children[0]->tensor_order();
          size_type dim1 = pnode->children[1]->tensor_order();
          if (dim0 > 2 || dim1 > 2)
            ga_syntax_error(expr, pnode->pos,
                            "Frobenius product acts only on matrices.");
          const bgeot::multi_index &size0 = pnode->children[0]->t.sizes();
          const bgeot::multi_index &size1 = pnode->children[1]->t.sizes();
          size_type c_size = std::min(dim0, dim1);
          bool compatible = true;
          for (size_type i = 0; i < c_size; ++i)
            if (size0[i] != size1[i]) compatible = false; 
          for (size_type i = c_size; i < dim0; ++i)
            if (size0[i] != 1) compatible = false;
          for (size_type i = c_size; i < dim1; ++i)
            if (size1[i] != 1) compatible = false;
          if (!compatible)
            ga_syntax_error(expr, pnode->pos, "Frobenius product "
                            "of expressions of different sizes");
          
          if (all_cte) {
            pnode->node_type = GA_NODE_CONSTANT;
            pnode->init_scalar_tensor
              (gmm::vect_sp(pnode->children[0]->t.as_vector(),
                            pnode->children[1]->t.as_vector()));
            tree.clear_children(pnode);
          }

        }
        break;

      case GA_TMULT:
        if (all_cte) {
          size_type dim0 = pnode->children[0]->tensor_order();
          size_type dim1 = pnode->children[1]->tensor_order();
          pnode->node_type = GA_NODE_CONSTANT;
          if (pnode->children[0]->t.size() == 1 &&
              pnode->children[0]->t.size() == 0) {
            pnode->init_scalar_tensor
              (pnode->children[0]->t[0] * pnode->children[1]->t[0]);
          } else if (pnode->children[0]->t.size() == 1) {
            pnode->t = pnode->children[1]->t;
            gmm::scale(pnode->t.as_vector(),
                       scalar_type(pnode->children[0]->t[0]));
          } else if (pnode->children[1]->t.size() == 1) {
            pnode->t = pnode->children[0]->t;
            gmm::scale(pnode->t.as_vector(),
                       scalar_type(pnode->children[1]->t[0]));
          } else {
            if (dim0+dim1 == 3 || dim0+dim1 > 4)
              ga_syntax_error(expr, pnode->pos, "Unauthorized "
                              "tensor multiplication.");
            bgeot::multi_index mi(dim0 + dim1);
            for (size_type i = 0; i < dim0; ++i)
              mi[i] = pnode->children[0]->t.size(i);
            for (size_type i = 0; i < dim1; ++i)
              mi[i+dim0] = pnode->children[1]->t.size(i);
            pnode->t.adjust_sizes(mi);
            size_type n0 = pnode->children[0]->t.size();
            size_type n1 = pnode->children[1]->t.size();
            for (size_type i = 0; i < n0; ++i)
              for (size_type j = 0; j < n1; ++j)
                pnode->t[i+j*n0] =
                  pnode->children[0]->t[i] * pnode->children[1]->t[j];
          }
          tree.clear_children(pnode);
        }
        break;

      case GA_MULT:
        if (all_cte) {
          pnode->node_type = GA_NODE_CONSTANT;
          if (pnode->children[0]->t.size() == 1 &&
              pnode->children[0]->t.size() == 0) {
            pnode->init_scalar_tensor
              (pnode->children[0]->t[0] * pnode->children[1]->t[0]);
          }           if (pnode->children[0]->t.size() == 1) {
            pnode->t = pnode->children[1]->t;
            gmm::scale(pnode->t.as_vector(), pnode->children[0]->t[0]);
          } else if (pnode->children[1]->t.size() == 1) {
            pnode->t = pnode->children[0]->t;
            gmm::scale(pnode->t.as_vector(), pnode->children[1]->t[0]);
          } else if (pnode->children[0]->tensor_order() == 2 &&
                     pnode->children[1]->tensor_order() == 1) {
            size_type m = pnode->children[0]->t.size(0);
            size_type n = pnode->children[0]->t.size(1);
            if (n != pnode->children[1]->t.size(0))
              ga_syntax_error(expr, pnode->pos, "Incompatible sizes in "
                              "matrix-vector multiplication.");
            pnode->init_vector_tensor(m);
            gmm::clear(pnode->t.as_vector());
            for (size_type i = 0; i < m; ++i)
              for (size_type j = 0; j < n; ++j)
                pnode->t[i] += pnode->children[0]->t(i,j)
                  * pnode->children[1]->t[j];
          } else if (pnode->children[0]->tensor_order() == 2 &&
                     pnode->children[1]->tensor_order() == 2) {
            size_type m = pnode->children[0]->t.size(0);
            size_type n = pnode->children[0]->t.size(1);
            size_type p = pnode->children[1]->t.size(1);
            if (n != pnode->children[1]->t.size(0))
              ga_syntax_error(expr, pnode->pos, "Incompatible sizes in "
                              "matrix-matrix multiplication.");
            pnode->init_matrix_tensor(m,p);
            gmm::clear(pnode->t.as_vector());
            for (size_type i = 0; i < m; ++i)
              for (size_type j = 0; j < n; ++j)
                for (size_type k = 0; k < p; ++k)
                  pnode->t(i,k) += pnode->children[0]->t(i,j)
                    * pnode->children[1]->t(j,k);
          }
          else if (pnode->children[0]->tensor_order() == 4 &&
                   pnode->children[1]->tensor_order() == 2) {
            size_type m = pnode->children[0]->t.size(0);
            size_type n = pnode->children[0]->t.size(1);
            size_type o = pnode->children[0]->t.size(2);
            size_type p = pnode->children[0]->t.size(3);
            if (o != pnode->children[1]->t.size(0) ||
                p != pnode->children[1]->t.size(1))
              ga_syntax_error(expr, pnode->pos, "Incompatible sizes in "
                              "tensor-matrix multiplication.");
            pnode->init_matrix_tensor(m,n);
            gmm::clear(pnode->t.as_vector());
            for (size_type i = 0; i < m; ++i)
              for (size_type j = 0; j < n; ++j)
                for (size_type k = 0; k < o; ++k)
                  for (size_type l = 0; l < p; ++l)
                    pnode->t(i,j) += pnode->children[0]->t(i,j,k,l)
                      * pnode->children[1]->t(k,l);
          } else ga_syntax_error(expr, pnode->pos,
                                 "Unauthorized multiplication.");
          tree.clear_children(pnode);
        }
        break;

      case GA_DIV:
        if (pnode->children[1]->tensor_order() > 0)
          ga_syntax_error(expr, pnode->pos, "Only the division by a scalar "
                          "is allowed.");
        if (pnode->children[1]->node_type == GA_NODE_CONSTANT &&
            pnode->children[1]->t[0] == scalar_type(0))
          ga_syntax_error(expr, pnode->children[1]->pos, "Division by zero");
        
        if (all_cte) {          
          pnode->node_type = GA_NODE_CONSTANT;
          pnode->t = pnode->children[0]->t;
          gmm::scale(pnode->t.as_vector(),
                     scalar_type(1) / pnode->children[1]->t[0]);
          tree.clear_children(pnode);
        }
        break;

      default:GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX:
      if (!all_sc) 
        ga_syntax_error(expr, pnode->pos, "Constant vector/matrix/tensor "
                        "components should be scalar valued.");
      if (all_cte) {
        size_type nbc1 = pnode->nbc1, nbc2 = pnode->nbc2, nbc3 = pnode->nbc3;
        size_type nbl = pnode->children.size() / (nbc1*nbc2*nbc3);
        pnode->node_type = GA_NODE_CONSTANT;
        if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1 && nbl == 1) {
          pnode->init_scalar_tensor(pnode->children[0]->t[0]);
        } else if (nbc1 == 1 && nbc2 == 1 && nbc3 == 1) {
          pnode->init_vector_tensor(nbl);
          for (size_type i = 0; i < nbl; ++i)
            pnode->t[i] = pnode->children[i]->t[0];
        } else if (nbc2 == 1 && nbc3 == 1) {
          pnode->init_matrix_tensor(nbl, nbc1);
          for (size_type i = 0; i < nbl; ++i)
            for (size_type j = 0; j < nbc1; ++j)
              pnode->t(i,j) = pnode->children[i*nbc1+j]->t[0];
        } else {
          pnode->init_order_four_tensor(nbl, nbc3, nbc2, nbc1);
          size_type n = 0;
          for (size_type i = 0; i < nbl; ++i)
            for (size_type j = 0; j < nbc3; ++j)
              for (size_type k = 0; k < nbc2; ++k)
                for (size_type l = 0; l < nbc1; ++l)
                  pnode->t(i,j,k,l) = pnode->children[n++]->t[0];
        }
        tree.clear_children(pnode);
      }
      break;

    case GA_NODE_PARAMS:
      // Should treat the params of a function and the index on
      // non constant terms
      
      if (pnode->children[0]->node_type == GA_NODE_CONSTANT) {
        if (pnode->children.size() != pnode->children[0]->tensor_order() + 1)
          ga_syntax_error(expr, pnode->pos, "Bad number of indices.");
        for (size_type i = 1; i < pnode->children.size(); ++i)
          if (pnode->children[i]->node_type != GA_NODE_ALLINDICES &&
              (pnode->children[i]->node_type != GA_NODE_CONSTANT ||
               pnode->children[i]->t.size() != 1))
            ga_syntax_error(expr, pnode->children[i]->pos,
                            "Indices should be constant integers or colon.");
        
        const bgeot::multi_index &mi0 = pnode->children[0]->t.sizes();
        bgeot::multi_index mi1(mi0.size()), mi2, indices;
        for (size_type i = 0; i < mi0.size(); ++i) {
          if (pnode->children[i+1]->node_type == GA_NODE_ALLINDICES)
            { mi2.push_back(mi0[i]); indices.push_back(i); mi1[i] = 0; }
          else {
            mi1[i] = size_type(::round(pnode->children[i+1]->t[0])-1);
            if (mi1[i] >= mi0[i])
              ga_syntax_error(expr, pnode->children[i+1]->pos,
                              "Index out of range.");
          }
        }
        if (mi2.size() == 3)
          ga_syntax_error(expr, pnode->pos,
                          "Sorry, order three tensors are not allowed.");
        
        pnode->node_type = GA_NODE_CONSTANT;
        pnode->t.adjust_sizes(mi2);

        for (bgeot::multi_index mi3(mi2.size()); !mi3.finished(mi2);
             mi3.incrementation(mi2)) {
          for (size_type j = 0; j < mi2.size(); ++j) {
            mi1[indices[j]] = mi3[j];
          }
          pnode->t(mi3) = pnode->children[0]->t(mi1);
        }
        tree.clear_children(pnode);
      }
      break;

    default:GMM_ASSERT1(false, "Unexpected node type. Internal error.");
    }
    
  }

  //=========================================================================
  // Small test for debug
  //=========================================================================

}

#include "getfem/getfem_regular_meshes.h"

namespace getfem {

  void lex_analysis(void) {

    // std::string expr="([1,2;3,4]@[1,2;1,2])(:,2,1,1)(1)+ [1,2;3,4](1,:)(2)"; // should give 4
    // std::string expr="[1,2;3,4]@[1,2;1,2]*[2,3;2,1]/4 + [1,2;3,1]*[1;1](1)"; // should give [[4, 8; 12, 13]]
    // std::string expr="[1,2;3,a](2,:) + b(:)"; // should give [[6, 9]]
    // std::string expr="[1,1;1,2,,1,1;1,2;;1,1;1,2,,1,1;1,3]";
    // std::string expr = "p*Trace(Grad_Test_u) + Test_u(1,2)+1.0E-1";
    std::string expr = "Grad_u.Grad_Test_u";
    // std::string expr = "-(4+(2*3)+2*(1+2))/-(-3+5)"; // should give 8
    
    ga_variables vars; 

    model_real_plain_vector a(1); a[0] = 3.0;
    vars.add_fixed_size_constant("a", a);
    model_real_plain_vector b(2); b[0] = 3.0; b[1] = 6.0;
    vars.add_fixed_size_constant("b", b);
    
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
