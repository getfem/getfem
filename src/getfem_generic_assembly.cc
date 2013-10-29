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

#include "getfem/getfem_config.h"
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
    GA_INVALID = 0,// invalid token
    GA_END,        // string end
    GA_NAME,       // A variable or user defined nonlinear function name
    GA_SCALAR,     // A real number
    GA_PLUS,       // '+'
    GA_MINUS,      // '-'
    GA_UNARY_MINUS,// '-'
    GA_MULT,       // '*'
    GA_DIV,        // '/'
    GA_COLON,      // ':'
    GA_DOT,        // '.'
    GA_TMULT,      // '@' tensor product
    GA_COMMA,      // ','
    GA_SEMICOLON,  // ';'
    GA_LPAR,       // '('
    GA_RPAR,       // ')'
    GA_LBRACKET,   // '['
    GA_RBRACKET,   // ']'
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
    default: return type;
    }
  }
  

  //=========================================================================
  // Tree structure for syntax analysis
  //=========================================================================

  enum GA_NODE_TYPE {
    GA_NODE_VOID = 0,
    GA_NODE_OP,
    GA_NODE_SCALAR,
    GA_NODE_VECTOR,
    GA_NODE_MATRIX,
    GA_NODE_TENSOR,
    GA_NODE_NAME,
    GA_NODE_PARAMS,
    GA_NODE_ALLINDICES,
    GA_NODE_C_MATRIX};

  struct ga_tree_node {
    GA_NODE_TYPE node_type;
    scalar_type val;
    base_vector vec;
    base_matrix mat;
    size_type nbc, pos;
    std::string name;
    GA_TOKEN_TYPE op_type;
    bool valid;
    ga_tree_node *parent; // only one for the moment ...
    std::vector<ga_tree_node *> children;

    ga_tree_node(void): node_type(GA_NODE_VOID), valid(true) {}
    ga_tree_node(GA_NODE_TYPE t, size_type p)
      : node_type(t), pos(p), valid(true) {}
    ga_tree_node(scalar_type v, size_type p)
      : node_type(GA_NODE_SCALAR), val(v), pos(p), valid(true) {}
    ga_tree_node(const char *n, size_type l, size_type p)
      : node_type(GA_NODE_NAME), pos(p), name(n, l), valid(true) {}
    ga_tree_node(GA_TOKEN_TYPE op, size_type p)
      : node_type(GA_NODE_OP), pos(p), op_type(op), valid(true) {}
    
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
      current_node->nbc = 0;
    }
    
    void add_op(GA_TOKEN_TYPE op_type, size_type pos) {
      while (current_node &&
             current_node->parent &&
             current_node->parent->node_type == GA_NODE_OP &&
             ga_operator_priorities[current_node->parent->op_type]
             >= ga_operator_priorities[op_type])
        current_node = current_node->parent;
      pga_tree_node new_node = new ga_tree_node(op_type, pos);
      if (current_node) {
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
      cout << "(";
      switch (pnode->op_type) {
      case GA_PLUS: cout << "+"; break;
      case GA_MINUS: case GA_UNARY_MINUS: cout << "-"; break;
      case GA_MULT: cout << "*"; break;
      case GA_DIV: cout << "/"; break;
      case GA_COLON: cout << ":"; break;
      case GA_DOT: cout << "."; break;
      case GA_TMULT: cout << "@"; break;
      default: cout << "Invalid or not taken into account operation"; break;
      }
      for (size_type i = 0; i < pnode->children.size(); ++i)
        { cout << " "; ga_print_node(pnode->children[i]); }
      cout << ")";
      break;

    case GA_NODE_SCALAR:
      cout << pnode->val;
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_VECTOR:
      cout << "[[";
      for (size_type i = 0; i < gmm::vect_size(pnode->vec); ++i)
        { if (i != 0) cout << ", "; cout << pnode->vec[i]; }
      cout << "]]";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_MATRIX:
      cout << "[[";
      for (size_type i = 0; i < gmm::mat_nrows(pnode->mat); ++i) {
        if (i != 0) cout << "; ";
        for (size_type j = 0; j < gmm::mat_ncols(pnode->mat); ++j)
          { if (j != 0) cout << ", "; cout << pnode->mat(i,j); }
      }
      cout << "]]";
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
        if (i%pnode->nbc != 0) cout << ", ";
        else if (i > 0 && (i%pnode->nbc) == 0) cout << "; ";
        ga_print_node(pnode->children[i]);
      }
      cout << "]";
      break;

    default: cout << "Invalid or not taken into account node type"; break;
    }
  }
 

  static void ga_print_tree(const ga_tree &tree) {
    verify_tree(tree.root, 0);
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

        case GA_LBRACKET: // Constant matrix
          {
            ga_tree sub_tree;
            GA_TOKEN_TYPE r_type;
            size_type nbc = 0, n = 0;
            bool foundsemi = false;
            tree.add_matrix(token_pos);
            for(;;) {
              r_type = ga_read_term(expr, pos, sub_tree);
              ++n; if (!foundsemi) ++nbc;
              if (r_type != GA_RBRACKET && r_type != GA_COMMA
                  && r_type != GA_SEMICOLON)
                ga_syntax_error(expr, pos-1, "The constant matrix components "
                                "should be separated by ',' and ';' and "
                                "the constant matrix should be ended by ']'.");
              if (((r_type == GA_SEMICOLON || r_type==GA_RBRACKET) && n != nbc)
                  || (n > nbc))
                ga_syntax_error(expr, pos-1, "Bad constant matrix format. "
                                "The number of components should be equal "
                                "on each line.");
              
              if (r_type == GA_SEMICOLON) { foundsemi = true; n = 0; }
              tree.add_sub_tree(sub_tree);
              if (r_type == GA_RBRACKET) { 
                tree.current_node->nbc = nbc;
                break;
              }
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
        case GA_END: case GA_RPAR: case GA_COMMA:
        case GA_RBRACKET: case GA_SEMICOLON:
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
  // Semantic analysis and ? compilation.
  //=========================================================================

  static void ga_node_analysis(const std::string &expr, ga_tree &tree,
                               pga_tree_node pnode) {

    bool all_sc = true;
    for (size_type i = 0; i < pnode->children.size(); ++i) {
      ga_node_analysis(expr, tree, pnode->children[i]);
      all_sc = all_sc && (pnode->children[i]->node_type == GA_NODE_SCALAR);
    }
    
    switch (pnode->node_type) {
    case GA_NODE_SCALAR: break;
    case GA_NODE_OP:
      switch(pnode->op_type) {
      case GA_PLUS: case GA_MINUS: // + traiter les cas qui ne marceh pas
        if (all_sc) {
          pnode->node_type = GA_NODE_SCALAR;
          pnode->val = pnode->children[0]->val + pnode->children[1]->val *
            ((pnode->op_type == GA_MINUS) ? scalar_type(-1) : scalar_type(1));
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_VECTOR &&
                   pnode->children[1]->node_type == GA_NODE_VECTOR) {
          pnode->node_type = GA_NODE_VECTOR;
          pnode->vec = pnode->children[0]->vec;
          if (gmm::vect_size(pnode->vec)
              != gmm::vect_size(pnode->children[1]->vec))
            ga_syntax_error(expr, pnode->pos, "Addition or substraction "
                            "of vectors of different sizes");
          gmm::add(pnode->children[1]->vec, pnode->vec);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_MATRIX &&
                   pnode->children[1]->node_type == GA_NODE_MATRIX) {
          pnode->node_type = GA_NODE_MATRIX;
          pnode->mat = pnode->children[0]->mat;
          if (gmm::mat_nrows(pnode->mat)
              != gmm::mat_nrows(pnode->children[1]->mat) ||
              gmm::mat_ncols(pnode->mat)
              != gmm::mat_ncols(pnode->children[1]->mat))
            ga_syntax_error(expr, pnode->pos, "Addition or substraction "
                            "of matrices of different sizes");
          gmm::add(pnode->children[1]->mat, pnode->mat);
          tree.clear_children(pnode);
        }
        break;

      case GA_DOT:
        if (pnode->children[0]->node_type == GA_NODE_VECTOR &&
            pnode->children[1]->node_type == GA_NODE_VECTOR) {
          pnode->node_type = GA_NODE_SCALAR;
          if (gmm::vect_size(pnode->children[0]->vec)
              != gmm::vect_size(pnode->children[1]->vec))
            ga_syntax_error(expr, pnode->pos, "Incompatible sizes in "
                            "scalar product.");
          pnode->val = gmm::vect_sp(pnode->children[0]->vec,
                                    pnode->children[1]->vec);
          tree.clear_children(pnode);
        }
        break;

      case GA_COLON:
        if (pnode->children[0]->node_type == GA_NODE_MATRIX &&
            pnode->children[1]->node_type == GA_NODE_MATRIX) {
          pnode->node_type = GA_NODE_SCALAR;
          if (gmm::mat_nrows(pnode->children[0]->mat)
              != gmm::mat_nrows(pnode->children[1]->mat) ||
              gmm::mat_ncols(pnode->children[0]->mat)
              != gmm::mat_ncols(pnode->children[1]->mat))
            ga_syntax_error(expr, pnode->pos, "Incompatible sizes in "
                            "Frobenius product.");
          pnode->val = gmm::vect_sp(pnode->children[0]->mat.as_vector(),
                                    pnode->children[1]->mat.as_vector());
          tree.clear_children(pnode);
        }
        break;

      case GA_MULT:
        if (all_sc) {
          pnode->node_type = GA_NODE_SCALAR;
          pnode->val = pnode->children[0]->val * pnode->children[1]->val;
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_VECTOR &&
                   pnode->children[1]->node_type == GA_NODE_SCALAR) {
          pnode->node_type = GA_NODE_VECTOR;
          pnode->vec = pnode->children[0]->vec;
          gmm::scale(pnode->vec, pnode->children[1]->val);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_MATRIX &&
                   pnode->children[1]->node_type == GA_NODE_SCALAR) {
          pnode->node_type = GA_NODE_MATRIX;
          pnode->mat = pnode->children[0]->mat;
          gmm::scale(pnode->mat, pnode->children[1]->val);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_SCALAR &&
                   pnode->children[1]->node_type == GA_NODE_VECTOR) {
          pnode->node_type = GA_NODE_VECTOR;
          pnode->vec = pnode->children[1]->vec;
          gmm::scale(pnode->vec, pnode->children[0]->val);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_SCALAR &&
                   pnode->children[1]->node_type == GA_NODE_MATRIX) {
          pnode->node_type = GA_NODE_MATRIX;
          pnode->mat = pnode->children[1]->mat;
          gmm::scale(pnode->mat, pnode->children[0]->val);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_MATRIX &&
                   pnode->children[1]->node_type == GA_NODE_VECTOR) {
          pnode->node_type = GA_NODE_VECTOR;
          gmm::resize(pnode->vec, gmm::mat_nrows(pnode->children[0]->mat));
          if (gmm::mat_ncols(pnode->children[0]->mat)
              != gmm::vect_size(pnode->children[1]->vec))
            ga_syntax_error(expr, pnode->pos, "Incompatible sizes in "
                            "matrix-vector multiplication.");
          gmm::mult(pnode->children[0]->mat, pnode->children[1]->vec,
                    pnode->vec);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_MATRIX &&
                   pnode->children[1]->node_type == GA_NODE_MATRIX) {
          pnode->node_type = GA_NODE_MATRIX;
          gmm::resize(pnode->mat, gmm::mat_nrows(pnode->children[0]->mat),
                      gmm::mat_ncols(pnode->children[1]->mat));
          if (gmm::mat_ncols(pnode->children[0]->mat)
              != gmm::mat_nrows(pnode->children[1]->mat))
            ga_syntax_error(expr, pnode->pos, "Incompatible sizes in "
                            "matrix-matrix multiplication.");
          gmm::mult(pnode->children[0]->mat, pnode->children[1]->mat,
                    pnode->mat);
          tree.clear_children(pnode);
        }
        break;

      case GA_DIV: // should control that the second argument is scalar valued
        if (all_sc) {
          pnode->node_type = GA_NODE_SCALAR;
          if (pnode->children[1]->val == scalar_type(0))
            ga_syntax_error(expr, pnode->children[1]->pos, "Division by zero");
          pnode->val = pnode->children[0]->val / pnode->children[1]->val;
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_VECTOR &&
                   pnode->children[1]->node_type == GA_NODE_SCALAR) {
          pnode->node_type = GA_NODE_VECTOR;
          pnode->vec = pnode->children[0]->vec;
          if (pnode->children[1]->val == scalar_type(0))
            ga_syntax_error(expr, pnode->children[1]->pos, "Division by zero");
          gmm::scale(pnode->vec, scalar_type(1)/pnode->children[1]->val);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_MATRIX &&
                   pnode->children[1]->node_type == GA_NODE_SCALAR) {
          pnode->node_type = GA_NODE_MATRIX;
          pnode->mat = pnode->children[0]->mat;
          if (pnode->children[1]->val == scalar_type(0))
            ga_syntax_error(expr, pnode->children[1]->pos, "Division by zero");
          gmm::scale(pnode->mat, scalar_type(1)/pnode->children[1]->val);
          tree.clear_children(pnode);
        }
        break;
        

      case GA_UNARY_MINUS:
        if (all_sc) {
          pnode->node_type = GA_NODE_SCALAR;
          pnode->val = -(pnode->children[0]->val);
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_VECTOR) {
          pnode->node_type = GA_NODE_VECTOR;
          pnode->vec = pnode->children[0]->vec;
          gmm::scale(pnode->vec, scalar_type(-1));
          tree.clear_children(pnode);
        } else if (pnode->children[0]->node_type == GA_NODE_MATRIX) {
          pnode->node_type = GA_NODE_MATRIX;
          pnode->mat = pnode->children[0]->mat;
          gmm::scale(pnode->mat, scalar_type(-1));
          tree.clear_children(pnode);
        }
        break;


      default:GMM_ASSERT1(false, "Unexpected operation. Internal error.");
      }
      break;

    case GA_NODE_C_MATRIX: // should control that all is scalar valued
      if (all_sc) {
        size_type nbc = pnode->nbc;
        size_type nbl = pnode->children.size() / nbc;
        if (nbc == 1 && nbl == 1) {
          pnode->node_type = GA_NODE_SCALAR;
          pnode->val = pnode->children[0]->val;
        } else if (nbc == 1) {
          pnode->node_type = GA_NODE_VECTOR;
          gmm::resize(pnode->vec, nbl);
          for (size_type i = 0; i < nbl; ++i)
            pnode->vec[i] = pnode->children[i]->val;
        } else {
          pnode->node_type = GA_NODE_MATRIX;
          gmm::resize(pnode->mat, nbl, nbc);
          for (size_type i = 0; i < nbl; ++i)
            for (size_type j = 0; j < nbc; ++j)
              pnode->mat(i,j) = pnode->children[i*nbc+j]->val;
        }
        tree.clear_children(pnode);
      }
      break;
      

    default:GMM_ASSERT1(false, "Unexpected node type. Internal error.");
    }

  }

  //=========================================================================
  // Tree simplification
  //=========================================================================
  
  // Reduce all constant operations, gather identical terms

  //=========================================================================
  // Small test for debug
  //=========================================================================

  static  std::string my_expr = "[1,1;1,2]:[1,2;2,1]+[6]";
  // static  std::string my_expr = "p*Trace(Grad_Test_u) + Test_u(1,2)+1.0E-1";
  // static  std::string my_expr = "-(at+(2*3)+2)/3";


  void lex_analysis(void) {
    ga_tree tree;
    ga_read_string(my_expr, tree);
    ga_print_tree(tree);
    ga_node_analysis(my_expr, tree, tree.root);
    ga_print_tree(tree);
  }
  

} /* end of namespace */
