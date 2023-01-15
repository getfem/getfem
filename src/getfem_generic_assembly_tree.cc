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

===========================================================================*/

#include "getfem/getfem_generic_assembly_tree.h"
#include "getfem/getfem_generic_assembly_functions_and_operators.h"
#include "getfem/getfem_generic_assembly_compile_and_exec.h"


namespace getfem {

  //=========================================================================
  // Lexical analysis for the generic assembly language
  //=========================================================================

  static GA_TOKEN_TYPE ga_char_type[256];
  static int ga_operator_priorities[GA_NB_TOKEN_TYPE];

  // Initialize ga_char_type and ga_operator_priorities arrays
  static bool init_ga_char_type() {
    for (int i = 0; i < 256; ++i) ga_char_type[i] = GA_INVALID;
    ga_char_type['+'] = GA_PLUS;        ga_char_type['-']  = GA_MINUS;
    ga_char_type['*'] = GA_MULT;        ga_char_type['/']  = GA_DIV;
    ga_char_type[':'] = GA_COLON;       ga_char_type['\''] = GA_QUOTE;
    ga_char_type['.'] = GA_DOT;         ga_char_type['@']  = GA_TMULT;
    ga_char_type[','] = GA_COMMA;       ga_char_type[';']  = GA_SEMICOLON;
    ga_char_type['('] = GA_LPAR;        ga_char_type[')']  = GA_RPAR;
    ga_char_type['['] = GA_LBRACKET;    ga_char_type[']']  = GA_RBRACKET;
    ga_char_type['_'] = GA_NAME;        ga_char_type['=']  = GA_COLON_EQ;
    for (unsigned i = 'a'; i <= 'z'; ++i)  ga_char_type[i] = GA_NAME;
    for (unsigned i = 'A'; i <= 'Z'; ++i)  ga_char_type[i] = GA_NAME;
    for (unsigned i = '0'; i <= '9'; ++i)  ga_char_type[i] = GA_SCALAR;

    for (unsigned i=0; i < GA_NB_TOKEN_TYPE; ++i) ga_operator_priorities[i] = 0;
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
    ga_operator_priorities[GA_SYM] = 4;
    ga_operator_priorities[GA_SKEW] = 4;
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
    if (type == GA_DOT) {
      if (pos >= expr.size()) return type;
      if (expr[pos] == '*') { ++pos; ++token_length; return GA_DOTMULT; }
      if (expr[pos] == '/') { ++pos; ++token_length; return GA_DOTDIV; }
      if (ga_char_type[unsigned(expr[pos])] != GA_SCALAR) return type;
      fdot = true; type = GA_SCALAR;
    }
    switch (type) {
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
      if (expr.compare(token_pos, token_length, "Sym") == 0)
        return GA_SYM;
      if (expr.compare(token_pos, token_length, "Def") == 0)
        return GA_DEF;
      if (expr.compare(token_pos, token_length, "Skew") == 0)
        return GA_SKEW;
      if (expr.compare(token_pos, token_length, "Trace") == 0)
        return GA_TRACE;
      if (expr.compare(token_pos, token_length, "Deviator") == 0)
        return GA_DEVIATOR;
      if (expr.compare(token_pos, token_length, "Interpolate") == 0)
        return GA_INTERPOLATE;
      if (expr.compare(token_pos, token_length, "Interpolate_derivative") == 0)
        return GA_INTERPOLATE_DERIVATIVE;
      if (expr.compare(token_pos, token_length, "Interpolate_filter") == 0)
        return GA_INTERPOLATE_FILTER;
      if (expr.compare(token_pos, token_length,
                       "Elementary_transformation") == 0)
        return GA_ELEMENTARY;
      if (expr.compare(token_pos, token_length, "Secondary_domain") == 0 ||
          expr.compare(token_pos, token_length, "Secondary_Domain") == 0)
        return GA_SECONDARY_DOMAIN;
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
    case GA_COLON:
      if (pos < expr.size() &&
          ga_char_type[unsigned(expr[pos])] == GA_COLON_EQ) {
        ++pos; return GA_COLON_EQ;
      }
      return type;
    case GA_COLON_EQ:
      return GA_INVALID;
    default: return type;
    }
  }

  //=========================================================================
  // Error handling
  //=========================================================================

  void ga_throw_error_msg(pstring expr, size_type pos,
                          const std::string &msg) {
    int length_before = 70, length_after = 70;
    if (expr && expr->size()) {
      int first = std::max(0, int(pos)-length_before);
      int last = std::min(int(pos)+length_after, int(expr->size()));
      if (last - first < length_before+length_after)
      first = std::max(0, int(pos)-length_before
                       -(length_before+length_after-last+first));
      if (last - first < length_before+length_after)
        last = std::min(int(pos)+length_after
                        +(length_before+length_after-last+first),
                        int(expr->size()));
      if (first > 0) cerr << "...";
      cerr << expr->substr(first, last-first);
      if (last < int(expr->size())) cerr << "...";
      cerr << endl;
      if (first > 0) cerr << "   ";
      if (int(pos) > first)
        cerr << std::setfill ('-') << std::setw(int(pos)-first) << '-'
             << std::setfill (' ');
      cerr << "^" << endl;
    }
    cerr << msg << endl;
  }

  //=========================================================================
  // Tree structure
  //=========================================================================

  void ga_tree_node::mult_test(const pga_tree_node n0, const pga_tree_node n1) {
    
    size_type test0 = n0->test_function_type, test1 = n1->test_function_type;
    if (test0 && test1 && (test0 == test1 || test0 >= 3 || test1 >= 3))
      ga_throw_error(expr, pos,
                     "Incompatibility of test functions in product.");
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

  void ga_tree::add_scalar(scalar_type val, size_type pos, pstring expr) {
    while (current_node && current_node->node_type != GA_NODE_OP)
      current_node = current_node->parent;
    if (current_node) {
      current_node->adopt_child(new ga_tree_node(val, pos, expr));
      current_node = current_node->children.back();
    }
    else {
      GMM_ASSERT1(root == nullptr, "Invalid tree operation");
      current_node = root = new ga_tree_node(val, pos, expr);
      root->parent = nullptr;
    }
  }

  void ga_tree::add_allindices(size_type pos, pstring expr) {
    while (current_node && current_node->node_type != GA_NODE_OP)
      current_node = current_node->parent;
    if (current_node) {
      current_node->adopt_child(new ga_tree_node(GA_NODE_ALLINDICES, pos,expr));
      current_node = current_node->children.back();
    }
    else {
      GMM_ASSERT1(root == nullptr, "Invalid tree operation");
      current_node = root = new ga_tree_node(GA_NODE_ALLINDICES, pos, expr);
      root->parent = nullptr;
    }
  }
  
  void ga_tree::add_name(const char *name, size_type length, size_type pos,
                         pstring expr) {
    while (current_node && current_node->node_type != GA_NODE_OP)
      current_node = current_node->parent;
    if (current_node) {
      current_node->adopt_child(new ga_tree_node(name, length, pos, expr));
      current_node = current_node->children.back();
    }
    else {
      GMM_ASSERT1(root == nullptr, "Invalid tree operation");
      current_node = root = new ga_tree_node(name, length, pos, expr);
      root->parent = nullptr;
    }
  }
  
  void ga_tree::add_sub_tree(ga_tree &sub_tree) {
    if (current_node &&
        (current_node->node_type == GA_NODE_PARAMS ||
         current_node->node_type == GA_NODE_INTERPOLATE_FILTER ||
         current_node->node_type == GA_NODE_C_MATRIX)) {
      GMM_ASSERT1(sub_tree.root, "Invalid tree operation");
      current_node->adopt_child(sub_tree.root);
    } else {
      GMM_ASSERT1(sub_tree.root, "Invalid tree operation");
      while (current_node && current_node->node_type != GA_NODE_OP)
        current_node = current_node->parent;
      if (current_node) {
        current_node->adopt_child(sub_tree.root);
        current_node = sub_tree.root;
      }
      else {
        GMM_ASSERT1(root == nullptr, "Invalid tree operation");
        current_node = root = sub_tree.root;
        root->parent = nullptr;
      }
    }
    sub_tree.root = sub_tree.current_node = nullptr;
  }
  
  void ga_tree::add_params(size_type pos, pstring expr) {
    GMM_ASSERT1(current_node, "internal error");
    while (current_node && current_node->parent &&
           current_node->parent->node_type == GA_NODE_OP &&
           ga_operator_priorities[current_node->parent->op_type] >= 4)
      current_node = current_node->parent;
    pga_tree_node new_node = new ga_tree_node(GA_NODE_PARAMS, pos, expr);
    new_node->parent = current_node->parent;
    if (current_node->parent)
      current_node->parent->replace_child(current_node, new_node);
    else
      root = new_node;
    new_node->adopt_child(current_node);
    current_node = new_node;
  }
  
  void ga_tree::add_matrix(size_type pos, pstring expr) {
    while (current_node && current_node->node_type != GA_NODE_OP)
      current_node = current_node->parent;
    if (current_node) {
      current_node->adopt_child(new ga_tree_node(GA_NODE_C_MATRIX, pos, expr));
      current_node = current_node->children.back();
    }
    else {
      GMM_ASSERT1(root == nullptr, "Invalid tree operation");
      current_node = root = new ga_tree_node(GA_NODE_C_MATRIX, pos, expr);
      root->parent = nullptr;
    }
    current_node->nbc1 = current_node->nbc2 = current_node->nbc3 = 0;
  }
  
  void ga_tree::add_op(GA_TOKEN_TYPE op_type, size_type pos,
                       pstring expr) {
    while (current_node && current_node->parent &&
           current_node->parent->node_type == GA_NODE_OP &&
           ga_operator_priorities[current_node->parent->op_type]
           >= ga_operator_priorities[op_type])
      current_node = current_node->parent;
    pga_tree_node new_node = new ga_tree_node(op_type, pos, expr);
    if (current_node) {
      if (op_type == GA_UNARY_MINUS
          || op_type == GA_SYM || op_type == GA_SKEW
          || op_type == GA_TRACE || op_type == GA_DEVIATOR
          || op_type == GA_PRINT) {
        current_node->adopt_child(new_node);
      } else {
        new_node->parent = current_node->parent;
        if (current_node->parent)
          current_node->parent->replace_child(current_node, new_node);
        else
          root = new_node;
        new_node->adopt_child(current_node);
      }
    } else {
      if (root) new_node->adopt_child(root);
      root = new_node;
      root->parent = nullptr;
    }
    current_node = new_node;
  }

  void ga_tree::clear_node_rec(pga_tree_node pnode) {
    if (pnode) {
      for (pga_tree_node &child : pnode->children)
        clear_node_rec(child);
      delete pnode;
      current_node = nullptr;
    }
  }
  
  void ga_tree::clear_node(pga_tree_node pnode) {
    if (pnode) {
      pga_tree_node parent = pnode->parent;
      if (parent) { // keep all siblings of pnode
        size_type j = 0;
        for (pga_tree_node &sibling : parent->children)
          if (sibling != pnode)
            parent->children[j++] = sibling;
        parent->children.resize(j);
      } else
        root = nullptr;
    }
    clear_node_rec(pnode);
  }
  
  void ga_tree::clear_children(pga_tree_node pnode) {
    for (pga_tree_node &child : pnode->children)
      clear_node_rec(child);
    pnode->children.resize(0);
  }

  void ga_tree::replace_node_by_child(pga_tree_node pnode, size_type i) {
    GMM_ASSERT1(i < pnode->children.size(), "Internal error");
    pga_tree_node child = pnode->children[i];
    child->parent = pnode->parent;
    if (pnode->parent)
      pnode->parent->replace_child(pnode, child);
    else
      root = child;
    current_node = nullptr;
    for (pga_tree_node &sibling : pnode->children)
      if (sibling != child) clear_node_rec(sibling);
    delete pnode;
  }
  
  void ga_tree::copy_node(pga_tree_node pnode, pga_tree_node parent,
                          pga_tree_node &child) {
    GMM_ASSERT1(child == nullptr, "Internal error");
    child = new ga_tree_node();
    *child = *pnode;
    child->parent = parent;
    for (pga_tree_node &grandchild : child->children)
      grandchild = nullptr;
    for (size_type j = 0; j < child->children.size(); ++j)
      copy_node(pnode->children[j], child, child->children[j]);
  }
  
  void ga_tree::duplicate_with_operation(pga_tree_node pnode,
                                         GA_TOKEN_TYPE op_type) {
    pga_tree_node newop = new ga_tree_node(op_type, pnode->pos, pnode->expr);
    newop->children.resize(2, nullptr);
    newop->children[0] = pnode;
    newop->pos = pnode->pos; newop->expr = pnode->expr;
    newop->parent = pnode->parent;
    if (pnode->parent)
      pnode->parent->replace_child(pnode, newop);
    else
      root = newop;
    pnode->parent = newop;
    copy_node(pnode, newop, newop->children[1]);
  }

  void ga_tree::add_child(pga_tree_node pnode, GA_NODE_TYPE node_type) {
    pga_tree_node newnode=new ga_tree_node();
    newnode->pos = pnode->pos; newnode->expr = pnode->expr;
    newnode->node_type = node_type; pnode->adopt_child(newnode);
  }

  void ga_tree::insert_node(pga_tree_node pnode, GA_NODE_TYPE node_type) {
    pga_tree_node newnode = new ga_tree_node();
    newnode->node_type = node_type;
    newnode->parent = pnode->parent;
    newnode->pos = pnode->pos; newnode->expr = pnode->expr;
    if (pnode->parent)
      pnode->parent->replace_child(pnode, newnode);
    else
      root = newnode;
    newnode->adopt_child(pnode);
  }
  
  bool sub_tree_are_equal
  (const pga_tree_node pnode1, const pga_tree_node pnode2,
   const ga_workspace &workspace, int version) {

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
      if (pnode1->name.compare(pnode2->name)) return false;
      break;
    case GA_NODE_PREDEF_FUNC: case GA_NODE_SPEC_FUNC:
      if (pnode1->name.compare(pnode2->name)) return false;
      break;
    case GA_NODE_CONSTANT: case GA_NODE_ZERO:
      if (pnode1->tensor().size() != pnode2->tensor().size()) return false;

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
      if (pnode1->tensor().size() != 1 &&
          pnode1->t.sizes().size() != pnode2->t.sizes().size()) return false;
      for (size_type i = 0; i < pnode1->t.sizes().size(); ++i)
        if (pnode1->t.sizes()[i] != pnode2->t.sizes()[i]) return false;
      for (size_type i = 0; i < pnode1->tensor().size(); ++i)
        if (gmm::abs(pnode1->tensor()[i] - pnode2->tensor()[i]) > 1E-25)
          return false;
      break;
    case GA_NODE_C_MATRIX:
      if (pnode1->children.size() != pnode2->children.size()) return false;
      if (pnode1->nb_test_functions() != pnode2->nb_test_functions())
        return false;
      if (pnode1->t.sizes().size() != pnode2->t.sizes().size()) return false;
      for (size_type i=pnode1->nb_test_functions();
           i<pnode1->t.sizes().size(); ++i)
        if (pnode1->t.sizes()[i] != pnode2->t.sizes()[i]) return false;
      if (pnode1->nbc1 != pnode2->nbc1) return false;
      break;
    case GA_NODE_INTERPOLATE_FILTER:
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name) ||
          pnode1->nbc1 != pnode2->nbc1)
        return false;
      break;
    case GA_NODE_INTERPOLATE_X:
    case GA_NODE_SECONDARY_DOMAIN_X:
    case GA_NODE_INTERPOLATE_ELT_K: case GA_NODE_INTERPOLATE_ELT_B:
    case GA_NODE_INTERPOLATE_NORMAL:
    case GA_NODE_SECONDARY_DOMAIN_NORMAL:
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name))
        return false;
      break;
    case GA_NODE_INTERPOLATE_DERIVATIVE:
      if (pnode1->interpolate_name_der.compare(pnode2->interpolate_name_der))
        return false;
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name) ||
          pnode1->elementary_name.compare(pnode2->elementary_name))
        return false;
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
    case GA_NODE_INTERPOLATE_VAL_TEST: case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST: case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_VAL_TEST: case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST: case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_SECONDARY_DOMAIN_VAL_TEST:
    case GA_NODE_SECONDARY_DOMAIN_GRAD_TEST:
    case GA_NODE_SECONDARY_DOMAIN_HESS_TEST:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST: case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST: case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST: case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST: case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      if (pnode1->interpolate_name.compare(pnode2->interpolate_name) ||
          pnode1->elementary_name.compare(pnode2->elementary_name))
        return false;
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
    case GA_NODE_SECONDARY_DOMAIN_VAL: case GA_NODE_SECONDARY_DOMAIN_GRAD:
    case GA_NODE_SECONDARY_DOMAIN_HESS:  case GA_NODE_SECONDARY_DOMAIN_DIVERG:
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

  static void verify_tree(const pga_tree_node pnode,
                          const pga_tree_node parent) {
    GMM_ASSERT1(pnode->parent == parent,
                "Invalid tree node " << pnode->node_type);
    for (pga_tree_node &child : pnode->children)
      verify_tree(child, pnode);
  }


  static void ga_print_constant_tensor(const pga_tree_node pnode,
                                       std::ostream &str) {
    size_type nt = pnode->nb_test_functions(); // for printing zero tensors
    switch (pnode->tensor_order()) {
    case 0:
      str << (nt ? scalar_type(0) : pnode->tensor()[0]);
      break;

    case 1:
      str << "[";
      for (size_type i = 0; i < pnode->tensor_proper_size(0); ++i) {
        if (i != 0) str << ",";
        str << (nt ? scalar_type(0) : pnode->tensor()[i]);
      }
      str << "]";
      break;

    case 2: case 3: case 4:
      {
        size_type ii(0);
        size_type n0 = pnode->tensor_proper_size(0);
        size_type n1 = pnode->tensor_proper_size(1);
        size_type n2 = ((pnode->tensor_order() > 2) ?
                        pnode->tensor_proper_size(2) : 1);
        size_type n3 = ((pnode->tensor_order() > 3) ?
                        pnode->tensor_proper_size(3) : 1);
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
                str << (nt ? scalar_type(0) : pnode->tensor()[ii++]);
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

    case 5: case 6: case 7: case 8:
      str << "Reshape([";
      for (size_type i = 0; i < pnode->tensor_proper_size(); ++i) {
        if (i != 0) str << ";";
        str << (nt ? scalar_type(0) : pnode->tensor()[i]);
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

  void ga_print_node(const pga_tree_node pnode,
                            std::ostream &str) {
    if (!pnode) return;
    long prec = str.precision(16);

    bool is_interpolate(false), is_elementary(false), is_secondary(false);
    bool is_xfem_plus(false), is_xfem_minus(false);
    switch(pnode->node_type) {
    case GA_NODE_INTERPOLATE:
    case GA_NODE_INTERPOLATE_X:
    case GA_NODE_INTERPOLATE_ELT_K: case GA_NODE_INTERPOLATE_ELT_B:
    case GA_NODE_INTERPOLATE_NORMAL:
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
    case GA_NODE_SECONDARY_DOMAIN:
    case GA_NODE_SECONDARY_DOMAIN_X:
    case GA_NODE_SECONDARY_DOMAIN_NORMAL:
    case GA_NODE_SECONDARY_DOMAIN_VAL:
    case GA_NODE_SECONDARY_DOMAIN_GRAD:
    case GA_NODE_SECONDARY_DOMAIN_HESS:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG:
    case GA_NODE_SECONDARY_DOMAIN_VAL_TEST:
    case GA_NODE_SECONDARY_DOMAIN_GRAD_TEST:
    case GA_NODE_SECONDARY_DOMAIN_HESS_TEST:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST:
      str << "Secondary_domain(";
      is_secondary = true;
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
    case GA_NODE_SECONDARY_DOMAIN_GRAD:
    case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_SECONDARY_DOMAIN_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_GRAD_TEST:
      str << "Grad_";
      break;
    case GA_NODE_HESS:
    case GA_NODE_INTERPOLATE_HESS:
    case GA_NODE_ELEMENTARY_HESS:
    case GA_NODE_SECONDARY_DOMAIN_HESS:
    case GA_NODE_XFEM_PLUS_HESS:
    case GA_NODE_XFEM_MINUS_HESS:
    case GA_NODE_HESS_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST:
    case GA_NODE_SECONDARY_DOMAIN_HESS_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST:
      str << "Hess_";
      break;
    case GA_NODE_DIVERG:
    case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG:
    case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_DIVERG:
    case GA_NODE_DIVERG_TEST:
    case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST:
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
        } else if (pnode->op_type == GA_SYM) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          str << "Sym("; ga_print_node(pnode->children[0], str); str << ")";
        } else if (pnode->op_type == GA_SKEW) {
          GMM_ASSERT1(pnode->children.size() == 1, "Invalid tree");
          str << "Skew("; ga_print_node(pnode->children[0], str); str << ")";
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
    case GA_NODE_INTERPOLATE_X: case GA_NODE_SECONDARY_DOMAIN_X:
      str << "X";
      break;
    case GA_NODE_INTERPOLATE_NORMAL: case GA_NODE_SECONDARY_DOMAIN_NORMAL:
      str << "Normal";
      break;
    case GA_NODE_INTERPOLATE_ELT_K:
      str << "element_K";
      break;
    case GA_NODE_INTERPOLATE_ELT_B:
      str << "element_B";
      break;
    case GA_NODE_INTERPOLATE_DERIVATIVE:
      str << (pnode->test_function_type == 1 ? "Test_" : "Test2_")
          << "Interpolate_derivative(" << pnode->interpolate_name_der << ","
          << pnode->name;
      if (pnode->interpolate_name.size())
        str << "," << pnode->interpolate_name;
      str <<  ")";
      break;
    case GA_NODE_INTERPOLATE:
    case GA_NODE_ELEMENTARY:
    case GA_NODE_SECONDARY_DOMAIN:
    case GA_NODE_XFEM_PLUS:
    case GA_NODE_XFEM_MINUS:
    case GA_NODE_VAL:
    case GA_NODE_INTERPOLATE_VAL:
    case GA_NODE_ELEMENTARY_VAL:
    case GA_NODE_SECONDARY_DOMAIN_VAL:
    case GA_NODE_XFEM_PLUS_VAL:
    case GA_NODE_XFEM_MINUS_VAL:
    case GA_NODE_GRAD:
    case GA_NODE_INTERPOLATE_GRAD:
    case GA_NODE_SECONDARY_DOMAIN_GRAD:
    case GA_NODE_ELEMENTARY_GRAD:
    case GA_NODE_XFEM_PLUS_GRAD:
    case GA_NODE_XFEM_MINUS_GRAD:
    case GA_NODE_HESS:
    case GA_NODE_INTERPOLATE_HESS:
    case GA_NODE_SECONDARY_DOMAIN_HESS:
    case GA_NODE_ELEMENTARY_HESS:
    case GA_NODE_XFEM_PLUS_HESS:
    case GA_NODE_XFEM_MINUS_HESS:
    case GA_NODE_DIVERG:
    case GA_NODE_INTERPOLATE_DIVERG:
    case GA_NODE_ELEMENTARY_DIVERG:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG:
    case GA_NODE_XFEM_PLUS_DIVERG:
    case GA_NODE_XFEM_MINUS_DIVERG:
      str << pnode->name;
      break;
    case GA_NODE_VAL_TEST:
    case GA_NODE_INTERPOLATE_VAL_TEST:
    case GA_NODE_ELEMENTARY_VAL_TEST:
    case GA_NODE_SECONDARY_DOMAIN_VAL_TEST:
    case GA_NODE_XFEM_PLUS_VAL_TEST:
    case GA_NODE_XFEM_MINUS_VAL_TEST:
    case GA_NODE_GRAD_TEST:
    case GA_NODE_INTERPOLATE_GRAD_TEST:
    case GA_NODE_ELEMENTARY_GRAD_TEST:
    case GA_NODE_SECONDARY_DOMAIN_GRAD_TEST:
    case GA_NODE_XFEM_PLUS_GRAD_TEST:
    case GA_NODE_XFEM_MINUS_GRAD_TEST:
    case GA_NODE_HESS_TEST:
    case GA_NODE_INTERPOLATE_HESS_TEST:
    case GA_NODE_ELEMENTARY_HESS_TEST:
    case GA_NODE_SECONDARY_DOMAIN_HESS_TEST:
    case GA_NODE_XFEM_PLUS_HESS_TEST:
    case GA_NODE_XFEM_MINUS_HESS_TEST:
    case GA_NODE_DIVERG_TEST:
    case GA_NODE_INTERPOLATE_DIVERG_TEST:
    case GA_NODE_ELEMENTARY_DIVERG_TEST:
    case GA_NODE_SECONDARY_DOMAIN_DIVERG_TEST:
    case GA_NODE_XFEM_PLUS_DIVERG_TEST:
    case GA_NODE_XFEM_MINUS_DIVERG_TEST:
      str << (pnode->test_function_type == 1 ? "Test_" : "Test2_")
          << pnode->name;
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
      for (size_type i = 1; i < pnode->children.size(); ++i) {
        if (i > 1) str << ", ";
        ga_print_node(pnode->children[i], str);
      }
      str << ")";
      break;

    case GA_NODE_NAME:
      str << pnode->name;
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;
      
    case GA_NODE_MACRO_PARAM:
      if (pnode->nbc2 == 1) str << "Grad_";
      if (pnode->nbc2 == 2) str << "Hess_";
      if (pnode->nbc2 == 3) str << "Div_";
      if (pnode->nbc3 == 1) str << "Test_";
      if (pnode->nbc3 == 2) str << "Test2_";
      str << "P" << pnode->nbc1;
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_RESHAPE:
      str << "Reshape";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;
      
    case GA_NODE_CROSS_PRODUCT:
      str << "Cross_product";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;
      
    case GA_NODE_SWAP_IND:
      str << "Swap_indices";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_IND_MOVE_LAST:
      str << "Index_move_last";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_CONTRACT:
      str << "Contract";
      GMM_ASSERT1(pnode->children.size() == 0, "Invalid tree");
      break;

    case GA_NODE_C_MATRIX:
      GMM_ASSERT1(pnode->children.size(), "Invalid tree");
      GMM_ASSERT1(pnode->nbc1 == pnode->tensor_order(), "Invalid C_MATRIX");
      switch (pnode->tensor_order()) {
      case 0:
        ga_print_node(pnode->children[0], str);
        break;

      case 1:
        str << "[";
        for (size_type i = 0; i < pnode->tensor_proper_size(0); ++i) {
          if (i != 0) str << ",";
          ga_print_node(pnode->children[i], str);
        }
        str << "]";
        break;

      case 2: case 3: case 4:
        {
          size_type ii(0);
          size_type n0 = pnode->tensor_proper_size(0);
          size_type n1 = pnode->tensor_proper_size(1);
          size_type n2 = ((pnode->tensor_order() > 2) ?
                          pnode->tensor_proper_size(2) : 1);
          size_type n3 = ((pnode->tensor_order() > 3) ?
                          pnode->tensor_proper_size(3) : 1);
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
                  ga_print_node(pnode->children[ii++], str);
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

      case 5: case 6: case 7: case 8:
        str << "Reshape([";
        for (size_type i = 0; i < pnode->tensor_proper_size(); ++i) {
          if (i != 0) str << ";";
          ga_print_node(pnode->children[i], str);
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
      break;

    default:
      str << "Invalid or not taken into account node type "
           << pnode->node_type;
      break;
    }

    if (is_interpolate)
      str << "," << pnode->interpolate_name << ")";
    else if (is_elementary) {
      str << "," << pnode->elementary_name;
      if (pnode->name.compare(pnode->elementary_target) != 0)
        str << "," << pnode->elementary_target;
      str << ")";
    } else if (is_secondary)
      str << ")";    
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

  // 0 : ok
  // 1 : function or operator name or "X"
  // 2 : reserved prefix Grad, Hess, Div, Derivative_ Test and Test2
  // 3 : reserved prefix Dot and Previous
  int ga_check_name_validity(const std::string &name) {
    if (name.compare(0, 11, "Derivative_") == 0)
      return 2;
    
    const ga_predef_operator_tab &PREDEF_OPERATORS
      = dal::singleton<ga_predef_operator_tab>::instance(0);
    const ga_spec_function_tab &SPEC_FUNCTIONS
      = dal::singleton<ga_spec_function_tab>::instance(0);
    const ga_spec_op_tab &SPEC_OP
      = dal::singleton<ga_spec_op_tab>::instance(0);
    const ga_predef_function_tab &PREDEF_FUNCTIONS
      = dal::singleton<ga_predef_function_tab>::instance(0);

    if (SPEC_OP.find(name) != SPEC_OP.end())
      return 1;

    if (PREDEF_FUNCTIONS.find(name) != PREDEF_FUNCTIONS.end())
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
  // Structure dealing with macros.
  //=========================================================================

  ga_macro::ga_macro() : ptree(new ga_tree), nbp(0) {}
  ga_macro::~ga_macro() { delete ptree; }
  ga_macro::ga_macro(const std::string &name, const ga_tree &t, size_type nbp_)
    : ptree(new ga_tree(t)), macro_name_(name), nbp(nbp_) {}
  ga_macro::ga_macro(const ga_macro &gam)
    : ptree(new ga_tree(gam.tree())), macro_name_(gam.name()),
      nbp(gam.nb_params()) {}
  ga_macro &ga_macro::operator =(const ga_macro &gam) {
    delete ptree; ptree = new ga_tree(gam.tree());
    macro_name_ = gam.name();
    nbp = gam.nb_params();
    return *this;
  }

  static void ga_replace_macro_params
  (ga_tree &tree, pga_tree_node pnode,
   const std::vector<pga_tree_node> &children) {
    if (!pnode) return;
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_replace_macro_params(tree, pnode->children[i], children);
    
    if (pnode->node_type == GA_NODE_MACRO_PARAM) {
      size_type po = pnode->nbc2;
      size_type pt = pnode->nbc3;
      GMM_ASSERT1(pnode->nbc1+1 < children.size(), "Internal error");
      pga_tree_node pchild = children[pnode->nbc1+1];

      if (po || pt || pnode->op_type != GA_NAME) {
        if (!(pchild->children.empty()) || pchild->node_type != GA_NODE_NAME)
          ga_throw_error(pchild->expr, pchild->pos, "Error in macro "
                         "expansion. Only variable name are allowed for macro "
                         "parameter preceded by Grad_ Hess_ Test_ or Test2_ "
                         "prefixes.");
        switch(pnode->op_type) {
        case GA_NAME : pnode->node_type = GA_NODE_NAME; break;
        case GA_INTERPOLATE : pnode->node_type = GA_NODE_INTERPOLATE; break;
        case GA_INTERPOLATE_DERIVATIVE :
          pnode->node_type = GA_NODE_INTERPOLATE_DERIVATIVE; break;
        case GA_ELEMENTARY : pnode->node_type = GA_NODE_ELEMENTARY; break;
        case GA_SECONDARY_DOMAIN :
          pnode->node_type = GA_NODE_SECONDARY_DOMAIN; break;
        case GA_XFEM_PLUS : pnode->node_type = GA_NODE_XFEM_PLUS; break;
        case GA_XFEM_MINUS: pnode->node_type = GA_NODE_XFEM_MINUS; break;
        default:break;
        }
        pnode->name = pchild->name;
        if (pt == 1) pnode->name = "Test_" + pnode->name;
        if (pt == 2) pnode->name = "Test2_" + pnode->name;
        if (po == 1) pnode->name = "Grad_" + pnode->name;
        if (po == 2) pnode->name = "Hess_" + pnode->name;
        if (po == 3) pnode->name = "Div_" + pnode->name;
      } else {
        pga_tree_node pnode_old = pnode;
        pnode = nullptr;
        tree.copy_node(pchild, pnode_old->parent, pnode);
        if (pnode_old->parent)
          pnode_old->parent->replace_child(pnode_old, pnode);
        else
          tree.root = pnode;
        GMM_ASSERT1(pnode_old->children.empty(), "Internal error");
        delete pnode_old;
      }
    }
  }
  
  static void ga_expand_macro(ga_tree &tree, pga_tree_node pnode,
                              const ga_macro_dictionary &macro_dict) {
    if (!pnode) return;
    
    if (pnode->node_type == GA_NODE_PARAMS) {
      
      for (size_type i = 1; i < pnode->children.size(); ++i)
        ga_expand_macro(tree, pnode->children[i], macro_dict);

      if (pnode->children[0]->node_type != GA_NODE_NAME) {
        ga_expand_macro(tree, pnode->children[0], macro_dict);
      } else {

        if (macro_dict.macro_exists(pnode->children[0]->name)) {

          const ga_macro &gam = macro_dict.get_macro(pnode->children[0]->name);

          if (gam.nb_params()==0) { // Macro without parameters
            pga_tree_node pnode_old = pnode->children[0];
            pnode->children[0] = nullptr;
            tree.copy_node(gam.tree().root,
                           pnode_old->parent,pnode->children[0]);
            GMM_ASSERT1(pnode_old->children.empty(), "Internal error");
            delete pnode_old;
          } else { // Macro with parameters
            if (gam.nb_params()+1 != pnode->children.size())
              ga_throw_error(pnode->expr, pnode->pos,
                             "Bad number of parameters in the use of macro '"
                             << gam.name() << "'. Expected " << gam.nb_params()
                             << " found " << pnode->children.size()-1 << ".");

            pga_tree_node pnode_old = pnode;
            pnode = nullptr;
            tree.copy_node(gam.tree().root, pnode_old->parent, pnode);
            if (pnode_old->parent)
              pnode_old->parent->replace_child(pnode_old, pnode);
            else
              tree.root = pnode;
            ga_replace_macro_params(tree, pnode, pnode_old->children);
            tree.clear_node_rec(pnode_old);
          }
        }
      }

    } else if (pnode->node_type == GA_NODE_NAME &&
               macro_dict.macro_exists(pnode->name)) {
      // Macro without parameters
      const ga_macro &gam = macro_dict.get_macro(pnode->name);
      if (gam.nb_params() != 0)
        ga_throw_error(pnode->expr, pnode->pos,
                       "Bad number of parameters in the use of macro '"
                       << gam.name() << "'. Expected " << gam.nb_params()
                       << " none found.");

      pga_tree_node pnode_old = pnode;
      pnode = nullptr;
      tree.copy_node(gam.tree().root, pnode_old->parent, pnode);
      if (pnode_old->parent)
        pnode_old->parent->replace_child(pnode_old, pnode);
      else
        tree.root = pnode;
      GMM_ASSERT1(pnode_old->children.empty(), "Internal error");
      delete pnode_old;
    } else {
      for (size_type i = 0; i < pnode->children.size(); ++i)
        ga_expand_macro(tree, pnode->children[i], macro_dict);
    }
  }

  static void ga_mark_macro_params_rec(const pga_tree_node pnode,
                                       const std::vector<std::string> &params) {
    if (!pnode) return;
    for (size_type i = 0; i < pnode->children.size(); ++i)
      ga_mark_macro_params_rec(pnode->children[i], params);
    
    if (pnode->node_type == GA_NODE_NAME ||
        pnode->node_type == GA_NODE_INTERPOLATE ||
        pnode->node_type == GA_NODE_ELEMENTARY ||
        pnode->node_type == GA_NODE_SECONDARY_DOMAIN ||
        pnode->node_type == GA_NODE_XFEM_PLUS ||
        pnode->node_type == GA_NODE_XFEM_MINUS) {
      std::string name = pnode->name;
      size_type po = ga_parse_prefix_operator(name);
      size_type pt = ga_parse_prefix_test(name);

      for (size_type i = 0; i < params.size(); ++i)
        if (name.compare(params[i]) == 0) {
          pnode->name = name;
          switch(pnode->node_type) {
          case GA_NODE_NAME : pnode->op_type = GA_NAME; break;
          case GA_NODE_INTERPOLATE : pnode->op_type = GA_INTERPOLATE; break;
          case GA_NODE_INTERPOLATE_DERIVATIVE :
            pnode->op_type = GA_INTERPOLATE_DERIVATIVE; break;
          case GA_NODE_ELEMENTARY : pnode->op_type = GA_ELEMENTARY; break;
          case GA_NODE_SECONDARY_DOMAIN :
            pnode->op_type = GA_SECONDARY_DOMAIN; break;
          case GA_NODE_XFEM_PLUS : pnode->op_type = GA_XFEM_PLUS; break;
          case GA_NODE_XFEM_MINUS: pnode->op_type = GA_XFEM_MINUS; break;
          default:break;
          }
          pnode->node_type = GA_NODE_MACRO_PARAM;
          pnode->nbc1 = i; pnode->nbc2 = po; pnode->nbc3 = pt;
        }
    }
  }
  
  static void ga_mark_macro_params(ga_macro &gam,
                                   const std::vector<std::string> &params,
                                   const ga_macro_dictionary &macro_dict) {
    if (gam.tree().root) {
      ga_mark_macro_params_rec(gam.tree().root, params);
      ga_expand_macro(gam.tree(), gam.tree().root, macro_dict);
    }
  }

  bool ga_macro_dictionary::macro_exists(const std::string &name) const {
    if (macros.find(name) != macros.end()) return true;
    if (parent && parent->macro_exists(name)) return true;
    return false;
  }

  const ga_macro &
  ga_macro_dictionary::get_macro(const std::string &name) const {
    auto it = macros.find(name);
    if (it != macros.end()) return it->second;
    if (parent) return parent->get_macro(name);
    GMM_ASSERT1(false, "Undefined macro");
  }

  void ga_macro_dictionary::add_macro(const ga_macro &gam)
  { macros[gam.name()] = gam; }

  void ga_macro_dictionary::add_macro(const std::string &name,
                                      const std::string &expr)
  { ga_tree tree; ga_read_string_reg("Def "+name+":="+expr, tree, *this); }

  void ga_macro_dictionary::del_macro(const std::string &name) {
    auto it = macros.find(name);
    GMM_ASSERT1(it != macros.end(), "Undefined macro (at this level)");
    macros.erase(it);
  }
  
  
  //=========================================================================
  // Syntax analysis for the generic assembly language
  //=========================================================================

  // Read a term with an (implicit) pushdown automaton.
  static GA_TOKEN_TYPE ga_read_term(pstring expr, size_type &pos,
                                    ga_tree &tree,
                                    ga_macro_dictionary &macro_dict) {
    size_type token_pos, token_length;
    GA_TOKEN_TYPE t_type;
    int state = 1; // 1 = reading term, 2 = reading after term

    for (;;) {

      t_type = ga_get_token(*expr, pos, token_pos, token_length);

      switch (state) {

      case 1:
        switch (t_type) {
        case GA_SCALAR:
          {
            char *endptr; const char *nptr = &((*expr)[token_pos]);
            scalar_type s_read = ::strtod(nptr, &endptr);
            if (endptr == nptr)
              ga_throw_error(expr, token_pos, "Bad numeric format.");
            tree.add_scalar(s_read, token_pos, expr);
          }
          state = 2; break;

        case GA_COLON:
          tree.add_allindices(token_pos, expr);
          state = 2; break;

        case GA_NAME:
          tree.add_name(&((*expr)[token_pos]), token_length, token_pos, expr);
          state = 2; break;

        case GA_MINUS: // unary -
          tree.add_op(GA_UNARY_MINUS, token_pos, expr);
          state = 1; break;

        case GA_PLUS:  // unary +
          state = 1; break;

        case GA_SYM:
          tree.add_op(GA_SYM, token_pos, expr);
          state = 1; break;

        case GA_SKEW:
          tree.add_op(GA_SKEW, token_pos, expr);
          state = 1; break;

        case GA_TRACE:
          tree.add_op(GA_TRACE, token_pos, expr);
          state = 1; break;

        case GA_DEVIATOR:
          tree.add_op(GA_DEVIATOR, token_pos, expr);
          state = 1; break;

        case GA_DEF:
          {
            ga_macro gam;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "Macro definition should begin with macro name");
            gam.name() = std::string(&((*expr)[token_pos]), token_length);
            if (ga_check_name_validity(gam.name()))
              ga_throw_error(expr, pos-1, "Invalid macro name.")
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            std::vector<std::string> params;
            if (t_type == GA_LPAR) {
              t_type = ga_get_token(*expr, pos, token_pos, token_length);
              while (t_type == GA_NAME) {
                params.push_back(std::string(&((*expr)[token_pos]),
                                             token_length));
                if (ga_check_name_validity(params.back()))
                  ga_throw_error(expr, pos-1, "Invalid macro parameter name.");
                for (size_type i = 0; i+1 < params.size(); ++i)
                  if (params.back().compare(params[i]) == 0)
                    ga_throw_error(expr, pos-1,
                                   "Invalid repeated macro parameter name.");
                t_type = ga_get_token(*expr, pos, token_pos, token_length);
                if (t_type == GA_COMMA)
                  t_type = ga_get_token(*expr, pos, token_pos, token_length);
              }
              if (t_type != GA_RPAR)
                ga_throw_error(expr, pos-1,
                              "Missing right parenthesis in macro definition.");
              t_type = ga_get_token(*expr, pos, token_pos, token_length);
            }
            if (t_type != GA_COLON_EQ)
              ga_throw_error(expr, pos-1, "Missing := for macro definition.");

            t_type = ga_read_term(expr, pos, gam.tree(), macro_dict);
            if (gam.tree().root)
              ga_expand_macro(gam.tree(), gam.tree().root, macro_dict);
            gam.nb_params() = params.size();
            if (params.size())
              ga_mark_macro_params(gam, params, macro_dict);
            macro_dict.add_macro(gam);
            
            // cout << "macro \"" << gam.name() << "\" registered with "
            //      << gam.nb_params() << " params  := "
            //      << ga_tree_to_string(gam.tree()) << endl;

            if (t_type == GA_END) return t_type;
            else if (t_type != GA_SEMICOLON)
              ga_throw_error(expr, pos-1,
                             "Syntax error at the end of macro definition.");
            state = 1;
          }
          break;

        case GA_INTERPOLATE:
          {
            tree.add_scalar(scalar_type(0), token_pos, expr);
            tree.current_node->node_type = GA_NODE_INTERPOLATE;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1, "Missing interpolate arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "First argument of Interpolate should be a "
                             "variable, test function, X or Normal.");
            tree.current_node->name = std::string(&((*expr)[token_pos]),
                                                  token_length);

            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_COMMA)
              ga_throw_error(expr, pos, "Bad format for Interpolate "
                             "arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "Second argument of Interpolate should be a "
                             "transformation name.");
            tree.current_node->interpolate_name
              = std::string(&((*expr)[token_pos]), token_length);
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "interpolate arguments.");
            state = 2;
          }
          break;

        case GA_INTERPOLATE_DERIVATIVE:
          {
            tree.add_scalar(scalar_type(0), token_pos, expr);
            tree.current_node->node_type = GA_NODE_INTERPOLATE_DERIVATIVE;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,
                             "Missing Interpolate_derivative arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "First argument of Interpolate should the "
                             "interpolate transformtion name ");
            tree.current_node->interpolate_name_der
              = std::string(&((*expr)[token_pos]), token_length);
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_COMMA)
              ga_throw_error(expr, pos, "Bad format for Interpolate_derivative "
                             "arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "Second argument of Interpolate should be a "
                             "variable name.");
            tree.current_node->name
              = std::string(&((*expr)[token_pos]), token_length);
            
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            tree.current_node->interpolate_name = "";
            if (t_type == GA_COMMA) {
              t_type = ga_get_token(*expr, pos, token_pos, token_length);
              if (t_type != GA_NAME)
                ga_throw_error(expr, pos,
                               "Third argument of Interpolate should be a "
                               "interpolate transformation name.");
              tree.current_node->interpolate_name
                = std::string(&((*expr)[token_pos]), token_length);
              t_type = ga_get_token(*expr, pos, token_pos, token_length);
            }
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Interpolate_derivative arguments.");
            state = 2;
          }
          break;

        case GA_ELEMENTARY:
          {
            tree.add_scalar(scalar_type(0), token_pos, expr);
            tree.current_node->node_type = GA_NODE_ELEMENTARY;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,
                             "Missing Elementary_transformation arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "First argument of Elementary_transformation "
                             "should be a variable or a test function.");
            tree.current_node->name = std::string(&((*expr)[token_pos]),
                                                  token_length);
            tree.current_node->elementary_target = tree.current_node->name;

            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_COMMA)
              ga_throw_error(expr, pos, "Bad format for "
                             "Elementary_transformation arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "Second argument of Elementary_transformation "
                             "should be a transformation name.");
            tree.current_node->elementary_name
              = std::string(&((*expr)[token_pos]), token_length);
            t_type = ga_get_token(*expr, pos, token_pos, token_length);

            if (t_type == GA_COMMA) {
              t_type = ga_get_token(*expr, pos, token_pos, token_length);
              if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "Third argument of Elementary_transformation "
                             "should be a variable or data name.");
              
              tree.current_node->elementary_target =
                std::string(&((*expr)[token_pos]), token_length);
              t_type = ga_get_token(*expr, pos, token_pos, token_length);
            }
            
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Elementary_transformation arguments.");
            state = 2;
          }
          break;

        case GA_SECONDARY_DOMAIN:
          {
            tree.add_scalar(scalar_type(0), token_pos, expr);
            tree.current_node->node_type = GA_NODE_SECONDARY_DOMAIN;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,"Missing Secondary_domain arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "First argument of Secondary_domain should be a "
                             "variable, test function, X or Normal.");
            tree.current_node->name = std::string(&((*expr)[token_pos]),
                                                  token_length);
            tree.current_node->interpolate_name =  tree.secondary_domain;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Secondary_domain arguments.");
            state = 2;
          }
          break;

        case GA_XFEM_PLUS:
          {
            tree.add_scalar(scalar_type(0), token_pos, expr);
            tree.current_node->node_type = GA_NODE_XFEM_PLUS;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,
                             "Missing Xfem_plus arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "The argument of Xfem_plus should be a "
                             "variable or a test function.");
            tree.current_node->name = std::string(&((*expr)[token_pos]),
                                                  token_length);
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Xfem_plus argument.");
            state = 2;
          }
          break;

        case GA_XFEM_MINUS:
          {
            tree.add_scalar(scalar_type(0), token_pos, expr);
            tree.current_node->node_type = GA_NODE_XFEM_MINUS;
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1,
                             "Missing Xfem_minus arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos,
                             "The argument of Xfem_minus should be a "
                             "variable or a test function.");
            tree.current_node->name = std::string(&((*expr)[token_pos]),
                                                  token_length);
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Missing a parenthesis after "
                             "Xfem_minus argument.");
            state = 2;
          }
          break;

        case GA_INTERPOLATE_FILTER:
          {
            tree.add_scalar(scalar_type(0), token_pos, expr);
            tree.current_node->node_type = GA_NODE_INTERPOLATE_FILTER;
            tree.current_node->nbc1 = size_type(-1);
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_LPAR)
              ga_throw_error(expr, pos-1, "Missing interpolate arguments.");
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_NAME)
              ga_throw_error(expr, pos, "First argument of Interpolate_filter "
                             "should be a transformation name.");
            tree.current_node->interpolate_name
              = std::string(&((*expr)[token_pos]), token_length);
            t_type = ga_get_token(*expr, pos, token_pos, token_length);
            if (t_type != GA_COMMA)
              ga_throw_error(expr, pos,
                             "Bad format for Interpolate_filter arguments.");
            ga_tree sub_tree;
            t_type = ga_read_term(expr, pos, sub_tree, macro_dict);
            if (t_type != GA_RPAR && t_type != GA_COMMA)
              ga_throw_error(expr, pos-1,
                             "Bad format for Interpolate_filter arguments.");
            tree.add_sub_tree(sub_tree);
            if (t_type == GA_COMMA) {
               ga_tree sub_tree2;
               t_type = ga_read_term(expr, pos, sub_tree2, macro_dict);
               tree.add_sub_tree(sub_tree2);
            }
            if (t_type != GA_RPAR)
              ga_throw_error(expr, pos-1, "Unbalanced parenthesis.");
            state = 2;
          }
          break;

        case GA_PRINT:
          tree.add_op(GA_PRINT, token_pos, expr);
          state = 1; break;

        case GA_LPAR: // Parenthesed expression
          {
            ga_tree sub_tree;
            GA_TOKEN_TYPE r_type;
            r_type = ga_read_term(expr, pos, sub_tree, macro_dict);
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
            bool foundcomma(false), foundsemi(false);

            r_type = ga_read_term(expr, pos, sub_tree, macro_dict);
            size_type nb_comp = 0;
            tree.add_matrix(token_pos, expr);

            if (sub_tree.root->node_type == GA_NODE_C_MATRIX) { // nested format
              bgeot::multi_index mii;
              do {
                if (nb_comp) {
                  sub_tree.clear();
                  r_type = ga_read_term(expr, pos, sub_tree, macro_dict);
                }
                // in the nested format only "," and "]" are expected
                if (sub_tree.root->node_type != GA_NODE_C_MATRIX ||
                    (r_type != GA_COMMA && r_type != GA_RBRACKET))
                  ga_throw_error(expr, pos-1, "Bad explicit "
                                 "vector/matrix/tensor format.");

                // convert a row vector [a,b] to a column vector [a;b]
                if (sub_tree.root->marked &&
                    sub_tree.root->tensor().sizes()[0] == 1 &&
                    sub_tree.root->tensor().size() != 1) {
                  bgeot::multi_index mi = sub_tree.root->tensor().sizes();
                  for (size_type i = mi.size()-1; i > 0; i--)
                    mi[i-1] = mi[i];
                  mi.pop_back();
                  sub_tree.root->tensor().adjust_sizes(mi);
                }
                if (!nb_comp) mii = sub_tree.root->tensor().sizes();
                else {
                  const bgeot::multi_index &mi=sub_tree.root->tensor().sizes();
                  bool cmp = true;
                  if (mii.size() == mi.size()) {
                     for (size_type i = 0; i < mi.size(); ++i)
                       if (mi[i] != mii[i]) cmp = false;
                  } else cmp = false;
                  if (!cmp)
                    ga_throw_error(expr, pos-1, "Bad explicit "
                                   "vector/matrix/tensor format.");
                }
                for (size_type i = 0; i < sub_tree.root->children.size(); ++i) {
                  sub_tree.root->children[i]->parent = tree.current_node;
                  tree.current_node->children.push_back
                    (sub_tree.root->children[i]);
                }
                sub_tree.root->children.resize(0);
                nb_comp++;
              } while (r_type != GA_RBRACKET);
              tree.current_node->marked = false;
              mii.push_back(nb_comp);
              tree.current_node->tensor().adjust_sizes(mii);
            } else { // non nested format
              do {
                if (nb_comp) {
                  sub_tree.clear();
                  r_type = ga_read_term(expr, pos, sub_tree, macro_dict);
                }
                nb_comp++;

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

              } while (r_type != GA_RBRACKET);
              bgeot::multi_index mi;
              nbc1 = tree.current_node->nbc1; nbc2 = tree.current_node->nbc2;
              nbc3 = tree.current_node->nbc3;

              size_type nbl = tree.current_node->children.size()
                / (nbc2 * nbc1 * nbc3);
              switch(tensor_order) {
              case 1:
                /* mi.push_back(1); */ mi.push_back(nbc1); break;
              case 2:
                mi.push_back(nbl); if (nbc1 > 1) mi.push_back(nbc1); break;
              case 3:
                mi.push_back(nbl); mi.push_back(nbc2);
                mi.push_back(nbc1);
                break;
              case 4:
                mi.push_back(nbl); mi.push_back(nbc3);
                mi.push_back(nbc2); mi.push_back(nbc1);
                break;
              default: GMM_ASSERT1(false, "Internal error");
              }
              tree.current_node->tensor().adjust_sizes(mi);
              std::vector<pga_tree_node> children = tree.current_node->children;
              auto it = tree.current_node->children.begin();
              for (size_type i = 0; i < nbc1; ++i)
                for (size_type j = 0; j < nbc2; ++j)
                  for (size_type k = 0; k < nbc3; ++k)
                    for (size_type l = 0; l < nbl; ++l, ++it)
                      *it = children[i+nbc1*(j+nbc2*(k+nbc3*l))];
              tree.current_node->marked = true;
            }
          }
          tree.current_node->nbc1 = tree.current_node->tensor().sizes().size();
          state = 2;
          break;

        default:
          ga_throw_error(expr, token_pos, "Unexpected token.");
        }
        break;

      case 2:
        switch (t_type) {
        case GA_PLUS: case GA_MINUS: case GA_MULT: case GA_DIV:
        case GA_COLON: case GA_DOT: case GA_DOTMULT: case GA_DOTDIV:
        case GA_TMULT:
          tree.add_op(t_type, token_pos, expr);
          state = 1; break;
        case GA_QUOTE:
          tree.add_op(t_type, token_pos, expr);
          state = 2; break;
        case GA_END: case GA_RPAR: case GA_COMMA: case GA_DCOMMA:
        case GA_RBRACKET: case GA_SEMICOLON: case GA_DSEMICOLON:
          return t_type;
        case GA_LPAR: // Parameter list
          {
            ga_tree sub_tree;
            GA_TOKEN_TYPE r_type;
            tree.add_params(token_pos, expr);
            do {
              r_type = ga_read_term(expr, pos, sub_tree, macro_dict);
              if (r_type != GA_RPAR && r_type != GA_COMMA)
                ga_throw_error(expr, pos-((r_type != GA_END)?1:0),
                               "Parameters should be separated "
                               "by ',' and parameter list ended by ')'.");
              tree.add_sub_tree(sub_tree);
            } while (r_type != GA_RPAR);
            state = 2;
          }
          break;

        default:
          ga_throw_error(expr, token_pos, "Unexpected token.");
        }
        break;
      }
    }

    return GA_INVALID;
  }

  // Syntax analysis of a string. Conversion to a tree. register the macros.
  void ga_read_string_reg(const std::string &expr, ga_tree &tree,
                          ga_macro_dictionary &macro_dict) {
    size_type pos = 0, token_pos, token_length;
    tree.clear();
    GA_TOKEN_TYPE t = ga_get_token(expr, pos, token_pos, token_length);
    if (t == GA_END) return;
    pos = 0;
    pstring nexpr(new std::string(expr));
    
    t = ga_read_term(nexpr, pos, tree, macro_dict);
    if (tree.root) ga_expand_macro(tree, tree.root, macro_dict);
    
    switch (t) {
    case GA_RPAR:
      ga_throw_error(nexpr, pos-1, "Unbalanced parenthesis.");
      break;
    case GA_RBRACKET:
      ga_throw_error(nexpr, pos-1, "Unbalanced bracket.");
      break;
    case GA_END:
      break;
    default:
      ga_throw_error(nexpr, pos-1, "Unexpected token.");
      break;
    }
  }
  
  // Syntax analysis of a string. Conversion to a tree.
  // Do not register the macros (but expand them).
  void ga_read_string(const std::string &expr, ga_tree &tree,
                      const ga_macro_dictionary &macro_dict) {
    ga_macro_dictionary macro_dict_loc(true, macro_dict);
    ga_read_string_reg(expr, tree, macro_dict_loc);
  }

  // Small tool to make basic substitutions into an assembly string
  // Should be replaced by macros now.
  std::string ga_substitute(const std::string &expr,
                            const std::map<std::string, std::string> &dict) {
    if (dict.size()) {
      size_type pos = 0, token_pos, token_length;
      std::stringstream exprs;

      while (true) {
        GA_TOKEN_TYPE t_type = ga_get_token(expr, pos, token_pos, token_length);
        if (t_type == GA_END) return exprs.str();
        std::string name(&(expr[token_pos]), token_length);
        if (t_type == GA_NAME && dict.find(name) != dict.end())
          exprs << dict.at(name); else exprs << name;
      }
    }
    return expr;
  }

} /* end of namespace */
