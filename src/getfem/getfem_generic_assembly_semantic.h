/*===========================================================================

 Copyright (C) 2013-2018 Yves Renard

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

/** @file   getfem_generic_assembly_semantic.h
    @author Yves Renard <Yves.Renard@insa-lyon.fr>
    @date   November 18, 2013.
    @brief  Semantic analysis of assembly trees and semantic manipulations.
 */


#ifndef GETFEM_GENERIC_ASSEMBLY_SEMANTIC_H__
#define GETFEM_GENERIC_ASSEMBLY_SEMANTIC_H__

#include "getfem/getfem_generic_assembly_tree.h"


namespace getfem {

  /* Performs the semantic analysis, tree simplification and tree enrichment
     - Control tensor sizes for operations, operator or function call
     - Compute all constant operations (i.e. non element dependent)
     - Build a ready to use tree for derivation/compilation
     option = 0 : strict analysis,
              1 : do not complain about incompatible test functions but
                  store them,
              2 : cut incompatible test function branches with respect to the
                  one in workspace.selected_pair
              3 : do not complain about incompatible test functions neither
                  store them.
  */
  void ga_semantic_analysis(ga_tree &tree,
                            const ga_workspace &workspace,
                            const mesh &m,
                            size_type ref_elt_dim,
                            bool eval_fixed_size,
                            bool ignore_X, int option = 0);

  /* Extract the variables used in a sub-tree, ignoring or not the data.
     Variable groups are taken into account. Return if at least one variable
     as been detected or not */
  bool ga_extract_variables(const pga_tree_node pnode,
                            const ga_workspace &workspace,
                            const mesh &m,
                            std::set<var_trans_pair> &vars,
                            bool ignore_data);
  
  /* Extract a sub tree which consists of the corresponding node and of
     the terms multiplying this term, but not the term in addition.
     The aim is to expand an expression is a sum of elementary factors.
     Complains if a nonlinear term is encountered. */
  void ga_extract_factor(ga_tree &result_tree, pga_tree_node pnode,
                         pga_tree_node &new_pnode);

  /* Extract the constant term of degree 1 expressions. */
  bool ga_node_extract_constant_term
  (ga_tree &tree, pga_tree_node pnode, const ga_workspace &workspace,
   const mesh &m);

  /* Extract the Neumann term of an assembly tree with respect to a variable. */
  void ga_extract_Neumann_term
  (ga_tree &tree, const std::string &varname,
   ga_workspace &workspace, pga_tree_node pnode, std::string &result);


  /* Derivation of the tree with respect to a variable.
     The tree is modified and should be copied first and passed to
     ga_semantic_analysis after for enrichment. */
  void ga_derivative(ga_tree &tree, const ga_workspace &workspace,
                     const mesh &m, const std::string &varname,
                     const std::string &interpolatename,
                     size_type order);

  std::string ga_derivative_scalar_function(const std::string &expr,
                                            const std::string &var);

  bool ga_is_affine(const ga_tree &tree, const ga_workspace &workspace,
                    const std::string &varname,
                    const std::string &interpolatename);
  
  // Function of internal use
  inline size_type ref_elt_dim_of_mesh(const mesh &m) {
    return m.convex_index().card() ?
      m.trans_of_convex(m.convex_index().first())->dim() : size_type(0);
  }
  

} /* end of namespace */


#endif /* GETFEM_GENERIC_ASSEMBLY_TREE_H__  */
