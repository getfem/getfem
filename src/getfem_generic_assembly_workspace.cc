/*===========================================================================

 Copyright (C) 2013-2020 Yves Renard

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

#include "getfem/getfem_generic_assembly_tree.h"
#include "getfem/getfem_generic_assembly_semantic.h"
#include "getfem/getfem_generic_assembly_compile_and_exec.h"
#include "getfem/getfem_generic_assembly_functions_and_operators.h"

namespace getfem {

  //=========================================================================
  // Structure dealing with user defined environment : constant, variables,
  // functions, operators.
  //=========================================================================


  void ga_workspace::init() {
    // allocate own storage for K an V to be used unless/until external
    // storage is provided with set_assembled_matrix/vector
    K = std::make_shared<model_real_sparse_matrix>(2,2);
    V = std::make_shared<base_vector>(2);
    KQJpr = std::make_shared<model_real_sparse_matrix>(2,2);
    // add default transformations
    add_interpolate_transformation // deprecated version
      ("neighbour_elt", interpolate_transformation_neighbor_instance());
    add_interpolate_transformation
      ("neighbor_element", interpolate_transformation_neighbor_instance());

    ga_tree tree1;
    pstring s1 = std::make_shared<std::string>("Hess_u");
    tree1.add_name(s1->c_str(), 6, 0, s1);
    tree1.root->name = "u";
    tree1.root->op_type = GA_NAME;
    tree1.root->node_type = GA_NODE_MACRO_PARAM;
    tree1.root->nbc1 = 0;
    tree1.root->nbc2 = ga_parse_prefix_operator(*s1);
    tree1.root->nbc3 = ga_parse_prefix_test(*s1);
    ga_macro gam1("Hess", tree1, 1);
    macro_dict.add_macro(gam1);

    ga_tree tree2;
    pstring s2 = std::make_shared<std::string>("Div_u");
    tree2.add_name(s2->c_str(), 5, 0, s2);
    tree2.root->name = "u";
    tree2.root->op_type = GA_NAME;
    tree2.root->node_type = GA_NODE_MACRO_PARAM;
    tree2.root->nbc1 = 0;
    tree2.root->nbc2 = ga_parse_prefix_operator(*s2);
    tree2.root->nbc3 = ga_parse_prefix_test(*s2);
    ga_macro gam2("Div", tree2, 1);
    macro_dict.add_macro(gam2);
  }

  // variables and variable groups
  void ga_workspace::add_fem_variable
  (const std::string &name, const mesh_fem &mf,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    GMM_ASSERT1(nb_intern_dof == 0 || I.last() < first_intern_dof,
                "The provided interval overlaps with internal dofs");
    nb_prim_dof = std::max(nb_prim_dof, I.last());
    variables.emplace(name, var_description(true, &mf, 0, I, &VV, 1));
  }

  void ga_workspace::add_im_variable
  (const std::string &name, const im_data &imd,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    GMM_ASSERT1(nb_intern_dof == 0 || I.last() <= first_intern_dof,
                "The provided interval overlaps with internal dofs");
    nb_prim_dof = std::max(nb_prim_dof, I.last());
    variables.emplace(name, var_description(true, 0, &imd, I, &VV, 1));
  }

  void ga_workspace::add_internal_im_variable
  (const std::string &name, const im_data &imd,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    GMM_ASSERT1(I.first() >= nb_prim_dof,
                "The provided interval overlaps with primary dofs");
    nb_intern_dof += first_intern_dof - std::min(first_intern_dof, I.first());
    first_intern_dof = std::min(first_intern_dof, I.first());
    nb_intern_dof += first_intern_dof + nb_intern_dof
                   - std::min(first_intern_dof + nb_intern_dof, I.last());
    variables.emplace(name, var_description(true, 0, &imd, I, &VV, 1, true));
  }

  void ga_workspace::add_fixed_size_variable
  (const std::string &name,
   const gmm::sub_interval &I, const model_real_plain_vector &VV) {
    GMM_ASSERT1(nb_intern_dof == 0 || I.last() <= first_intern_dof,
                "The provided interval overlaps with internal dofs");
    nb_prim_dof = std::max(nb_prim_dof, I.last());
    variables.emplace(name, var_description(true, 0, 0, I, &VV,
                                            dim_type(gmm::vect_size(VV))));
  }

  void ga_workspace::add_fem_constant
  (const std::string &name, const mesh_fem &mf,
   const model_real_plain_vector &VV) {
    GMM_ASSERT1(mf.nb_dof(), "The provided mesh_fem of variable" << name
                             << "has zero degrees of freedom");
    size_type Q = gmm::vect_size(VV)/mf.nb_dof();
    if (Q == 0) Q = size_type(1);
    variables.emplace(name, var_description(false, &mf, 0,
                                            gmm::sub_interval(), &VV, Q));
  }

  void ga_workspace::add_fixed_size_constant
  (const std::string &name, const model_real_plain_vector &VV) {
    variables.emplace(name, var_description(false, 0, 0,
                                            gmm::sub_interval(), &VV,
                                            gmm::vect_size(VV)));
  }

  void ga_workspace::add_im_data(const std::string &name, const im_data &imd,
                                 const model_real_plain_vector &VV) {
    variables.emplace(name, var_description
      (false, 0, &imd, gmm::sub_interval(), &VV,
       gmm::vect_size(VV)/(imd.nb_filtered_index() * imd.nb_tensor_elem())));
  }

  bool ga_workspace::is_internal_variable(const std::string &name) const {

    if ((md && md->variable_exists(name) && md->is_internal_variable(name)) ||
        (parent_workspace && parent_workspace->variable_exists(name)
                          && parent_workspace->is_internal_variable(name)))
      return true;
    else {
      VAR_SET::const_iterator it = variables.find(name);
      return it == variables.end() ? false : it->second.is_internal;
    }
  }

  bool ga_workspace::variable_exists(const std::string &name) const {
    return (md && md->variable_exists(name)) ||
      (parent_workspace && parent_workspace->variable_exists(name)) ||
      (variables.find(name) != variables.end());
  }

  bool ga_workspace::variable_group_exists(const std::string &name) const {
    return (variable_groups.find(name) != variable_groups.end()) ||
      (md && md->variable_group_exists(name)) ||
      (parent_workspace && parent_workspace->variable_group_exists(name));
  }

  const std::vector<std::string>&
  ga_workspace::variable_group(const std::string &group_name) const {
    std::map<std::string, std::vector<std::string> >::const_iterator
      it = variable_groups.find(group_name);
    if (it != variable_groups.end())
      return (variable_groups.find(group_name))->second;
    if (md && md->variable_group_exists(group_name))
      return md->variable_group(group_name);
    if (parent_workspace &&
        parent_workspace->variable_group_exists(group_name))
      return parent_workspace->variable_group(group_name);
    GMM_ASSERT1(false, "Undefined variable group " << group_name);
  }

  const std::string&
  ga_workspace::first_variable_of_group(const std::string &name) const {
    const std::vector<std::string> &t = variable_group(name);
    GMM_ASSERT1(t.size(), "Variable group " << name << " is empty");
    return t[0];
  }

  bool ga_workspace::is_constant(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end())
      return !(it->second.is_variable);
    else if (variable_group_exists(name))
      return is_constant(first_variable_of_group(name));
    else if (reenabled_var_intervals.count(name))
      return false;
    else if (md && md->variable_exists(name))
      return md->is_data(name);
    else if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->is_constant(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  bool ga_workspace::is_disabled_variable(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end())
      return false;
    else if (variable_group_exists(name))
      return is_disabled_variable(first_variable_of_group(name));
    else if (reenabled_var_intervals.count(name))
      return false;
    else if (md && md->variable_exists(name))
      return md->is_disabled_variable(name);
    else if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->is_disabled_variable(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  const scalar_type &
  ga_workspace::factor_of_variable(const std::string &name) const {
    static const scalar_type one(1);
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return one;
    if (variable_group_exists(name))
      return one;
    if (md && md->variable_exists(name)) return md->factor_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->factor_of_variable(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  const gmm::sub_interval &
  ga_workspace::interval_of_variable(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return it->second.I;
    const auto it2 = reenabled_var_intervals.find(name);
    if (it2 != reenabled_var_intervals.end()) return it2->second;
    if (with_parent_variables && md && md->variable_exists(name))
      return md->interval_of_variable(name);
    else if (with_parent_variables &&
             parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->interval_of_variable(name);
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  const mesh_fem *
  ga_workspace::associated_mf(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end())
      return it->second.is_fem_dofs ? it->second.mf : 0;
    if (md && md->variable_exists(name))
      return md->pmesh_fem_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->associated_mf(name);
    if (variable_group_exists(name))
      return associated_mf(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  const im_data *
  ga_workspace::associated_im_data(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) return it->second.imd;
    if (md && md->variable_exists(name))
      return md->pim_data_of_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->associated_im_data(name);
    if (variable_group_exists(name)) return 0;
    GMM_ASSERT1(false, "Undefined variable " << name);
  }

  size_type ga_workspace::qdim(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) {
      const mesh_fem *mf = it->second.is_fem_dofs ? it->second.mf : 0;
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

  bgeot::multi_index
  ga_workspace::qdims(const std::string &name) const {
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

  const model_real_plain_vector &
  ga_workspace::value(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end())
      return *(it->second.V);
    if (md && md->variable_exists(name))
      return md->real_variable(name);
    if (parent_workspace && parent_workspace->variable_exists(name))
      return parent_workspace->value(name);
    if (variable_group_exists(name))
      return value(first_variable_of_group(name));
    GMM_ASSERT1(false, "Undefined variable or group " << name);
  }

  scalar_type ga_workspace::get_time_step() const {
    if (md) return md->get_time_step();
    if (parent_workspace) return parent_workspace->get_time_step();
    GMM_ASSERT1(false, "No time step defined here");
  }

  void ga_workspace::add_interpolate_transformation
  (const std::string &name, pinterpolate_transformation ptrans) {
     if (secondary_domain_exists(name))
      GMM_ASSERT1(false, "A secondary domain with the same "
                  "name already exists");
    if (transformations.find(name) != transformations.end())
      GMM_ASSERT1(name != "neighbor_element", "neighbor_element is a "
                  "reserved interpolate transformation name");
    transformations[name] = ptrans;
  }

  bool ga_workspace::interpolate_transformation_exists
  (const std::string &name) const {
    return (md && md->interpolate_transformation_exists(name)) ||
      (parent_workspace &&
       parent_workspace->interpolate_transformation_exists(name)) ||
      (transformations.find(name) != transformations.end());
  }

  pinterpolate_transformation
  ga_workspace::interpolate_transformation(const std::string &name) const {
    auto it = transformations.find(name);
    if (it != transformations.end()) return it->second;
    if (md && md->interpolate_transformation_exists(name))
      return md->interpolate_transformation(name);
    if (parent_workspace &&
       parent_workspace->interpolate_transformation_exists(name))
      return parent_workspace->interpolate_transformation(name);
    GMM_ASSERT1(false, "Inexistent transformation " << name);
  }

  bool ga_workspace::elementary_transformation_exists
  (const std::string &name) const {
    return (md && md->elementary_transformation_exists(name)) ||
      (parent_workspace &&
       parent_workspace->elementary_transformation_exists(name)) ||
      (elem_transformations.find(name) != elem_transformations.end());
  }

  pelementary_transformation
  ga_workspace::elementary_transformation(const std::string &name) const {
    auto it = elem_transformations.find(name);
    if (it != elem_transformations.end()) return it->second;
    if (md && md->elementary_transformation_exists(name))
      return md->elementary_transformation(name);
    if (parent_workspace &&
       parent_workspace->elementary_transformation_exists(name))
      return parent_workspace->elementary_transformation(name);
    GMM_ASSERT1(false, "Inexistent elementary transformation " << name);
  }

  void ga_workspace::add_secondary_domain(const std::string &name,
                                          psecondary_domain psecdom) {
    if (interpolate_transformation_exists(name))
      GMM_ASSERT1(false, "An interpolate transformation with the same "
                  "name already exists");
    secondary_domains[name] = psecdom;
  }

  bool ga_workspace::secondary_domain_exists
  (const std::string &name) const {
    return (md && md->secondary_domain_exists(name)) ||
      (parent_workspace &&
       parent_workspace->secondary_domain_exists(name)) ||
      (secondary_domains.find(name) != secondary_domains.end());
  }

  psecondary_domain
  ga_workspace::secondary_domain(const std::string &name) const {
    auto it = secondary_domains.find(name);
    if (it != secondary_domains.end()) return it->second;
    if (md && md->secondary_domain_exists(name))
      return md->secondary_domain(name);
    if (parent_workspace &&
       parent_workspace->secondary_domain_exists(name))
      return parent_workspace->secondary_domain(name);
    GMM_ASSERT1(false, "Inexistent secondary domain " << name);
  }


  const mesh_region &
  ga_workspace::register_region(const mesh &m, const mesh_region &region) {
    if (&m == &dummy_mesh())
      return dummy_mesh_region();

    std::list<mesh_region> &lmr = registred_mesh_regions[&m];
    for (const mesh_region &rg : lmr)
      if (rg.compare(m, region, m)) return rg;
    lmr.push_back(region);
    return lmr.back();
  }

  void ga_workspace::add_tree(ga_tree &tree, const mesh &m,
                              const mesh_im &mim, const mesh_region &rg,
                              const std::string &expr,
                              size_type add_derivative_order,
                              bool function_expr, operation_type op_type,
                              const std::string varname_interpolation) {
    if (tree.root) {
      // cout << "add tree with tests functions of " <<  tree.root->name_test1
      //     << " and " << tree.root->name_test2 << endl;
      //     ga_print_node(tree.root, cout); cout << endl;

      // Eliminate the term if it corresponds to disabled variables
      if ((tree.root->test_function_type >= 1 &&
           is_disabled_variable(tree.root->name_test1)) ||
          (tree.root->test_function_type >= 2 &&
           is_disabled_variable(tree.root->name_test2))) {
        // cout<<"disabling term ";  ga_print_node(tree.root, cout); cout<<endl;
        return;
      }

      bool remain = true;
      size_type order = 0, ind_tree = 0;

      if (op_type != ga_workspace::ASSEMBLY)
        order = add_derivative_order;
      else {
        switch(tree.root->test_function_type) {
        case 0: order = 0; break;
        case 1: order = 1; break;
        case 3: order = 2; break;
        default: GMM_ASSERT1(false, "Inconsistent term "
                             << tree.root->test_function_type);
        }
      }

      bool found = false;
      for (const ga_workspace::tree_description &td : trees) {
        if (td.mim == &mim &&
            td.m == &m &&
            td.secondary_domain == tree.secondary_domain &&
            td.order == order &&
            td.name_test1 == tree.root->name_test1 &&
            td.interpolate_name_test1 == tree.root->interpolate_name_test1 &&
            td.name_test2 == tree.root->name_test2 &&
            td.interpolate_name_test2 == tree.root->interpolate_name_test2 &&
            td.rg == &rg &&
            td.operation == op_type &&
            td.varname_interpolation == varname_interpolation) {
          ga_tree &ftree = *(td.ptree);

          ftree.insert_node(ftree.root, GA_NODE_OP);
          ftree.root->op_type = GA_PLUS;
          ftree.root->children.resize(2, nullptr);
          ftree.copy_node(tree.root, ftree.root, ftree.root->children[1]);
          ga_semantic_analysis(ftree, *this, m,
                               ref_elt_dim_of_mesh(m), false, function_expr);
          found = true;
          break;
        }
      }

      if (!found) {
        ind_tree = trees.size();
        remain = false;
        trees.push_back(tree_description());
        trees.back().mim = &mim;
        trees.back().m = &m;
        trees.back().rg = &rg;
        trees.back().secondary_domain = tree.secondary_domain;
        trees.back().ptree = new ga_tree;
        trees.back().ptree->swap(tree);
        pga_tree_node root = trees.back().ptree->root;
        trees.back().name_test1 = root->name_test1;
        trees.back().name_test2 = root->name_test2;
        trees.back().interpolate_name_test1 = root->interpolate_name_test1;
        trees.back().interpolate_name_test2 = root->interpolate_name_test2;
        trees.back().order = order;
        trees.back().operation = op_type;
        trees.back().varname_interpolation = varname_interpolation;
       }

      if (op_type == ga_workspace::ASSEMBLY && order < add_derivative_order) {
        std::set<var_trans_pair> expr_variables;
        ga_extract_variables((remain ? tree : *(trees[ind_tree].ptree)).root,
                             *this, m, expr_variables, true);
        for (const var_trans_pair &var : expr_variables) {
          if (!(is_constant(var.varname))) {
            ga_tree dtree = (remain ? tree : *(trees[ind_tree].ptree));
            // cout << "Derivation with respect to " << var.varname << " : "
            //   << var.transname << " of " << ga_tree_to_string(dtree) << endl;
            // GA_TIC;
            ga_derivative(dtree, *this, m, var.varname, var.transname, 1+order);
            // cout << "Result : " << ga_tree_to_string(dtree) << endl;
            // GA_TOCTIC("Derivative time");
            ga_semantic_analysis(dtree, *this, m,
                                 ref_elt_dim_of_mesh(m), false, function_expr);
            // GA_TOCTIC("Analysis after Derivative time");
            // cout << "after analysis "  << ga_tree_to_string(dtree) << endl;
            add_tree(dtree, m, mim, rg, expr, add_derivative_order,
                     function_expr, op_type, varname_interpolation);
          }
        }
      }
    }
  }

  size_type ga_workspace::add_expression(const std::string &expr,
                                         const mesh_im &mim,
                                         const mesh_region &rg_,
                                         size_type add_derivative_order,
                                         const std::string &secondary_dom) {
    const mesh_region &rg = register_region(mim.linked_mesh(), rg_);
    // cout << "adding expression " << expr << endl;
    GA_TIC;
    size_type max_order = 0;
    std::vector<ga_tree> ltrees(1);
    ga_read_string(expr, ltrees[0], macro_dictionary());
    if (secondary_dom.size()) {
      GMM_ASSERT1(secondary_domain_exists(secondary_dom),
                  "Unknown secondary domain " << secondary_dom);
      ltrees[0].secondary_domain = secondary_dom;
    }
    // cout << "read : " << ga_tree_to_string(ltrees[0])  << endl;
    ga_semantic_analysis(ltrees[0], *this, mim.linked_mesh(),
                         ref_elt_dim_of_mesh(mim.linked_mesh()),
                         false, false, 1);
    // cout << "analysed : " << ga_tree_to_string(ltrees[0]) << endl;
    GA_TOC("First analysis time");
    if (ltrees[0].root) {
      if (test1.size() > 1 || test2.size() > 1) {
        size_type ntest2 = test2.size();
        if (ntest2 == 0)                  // temporarily add an element to
          test2.insert(var_trans_pair()); // allow entering the inner loop
        ltrees.resize(test1.size()*test2.size(), ltrees[0]);
        auto ltree = ltrees.begin();
        for (const auto &t1 : test1) {
          for (const auto &t2 : test2) {
            selected_test1 = t1;
            if (ntest2 > 0) selected_test2 = t2;
            // cout << "analysis with " << selected_test1.first << endl;
            ga_semantic_analysis(*ltree, *this, mim.linked_mesh(),
                                 ref_elt_dim_of_mesh(mim.linked_mesh()),
                                 false, false, 2);
            // cout <<"split: "<< ga_tree_to_string(*ltree) << endl;
            if (ltree != ltrees.end()) ++ltree;
          }
        }
        if (ntest2 == 0) test2.clear(); // remove temporarily added element
      }

      for (ga_tree &ltree : ltrees) {
        if (ltree.root) {
          // cout << "adding tree " << ga_tree_to_string(ltree) << endl;
          max_order = std::max(ltree.root->nb_test_functions(), max_order);
          add_tree(ltree, mim.linked_mesh(), mim, rg, expr,
                   add_derivative_order, true);
        }
      }
    }
    GA_TOC("Time for add expression");
    return max_order;
  }

  void ga_workspace::add_function_expression(const std::string &expr) {
    ga_tree tree;
    ga_read_string(expr, tree, macro_dictionary());
    ga_semantic_analysis(tree, *this, dummy_mesh(), 1, false, true);
    if (tree.root) {
      // GMM_ASSERT1(tree.root->nb_test_functions() == 0,
      //            "Invalid function expression");
      add_tree(tree, dummy_mesh(), dummy_mesh_im(), dummy_mesh_region(),
               expr, 0, true);
    }
  }

  void ga_workspace::add_interpolation_expression(const std::string &expr,
                                                  const mesh &m,
                                                  const mesh_region &rg_) {
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree, macro_dictionary());
    ga_semantic_analysis(tree, *this, m, ref_elt_dim_of_mesh(m),
                         false, false);
    if (tree.root) {
      // GMM_ASSERT1(tree.root->nb_test_functions() == 0,
      //            "Invalid expression containing test functions");
      add_tree(tree, m, dummy_mesh_im(), rg, expr, 0, false,
               ga_workspace::PRE_ASSIGNMENT);
    }
  }

  void ga_workspace::add_interpolation_expression(const std::string &expr,
                                                  const mesh_im &mim,
                                                  const mesh_region &rg_) {
    const mesh &m = mim.linked_mesh();
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree, macro_dictionary());
    ga_semantic_analysis(tree, *this, m, ref_elt_dim_of_mesh(m),
                         false, false);
    if (tree.root) {
      GMM_ASSERT1(tree.root->nb_test_functions() == 0,
                  "Invalid expression containing test functions");
      add_tree(tree, m, mim, rg, expr, 0, false,
               ga_workspace::PRE_ASSIGNMENT);
    }
  }

  void ga_workspace::add_assignment_expression
  (const std::string &varname, const std::string &expr, const mesh_region &rg_,
   size_type order, bool before) {
    const im_data *imd = associated_im_data(varname);
    GMM_ASSERT1(imd != 0, "Only applicable to im_data");
    const mesh_im &mim = imd->linked_mesh_im();
    const mesh &m = mim.linked_mesh();
    const mesh_region &rg = register_region(m, rg_);
    ga_tree tree;
    ga_read_string(expr, tree, macro_dictionary());
    ga_semantic_analysis(tree, *this, m, ref_elt_dim_of_mesh(m), false, false);
    if (tree.root) {
      GMM_ASSERT1(tree.root->nb_test_functions() == 0,
                  "Invalid expression containing test functions");
      add_tree(tree, m, mim, rg, expr, order, false,
               before ? ga_workspace::PRE_ASSIGNMENT
                      : ga_workspace::POST_ASSIGNMENT,
               varname);
    }
  }

  size_type ga_workspace::nb_trees() const { return trees.size(); }

  ga_workspace::tree_description &ga_workspace::tree_info(size_type i)
  { return trees[i]; }

  bool ga_workspace::used_variables(std::vector<std::string> &vl,
                                    std::vector<std::string> &vl_test1,
                                    std::vector<std::string> &vl_test2,
                                    std::vector<std::string> &dl,
                                    size_type order) {
    bool islin = true;
    std::set<var_trans_pair> vll, dll;
    for (const std::string &v : vl) vll.insert(var_trans_pair(v, ""));
    for (const std::string &d : dl) dll.insert(var_trans_pair(d, ""));

    for (const ga_workspace::tree_description &td : trees) {
      std::set<var_trans_pair> dllaux;
      bool fv = ga_extract_variables(td.ptree->root, *this, *(td.m),
                                     dllaux, false);

      if (td.order == order)
        for (const auto &t : dllaux)
          dll.insert(t);

      switch (td.order) {
      case 0:  break;
      case 1:
        if (td.order == order) {
          if (variable_group_exists(td.name_test1)) {
            for (const std::string &t : variable_group(td.name_test1))
              vll.insert(var_trans_pair(t, td.interpolate_name_test1));
          } else {
            vll.insert(var_trans_pair(td.name_test1,
                                      td.interpolate_name_test1));
          }
          bool found = false;
          for (const std::string &t : vl_test1)
            if (td.name_test1 == t)
              found = true;
          if (!found)
            vl_test1.push_back(td.name_test1);
        }
        break;
      case 2:
        if (td.order == order) {
          if (variable_group_exists(td.name_test1)) {
            for (const std::string &t : variable_group(td.name_test1))
              vll.insert(var_trans_pair(t, td.interpolate_name_test1));
          } else {
            vll.insert(var_trans_pair(td.name_test1,
                                      td.interpolate_name_test1));
          }
          if (variable_group_exists(td.name_test2)) {
            for (const std::string &t : variable_group(td.name_test2))
              vll.insert(var_trans_pair(t, td.interpolate_name_test2));
          } else {
            vll.insert(var_trans_pair(td.name_test2,
                                      td.interpolate_name_test2));
          }
          bool found = false;
          for (size_type j = 0; j < vl_test1.size(); ++j)
            if ((td.name_test1 == vl_test1[j]) &&
                (td.name_test2 == vl_test2[j]))
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
    for (const auto &var : vll)
      if (vl.size() == 0 || var.varname != vl.back())
        vl.push_back(var.varname);
    dl.clear();
    for (const auto &var : dll)
      if (dl.size() == 0 || var.varname != dl.back())
        dl.push_back(var.varname);

    return islin;
  }

  bool ga_workspace::is_linear(size_type order) {
    std::vector<std::string> vl, vl_test1, vl_test2, dl;
    return used_variables(vl, vl_test1, vl_test2, dl, order);
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
      for (const std::string &t : variable_group(group_name))
        if (&(associated_mf(t)->linked_mesh()) == &m)
          return t;
      GMM_ASSERT1(false, "No variable in this group for the given mesh");
    } else
      return group_name;
  }


  void ga_workspace::assembly(size_type order, bool condensation) {

    const ga_workspace *w = this;
    while (w->parent_workspace) w = w->parent_workspace;
    if (w->md) w->md->nb_dof(); // To eventually call actualize_sizes()

    GA_TIC;
    ga_instruction_set gis;
    ga_compile(*this, gis, order, condensation);
    GA_TOCTIC("Compile time");

    size_type nb_tot_dof = condensation ? nb_prim_dof + nb_intern_dof
                                        : nb_prim_dof;
    if (order == 2) {
      if (K.use_count()) {
        gmm::clear(*K);
        gmm::resize(*K, nb_prim_dof, nb_prim_dof);
      } // else
      // We trust that the caller has provided a matrix large enough for the
      // terms to be assembled (can actually be smaller than the full matrix)
      //  GMM_ASSERT1(gmm::mat_nrows(*K) == nb_prim_dof &&
      //              gmm::mat_ncols(*K) == nb_prim_dof, "Wrong sizes");
      if (KQJpr.use_count()) {
        gmm::clear(*KQJpr);
        gmm::resize(*KQJpr, nb_prim_dof+nb_intern_dof, nb_prim_dof); // redundant if condensation == false
      } else if (condensation)
        GMM_ASSERT1(gmm::mat_nrows(*KQJpr) == nb_prim_dof+nb_intern_dof &&
                    gmm::mat_ncols(*KQJpr) == nb_prim_dof, "Wrong sizes");
      gmm::clear(col_unreduced_K);
      gmm::clear(row_unreduced_K);
      gmm::clear(row_col_unreduced_K);
      gmm::resize(col_unreduced_K, nb_tot_dof, nb_tmp_dof);
      gmm::resize(row_unreduced_K, nb_tmp_dof, nb_tot_dof);
      gmm::resize(row_col_unreduced_K, nb_tmp_dof, nb_tmp_dof);
      if (condensation) {
        gmm::clear(unreduced_V);
        gmm::resize(unreduced_V, nb_tmp_dof);
      }
    } else if (order == 1) {
      if (V.use_count()) {
        gmm::clear(*V);
        gmm::resize(*V, nb_tot_dof);
      } else
        GMM_ASSERT1(V->size() == nb_tot_dof, "Wrong size");
      gmm::clear(unreduced_V);
      gmm::resize(unreduced_V, nb_tmp_dof);
    }
    gmm::clear(assembled_tensor().as_vector());

    GA_TOCTIC("Init time");
    ga_exec(gis, *this);     // --> unreduced_V, *V,
    GA_TOCTIC("Exec time");  //     unreduced_K, *K

    if (order == 0) {
      MPI_SUM_VECTOR(assemb_t.as_vector());
    } else if (order == 1 || (order == 2 && condensation)) {
      MPI_SUM_VECTOR(*V);
      MPI_SUM_VECTOR(unreduced_V);
    }

    // Deal with reduced fems, unreduced_K --> *K, *KQJpr,
    //                         unreduced_V --> *V
    if (order > 0) {
      std::set<std::string> vars_vec_done;
      std::set<std::pair<std::string, std::string> > vars_mat_done;
      for (const auto &term : gis.unreduced_terms) {
        const std::string &name1 = term.first;
        const std::string &name2 = term.second;
        const std::vector<std::string>
          vg1_(1,name1), vg2_(1,name2),
          &vg1 = variable_group_exists(name1) ? variable_group(name1) : vg1_,
          &vg2 = variable_group_exists(name2) ? variable_group(name2) : vg2_;
        if (order == 1) {
          for (const std::string &vname1 : vg1) {
            const mesh_fem *mf1 = associated_mf(vname1);
            if (mf1 && mf1->is_reduced() && vars_vec_done.count(vname1) == 0) {
              gmm::sub_interval uI1 = temporary_interval_of_variable(vname1),
                                I1 = interval_of_variable(vname1);
              gmm::mult_add(gmm::transposed(mf1->extension_matrix()),
                            gmm::sub_vector(unreduced_V, uI1),
                            gmm::sub_vector(*V, I1));
              vars_vec_done.insert(vname1);
            }
          }
        } else {
          for (const std::string &vname1 : vg1) {
            for (const std::string &vname2 : vg2) {
              const mesh_fem *mf1 = associated_mf(vname1),
                             *mf2 = associated_mf(vname2);
              if (((mf1 && mf1->is_reduced()) || (mf2 && mf2->is_reduced()))
                  && vars_mat_done.count(std::make_pair(vname1,vname2)) == 0) {
                gmm::sub_interval
                  uI1 = temporary_interval_of_variable(vname1),
                  uI2 = temporary_interval_of_variable(vname2),
                  I1 = interval_of_variable(vname1),
                  I2 = interval_of_variable(vname2);
                if (mf1 && mf1->is_reduced() && mf2 && mf2->is_reduced()) {
                  model_real_sparse_matrix aux(I1.size(), uI2.size());
                  model_real_row_sparse_matrix M(I1.size(), I2.size());
                  gmm::mult(gmm::transposed(mf1->extension_matrix()),
                            gmm::sub_matrix(row_col_unreduced_K, uI1, uI2),
                            aux);
                  gmm::mult(aux, mf2->extension_matrix(), M);
                  gmm::add(M, gmm::sub_matrix(*K, I1, I2));
                } else if (mf1 && mf1->is_reduced()) {
                  if (condensation && vars_vec_done.count(vname1) == 0) {
                    gmm::mult_add(gmm::transposed(mf1->extension_matrix()),
                                  gmm::sub_vector(unreduced_V, uI1),
                                  gmm::sub_vector(*V, I1));
                    vars_vec_done.insert(vname1);
                  }
                  model_real_sparse_matrix M(I1.size(), I2.size());
                  gmm::mult(gmm::transposed(mf1->extension_matrix()),
                            gmm::sub_matrix(row_unreduced_K, uI1, I2), M);
                  gmm::add(M, gmm::sub_matrix(*K, I1, I2));
                } else {
                  model_real_row_sparse_matrix M(I1.size(), I2.size());
                  gmm::mult(gmm::sub_matrix(col_unreduced_K, I1, uI2),
                            mf2->extension_matrix(), M);
                  if (I1.first() < nb_prim_dof) {
                    GMM_ASSERT1(I1.last() <= nb_prim_dof, "Internal error");
                    gmm::add(M, gmm::sub_matrix(*K, I1, I2)); // -> *K
                  } else { // vname1 is an internal variable
                    gmm::add(M, gmm::sub_matrix(*KQJpr, I1, I2)); // -> *KQJpr
                  }
                }
                vars_mat_done.insert(std::make_pair(vname1,vname2));
              }
            }
          }
        }
      }
    }
  }

  void ga_workspace::set_include_empty_int_points(bool include) {
    include_empty_int_pts = include;
  }

  bool ga_workspace::include_empty_int_points() const {
    return include_empty_int_pts;
  }

  void ga_workspace::add_temporary_interval_for_unreduced_variable
    (const std::string &name)
  {
    if (variable_group_exists(name)) {
      for (const std::string &v : variable_group(name))
        add_temporary_interval_for_unreduced_variable(v);
    } else if (tmp_var_intervals.count(name) == 0) {
      const mesh_fem *mf = associated_mf(name);
      if (mf && mf->is_reduced()) {
        size_type nd = mf->nb_basic_dof();
        tmp_var_intervals[name] = gmm::sub_interval(nb_tmp_dof, nd);
        nb_tmp_dof += nd;
      }
    }
  }

  void ga_workspace::clear_expressions() { trees.clear(); }

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
    operation = td.operation;
    varname_interpolation = td.varname_interpolation;
    name_test1 = td.name_test1;
    name_test2 = td.name_test2;
    interpolate_name_test1 = td.interpolate_name_test1;
    interpolate_name_test2 = td.interpolate_name_test2;
    mim = td.mim;
    m = td.m;
    rg = td.rg;
    ptree = 0;
    if (td.ptree) ptree = new ga_tree(*(td.ptree));
  }

  ga_workspace::tree_description &ga_workspace::tree_description::operator =
  (const ga_workspace::tree_description& td)
  { if (ptree) delete ptree; ptree = 0; copy(td); return *this; }
  ga_workspace::tree_description::~tree_description()
  { if (ptree) delete ptree; ptree = 0; }

  ga_workspace::ga_workspace(const getfem::model &md_,
                             const inherit var_inherit)
    : md(&md_), parent_workspace(0),
      with_parent_variables(var_inherit == inherit::ENABLED ||
                            var_inherit == inherit::ALL),
      nb_tmp_dof(0), macro_dict(md_.macro_dictionary())
  {
    init();
    nb_prim_dof = with_parent_variables ? md->nb_primary_dof() : 0;
    nb_intern_dof = with_parent_variables ? md->nb_internal_dof() : 0;
    if (var_inherit == inherit::ALL) { // enable model's disabled variables
      model::varnamelist vlmd;
      md->variable_list(vlmd);
      for (const auto &varname : vlmd)
        if (md->is_disabled_variable(varname)) {
          if (md->is_affine_dependent_variable(varname)) {
            std::string orgvarname = md->org_variable(varname);
            if (reenabled_var_intervals.count(orgvarname) == 0) {
              size_type varsize = gmm::vect_size(md->real_variable(orgvarname));
              reenabled_var_intervals[orgvarname]
                = gmm::sub_interval (nb_prim_dof, varsize);
              nb_prim_dof += varsize;
            }
            reenabled_var_intervals[varname]
              = reenabled_var_intervals[orgvarname];
          } else  if (reenabled_var_intervals.count(varname) == 0) {
            size_type varsize = gmm::vect_size(md->real_variable(varname));
            reenabled_var_intervals[varname]
              = gmm::sub_interval(nb_prim_dof, varsize);
            nb_prim_dof += varsize;
          }
        }
    }
    first_intern_dof = nb_prim_dof; // dofs are contiguous in getfem::model
  }
  ga_workspace::ga_workspace(const ga_workspace &gaw,
                             const inherit var_inherit)
    : md(0), parent_workspace(&gaw),
      with_parent_variables(var_inherit == inherit::ENABLED ||
                            var_inherit == inherit::ALL),
      nb_tmp_dof(0), macro_dict(gaw.macro_dictionary())
  {
    init();
    nb_prim_dof = with_parent_variables ? gaw.nb_primary_dof() : 0;
    nb_intern_dof = with_parent_variables ? gaw.nb_internal_dof() : 0;
    first_intern_dof = with_parent_variables ? gaw.first_internal_dof() : 0;
  }
  ga_workspace::ga_workspace()
    : md(0), parent_workspace(0), with_parent_variables(false),
      nb_prim_dof(0), nb_intern_dof(0), first_intern_dof(0), nb_tmp_dof(0)
  { init(); }
  ga_workspace::~ga_workspace() { clear_expressions(); }

  //=========================================================================
  // Extract the constant term of degree 1 expressions
  //=========================================================================

  std::string ga_workspace::extract_constant_term(const mesh &m) {
    std::string constant_term;
    for (const ga_workspace::tree_description &td : trees) {
      if (td.order == 1) {
        ga_tree local_tree = *(td.ptree);
        if (local_tree.root)
          ga_node_extract_constant_term(local_tree, local_tree.root, *this, m);
        if (local_tree.root)
          ga_semantic_analysis(local_tree, *this, m,
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

  std::string ga_workspace::extract_order0_term() {
    std::string term;
    for (const ga_workspace::tree_description &td : trees) {
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
    for (const ga_workspace::tree_description &td : trees) {
      if (td.order == 1 && td.name_test1 == varname) {
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

  std::string ga_workspace::extract_Neumann_term(const std::string &varname) {
    std::string result;
    for (const ga_workspace::tree_description &td : trees) {
      if (td.order == 1 && td.name_test1 == varname) {
        ga_tree &local_tree = *(td.ptree);
        if (local_tree.root)
          ga_extract_Neumann_term(local_tree, varname, *this,
                                  local_tree.root, result);
      }
    }
    return result;
  }

} /* end of namespace */
