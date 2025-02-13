/*===========================================================================

 Copyright (C) 2009-2020 Yves Renard

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

#include <iomanip>
#include "gmm/gmm_range_basis.h"
#include "gmm/gmm_solver_cg.h"
#include "gmm/gmm_condition_number.h"
#include "getfem/getfem_models.h"
#include "getfem/getfem_accumulated_distro.h"
#include "getfem/getfem_assembling.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_generic_assembly_tree.h"


namespace getfem {

  model::model(bool comp_version) {
    init(); complex_version = comp_version;
    is_linear_ = is_symmetric_ = is_coercive_ = true;
    leading_dim = 0;
    time_integration = 0; init_step = false; time_step = scalar_type(1);
    add_interpolate_transformation
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

  void model::var_description::set_size() {
    clear_temporaries();
    v_num_var_iter.resize(n_iter);
    v_num_iter.resize(n_iter);
    size_type s = mf ? passociated_mf()->nb_dof()
                     : (imd ? imd->nb_filtered_index()
                              *imd->nb_tensor_elem()
                            : 1);
    s *= qdim();
    for (size_type i = 0; i < n_iter; ++i)
      if (is_complex)
        complex_value[i].resize(s);
      else
        real_value[i].resize(s);
    if (is_affine_dependent) {
      if (is_complex)
        affine_complex_value.resize(s);
      else
        affine_real_value.resize(s);
    }
  }

  size_type model::var_description::add_temporary(gmm::uint64_type id_num) {
    size_type nit = n_iter;
    for (; nit < n_iter + n_temp_iter ; ++nit)
      if (v_num_var_iter[nit] == id_num) break;
    if (nit >=  n_iter + n_temp_iter) {
      ++n_temp_iter;
      v_num_var_iter.resize(nit+1);
      v_num_var_iter[nit] = id_num;
      v_num_iter.resize(nit+1);
      v_num_iter[nit] = 0;
      if (is_complex) {
        size_type s = complex_value[0].size();
        complex_value.resize(n_iter + n_temp_iter);
        complex_value[nit].resize(s);
      } else {
        size_type s = real_value[0].size();
        real_value.resize(n_iter + n_temp_iter);
        real_value[nit].resize(s);
      }
    }
    return nit;
  }

  void model::var_description::clear_temporaries() {
    n_temp_iter = 0;
    default_iter = 0;
    if (is_complex)
      complex_value.resize(n_iter);
    else
      real_value.resize(n_iter);
  }

  bool model::check_name_validity(const std::string &name, bool assert) const {

    if (variables.count(name) != 0) {
      GMM_ASSERT1(!assert, "Variable " << name << " already exists");
      return false;
    } else if (variable_groups.count(name) != 0) {
      GMM_ASSERT1(!assert,
                  name << " corresponds to an already existing group of "
                  "variables name");
      return false;
    } else if (macro_exists(name)) {
      GMM_ASSERT1(!assert,
                  name << " corresponds to an already existing macro");
      return false;
    } else if (name.compare("X") == 0) {
      GMM_ASSERT1(!assert, "X is a reserved keyword of the generic "
                  "assembly language");
      return false;
    }

    int ga_valid = ga_check_name_validity(name);
    if (ga_valid == 1) {
      GMM_ASSERT1(!assert, "Invalid variable name, corresponds to an "
                "operator or function name of the generic assembly language");
      return false;
    } else if (ga_valid == 2) {
      GMM_ASSERT1(!assert, "Invalid variable name having a reserved "
                  "prefix used by the generic assembly language");
      return false;
    } else if (ga_valid == 3) {
      std::string org_name = sup_previous_and_dot_to_varname(name);
      if (org_name.size() < name.size() &&
          variables.find(org_name) != variables.end()) {
        GMM_ASSERT1(!assert,
                    "Dot and Previous are reserved prefix used for time "
                    "integration schemes");
        return false;
      }
    }

    bool valid = !name.empty() && isalpha(name[0]);
    if (valid)
      for (size_type i = 1; i < name.size(); ++i)
        if (!(isalnum(name[i]) || name[i] == '_')) valid = false;
    GMM_ASSERT1(!assert || valid,
                "Illegal variable name : \"" << name << "\"");
    return valid;
  }

  std::string model::new_name(const std::string &name) {
    std::string res_name = name;
    bool valid = check_name_validity(res_name, false);
    for (size_type i = 2; !valid && i < 50; ++i) {
      std::stringstream m;
      m << name << '_' << i;
      res_name = m.str();
      valid = check_name_validity(res_name, false);
    }
    for (size_type i = 2; !valid && i < 1000; ++i) {
      std::stringstream m;
      m << "new_" << name << '_' << i;
      res_name = m.str();
      valid = check_name_validity(res_name, false);
    }
    GMM_ASSERT1(valid, "Illegal variable name: " << name);
    return res_name;
  }


  model::VAR_SET::const_iterator
  model::find_variable(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
    return it;
  }

  const model::var_description &
  model::variable_description(const std::string &name) const {
    return find_variable(name)->second;
  }

  std::string sup_previous_and_dot_to_varname(std::string v) {
    if (!(v.compare(0, 8, "Previous")) && (v[8] == '_' || v[9] == '_')) {
      v = v.substr((v[8] == '_') ? 9 : 10);
    }
    if (!(v.compare(0, 3, "Dot")) && (v[3] == '_' || v[4] == '_')) {
      v = v.substr((v[3] == '_') ? 4 : 5);
    }
    if (is_old(v)) v = no_old_prefix_name(v);
    return v;
  }

  bool is_old(const std::string &name){
    return name.substr(0, PREFIX_OLD_LENGTH) == PREFIX_OLD;
  }

  std::string no_old_prefix_name(const std::string &name){
    return is_old(name) ? name.substr(PREFIX_OLD_LENGTH) : name;
  }

  bool model::is_disabled_variable(const std::string &name) const {
    if (is_old(name)) return false;
    VAR_SET::const_iterator it = find_variable(name);
    if (!(it->second.is_variable)) return false;
    if (it->second.is_affine_dependent)
      it = variables.find(it->second.org_name);
    return it->second.is_disabled;
  }

  bool model::is_data(const std::string &name) const {
    if (is_old(name)) return true;
    VAR_SET::const_iterator it = find_variable(name);
    if (it->second.is_affine_dependent)
      it = variables.find(it->second.org_name);
    return !(it->second.is_variable) || it->second.is_disabled;
  }

  bool model::is_true_data(const std::string &name) const {
    return is_old(name) || !(variable_description(name).is_variable);
  }

  bool model::is_internal_variable(const std::string &name) const {
    if (is_old(name)) return false;
    const auto &var_descr = variable_description(name);
    return var_descr.is_internal && var_descr.is_enabled();
  }

  bool model::is_affine_dependent_variable(const std::string &name) const {
    return !(is_old(name)) && variable_description(name).is_affine_dependent;
  }

  const std::string &
  model::org_variable(const std::string &name) const {
    GMM_ASSERT1(is_affine_dependent_variable(name),
                "For affine dependent variables only");
    return variable_description(name).org_name;
  }

  const scalar_type &
  model::factor_of_variable(const std::string &name) const {
    return variable_description(name).alpha;
  }

  void model::set_factor_of_variable(const std::string &name, scalar_type a) {
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
    if (it->second.alpha != a) {
      it->second.alpha = a;
      for (auto &v_num : it->second.v_num_data) v_num = act_counter();
    }
  }

  bool model::is_im_data(const std::string &name) const {
    return variable_description(no_old_prefix_name(name)).imd != 0;
  }

  const im_data *
  model::pim_data_of_variable(const std::string &name) const {
    return variable_description(no_old_prefix_name(name)).imd;
  }

  const gmm::uint64_type &
  model::version_number_of_data_variable(const std::string &name,
                                         size_type niter) const {
    VAR_SET::const_iterator it = find_variable(name);
    if (niter == size_type(-1)) niter = it->second.default_iter;
    return it->second.v_num_data[niter];
  }

  size_type model::nb_dof(bool with_internal) const {
    context_check();
    if (act_size_to_be_done) actualize_sizes();
    if (complex_version)
      return gmm::vect_size(crhs);
    else if (with_internal && gmm::vect_size(full_rrhs))
      return gmm::vect_size(full_rrhs);
    else
      return gmm::vect_size(rrhs);
  }

  void model::resize_global_system() const {

    size_type full_size = 0;
    for (auto &&v : variables)
      if (v.second.is_variable) {
        if (v.second.is_disabled)
          v.second.I  = gmm::sub_interval(0,0);
        else if (!v.second.is_affine_dependent && !v.second.is_internal) {
          v.second.I = gmm::sub_interval(full_size, v.second.size());
          full_size += v.second.size();
        }
      }
    size_type primary_size = full_size;

    for (auto &&v : variables)
      if (v.second.is_internal && v.second.is_enabled()) { // is_internal_variable()
        v.second.I = gmm::sub_interval(full_size, v.second.size());
        full_size += v.second.size();
      }

    for (auto &&v : variables)
      if (v.second.is_affine_dependent) {
        v.second.I = variables.find(v.second.org_name)->second.I;
        v.second.set_size();
      }

    if (complex_version) {
      gmm::resize(cTM, primary_size, primary_size);
      gmm::resize(crhs, primary_size);
    }
    else {
      gmm::resize(rTM, primary_size, primary_size);
      gmm::resize(rrhs, primary_size);
    }

    if (full_size > primary_size) {
      GMM_ASSERT1(has_internal_variables(), "Internal error");
      gmm::resize(internal_rTM, full_size-primary_size, primary_size);
      gmm::resize(full_rrhs, full_size);
      gmm::resize(internal_sol, full_size-primary_size);
    } else {
      GMM_ASSERT1(!(has_internal_variables()), "Internal error");
      gmm::resize(internal_rTM, 0, 0);
      full_rrhs.clear();
    }

    for (dal::bv_visitor ib(valid_bricks); !ib.finished(); ++ib)
      for (const term_description &term : bricks[ib].tlist)
        if (term.is_global) {
          bricks[ib].terms_to_be_computed = true;
          break;
        }
  }

  void model::actualize_sizes() const {
    // cout << "begin act size" << endl;
    bool actualized = false;
    getfem::local_guard lock = locks_.get_lock();
    if (actualized) return; // If multiple threads are calling the method

    act_size_to_be_done = false;

    std::map<std::string, std::vector<std::string> > multipliers;
    std::set<std::string> tobedone;

//     #if GETFEM_PARA_LEVEL > 1
//     int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk);
//     double t_ref = MPI_Wtime();
//     cout << "Actualize size called from thread " << rk << endl;
//     #endif


    // In case of change in fems or mims, linear terms have to be recomputed
    // We could select which brick is to be recomputed if we would be able
    // to know which fem or mim is changed
    // -> already taken into account in update_brick()
    // for (dal::bv_visitor ib(valid_bricks); !ib.finished(); ++ib)
    //  bricks[ib].terms_to_be_computed = true;

    for (auto &&v : variables) {
      const std::string &vname = v.first;
      var_description &vdescr = v.second;
      if (vdescr.mf && !vdescr.is_affine_dependent) {
        if ((vdescr.filter & VDESCRFILTER_CTERM)
            || (vdescr.filter & VDESCRFILTER_INFSUP)) {
          VAR_SET::iterator vfilt = variables.find(vdescr.filter_var);
          GMM_ASSERT1(vfilt != variables.end(), "The primal variable of the"
                      " multiplier does not exist : " << vdescr.filter_var);
          GMM_ASSERT1(vfilt->second.mf, "The primal variable of the "
                      "multiplier should be a fem variable");
          multipliers[vdescr.filter_var].push_back(vname);
          if (vdescr.v_num < vdescr.mf->version_number() ||
              vdescr.v_num < vfilt->second.mf->version_number()) {
            tobedone.insert(vdescr.filter_var);
          }
        }
        switch (vdescr.filter) {
        case VDESCRFILTER_NO:
          if (vdescr.v_num < vdescr.mf->version_number()) {
            vdescr.set_size();
            vdescr.v_num = act_counter();
          }
          break;
        case VDESCRFILTER_REGION:
          if (vdescr.v_num < vdescr.mf->version_number()) {
            dal::bit_vector
              dor = vdescr.mf->dof_on_region(vdescr.filter_region);
            vdescr.partial_mf->adapt(dor);
            vdescr.set_size();
            vdescr.v_num = act_counter();
          }
          break;
        default : break;
        }
      }

      if (vdescr.imd != 0
          && vdescr.v_num < vdescr.imd->version_number()) {
        vdescr.set_size();
        vdescr.v_num = act_counter();
      }
    }

    for (auto &&v : variables) {
      var_description &vdescr = v.second;
      if (vdescr.mf && !(vdescr.is_affine_dependent) &&
          ((vdescr.filter & VDESCRFILTER_CTERM)
           || (vdescr.filter & VDESCRFILTER_INFSUP))) {
        if (tobedone.count(vdescr.filter_var)) {
          // This step forces the recomputation of corresponding bricks.
          // A test to check if a modification is really necessary could
          // be done first ... (difficult to coordinate with other
          // multipliers)
          dal::bit_vector alldof; alldof.add(0, vdescr.mf->nb_dof());
          vdescr.partial_mf->adapt(alldof);
          vdescr.set_size();
          vdescr.v_num = act_counter();
        }
      }
    }

    resize_global_system();

    for (const std::string &vname : tobedone) {
//       #if GETFEM_PARA_LEVEL > 1
//       double tt_ref = MPI_Wtime();
//       if (!rk) cout << "compute size of multipliers for " << vname
//                     << endl;
//       #endif

      const std::vector<std::string> &mults = multipliers[vname];
      const var_description &vdescr = variable_description(vname);

      gmm::col_matrix< gmm::rsvector<scalar_type> > MGLOB;
      if (mults.size() > 1) {
        size_type s = 0;
        for (const std::string &mult : mults)
          s += variable_description(mult).mf->nb_dof();
        gmm::resize(MGLOB, vdescr.mf->nb_dof(), s);
      }
      size_type s = 0;
      std::set<size_type> glob_columns;
      for (const std::string &multname : mults) {
        var_description &multdescr = variables.find(multname)->second;

        // Obtaining the coupling matrix between the multipier and
        // the primal variable. A search is done on all the terms of the
        // model. Only the corresponding linear terms are added.
        // If no linear term is available, a mass matrix is used.
        gmm::col_matrix< gmm::rsvector<scalar_type> >
          MM(vdescr.associated_mf().nb_dof(), multdescr.mf->nb_dof());
        bool termadded = false;

        if (multdescr.filter & VDESCRFILTER_CTERM) {

          for (dal::bv_visitor ib(valid_bricks); !ib.finished(); ++ib) {
            const brick_description &brick = bricks[ib];
            bool bupd = false;
            bool cplx = is_complex() && brick.pbr->is_complex();

            if (brick.tlist.size() == 0) {
              bool varc = false, multc = false;
              for (const std::string &var : brick.vlist) {
                if (multname.compare(var) == 0) multc = true;
                if (vname.compare(var) == 0) varc = true;
              }
              if (multc && varc) {
                GMM_ASSERT1(!cplx, "Sorry, not taken into account");
                generic_expressions.clear();
                brick.terms_to_be_computed = true;
                update_brick(ib, BUILD_MATRIX);
                if (generic_expressions.size()) {
                  GMM_TRACE2("Generic assembly for actualize sizes");
                  {
                    gmm::clear(rTM);
                    accumulated_distro<decltype(rTM)>  distro_rTM(rTM);
                    GETFEM_OMP_PARALLEL(
                        ga_workspace workspace(*this);
                        for (const auto &ge : generic_expressions)
                          workspace.add_expression(ge.expr, ge.mim, ge.region,
                                                   2, ge.secondary_domain);
                        workspace.set_assembled_matrix(distro_rTM);
                        workspace.assembly(2);
                    );
                  } //accumulated_distro scope
                  gmm::add(gmm::sub_matrix(rTM, vdescr.I, multdescr.I), MM);
                  gmm::add(gmm::transposed
                           (gmm::sub_matrix(rTM, multdescr.I, vdescr.I)), MM);
                  bupd = false;
                }
              }
            }

            for (size_type j = 0; j < brick.tlist.size(); ++j) {
              const term_description &term = brick.tlist[j];

              if (term.is_matrix_term) {
                if (term.is_global) {
                  bool varc = false, multc = false;
                  for (const std::string &var : brick.vlist) {
                    if (multname.compare(var) == 0) multc = true;
                    if (vname.compare(var) == 0) varc = true;
                  }
                  if (multc && varc) {
                    GMM_ASSERT1(!cplx, "Sorry, not taken into account");
                    generic_expressions.clear();
                    if (!bupd) {
                      brick.terms_to_be_computed = true;
                      update_brick(ib, BUILD_MATRIX);
                      bupd = true;
                    }
                    gmm::add(gmm::sub_matrix(brick.rmatlist[j],
                                             multdescr.I, vdescr.I),
                             MM);
                    gmm::add(gmm::transposed(gmm::sub_matrix
                                             (brick.rmatlist[j],
                                              vdescr.I, multdescr.I)),
                             MM);
                    termadded = true;
                  }
                } else if (multname.compare(term.var1) == 0 &&
                           vname.compare(term.var2) == 0) {
                  if (!bupd) {
                    brick.terms_to_be_computed = true;
                    update_brick(ib, BUILD_MATRIX);
                    bupd = true;
                  }
                  if (cplx)
                    gmm::add(gmm::transposed(gmm::real_part(brick.cmatlist[j])),
                             MM);
                  else
                    gmm::add(gmm::transposed(brick.rmatlist[j]), MM);
                  termadded = true;

                } else if (multname.compare(term.var2) == 0 &&
                           vname.compare(term.var1) == 0) {
                  if (!bupd) {
                    brick.terms_to_be_computed = true;
                    update_brick(ib, BUILD_MATRIX);
                    bupd = true;
                  }
                  if (cplx)
                    gmm::add(gmm::real_part(brick.cmatlist[j]), MM);
                  else
                    gmm::add(brick.rmatlist[j], MM);
                  termadded = true;
                }
              }
            }
          }

          if (!termadded)
            GMM_WARNING1("No term found to filter multiplier " << multname
                         << ". Variable is cancelled");
        } else if (multdescr.filter & VDESCRFILTER_INFSUP) {
          mesh_region rg(multdescr.filter_region);
          multdescr.filter_mim->linked_mesh().intersect_with_mpi_region(rg);
          asm_mass_matrix(MM, *(multdescr.filter_mim), vdescr.associated_mf(),
                          *(multdescr.mf), rg);
        }

        MPI_SUM_SPARSE_MATRIX(MM);

        //
        // filtering
        //
        std::set<size_type> columns;
        gmm::range_basis(MM, columns);
        if (columns.size() == 0)
          GMM_WARNING1("Empty basis found for multiplier " << multname);

        if (mults.size() > 1) {
          gmm::copy(MM, gmm::sub_matrix
                    (MGLOB,
                     gmm::sub_interval(0, vdescr.associated_mf().nb_dof()),
                     gmm::sub_interval(s, multdescr.mf->nb_dof())));
          for (const size_type &icol : columns)
            glob_columns.insert(s + icol);
          s += multdescr.mf->nb_dof();
        } else {
          dal::bit_vector kept;
          for (const size_type &icol : columns)
            kept.add(icol);
          if (multdescr.filter & VDESCRFILTER_REGION)
            kept &= multdescr.mf->dof_on_region(multdescr.filter_region);
          multdescr.partial_mf->adapt(kept);
          multdescr.set_size();
          multdescr.v_num = act_counter();
        }
      }

//         #if GETFEM_PARA_LEVEL > 1
//         if (!rk) cout << "Range basis for  multipliers for " << vname << " time : " << MPI_Wtime()-tt_ref << endl;
//         #endif

      if (mults.size() > 1) {
        gmm::range_basis(MGLOB, glob_columns, 1E-12, gmm::col_major(), true);

//         #if GETFEM_PARA_LEVEL > 1
//         if (!rk) cout << "Producing partial mf for  multipliers for " << vname << " time : " << MPI_Wtime()-tt_ref << endl;
//         #endif

        s = 0;
        for (const std::string &multname : mults) {
          var_description &multdescr = variables.find(multname)->second;
          dal::bit_vector kept;
          size_type nbdof = multdescr.mf->nb_dof();
          for (const size_type &icol : glob_columns)
            if (icol >= s && icol < s + nbdof) kept.add(icol-s);
          if (multdescr.filter & VDESCRFILTER_REGION)
            kept &= multdescr.mf->dof_on_region(multdescr.filter_region);
          multdescr.partial_mf->adapt(kept);
          multdescr.set_size();
          multdescr.v_num = act_counter();
          s += multdescr.mf->nb_dof();
        }
      }
//       #if GETFEM_PARA_LEVEL > 1
//       if (!rk) cout << "End compute size of  multipliers for " << vname << " time : " << MPI_Wtime()-tt_ref << endl;
//       #endif
    }

    resize_global_system();
    actualized = true;
//     #if GETFEM_PARA_LEVEL > 1
//     cout << "Actualize sizes time from thread " << rk << " : " << MPI_Wtime()-t_ref << endl;

//     #endif

    // cout << "end act size" << endl;
  }


  void model::listvar(std::ostream &ost) const {
    if (variables.size() == 0)
      ost << "Model with no variable nor data" << endl;
    else {
      ost << "List of model variables and data:" << endl;
      for (int vartype=0; vartype < 3; ++vartype)
        for (const auto &v : variables) {
          const var_description &vdescr = v.second;
          bool is_variable = vdescr.is_variable;
          bool is_disabled = is_variable && is_disabled_variable(v.first);
          if (vartype == 0) {      // Only enabled variables
            if (!is_variable || is_disabled) continue;
          } else if (vartype == 1) { // Only disabled variables
            if (!is_disabled) continue;
          } else if (vartype == 2) { // Only data
            if (is_variable) continue;
          }
          ost << (is_variable ? "Variable       " : "Data           ");
          ost << std::setw(30) << std::left << v.first;
          ost << std::setw(2) << std::right << vdescr.n_iter;
          ost << ((vdescr.n_iter == 1) ? " copy   " : " copies ");
          ost << (vdescr.mf ? "fem dependant " : "constant size ");
          ost << std::setw(8) << std::right << vdescr.size();
          if (is_complex()) ost << " complex";
          ost << ((vdescr.size() > 1) ? " doubles." : " double.");
          ost << (is_disabled ? "\t (disabled)" : "\t           ");
          if (vdescr.imd != 0) ost << "\t (is im_data)";
          if (vdescr.is_affine_dependent) ost << "\t (is affine dependent)";
          ost << endl;
        }
      for (const auto &vargroup : variable_groups) {
        ost << "Variable group " << std::setw(30) << std::left
            << vargroup.first;
        if (vargroup.second.size()) {
          bool first(true);
          for (const std::string &vname : vargroup.second) {
            ost << (first ? " " : ", ") << vname;
            first = false;
          }
          ost << endl;
        } else
          ost << " empty" << endl;
      }
    }
  }

  void model::listresiduals(std::ostream &ost) const {
    context_check(); if (act_size_to_be_done) actualize_sizes();
    if (variables.size() == 0)
      ost << "Model with no variable nor data" << endl;
    else {
      bool firstvar(true);
      for (const auto &v : variables) {
        if (v.second.is_variable) {
          const model_real_plain_vector &rhs = v.second.is_internal
                                             ? full_rrhs : rrhs;
          const gmm::sub_interval &II = interval_of_variable(v.first);
          scalar_type res = gmm::vect_norm2(gmm::sub_vector(rhs, II));
          if (!firstvar) cout << ", ";
          ost << "res_" << v.first << "= " << std::setw(11) << res;
          firstvar = false;
        }
      }
      ost << endl;
    }
  }

  void model::add_fixed_size_variable(const std::string &name, size_type size,
                                      size_type niter) {
    bgeot::multi_index sizes(1);
    sizes[0] = size;
    add_fixed_size_variable(name, sizes, niter);
  }

  void model::add_fixed_size_variable(const std::string &name,
                                      const bgeot::multi_index &sizes,
                                      size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(true, is_complex(), 0, 0, niter));
    variables[name].qdims = sizes;
    act_size_to_be_done = true;
    variables[name].set_size();
    GMM_ASSERT1(variables[name].qdim(),
                "Variables of null size are not allowed");
  }

  void model::resize_fixed_size_variable(const std::string &name,
                                         size_type size) {
    bgeot::multi_index sizes(1);
    sizes[0] = size;
    resize_fixed_size_variable(name, sizes);
  }

  void model::resize_fixed_size_variable(const std::string &name,
                                         const bgeot::multi_index &sizes) {
    GMM_ASSERT1(variables[name].mf == 0,
                "Cannot explicitly resize a fem variable or data");
    GMM_ASSERT1(variables[name].imd == 0,
                "Cannot explicitly resize an im variable or data");
    variables[name].qdims = sizes;
    variables[name].set_size();
  }

  void model::add_fixed_size_data(const std::string &name, size_type size,
                                  size_type niter) {
    bgeot::multi_index sizes(1);
    sizes[0] = size;
    add_fixed_size_data(name, sizes, niter);
  }

  void model::add_fixed_size_data(const std::string &name,
                                  const bgeot::multi_index &sizes,
                                  size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(false, is_complex(), 0, 0, niter));
    variables[name].qdims = sizes;
    GMM_ASSERT1(variables[name].qdim(), "Data of null size are not allowed");
    variables[name].set_size();
  }

  void model::add_initialized_matrix_data(const std::string &name,
                                          const base_matrix &M) {
    add_fixed_size_data(name, bgeot::multi_index(gmm::mat_nrows(M),
                                                 gmm::mat_ncols(M)));
    GMM_ASSERT1(!(is_complex()), "Sorry, complex version to be done");
    gmm::copy(M.as_vector(), set_real_variable(name));
  }

  void model::add_initialized_matrix_data(const std::string &name,
                                          const base_complex_matrix &M) {
    add_fixed_size_data(name, bgeot::multi_index(gmm::mat_nrows(M),
                                                 gmm::mat_ncols(M)));
    GMM_ASSERT1(!(is_complex()), "Sorry, complex version to be done");
    gmm::copy(M.as_vector(), set_complex_variable(name));
  }

  void model::add_initialized_tensor_data(const std::string &name,
                                         const base_tensor &t) {
    add_fixed_size_data(name, t.sizes(), 1);
    GMM_ASSERT1(!(is_complex()), "Sorry, complex version to be done");
    gmm::copy(t.as_vector(), set_real_variable(name));
  }

  void model::add_initialized_tensor_data(const std::string &name,
                                          const base_complex_tensor &t) {
    add_fixed_size_data(name, t.sizes(), 1);
    GMM_ASSERT1(!(is_complex()), "Sorry, complex version to be done");
    gmm::copy(t.as_vector(), set_complex_variable(name));
  }

  void model::add_im_variable(const std::string &name, const im_data &imd,
                              size_type niter) {
    check_name_validity(name);
    variables.emplace(name,
                      var_description(true, is_complex(), 0, &imd, niter));
    variables[name].set_size();
    add_dependency(imd);
    act_size_to_be_done = true;
  }

  void model::add_internal_im_variable(const std::string &name,
                                       const im_data &imd) {
    add_im_variable(name, imd);
    variables[name].is_internal = true;
  }

  void model::add_im_data(const std::string &name, const im_data &imd,
                          size_type niter) {
    check_name_validity(name);
    variables.emplace(name,
                      var_description(false, is_complex(), 0, &imd, niter));
    variables[name].set_size();
    add_dependency(imd);
  }

  void model::add_fem_variable(const std::string &name, const mesh_fem &mf,
                               size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(true, is_complex(), &mf, 0, niter,
                                            VDESCRFILTER_NO));
    variables[name].set_size();
    add_dependency(mf);
    act_size_to_be_done = true;
    leading_dim = std::max(leading_dim, mf.linked_mesh().dim());
  }

  void model::add_filtered_fem_variable(const std::string &name,
                                        const mesh_fem &mf,
                                        size_type region, size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(true, is_complex(), &mf, 0, niter,
                                            VDESCRFILTER_REGION, region));
    variables[name].set_size();
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_affine_dependent_variable(const std::string &name,
                                            const std::string &org_name,
                                            scalar_type alpha) {
    check_name_validity(name);
    VAR_SET::const_iterator it = find_variable(org_name);
    GMM_ASSERT1(it->second.is_variable && !(it->second.is_affine_dependent),
                "The original variable should be a variable");
    variables.emplace(name, variables[org_name]);
    variables[name].is_affine_dependent = true;
    variables[name].org_name = org_name;
    variables[name].alpha = alpha;
    variables[name].set_size();
  }

  void model::add_fem_data(const std::string &name, const mesh_fem &mf,
                           dim_type qdim, size_type niter) {
    bgeot::multi_index sizes(1);
    sizes[0] = qdim;
    add_fem_data(name, mf, sizes, niter);
  }

  void model::add_fem_data(const std::string &name, const mesh_fem &mf,
                           const bgeot::multi_index &sizes, size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(false, is_complex(), &mf, 0, niter,
                                            VDESCRFILTER_NO));
    variables[name].qdims = sizes;
    GMM_ASSERT1(variables[name].qdim(), "Data of null size are not allowed");
    variables[name].set_size();
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             const std::string &primal_name,
                             size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(true, is_complex(), &mf, 0, niter,
                                            VDESCRFILTER_CTERM, size_type(-1),
                                            primal_name));
    variables[name].set_size();
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             size_type region, const std::string &primal_name,
                             size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(true, is_complex(), &mf, 0, niter,
                                            VDESCRFILTER_REGION_CTERM, region,
                                            primal_name));
    variables[name].set_size();
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             const std::string &primal_name,
                             const mesh_im &mim,
                             size_type region, size_type niter) {
    check_name_validity(name);
    variables.emplace(name, var_description(true, is_complex(), &mf, 0, niter,
                                            VDESCRFILTER_INFSUP, region,
                                            primal_name, &mim));
    variables[name].set_size();
    act_size_to_be_done = true;
    add_dependency(mf);
    add_dependency(mim);
  }

  void model::disable_variable(const std::string &name) {
    enable_variable(name, false);
  }

  void model::enable_variable(const std::string &name, bool enabled) {
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
    it->second.is_disabled = !enabled;
    for (auto &&v : variables) {
      if (((v.second.filter & VDESCRFILTER_INFSUP) ||
           (v.second.filter & VDESCRFILTER_CTERM))
          && name.compare(v.second.filter_var) == 0) {
        v.second.is_disabled = !enabled;
      }
      if (v.second.is_variable && v.second.is_affine_dependent
          && name.compare(v.second.org_name) == 0)
        v.second.is_disabled = !enabled;
    }
    if (!act_size_to_be_done) resize_global_system();
  }

  bool model::variable_exists(const std::string &name) const {
    return variables.count(no_old_prefix_name(name)) > 0;
  }

  void model::add_macro(const std::string &name, const std::string &expr) {
    check_name_validity(name.substr(0, name.find("(")));
    macro_dict.add_macro(name, expr);
  }

  void model::del_macro(const std::string &name)
  { macro_dict.del_macro(name); }

  void model::delete_brick(size_type ib) {
     GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
     valid_bricks.del(ib);
     active_bricks.del(ib);

     for (size_type i = 0; i < bricks[ib].mims.size(); ++i) {
       const mesh_im *mim = bricks[ib].mims[i];
       bool found = false;
       for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
         for  (size_type j = 0; j < bricks[ibb].mims.size(); ++j)
           if (bricks[ibb].mims[j] == mim) found = true;
       }
       for (const auto &v : variables) {
         if (v.second.mf && (v.second.filter & VDESCRFILTER_INFSUP) &&
             mim == v.second.filter_mim) found = true;
        }
       if (!found) sup_dependency(*mim);
     }

     is_linear_ = is_symmetric_ = is_coercive_ = true;
     for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
       is_linear_ = is_linear_ && bricks[ibb].pbr->is_linear();
       is_symmetric_ = is_symmetric_ &&  bricks[ibb].pbr->is_symmetric();
       is_coercive_ = is_coercive_ &&  bricks[ibb].pbr->is_coercive();
     }
     bricks[ib] = brick_description();
  }

  void model::delete_variable(const std::string &varname) {
    for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
      for (const auto &vname : bricks[ibb].vlist)
        GMM_ASSERT1(varname.compare(vname),
                    "Cannot delete a variable which is still used by a brick");
      for (const auto &dname : bricks[ibb].dlist)
        GMM_ASSERT1(varname.compare(dname),
                    "Cannot delete a data which is still used by a brick");
    }

    VAR_SET::const_iterator it = find_variable(varname);

    if (it->second.mf) {
      const mesh_fem *mf = it->second.mf;
      bool found = false;
      for(VAR_SET::iterator it2 = variables.begin();
          it2 != variables.end(); ++it2) {
        if (it != it2 && it2->second.mf && mf == it2->second.mf)
          found = true;
      }
      if (!found) sup_dependency(*mf);

      if (it->second.filter & VDESCRFILTER_INFSUP) {
        const mesh_im *mim = it->second.filter_mim;
        found = false;
        for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
          for  (size_type j = 0; j < bricks[ibb].mims.size(); ++j)
            if (bricks[ibb].mims[j] == mim) found = true;
        }
        for (VAR_SET::iterator it2 = variables.begin();
             it2 != variables.end(); ++it2) {
          if (it != it2 && it2->second.mf &&
              (it2->second.filter & VDESCRFILTER_INFSUP) &&
              mim == it2->second.filter_mim) found = true;
        }
        if (!found) sup_dependency(*mim);
      }
    }

    if (it->second.imd != 0) sup_dependency(*(it->second.imd));

    variables.erase(varname);
    act_size_to_be_done = true;
  }

  size_type model::add_brick(pbrick pbr, const varnamelist &varnames,
                             const varnamelist &datanames,
                             const termlist &terms,
                             const mimlist &mims, size_type region) {
    size_type ib = valid_bricks.first_false();

    for (size_type i = 0; i < terms.size(); ++i)
      if (terms[i].is_global && terms[i].is_matrix_term && pbr->is_linear())
        GMM_ASSERT1(false, "Global linear matrix terms are not allowed");

    if (ib == bricks.size())
      bricks.push_back(brick_description(pbr, varnames, datanames, terms,
                                         mims, region));
    else
      bricks[ib] = brick_description(pbr, varnames, datanames, terms,
                                     mims, region);
    active_bricks.add(ib);
    valid_bricks.add(ib);

    // The brick itself already reacts to a mesh_im change in update_brick()
    // for (size_type i = 0; i < bricks[ib].mims.size(); ++i)
    //   add_dependency(*(bricks[ib].mims[i]));

    GMM_ASSERT1(pbr->is_real() || is_complex(),
                "Impossible to add a complex brick to a real model");
    if (is_complex() && pbr->is_complex()) {
      bricks[ib].cmatlist.resize(terms.size());
      bricks[ib].cveclist[0].resize(terms.size());
      bricks[ib].cveclist_sym[0].resize(terms.size());
    } else {
      bricks[ib].rmatlist.resize(terms.size());
      bricks[ib].rveclist[0].resize(terms.size());
      bricks[ib].rveclist_sym[0].resize(terms.size());
    }
    is_linear_ = is_linear_ && pbr->is_linear();
    is_symmetric_ = is_symmetric_ && pbr->is_symmetric();
    is_coercive_ = is_coercive_ && pbr->is_coercive();

    for (const auto &vname : varnames)
      GMM_ASSERT1(variables.count(vname),
                  "Undefined model variable " << vname);
    // cout << "dl == " << datanames << endl;
    for (const auto &dname : datanames)
      GMM_ASSERT1(variables.count(dname),
                  "Undefined model data or variable " << dname);

    return ib;
  }

  void model::add_mim_to_brick(size_type ib, const mesh_im &mim) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].mims.push_back(&mim);
    add_dependency(mim);
  }

  void model::change_terms_of_brick(size_type ib, const termlist &terms) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].tlist = terms;
    if (is_complex() && bricks[ib].pbr->is_complex()) {
      bricks.back().cmatlist.resize(terms.size());
      bricks.back().cveclist[0].resize(terms.size());
      bricks.back().cveclist_sym[0].resize(terms.size());
    } else {
      bricks.back().rmatlist.resize(terms.size());
      bricks.back().rveclist[0].resize(terms.size());
      bricks.back().rveclist_sym[0].resize(terms.size());
    }
  }

  void model::change_variables_of_brick(size_type ib, const varnamelist &vl) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].vlist = vl;
    for (const auto &v : vl)
      GMM_ASSERT1(variables.count(v), "Undefined model variable " << v);
  }

  void model::change_data_of_brick(size_type ib, const varnamelist &dl) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].dlist = dl;
    for (const auto &v : dl)
      GMM_ASSERT1(variables.count(v), "Undefined model variable " << v);
  }

  void model::change_mims_of_brick(size_type ib, const mimlist &ml) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].mims = ml;
    for (const auto &mim : ml) add_dependency(*mim);
  }

  void model::change_update_flag_of_brick(size_type ib, bool flag) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].is_update_brick = flag;
  }

  void model::set_time(scalar_type t, bool to_init) {
    static const std::string varname("t");
    VAR_SET::iterator it = variables.find(varname);
    if (it == variables.end()) {
      add_fixed_size_data(varname, 1);
    } else {
      GMM_ASSERT1(it->second.size() == 1, "Time data should be of size 1");
    }
    if (it == variables.end() || to_init) {
      if (is_complex())
        set_complex_variable(varname)[0] = complex_type(t);
      else
        set_real_variable(varname)[0] = t;
    }
  }

  scalar_type model::get_time() {
    static const std::string varname("t");
    set_time(scalar_type(0), false);
    if (is_complex())
      return gmm::real(complex_variable(varname)[0]);
    else
      return real_variable(varname)[0];
  }

  void model::call_init_affine_dependent_variables(int version) {
    for (VAR_SET::iterator it = variables.begin();
         it != variables.end(); ++it) {
      var_description &vdescr = it->second;
      if (vdescr.is_variable && vdescr.ptsc) {
        if (version == 2)
          vdescr.ptsc->init_affine_dependent_variables_precomputation(*this);
        else
          vdescr.ptsc->init_affine_dependent_variables(*this);
      }
    }
  }

  void model::shift_variables_for_time_integration() {
    for (VAR_SET::iterator it = variables.begin();
         it != variables.end(); ++it)
      if (it->second.is_variable && it->second.ptsc)
        it->second.ptsc->shift_variables(*this);
  }

  void model::add_time_integration_scheme(const std::string &varname,
                                          ptime_scheme ptsc) {
    VAR_SET::iterator it = variables.find(varname);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << varname);
    GMM_ASSERT1(it->second.is_variable && !(it->second.is_affine_dependent),
                "Cannot apply an integration scheme to " << varname);
    it->second.ptsc = ptsc;
//     if (!first_step)
//       GMM_WARNING2("When you apply a scheme to new variable or change the "
//                    "scheme of a variable after the first time step, "
//                    "the precomputation of time derivative will not be "
//                    "executed. Caution: You have to care by yourself of "
//                    "the compatbility of the operation");
    time_integration = 1;
  }

  void model::copy_init_time_derivative() {

    for (VAR_SET::iterator it = variables.begin();
         it != variables.end(); ++it)
      if (it->second.is_variable && it->second.ptsc) {

        std::string name_v, name_previous_v;
        it->second.ptsc->time_derivative_to_be_initialized(name_v,
                                                          name_previous_v);

        if (name_v.size()) {
          if (is_complex()) {
            model_complex_plain_vector v0 = complex_variable(name_v);
            gmm::copy(v0, set_complex_variable(name_previous_v));
          } else {
            const model_real_plain_vector &v0 = real_variable(name_v);
            gmm::copy(v0, set_real_variable(name_previous_v));
          }
        }
      }
  }

  // ----------------------------------------------------------------------
  //
  // Theta-method scheme for first order problems
  //
  // ----------------------------------------------------------------------

    class APIDECL first_order_theta_method_scheme
      : public virtual_time_scheme {

      std::string U, U0, V, V0;
      scalar_type theta;

    public:
      // V = (U-U0)/(theta*dt) - ((1-theta)/theta)*V0
      virtual void init_affine_dependent_variables(model &md) const {
        scalar_type dt = md.get_time_step();
        scalar_type a = scalar_type(1)/(theta*dt);
        scalar_type b = (scalar_type(1)-theta)/theta;
        md.set_factor_of_variable(V, a);
        if (md.is_complex()) {
           gmm::add(gmm::scaled(md.complex_variable(U0), -complex_type(a)),
                    gmm::scaled(md.complex_variable(V0), -complex_type(b)),
                    md.set_complex_constant_part(V));

        } else {
          gmm::add(gmm::scaled(md.real_variable(U0), -a),
                   gmm::scaled(md.real_variable(V0), -b),
                   md.set_real_constant_part(V));
        }
      }

      // V = (U-U0)/dt (backward Euler for precomputation)
      virtual void init_affine_dependent_variables_precomputation(model &md)
        const {
        scalar_type dt = md.get_time_step();
        md.set_factor_of_variable(V, scalar_type(1)/dt);
        if (md.is_complex()) {
          gmm::copy(gmm::scaled(md.complex_variable(U0), -complex_type(1)/dt),
                    md.set_complex_constant_part(V));

        } else {
          gmm::copy(gmm::scaled(md.real_variable(U0), -scalar_type(1)/dt),
                    md.set_real_constant_part(V));
        }
      }

      virtual void time_derivative_to_be_initialized
      (std::string &name_v, std::string &name_previous_v) const
      { if (theta != scalar_type(1)) { name_v = V; name_previous_v = V0; } }

      virtual void shift_variables(model &md) const {
        if (md.is_complex()) {
          gmm::copy(md.complex_variable(U), md.set_complex_variable(U0));
          gmm::copy(md.complex_variable(V), md.set_complex_variable(V0));
        } else {
          gmm::copy(md.real_variable(U), md.set_real_variable(U0));
          gmm::copy(md.real_variable(V), md.set_real_variable(V0));
        }
      }


      first_order_theta_method_scheme(model &md, std::string varname,
                                      scalar_type th) {
        U = varname;
        U0 = "Previous_" + U;
        V = "Dot_" + U;
        V0 = "Previous_Dot_" + U;
        theta = th;
        GMM_ASSERT1(theta > scalar_type(0) && theta <=  scalar_type(1),
                    "Invalid value of theta parameter for the theta-method");

        if (!(md.variable_exists(V)))
          md.add_affine_dependent_variable(V, U);
        const mesh_fem *mf =  md.pmesh_fem_of_variable(U);
        size_type s = md.is_complex() ? gmm::vect_size(md.complex_variable(U))
          : gmm::vect_size(md.real_variable(U));

        if (mf) {
          if (!(md.variable_exists(U0))) md.add_fem_data(U0, *mf);
          if (!(md.variable_exists(V0))) md.add_fem_data(V0, *mf);
        } else {
          if (!(md.variable_exists(U0))) md.add_fixed_size_data(U0, s);
          if (!(md.variable_exists(V0))) md.add_fixed_size_data(V0, s);
        }
      }


    };

  void add_theta_method_for_first_order(model &md, const std::string &varname,
                                        scalar_type theta) {
    ptime_scheme ptsc
      = std::make_shared<first_order_theta_method_scheme>(md, varname,theta);
    md.add_time_integration_scheme(varname, ptsc);
  }

  // ----------------------------------------------------------------------
  //
  // Theta-method for second order problems
  //
  // ----------------------------------------------------------------------

  class APIDECL second_order_theta_method_scheme
    : public virtual_time_scheme {

    std::string U, U0, V, V0, A, A0;
    scalar_type theta;

  public:
    // V = (U-U0)/(theta*dt) - dt*(1-theta)*V0
    // A = (U-U0)/(theta^2*dt^2) - V0/(theta^2*dt) - dt*(1-theta)*A0
    virtual void init_affine_dependent_variables(model &md) const {
      scalar_type dt = md.get_time_step();
      md.set_factor_of_variable(V, scalar_type(1)/(theta*dt));
      md.set_factor_of_variable(A, scalar_type(1)/(theta*theta*dt*dt));
      if (md.is_complex()) {
        gmm::add(gmm::scaled(md.complex_variable(U0),
                             -complex_type(1)/(theta*dt)),
                 gmm::scaled(md.complex_variable(V0),
                             -(complex_type(1)-complex_type(theta))/theta),
                 md.set_complex_constant_part(V));
        gmm::add(gmm::scaled(md.complex_variable(U0),
                             -complex_type(1)/(theta*theta*dt*dt)),
                 gmm::scaled(md.complex_variable(A0),
                             -(complex_type(1)-complex_type(theta))/theta),
                 md.set_complex_constant_part(A));
        gmm::add(gmm::scaled(md.complex_variable(V0),
                             -complex_type(1)/(theta*theta*dt)),
                 md.set_complex_constant_part(A));


      } else {
        gmm::add(gmm::scaled(md.real_variable(U0),
                             -scalar_type(1)/(theta*dt)),
                 gmm::scaled(md.real_variable(V0),
                             -(scalar_type(1)-theta)/theta),
                 md.set_real_constant_part(V));
        gmm::add(gmm::scaled(md.real_variable(U0),
                             -scalar_type(1)/(theta*theta*dt*dt)),
                 gmm::scaled(md.real_variable(A0),
                             -(scalar_type(1)-theta)/theta),
                 md.set_real_constant_part(A));
        gmm::add(gmm::scaled(md.real_variable(V0),
                             -scalar_type(1)/(theta*theta*dt)),
                 md.set_real_constant_part(A));

      }
    }

    // V = (U-U0)/dt (backward Euler for precomputation)
    // A = (U-U0)/(dt^2) - V0/dt
    virtual void init_affine_dependent_variables_precomputation(model &md)
      const {
      scalar_type dt = md.get_time_step();
      md.set_factor_of_variable(V, scalar_type(1)/dt);
      md.set_factor_of_variable(A, scalar_type(1)/(dt*dt));
      if (md.is_complex()) {
        gmm::copy(gmm::scaled(md.complex_variable(U0),
                              -complex_type(1)/dt),
                  md.set_complex_constant_part(V));
        gmm::add(gmm::scaled(md.complex_variable(U0),
                             -complex_type(1)/(dt*dt)),
                 gmm::scaled(md.complex_variable(V0),
                             -complex_type(1)/dt),
                 md.set_complex_constant_part(A));
      } else {
        gmm::copy(gmm::scaled(md.real_variable(U0),
                              -scalar_type(1)/dt),
                  md.set_real_constant_part(V));
        gmm::add(gmm::scaled(md.real_variable(U0),
                             -scalar_type(1)/(dt*dt)),
                 gmm::scaled(md.real_variable(V0),
                             -scalar_type(1)/dt),
                 md.set_real_constant_part(A));
      }
    }

    virtual void time_derivative_to_be_initialized
    (std::string &name_v, std::string &name_previous_v) const
    { if (theta != scalar_type(1)) { name_v = A; name_previous_v = A0; } }

    virtual void shift_variables(model &md) const {
      if (md.is_complex()) {
        gmm::copy(md.complex_variable(U), md.set_complex_variable(U0));
        gmm::copy(md.complex_variable(V), md.set_complex_variable(V0));
        gmm::copy(md.complex_variable(A), md.set_complex_variable(A0));
      } else {
        gmm::copy(md.real_variable(U), md.set_real_variable(U0));
        gmm::copy(md.real_variable(V), md.set_real_variable(V0));
        gmm::copy(md.real_variable(A), md.set_real_variable(A0));
      }
    }


    second_order_theta_method_scheme(model &md, std::string varname,
                                     scalar_type th) {
      U = varname;
      U0 = "Previous_" + U;
      V = "Dot_" + U;
      V0 = "Previous_Dot_" + U;
      A = "Dot2_" + U;
      A0 = "Previous_Dot2_" + U;
      theta = th;
      GMM_ASSERT1(theta > scalar_type(0) && theta <=  scalar_type(1),
                  "Invalid value of theta parameter for the theta-method");

      if (!(md.variable_exists(V)))
        md.add_affine_dependent_variable(V, U);
      if (!(md.variable_exists(A)))
        md.add_affine_dependent_variable(A, U);

      const mesh_fem *mf =  md.pmesh_fem_of_variable(U);
      size_type s = md.is_complex() ? gmm::vect_size(md.complex_variable(U))
        : gmm::vect_size(md.real_variable(U));

      if (mf) {
        if (!(md.variable_exists(U0))) md.add_fem_data(U0, *mf);
        if (!(md.variable_exists(V0))) md.add_fem_data(V0, *mf);
        if (!(md.variable_exists(A0))) md.add_fem_data(A0, *mf);
      } else {
        if (!(md.variable_exists(U0))) md.add_fixed_size_data(U0, s);
        if (!(md.variable_exists(V0))) md.add_fixed_size_data(V0, s);
        if (!(md.variable_exists(A0))) md.add_fixed_size_data(A0, s);
      }
    }


  };

  void add_theta_method_for_second_order(model &md, const std::string &varname,
                                         scalar_type theta) {
    ptime_scheme ptsc = std::make_shared<second_order_theta_method_scheme>
      (md,varname,theta);
    md.add_time_integration_scheme(varname, ptsc);
  }


  // ----------------------------------------------------------------------
  //
  // Newmark method for second order problems
  //
  // ----------------------------------------------------------------------

    class APIDECL Newmark_scheme
      : public virtual_time_scheme {

      std::string U, U0, V, V0, A, A0;
      scalar_type beta, gamma;

    public:
      // V = (U-U0)/(theta*dt) - dt*(1-theta)*V0
      // A = (U-U0)/(theta^2*dt^2) - V0/(theta^2*dt) - dt*(1-theta)*A0
      virtual void init_affine_dependent_variables(model &md) const {
        scalar_type dt = md.get_time_step();
        scalar_type a0 = scalar_type(1)/(beta*dt*dt), a1 =  dt*a0;
        scalar_type a2 = (scalar_type(1) - scalar_type(2)*beta)
          / (scalar_type(2)*beta);
        scalar_type b0 = gamma/(beta*dt), b1 = (beta-gamma)/beta;
        scalar_type b2 = dt*(scalar_type(1)-gamma/(scalar_type(2)*beta));

        md.set_factor_of_variable(V, b0);
        md.set_factor_of_variable(A, a0);
        if (md.is_complex()) {
          gmm::add(gmm::scaled(md.complex_variable(U0), -complex_type(b0)),
                   gmm::scaled(md.complex_variable(V0), complex_type(b1)),
                   md.set_complex_constant_part(V));
          gmm::add(gmm::scaled(md.complex_variable(A0), complex_type(b2)),
                   md.set_complex_constant_part(V));
          gmm::add(gmm::scaled(md.complex_variable(U0), -complex_type(a0)),
                   gmm::scaled(md.complex_variable(V0), -complex_type(a1)),
                   md.set_complex_constant_part(A));
          gmm::add(gmm::scaled(md.complex_variable(A0), -complex_type(a2)),
                   md.set_complex_constant_part(A));
        } else {
          gmm::add(gmm::scaled(md.real_variable(U0), -b0),
                   gmm::scaled(md.real_variable(V0), b1),
                   md.set_real_constant_part(V));
          gmm::add(gmm::scaled(md.real_variable(A0), b2),
                   md.set_real_constant_part(V));
          gmm::add(gmm::scaled(md.real_variable(U0), -a0),
                   gmm::scaled(md.real_variable(V0), -a1),
                   md.set_real_constant_part(A));
          gmm::add(gmm::scaled(md.real_variable(A0), -a2),
                   md.set_real_constant_part(A));

        }
      }

      // V = (U-U0)/dt (backward Euler for precomputation)
      // A = (U-U0)/(dt^2) - V0/dt
      virtual void init_affine_dependent_variables_precomputation(model &md)
        const {
        scalar_type dt = md.get_time_step();
        md.set_factor_of_variable(V, scalar_type(1)/dt);
        md.set_factor_of_variable(A, scalar_type(1)/(dt*dt));
        if (md.is_complex()) {
          gmm::copy(gmm::scaled(md.complex_variable(U0),
                                -complex_type(1)/dt),
                    md.set_complex_constant_part(V));
          gmm::add(gmm::scaled(md.complex_variable(U0),
                               -complex_type(1)/(dt*dt)),
                   gmm::scaled(md.complex_variable(V0),
                               -complex_type(1)/dt),
                   md.set_complex_constant_part(A));
        } else {
          gmm::copy(gmm::scaled(md.real_variable(U0),
                                -scalar_type(1)/dt),
                    md.set_real_constant_part(V));
          gmm::add(gmm::scaled(md.real_variable(U0),
                               -scalar_type(1)/(dt*dt)),
                   gmm::scaled(md.real_variable(V0),
                               -scalar_type(1)/dt),
                   md.set_real_constant_part(A));
        }
      }

      virtual void time_derivative_to_be_initialized
      (std::string &name_v, std::string &name_previous_v) const {
        if (beta != scalar_type(0.5) || gamma != scalar_type(1))
          { name_v = A; name_previous_v = A0; }
      }

      virtual void shift_variables(model &md) const {
        if (md.is_complex()) {
          gmm::copy(md.complex_variable(U), md.set_complex_variable(U0));
          gmm::copy(md.complex_variable(V), md.set_complex_variable(V0));
          gmm::copy(md.complex_variable(A), md.set_complex_variable(A0));
        } else {
          gmm::copy(md.real_variable(U), md.set_real_variable(U0));
          gmm::copy(md.real_variable(V), md.set_real_variable(V0));
          gmm::copy(md.real_variable(A), md.set_real_variable(A0));
        }
      }


      Newmark_scheme(model &md, std::string varname,
                     scalar_type be, scalar_type ga) {
        U = varname;
        U0 = "Previous_" + U;
        V = "Dot_" + U;
        V0 = "Previous_Dot_" + U;
        A = "Dot2_" + U;
        A0 = "Previous_Dot2_" + U;
        beta = be; gamma = ga;
        GMM_ASSERT1(beta > scalar_type(0) && beta <=  scalar_type(1)
                    && gamma >= scalar_type(0.5) && gamma <=  scalar_type(1),
                    "Invalid parameter values for the Newmark scheme");

        if (!(md.variable_exists(V)))
          md.add_affine_dependent_variable(V, U);
        if (!(md.variable_exists(A)))
          md.add_affine_dependent_variable(A, U);

        const mesh_fem *mf =  md.pmesh_fem_of_variable(U);
        size_type s = md.is_complex() ? gmm::vect_size(md.complex_variable(U))
          : gmm::vect_size(md.real_variable(U));

        if (mf) {
          if (!(md.variable_exists(U0))) md.add_fem_data(U0, *mf);
          if (!(md.variable_exists(V0))) md.add_fem_data(V0, *mf);
          if (!(md.variable_exists(A0))) md.add_fem_data(A0, *mf);
        } else {
          if (!(md.variable_exists(U0))) md.add_fixed_size_data(U0, s);
          if (!(md.variable_exists(V0))) md.add_fixed_size_data(V0, s);
          if (!(md.variable_exists(A0))) md.add_fixed_size_data(A0, s);
        }
      }


    };

  void add_Newmark_scheme(model &md, const std::string &varname,
                          scalar_type beta, scalar_type gamma) {
    ptime_scheme ptsc = std::make_shared<Newmark_scheme>
      (md, varname, beta, gamma);
    md.add_time_integration_scheme(varname, ptsc);
  }

  // ----------------------------------------------------------------------
  //
  // Houbolt method
  //
  // ----------------------------------------------------------------------

    class APIDECL Houbolt_scheme
      : public virtual_time_scheme {

      std::string U, U01, U02, U03, V, A;

    public:
      // V = 1/(6*dt)*(11*U-18*U01+9*U02-2*U03)
      // A = 1/(dt**2)*(2*U-5*U01+4*U02-U03)
      virtual void init_affine_dependent_variables(model &md) const {
        scalar_type dt = md.get_time_step();
        scalar_type a0 = scalar_type(2)/(dt*dt);
        scalar_type a1 = scalar_type(5)/(dt*dt);
        scalar_type a2 = scalar_type(4)/(dt*dt);
        scalar_type a3 = scalar_type(1)/(dt*dt);
        scalar_type b0 = scalar_type(11)/(scalar_type(6)*dt);
        scalar_type b1 = scalar_type(18)/(scalar_type(6)*dt);
        scalar_type b2 = scalar_type(9)/(scalar_type(6)*dt);
        scalar_type b3 = scalar_type(2)/(scalar_type(6)*dt);

        md.set_factor_of_variable(V, b0);
        md.set_factor_of_variable(A, a0);
        if (md.is_complex()) {
          gmm::add(gmm::scaled(md.complex_variable(U01), -complex_type(b1)),
                   gmm::scaled(md.complex_variable(U02), complex_type(b2)),
                   md.set_complex_constant_part(V));
          gmm::add(gmm::scaled(md.complex_variable(U03), -complex_type(b3)),
                   md.set_complex_constant_part(V));
          gmm::add(gmm::scaled(md.complex_variable(U01), -complex_type(a1)),
                   gmm::scaled(md.complex_variable(U02), complex_type(a2)),
                   md.set_complex_constant_part(A));
          gmm::add(gmm::scaled(md.complex_variable(U03), -complex_type(a3)),
                   md.set_complex_constant_part(A));
        } else {
          gmm::add(gmm::scaled(md.real_variable(U01), -b1),
                   gmm::scaled(md.real_variable(U02), b2),
                   md.set_real_constant_part(V));
          gmm::add(gmm::scaled(md.real_variable(U03), -b3),
                   md.set_real_constant_part(V));
          gmm::add(gmm::scaled(md.real_variable(U01), -a1),
                   gmm::scaled(md.real_variable(U02), a2),
                   md.set_real_constant_part(A));
          gmm::add(gmm::scaled(md.real_variable(U03), -a3),
                   md.set_real_constant_part(A));
        }
      }

      virtual void init_affine_dependent_variables_precomputation(model &md)
        const {
        (void) md;
      }

      virtual void time_derivative_to_be_initialized
      (std::string &name_v, std::string &name_previous_v) const {
        (void) name_v;
        (void) name_previous_v;
      }

      virtual void shift_variables(model &md) const {
        if (md.is_complex()) {
          gmm::copy(md.complex_variable(U02), md.set_complex_variable(U03));
          gmm::copy(md.complex_variable(U01), md.set_complex_variable(U02));
          gmm::copy(md.complex_variable(U), md.set_complex_variable(U01));
        } else {
          gmm::copy(md.real_variable(U02), md.set_real_variable(U03));
          gmm::copy(md.real_variable(U01), md.set_real_variable(U02));
          gmm::copy(md.real_variable(U), md.set_real_variable(U01));
        }
      }


      Houbolt_scheme(model &md, std::string varname) {
        U = varname;
        U01 = "Previous_" + U;
        U02 = "Previous2_" + U;
        U03 = "Previous3_" + U;
        V = "Dot_" + U;
        A = "Dot2_" + U;

        if (!(md.variable_exists(V)))
          md.add_affine_dependent_variable(V, U);
        if (!(md.variable_exists(A)))
          md.add_affine_dependent_variable(A, U);

        const mesh_fem *mf = md.pmesh_fem_of_variable(U);
        size_type s = md.is_complex() ? gmm::vect_size(md.complex_variable(U))
          : gmm::vect_size(md.real_variable(U));

        if (mf) {
          if (!(md.variable_exists(U01))) md.add_fem_data(U01, *mf);
          if (!(md.variable_exists(U02))) md.add_fem_data(U02, *mf);
          if (!(md.variable_exists(U03))) md.add_fem_data(U03, *mf);
        } else {
          if (!(md.variable_exists(U01))) md.add_fixed_size_data(U01, s);
          if (!(md.variable_exists(U02))) md.add_fixed_size_data(U02, s);
          if (!(md.variable_exists(U03))) md.add_fixed_size_data(U03, s);
        }

      }

    };

  void add_Houbolt_scheme(model &md, const std::string &varname) {
    ptime_scheme ptsc = std::make_shared<Houbolt_scheme>
      (md, varname);
    md.add_time_integration_scheme(varname, ptsc);
  }









  void model::add_time_dispatcher(size_type ibrick, pdispatcher pdispatch) {
    GMM_ASSERT1(valid_bricks[ibrick], "Inexistent brick");
    pbrick pbr = bricks[ibrick].pbr;

    bricks[ibrick].pdispatch = pdispatch;

    size_type nbrhs = bricks[ibrick].nbrhs
      = std::max(size_type(1), pdispatch->nbrhs());

    gmm::resize(bricks[ibrick].coeffs, nbrhs);

    if (is_complex() && pbr->is_complex()) {
      bricks[ibrick].cveclist.resize(nbrhs);
      bricks[ibrick].cveclist_sym.resize(nbrhs);
      for (size_type k = 1; k < nbrhs; ++k) {
        bricks[ibrick].cveclist[k] = bricks[ibrick].cveclist[0];
        bricks[ibrick].cveclist_sym[k] = bricks[ibrick].cveclist_sym[0];
      }
    } else {
      bricks[ibrick].rveclist.resize(nbrhs);
      bricks[ibrick].rveclist_sym.resize(nbrhs);
      for (size_type k = 1; k < nbrhs; ++k) {
        bricks[ibrick].rveclist[k] = bricks[ibrick].rveclist[0];
        bricks[ibrick].rveclist_sym[k] = bricks[ibrick].rveclist_sym[0];
      }
    }
  }

  const std::string &model::varname_of_brick(size_type ind_brick,
                                             size_type ind_var) {
    GMM_ASSERT1(valid_bricks[ind_brick], "Inexistent brick");
    GMM_ASSERT1(ind_var < bricks[ind_brick].vlist.size(),
               "Inexistent brick variable");
    return bricks[ind_brick].vlist[ind_var];
  }

  const std::string &model::dataname_of_brick(size_type ind_brick,
                                              size_type ind_data) {
    GMM_ASSERT1(valid_bricks[ind_brick], "Inexistent brick");
    GMM_ASSERT1(ind_data < bricks[ind_brick].dlist.size(),
                "Inexistent brick data");
    return bricks[ind_brick].dlist[ind_data];
  }

  void model::listbricks(std::ostream &ost, size_type base_id) const {
    if (valid_bricks.card() == 0)
      ost << "Model with no bricks" << endl;
    else {
      ost << "List of model bricks:" << endl;
      for (dal::bv_visitor i(valid_bricks); !i.finished(); ++i) {
        ost << "Brick " << std::setw(3) << std::right << i + base_id
            << " " << std::setw(20) << std::right
            << bricks[i].pbr->brick_name();
        if (!(active_bricks[i])) ost << " (deactivated)";
        if (bricks[i].pdispatch) ost << " (dispatched)";
        ost << endl << "  concerned variables: " << bricks[i].vlist[0];
        for (size_type j = 1; j < bricks[i].vlist.size(); ++j)
          ost << ", " << bricks[i].vlist[j];
        ost << "." << endl;
        ost << "  brick with " << bricks[i].tlist.size() << " term";
        if (bricks[i].tlist.size() > 1) ost << "s";
        ost << endl;
        // + lister les termes
      }
    }
  }

  // before call to asm_real_tangent_terms or asm_complex_tangent_terms
  // from the assembly procedure or a time dispatcher
  void model::brick_init(size_type ib, build_version version,
                         size_type rhs_ind) const {
    const brick_description &brick = bricks[ib];
    bool cplx = is_complex() && brick.pbr->is_complex();

    // Initialization of vector and matrices.
    for (size_type j = 0; j < brick.tlist.size(); ++j) {
      const term_description &term = brick.tlist[j];
      bool isg = term.is_global;
      size_type nbgdof = is_complex() ?
        gmm::vect_size(crhs) : gmm::vect_size(rrhs);
      size_type nbd1 = isg ? nbgdof : variables[term.var1].size();
      size_type nbd2 = isg ? nbgdof : (term.is_matrix_term ?
                                      variables[term.var2].size() : 0);
      if (term.is_matrix_term &&
          (brick.pbr->is_linear() || (version | BUILD_MATRIX))) {
        if (version | BUILD_ON_DATA_CHANGE) {
          if (cplx)
            gmm::resize(brick.cmatlist[j], nbd1, nbd2);
          else
            gmm::resize(brick.rmatlist[j], nbd1, nbd2);
        } else {
          if (cplx)
            brick.cmatlist[j] = model_complex_sparse_matrix(nbd1, nbd2);
          else
            brick.rmatlist[j] = model_real_sparse_matrix(nbd1, nbd2);
        }
      }
      if (brick.pbr->is_linear() || (version | BUILD_RHS)) {
        for (size_type k = 0; k < brick.nbrhs; ++k) {
          if (cplx) {
            if (k == rhs_ind) gmm::clear(brick.cveclist[k][j]);
            gmm::resize(brick.cveclist[k][j], nbd1);
            if (term.is_symmetric && term.var1.compare(term.var2)) {
              if (k == rhs_ind) gmm::clear(brick.cveclist_sym[k][j]);
              gmm::resize(brick.cveclist_sym[k][j], nbd2);
            }
          } else {
            if (k == rhs_ind) gmm::clear(brick.rveclist[k][j]);
            gmm::resize(brick.rveclist[k][j], nbd1);
            if (term.is_symmetric && term.var1.compare(term.var2)) {
              if (k == rhs_ind) gmm::clear(brick.rveclist_sym[k][j]);
              gmm::resize(brick.rveclist_sym[k][j], nbd2);
            }
          }
        }
      }
    }
  }

  void model::post_to_variables_step(){}

  void model::brick_call(size_type ib, build_version version,
                         size_type rhs_ind) const
  {
    const brick_description &brick = bricks[ib];
    bool cplx = is_complex() && brick.pbr->is_complex();

    brick_init(ib, version, rhs_ind);

    if (cplx)
    {
      brick.pbr->complex_pre_assembly_in_serial(*this, ib, brick.vlist,
                                                brick.dlist, brick.mims,
                                                brick.cmatlist,
                                                brick.cveclist[rhs_ind],
                                                brick.cveclist_sym[rhs_ind],
                                                brick.region, version);

      /*distributing the resulting vectors and matrices for individual threads.*/
      { //brackets are needed because accumulated_distro has constructor/destructor
        //semantics (as in RAII)
        accumulated_distro<complex_matlist> cmatlist(brick.cmatlist);
        accumulated_distro<complex_veclist> cveclist(brick.cveclist[rhs_ind]);
        accumulated_distro<complex_veclist> cveclist_sym(brick.cveclist_sym[rhs_ind]);

        /*running the assembly in parallel*/
        GETFEM_OMP_PARALLEL(
          brick.pbr->asm_complex_tangent_terms(*this, ib, brick.vlist,
                                                brick.dlist, brick.mims,
                                                cmatlist,
                                                cveclist,
                                                cveclist_sym,
                                                brick.region, version);
        )
      }
      brick.pbr->complex_post_assembly_in_serial(*this, ib, brick.vlist,
                                                 brick.dlist, brick.mims,
                                                 brick.cmatlist,
                                                 brick.cveclist[rhs_ind],
                                                 brick.cveclist_sym[rhs_ind],
                                                 brick.region, version);

      if (brick.is_update_brick) //contributions of pure update bricks must be deleted
      {
        for (auto &&mat : brick.cmatlist)
          gmm::clear(mat);

        for (auto &&vecs : brick.cveclist)
          for (auto &&vec : vecs)
            gmm::clear(vec);

        for (auto &&vecs : brick.cveclist_sym)
          for (auto &&vec : vecs)
            gmm::clear(vec);
      }
    }
    else //not cplx
    {
      brick.pbr->real_pre_assembly_in_serial(*this, ib, brick.vlist,
                                             brick.dlist, brick.mims,
                                             brick.rmatlist,
                                             brick.rveclist[rhs_ind],
                                             brick.rveclist_sym[rhs_ind],
                                             brick.region, version);
      {
        /*distributing the resulting vectors and matrices for individual threads.*/
        accumulated_distro<real_matlist> rmatlist(brick.rmatlist);
        accumulated_distro<real_veclist> rveclist(brick.rveclist[rhs_ind]);
        accumulated_distro<real_veclist> rveclist_sym(brick.rveclist_sym[rhs_ind]);

        /*running the assembly in parallel*/
        GETFEM_OMP_PARALLEL(
            brick.pbr->asm_real_tangent_terms(*this, ib, brick.vlist,
                                              brick.dlist, brick.mims,
                                              rmatlist,
                                              rveclist,
                                              rveclist_sym,
                                              brick.region,
                                              version);
           );
      }
      brick.pbr->real_post_assembly_in_serial(*this, ib, brick.vlist,
                                              brick.dlist, brick.mims,
                                              brick.rmatlist,
                                              brick.rveclist[rhs_ind],
                                              brick.rveclist_sym[rhs_ind],
                                              brick.region, version);

      if (brick.is_update_brick) //contributions of pure update bricks must be deleted
      {
        for (auto &&mat : brick.rmatlist)
          gmm::clear(mat);

        for (auto &&vecs : brick.rveclist)
          for (auto &&vec : vecs)
            gmm::clear(vec);

        for (auto &&vecs : brick.rveclist_sym)
          for (auto &&vec : vecs)
            gmm::clear(vec);
      }
    }
  }


  void model::set_dispatch_coeff() {
    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];
      if (brick.pdispatch)
        brick.pdispatch->set_dispatch_coeff(*this, ib);

    }
  }

  void model::first_iter() {
    context_check(); if (act_size_to_be_done) actualize_sizes();
    for (auto && v : variables) v.second.clear_temporaries();

    set_dispatch_coeff();

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];
      if (brick.pdispatch) {
        if (is_complex() && brick.pbr->is_complex())
          brick.pdispatch->next_complex_iter(*this, ib, brick.vlist,
                                             brick.dlist,
                                             brick.cmatlist, brick.cveclist,
                                             brick.cveclist_sym, true);
        else
          brick.pdispatch->next_real_iter(*this, ib, brick.vlist, brick.dlist,
                                             brick.rmatlist, brick.rveclist,
                                             brick.rveclist_sym, true);
      }
    }
  }

  void model::next_iter() {
    context_check(); if (act_size_to_be_done) actualize_sizes();
    set_dispatch_coeff();

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];
      if (brick.pdispatch) {
        if (is_complex() && brick.pbr->is_complex())
          brick.pdispatch->next_complex_iter(*this, ib, brick.vlist,
                                             brick.dlist,
                                             brick.cmatlist, brick.cveclist,
                                             brick.cveclist_sym, false);
        else
          brick.pdispatch->next_real_iter(*this, ib, brick.vlist, brick.dlist,
                                          brick.rmatlist, brick.rveclist,
                                          brick.rveclist_sym, false);
      }
    }

    for (auto &&v : variables)
      for (size_type i = 1; i < v.second.n_iter; ++i) {
        if (is_complex())
          gmm::copy(v.second.complex_value[i-1], v.second.complex_value[i]);
        else
          gmm::copy(v.second.real_value[i-1], v.second.real_value[i]);
        v.second.v_num_data[i] = act_counter();
      }
  }

  bool model::is_var_newer_than_brick(const std::string &varname,
                                      size_type ib, size_type niter) const {
    const brick_description &brick = bricks[ib];
    var_description &vd = variables[varname];
    if (niter == size_type(-1)) niter = vd.default_iter;
    return (vd.v_num > brick.v_num || vd.v_num_data[niter] > brick.v_num);
  }

  bool model::is_var_mf_newer_than_brick(const std::string &varname,
                                      size_type ib) const {
    const brick_description &brick = bricks[ib];
    var_description &vd = variables[varname];
    return (vd.v_num > brick.v_num);
  }

  bool model::is_mim_newer_than_brick(const mesh_im &im,
                                      size_type ib) const {
    const brick_description &brick = bricks[ib];
    return (im.version_number() > brick.v_num);
  }

  void model::define_variable_group(const std::string &group_name,
                                    const std::vector<std::string> &nl) {
    GMM_ASSERT1(!(variable_exists(group_name)), "The name of a group of "
                "variables cannot be the same as a variable name");

    std::set<const mesh *> ms;
    bool is_data_ = false;
    for (size_type i = 0; i < nl.size(); ++i) {
      if (i == 0)
        is_data_ = is_true_data(nl[i]);
      else {
        GMM_ASSERT1(is_data_ == is_true_data(nl[i]),
                    "It is not possible to mix variables and data in a group");
      }
      GMM_ASSERT1(variable_exists(nl[i]),
                  "All variables in a group have to exist in the model");
      const mesh_fem *mf = pmesh_fem_of_variable(nl[i]);
      GMM_ASSERT1(mf, "Variables in a group should be fem variables");
      GMM_ASSERT1(ms.find(&(mf->linked_mesh())) == ms.end(),
                  "Two variables in a group cannot share the same mesh");
      ms.insert(&(mf->linked_mesh()));
    }
    variable_groups[group_name] = nl;
  }

  void model::add_assembly_assignments(const std::string &varname,
                                       const std::string &expr, size_type rg,
                                       size_type order, bool before) {
    GMM_ASSERT1(order < 3 || order == size_type(-1), "Bad order value");
    const im_data *imd = pim_data_of_variable(varname);
    GMM_ASSERT1(imd != 0, "Only applicable to im_data");
    assignement_desc as;
    as.varname = varname; as.expr = expr; as.region = rg; as.order = order;
    as.before = before;
    assignments.push_back(as);
  }

  void model::add_temporaries(const varnamelist &vl,
                              gmm::uint64_type id_num) const {
    for (size_type i = 0; i < vl.size(); ++i) {
      var_description &vd = variables[vl[i]];
      if (vd.n_iter > 1) {
        vd.add_temporary(id_num);
      }
    }
  }

  bool model::temporary_uptodate(const std::string &varname,
                                 gmm::uint64_type  id_num,
                                 size_type &ind) const {
    var_description &vd = variables[varname];
    ind = vd.n_iter;
    for (; ind < vd.n_iter + vd.n_temp_iter ; ++ind) {
      if (vd.v_num_var_iter[ind] == id_num) break;
    }
    if (ind <  vd.n_iter + vd.n_temp_iter) {
      if (vd.v_num_iter[ind] <= vd.v_num_data[vd.default_iter]) {
        vd.v_num_iter[ind] = act_counter();
        return false;
      }
      return true;
    }
    ind = size_type(-1);
    return true;
  }

  void model::set_default_iter_of_variable(const std::string &varname,
                                    size_type ind) const {
    if (ind != size_type(-1)) {
      var_description &vd = variables[varname];
      GMM_ASSERT1(ind < vd.n_iter + vd.n_temp_iter,
                  "Inexistent iteration " << ind);
      vd.default_iter = ind;
    }
  }

  void model::reset_default_iter_of_variables(const varnamelist &vl) const {
    for (size_type i = 0; i < vl.size(); ++i)
      variables[vl[i]].default_iter = 0;
  }

  const model_real_sparse_matrix &
  model::linear_real_matrix_term(size_type ib, size_type iterm) {
    GMM_ASSERT1(bricks[ib].tlist[iterm].is_matrix_term,
                "Not a matrix term !");
    GMM_ASSERT1(bricks[ib].pbr->is_linear(), "Nonlinear term !");
    return bricks[ib].rmatlist[iterm];
  }

  const model_complex_sparse_matrix &
  model::linear_complex_matrix_term(size_type ib, size_type iterm) {
    GMM_ASSERT1(bricks[ib].tlist[iterm].is_matrix_term,
                "Not a matrix term !");
    GMM_ASSERT1(bricks[ib].pbr->is_linear(), "Nonlinear term !");
    return bricks[ib].cmatlist[iterm];
  }

  // Call the brick to compute the terms
  void model::update_brick(size_type ib, build_version version) const {
    const brick_description &brick = bricks[ib];
    bool cplx = is_complex() && brick.pbr->is_complex();
    bool tobecomputed = brick.terms_to_be_computed
      || brick.pbr->is_to_be_computed_each_time()
      || !(brick.pbr->is_linear());

    // check variable list to test if a mesh_fem has changed.
    if (!tobecomputed ) {
      for (size_type i = 0; i < brick.vlist.size() && !tobecomputed; ++i) {
        var_description &vd = variables[brick.vlist[i]];
        if (vd.v_num > brick.v_num) tobecomputed = true;
      }
    }

    // check data list to test if a vector value of a data has changed.
    for (size_type i = 0; i < brick.dlist.size() && !tobecomputed; ++i) {
      var_description &vd = variables[brick.dlist[i]];
      if (vd.v_num > brick.v_num || vd.v_num_data[vd.default_iter] > brick.v_num) {
        tobecomputed = true;
        version = build_version(version | BUILD_ON_DATA_CHANGE);
      }
    }

    // Check if a mesh_im has changed
    if (!tobecomputed ) {
      for (size_type i = 0; i < brick.mims.size() && !tobecomputed; ++i) {
        if (brick.mims[i]->version_number() > brick.v_num) tobecomputed = true;
      }
    }

    if (tobecomputed) {
      brick.external_load = scalar_type(0);

      if (!(brick.pdispatch))
        { brick_call(ib, version, 0); }
      else {
        if (cplx)
          brick.pdispatch->asm_complex_tangent_terms
            (*this, ib, brick.cmatlist, brick.cveclist, brick.cveclist_sym,
             version);
        else
          brick.pdispatch->asm_real_tangent_terms
            (*this, ib, brick.rmatlist, brick.rveclist, brick.rveclist_sym,
             version);
      }
      brick.v_num = act_counter();
    }

    if (brick.pbr->is_linear()) brick.terms_to_be_computed = false;
  }

  // OBSOLETE (linked to time dispatchers) or to be corrected to take
  // into account global matrices
  void model::linear_brick_add_to_rhs(size_type ib, size_type ind_data,
                                      size_type n_iter) const {
    const brick_description &brick = bricks[ib];
    if (brick.pbr->is_linear()) {

      bool cplx = is_complex() && brick.pbr->is_complex();

      for (size_type j = 0; j < brick.tlist.size(); ++j) {
        const term_description &term = brick.tlist[j];
        bool isg = term.is_global;
        size_type nbgdof = nb_dof();

        size_type n_iter_1 = n_iter, n_iter_2 = n_iter;
        if (!isg && n_iter == size_type(-1)) {
          n_iter_1 = variables[term.var1].default_iter;
          if (term.is_matrix_term)
            n_iter_2 = variables[term.var2].default_iter;
        }



        if (term.is_matrix_term) {
          if (cplx) {
            if (isg) {
              model_complex_plain_vector V(nbgdof);
              for (VAR_SET::iterator it = variables.begin();
                   it != variables.end(); ++it)
                if (it->second.is_variable) {
                  size_type n_iter_i = (n_iter == size_type(-1))
                    ? it->second.default_iter : n_iter;
                  gmm::copy(it->second.complex_value[n_iter_i],
                            gmm::sub_vector(V, it->second.I));
                }
              gmm::mult_add
                (brick.cmatlist[j],
                 gmm::scaled(V,  complex_type(-1)),
                 brick.cveclist[ind_data][j]);
            } else
              gmm::mult_add
                (brick.cmatlist[j],
                 gmm::scaled(variables[term.var2].complex_value[n_iter_2],
                             complex_type(-1)),
                 brick.cveclist[ind_data][j]);
          }
          else {
            if (isg) {
              model_real_plain_vector V(nbgdof);
              for (VAR_SET::iterator it = variables.begin();
                   it != variables.end(); ++it)
                if (it->second.is_variable) {
                  size_type n_iter_i = (n_iter == size_type(-1))
                    ? it->second.default_iter : n_iter;
                  gmm::copy(it->second.real_value[n_iter_i],
                            gmm::sub_vector(V, it->second.I));
                }
              gmm::mult_add
                (brick.rmatlist[j], gmm::scaled(V,  scalar_type(-1)),
                 brick.rveclist[ind_data][j]);
            } else
              gmm::mult_add
                (brick.rmatlist[j],
                 gmm::scaled(variables[term.var2].real_value[n_iter_2],
                             scalar_type(-1)), brick.rveclist[ind_data][j]);
          }

          if (term.is_symmetric && term.var1.compare(term.var2)) {
            if (cplx)
              gmm::mult_add
                (gmm::conjugated(brick.cmatlist[j]),
                 gmm::scaled(variables[term.var1].complex_value[n_iter_1],
                             complex_type(-1)),
                 brick.cveclist_sym[ind_data][j]);
            else
              gmm::mult_add
                (gmm::transposed(brick.rmatlist[j]),
                 gmm::scaled(variables[term.var1].real_value[n_iter_1],
                             scalar_type(-1)),
                 brick.rveclist_sym[ind_data][j]);
          }
        }
      }
    }
  }

  void model::update_affine_dependent_variables() {
    for (VAR_SET::iterator it = variables.begin(); it != variables.end(); ++it)
      if (it->second.is_affine_dependent) {
        VAR_SET::iterator it2 = variables.find(it->second.org_name);
        if (it->second.size() != it2->second.size())
          it->second.set_size();
         if (it->second.is_complex) {
           gmm::add(gmm::scaled(it2->second.complex_value[0],
                                complex_type(it->second.alpha)),
                    it->second.affine_complex_value,
                    it->second.complex_value[0]);
         } else {
           gmm::add(gmm::scaled(it2->second.real_value[0], it->second.alpha),
                    it->second.affine_real_value, it->second.real_value[0]);
         }
        it->second.v_num = std::max(it->second.v_num, it2->second.v_num);
        for (size_type i = 0; i < it->second.n_iter; ++i)
        {
          it->second.v_num_data[i] = std::max(it->second.v_num_data[i],
                                              it2->second.v_num_data[i]);
        }
      }
  }



  std::string model::Neumann_term(const std::string &varname,
                                  size_type region) {
    std::string result;
    const mesh_fem *mf = pmesh_fem_of_variable(varname);
    GMM_ASSERT1(mf, "Works only with fem variables.");
    mesh &m = const_cast<mesh &>(mf->linked_mesh());
    mesh_im dummy_mim(m);

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];

      bool detected = false;
      for (size_type i = 0; i <  brick.vlist.size(); ++i)
        if (brick.vlist[i].compare(varname) == 0)
          { detected = true; break; }

      if (detected && brick.mims.size()) {
        int ifo = -1;
        for (auto &pmim :  brick.mims)
          ifo = std::max(ifo, mf->linked_mesh().region(region)
                              .region_is_faces_of(m, brick.region,
                                                  pmim->linked_mesh()));
        GMM_ASSERT1(ifo >= 0,
                    "The given region is only partially covered by "
                    "region of brick \"" << brick.pbr->brick_name()
                    << "\". Please subdivise the region");
        if (ifo == 1) {
          std::string expr = brick.pbr->declare_volume_assembly_string
            (*this, ib, brick.vlist, brick.dlist);

          ga_workspace workspace(*this);
          size_type order = workspace.add_expression
            (expr, dummy_mim, region);
          GMM_ASSERT1(order <= 1, "Wrong order for a Neumann term");
          expr = workspace.extract_Neumann_term(varname);
          if (expr.size()) {
            if (result.size())
              result += " + " + workspace.extract_Neumann_term(varname);
            else
              result = workspace.extract_Neumann_term(varname);
          }
        }
      }
    }
    return result;
  }



  void model::assembly(build_version version) {

    GMM_ASSERT1(version != BUILD_ON_DATA_CHANGE,
                "Invalid assembly version BUILD_ON_DATA_CHANGE");
    GMM_ASSERT1(version != BUILD_WITH_LIN,
                "Invalid assembly version BUILD_WITH_LIN");
    GMM_ASSERT1(version != BUILD_WITH_INTERNAL,
                "Invalid assembly version BUILD_WITH_INTERNAL");
    int nbp=1;
#if GETFEM_PARA_LEVEL > 0
    double t_ref = MPI_Wtime();
    int rk=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &nbp);
#endif

    context_check(); if (act_size_to_be_done) actualize_sizes();
    if (is_complex()) {
      if (version & BUILD_MATRIX) gmm::clear(cTM);
      if (version & BUILD_RHS) gmm::clear(crhs);
    }
    else {
      if (version & BUILD_MATRIX) gmm::clear(rTM);
      if (version & BUILD_RHS) gmm::clear(rrhs);
    }
    clear_dof_constraints();
    generic_expressions.clear();
    update_affine_dependent_variables();

    if (version & BUILD_RHS) approx_external_load_ = scalar_type(0);

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {

      brick_description &brick = bricks[ib];

      // Disables the brick if all its variables are disabled.
      bool auto_disabled_brick = true;
      for (size_type j = 0; j < brick.vlist.size(); ++j) {
        if (!(is_disabled_variable(brick.vlist[j])))
          auto_disabled_brick = false;
      }
      if (auto_disabled_brick) continue;

      update_brick(ib, version);

      bool cplx = is_complex() && brick.pbr->is_complex();

      scalar_type coeff0 = scalar_type(1);
      if (brick.pdispatch) coeff0 = brick.matrix_coeff;

      // Assembly of terms

      for (size_type j = 0; j < brick.tlist.size(); ++j) {
        term_description &term = brick.tlist[j];
        bool isg = term.is_global, isprevious = false;
        size_type nbgdof = nb_dof();
        scalar_type alpha = coeff0, alpha1 = coeff0, alpha2 = coeff0;
        gmm::sub_interval I1(0,nbgdof), I2(0,nbgdof);
        var_description *var1=nullptr, *var2=nullptr;
        if (!isg) {
          VAR_SET::iterator it1 = variables.find(term.var1);
          GMM_ASSERT1(it1 != variables.end(), "Internal error");
          var1 = &(it1->second);
          GMM_ASSERT1(var1->is_variable, "Assembly of data not allowed");
          I1 = var1->I;
          if (term.is_matrix_term) {
            VAR_SET::iterator it2 = variables.find(term.var2);
            GMM_ASSERT1(it2 != variables.end(), "Internal error");
            var2 = &(it2->second);
            I2 = var2->I;
            if (!(var2->is_variable)) {
              std::string vorgname = sup_previous_and_dot_to_varname(term.var2);
              VAR_SET::iterator it3 = variables.find(vorgname);
              GMM_ASSERT1(it3->second.is_variable,
                          "Assembly of data not allowed");
              I2 = it3->second.I;
              isprevious = true;
            }
            alpha *= var1->alpha * var2->alpha;
            alpha1 *= var1->alpha;
            alpha2 *= var2->alpha;
          }
        }

        if (cplx) { // complex term in complex model
          if (term.is_matrix_term && (version & BUILD_MATRIX) && !isprevious
              && (isg || (var1->is_enabled() && var2->is_enabled()))) {
            gmm::add(gmm::scaled(brick.cmatlist[j], alpha),
                     gmm::sub_matrix(cTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first())
              gmm::add(gmm::scaled(gmm::transposed(brick.cmatlist[j]), alpha),
                       gmm::sub_matrix(cTM, I2, I1));
          }
          if (version & BUILD_RHS) {
            //FIXME MPI_SUM_VECTOR(crhs)
            if (isg || var1->is_enabled()) {
              if (brick.pdispatch)
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.cveclist[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I1));
              else
                gmm::add(gmm::scaled(brick.cveclist[0][j],
                                     complex_type(alpha1)),
                         gmm::sub_vector(crhs, I1));
            }
            if (term.is_matrix_term && brick.pbr->is_linear() && is_linear()) {
              if (var2->is_affine_dependent && var1->is_enabled())
                gmm::mult_add(brick.cmatlist[j],
                              gmm::scaled(var2->affine_complex_value,
                                          complex_type(-alpha1)),
                              gmm::sub_vector(crhs, I1));
              if (term.is_symmetric && I1.first() != I2.first()
                  && var1->is_affine_dependent && var2->is_enabled())
                gmm::mult_add(gmm::conjugated(brick.cmatlist[j]),
                              gmm::scaled(var1->affine_complex_value,
                                          complex_type(-alpha2)),
                              gmm::sub_vector(crhs, I2));
            }
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_LIN))
                && var1->is_enabled())
              gmm::mult_add(brick.cmatlist[j],
                            gmm::scaled(var2->complex_value[0],
                                        complex_type(-alpha1)),
                            gmm::sub_vector(crhs, I1));
            if (term.is_symmetric && I1.first() != I2.first()
                && var2->is_enabled()) {
              if (brick.pdispatch)
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.cveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I2));
              else
                gmm::add(gmm::scaled(brick.cveclist_sym[0][j],
                                     complex_type(alpha2)),
                         gmm::sub_vector(crhs, I2));
              if (brick.pbr->is_linear()
                  && (!is_linear() || (version & BUILD_WITH_LIN)))
                 gmm::mult_add(gmm::conjugated(brick.cmatlist[j]),
                               gmm::scaled(var1->complex_value[0],
                                           complex_type(-alpha2)),
                               gmm::sub_vector(crhs, I2));
            }
          }
        } else if (is_complex()) { // real term in complex model
          if (term.is_matrix_term && (version & BUILD_MATRIX) && !isprevious
              && (isg || (var1->is_enabled() && var2->is_enabled()))) {
            gmm::add(gmm::scaled(brick.rmatlist[j], alpha),
                     gmm::sub_matrix(cTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first())
              gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), alpha),
                       gmm::sub_matrix(cTM, I2, I1));
          }
          if (version & BUILD_RHS) {
            //FIXME MPI_SUM_VECTOR(crhs)
            if (isg || var1->is_enabled()) {
              if (brick.pdispatch)
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I1));
              else
                gmm::add(gmm::scaled(brick.rveclist[0][j], alpha1),
                         gmm::sub_vector(crhs, I1));
            }
            if (term.is_matrix_term && brick.pbr->is_linear() && is_linear()) {
              if (var2->is_affine_dependent && var1->is_enabled())
                gmm::mult_add(brick.rmatlist[j],
                              gmm::scaled(var2->affine_complex_value,
                                          complex_type(-alpha1)),
                              gmm::sub_vector(crhs, I1));
              if (term.is_symmetric && I1.first() != I2.first()
                  && var1->is_affine_dependent && var2->is_enabled())
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                              gmm::scaled(var1->affine_complex_value,
                                          complex_type(-alpha2)),
                              gmm::sub_vector(crhs, I2));
            }
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_LIN))
                && var1->is_enabled())
              gmm::mult_add(brick.rmatlist[j],
                            gmm::scaled(var2->complex_value[0],
                                        complex_type(-alpha1)),
                            gmm::sub_vector(crhs, I1));
            if (term.is_symmetric && I1.first() != I2.first()
                && var2->is_enabled()) {
              if (brick.pdispatch)
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I2));
              else
                gmm::add(gmm::scaled(brick.rveclist_sym[0][j], alpha2),
                         gmm::sub_vector(crhs, I2));

              if (brick.pbr->is_linear()
                  && (!is_linear() || (version & BUILD_WITH_LIN)))
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                              gmm::scaled(var1->complex_value[0],
                                          complex_type(-alpha2)),
                              gmm::sub_vector(crhs, I2));
            }
          }
        } else { // real term in real model
          if (term.is_matrix_term && (version & BUILD_MATRIX) && !isprevious
              && (isg || (var1->is_enabled() && var2->is_enabled()))) {
            gmm::add(gmm::scaled(brick.rmatlist[j], alpha),
                     gmm::sub_matrix(rTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first())
              gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), alpha),
                       gmm::sub_matrix(rTM, I2, I1));
          }
          if (version & BUILD_RHS) {
            // Contributions to interval I1 of var1
            auto vec_out1 = gmm::sub_vector(rrhs, I1);
            if (isg || var1->is_enabled()) {
              if (brick.pdispatch)
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist[k][j],
                                       brick.coeffs[k]),
                           vec_out1);
              else
                gmm::add(gmm::scaled(brick.rveclist[0][j], alpha1),
                         vec_out1);
            }
            if (var1->is_enabled()
                && term.is_matrix_term && brick.pbr->is_linear()) {
              bool affine_contrib(is_linear() && var2->is_affine_dependent);
              bool linear_contrib(!is_linear() || (version & BUILD_WITH_LIN));
              const auto &matj = brick.rmatlist[j];
              const auto vec_affine2 = gmm::scaled(var2->affine_real_value,
                                                   -alpha1);
              const auto vec_linear2 = gmm::scaled(var2->real_value[0],
                                                   -alpha1);
              if (nbp > 1) {
                model_real_plain_vector vec_tmp1(I1.size(), 0.);
                if (affine_contrib) // Affine dependent variable contribution
                  gmm::mult(matj, vec_affine2, vec_tmp1);
                if (linear_contrib) // Linear term contribution
                  gmm::mult_add(matj, vec_linear2, vec_tmp1);
                MPI_SUM_VECTOR(vec_tmp1);
                gmm::add(vec_tmp1, vec_out1);
              } else { // nbp == 1
                if (affine_contrib) // Affine dependent variable contribution
                  gmm::mult_add(matj, vec_affine2, vec_out1);
                if (linear_contrib) // Linear term contribution
                  gmm::mult_add(matj, vec_linear2, vec_out1);
              }
            }

            // Contributions to interval I2 of var2 due to symmetric terms
            if (term.is_symmetric && I1.first() != I2.first() &&
                var2->is_enabled()) {
              auto vec_out2 = gmm::sub_vector(rrhs, I2);
              if (brick.pdispatch)
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist_sym[k][j],
                                       brick.coeffs[k]),
                           vec_out2);
              else
                gmm::add(gmm::scaled(brick.rveclist_sym[0][j], alpha2),
                         vec_out2);
              if (term.is_matrix_term && brick.pbr->is_linear()) {
                bool affine_contrib(is_linear() && var1->is_affine_dependent);
                bool linear_contrib(!is_linear() || (version & BUILD_WITH_LIN));
                const auto matj_trans = gmm::transposed(brick.rmatlist[j]);
                const auto vec_affine1 = gmm::scaled(var1->affine_real_value,
                                                     -alpha2);
                const auto vec_linear1 = gmm::scaled(var1->real_value[0],
                                                     -alpha2);
                if (nbp > 1) {
                  model_real_plain_vector vec_tmp2(I2.size(),0.);
                  if (affine_contrib) // Affine dependent variable contribution
                    gmm::mult(matj_trans, vec_affine1, vec_tmp2);
                  if (linear_contrib) // Linear term contribution
                    gmm::mult_add(matj_trans, vec_linear1, vec_tmp2);
                  MPI_SUM_VECTOR(vec_tmp2);
                  gmm::add(vec_tmp2, vec_out2);
                } else { // nbp == 1
                  if (affine_contrib) // Affine dependent variable contribution
                    gmm::mult_add(matj_trans, vec_affine1, vec_out2);
                  if (linear_contrib) // Linear term contribution
                    gmm::mult_add(matj_trans, vec_linear1, vec_out2);
                }
              }
            }
          }
        }
      }

      if (brick.pbr->is_linear())
        brick.terms_to_be_computed = false;
      // Commented to allow to get the information after assembly. Used in
      // some aplications. Should be optional ?
//       else
//         if (cplx) {
//           brick.cmatlist = complex_matlist(brick.tlist.size());
//           brick.cveclist[0] = complex_veclist(brick.tlist.size());
//         } else {
//           brick.rmatlist = real_matlist(brick.tlist.size());
//           brick.rveclist[0] = real_veclist(brick.tlist.size());
//         }

      if (version & BUILD_RHS) approx_external_load_ += brick.external_load;
    }

    if (version & BUILD_RHS && version & BUILD_WITH_INTERNAL) {
      GMM_ASSERT1(gmm::vect_size(full_rrhs) > 0 && has_internal_variables(),
                  "Internal error");
      gmm::sub_interval IP(0,gmm::vect_size(rrhs));
      gmm::fill(full_rrhs, 0.);
      gmm::copy(rrhs, gmm::sub_vector(full_rrhs, IP)); // TICTIC
    }

    // Generic expressions
    if (generic_expressions.size()) {
      GMM_ASSERT1(!is_complex(), "to be done");

      if (version & BUILD_RHS)
        GMM_TRACE2("Global generic assembly RHS");
      if (version & BUILD_MATRIX)
        GMM_TRACE2("Global generic assembly tangent term");

      // auxilliary lambda function
      auto add_assignments_and_expressions_to_workspace =
      [&](ga_workspace &workspace)
      {
        for (const auto &ad : assignments)
          workspace.add_assignment_expression
            (ad.varname, ad.expr, ad.region, ad.order, ad.before);
        for (const auto &ge : generic_expressions)
          workspace.add_expression(ge.expr, ge.mim, ge.region,
                                   2, ge.secondary_domain);
      };

      const bool with_internal = version & BUILD_WITH_INTERNAL
                                 && has_internal_variables();
      model_real_sparse_matrix intern_mat; // temp for extracting condensation info
      model_real_plain_vector res0, // holds the original RHS
                              res1; // holds the condensed RHS

      size_type full_size = gmm::vect_size(full_rrhs),
                primary_size = gmm::vect_size(rrhs);

      if ((version & BUILD_RHS) || (version & BUILD_MATRIX && with_internal))
        gmm::resize(res0, with_internal ? full_size : primary_size);
      if (version & BUILD_MATRIX && with_internal)
        gmm::resize(res1, full_size);

      if (version & BUILD_MATRIX) {
        if (with_internal) {
          gmm::resize(intern_mat, full_size, primary_size);
          gmm::resize(res1, full_size);
        }
        accumulated_distro<decltype(rTM)> tangent_matrix_distro(rTM);
        accumulated_distro<decltype(intern_mat)> intern_mat_distro(intern_mat);
        accumulated_distro<model_real_plain_vector> res1_distro(res1);

        if (version & BUILD_RHS) { // both BUILD_RHS & BUILD_MATRIX
          accumulated_distro<model_real_plain_vector> res0_distro(res0);
          GETFEM_OMP_PARALLEL( // running the assembly in parallel
            ga_workspace workspace(*this);
            add_assignments_and_expressions_to_workspace(workspace);
            workspace.set_assembled_vector(res0_distro);
            workspace.assembly(1, with_internal);
            if (with_internal) { // Condensation reads from/writes to rhs
              gmm::copy(res0_distro.get(), res1_distro.get());
              gmm::add(gmm::scaled(full_rrhs, scalar_type(-1)),
                       res1_distro.get()); // initial value residual=-rhs (actually only the internal variables residual is needed)
              workspace.set_assembled_vector(res1_distro);
              workspace.set_internal_coupling_matrix(intern_mat_distro);
            }
            workspace.set_assembled_matrix(tangent_matrix_distro);
            workspace.assembly(2, with_internal);
          ) // end GETFEM_OMP_PARALLEL
        } // end of res0_distro scope
        else { // only BUILD_MATRIX
          GETFEM_OMP_PARALLEL( // running the assembly in parallel
            ga_workspace workspace(*this);
            add_assignments_and_expressions_to_workspace(workspace);
            if (with_internal) { // Condensation reads from/writes to rhs
              gmm::copy(gmm::scaled(full_rrhs, scalar_type(-1)),
                        res1_distro.get()); // initial value residual=-rhs (actually only the internal variables residual is needed)
              workspace.set_assembled_vector(res1_distro);
              workspace.set_internal_coupling_matrix(intern_mat_distro);
            }
            workspace.set_assembled_matrix(tangent_matrix_distro);
            workspace.assembly(2, with_internal);
          ) // end GETFEM_OMP_PARALLEL
        }
      } // end of tangent_matrix_distro, intern_mat_distro, res1_distro scope
      else if (version & BUILD_RHS) {
        accumulated_distro<model_real_plain_vector> res0_distro(res0);
        GETFEM_OMP_PARALLEL( // running the assembly in parallel
          ga_workspace workspace(*this);
          add_assignments_and_expressions_to_workspace(workspace);
          workspace.set_assembled_vector(res0_distro);
          workspace.assembly(1, with_internal);
        ) // end GETFEM_OMP_PARALLEL
      } // end of res0_distro scope

      if (version & BUILD_RHS) {
        gmm::scale(res0, scalar_type(-1)); // from residual to rhs
        if (with_internal) {
          gmm::sub_interval IP(0,gmm::vect_size(rrhs));
          gmm::add(gmm::sub_vector(res0, IP), rrhs); // TOCTOC
          gmm::add(res0, full_rrhs);
        } else
          gmm::add(res0, rrhs);
      }

      if (version & BUILD_MATRIX && with_internal) {
        gmm::scale(res1, scalar_type(-1)); // from residual to rhs
        gmm::sub_interval IP(0, primary_size),
                          II(primary_size, full_size-primary_size);
        gmm::copy(gmm::sub_matrix(intern_mat, II, IP), internal_rTM); // --> internal_rTM
        gmm::add(gmm::sub_vector(res1, IP), rrhs);                    // --> rrhs
        gmm::copy(gmm::sub_vector(res1, II), internal_sol);           // --> internal_sol
      }
    }

    // Post simplification for dof constraints
    if ((version & BUILD_RHS) || (version & BUILD_MATRIX)) {
      if (is_complex()) {
        std::vector<size_type> dof_indices;
        std::vector<complex_type> dof_pr_values;
        std::vector<complex_type> dof_go_values;
        for (const auto &keyval : complex_dof_constraints) {
          const gmm::sub_interval &I = interval_of_variable(keyval.first);
          const model_complex_plain_vector &V = complex_variable(keyval.first);
          for (const auto &val : keyval.second) {
            dof_indices.push_back(val.first + I.first());
            dof_go_values.push_back(val.second);
            dof_pr_values.push_back(V[val.first]);
          }
        }

        if (dof_indices.size()) {
          gmm::sub_index SI(dof_indices);
          gmm::sub_interval II(0, nb_dof());

          if (version & BUILD_RHS) {
            if (MPI_IS_MASTER())
              approx_external_load_ += gmm::vect_norm1(dof_go_values);
            if (is_linear_) {
              if (is_symmetric_) {
                scalar_type valnorm = gmm::vect_norm2(dof_go_values);
                if (valnorm > scalar_type(0)) {
                  GMM_ASSERT1(version & BUILD_MATRIX, "Rhs only for a "
                              "symmetric linear problem with dof "
                              "constraint not allowed");
                  model_complex_plain_vector vv(gmm::vect_size(crhs));
                  gmm::mult(gmm::sub_matrix(cTM, II, SI), dof_go_values, vv);
                  MPI_SUM_VECTOR(vv);
                  gmm::add(gmm::scaled(vv, scalar_type(-1)), crhs);
                }
              }
              gmm::copy(dof_go_values, gmm::sub_vector(crhs, SI));
            } else {
              gmm::add(dof_go_values,
                       gmm::scaled(dof_pr_values, complex_type(-1)),
                       gmm::sub_vector(crhs, SI));
            }
          }
          if (version & BUILD_MATRIX) {
            gmm::clear(gmm::sub_matrix(cTM, SI, II));
            if (is_symmetric_) gmm::clear(gmm::sub_matrix(cTM, II, SI));

            if (MPI_IS_MASTER()) {
              for (size_type i = 0; i < dof_indices.size(); ++i)
                cTM(dof_indices[i], dof_indices[i]) = complex_type(1);
            }
          }
        }
      } else { // !is_complex()
        std::vector<size_type> dof_indices;
        std::vector<scalar_type> dof_pr_values;
        std::vector<scalar_type> dof_go_values;
        for (const auto &keyval : real_dof_constraints) {
          const gmm::sub_interval &I = interval_of_variable(keyval.first);
          const model_real_plain_vector &V = real_variable(keyval.first);
          for (const auto &val : keyval.second) {
            dof_indices.push_back(val.first + I.first());
            dof_go_values.push_back(val.second);
            dof_pr_values.push_back(V[val.first]);
          }
        }

        #if GETFEM_PARA_LEVEL > 1
        GMM_ASSERT1(MPI_IS_MASTER() || (dof_indices.size() == 0),
                    "Sorry, for the moment, the dof constraints have to be "
                    "added on the master process only");
        size_type dof_indices_size = dof_indices.size();
        MPI_BCAST0_SCALAR(dof_indices_size);
        dof_indices.resize(dof_indices_size);
        MPI_BCAST0_VECTOR(dof_indices);
        dof_pr_values.resize(dof_indices_size);
        MPI_BCAST0_VECTOR(dof_pr_values);
        dof_go_values.resize(dof_indices_size);
        MPI_BCAST0_VECTOR(dof_go_values);
        #endif

        if (dof_indices.size()) {
          gmm::sub_index SI(dof_indices);
          gmm::sub_interval II(0, nb_dof());

          if (version & BUILD_RHS) {
            if (MPI_IS_MASTER())
              approx_external_load_ += gmm::vect_norm1(dof_go_values);
            if (is_linear_) {
              if (is_symmetric_) {
                scalar_type valnorm = gmm::vect_norm2(dof_go_values);
                if (valnorm > scalar_type(0)) {
                  GMM_ASSERT1(version & BUILD_MATRIX, "Rhs only for a "
                              "symmetric linear problem with dof "
                              "constraint not allowed");
                  model_real_plain_vector vv(gmm::vect_size(rrhs));
                  gmm::mult(gmm::sub_matrix(rTM, II, SI), dof_go_values, vv);
                  MPI_SUM_VECTOR(vv);
                  gmm::add(gmm::scaled(vv, scalar_type(-1)), rrhs);
                }
              }
              gmm::copy(dof_go_values, gmm::sub_vector(rrhs, SI));
            } else {
              gmm::add(dof_go_values,
                       gmm::scaled(dof_pr_values, scalar_type(-1)),
                       gmm::sub_vector(rrhs, SI));
            }
          }
          if (version & BUILD_MATRIX) {
            gmm::clear(gmm::sub_matrix(rTM, SI, II));
            if (is_symmetric_) gmm::clear(gmm::sub_matrix(rTM, II, SI));

            if (MPI_IS_MASTER()) {
              for (size_type i = 0; i < dof_indices.size(); ++i)
                rTM(dof_indices[i], dof_indices[i]) = scalar_type(1);
            }
          }
        }
      }
    }

    if (version & BUILD_RHS) {
       // some contributions are added only in the master process
       // send the correct result to all other processes
       MPI_BCAST0_SCALAR(approx_external_load_);
    }

    #if GETFEM_PARA_LEVEL > 0
    if (MPI_IS_MASTER()) cout << "Assembly time " << MPI_Wtime()-t_ref << endl;
    #endif

  }


  const mesh_fem &
  model::mesh_fem_of_variable(const std::string &name) const {
    auto it = find_variable(no_old_prefix_name(name));
    return it->second.associated_mf();
  }

  const mesh_fem *
  model::pmesh_fem_of_variable(const std::string &name) const {
    auto it = find_variable(no_old_prefix_name(name));
    return it->second.passociated_mf();
  }

  bgeot::multi_index
  model::qdims_of_variable(const std::string &name) const {
    auto it = find_variable(no_old_prefix_name(name));
    const mesh_fem *mf = it->second.passociated_mf();
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

  size_type model::qdim_of_variable(const std::string &name) const {
    auto it = find_variable(no_old_prefix_name(name));
    const mesh_fem *mf = it->second.passociated_mf();
    const im_data *imd = it->second.imd;
    size_type n = it->second.qdim();
    if (mf) {
      return mf->get_qdim() * n;
    } else if (imd) {
      return imd->tensor_size().total_size() * n;
    }
    return n;
  }


  const gmm::sub_interval &
  model::interval_of_variable(const std::string &name) const {
    context_check();
    if (act_size_to_be_done) actualize_sizes();
    VAR_SET::const_iterator it = find_variable(name);
    return it->second.I;
  }

  const model_real_plain_vector &
  model::real_variable(const std::string &name) const {
    return is_old(name) ? real_variable(no_old_prefix_name(name), 1)
                        : real_variable(name, size_type(-1));
  }

  const model_real_plain_vector &
  model::real_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(!complex_version, "This model is a complex one");
    GMM_ASSERT1(!is_old(name), "Please don't use Old_ prefix in combination "
                               "with variable version");
    context_check();
    auto it = variables.find(name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.mf) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size();
    }
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter > niter,
                "Invalid iteration number " << niter << " for " << name);
    return it->second.real_value[niter];
  }

  const model_complex_plain_vector &
  model::complex_variable(const std::string &name) const {
    return is_old(name) ? complex_variable(no_old_prefix_name(name), 1)
                        : complex_variable(name, size_type(-1));
  }

  const model_complex_plain_vector &
  model::complex_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    GMM_ASSERT1(!is_old(name), "Please don't use Old_ prefix in combination with"
                               " variable version");
    context_check();
    auto it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.mf) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size();
    }
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter  > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.complex_value[niter];
  }

  model_real_plain_vector &
  model::set_real_variable(const std::string &name) const {
    return is_old(name) ? set_real_variable(no_old_prefix_name(name), 1)
                        : set_real_variable(name, size_type(-1));
  }


  model_real_plain_vector &
  model::set_real_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(!complex_version, "This model is a complex one");
    GMM_ASSERT1(!is_old(name), "Please don't use Old_ prefix in combination with"
                               " variable version");
    context_check();
    auto it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.mf) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size();
    }
    if (niter == size_type(-1)) niter = it->second.default_iter;
    it->second.v_num_data[niter] = act_counter();
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.real_value[niter];
  }

  model_complex_plain_vector &
  model::set_complex_variable(const std::string &name) const {
    return is_old(name) ? set_complex_variable(no_old_prefix_name(name), 1)
                        : set_complex_variable(name, size_type(-1));
  }

  model_complex_plain_vector &
  model::set_complex_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    GMM_ASSERT1(!is_old(name), "Please don't use Old_ prefix in combination with"
                               " variable version");
    context_check();
    auto it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.mf) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size();
    }
    if (niter == size_type(-1)) niter = it->second.default_iter;
    it->second.v_num_data[niter] = act_counter();
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.complex_value[niter];
  }

  model_real_plain_vector &
  model::set_real_constant_part(const std::string &name) const {
    GMM_ASSERT1(!complex_version, "This model is a complex one");
    context_check();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
    GMM_ASSERT1(it->second.is_affine_dependent,
                "Only for affine dependent variables");
    if (act_size_to_be_done && it->second.mf) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size();
    }
    for (auto &v_num : it->second.v_num_data) v_num = act_counter();
    return it->second.affine_real_value;
  }

  model_complex_plain_vector &
  model::set_complex_constant_part(const std::string &name) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    context_check();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.mf) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size();
    }
    for (auto &v_num : it->second.v_num_data) v_num = act_counter();
    return it->second.affine_complex_value;
  }

  void model::check_brick_stiffness_rhs(size_type ind_brick) const {

    const brick_description &brick = bricks[ind_brick];
    update_brick(ind_brick, model::BUILD_ALL);

    brick.pbr->check_stiffness_matrix_and_rhs(*this, ind_brick, brick.tlist,
                      brick.vlist, brick.dlist, brick.mims, brick.rmatlist,
                      brick.rveclist[0], brick.rveclist_sym[0], brick.region);
  }

  void model::clear() {
    variables.clear();
    active_bricks.clear();
    valid_bricks.clear();
    real_dof_constraints.clear();
    complex_dof_constraints.clear();
    bricks.resize(0);
    rTM = model_real_sparse_matrix();
    cTM = model_complex_sparse_matrix();
    rrhs = model_real_plain_vector();
    crhs = model_complex_plain_vector();
  }



  void virtual_brick::full_asm_real_tangent_terms_(const model &md, size_type ind_brick,
    const model::varnamelist &term_list,
    const model::varnamelist &var_list,
    const model::mimlist &mim_list,
    model::real_matlist &rmatlist,
    model::real_veclist &rveclist,
    model::real_veclist &rveclist_sym,
    size_type region, build_version version) const
  {
    real_pre_assembly_in_serial(md, ind_brick, term_list, var_list, mim_list, rmatlist,
      rveclist, rveclist_sym, region, version);
    asm_real_tangent_terms(md, ind_brick, term_list, var_list, mim_list, rmatlist,
      rveclist, rveclist_sym, region, version);
    real_post_assembly_in_serial(md, ind_brick, term_list, var_list, mim_list, rmatlist,
      rveclist, rveclist_sym, region, version);
  }

  void virtual_brick::check_stiffness_matrix_and_rhs
  (const model &md, size_type s,
   const model::termlist& tlist,
   const model::varnamelist &vl,
   const model::varnamelist &dl,
   const model::mimlist &mims,
   model::real_matlist &matl,
   model::real_veclist &rvc1,
   model::real_veclist &rvc2,
   size_type rg,
   const scalar_type TINY) const
  {
    std::cout<<"******Verifying stiffnesses of *******"<<std::endl;
    std::cout<<"*** "<<brick_name()<<std::endl;

    //Build the index for the corresponding RHS
    std::map<std::string,size_type> rhs_index;
    for(size_type iterm=0;iterm<matl.size();iterm++)
      if (tlist[iterm].var1==tlist[iterm].var2) rhs_index[tlist[iterm].var1]=iterm;

    if (rhs_index.size()==0){
      GMM_WARNING0("*** cannot verify stiffness for this brick***");
      return;
    }
    full_asm_real_tangent_terms_(md, s, vl, dl, mims, matl, rvc1, rvc2,
                           rg, model::BUILD_MATRIX);
    for(size_type iterm=0;iterm<matl.size();iterm++)
    {

      std::cout<<std::endl;
      std::cout<<"    Stiffness["<<tlist[iterm].var1
               <<","<<tlist[iterm].var2<<"]:"<<std::endl;
      if (md.real_variable(tlist[iterm].var1).size()==0)
        {
          std::cout<<"    "<<tlist[iterm].var1<<" has zero size. Skipping this term"<<std::endl;
          continue;
        }
      if (md.real_variable(tlist[iterm].var2).size()==0)
        {
          std::cout<<"    "<<tlist[iterm].var2<<" has zero size. Skipping this term"<<std::endl;
          continue;
        }

      model_real_sparse_matrix SM(matl[iterm]);
      gmm::fill(rvc1[rhs_index[tlist[iterm].var1]], 0.0);
      full_asm_real_tangent_terms_(md, s, vl, dl, mims, matl, rvc1, rvc2,
                             rg, model::BUILD_RHS);
      if (gmm::mat_euclidean_norm(matl[iterm])<1e-12){
        std::cout<<"    The assembled matrix is nearly zero, skipping."<<std::endl;
        continue;
      }
      model_real_plain_vector RHS0(rvc1[rhs_index[tlist[iterm].var1]]);

      //finite difference stiffness
      model_real_sparse_matrix fdSM(matl[iterm].nrows(), matl[iterm].ncols());
      model_real_plain_vector&U = md.set_real_variable(tlist[iterm].var2);
      model_real_plain_vector& RHS1 = rvc1[rhs_index[tlist[iterm].var1]];

      scalar_type relative_tiny = gmm::vect_norminf(RHS1)*TINY;
      if (relative_tiny < TINY) relative_tiny = TINY;

      for (size_type j = 0; j < matl[iterm].ncols(); j++){
        U[j] += relative_tiny;
        gmm::fill(RHS1, 0.0);
        full_asm_real_tangent_terms_(md, s, vl, dl, mims, matl, rvc1, rvc2,
          rg, model::BUILD_RHS);
        for (size_type i = 0; i<matl[iterm].nrows(); i++)
          fdSM(i, j) = (RHS0[i] - RHS1[i]) / relative_tiny;
        U[j] -= relative_tiny;
      }

      model_real_sparse_matrix diffSM(matl[iterm].nrows(),matl[iterm].ncols());
      gmm::add(SM,gmm::scaled(fdSM,-1.0),diffSM);
      scalar_type norm_error_euc
        = gmm::mat_euclidean_norm(diffSM)/gmm::mat_euclidean_norm(SM)*100;
      scalar_type norm_error_1
        = gmm::mat_norm1(diffSM)/gmm::mat_norm1(SM)*100;
      scalar_type norm_error_max
        = gmm::mat_maxnorm(diffSM)/gmm::mat_maxnorm(SM)*100;

      //checking symmetry of diagonal terms
      scalar_type nsym_norm_error_euc=0.0;
      scalar_type nsym_norm_error_1=0.0;
      scalar_type nsym_norm_error_max=0.0;
      if (tlist[iterm].var1==tlist[iterm].var2){
        model_real_sparse_matrix diffSMtransposed(matl[iterm].nrows(),matl[iterm].ncols());
        gmm::add(gmm::transposed(fdSM),gmm::scaled(fdSM,-1.0),diffSMtransposed);
        nsym_norm_error_euc
          = gmm::mat_euclidean_norm(diffSMtransposed)/gmm::mat_euclidean_norm(fdSM)*100;
        nsym_norm_error_1
          = gmm::mat_norm1(diffSMtransposed)/gmm::mat_norm1(fdSM)*100;
        nsym_norm_error_max
          = gmm::mat_maxnorm(diffSMtransposed)/gmm::mat_maxnorm(fdSM)*100;
      }

      //print matrix if the size is small
      if(rvc1[0].size()<8){
        std::cout << "RHS Stiffness Matrix: \n";
        std::cout << "------------------------\n";
        for(size_type i=0; i < rvc1[iterm].size(); ++i){
          std::cout << "[";
          for(size_type j=0; j < rvc1[iterm].size(); ++j){
            std::cout << fdSM(i,j) << "  ";
          }
          std::cout << "]\n";
        }
        std::cout << "Analytical Stiffness Matrix: \n";
        std::cout << "------------------------\n";
        for(size_type i=0; i < rvc1[iterm].size(); ++i){
          std::cout << "[";
          for(size_type j=0; j < rvc1[iterm].size(); ++j){
            std::cout << matl[iterm](i,j) << "  ";
          }
          std::cout << "]\n";
        }
        std::cout << "Vector U: \n";
        std::cout << "------------------------\n";
        for(size_type i=0; i < rvc1[iterm].size(); ++i){
          std::cout << "[";
          std::cout << md.real_variable(tlist[iterm].var2)[i] << "  ";
          std::cout << "]\n";
        }
      }
      std::cout
        << "\n\nfinite diff test error_norm_eucl: " << norm_error_euc << "%\n"
        << "finite diff test error_norm1: " << norm_error_1 << "%\n"
        << "finite diff test error_max_norm: " << norm_error_max << "%\n\n\n";

      if (tlist[iterm].var1==tlist[iterm].var2){
        std::cout
          << "Nonsymmetrical test error_norm_eucl: "<< nsym_norm_error_euc<< "%\n"
          << "Nonsymmetrical test error_norm1: " << nsym_norm_error_1 << "%\n"
          << "Nonsymmetrical test error_max_norm: " << nsym_norm_error_max << "%"
          << std::endl;
      }
    }
  }

  // ----------------------------------------------------------------------
  //
  //
  // Standard bricks
  //
  //
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  //
  // Generic assembly source term brick
  //
  // ----------------------------------------------------------------------

  struct gen_source_term_assembly_brick : public virtual_brick {

    std::string expr, directvarname, directdataname;
    model::varnamelist vl_test1;
    std::string secondary_domain;

    void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                const model::varnamelist &,
                                const model::varnamelist &,
                                const model::mimlist &mims,
                                model::real_matlist &,
                                model::real_veclist &vecl,
                                model::real_veclist &,
                                size_type region,
                                build_version) const override {
      GMM_ASSERT1(vecl.size() ==  vl_test1.size()
                  + ((directdataname.size() == 0) ? 0 : 1), "Wrong number "
                  "of terms for Generic source term assembly brick ");
      GMM_ASSERT1(mims.size() == 1, "Generic source term assembly brick "
                  "needs one and only one mesh_im");
      GMM_TRACE2("Generic source term assembly");

      gmm::clear(vecl[0]);

      if (expr.size()) {
        const mesh_im &mim = *mims[0];
        mesh_region rg(region);
        mim.linked_mesh().intersect_with_mpi_region(rg);

        // reenables disabled variables
        ga_workspace workspace(md, ga_workspace::inherit::ALL);
        GMM_TRACE2(name << ": generic source term assembly");
        workspace.add_expression(expr, mim, rg, 1, secondary_domain);
        workspace.assembly(1);
        const auto &V=workspace.assembled_vector();
        for (size_type i = 0; i < vl_test1.size(); ++i) {
          const auto &I=workspace.interval_of_variable(vl_test1[i]);
          gmm::copy(gmm::sub_vector(V, I), vecl[i]);
        }
      }

      if (directvarname.size()) {
        gmm::copy(md.real_variable(directdataname), vecl.back());
      }
    }

    void real_post_assembly_in_serial(const model &md, size_type ib,
                                      const model::varnamelist &/* vl */,
                                      const model::varnamelist &/* dl */,
                                      const model::mimlist &/* mims */,
                                      model::real_matlist &/*matl*/,
                                      model::real_veclist &vecl,
                                      model::real_veclist &,
                                      size_type /*region*/,
                                      build_version) const override {
      scalar_type el = scalar_type(0);
      for (size_type i=0; i < vecl.size(); ++i) el += gmm::vect_norm1(vecl[i]);
      md.add_external_load(ib, el);
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    gen_source_term_assembly_brick(const std::string &expr_,
                                   std::string brickname,
                                   const model::varnamelist &vl_test1_,
                                   const std::string &directvarname_,
                                   const std::string &directdataname_,
                                   const std::string &secdom)
      : vl_test1(vl_test1_), secondary_domain(secdom) {
      if (brickname.size() == 0)
        brickname = "Generic source term assembly brick";
      expr = expr_;
      set_flags(brickname, true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, false /* is complex */,
                false /* compute each time */);
      directvarname = directvarname_; directdataname = directdataname_;
    }

  };

  static bool check_compatibility_vl_test(model &md,
                                          const model::varnamelist vl_test) {
    model::varnamelist org;
    for (size_type i = 0; i < vl_test.size(); ++i) {
      if (md.is_affine_dependent_variable(vl_test[i]))
        org.push_back(md.org_variable(vl_test[i]));
    }
    for (size_type i = 0; i < vl_test.size(); ++i)
      for (size_type j = 0; j < org.size(); ++j)
        if (vl_test[i].compare(org[j]) == 0) return false;
    return true;
  }

  size_type add_source_term_
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   const std::string &brickname, std::string directvarname,
   const std::string &directdataname, bool return_if_nonlin,
   const std::string &secondary_domain) {

    ga_workspace workspace(md);
    size_type order = workspace.add_expression(expr, mim, region, 1,
                                               secondary_domain);
    GMM_ASSERT1(order <= 1, "Wrong order for a source term");
    model::varnamelist vl, vl_test1, vl_test2, dl;
    bool is_lin = workspace.used_variables(vl, vl_test1, vl_test2, dl, 1);
    if (!is_lin && return_if_nonlin) return size_type(-1);
    GMM_ASSERT1(is_lin, "Nonlinear term");
    GMM_ASSERT1(check_compatibility_vl_test(md, vl_test1),
                "This brick do not support the assembly on both an affine "
                "dependent variable and its original variable. "
                "Split the brick.");

    if (directdataname.size()) {
      vl.push_back(directvarname);
      dl.push_back(directdataname);
    } else directvarname = "";

    pbrick pbr = std::make_shared<gen_source_term_assembly_brick>
      (expr, brickname, vl_test1, directvarname, directdataname,
       secondary_domain);
    model::termlist tl;

    for (size_type i = 0; i < vl_test1.size(); ++i)
      tl.push_back(model::term_description(vl_test1[i]));
    if (directdataname.size())
      tl.push_back(model::term_description(directvarname));

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  size_type add_source_term
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   const std::string &brickname, const std::string &directvarname,
   const std::string &directdataname, bool return_if_nonlin) {
    return add_source_term_(md, mim, expr, region, brickname, directvarname,
                            directdataname, return_if_nonlin, "");
  }

  size_type add_twodomain_source_term
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   const std::string &secondary_domain,
   const std::string &brickname, const std::string &directvarname,
   const std::string &directdataname, bool return_if_nonlin) {
    return add_source_term_(md, mim, expr, region, brickname, directvarname,
                            directdataname, return_if_nonlin, secondary_domain);
  }

  // ----------------------------------------------------------------------
  //
  // Linear generic assembly brick
  //
  // ----------------------------------------------------------------------

  struct gen_linear_assembly_brick : public virtual_brick {

    std::string expr;
    bool is_lower_dim;
    model::varnamelist vl_test1, vl_test2;
    std::string secondary_domain;

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &/* vl */,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &/* vecl */,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(matl.size() == vl_test1.size(),
                  "Wrong number of terms for Generic linear assembly brick");
      GMM_ASSERT1(mims.size() == 1,
                  "Generic linear assembly brick needs one and only one "
                  "mesh_im");
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0);
      for (size_type i = 0; i < dl.size(); ++i)
        recompute_matrix = recompute_matrix ||
          md.is_var_newer_than_brick(dl[i], ib);

      if (recompute_matrix) {
        // reenables disabled variables
        ga_workspace workspace(md, ga_workspace::inherit::ALL);
        workspace.add_expression(expr, *(mims[0]), region, 2, secondary_domain);
        GMM_TRACE2(name << ": generic matrix assembly");
        workspace.assembly(2);
        const auto &R=workspace.assembled_matrix();
        for (size_type i = 0; i < vl_test1.size(); ++i) {
          scalar_type alpha = scalar_type(1)
            / ( workspace.factor_of_variable(vl_test1[i]) *
                workspace.factor_of_variable(vl_test2[i]));
          const auto &I1=workspace.interval_of_variable(vl_test1[i]),
                     &I2=workspace.interval_of_variable(vl_test2[i]);
          gmm::copy(gmm::scaled(gmm::sub_matrix(R, I1, I2), alpha),
                    matl[i]);
        }
      }
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return is_lower_dim ? std::string() : expr;
    }

    gen_linear_assembly_brick(const std::string &expr_, const mesh_im &mim,
                              bool is_sym,
                              bool is_coer, std::string brickname,
                              const model::varnamelist &vl_test1_,
                              const model::varnamelist &vl_test2_,
                              const std::string &secdom)
      : vl_test1(vl_test1_), vl_test2(vl_test2_), secondary_domain(secdom) {
      if (brickname.size() == 0) brickname = "Generic linear assembly brick";
      expr = expr_;
      is_lower_dim = mim.is_lower_dimensional();
      set_flags(brickname, true /* is linear*/,
                is_sym /* is symmetric */, is_coer /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };

  static bool check_compatibility_vl_test(model &md,
                                          const model::varnamelist vl_test1,
                                          const model::varnamelist vl_test2) {
    for (size_type i = 0; i < vl_test1.size(); ++i)
      for (size_type j = 0; j < vl_test1.size(); ++j) {
        bool is1 = md.is_affine_dependent_variable(vl_test1[i]);
        bool is2 = md.is_affine_dependent_variable(vl_test2[i]);
        if (is1 || is2) {
          const std::string &org1
            = is1 ? md.org_variable(vl_test1[i]) : vl_test1[i];
          const std::string &org2
            = is2 ? md.org_variable(vl_test2[i]) : vl_test2[i];
          bool is1_bis = md.is_affine_dependent_variable(vl_test1[j]);
          bool is2_bis = md.is_affine_dependent_variable(vl_test2[j]);
          const std::string &org1_bis = is1_bis ? md.org_variable(vl_test1[j])
            : vl_test1[j];
          const std::string &org2_bis = is2_bis ? md.org_variable(vl_test2[j])
            : vl_test2[j];
          if (org1.compare(org1_bis) == 0 && org2.compare(org2_bis))
            return false;
        }
      }
    return true;
  }



  size_type add_linear_term_
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   bool is_sym, bool is_coercive, const std::string &brickname,
   bool return_if_nonlin, const std::string &secondary_domain) {
    // reenables disabled variables
    ga_workspace workspace(md, ga_workspace::inherit::ALL);
    size_type order = workspace.add_expression(expr, mim, region,
                                               2, secondary_domain);
    model::varnamelist vl, vl_test1, vl_test2, dl;
    bool is_lin = workspace.used_variables(vl, vl_test1, vl_test2, dl, 2);

    if (!is_lin && return_if_nonlin) return size_type(-1);
    GMM_ASSERT1(is_lin, "Nonlinear term");
    if (order == 0) { is_coercive = is_sym = true; }

    size_type brick_id(-1);
    std::string const_expr= workspace.extract_constant_term(mim.linked_mesh());
    if (const_expr.size())
      brick_id = add_source_term_(md, mim, const_expr, region,
                                  brickname+" (source term)",
                                  "", "", false, secondary_domain);

    // GMM_ASSERT1(order <= 1,
    //             "This brick does not support a second order term");
    GMM_ASSERT1(check_compatibility_vl_test(md, vl_test1, vl_test2),
                "This brick do not support the assembly on both an affine "
                "dependent variable and its original variable. "
                "Split the brick.");

    if (vl_test1.size()) {
      pbrick pbr = std::make_shared<gen_linear_assembly_brick>
        (expr, mim, is_sym, is_coercive, brickname, vl_test1, vl_test2,
         secondary_domain);
      model::termlist tl;
      for (size_type i = 0; i < vl_test1.size(); ++i)
        tl.push_back(model::term_description(vl_test1[i], vl_test2[i], false));

      return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
    }
    return brick_id;
  }

  size_type add_linear_term
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   bool is_sym, bool is_coercive, const std::string &brickname,
   bool return_if_nonlin) {
    return add_linear_term_(md, mim, expr, region, is_sym, is_coercive,
                            brickname, return_if_nonlin, "");
  }

  size_type add_linear_twodomain_term
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   const std::string &secondary_domain, bool is_sym, bool is_coercive,
   const std::string &brickname, bool return_if_nonlin) {
    return add_linear_term_(md, mim, expr, region, is_sym, is_coercive,
                            brickname, return_if_nonlin, secondary_domain);
  }


  // ----------------------------------------------------------------------
  //
  // Nonlinear generic assembly brick
  //
  // ----------------------------------------------------------------------

  struct gen_nonlinear_assembly_brick : public virtual_brick {

    std::string expr;
    bool is_lower_dim;
    std::string secondary_domain;

    virtual void real_post_assembly_in_serial(const model &md, size_type ,
                                              const model::varnamelist &,
                                              const model::varnamelist &,
                                              const model::mimlist &mims,
                                              model::real_matlist &,
                                              model::real_veclist &,
                                              model::real_veclist &,
                                              size_type region,
                                              build_version) const {
      GMM_ASSERT1(mims.size() == 1,
                  "Generic linear assembly brick needs one and only one "
                  "mesh_im");
      md.add_generic_expression(expr, *(mims[0]), region, secondary_domain);
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return expr;
    }


    gen_nonlinear_assembly_brick(const std::string &expr_, const mesh_im &mim,
                                 bool is_sym,
                                 bool is_coer,
                                 std::string brickname,
                                 const std::string &secdom) {
      if (brickname.size() == 0) brickname = "Generic linear assembly brick";
      expr = expr_;
      secondary_domain = secdom;
      is_lower_dim = mim.is_lower_dimensional();
      set_flags(brickname, false /* is linear*/,
                is_sym /* is symmetric */, is_coer /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };

  size_type add_nonlinear_term_
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   bool is_sym, bool is_coercive, const std::string &brickname,
   const std::string &secondary_domain) {

    ga_workspace workspace(md);
    size_type order = workspace.add_expression(expr, mim, region, 2,
                                               secondary_domain);
    GMM_ASSERT1(order < 2, "Order two test functions (Test2) are not allowed"
                " in assembly string for nonlinear terms");
    model::varnamelist vl, vl_test1, vl_test2, ddl, dl;
    workspace.used_variables(vl, vl_test1, vl_test2, ddl, order);
    for (size_type i = 0; i < ddl.size(); ++i)
      if (md.is_true_data(ddl[i])) dl.push_back(ddl[i]);
      else vl.push_back(ddl[i]);
    if (order == 0) { is_coercive = is_sym = true; }
    pbrick pbr = std::make_shared<gen_nonlinear_assembly_brick>
      (expr, mim, is_sym, is_coercive, brickname, secondary_domain);
    model::termlist tl; // No term
    // tl.push_back(model::term_description(true, is_sym));
    // TODO to be changed.
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  size_type add_nonlinear_term
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   bool is_sym, bool is_coercive, const std::string &brickname) {
    return add_nonlinear_term_(md, mim, expr, region, is_sym, is_coercive,
                               brickname, "");
  }

  size_type add_nonlinear_twodomain_term
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   const std::string &secondary_domain, bool is_sym, bool is_coercive,
   const std::string &brickname) {
    return add_nonlinear_term_(md, mim, expr, region, is_sym, is_coercive,
                               brickname, secondary_domain);
  }


  // ----------------------------------------------------------------------
  //
  // Generic elliptic brick
  //
  // ----------------------------------------------------------------------

  // Kept only for the complex version
  struct generic_elliptic_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type /*ib*/,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Generic elliptic brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Generic elliptic brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() <= 1,
                  "Wrong number of variables for generic elliptic brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      size_type N = m.dim(), Q = mf_u.get_qdim(), s = 1;
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector *A = 0;
      const mesh_fem *mf_a = 0;
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);
      if (dl.size() > 0) {
        A = &(md.real_variable(dl[0]));
        mf_a = md.pmesh_fem_of_variable(dl[0]);
        s = gmm::vect_size(*A);
        if (mf_a) s = s * mf_a->get_qdim() / mf_a->nb_dof();
      }

      gmm::clear(matl[0]);
      GMM_TRACE2("Generic elliptic term assembly");
      if (s == 1) {
        if (mf_a) {
          if (Q > 1)
            asm_stiffness_matrix_for_laplacian_componentwise
              (matl[0], mim, mf_u, *mf_a, *A, rg);
          else
            asm_stiffness_matrix_for_laplacian
              (matl[0], mim, mf_u, *mf_a, *A, rg);

        } else {
          if (Q > 1)
            asm_stiffness_matrix_for_homogeneous_laplacian_componentwise
              (matl[0], mim, mf_u, rg);
          else
            asm_stiffness_matrix_for_homogeneous_laplacian
              (matl[0], mim, mf_u, rg);
          if (A) gmm::scale(matl[0], (*A)[0]);
        }
      } else if (s == N*N) {
        if (mf_a) {
          if (Q > 1)
            asm_stiffness_matrix_for_scalar_elliptic_componentwise
              (matl[0], mim, mf_u, *mf_a, *A, rg);
          else
            asm_stiffness_matrix_for_scalar_elliptic
              (matl[0], mim, mf_u, *mf_a, *A, rg);
        } else {
          if (Q > 1)
            asm_stiffness_matrix_for_homogeneous_scalar_elliptic_componentwise
              (matl[0], mim, mf_u, *A, rg);
          else
            asm_stiffness_matrix_for_homogeneous_scalar_elliptic
              (matl[0], mim, mf_u, *A, rg);
        }
      } else if (s == N*N*Q*Q) {
        if (mf_a)
          asm_stiffness_matrix_for_vector_elliptic
            (matl[0], mim, mf_u, *mf_a, *A, rg);
        else
          asm_stiffness_matrix_for_homogeneous_vector_elliptic
            (matl[0], mim, mf_u, *A, rg);
      } else
        GMM_ASSERT1(false, "Bad format generic elliptic brick coefficient");

    }

    virtual void real_post_assembly_in_serial(const model &md, size_type,
                                              const model::varnamelist &,
                                              const model::varnamelist &dl,
                                              const model::mimlist &/* mims */,
                                              model::real_matlist &/*matl*/,
                                              model::real_veclist &,
                                              model::real_veclist &,
                                              size_type /*region*/,
                                              build_version) const
    {
      const model_real_plain_vector *A = 0;
      const mesh_fem *mf_a = 0;
      if (dl.size() > 0)
      {
        A = &(md.real_variable(dl[0]));
        mf_a = md.pmesh_fem_of_variable(dl[0]);
      }
    }

    virtual void complex_post_assembly_in_serial(const model &md, size_type,
                                              const model::varnamelist &,
                                              const model::varnamelist &dl,
                                              const model::mimlist &/*mims*/,
                                              model::complex_matlist &/*matl*/,
                                              model::complex_veclist &,
                                              model::complex_veclist &,
                                              size_type /* region */,
                                              build_version) const
    {
      const model_real_plain_vector *A = 0;
      const mesh_fem *mf_a = 0;
      if (dl.size() > 0)
      {
        A = &(md.real_variable(dl[0]));
        mf_a = md.pmesh_fem_of_variable(dl[0]);
      }
    }

    virtual scalar_type asm_complex_tangent_terms(const model &md, size_type,
                                                  const model::varnamelist &vl,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::complex_matlist &matl,
                                                  model::complex_veclist &,
                                                  model::complex_veclist &,
                                                  size_type) const {
      const model_complex_plain_vector &U = md.complex_variable(vl[0]);
      return gmm::abs(gmm::vect_hp(matl[0], U, U)) / scalar_type(2);
    }


    virtual void asm_complex_tangent_terms(const model &md, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Generic elliptic brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Generic elliptic brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() <= 1,
                  "Wrong number of variables for generic elliptic brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      size_type N = m.dim(), Q = mf_u.get_qdim(), s = 1;
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector *A = 0;
      const mesh_fem *mf_a = 0;
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);


      if (dl.size() > 0) {
        A = &(md.real_variable(dl[0]));
        mf_a = md.pmesh_fem_of_variable(dl[0]);
        s = gmm::vect_size(*A);
        if (mf_a) s = s * mf_a->get_qdim() / mf_a->nb_dof();
      }

      gmm::clear(matl[0]);
      GMM_TRACE2("Generic elliptic term assembly");
      if (s == 1) {
        if (mf_a) {
          if (Q > 1)
            asm_stiffness_matrix_for_laplacian_componentwise
              (matl[0], mim, mf_u, *mf_a, *A, rg);
          else
            asm_stiffness_matrix_for_laplacian
              (matl[0], mim, mf_u, *mf_a, *A, rg);

        } else {
          if (Q > 1)
            asm_stiffness_matrix_for_homogeneous_laplacian_componentwise
              (gmm::real_part(matl[0]), mim, mf_u, rg);
          else
            asm_stiffness_matrix_for_homogeneous_laplacian
              (gmm::real_part(matl[0]), mim, mf_u, rg);
          if (A) gmm::scale(matl[0], (*A)[0]);
        }
      } else if (s == N*N) {
        if (mf_a) {
          if (Q > 1)
            asm_stiffness_matrix_for_scalar_elliptic_componentwise
              (matl[0], mim, mf_u, *mf_a, *A, rg);
          else
            asm_stiffness_matrix_for_scalar_elliptic
              (matl[0], mim, mf_u, *mf_a, *A, rg);
        } else {
          if (Q > 1)
            asm_stiffness_matrix_for_homogeneous_scalar_elliptic_componentwise
              (matl[0], mim, mf_u, *A, rg);
          else
            asm_stiffness_matrix_for_homogeneous_scalar_elliptic
              (matl[0], mim, mf_u, *A, rg);
        }
      } else if (s == N*N*Q*Q) {
        if (mf_a)
          asm_stiffness_matrix_for_vector_elliptic
            (matl[0], mim, mf_u, *mf_a, *A, rg);
        else
          asm_stiffness_matrix_for_homogeneous_vector_elliptic
            (matl[0], mim, mf_u, *A, rg);
      } else
        GMM_ASSERT1(false,
                    "Bad format generic elliptic brick coefficient");
    }

    generic_elliptic_brick() {
      set_flags("Generic elliptic", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */);
    }

  };

  size_type add_Laplacian_brick(model &md, const mesh_im &mim,
                                const std::string &varname,
                                size_type region) {
    if (md.is_complex()) {
      pbrick pbr = std::make_shared<generic_elliptic_brick>();
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          model::varnamelist(), tl, model::mimlist(1, &mim),
                          region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
      size_type qdim = mf_u.get_qdim();
      std::string expr;
      if (qdim == 1)
        expr = "Grad_"+varname+".Grad_"+test_varname;
      else
        expr = "Grad_"+varname+":Grad_"+test_varname;
      return add_linear_term(md, mim, expr, region, true, true,
                             "Laplacian", false);
    }
  }

  size_type add_generic_elliptic_brick(model &md, const mesh_im &mim,
                                       const std::string &varname,
                                       const std::string &dataname,
                                       size_type region) {
    if (md.is_complex()) {
      pbrick pbr = std::make_shared<generic_elliptic_brick>();
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          model::varnamelist(1, dataname), tl,
                          model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
      size_type dim = mf_u.linked_mesh().dim(), qdim = mf_u.get_qdim(), qdim_data = 1;
      std::string expr;

      if (md.variable_exists(dataname)) {
        const mesh_fem *mf = md.pmesh_fem_of_variable(dataname);
        size_type n = gmm::vect_size(md.real_variable(dataname));
        if (mf) qdim_data = mf->get_qdim() * (n / mf->nb_dof());
        else  qdim_data = n;
      }

      if (qdim == 1) {
        if (qdim_data != 1) {
          GMM_ASSERT1(qdim_data == gmm::sqr(dim),
                      "Wrong data size for generic elliptic brick");
          expr = "((Reshape("+dataname+",meshdim,meshdim))*Grad_"+varname+").Grad_"
            + test_varname;
        } else {
          expr = "(("+dataname+")*Grad_"+varname+").Grad_"+test_varname;
        }
      } else {
        if (qdim_data != 1) {
          if (qdim_data == gmm::sqr(dim))
            expr = "((Reshape("+dataname+",meshdim,meshdim))*Grad_"+varname+"):Grad_"
              +test_varname;
          else if (qdim_data == gmm::sqr(gmm::sqr(dim)))
            expr = "((Reshape("+dataname+",meshdim,meshdim,meshdim,meshdim))*Grad_"
              +varname+"):Grad_"+test_varname;
          else GMM_ASSERT1(false, "Wrong data size for generic elliptic brick");
        } else {
          expr = "(("+dataname+")*Grad_"+varname+"):Grad_"+test_varname;
        }
      }
      size_type ib = add_linear_term
        (md, mim, expr, region, true, true, "Generic elliptic", true);
      if (ib == size_type(-1))
        ib = add_nonlinear_term(md, mim, expr, region, false, false,
                                "Generic elliptic (nonlinear)");
      return ib;
    }
  }

  // ----------------------------------------------------------------------
  //
  // Source term brick
  //
  // ----------------------------------------------------------------------

  // Kept only for the complex version
  struct source_term_brick : public virtual_brick {

    void asm_real_tangent_terms(const model &md, size_type /*ib*/,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const override {
      GMM_ASSERT1(vecl.size() == 1,
                  "Source term brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Source term brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() > 0 && dl.size() <= 2,
                  "Wrong number of variables for source term brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector &A = md.real_variable(dl[0]);
      const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      size_type s = gmm::vect_size(A);
      if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();

      GMM_ASSERT1(mf_u.get_qdim() == s,
                  dl[0] << ": bad format of source term data. "
                  "Detected dimension is " << s << " should be "
                  << size_type(mf_u.get_qdim()));

      GMM_TRACE2("Source term assembly");
      if (mf_data)
        asm_source_term(vecl[0], mim, mf_u, *mf_data, A, rg);
      else
        asm_homogeneous_source_term(vecl[0], mim, mf_u, A, rg);

      if (dl.size() > 1) gmm::add(md.real_variable(dl[1]), vecl[0]);
    }

    void real_post_assembly_in_serial(const model &md, size_type ib,
                                      const model::varnamelist &/* vl */,
                                      const model::varnamelist &/* dl */,
                                      const model::mimlist &/* mims */,
                                      model::real_matlist &/*matl*/,
                                      model::real_veclist &vecl,
                                      model::real_veclist &,
                                      size_type /*region*/,
                                      build_version) const override
    {
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }


    void asm_complex_tangent_terms(const model &md, size_type /*ib*/,
                                   const model::varnamelist &vl,
                                   const model::varnamelist &dl,
                                   const model::mimlist &mims,
                                   model::complex_matlist &,
                                   model::complex_veclist &vecl,
                                   model::complex_veclist &,
                                   size_type region,
                                   build_version) const override {
      GMM_ASSERT1(vecl.size() == 1,
                  "Source term brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Source term brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() > 0 && dl.size() <= 2,
                  "Wrong number of variables for source term brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_im &mim = *mims[0];
      const model_complex_plain_vector &A = md.complex_variable(dl[0]);
      const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      size_type s = gmm::vect_size(A);
      if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();

      GMM_ASSERT1(mf_u.get_qdim() == s, "Bad format of source term data");

      GMM_TRACE2("Source term assembly");
      if (mf_data)
        asm_source_term(vecl[0], mim, mf_u, *mf_data, A, rg);
      else
        asm_homogeneous_source_term(vecl[0], mim, mf_u, A, rg);

      if (dl.size() > 1) gmm::add(md.complex_variable(dl[1]), vecl[0]);
    }

    void complex_post_assembly_in_serial(const model &md,
                                         size_type ib,
                                         const model::varnamelist &,
                                         const model::varnamelist &,
                                         const model::mimlist &,
                                         model::complex_matlist &,
                                         model::complex_veclist &vecl,
                                         model::complex_veclist &,
                                         size_type, build_version) const override
    {
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }



    source_term_brick() {
      set_flags("Source term", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }


  };

  size_type add_source_term_brick(model &md, const mesh_im &mim,
                                  const std::string &varname,
                                  const std::string &dataexpr,
                                  size_type region,
                                  const std::string &directdataname) {
    if (md.is_complex()) {
      pbrick pbr = std::make_shared<source_term_brick>();
      model::termlist tl;
      tl.push_back(model::term_description(varname));
      model::varnamelist vdata(1, dataexpr);
      if (directdataname.size()) vdata.push_back(directdataname);
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          vdata, tl, model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
      size_type qdim = mf_u.get_qdim();
      std::string expr;
      if (qdim == 1)
        expr = "("+dataexpr+")*"+test_varname;
      else
        expr = "("+dataexpr+")."+test_varname;
      size_type ib = add_source_term_generic_assembly_brick
        (md, mim, expr, region, "Source term", varname, directdataname, true);
      if (ib == size_type(-1)) {
        ib = add_nonlinear_term(md, mim, "-("+expr+")", region, false, false,
                                "Source term (nonlinear)");
        if (directdataname.size())
          add_source_term_generic_assembly_brick
            (md, mim, "", region, "Source term", varname, directdataname);
      }
      return ib;
    }
  }

  // ----------------------------------------------------------------------
  //
  // Normal source term brick
  //
  // ----------------------------------------------------------------------

  struct normal_source_term_brick : public virtual_brick {

    void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                const model::varnamelist &vl,
                                const model::varnamelist &dl,
                                const model::mimlist &mims,
                                model::real_matlist &,
                                model::real_veclist &vecl,
                                model::real_veclist &,
                                size_type region,
                                build_version) const override {
      GMM_ASSERT1(vecl.size() == 1,
                  "Source term brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Source term brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
                  "Wrong number of variables for source term brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector &A = md.real_variable(dl[0]);
      const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      size_type s = gmm::vect_size(A), N = mf_u.linked_mesh().dim();
      if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();

      GMM_ASSERT1(mf_u.get_qdim()*N == s,
                  dl[0] << ": bad format of normal source term data. "
                  "Detected dimension is " << s << " should be "
                  << size_type(mf_u.get_qdim()*N));

      GMM_TRACE2("source term assembly");
      if (mf_data)
        asm_normal_source_term(vecl[0], mim, mf_u, *mf_data, A, rg);
      else
        asm_homogeneous_normal_source_term(vecl[0], mim, mf_u, A, rg);
    }

    void real_post_assembly_in_serial(const model &md, size_type ib,
                                      const model::varnamelist &/* vl */,
                                      const model::varnamelist &/* dl */,
                                      const model::mimlist &/* mims */,
                                      model::real_matlist &/*matl*/,
                                      model::real_veclist &vecl,
                                      model::real_veclist &,
                                      size_type /*region*/,
                                      build_version) const override {
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }


    virtual void asm_complex_tangent_terms(const model &md, size_type /* ib */,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version) const {
      GMM_ASSERT1(vecl.size() == 1,
                  "Source term brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Source term brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
                  "Wrong number of variables for source term brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_im &mim = *mims[0];
      const model_complex_plain_vector &A = md.complex_variable(dl[0]);
      const mesh_fem *mf_data = md.pmesh_fem_of_variable(dl[0]);
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      size_type s = gmm::vect_size(A), N = mf_u.linked_mesh().dim();
      if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();

      GMM_ASSERT1(s == mf_u.get_qdim()*N, "Bad format of source term data");

      GMM_TRACE2("source term assembly");
      if (mf_data)
        asm_normal_source_term(vecl[0], mim, mf_u, *mf_data, A, rg);
      else
        asm_homogeneous_normal_source_term(vecl[0], mim, mf_u, A, rg);
    }

    void complex_post_assembly_in_serial(const model &md, size_type ib,
                                         const model::varnamelist &/* vl */,
                                         const model::varnamelist &/* dl */,
                                         const model::mimlist &/* mims */,
                                         model::complex_matlist &/*matl*/,
                                         model::complex_veclist &vecl,
                                         model::complex_veclist &,
                                         size_type /*region*/,
                                         build_version) const override {
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }

    normal_source_term_brick() {
      set_flags("Normal source term", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }


  };

  size_type add_normal_source_term_brick(model &md, const mesh_im &mim,
                                         const std::string &varname,
                                         const std::string &dataexpr,
                                         size_type region) {
    if (md.is_complex()) {
      pbrick pbr = std::make_shared<normal_source_term_brick>();
      model::termlist tl;
      tl.push_back(model::term_description(varname));
      model::varnamelist vdata(1, dataexpr);
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          vdata, tl, model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
      size_type qdim = mf_u.get_qdim();
      std::string expr;
      if (qdim == 1)
        expr = "(("+dataexpr+").Normal)*"+test_varname;
      else
        expr = "(Reshape("+dataexpr+",qdim("+varname
          + "),meshdim)*Normal)."+test_varname;
      return add_source_term_generic_assembly_brick
        (md, mim, expr, region, "Source term");
    }
  }


  // ----------------------------------------------------------------------
  //
  // Dirichlet condition brick
  //
  // ----------------------------------------------------------------------
  // Two variables : with multipliers
  // One variable : penalization

  struct Dirichlet_condition_brick : public virtual_brick {

    bool H_version; // The version hu = r for vector fields.
    bool normal_component; // Dirichlet on normal component for vector field.
    const mesh_fem *mf_mult_;
    mutable getfem::omp_distribute<model_real_sparse_matrix> rB_th;
    mutable getfem::omp_distribute<model_real_plain_vector> rV_th;
    mutable getfem::omp_distribute<model_complex_sparse_matrix> cB_th;
    mutable getfem::omp_distribute<model_complex_plain_vector> cV_th;

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Dirichlet condition brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Dirichlet condition brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 3,
                  "Wrong number of variables for Dirichlet condition brick");

      model_real_sparse_matrix& rB = rB_th;
      model_real_plain_vector&  rV = rV_th;

      bool penalized = (vl.size() == 1);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_mult = penalized ? (mf_mult_ ? *mf_mult_ : mf_u)
        : md.mesh_fem_of_variable(vl[1]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector *A = 0, *COEFF = 0, *H = 0;
      const mesh_fem *mf_data = 0, *mf_H = 0;
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
        || (penalized && md.is_var_newer_than_brick(dl[0], ib));

      if (penalized) {
        COEFF = &(md.real_variable(dl[0]));
        GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                    "Data for coefficient should be a scalar");
      }

      size_type s = 0, ind = (penalized ? 1 : 0);
      if (dl.size() > ind) {
        A = &(md.real_variable(dl[ind]));
        mf_data = md.pmesh_fem_of_variable(dl[ind]);
        s = gmm::vect_size(*A);
        if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();
        size_type ss = ((normal_component) ? 1 :  mf_u.get_qdim());
        GMM_ASSERT1(s == ss, dl[ind] << ": bad format of "
                    "Dirichlet data. Detected dimension is " << s
                    << " should be " << ss);
      }

      if (dl.size() > ind + 1) {
        GMM_ASSERT1(H_version,
                    "Wrong number of data for Dirichlet condition brick");
        H = &(md.real_variable(dl[ind+1]));
        mf_H = md.pmesh_fem_of_variable(dl[ind+1]);
        s = gmm::vect_size(*A);
        if (mf_H) {
          s = s * mf_H->get_qdim() / mf_H->nb_dof();
          GMM_ASSERT1(mf_H->get_qdim() == 1,  "Implemented only for mf_H "
                      "a scalar finite element method");
        }
        GMM_ASSERT1(s = gmm::sqr(mf_u.get_qdim()),
                    dl[ind+1] << ": bad format of Dirichlet data. "
                    "Detected dimension is " << s << " should be "
                    << size_type(gmm::sqr(mf_u.get_qdim())));
      }

      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (recompute_matrix) {
        model_real_sparse_matrix *B = &(matl[0]);
        if (penalized && (&mf_mult != &mf_u)) {
          gmm::resize(rB, mf_mult.nb_dof(), mf_u.nb_dof());
          gmm::clear(rB);
          B = &rB;
        } else {
          gmm::clear(matl[0]);
        }
        GMM_TRACE2("Mass term assembly for Dirichlet condition");
        if (H_version || normal_component) {
          ga_workspace workspace;
          gmm::sub_interval Imult(0, mf_mult.nb_dof()), Iu(0, mf_u.nb_dof());
          base_vector u(mf_u.nb_dof());
          base_vector mult(mf_mult.nb_dof());
          workspace.add_fem_variable("u", mf_u, Iu, u);
          workspace.add_fem_variable("mult", mf_mult, Imult, mult);
          auto expression = std::string{};
          if (H_version){
            if (mf_H) workspace.add_fem_constant("A", *mf_H, *H);
            else workspace.add_fixed_size_constant("A", *H);
            expression = (mf_u.get_qdim() == 1) ?
                          "Test_mult . (A . Test2_u)"
                          :
                          "Test_mult. (Reshape(A, qdim(u), qdim(u)) . Test2_u)";
          } else if (normal_component) {
            expression = "Test_mult . (Test2_u . Normal)";
          }
          workspace.add_expression(expression, mim, rg);
          workspace.set_assembled_matrix(*B);
          workspace.assembly(2);
        } else {
          asm_mass_matrix(*B, mim, mf_mult, mf_u, rg);
        }

        if (penalized && (&mf_mult != &mf_u)) {
          GMM_ASSERT1(MPI_IS_MASTER(), "Sorry, the penalized option "
                      "filtered by a multiplier space is not parallelized");
          gmm::mult(gmm::transposed(rB), rB, matl[0]);
          gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
        } else if (penalized) {
          gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
        }
      }

      if (dl.size() > ind) {
        GMM_TRACE2("Source term assembly for Dirichlet condition");

        if (penalized && (&mf_mult != &mf_u)) {
          gmm::resize(rV, mf_mult.nb_dof());
          gmm::clear(rV);
          if (mf_data)
            asm_source_term(rV, mim, mf_mult, *mf_data, *A, rg);
          else
            asm_homogeneous_source_term(rV, mim, mf_mult, *A, rg);
        } else {
          if (mf_data)
            asm_source_term(vecl[0], mim, mf_mult, *mf_data, *A, rg);
          else
            asm_homogeneous_source_term(vecl[0], mim, mf_mult, *A, rg);
        }

        if (penalized && (&mf_mult != &mf_u))  {
          gmm::mult(gmm::transposed(rB), rV, vecl[0]);
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
          rV = model_real_plain_vector();
        } else if (penalized)
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }

    }

    void real_post_assembly_in_serial(const model &md, size_type ib,
                                      const model::varnamelist &/* vl */,
                                      const model::varnamelist &/* dl */,
                                      const model::mimlist &/* mims */,
                                      model::real_matlist &/*matl*/,
                                      model::real_veclist &vecl,
                                      model::real_veclist &,
                                      size_type /*region*/,
                                      build_version) const override {
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type ib,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Dirichlet condition brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Dirichlet condition brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 3,
                  "Wrong number of variables for Dirichlet condition brick");

      model_complex_sparse_matrix& cB = cB_th;
      model_complex_plain_vector&  cV = cV_th;

      bool penalized = (vl.size() == 1);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_mult = penalized ? (mf_mult_ ? *mf_mult_ : mf_u)
        : md.mesh_fem_of_variable(vl[1]);
      const mesh_im &mim = *mims[0];
      const model_complex_plain_vector *A = 0, *COEFF = 0, *H = 0;
      const mesh_fem *mf_data = 0, *mf_H = 0;
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
        || (penalized && md.is_var_newer_than_brick(dl[0], ib));

      if (penalized) {
        COEFF = &(md.complex_variable(dl[0]));
        GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                    "Data for coefficient should be a scalar");
      }

      size_type s = 0, ind = (penalized ? 1 : 0);
      if (dl.size() > ind) {
        A = &(md.complex_variable(dl[ind]));
        mf_data = md.pmesh_fem_of_variable(dl[ind]);
        s = gmm::vect_size(*A);
        if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();
        size_type ss = s * ((normal_component) ? mf_u.linked_mesh().dim() : 1);
        GMM_ASSERT1(mf_u.get_qdim() == ss,
                    dl[ind] << ": bad format of Dirichlet data. "
                    "Detected dimension is " << ss << " should be "
                    << size_type(mf_u.get_qdim()));
      }

      if (dl.size() > ind + 1) {
        GMM_ASSERT1(H_version,
                    "Wrong number of data for Dirichlet condition brick");
        H = &(md.complex_variable(dl[ind+1]));
        mf_H = md.pmesh_fem_of_variable(dl[ind+1]);
        s = gmm::vect_size(*A);
        if (mf_H) {
          s = s * mf_H->get_qdim() / mf_H->nb_dof();
          GMM_ASSERT1(mf_H->get_qdim() == 1,  "Implemented only for mf_H "
                      "a scalar finite element method");
        }
        GMM_ASSERT1(s = gmm::sqr(mf_u.get_qdim()),
                    dl[ind+1] << ": bad format of Dirichlet data. "
                    "Detected dimension is " << s << " should be "
                    << size_type(gmm::sqr(mf_u.get_qdim())));
      }

      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (recompute_matrix) {
        model_complex_sparse_matrix *B = &(matl[0]);
        if (penalized && (&mf_mult != &mf_u)) {
          gmm::resize(cB, mf_mult.nb_dof(), mf_u.nb_dof());
          gmm::clear(cB);
          B = &cB;
        } else {
          gmm::clear(matl[0]);
        }
        GMM_TRACE2("Mass term assembly for Dirichlet condition");
        if (H_version) {
          if (mf_u.get_qdim() == 1)
            asm_real_or_complex_1_param_mat(*B, mim, mf_mult, mf_H, *H, rg,
                                            "(A*Test_u).Test2_u");
          else
            asm_real_or_complex_1_param_mat(*B, mim, mf_mult, mf_H, *H, rg,
                            "(Reshape(A,qdim(u),qdim(u))*Test2_u).Test_u");
          // if (mf_H)
          //   asm_real_or_complex_1_param
          //     (*B, mim, mf_mult, *mf_H, *H, rg, (mf_u.get_qdim() == 1) ?
          //      "F=data(#2);"
          //      "M(#1,#3)+=sym(comp(Base(#1).Base(#3).Base(#2))(:,:,i).F(i))"
          //      : "F=data(qdim(#1),qdim(#1),#2);"
          //      "M(#1,#3)+=sym(comp(vBase(#1).vBase(#3).Base(#2))(:,i,:,j,k).F(i,j,k));", &mf_u);
          // else
          //    asm_real_or_complex_1_param
          //     (*B, mim, mf_mult, mf_u, *H, rg, (mf_u.get_qdim() == 1) ?
          //      "F=data(1);"
          //      "M(#1,#2)+=sym(comp(Base(#1).Base(#2)).F(1))"
          //      : "F=data(qdim(#1),qdim(#1));"
          //      "M(#1,#2)+=sym(comp(vBase(#1).vBase(#2))(:,i,:,j).F(i,j));");
        }
        else if (normal_component) {
          ga_workspace workspace;
          gmm::sub_interval Imult(0, mf_mult.nb_dof()), Iu(0, mf_u.nb_dof());
          base_vector mult(mf_mult.nb_dof()), u(mf_u.nb_dof());
          workspace.add_fem_variable("mult", mf_mult, Imult, mult);
          workspace.add_fem_variable("u", mf_u, Iu, u);
          workspace.add_expression("Test_mult.(Test2_u.Normal)", mim, rg);
          model_real_sparse_matrix BB(mf_mult.nb_dof(), mf_u.nb_dof());
          workspace.set_assembled_matrix(BB);
          workspace.assembly(2);
          gmm::add(BB, *B);

          // generic_assembly assem;
          // if (mf_mult.get_qdim() == 1)
          //   assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1).Normal())(:,:,i,i);");
          // else
          //   assem.set("M(#2,#1)+=comp(vBase(#2).mBase(#1).Normal())(:,i,:,i,j,j);");
          // assem.push_mi(mim);
          // assem.push_mf(mf_u);
          // assem.push_mf(mf_mult);
          // assem.push_mat(gmm::real_part(*B));
          // assem.assembly(rg);
        } else {
          asm_mass_matrix(*B, mim, mf_mult, mf_u, rg);
        }
        if (penalized && (&mf_mult != &mf_u)) {
          gmm::mult(gmm::transposed(cB), cB, matl[0]);
          gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
        } else if (penalized) {
          gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
        }
      }

      if (dl.size() > ind) {
        GMM_TRACE2("Source term assembly for Dirichlet condition");

        if (penalized && (&mf_mult != &mf_u)) {
          gmm::resize(cV, mf_mult.nb_dof());
          gmm::clear(cV);
          if (mf_data)
            asm_source_term(cV, mim, mf_mult, *mf_data, *A, rg);
          else
            asm_homogeneous_source_term(cV, mim, mf_mult, *A, rg);
        } else {
          if (mf_data)
            asm_source_term(vecl[0], mim, mf_mult, *mf_data, *A, rg);
          else
            asm_homogeneous_source_term(vecl[0], mim, mf_mult, *A, rg);
        }

        if (penalized && (&mf_mult != &mf_u))  {
          gmm::mult(gmm::transposed(cB), cV, vecl[0]);
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
          cV = model_complex_plain_vector();
        } else if (penalized)
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }
    }

    void complex_post_assembly_in_serial(const model &md, size_type ib,
                                         const model::varnamelist &/* vl */,
                                         const model::varnamelist &/* dl */,
                                         const model::mimlist &/* mims */,
                                         model::complex_matlist &/*matl*/,
                                         model::complex_veclist &vecl,
                                         model::complex_veclist &,
                                         size_type /*region*/,
                                         build_version) const override {
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }


    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    Dirichlet_condition_brick(bool penalized, bool H_version_,
                              bool normal_component_,
                              const mesh_fem *mf_mult__ = 0) {
      mf_mult_ = mf_mult__;
      H_version = H_version_;
      normal_component = normal_component_;
      GMM_ASSERT1(!(H_version && normal_component), "Bad Dirichlet version");
      set_flags(penalized ? "Dirichlet with penalization brick"
                          : "Dirichlet with multipliers brick",
                true /* is linear*/,
                true /* is symmetric */, penalized /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }
  };

  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname) {
    pbrick pbr = std::make_shared<Dirichlet_condition_brick>(false,false,false);
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const mesh_fem &mf_mult, size_type region,
   const std::string &dataname) {
    std::string multname = md.new_name("mult_on_" + varname);
    md.add_multiplier(multname, mf_mult, varname);
    return add_Dirichlet_condition_with_multipliers
      (md, mim, varname, multname, region, dataname);
  }

  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   dim_type degree, size_type region,
   const std::string &dataname) {
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    const mesh_fem &mf_mult = classical_mesh_fem(mf_u.linked_mesh(),
                                                 degree, mf_u.get_qdim());
    return add_Dirichlet_condition_with_multipliers
      (md, mim, varname, mf_mult, region, dataname);
  }

  const std::string &mult_varname_Dirichlet(model &md, size_type ind_brick) {
    return md.varname_of_brick(ind_brick, 1);
  }

  size_type add_Dirichlet_condition_with_penalization
  (model &md, const mesh_im &mim, const std::string &varname,
   scalar_type penalisation_coeff, size_type region,
   const std::string &dataname, const mesh_fem *mf_mult) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = std::make_shared<Dirichlet_condition_brick>
      (true, false, false, mf_mult);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  size_type add_normal_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname) {
    pbrick pbr = std::make_shared<Dirichlet_condition_brick>(false,false,true);
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  size_type add_normal_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const mesh_fem &mf_mult, size_type region,
   const std::string &dataname) {
    std::string multname = md.new_name("mult_on_" + varname);
    md.add_multiplier(multname, mf_mult, varname);
    return add_normal_Dirichlet_condition_with_multipliers
      (md, mim, varname, multname, region, dataname);
  }

  size_type add_normal_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   dim_type degree, size_type region,
   const std::string &dataname) {
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    const mesh_fem &mf_mult = classical_mesh_fem(mf_u.linked_mesh(),degree, 1);
    return add_normal_Dirichlet_condition_with_multipliers
      (md, mim, varname, mf_mult, region, dataname);
  }

  size_type add_normal_Dirichlet_condition_with_penalization
  (model &md, const mesh_im &mim, const std::string &varname,
   scalar_type penalisation_coeff, size_type region,
   const std::string &dataname, const mesh_fem *mf_mult) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = std::make_shared<Dirichlet_condition_brick>
      (true, false, true, mf_mult);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  size_type add_generalized_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname, const std::string &Hname) {
    pbrick pbr = std::make_shared<Dirichlet_condition_brick>(false,true,false);
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    dl.push_back(dataname);
    dl.push_back(Hname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  size_type add_generalized_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const mesh_fem &mf_mult, size_type region,
   const std::string &dataname, const std::string &Hname) {
    std::string multname = md.new_name("mult_on_" + varname);
    md.add_multiplier(multname, mf_mult, varname);
    return add_generalized_Dirichlet_condition_with_multipliers
      (md, mim, varname, multname, region, dataname, Hname);
  }

  size_type add_generalized_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   dim_type degree, size_type region,
   const std::string &dataname, const std::string &Hname) {
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    const mesh_fem &mf_mult = classical_mesh_fem(mf_u.linked_mesh(),
                                                 degree, mf_u.get_qdim());
    return add_generalized_Dirichlet_condition_with_multipliers
      (md, mim, varname, mf_mult, region, dataname, Hname);
  }

  size_type add_generalized_Dirichlet_condition_with_penalization
  (model &md, const mesh_im &mim, const std::string &varname,
   scalar_type penalisation_coeff, size_type region,
   const std::string &dataname, const std::string &Hname,
   const mesh_fem *mf_mult) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = std::make_shared<Dirichlet_condition_brick>
      (true, true, false, mf_mult);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    dl.push_back(dataname);
    dl.push_back(Hname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  void change_penalization_coeff(model &md, size_type ind_brick,
                                 scalar_type penalisation_coeff) {
    const std::string &coeffname = md.dataname_of_brick(ind_brick, 0);
    if (!md.is_complex()) {
      model_real_plain_vector &d = md.set_real_variable(coeffname);
      GMM_ASSERT1(gmm::vect_size(d)==1,
                  "Wrong coefficient size, may be not a Dirichlet brick "
                  "with penalization");
      d[0] = penalisation_coeff;
    }
    else {
      model_complex_plain_vector &d = md.set_complex_variable(coeffname);
      GMM_ASSERT1(gmm::vect_size(d)==1,
                  "Wrong coefficient size, may be not a Dirichlet brick "
                  "with penalization");
      d[0] = penalisation_coeff;
    }
  }

  // ----------------------------------------------------------------------
  //
  // Dirichlet condition brick with simplification
  //
  // ----------------------------------------------------------------------

  struct simplification_Dirichlet_condition_brick : public virtual_brick {

    virtual void real_pre_assembly_in_serial(const model &md, size_type /*ib*/,
                                             const model::varnamelist &vl,
                                             const model::varnamelist &dl,
                                             const model::mimlist &mims,
                                             model::real_matlist &matl,
                                             model::real_veclist &vecl,
                                             model::real_veclist &,
                                             size_type region,
                                             build_version /*version*/) const {
      if (MPI_IS_MASTER()) {

        GMM_ASSERT1(vecl.size() == 0 && matl.size() == 0,
                    "Dirichlet condition brick by simplification has no term");
        GMM_ASSERT1(mims.size() == 0,
                    "Dirichlet condition brick by simplification need no "
                    "mesh_im");
        GMM_ASSERT1(vl.size() == 1 && dl.size() <= 1,
                    "Wrong number of variables for Dirichlet condition brick "
                    "by simplification");

        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        const model_real_plain_vector *A = 0;
        const mesh_fem *mf_data = 0;
        size_type s = 0;

        if (dl.size() == 1) {
          A = &(md.real_variable(dl[0]));
          mf_data = md.pmesh_fem_of_variable(dl[0]);

          if (mf_data) {
            GMM_ASSERT1(mf_data == &mf_u, "Sorry, for this brick, the data has"
                        " to be defined on the same f.e.m. as the unknown");
          } else {
            s = gmm::vect_size(*A);
            GMM_ASSERT1(mf_u.get_qdim() == s, ": bad format of "
                        "Dirichlet data. Detected dimension is " << s
                        << " should be " << size_type(mf_u.get_qdim()));
          }
        }

        mesh_region rg(region);
        // mf_u.linked_mesh().intersect_with_mpi_region(rg); // Not distributed

        if (mf_u.get_qdim() > 1 || (!mf_data && A)) {
          for (mr_visitor i(rg, mf_u.linked_mesh()); !i.finished(); ++i) {
            pfem pf = mf_u.fem_of_element(i.cv());
            if (pf) {
              GMM_ASSERT1(pf->target_dim() == 1,
                          "Intrinsically vectorial fems are not allowed");
              GMM_ASSERT1(mf_data || pf->is_lagrange(), "Constant Dirichlet "
                          "data allowed for lagrange fems only");
            }
          }
        }

        dal::bit_vector dofs = mf_u.dof_on_region(rg);

        if (A && !mf_data) {
          GMM_ASSERT1(dofs.card() % s == 0, "Problem with dof vectorization");
        }

        for (dal::bv_visitor i(dofs); !i.finished(); ++i) {
          scalar_type val(0);
          if (A) val = (mf_data ? (*A)[i] :  (*A)[i%s]);
          md.add_real_dof_constraint(vl[0], i, val);
        }
      }
    }

    virtual void complex_pre_assembly_in_serial(const model &md, size_type /*ib*/,
                                                const model::varnamelist &vl,
                                                const model::varnamelist &dl,
                                                const model::mimlist &mims,
                                                model::complex_matlist &matl,
                                                model::complex_veclist &vecl,
                                                model::complex_veclist &,
                                                size_type region,
                                                build_version /*version*/) const {
      if (MPI_IS_MASTER()) {
        GMM_ASSERT1(vecl.size() == 0 && matl.size() == 0,
                    "Dirichlet condition brick by simplification has no term");
        GMM_ASSERT1(mims.size() == 0,
                    "Dirichlet condition brick by simplification need no "
                    "mesh_im");
        GMM_ASSERT1(vl.size() == 1 && dl.size() <= 1,
                    "Wrong number of variables for Dirichlet condition brick "
                    "by simplification");

        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        const model_complex_plain_vector *A = 0;
        const mesh_fem *mf_data = 0;
        size_type s = 0;

        if (dl.size() == 1) {
          A = &(md.complex_variable(dl[0]));
          mf_data = md.pmesh_fem_of_variable(dl[0]);

          if (mf_data) {
            GMM_ASSERT1(mf_data == &mf_u, "Sorry, for this brick, the data has"
                        " to be define on the same f.e.m. than the unknown");
          } else {
            s = gmm::vect_size(*A);
            GMM_ASSERT1(mf_u.get_qdim() == s, ": bad format of "
                        "Dirichlet data. Detected dimension is " << s
                        << " should be " << size_type(mf_u.get_qdim()));
          }
        }

        mesh_region rg(region);
        // mf_u.linked_mesh().intersect_with_mpi_region(rg); // Not distributed

        if (mf_u.get_qdim() > 1 || (!mf_data && A)) {
          for (mr_visitor i(rg, mf_u.linked_mesh()); !i.finished(); ++i) {
            pfem pf = mf_u.fem_of_element(i.cv());
            if (pf) {
              GMM_ASSERT1(pf->target_dim() == 1,
                          "Intrinsically vectorial fems are not allowed");
              GMM_ASSERT1(mf_data || pf->is_lagrange(), "Constant Dirichlet "
                          "data allowed for lagrange fems only");
            }
          }
        }

        dal::bit_vector dofs = mf_u.dof_on_region(rg);

        if (A && !mf_data) {
          GMM_ASSERT1(dofs.card() % s == 0, "Problem with dof vectorization");
        }

        for (dal::bv_visitor i(dofs); !i.finished(); ++i) {
          complex_type val(0);
          if (A) val = (mf_data ? (*A)[i] :  (*A)[i%s]);
          md.add_complex_dof_constraint(vl[0], i, val);
        }
      }
    }


    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    simplification_Dirichlet_condition_brick() {
      set_flags("Dirichlet with simplification brick",
                true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                true /* compute each time */);
    }
  };

  size_type add_Dirichlet_condition_with_simplification
  (model &md, const std::string &varname,
   size_type region, const std::string &dataname) {
    pbrick pbr = std::make_shared<simplification_Dirichlet_condition_brick>();
    model::termlist tl;
    model::varnamelist vl(1, varname);
    model::varnamelist dl;
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), region);
  }

  // ----------------------------------------------------------------------
  //
  // Dirichlet condition brick with Nitsche's method
  //
  // ----------------------------------------------------------------------

  size_type add_Dirichlet_condition_with_Nitsche_method
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &Neumannterm,
   const std::string &datagamma0, size_type region, scalar_type theta_,
   const std::string &datag) {
    std::string theta = std::to_string(theta_);
    // reenables disabled variables
    ga_workspace workspace(md, ga_workspace::inherit::ALL);
    size_type order = workspace.add_expression(Neumannterm, mim, region, 1);
    GMM_ASSERT1(order == 0, "Wrong expression of the Neumann term");
    bool is_lin = workspace.is_linear(1);

    std::string condition = "("+varname + (datag.size() ? "-("+datag+"))":")");
    std::string gamma = "(("+datagamma0+")*element_size)";
    std::string r = "(1/"+gamma+")";
    std::string expr = "("+r+"*"+condition+"-("+Neumannterm+")).Test_"+varname;
    if (theta_ != scalar_type(0)) {
      std::string derivative_Neumann = workspace.extract_order1_term(varname);
      if (derivative_Neumann.size())
        expr+="-"+theta+"*"+condition+".("+derivative_Neumann+")";
    }

    // cout << "Nitsche expression : " << expr << endl;
    // cout << "is_lin : " << int(is_lin) << endl;

    if (is_lin) {
      return add_linear_term(md, mim, expr, region, false, false,
                             "Dirichlet condition with Nitsche's method");
    } else {
      return add_nonlinear_term(md, mim, expr, region, false, false,
                                "Dirichlet condition with Nitsche's method");
    }
  }

  size_type add_normal_Dirichlet_condition_with_Nitsche_method
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &Neumannterm,
   const std::string &datagamma0, size_type region, scalar_type theta_,
   const std::string &datag) {
    std::string theta = std::to_string(theta_);
    // reenables disabled variables
    ga_workspace workspace(md, ga_workspace::inherit::ALL);
    size_type order = workspace.add_expression(Neumannterm, mim, region, 1);
    GMM_ASSERT1(order == 0, "Wrong expression of the Neumann term");
    bool is_lin = workspace.is_linear(1);

    std::string condition = "("+varname+".Normal"
      + (datag.size() ? "-("+datag+"))":")");
    std::string gamma = "(("+datagamma0+")*element_size)";
    std::string r = "(1/"+gamma+")";
    std::string expr = "("+r+"*"+condition+"-Normal.("+Neumannterm
      +"))*(Normal.Test_"+varname+")";
    if (theta_ != scalar_type(0)) {
      std::string derivative_Neumann = workspace.extract_order1_term(varname);
      if (derivative_Neumann.size())
        expr+="-"+theta+"*"+condition+"*Normal.("
          +derivative_Neumann+")";
    }
    if (is_lin) {
      return add_linear_term(md, mim, expr, region, false, false,
                             "Dirichlet condition with Nitsche's method");
    } else {
      return add_nonlinear_term(md, mim, expr, region, false, false,
                                "Dirichlet condition with Nitsche's method");
    }
  }

  size_type add_generalized_Dirichlet_condition_with_Nitsche_method
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &Neumannterm,
   const std::string &datagamma0, size_type region, scalar_type theta_,
   const std::string &datag, const std::string &dataH) {
    std::string theta = std::to_string(theta_);
    // reenables disabled variables
    ga_workspace workspace(md, ga_workspace::inherit::ALL);
    size_type order = workspace.add_expression(Neumannterm, mim, region, 1);
    GMM_ASSERT1(order == 0, "Wrong expression of the Neumann term");
    bool is_lin = workspace.is_linear(1);

    std::string condition = "(("+dataH+")*"+varname
      + (datag.size() ? "-("+datag+"))":")");
    std::string gamma = "(("+datagamma0+")*element_size)";
    std::string r = "(1/"+gamma+")";
    std::string expr = "("+r+"*"+condition+"-("+dataH+")*("+Neumannterm
      +"))*(("+dataH+")*Test_"+varname+")";
    if (theta_ != scalar_type(0)) {
      std::string derivative_Neumann = workspace.extract_order1_term(varname);
      if (derivative_Neumann.size())
        expr+="-"+theta+"*"+condition+"*(("+dataH+")*("
          +derivative_Neumann+"))";
    }
    if (is_lin) {
      return add_linear_term(md, mim, expr, region, false, false,
                             "Dirichlet condition with Nitsche's method");
    } else {
      return add_nonlinear_term(md, mim, expr, region, false, false,
                                "Dirichlet condition with Nitsche's method");
    }
  }

  // ----------------------------------------------------------------------
  //
  // Pointwise constraints brick
  //
  // ----------------------------------------------------------------------
  // Two variables : with multipliers
  // One variable : penalization

  struct pointwise_constraints_brick : public virtual_brick {

    mutable gmm::row_matrix<model_real_sparse_vector> rB;
    mutable gmm::row_matrix<model_complex_sparse_vector> cB;

    virtual void  real_pre_assembly_in_serial(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &/*rvecl*/,
                                        size_type,
                                        build_version version) const {
      if (MPI_IS_MASTER()) {

        GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                    "Pointwize constraints brick has only one term");
        GMM_ASSERT1(mims.size() == 0,
                    "Pointwize constraints brick does not need a mesh_im");
        GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2,
                  "Wrong number of variables for pointwize constraints brick");
        bool penalized = (vl.size() == 1);
        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        dim_type N = mf_u.linked_mesh().dim(), Q = mf_u.get_qdim(), ind_pt = 0;
        size_type dlsize = size_type((penalized ? 1 : 0) + 1 + (Q > 1 ? 1 : 0));
        GMM_ASSERT1(dl.size() == dlsize || dl.size() == dlsize+1,
                    "Wrong number of data for pointwize constraints brick");


        const model_real_plain_vector *COEFF = 0;
        if (penalized) {
          COEFF = &(md.real_variable(dl[0]));
          ind_pt = 1;
          GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                      "Data for coefficient should be a scalar");
        }

        const model_real_plain_vector &PT = md.real_variable(dl[ind_pt]);
        size_type nb_co = gmm::vect_size(PT) / N;

        dim_type ind_unitv = dim_type((Q > 1) ? ind_pt+1 : 0);
        const model_real_plain_vector &unitv =md.real_variable(dl[ind_unitv]);
        GMM_ASSERT1((!ind_unitv || gmm::vect_size(unitv) == nb_co * Q),
                    "Wrong size for vector of unit vectors");

        dim_type ind_rhs = dim_type((Q > 1) ? ind_pt+2 : ind_pt+1);
        if (dl.size() < size_type(ind_rhs + 1)) ind_rhs = 0;
        const model_real_plain_vector &rhs =  md.real_variable(dl[ind_rhs]);
        GMM_ASSERT1((!ind_rhs || gmm::vect_size(rhs) == nb_co),
                    "Wrong size for vector of rhs");

        bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
          || (penalized && (md.is_var_newer_than_brick(dl[ind_pt], ib)
                            || md.is_var_newer_than_brick(dl[ind_unitv], ib)
                            || md.is_var_newer_than_brick(dl[ind_rhs], ib)));

        if (recompute_matrix) {
          gmm::row_matrix<model_real_sparse_vector> BB(nb_co*Q, mf_u.nb_dof());
          gmm::clear(rB); gmm::resize(rB, nb_co, mf_u.nb_dof());

          dal::bit_vector dof_untouched;
          getfem::mesh_trans_inv mti(mf_u.linked_mesh());
          base_node pt(N);
          for (size_type i = 0; i < nb_co; ++i) {
            gmm::copy(gmm::sub_vector(PT, gmm::sub_interval(i*N, N)), pt);
            mti.add_point(pt);
          }
          gmm::row_matrix<model_real_sparse_vector> &BBB = ((Q > 1) ? BB : rB);
          model_real_plain_vector vv;
          interpolation(mf_u, mti, vv, vv, BBB,  1, 1, &dof_untouched);
          GMM_ASSERT1(dof_untouched.card() == 0,
                      "Pointwize constraints : some of the points are outside "
                      "the mesh: " << dof_untouched);

          if (Q > 1) {
            for (size_type i = 0; i < nb_co; ++i)
              for (size_type q = 0; q < Q; ++q)
                gmm::add(gmm::scaled(gmm::mat_row(BB, i*Q+q), unitv[i*Q+q]),
                         gmm::mat_row(rB, i));
          }
          if (penalized) {
          gmm::mult(gmm::transposed(rB), rB, matl[0]);
          gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
          } else
            gmm::copy(rB, matl[0]);
        }

        if (ind_rhs) {
          if (penalized) {
            gmm::mult(gmm::transposed(rB), rhs, vecl[0]);
            gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
          }
          else gmm::copy(rhs, vecl[0]);
        }
        else gmm::clear(vecl[0]);
      }
    }

    virtual void complex_pre_assembly_in_serial(const model &md, size_type ib,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type,
                                           build_version version) const {
      if (MPI_IS_MASTER()) {
        GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                    "Pointwize constraints brick only one term");
        GMM_ASSERT1(mims.size() == 0,
                    "Pointwize constraints brick does not need a mesh_im");
        GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2,
                  "Wrong number of variables for pointwize constraints brick");
        bool penalized = (vl.size() == 1);
        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        dim_type N = mf_u.linked_mesh().dim(), Q = mf_u.get_qdim(), ind_pt = 0;
        size_type dlsize = size_type((penalized ? 1 : 0) + 1 + (Q > 1 ? 1 :0));
        GMM_ASSERT1(dl.size() == dlsize || dl.size() == dlsize+1,
                    "Wrong number of data for pointwize constraints brick");


        const model_complex_plain_vector *COEFF = 0;
        if (penalized) {
          COEFF = &(md.complex_variable(dl[0]));
          ind_pt = 1;
          GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                      "Data for coefficient should be a scalar");
        }

        const model_complex_plain_vector &PT = md.complex_variable(dl[ind_pt]);
        size_type nb_co = gmm::vect_size(PT) / N;

        dim_type ind_unitv = dim_type((Q > 1) ? ind_pt+1 : 0);
        const model_complex_plain_vector &unitv
          = md.complex_variable(dl[ind_unitv]);
        GMM_ASSERT1((!ind_unitv || gmm::vect_size(unitv) == nb_co * Q),
                    "Wrong size for vector of unit vectors");

        dim_type ind_rhs = dim_type((Q > 1) ? ind_pt+2 : ind_pt+1);
        if (dl.size() < size_type(ind_rhs + 1)) ind_rhs = 0;
        const model_complex_plain_vector &rhs
          = md.complex_variable(dl[ind_rhs]);
        GMM_ASSERT1((!ind_rhs || gmm::vect_size(rhs) == nb_co),
                    "Wrong size for vector of rhs");

        bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
          || (penalized && (md.is_var_newer_than_brick(dl[ind_pt], ib)
                            || md.is_var_newer_than_brick(dl[ind_unitv], ib)
                            || md.is_var_newer_than_brick(dl[ind_rhs], ib)));

        if (recompute_matrix) {
          gmm::row_matrix<model_complex_sparse_vector>
            BB(nb_co*Q,mf_u.nb_dof());
          gmm::clear(cB); gmm::resize(cB, nb_co, mf_u.nb_dof());
          dal::bit_vector dof_untouched;
          getfem::mesh_trans_inv mti(mf_u.linked_mesh());
          base_node pt(N);
          for (size_type i = 0; i < nb_co; ++i) {
            gmm::copy(gmm::real_part
                      (gmm::sub_vector(PT, gmm::sub_interval(i*N, N))), pt);
            mti.add_point(pt);
          }
          gmm::row_matrix<model_complex_sparse_vector> &BBB = ((Q > 1) ? BB :cB);
          model_complex_plain_vector vv;
          interpolation(mf_u, mti, vv, vv, BBB,  1, 1, &dof_untouched);
          GMM_ASSERT1(dof_untouched.card() == 0,
                      "Pointwize constraints : some of the points are outside "
                      "the mesh: " << dof_untouched);

          if (Q > 1) {
            for (size_type i = 0; i < nb_co; ++i)
              for (size_type q = 0; q < Q; ++q)
                gmm::add(gmm::scaled(gmm::mat_row(BB, i*Q+q), unitv[i*Q+q]),
                         gmm::mat_row(cB, i));
          }

          if (penalized) {
            gmm::mult(gmm::transposed(cB), cB, matl[0]);
            gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
          } else
            gmm::copy(cB, matl[0]);
        }


        if (ind_rhs) {
          if (penalized) {
            gmm::mult(gmm::transposed(cB), rhs, vecl[0]);
            gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
          }
          else gmm::copy(rhs, vecl[0]);
        }
        else gmm::clear(vecl[0]);
      }
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    pointwise_constraints_brick(bool penalized) {
      set_flags(penalized ? "Pointwise cosntraints with penalization brick"
                          : "Pointwise cosntraints with multipliers brick",
                true /* is linear*/,
                true /* is symmetric */, penalized /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }
  };


  size_type add_pointwise_constraints_with_penalization
  (model &md, const std::string &varname,
   scalar_type penalisation_coeff, const std::string &dataname_pt,
   const std::string &dataname_unitv, const std::string &dataname_val) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = std::make_shared<pointwise_constraints_brick>(true);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    dl.push_back(dataname_pt);
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    if (mf_u.get_qdim() > 1) dl.push_back(dataname_unitv);
    if (dataname_val.size() > 0) dl.push_back(dataname_val);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }

  size_type add_pointwise_constraints_with_given_multipliers
  (model &md, const std::string &varname,
   const std::string &multname, const std::string &dataname_pt,
   const std::string &dataname_unitv, const std::string &dataname_val) {
    pbrick pbr = std::make_shared<pointwise_constraints_brick>(false);
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl(1, dataname_pt);
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    if (mf_u.get_qdim() > 1) dl.push_back(dataname_unitv);
    if (dataname_val.size() > 0) dl.push_back(dataname_val);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }

  size_type add_pointwise_constraints_with_multipliers
  (model &md, const std::string &varname, const std::string &dataname_pt,
   const std::string &dataname_unitv, const std::string &dataname_val) {
    std::string multname = md.new_name("mult_on_" + varname);
    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    size_type nb_co =
      ((md.is_complex()) ? gmm::vect_size(md.complex_variable(dataname_pt))
       : gmm::vect_size(md.real_variable(dataname_pt)))
      / mf_u.linked_mesh().dim();
    md.add_fixed_size_variable(multname, nb_co);
    return add_pointwise_constraints_with_given_multipliers
      (md, varname, multname, dataname_pt, dataname_unitv, dataname_val);
  }


  // ----------------------------------------------------------------------
  //
  // Helmholtz brick
  //
  // ----------------------------------------------------------------------

  struct Helmholtz_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Helmholtz brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Helmholtz brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
                  "Wrong number of variables for Helmholtz brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      size_type Q = mf_u.get_qdim(), s = 1;
      GMM_ASSERT1(Q == 1, "Helmholtz brick is only for scalar field, sorry.");
      const mesh_im &mim = *mims[0];
      const mesh_fem *mf_a = 0;
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);
      const model_real_plain_vector *A = &(md.real_variable(dl[0]));
      mf_a = md.pmesh_fem_of_variable(dl[0]);
      s = gmm::vect_size(*A);
      if (mf_a) s = s * mf_a->get_qdim() / mf_a->nb_dof();

      if (s == 1) {
        GMM_TRACE2("Stiffness matrix assembly for Helmholtz problem");
        gmm::clear(matl[0]);
        model_real_plain_vector A2(gmm::vect_size(*A));
        for (size_type i=0; i < gmm::vect_size(*A); ++i) // Not valid for
          A2[i] = gmm::sqr((*A)[i]); // non lagrangian fem ...
        if (mf_a)
          asm_Helmholtz(matl[0], mim, mf_u, *mf_a, A2, rg);
        else
          asm_homogeneous_Helmholtz(matl[0], mim, mf_u, A2, rg);
      } else
        GMM_ASSERT1(false, "Bad format Helmholtz brick coefficient");
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Helmholtz brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Helmholtz brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
                  "Wrong number of variables for Helmholtz brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      size_type Q = mf_u.get_qdim(), s = 1;
      GMM_ASSERT1(Q == 1, "Helmholtz brick is only for scalar field, sorry.");
      const mesh_im &mim = *mims[0];
      const mesh_fem *mf_a = 0;
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);
      const model_complex_plain_vector *A = &(md.complex_variable(dl[0]));
      mf_a = md.pmesh_fem_of_variable(dl[0]);
      s = gmm::vect_size(*A);
      if (mf_a) s = s * mf_a->get_qdim() / mf_a->nb_dof();

      if (s == 1) {
        GMM_TRACE2("Stiffness matrix assembly for Helmholtz problem");
        gmm::clear(matl[0]);
        model_complex_plain_vector A2(gmm::vect_size(*A));
        for (size_type i=0; i < gmm::vect_size(*A); ++i) // Not valid for
          A2[i] = gmm::sqr((*A)[i]); // non lagrangian fem ...
        if (mf_a)
          asm_Helmholtz(matl[0], mim, mf_u, *mf_a, A2, rg);
        else
          asm_homogeneous_Helmholtz(matl[0], mim, mf_u, A2, rg);
      } else
        GMM_ASSERT1(false, "Bad format Helmholtz brick coefficient");
    }

    Helmholtz_brick() {
      set_flags("Helmholtz", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */);
    }

  };

  size_type add_Helmholtz_brick(model &md, const mesh_im &mim,
                                const std::string &varname,
                                const std::string &dataexpr,
                                size_type region) {
    if (md.is_complex()) {
      pbrick pbr = std::make_shared<Helmholtz_brick>();
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          model::varnamelist(1, dataexpr), tl,
                          model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      std::string expr = "Grad_"+varname+".Grad_"+test_varname
        +" + sqr("+dataexpr+")*"+varname+"*"+test_varname;

       size_type ib = add_linear_term(md, mim, expr, region, true, true,
                                      "Helmholtz", true);
       if (ib == size_type(-1))
         ib = add_nonlinear_term(md, mim, expr, region, false, false,
                                 "Helmholtz (nonlinear)");
       return ib;
    }
  }



  // ----------------------------------------------------------------------
  //
  // Fourier-Robin brick
  //
  // ----------------------------------------------------------------------

  struct Fourier_Robin_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Fourier-Robin brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Fourier-Robin brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
                  "Wrong number of variables for Fourier-Robin brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      size_type Q = mf_u.get_qdim(), s = 1;
      const mesh_im &mim = *mims[0];
      const mesh_fem *mf_a = 0;
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);
      const model_real_plain_vector *A = &(md.real_variable(dl[0]));
      mf_a = md.pmesh_fem_of_variable(dl[0]);
      s = gmm::vect_size(*A);
      if (mf_a) s = s * mf_a->get_qdim() / mf_a->nb_dof();
      GMM_ASSERT1(s == Q*Q,
                  "Bad format Fourier-Robin brick coefficient");

      GMM_TRACE2("Fourier-Robin term assembly");
      gmm::clear(matl[0]);
      if (mf_a)
        asm_qu_term(matl[0], mim, mf_u, *mf_a, *A, rg);
      else
        asm_homogeneous_qu_term(matl[0], mim, mf_u, *A, rg);
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Fourier-Robin brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Fourier-Robin brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 1,
                  "Wrong number of variables for Fourier-Robin brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      size_type Q = mf_u.get_qdim(), s = 1;
      const mesh_im &mim = *mims[0];
      const mesh_fem *mf_a = 0;
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);
      const model_complex_plain_vector *A = &(md.complex_variable(dl[0]));
      mf_a = md.pmesh_fem_of_variable(dl[0]);
      s = gmm::vect_size(*A);
      if (mf_a) s = s * mf_a->get_qdim() / mf_a->nb_dof();
      GMM_ASSERT1(s == Q*Q,
                  "Bad format Fourier-Robin brick coefficient");

      GMM_TRACE2("Fourier-Robin term assembly");
      gmm::clear(matl[0]);
      if (mf_a)
        asm_qu_term(matl[0], mim, mf_u, *mf_a, *A, rg);
      else
        asm_homogeneous_qu_term(matl[0], mim, mf_u, *A, rg);
    }

    Fourier_Robin_brick() {
      set_flags("Fourier Robin condition", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }

  };

  size_type add_Fourier_Robin_brick(model &md, const mesh_im &mim,
                                    const std::string &varname,
                                    const std::string &dataexpr,
                                    size_type region) {
    if (md.is_complex()) {
      pbrick pbr = std::make_shared<Fourier_Robin_brick>();
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          model::varnamelist(1, dataexpr), tl,
                          model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      std::string expr = "(("+dataexpr+")*"+varname+")."+test_varname;
      size_type ib = add_linear_term(md, mim, expr, region, true, true,
                                     "Fourier-Robin", true);
      if (ib == size_type(-1))
        ib = add_nonlinear_term(md, mim, expr, region, false, false,
                                "Fourier-Robin (nonlinear)");
      return ib;
    }
  }

  // ----------------------------------------------------------------------
  //
  // Constraint brick
  //
  // ----------------------------------------------------------------------

  struct have_private_data_brick : public virtual_brick {

    model_real_sparse_matrix rB;
    model_complex_sparse_matrix cB;
    model_real_plain_vector rL;
    model_complex_plain_vector cL;
    std::string nameL;
  };

  struct constraint_brick : public have_private_data_brick {

    virtual void real_pre_assembly_in_serial(const model &md, size_type,
                                             const model::varnamelist &vl,
                                             const model::varnamelist &dl,
                                             const model::mimlist &mims,
                                             model::real_matlist &matl,
                                             model::real_veclist &vecl,
                                             model::real_veclist &,
                                             size_type, build_version) const {
      if (MPI_IS_MASTER()) {

        GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                    "Constraint brick has one and only one term");
        GMM_ASSERT1(mims.size() == 0,
                    "Constraint brick need no mesh_im");
        GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 1,
                    "Wrong number of variables for constraint brick");

        bool penalized = (vl.size() == 1);
        const model_real_plain_vector *COEFF = 0;

        bool has_data = (nameL.compare("") != 0);
        if (has_data)
          GMM_ASSERT1(nameL.compare(dl.back()) == 0 &&
                      md.variable_exists(nameL) && md.is_data(nameL),
                      "Internal error");
        const model_real_plain_vector &
          rrL = has_data ? md.real_variable(nameL) : rL;

        if (penalized) {
          COEFF = &(md.real_variable(dl[0]));
          GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                      "Data for coefficient should be a scalar");

          gmm::mult(gmm::transposed(rB),
                    gmm::scaled(rrL, gmm::abs((*COEFF)[0])), vecl[0]);
          gmm::mult(gmm::transposed(rB),
                    gmm::scaled(rB, gmm::abs((*COEFF)[0])), matl[0]);
        } else {
          gmm::copy(rrL, vecl[0]);
          gmm::copy(rB, matl[0]);
        }
      }
    }

    virtual void complex_pre_assembly_in_serial(const model &md, size_type,
                                                const model::varnamelist &vl,
                                                const model::varnamelist &dl,
                                                const model::mimlist &mims,
                                                model::complex_matlist &matl,
                                                model::complex_veclist &vecl,
                                                model::complex_veclist &,
                                                size_type,
                                                build_version) const {
      if (MPI_IS_MASTER()) {

        GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                    "Constraint brick has one and only one term");
        GMM_ASSERT1(mims.size() == 0,
                    "Constraint brick need no mesh_im");
        GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 1,
                    "Wrong number of variables for constraint brick");

        bool penalized = (vl.size() == 1);
        const model_complex_plain_vector *COEFF = 0;

        bool has_data = (nameL.compare("") != 0);
        if (has_data)
          GMM_ASSERT1(nameL.compare(dl.back()) == 0 &&
                      md.variable_exists(nameL) && md.is_data(nameL),
                      "Internal error");
        const model_complex_plain_vector &
          ccL = has_data ? md.complex_variable(nameL) : cL;

        if (penalized) {
          COEFF = &(md.complex_variable(dl[0]));
          GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                      "Data for coefficient should be a scalar");

          gmm::mult(gmm::transposed(cB),
                    gmm::scaled(ccL, gmm::abs((*COEFF)[0])), vecl[0]);
          gmm::mult(gmm::transposed(cB),
                    gmm::scaled(cB, gmm::abs((*COEFF)[0])), matl[0]);
        } else {
          gmm::copy(ccL, vecl[0]);
          gmm::copy(cB, matl[0]);
        }
      }
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    constraint_brick(bool penalized) {
      set_flags(penalized ? "Constraint with penalization brick"
                          : "Constraint with multipliers brick",
                true /* is linear*/,
                true /* is symmetric */, penalized /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }

  };

  model_real_sparse_matrix &set_private_data_brick_real_matrix
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    have_private_data_brick *p = dynamic_cast<have_private_data_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->rB;
  }

  model_real_plain_vector &set_private_data_brick_real_rhs
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    have_private_data_brick *p = dynamic_cast<have_private_data_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    if (p->nameL.compare("") != 0) GMM_WARNING1("Rhs already set by data name");
    return p->rL;
  }

  model_complex_sparse_matrix &set_private_data_brick_complex_matrix
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    have_private_data_brick *p = dynamic_cast<have_private_data_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    return p->cB;
  }

  model_complex_plain_vector &set_private_data_brick_complex_rhs
  (model &md, size_type indbrick) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    have_private_data_brick *p = dynamic_cast<have_private_data_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    if (p->nameL.compare("") != 0) GMM_WARNING1("Rhs already set by data name");
    return p->cL;
  }

  void set_private_data_rhs
  (model &md, size_type indbrick, const std::string &varname) {
    pbrick pbr = md.brick_pointer(indbrick);
    md.touch_brick(indbrick);
    have_private_data_brick *p = dynamic_cast<have_private_data_brick *>
      (const_cast<virtual_brick *>(pbr.get()));
    GMM_ASSERT1(p, "Wrong type of brick");
    if (p->nameL.compare(varname) != 0) {
      model::varnamelist dl = md.datanamelist_of_brick(indbrick);
      if (p->nameL.compare("") == 0) dl.push_back(varname);
      else dl.back() = varname;
      md.change_data_of_brick(indbrick, dl);
      p->nameL = varname;
    }
  }

  size_type add_constraint_with_penalization
  (model &md, const std::string &varname, scalar_type penalisation_coeff) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = std::make_shared<constraint_brick>(true);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }

  size_type add_constraint_with_multipliers
  (model &md, const std::string &varname, const std::string &multname) {
    pbrick pbr = std::make_shared<constraint_brick>(false);
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }


  // ----------------------------------------------------------------------
  //
  // Explicit matrix brick
  //
  // ----------------------------------------------------------------------

  struct explicit_matrix_brick : public have_private_data_brick {

    virtual void real_pre_assembly_in_serial(const model &, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type, build_version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Explicit matrix has one and only one term");
      GMM_ASSERT1(mims.size() == 0, "Explicit matrix need no mesh_im");
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() == 0,
                  "Wrong number of variables for explicit matrix brick");
      GMM_ASSERT1(gmm::mat_ncols(rB) == gmm::mat_ncols(matl[0]) &&
                  gmm::mat_nrows(rB) == gmm::mat_nrows(matl[0]),
                  "Explicit matrix brick dimension mismatch ("<<
                  gmm::mat_ncols(rB)<<"x"<<gmm::mat_nrows(rB)<<") != ("<<
                  gmm::mat_ncols(matl[0])<<"x"<<gmm::mat_nrows(matl[0])<<")");
      gmm::copy(rB, matl[0]);
    }

    virtual void complex_pre_assembly_in_serial(const model &, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type,
                                           build_version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Explicit matrix has one and only one term");
      GMM_ASSERT1(mims.size() == 0, "Explicit matrix need no mesh_im");
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() == 0,
                  "Wrong number of variables for explicit matrix brick");
      gmm::copy(cB, matl[0]);
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    explicit_matrix_brick(bool symmetric_, bool coercive_) {
      set_flags("Explicit matrix brick",
                true /* is linear*/,
                symmetric_ /* is symmetric */, coercive_ /* is coercive */,
                true /* is real */, true /* is complex */,
                true /* is to be computed each time */);
    }
  };

  size_type add_explicit_matrix
  (model &md, const std::string &varname1, const std::string &varname2,
   bool issymmetric, bool iscoercive) {
    pbrick pbr = std::make_shared<explicit_matrix_brick>(issymmetric,
                                                         iscoercive);
    model::termlist tl;
    tl.push_back(model::term_description(varname1, varname2, issymmetric));
    model::varnamelist vl(1, varname1);
    vl.push_back(varname2);
    model::varnamelist dl;
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }

  // ----------------------------------------------------------------------
  //
  // Explicit rhs brick
  //
  // ----------------------------------------------------------------------

  struct explicit_rhs_brick : public have_private_data_brick {

    virtual void real_pre_assembly_in_serial(const model &, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type, build_version) const {
      if (MPI_IS_MASTER()) {
        GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                    "Explicit rhs has one and only one term");
        GMM_ASSERT1(mims.size() == 0, "Explicit rhs need no mesh_im");
        GMM_ASSERT1(vl.size() == 1 && dl.size() == 0,
                    "Wrong number of variables for explicit rhs brick");
        gmm::copy(rL, vecl[0]);
      }
    }

    virtual void complex_pre_assembly_in_serial(const model &, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type,
                                           build_version) const {
      if (MPI_IS_MASTER()) {
        GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                    "Explicit rhs has one and only one term");
        GMM_ASSERT1(mims.size() == 0, "Explicit rhs need no mesh_im");
        GMM_ASSERT1(vl.size() == 1 && dl.size() == 0,
                    "Wrong number of variables for explicit rhs brick");
        gmm::copy(cL, vecl[0]);
      }

    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    explicit_rhs_brick() {
      set_flags("Explicit rhs brick",
                true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                true /* is to be computed each time */);
    }

  };

  size_type add_explicit_rhs
  (model &md, const std::string &varname) {
    pbrick pbr = std::make_shared<explicit_rhs_brick>();
    model::termlist tl;
    tl.push_back(model::term_description(varname));
    model::varnamelist vl(1, varname);
    model::varnamelist dl;
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }


  // ----------------------------------------------------------------------
  //
  // Isotropic linearized elasticity brick
  //
  // ----------------------------------------------------------------------

  struct iso_lin_elasticity_new_brick : public virtual_brick {

    std::string expr, dataname3;

    void asm_real_tangent_terms(const model &md, size_type ib,
                                const model::varnamelist &vl,
                                const model::varnamelist &dl,
                                const model::mimlist &mims,
                                model::real_matlist &matl,
                                model::real_veclist &vecl,
                                model::real_veclist &,
                                size_type region,
                                build_version version) const override {
      GMM_ASSERT1(vl.size() == 1, "Linearized isotropic elasticity brick "
                  "has one and only one variable");
      GMM_ASSERT1(matl.size() == 1, "Linearized isotropic elasticity brick "
                  "has one and only one term");
      GMM_ASSERT1(mims.size() == 1, "Linearized isotropic elasticity brick "
                  "needs one and only one mesh_im");

      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0);
      for (size_type i = 0; i < dl.size(); ++i) {
        recompute_matrix = recompute_matrix ||
          md.is_var_newer_than_brick(dl[i], ib);
      }

      if (recompute_matrix) {
        // reenables disabled variables
        ga_workspace workspace(md, ga_workspace::inherit::ALL);
        workspace.add_expression(expr, *(mims[0]), region);
        GMM_TRACE2(name << ": generic matrix assembly");
        workspace.assembly(2);
        scalar_type alpha = scalar_type(1)
          / (workspace.factor_of_variable(vl[0]));
        const auto &R=workspace.assembled_matrix();
        gmm::sub_interval I = workspace.interval_of_variable(vl[0]);
        gmm::copy(gmm::scaled(gmm::sub_matrix(R, I, I), alpha),
                  matl[0]);
      }

      if  (dataname3.size()) { // Pre-constraints given by an "initial"
        // displacement u0. Means that the computed displacement will be u - u0
        // The displacement u0 should be discribed on the same fem as the
        // variable.
        gmm::clear(vecl[0]);
        gmm::mult(matl[0],
                  gmm::scaled(md.real_variable(dataname3), scalar_type(-1)),
                  vecl[0]);
      }

    }

    void real_post_assembly_in_serial(const model &md, size_type ib,
                                      const model::varnamelist &/* vl */,
                                      const model::varnamelist &/* dl */,
                                      const model::mimlist &/* mims */,
                                      model::real_matlist &/*matl*/,
                                      model::real_veclist &vecl,
                                      model::real_veclist &,
                                      size_type /*region*/,
                                      build_version) const override {
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }


    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return expr;
    }

    iso_lin_elasticity_new_brick(const std::string &expr_,
                                 const std::string &dataname3_) {
      expr = expr_; dataname3 = dataname3_;
      set_flags("Linearized isotropic elasticity", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };


  size_type add_isotropic_linearized_elasticity_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataexpr1, const std::string &dataexpr2,
   size_type region, const std::string &dataname3) {
    std::string test_varname
      = "Test_" + sup_previous_and_dot_to_varname(varname);

    std::string expr1 = "((("+dataexpr1+")*(Div_"+varname+"-Div_"+dataname3
      +"))*Id(meshdim)+(2*("+dataexpr2+"))*(Sym(Grad_"+varname
      +")-Sym(Grad_"+dataname3+"))):Grad_" +test_varname;
    std::string expr2 = "(Div_"+varname+"*(("+dataexpr1+")*Id(meshdim))"
      +"+(2*("+dataexpr2+"))*Sym(Grad_"+varname+")):Grad_"+test_varname;

    bool is_lin;
    model::varnamelist vl, dl;
    { // reenables disabled variables
      ga_workspace workspace(md, ga_workspace::inherit::ALL);
      workspace.add_expression(expr2, mim, region);
      model::varnamelist vl_test1, vl_test2;
      is_lin = workspace.used_variables(vl, vl_test1, vl_test2, dl, 2);
    }
    if (is_lin) {
      pbrick pbr = std::make_shared<iso_lin_elasticity_new_brick>
        (expr2, dataname3);
      model::termlist tl;
      tl.push_back(model::term_description(varname,
                           sup_previous_and_dot_to_varname(varname), true));
      if (dataname3.size()) dl.push_back(dataname3);
      return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
    } else {
      return add_nonlinear_term
        (md, mim, dataname3.size() ? expr1 : expr2, region, false, false,
         "Linearized isotropic elasticity (with nonlinear dependance)");
    }
  }

  size_type add_isotropic_linearized_elasticity_pstrain_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &data_E, const std::string &data_nu,
   size_type region) {
    std::string test_varname
      = "Test_" + sup_previous_and_dot_to_varname(varname);

    std::string mu = "(("+data_E+")/(2*(1+("+data_nu+"))))";
    std::string lambda = "(("+data_E+")*("+data_nu+")/((1+("+data_nu+"))*(1-2*("
      +data_nu+"))))";
    std::string expr = lambda+"*Div_"+varname+"*Div_"+test_varname
      + "+"+mu+"*(Grad_"+varname+"+Grad_"+varname+"'):Grad_"+test_varname;

    bool is_lin;
    { // reenables disabled variables
      ga_workspace workspace(md, ga_workspace::inherit::ALL);
      workspace.add_expression(expr, mim, region);
      is_lin = workspace.is_linear(2);
    }
    if (is_lin) {
      return add_linear_term(md, mim, expr, region, false, false,
                             "Linearized isotropic elasticity");
    } else {
      return add_nonlinear_term
        (md, mim, expr, region, false, false,
         "Linearized isotropic elasticity (with nonlinear dependance)");
    }
  }

  size_type add_isotropic_linearized_elasticity_pstress_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &data_E, const std::string &data_nu,
   size_type region) {
    std::string test_varname
      = "Test_" + sup_previous_and_dot_to_varname(varname);

    const mesh_fem *mfu = md.pmesh_fem_of_variable(varname);
    GMM_ASSERT1(mfu, "The variable should be a fem variable");
    size_type N = mfu->linked_mesh().dim();

    std::string mu = "(("+data_E+")/(2*(1+("+data_nu+"))))";
    std::string lambda =  "(("+data_E+")*("+data_nu+")/((1+("+data_nu+"))*(1-2*("
      +data_nu+"))))";
    if (N == 2)
      lambda = "(("+data_E+")*("+data_nu+")/((1-sqr("+data_nu+"))))";
    std::string expr = lambda+"*Div_"+varname+"*Div_"+test_varname
      + "+"+mu+"*(Grad_"+varname+"+Grad_"+varname+"'):Grad_"+test_varname;

    bool is_lin;
    { // reenables disabled variables
      ga_workspace workspace(md, ga_workspace::inherit::ALL);
      workspace.add_expression(expr, mim, region);
      is_lin = workspace.is_linear(2);
    }
    if (is_lin) {
      return add_linear_term(md, mim, expr, region, false, false,
                             "Linearized isotropic elasticity");
    } else {
      return add_nonlinear_term
        (md, mim, expr, region, false, false,
         "Linearized isotropic elasticity (with nonlinear dependance)");
    }
  }

  // Tresca to be implemented with generic interpolation
  void compute_isotropic_linearized_Von_Mises_or_Tresca
  (model &md, const std::string &varname, const std::string &data_lambda,
   const std::string &data_mu, const mesh_fem &mf_vm,
   model_real_plain_vector &VM, bool tresca) {

    if (tresca) {
      const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
      const mesh_fem *mf_lambda = md.pmesh_fem_of_variable(data_lambda);
      const model_real_plain_vector *lambda=&(md.real_variable(data_lambda));
      const mesh_fem *mf_mu = md.pmesh_fem_of_variable(data_mu);
      const model_real_plain_vector *mu = &(md.real_variable(data_mu));

      size_type sl = gmm::vect_size(*lambda);
      if (mf_lambda) sl = sl * mf_lambda->get_qdim() / mf_lambda->nb_dof();
      size_type sm = gmm::vect_size(*mu);
      if (mf_mu) sm = sm * mf_mu->get_qdim() / mf_mu->nb_dof();

      GMM_ASSERT1(sl == 1 && sm == 1, "Bad format for Lame coefficients");
      GMM_ASSERT1(mf_lambda == mf_mu,
                  "The two Lame coefficients should be described on the same "
                  "finite element method.");

      if (mf_lambda) {
        getfem::interpolation_von_mises_or_tresca(mf_u, mf_vm,
                                                  md.real_variable(varname), VM,
                                                  *mf_lambda, *lambda,
                                                  *mf_lambda, *mu,
                                                  tresca);
      } else {
        mf_lambda = &(classical_mesh_fem(mf_u.linked_mesh(), 0));
        model_real_plain_vector LAMBDA(mf_lambda->nb_dof(), (*lambda)[0]);
        model_real_plain_vector MU(mf_lambda->nb_dof(), (*mu)[0]);
        getfem::interpolation_von_mises_or_tresca(mf_u, mf_vm,
                                                  md.real_variable(varname), VM,
                                                  *mf_lambda, LAMBDA,
                                                  *mf_lambda, MU,
                                                  tresca);
      }
    } else {
      // The Lambda part is not necessary for Von Mises stress ...
      // std::string sigma = "("+data_lambda+")*Div_"+varname+"*Id(meshdim)+("
      //   + data_mu+")*(Grad_"+varname+"+Grad_"+varname+"')";
      std::string sigma_d="("+data_mu+")*(Grad_"+varname+"+Grad_"+varname+"')";
      std::string expr = "sqrt(3/2)*Norm(Deviator("+sigma_d+"))";
      ga_interpolation_Lagrange_fem(md, expr, mf_vm, VM);
    }
  }


  void compute_isotropic_linearized_Von_Mises_pstrain
  (model &md, const std::string &varname, const std::string &data_E,
   const std::string &data_nu, const mesh_fem &mf_vm,
   model_real_plain_vector &VM) {
    // The Lambda part is not necessary for Von Mises stress ...
    // std::string lambda = "(("+data_E+")*("+data_nu+")/((1+("+data_nu+"))*(1-2*("
    //  +data_nu+"))))";
    std::string mu = "(("+data_E+")/(2*(1+("+data_nu+"))))";
    // std::string sigma = lambda+"*Div_"+varname+"*Id(meshdim)+"
    //  + mu+"*(Grad_"+varname+"+Grad_"+varname+"')";
    std::string sigma_d = mu+"*(Grad_"+varname+"+Grad_"+varname+"')";
    std::string expr = "sqrt(3/2)*Norm(Deviator("+sigma_d+"))";
    ga_interpolation_Lagrange_fem(md, expr, mf_vm, VM);
  }

  void compute_isotropic_linearized_Von_Mises_pstress
  (model &md, const std::string &varname, const std::string &data_E,
   const std::string &data_nu, const mesh_fem &mf_vm,
   model_real_plain_vector &VM) {
    // The Lambda part is not necessary for Von Mises stress ...
    // const mesh_fem *mfu = md.pmesh_fem_of_variable(varname);
    // GMM_ASSERT1(mfu, "The variable should be a fem variable");
    // size_type N = mfu->linked_mesh().dim();
    // std::string lambda =  "(("+data_E+")*("+data_nu+")/((1+("+data_nu
    //  +"))*(1-2*("+data_nu+"))))";
    // if (N == 2)
    //  lambda = "(("+data_E+")*("+data_nu+")/((1-sqr("+data_nu+"))))";
    std::string mu = "(("+data_E+")/(2*(1+("+data_nu+"))))";
    // std::string sigma = lambda+"*Div_"+varname+"*Id(meshdim)+"
    //   + mu+"*(Grad_"+varname+"+Grad_"+varname+"')";
    std::string sigma_d = mu+"*(Grad_"+varname+"+Grad_"+varname+"')";
    std::string expr = "sqrt(3/2)*Norm(Deviator("+sigma_d+"))";
    ga_interpolation_Lagrange_fem(md, expr, mf_vm, VM);
  }


  // --------------------------------------------------------------------
  //
  // linearized incompressibility brick  (div u = 0)
  //
  // ----------------------------------------------------------------------

  struct linear_incompressibility_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type /*ib*/,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {

      GMM_ASSERT1((matl.size() == 1 && dl.size() == 0)
                  || (matl.size() == 2 && dl.size() == 1),
                  "Wrong term and/or data number for Linear incompressibility "
                  "brick.");
      GMM_ASSERT1(mims.size() == 1, "Linear incompressibility brick need one "
                  "and only one mesh_im");
      GMM_ASSERT1(vl.size() == 2, "Wrong number of variables for linear "
                  "incompressibility brick");

      bool penalized = (dl.size() == 1);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_p = md.mesh_fem_of_variable(vl[1]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector *COEFF = 0;
      const mesh_fem *mf_data = 0;

      if (penalized) {
        COEFF = &(md.real_variable(dl[0]));
        mf_data = md.pmesh_fem_of_variable(dl[0]);
        size_type s = gmm::vect_size(*COEFF);
        if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();
        GMM_ASSERT1(s == 1, "Bad format for the penalization parameter");
      }

      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      GMM_TRACE2("Stokes term assembly");
      gmm::clear(matl[0]);
      asm_stokes_B(matl[0], mim, mf_u, mf_p, rg);

      if (penalized) {
        gmm::clear(matl[1]);
        if (mf_data) {
          asm_mass_matrix_param(matl[1], mim, mf_p, *mf_data, *COEFF, rg);
          gmm::scale(matl[1], scalar_type(-1));
        }
        else {
          asm_mass_matrix(matl[1], mim, mf_p, rg);
          gmm::scale(matl[1], -(*COEFF)[0]);
        }
      }

    }


    virtual void real_post_assembly_in_serial(const model &, size_type,
                                              const model::varnamelist &,
                                              const model::varnamelist &/*dl*/,
                                              const model::mimlist &/*mims*/,
                                              model::real_matlist &/*matl*/,
                                              model::real_veclist &,
                                              model::real_veclist &,
                                              size_type /*region*/,
                                              build_version) const
    { }


    linear_incompressibility_brick() {
      set_flags("Linear incompressibility brick",
                true /* is linear*/,
                true /* is symmetric */, false /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };

  size_type add_linear_incompressibility
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataexpr) {
#if 0
    pbrick pbr = std::make_shared<linear_incompressibility_brick>();
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    if (dataname.size()) {
      dl.push_back(dataexpr);
      tl.push_back(model::term_description(multname, multname, true));
    }
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
#else
    std::string test_varname
      = "Test_" + sup_previous_and_dot_to_varname(varname);
    std::string test_multname
      = "Test_" + sup_previous_and_dot_to_varname(multname);
    std::string expr;
    if (dataexpr.size())
      expr = "-"+multname+"*Div_"+test_varname + "-"+test_multname
        +"*Div_"+varname+"+(("+dataexpr+")*"+multname+")*"+test_multname;
    else
      expr = "-"+multname+"*Div_"+test_varname + "-"+test_multname
        +"*Div_"+varname;
    size_type ib = add_linear_term(md, mim, expr, region, true, true,
                                   "Linear incompressibility", true);
    if (ib == size_type(-1))
      ib = add_nonlinear_term
        (md, mim, expr, region, false, false,
         "Linear incompressibility (with nonlinear dependance)");
    return ib;
#endif
  }



  // ----------------------------------------------------------------------
  //
  // Mass brick
  //
  // ----------------------------------------------------------------------

  struct mass_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Mass brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Mass brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() <= 1,
                  "Wrong number of variables for mass brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      const mesh_im &mim = *mims[0];
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);

      const mesh_fem *mf_rho = 0;
      const model_real_plain_vector *rho = 0;

      if (dl.size()) {
        mf_rho = md.pmesh_fem_of_variable(dl[0]);
        rho = &(md.real_variable(dl[0]));
        size_type sl = gmm::vect_size(*rho);
        if (mf_rho) sl = sl * mf_rho->get_qdim() / mf_rho->nb_dof();
        GMM_ASSERT1(sl == 1, "Bad format of mass brick coefficient");
      }

      GMM_TRACE2("Mass matrix assembly");
      gmm::clear(matl[0]);
      if (dl.size() && mf_rho) {
        asm_mass_matrix_param(matl[0], mim, mf_u, *mf_rho, *rho, rg);
      } else {
        asm_mass_matrix(matl[0], mim, mf_u, rg);
        if (dl.size()) gmm::scale(matl[0], (*rho)[0]);
      }
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Mass brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Mass brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() <= 1,
                  "Wrong number of variables for mass brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      const mesh_im &mim = *mims[0];
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);

      const mesh_fem *mf_rho = 0;
      const model_complex_plain_vector *rho = 0;

      if (dl.size()) {
        mf_rho = md.pmesh_fem_of_variable(dl[0]);
        rho = &(md.complex_variable(dl[0]));
        size_type sl = gmm::vect_size(*rho);
        if (mf_rho) sl = sl * mf_rho->get_qdim() / mf_rho->nb_dof();
        GMM_ASSERT1(sl == 1, "Bad format of mass brick coefficient");
      }

      GMM_TRACE2("Mass matrix assembly");
      gmm::clear(matl[0]);
      if (dl.size() && mf_rho) {
        asm_mass_matrix_param(matl[0], mim, mf_u, *mf_rho, *rho, rg);
      } else {
        asm_mass_matrix(matl[0], mim, mf_u, rg);
        if (dl.size()) gmm::scale(matl[0], (*rho)[0]);
      }
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    mass_brick() {
      set_flags("Mass brick", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }

  };

  size_type add_mass_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataexpr_rho,  size_type region) {
    if (md.is_complex()) {
      pbrick pbr = std::make_shared<mass_brick>();
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      model::varnamelist dl;
      if (dataexpr_rho.size())
        dl.push_back(dataexpr_rho);
      return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
                          model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
         = "Test_" + sup_previous_and_dot_to_varname(varname);
      std::string expr;
      if (dataexpr_rho.size())
        expr ="(("+dataexpr_rho+")*"+varname+")."+test_varname;
      else
        expr = varname+"."+test_varname;
      size_type ib = add_linear_term(md, mim, expr, region, true, true,
                                     "Mass matrix", true);
      if (ib == size_type(-1))
        ib = add_nonlinear_term(md, mim, expr, region, false, false,
                                "Mass matrix (nonlinear)");
      return ib;
    }
  }

  // ----------------------------------------------------------------------
  //
  // Lumped Mass brick for first order
  //
  // ----------------------------------------------------------------------

  struct lumped_mass_for_first_order_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Lumped Mass brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Lumped Mass brick needs one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() <= 1,
                  "Wrong number of variables for lumped mass brick");

      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh &m = mf_u.linked_mesh();
      const mesh_im &mim = *mims[0];
      mesh_region rg(region);
      m.intersect_with_mpi_region(rg);

      const mesh_fem *mf_rho = 0;
      const model_real_plain_vector *rho = 0;

      if (dl.size()) {
        mf_rho = md.pmesh_fem_of_variable(dl[0]);
        rho = &(md.real_variable(dl[0]));
        size_type sl = gmm::vect_size(*rho);
        if (mf_rho) sl = sl * mf_rho->get_qdim() / mf_rho->nb_dof();
        GMM_ASSERT1(sl == 1, "Bad format of mass brick coefficient");
      }

      GMM_TRACE2("Lumped mass matrix assembly (please check that integration is 1st order.)");
      gmm::clear(matl[0]);
      if (dl.size() && mf_rho) {
        asm_lumped_mass_matrix_for_first_order_param(matl[0], mim, mf_u, *mf_rho, *rho, rg);
      } else {
        asm_lumped_mass_matrix_for_first_order(matl[0], mim, mf_u, rg);
        if (dl.size()) gmm::scale(matl[0], (*rho)[0]);
      }

    }

    lumped_mass_for_first_order_brick() {
      set_flags("Lumped mass brick", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, false /* no complex version */,
                false /* compute each time */);
    }

  };

  size_type add_lumped_mass_for_first_order_brick
  (model & md, const mesh_im &mim, const std::string &varname,
   const std::string &dataexpr_rho, size_type region) {
    pbrick pbr = std::make_shared<lumped_mass_for_first_order_brick>();
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl;
    if (dataexpr_rho.size())
      dl.push_back(dataexpr_rho);
    return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
                        model::mimlist(1, &mim), region);
  }

  // ----------------------------------------------------------------------
  //
  // From now on, DEPRECATED PART
  //
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  //
  // Generic first order time derivative brick.
  // Represents M(U^{n+1} - U^n) / dt
  //
  // ----------------------------------------------------------------------

  struct basic_d_on_dt_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Basic d/dt brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Basic d/dt brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() >= 2 && dl.size() <= 3,
                  "Wrong number of variables for basic d/dt brick");

      // It should me more convenient not to recompute the matrix if only
      // dt is modified
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
        || (md.is_var_newer_than_brick(dl[1], ib));
      if (dl.size() > 2)
        recompute_matrix = recompute_matrix ||
          md.is_var_newer_than_brick(dl[2], ib);


      if (recompute_matrix) {
        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        const mesh &m = mf_u.linked_mesh();
        const mesh_im &mim = *mims[0];
        mesh_region rg(region);
        m.intersect_with_mpi_region(rg);

        const model_real_plain_vector &dt = md.real_variable(dl[1]);
        GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for time step");

        const mesh_fem *mf_rho = 0;
        const model_real_plain_vector *rho = 0;

        if (dl.size() > 2) {
          mf_rho = md.pmesh_fem_of_variable(dl[2]);
          rho = &(md.real_variable(dl[2]));
          size_type sl = gmm::vect_size(*rho);
          if (mf_rho) sl = sl * mf_rho->get_qdim() / mf_rho->nb_dof();
          GMM_ASSERT1(sl == 1, "Bad format for density");
        }

        GMM_TRACE2("Mass matrix assembly for d_on_dt brick");
        if (dl.size() > 2 && mf_rho) {
          gmm::clear(matl[0]);
          asm_mass_matrix_param(matl[0], mim, mf_u, *mf_rho, *rho, rg);
          gmm::scale(matl[0], scalar_type(1) / dt[0]);
        } else {
          gmm::clear(matl[0]);
          asm_mass_matrix(matl[0], mim, mf_u, rg);
          if (dl.size() > 2) gmm::scale(matl[0], (*rho)[0] / dt[0]);
          else gmm::scale(matl[0], scalar_type(1) / dt[0]);
        }
      }
      gmm::mult(matl[0], md.real_variable(dl[0], 1), vecl[0]);
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type ib,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Basic d/dt brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Basic d/dt brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() >= 2 && dl.size() <= 3,
                  "Wrong number of variables for basic d/dt brick");

      // It should me more convenient not to recompute the matrix if only
      // dt is modified
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
        || (md.is_var_newer_than_brick(dl[1], ib));
      if (dl.size() > 2)
        recompute_matrix = recompute_matrix ||
          md.is_var_newer_than_brick(dl[2], ib);

      if (recompute_matrix) {
        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        const mesh &m = mf_u.linked_mesh();
        const mesh_im &mim = *mims[0];

        mesh_region rg(region);
        m.intersect_with_mpi_region(rg);

        const model_complex_plain_vector &dt = md.complex_variable(dl[1]);
        GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for time step");

        const mesh_fem *mf_rho = 0;
        const model_complex_plain_vector *rho = 0;

        if (dl.size() > 2) {
          mf_rho = md.pmesh_fem_of_variable(dl[2]);
          rho = &(md.complex_variable(dl[2]));
          size_type sl = gmm::vect_size(*rho);
          if (mf_rho) sl = sl * mf_rho->get_qdim() / mf_rho->nb_dof();
          GMM_ASSERT1(sl == 1, "Bad format for density");
        }

        GMM_TRACE2("Mass matrix assembly for d_on_dt brick");
        if (dl.size() > 2 && mf_rho) {
          gmm::clear(matl[0]);
          asm_mass_matrix_param(matl[0], mim, mf_u, *mf_rho, *rho, rg);
          gmm::scale(matl[0], scalar_type(1) / dt[0]);
        } else {
          gmm::clear(matl[0]);
          asm_mass_matrix(matl[0], mim, mf_u, rg);
          if (dl.size() > 2) gmm::scale(matl[0], (*rho)[0] / dt[0]);
          else gmm::scale(matl[0], scalar_type(1) / dt[0]);
        }
      }
      gmm::mult(matl[0], md.complex_variable(dl[0], 1), vecl[0]);
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    basic_d_on_dt_brick() {
      set_flags("Basic d/dt brick", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }

  };

  size_type add_basic_d_on_dt_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname_dt, const std::string &dataname_rho,
   size_type region) {
    pbrick pbr = std::make_shared<basic_d_on_dt_brick>();
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, varname);
    dl.push_back(dataname_dt);
    if (dataname_rho.size())
      dl.push_back(dataname_rho);
    return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
                        model::mimlist(1, &mim), region);
  }

  // ----------------------------------------------------------------------
  //
  // Generic second order time derivative brick. The velocity is considered
  // as a separate data.
  // Represents M(U^{n+1} - U^n) / (\alpha dt^2) - M V^n / (\alpha dt)
  //
  // ----------------------------------------------------------------------

  struct basic_d2_on_dt2_brick : public virtual_brick {

    mutable scalar_type old_alphadt2;

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Basic d2/dt2 brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Basic d2/dt2 brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() >= 4 && dl.size() <= 5,
                  "Wrong number of variables for basic d2/dt2 brick");

      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0);

      recompute_matrix = recompute_matrix || md.is_var_newer_than_brick(dl[2], ib);
      if (dl.size() > 4) recompute_matrix || md.is_var_newer_than_brick(dl[4], ib);

      const model_real_plain_vector &dt = md.real_variable(dl[2]);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for time step");
      const model_real_plain_vector &alpha = md.real_variable(dl[3]);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for parameter alpha");
      scalar_type alphadt2 = gmm::sqr(dt[0]) * alpha[0];

      if (!recompute_matrix && alphadt2 != old_alphadt2)
        gmm::scale(matl[0], old_alphadt2/alphadt2);
      old_alphadt2 = alphadt2;

      if (recompute_matrix) {
        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        const mesh &m = mf_u.linked_mesh();
        const mesh_im &mim = *mims[0];
        mesh_region rg(region);
        m.intersect_with_mpi_region(rg);

        const mesh_fem *mf_rho = 0;
        const model_real_plain_vector *rho = 0;

        if (dl.size() > 4) {
          mf_rho = md.pmesh_fem_of_variable(dl[4]);
          rho = &(md.real_variable(dl[4]));
          size_type sl = gmm::vect_size(*rho);
          if (mf_rho) sl = sl * mf_rho->get_qdim() / mf_rho->nb_dof();
          GMM_ASSERT1(sl == 1, "Bad format for density");
        }

        GMM_TRACE2("Mass matrix assembly for d2_on_dt2 brick");
        if (dl.size() > 4 && mf_rho) {
          gmm::clear(matl[0]);
          asm_mass_matrix_param(matl[0], mim, mf_u, *mf_rho, *rho, rg);
          gmm::scale(matl[0], scalar_type(1) / alphadt2);
        } else {
          gmm::clear(matl[0]);
          asm_mass_matrix(matl[0], mim, mf_u, rg);
          if (dl.size() > 4) gmm::scale(matl[0], (*rho)[0] / alphadt2);
          else gmm::scale(matl[0], scalar_type(1) / alphadt2);
        }
      }
      gmm::mult(matl[0], md.real_variable(dl[0], 1), vecl[0]);
      gmm::mult_add(matl[0], gmm::scaled(md.real_variable(dl[1], 1), dt[0]),
                    vecl[0]);
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type ib,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version version) const {
      GMM_ASSERT1(matl.size() == 1,
                  "Basic d2/dt2 brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "Basic d2/dt2 brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() >= 4 && dl.size() <= 5,
                  "Wrong number of variables for basic d2/dt2 brick");

      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0);

      recompute_matrix = recompute_matrix || md.is_var_newer_than_brick(dl[2], ib);
      if (dl.size() > 4) recompute_matrix || md.is_var_newer_than_brick(dl[4], ib);


      const model_complex_plain_vector &dt = md.complex_variable(dl[2]);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for time step");
      const model_complex_plain_vector &alpha = md.complex_variable(dl[3]);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for parameter alpha");
      scalar_type alphadt2 = gmm::real(gmm::sqr(dt[0]) * alpha[0]);

      if (!recompute_matrix && alphadt2 != old_alphadt2)
        gmm::scale(matl[0], old_alphadt2/alphadt2);
      old_alphadt2 = alphadt2;

      if (recompute_matrix) {
        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        const mesh &m = mf_u.linked_mesh();
        const mesh_im &mim = *mims[0];
        mesh_region rg(region);
        m.intersect_with_mpi_region(rg);

        const mesh_fem *mf_rho = 0;
        const model_complex_plain_vector *rho = 0;

        if (dl.size() > 4) {
          mf_rho = md.pmesh_fem_of_variable(dl[4]);
          rho = &(md.complex_variable(dl[4]));
          size_type sl = gmm::vect_size(*rho);
          if (mf_rho) sl = sl * mf_rho->get_qdim() / mf_rho->nb_dof();
          GMM_ASSERT1(sl == 1, "Bad format for density");
        }

        GMM_TRACE2("Mass matrix assembly for d2_on_dt2 brick");
        if (dl.size() > 4 && mf_rho) {
          gmm::clear(matl[0]);
          asm_mass_matrix_param(matl[0], mim, mf_u, *mf_rho, *rho, rg);
          gmm::scale(matl[0], scalar_type(1) / alphadt2);
        } else {
          gmm::clear(matl[0]);
          asm_mass_matrix(matl[0], mim, mf_u, rg);
          if (dl.size() > 4) gmm::scale(matl[0], (*rho)[0] / alphadt2);
          else gmm::scale(matl[0], scalar_type(1) / alphadt2);
        }
      }
      gmm::mult(matl[0], md.complex_variable(dl[0], 1), vecl[0]);
      gmm::mult_add(matl[0], gmm::scaled(md.complex_variable(dl[1], 1), dt[0]),
                    vecl[0]);
    }

    virtual std::string declare_volume_assembly_string
    (const model &, size_type, const model::varnamelist &,
     const model::varnamelist &) const {
      return std::string();
    }

    basic_d2_on_dt2_brick() {
      set_flags("Basic d2/dt2 brick", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */);
    }

  };

  size_type add_basic_d2_on_dt2_brick
  (model &md, const mesh_im &mim, const std::string &varnameU,
   const std::string &datanameV,
   const std::string &dataname_dt,
   const std::string &dataname_alpha,
   const std::string &dataname_rho,
   size_type region) {
    pbrick pbr = std::make_shared<basic_d2_on_dt2_brick>();
    model::termlist tl;
    tl.push_back(model::term_description(varnameU, varnameU, true));
    model::varnamelist dl(1, varnameU);
    dl.push_back(datanameV);
    dl.push_back(dataname_dt);
    dl.push_back(dataname_alpha);
    if (dataname_rho.size())
      dl.push_back(dataname_rho);
    return md.add_brick(pbr, model::varnamelist(1, varnameU), dl, tl,
                        model::mimlist(1, &mim), region);
  }



  // ----------------------------------------------------------------------
  //
  //
  // Standard time dispatchers
  //
  //
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  //
  // theta-method dispatcher
  //
  // ----------------------------------------------------------------------

    void theta_method_dispatcher::set_dispatch_coeff(const model &md, size_type ib) const {
      scalar_type theta;
      if (md.is_complex())
        theta = gmm::real(md.complex_variable(param_names[0])[0]);
      else
        theta = md.real_variable(param_names[0])[0];
      // coefficient for the matrix term
      md.matrix_coeff_of_brick(ib) = theta;
      // coefficient for the standard rhs
      md.rhs_coeffs_of_brick(ib)[0] = theta;
      // coefficient for the additional rhs
      md.rhs_coeffs_of_brick(ib)[1] = (scalar_type(1) - theta);
    }


    void theta_method_dispatcher::next_real_iter
    (const model &md, size_type ib, const model::varnamelist &vl,
     const model::varnamelist &dl, model::real_matlist &matl,
     std::vector<model::real_veclist> &vectl,
     std::vector<model::real_veclist> &vectl_sym, bool first_iter) const {
      next_iter(md, ib, vl, dl, matl, vectl, vectl_sym, first_iter);
    }

    void theta_method_dispatcher::next_complex_iter
    (const model &md, size_type ib, const model::varnamelist &vl,
     const model::varnamelist &dl,
     model::complex_matlist &matl,
     std::vector<model::complex_veclist> &vectl,
     std::vector<model::complex_veclist> &vectl_sym,
     bool first_iter) const {
      next_iter(md, ib, vl, dl, matl, vectl, vectl_sym, first_iter);
    }

    void theta_method_dispatcher::asm_real_tangent_terms
    (const model &md, size_type ib, model::real_matlist &/* matl */,
     std::vector<model::real_veclist> &/* vectl */,
     std::vector<model::real_veclist> &/* vectl_sym */,
     build_version version) const
    { md.brick_call(ib, version, 0); }

    void theta_method_dispatcher::asm_complex_tangent_terms
    (const model &md, size_type ib, model::complex_matlist &/* matl */,
     std::vector<model::complex_veclist> &/* vectl */,
     std::vector<model::complex_veclist> &/* vectl_sym */,
     build_version version) const
    { md.brick_call(ib, version, 0); }

    theta_method_dispatcher::theta_method_dispatcher(const std::string &THETA)
      : virtual_dispatcher(2) {
      param_names.push_back(THETA);
    }

  void add_theta_method_dispatcher
  (model &md, dal::bit_vector ibricks, const std::string &THETA) {
    pdispatcher pdispatch = std::make_shared<theta_method_dispatcher>(THETA);
    for (dal::bv_visitor i(ibricks); !i.finished(); ++i)
      md.add_time_dispatcher(i, pdispatch);
  }

  void velocity_update_for_order_two_theta_method
  (model &md, const std::string &U, const std::string &V,
   const std::string &pdt, const std::string &ptheta) {

    // V^{n+1} = (1-1/theta)*V^n + (1/theta)*(U^{n+1} - U^n)/dt

    if (md.is_complex()) {
      const model_complex_plain_vector &dt = md.complex_variable(pdt);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for time step");
      const model_complex_plain_vector &theta = md.complex_variable(ptheta);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for parameter theta");

      gmm::copy(gmm::scaled(md.complex_variable(V, 1),
                            scalar_type(1) - scalar_type(1) / theta[0]),
                md.set_complex_variable(V, 0));
      gmm::add(gmm::scaled(md.complex_variable(U, 0),
                           scalar_type(1) / (theta[0]*dt[0])),
               md.set_complex_variable(V, 0));
      gmm::add(gmm::scaled(md.complex_variable(U, 1),
                           -scalar_type(1) / (theta[0]*dt[0])),
               md.set_complex_variable(V, 0));
    } else {
      const model_real_plain_vector &dt = md.real_variable(pdt);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for time step");
      const model_real_plain_vector &theta = md.real_variable(ptheta);
      GMM_ASSERT1(gmm::vect_size(dt) == 1, "Bad format for parameter theta");

      gmm::copy(gmm::scaled(md.real_variable(V, 1),
                            scalar_type(1) - scalar_type(1) / theta[0]),
                md.set_real_variable(V, 0));
      gmm::add(gmm::scaled(md.real_variable(U, 0),
                           scalar_type(1) / (theta[0]*dt[0])),
               md.set_real_variable(V, 0));
      gmm::add(gmm::scaled(md.real_variable(U, 1),
                           -scalar_type(1) / (theta[0]*dt[0])),
               md.set_real_variable(V, 0));
    }
  }

  // ----------------------------------------------------------------------
  //
  // Newmark scheme dispatcher
  //
  // ----------------------------------------------------------------------

  void velocity_update_for_Newmark_scheme
  (model &md, size_type id2dt2b, const std::string &U, const std::string &V,
   const std::string &pdt, const std::string &ptwobeta,
   const std::string &pgamma) {

    md.disable_brick(id2dt2b);

    if (md.is_complex()) {
      complex_type twobeta = md.complex_variable(ptwobeta)[0];
      complex_type gamma = md.complex_variable(pgamma)[0];
      complex_type dt = md.complex_variable(pdt)[0];

      // Modification of the parameter for the theta-method.
      if (twobeta != gamma) {
        md.set_complex_variable(ptwobeta)[0] = gamma;
        md.set_dispatch_coeff();  // valid the change of coefficients.
      }

      // Computation of the residual (including the linear parts).
      md.assembly(model::BUILD_RHS_WITH_LIN);

      size_type nbdof = gmm::vect_size(md.complex_variable(U));
      model_complex_plain_vector W(nbdof), RHS(nbdof);
      gmm::copy(gmm::sub_vector(md.complex_rhs(), md.interval_of_variable(U)),
                RHS);

      // Compute the velocity. Inversion with CG.
      gmm::iteration iter(1e-12, 0, 100000);
      gmm::cg(md.linear_complex_matrix_term(id2dt2b, 0),
              W, RHS, gmm::identity_matrix(), iter);
      GMM_ASSERT1(iter.converged(), "Velocity not well computed");
      gmm::add(md.complex_variable(V, 1),
               gmm::scaled(W, complex_type(1)/(twobeta*dt)),
               md.set_complex_variable(V, 0));

      // Cancel the modification of the parameter for the theta-method.
      if (twobeta != gamma) {
        md.set_complex_variable(ptwobeta)[0] = twobeta;
        md.set_dispatch_coeff();  // valid the change of coefficients.
      }


      GMM_ASSERT1(false, "to be done");
    } else {
      scalar_type twobeta = md.real_variable(ptwobeta)[0];
      scalar_type gamma = md.real_variable(pgamma)[0];
      scalar_type dt = md.real_variable(pdt)[0];



      // Modification of the parameter for the theta-method.
      if (twobeta != gamma) {
        md.set_real_variable(ptwobeta)[0] = gamma;
        md.set_dispatch_coeff();  // valid the change of coefficients.
      }

      // Computation of the residual (including the linear parts).
      md.assembly(model::BUILD_RHS_WITH_LIN);

      size_type nbdof = gmm::vect_size(md.real_variable(U));
      model_real_plain_vector W(nbdof), RHS(nbdof);
      gmm::copy(gmm::sub_vector(md.real_rhs(), md.interval_of_variable(U)),
                RHS);

      // Compute the velocity. Inversion with CG.
      gmm::iteration iter(1e-12, 0, 100000);
      gmm::cg(md.linear_real_matrix_term(id2dt2b, 0),
              W, RHS, gmm::identity_matrix(), iter);
      GMM_ASSERT1(iter.converged(), "Velocity not well computed");
      gmm::add(md.real_variable(V, 1),
               gmm::scaled(W, scalar_type(1)/(twobeta*dt)),
               md.set_real_variable(V, 0));

      // Cancel the modification of the parameter for the theta-method.
      if (twobeta != gamma) {
        md.set_real_variable(ptwobeta)[0] = twobeta;
        md.set_dispatch_coeff();  // valid the change of coefficients.
      }

    }
    md.enable_brick(id2dt2b);
  }


  // ----------------------------------------------------------------------
  //
  // midpoint dispatcher
  //
  // ----------------------------------------------------------------------


  class midpoint_dispatcher : public virtual_dispatcher {

    gmm::uint64_type id_num;

  public :

    typedef model::build_version build_version;

    void set_dispatch_coeff(const model &md, size_type ib) const {
      md.matrix_coeff_of_brick(ib) = scalar_type(1)/scalar_type(2);
      md.rhs_coeffs_of_brick(ib)[0] = scalar_type(1);
      md.rhs_coeffs_of_brick(ib)[1] = scalar_type(1)/scalar_type(2);
    }

    template <typename MATLIST, typename VECTLIST>
    inline void next_iter(const model &md, size_type ib,
                          const model::varnamelist &vl,
                          const model::varnamelist &dl,
                          MATLIST &/* matl */,
                          VECTLIST &vectl, VECTLIST &vectl_sym,
                          bool first_iter) const {

      pbrick pbr = md.brick_pointer(ib);

      if (first_iter) { // For the moment, temporaries are deleted by
        // model::first_iter before the call to virtual_dispatcher::next_iter
        if (!(pbr->is_linear()))
          md.add_temporaries(vl, id_num); // add temporaries for all variables
        md.add_temporaries(dl, id_num); // add temporaries for versionned data
        for (auto &&v : vectl[1]) gmm::clear(v);
        for (auto &&v : vectl_sym[1]) gmm::clear(v);
      }

      if (pbr->is_linear()) { // If the problem is linear, add the term
        // coming from the previous iteration as a second rhs.
        // This rhs is only used for this.
        if (first_iter) md.update_brick(ib, model::BUILD_RHS);
        for (auto &&v : vectl[1]) gmm::clear(v);
        for (auto &&v : vectl_sym[1]) gmm::clear(v);
        md.linear_brick_add_to_rhs(ib, 1, 0);
      }
    }

    void next_real_iter
    (const model &md, size_type ib, const model::varnamelist &vl,
     const model::varnamelist &dl, model::real_matlist &matl,
     std::vector<model::real_veclist> &vectl,
     std::vector<model::real_veclist> &vectl_sym, bool first_iter) const {
      next_iter(md, ib, vl, dl, matl, vectl, vectl_sym, first_iter);
    }

    void next_complex_iter
    (const model &md, size_type ib, const model::varnamelist &vl,
     const model::varnamelist &dl,
     model::complex_matlist &matl,
     std::vector<model::complex_veclist> &vectl,
     std::vector<model::complex_veclist> &vectl_sym,
     bool first_iter) const {
      next_iter(md, ib, vl, dl, matl, vectl, vectl_sym, first_iter);
    }

    void asm_real_tangent_terms
    (const model &md, size_type ib, model::real_matlist &/* matl */,
     std::vector<model::real_veclist> &vectl,
     std::vector<model::real_veclist> &vectl_sym,
     build_version version) const {

      scalar_type half = scalar_type(1)/scalar_type(2);
      pbrick pbr = md.brick_pointer(ib);
      size_type ind;

      const model::varnamelist &vl = md.varnamelist_of_brick(ib);
      const model::varnamelist &dl = md.datanamelist_of_brick(ib);

      if (!(pbr->is_linear())) { // compute the mean variables
        for (size_type i = 0; i < vl.size(); ++i) {
          bool is_uptodate = md.temporary_uptodate(vl[i], id_num, ind);
          if (!is_uptodate && ind != size_type(-1))
            gmm::add(gmm::scaled(md.real_variable(vl[i], 0), half),
                     gmm::scaled(md.real_variable(vl[i], 1), half),
                     md.set_real_variable(vl[i], ind));
          md.set_default_iter_of_variable(vl[i], ind);
        }
      }

      // compute the mean data
      for (size_type i = 0; i < dl.size(); ++i) {
        bool is_uptodate = md.temporary_uptodate(dl[i], id_num, ind);
        if (!is_uptodate && ind != size_type(-1)) {
          gmm::add(gmm::scaled(md.real_variable(dl[i], 0), half),
                   gmm::scaled(md.real_variable(dl[i], 1), half),
                   md.set_real_variable(dl[i], ind));
        }
        md.set_default_iter_of_variable(dl[i], ind);
      }

      // call the brick for the mid-time step.
      md.brick_call(ib, version, 0);
      if (pbr->is_linear()) { // update second rhs (is updated by next_iter
        // but the call to the brick may have changed the matrices.
        for (auto &&v : vectl[1]) gmm::clear(v);
        for (auto &&v : vectl_sym[1]) gmm::clear(v);
        md.linear_brick_add_to_rhs(ib, 1, 1);
      }

      md.reset_default_iter_of_variables(dl);
      if (!(pbr->is_linear()))
        md.reset_default_iter_of_variables(vl);
    }

    virtual void asm_complex_tangent_terms
    (const model &md, size_type ib, model::complex_matlist &/* matl */,
     std::vector<model::complex_veclist> &vectl,
     std::vector<model::complex_veclist> &vectl_sym,
     build_version version) const {

      scalar_type half = scalar_type(1)/scalar_type(2);
      pbrick pbr = md.brick_pointer(ib);
      size_type ind;

      const model::varnamelist &vl = md.varnamelist_of_brick(ib);
      const model::varnamelist &dl = md.datanamelist_of_brick(ib);

      if (!(pbr->is_linear())) { // compute the mean variables
        for (size_type i = 0; i < vl.size(); ++i) {
          bool is_uptodate = md.temporary_uptodate(vl[i], id_num, ind);
          if (!is_uptodate && ind != size_type(-1))
            gmm::add(gmm::scaled(md.complex_variable(vl[i], 0), half),
                     gmm::scaled(md.complex_variable(vl[i], 1), half),
                     md.set_complex_variable(vl[i], ind));
          md.set_default_iter_of_variable(vl[i], ind);
        }
      }

      // compute the mean data
      for (size_type i = 0; i < dl.size(); ++i) {
        bool is_uptodate = md.temporary_uptodate(dl[i], id_num, ind);
        if (!is_uptodate && ind != size_type(-1)) {
          gmm::add(gmm::scaled(md.complex_variable(dl[i], 0), half),
                   gmm::scaled(md.complex_variable(dl[i], 1), half),
                   md.set_complex_variable(dl[i], ind));
        }
        md.set_default_iter_of_variable(dl[i], ind);
      }

      // call the brick for the mid-time step.
      md.brick_call(ib, version, 0);
      if (pbr->is_linear()) { // update second rhs (is updated by next_iter
        // but the call to the brick may have changed the matrices.
        for (auto &&v : vectl[1]) gmm::clear(v);
        for (auto &&v : vectl_sym[1]) gmm::clear(v);
        md.linear_brick_add_to_rhs(ib, 1, 1);
      }

      md.reset_default_iter_of_variables(dl);
      if (!(pbr->is_linear()))
        md.reset_default_iter_of_variables(vl);
    }

    midpoint_dispatcher() : virtual_dispatcher(2)
    { id_num = act_counter(); }

  };

  void add_midpoint_dispatcher(model &md, dal::bit_vector ibricks) {
    pdispatcher pdispatch = std::make_shared<midpoint_dispatcher>();
    for (dal::bv_visitor i(ibricks); !i.finished(); ++i)
      md.add_time_dispatcher(i, pdispatch);
  }


}  /* end of namespace getfem.                                             */

