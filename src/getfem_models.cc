/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================

 Copyright (C) 2009-2015 Yves Renard

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

#include <iomanip>
#include "gmm/gmm_range_basis.h"
#include "gmm/gmm_solver_cg.h"
#include "gmm/gmm_condition_number.h"
#include "getfem/getfem_models.h"
#include "getfem/getfem_assembling.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_generic_assembly.h"


namespace getfem {

  void model::var_description::set_size(size_type s) {
    n_temp_iter = 0;
    default_iter = 0;
    if (is_complex)
      complex_value.resize(n_iter);
    else
      real_value.resize(n_iter);
    v_num_var_iter.resize(n_iter);
    v_num_iter.resize(n_iter);
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

  void model::var_description::clear_temporaries(void) {
    n_temp_iter = 0;
    default_iter = 0;
    if (is_complex)
      complex_value.resize(n_iter);
    else
      real_value.resize(n_iter);
  }

  bool model::check_name_validity(const std::string &name, bool assert) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) {
      GMM_ASSERT1(!assert, "Variable " << name << " already exists");
      return false;
    }

    if (variable_groups.find(name) != variable_groups.end()) {
      GMM_ASSERT1(!assert,
                  name << " corresponds to an already existing group of "
                  "variables name");
      return false;
    }

    if (macros.find(name) != macros.end()) {
      GMM_ASSERT1(!assert,
                  name << " corresponds to an already existing macro");
      return false;
    }

    if (name.compare("X") == 0) {
      GMM_ASSERT1(!assert, "X is a reserved keyword of the generic "
                  "assembly language");
      return false;
    }

    int ga_valid = ga_check_name_validity(name);
    if (ga_valid == 1) {
      GMM_ASSERT1(!assert, "Invalid variable name, corresponds to an "
                "operator or function name of the generic assembly language");
      return false;
    }

    if (ga_valid == 2) {
      GMM_ASSERT1(!assert, "Invalid variable name having a reserved "
                  "prefix used by the generic assembly language");
      return false;
    }

    if (ga_valid == 3) {
      std::string org_name = sup_previous_and_dot_to_varname(name);
      if (org_name.size() < name.size() &&
          variables.find(org_name) != variables.end()) {
        GMM_ASSERT1(!assert,
                    "Dot and Previous are reserved prefix used for time "
                    "integration schemes");
        return false;
      }
    }

    bool valid = true;
    if (name.size() == 0) valid = false;
    else {
      if (!isalpha(name[0])) valid = false;
      for (size_type i = 1; i < name.size(); ++i)
        if (!(isalnum(name[i]) || name[i] == '_')) valid = false;
    }
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

  void model::resize_global_system() const {
    size_type tot_size = 0;

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
         ++it) {
      if (it->second.is_variable && it->second.is_disabled)
        it->second.I  = gmm::sub_interval(0,0);
      if (it->second.is_variable && !(it->second.is_affine_dependent)
          && !(it->second.is_disabled)) {
        it->second.I = gmm::sub_interval(tot_size, it->second.size());
        tot_size += it->second.size();
      }
    }

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
         ++it) {
      if (it->second.is_affine_dependent) {
        VAR_SET::iterator it2 = variables.find(it->second.org_name);
        it->second.I = it2->second.I;
        it->second.set_size(it2->second.size());
      }
    }

    if (complex_version) {
      gmm::resize(cTM, tot_size, tot_size);
      gmm::resize(crhs, tot_size);
    }
    else {
      gmm::resize(rTM, tot_size, tot_size);
      gmm::resize(rrhs, tot_size);
    }

    for (dal::bv_visitor ib(valid_bricks); !ib.finished(); ++ib) {
      const brick_description &brick = bricks[ib];
      for (size_type j = 0; j < brick.tlist.size(); ++j) {
        if (brick.tlist[j].is_global) {
          brick.terms_to_be_computed = true;
          break;
        }
      }
    }
  }

  void model::actualize_sizes(void) const {
    // cout << "begin act size" << endl;
    bool actualized = false;
    getfem::local_guard lock = locks_.get_lock();
    if (actualized) return; // If multiple threads are calling the method

    act_size_to_be_done = false;

    std::map<std::string, std::vector<std::string> > multipliers;
    std::map<std::string, bool > tobedone;

//     #if GETFEM_PARA_LEVEL > 1
//     int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk);
//     double t_ref = MPI_Wtime();
//     cout << "Actualize size called from thread " << rk << endl;
//     #endif


    // In case of change in fems or mims, linear terms have to be recomputed
    // We couls select which brick is to be recomputed if we would be able
    // to know which fem or mim is changed.
    for (dal::bv_visitor ib(valid_bricks); !ib.finished(); ++ib)
      bricks[ib].terms_to_be_computed = true;

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
         ++it) {
      if (it->second.is_fem_dofs && !(it->second.is_affine_dependent)) {
        if ((it->second.filter & VDESCRFILTER_CTERM)
            || (it->second.filter & VDESCRFILTER_INFSUP)) {
          VAR_SET::iterator it2 = variables.find(it->second.filter_var);
          GMM_ASSERT1(it2 != variables.end(), "The primal variable of the "
                      "multiplier does not exist");
          GMM_ASSERT1(it2->second.is_fem_dofs, "The primal variable of the "
                      "multiplier is not a fem variable");
          multipliers[it->second.filter_var].push_back(it->first);
          if (it->second.v_num < it->second.mf->version_number() ||
              it->second.v_num < it2->second.mf->version_number()) {
            tobedone[it->second.filter_var] = true;
          }
        }
        switch (it->second.filter) {
        case VDESCRFILTER_NO:
          if (it->second.v_num < it->second.mf->version_number()) {
            size_type s = it->second.mf->nb_dof();
            if (!it->second.is_variable) s *= it->second.qdim;
            it->second.set_size(s);
            it->second.v_num = act_counter();
          }
          break;
        case VDESCRFILTER_REGION:
          if (it->second.v_num < it->second.mf->version_number()) {
            dal::bit_vector dor
              = it->second.mf->dof_on_region(it->second.m_region);
            it->second.partial_mf->adapt(dor);
            it->second.set_size(it->second.partial_mf->nb_dof());
            it->second.v_num = act_counter();
          }
          break;
        default : break;
        }
      }
    
      if (it->second.pim_data != 0
          && it->second.v_num < it->second.pim_data->version_number()) {
        const im_data *pimd = it->second.pim_data;
        size_type data_size = pimd->nb_filtered_index()*pimd->nb_tensor_elem();
        it->second.set_size(data_size);
        it->second.v_num = act_counter();
      }
    }

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
         ++it) {
      if (it->second.is_fem_dofs && !(it->second.is_affine_dependent) &&
          ((it->second.filter & VDESCRFILTER_CTERM)
           || (it->second.filter & VDESCRFILTER_INFSUP))) {
        if (tobedone.find(it->second.filter_var) != tobedone.end()) {
          // This step forces the recomputation of corresponding bricks.
          // A test to check if a modification is really necessary could
          // be done first ... (difficult to coordinate with other
          // multipliers)
          dal::bit_vector alldof; alldof.add(0, it->second.mf->nb_dof());
          it->second.partial_mf->adapt(alldof);
          it->second.set_size(it->second.partial_mf->nb_dof());
          it->second.v_num = act_counter();
        }
      }
    }

    resize_global_system();

    for (std::map<std::string, bool>::iterator itbd = tobedone.begin();
         itbd != tobedone.end(); ++itbd) {
//       #if GETFEM_PARA_LEVEL > 1
//       double tt_ref = MPI_Wtime();
//       if (!rk) cout << "compute size of multipliers for " << itbd->first
//                     << endl;
//       #endif

      std::vector<std::string> &mults = multipliers[itbd->first];
      VAR_SET::iterator it2 = variables.find(itbd->first);

      gmm::col_matrix< gmm::rsvector<scalar_type> > MGLOB;
      if (mults.size() > 1) {
        size_type s = 0;
        for (size_type k = 0; k < mults.size(); ++k) {
          VAR_SET::iterator it = variables.find(mults[k]);
          s += it->second.mf->nb_dof();
        }
        gmm::resize(MGLOB, it2->second.mf->nb_dof(), s);
      }
      size_type s = 0;
      std::set<size_type> glob_columns;
      // std::vector<dal::bit_vector> mult_kept_dofs;
      for (size_type k = 0; k < mults.size(); ++k) {
        VAR_SET::iterator it = variables.find(mults[k]);

        // Obtaining the coupling matrix between the multipier and
        // the primal variable. A search is done on all the terms of the
        // model. Only the the corresponding linear terms are added.
        // If no linear term is available, a mass matrix is used.
        gmm::col_matrix< gmm::rsvector<scalar_type> >
          MM(it2->second.associated_mf().nb_dof(), it->second.mf->nb_dof());
        bool termadded = false;

        if (it->second.filter & VDESCRFILTER_CTERM) {

          for (dal::bv_visitor ib(valid_bricks); !ib.finished(); ++ib) {
            const brick_description &brick = bricks[ib];
            bool bupd = false;
            bool cplx = is_complex() && brick.pbr->is_complex();

            for (size_type j = 0; j < brick.tlist.size(); ++j) {

              const term_description &term = brick.tlist[j];

              if (term.is_matrix_term) {
                if (term.is_global) {
                  bool varc = false, multc = false;
                  for (size_type iv = 0; iv < brick.vlist.size(); ++iv) {
                    if (!(mults[k].compare(brick.vlist[iv]))) multc = true;
                    if (!(it2->first.compare(brick.vlist[iv]))) varc = true;
                  }
                  if (multc && varc) {
                    GMM_ASSERT1(!cplx, "Sorry, not taken into account");
                    generic_expressions.clear();
                    if (!bupd) {
                      brick.terms_to_be_computed = true;
                      update_brick(ib, BUILD_MATRIX);
                      bupd = true;
                    }
                    if (generic_expressions.size()) {
                      GMM_TRACE2("Generic assembly for actualize sizes");
                      ga_workspace workspace(*this);
                      for (std::list<gen_expr>::iterator ig
                             = generic_expressions.begin();
                           ig != generic_expressions.end(); ++ig) {
                        workspace.add_expression(ig->expr,ig->mim,ig->region);
                      }
                      gmm::clear(rTM);
                      workspace.set_assembled_matrix(rTM);
                      workspace.assembly(2);
                      gmm::add
                        (gmm::sub_matrix(rTM, it->second.I, it2->second.I),MM);
                      gmm::add(gmm::transposed
                               (gmm::sub_matrix(rTM, it2->second.I,
                                                it->second.I)), MM);
                      bupd = false;
                    } else {
                      gmm::add(gmm::sub_matrix(brick.rmatlist[j],
                                               it->second.I, it2->second.I),
                               MM);
                      gmm::add(gmm::transposed(gmm::sub_matrix
                                               (brick.rmatlist[j],
                                                it2->second.I, it->second.I)),
                               MM);
                    }
                  }
                } else if (!mults[k].compare(term.var1) && 
                    !it2->first.compare(term.var2)) {
                  if (!bupd) {
                    brick.terms_to_be_computed = true;
                    update_brick(ib, BUILD_MATRIX);
                    bupd = true;
                  }
                  if (cplx)
                    gmm::add
                      (gmm::transposed(gmm::real_part(brick.cmatlist[j])), MM);
                  else
                    gmm::add(gmm::transposed(brick.rmatlist[j]), MM);
                  termadded = true;
                  
                } else if (!mults[k].compare(term.var2) &&
                           !it2->first.compare(term.var1)) {
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
            GMM_WARNING1("No term found to filter multiplier " << it->first
                         << ". Variable is cancelled");
        } else if (it->second.filter & VDESCRFILTER_INFSUP) {
          mesh_region rg(it->second.m_region);
          it->second.mim->linked_mesh().intersect_with_mpi_region(rg);
          asm_mass_matrix(MM, *(it->second.mim), it2->second.associated_mf(),
                          *(it->second.mf), rg);
        }

        MPI_SUM_SPARSE_MATRIX(MM);

        //
        // filtering
        //
        std::set<size_type> columns;
        gmm::range_basis(MM, columns);

        if (mults.size() > 1) {
          gmm::copy(MM, gmm::sub_matrix
                    (MGLOB, gmm::sub_interval(0,
                                        it2->second.associated_mf().nb_dof()),
                     gmm::sub_interval(s, it->second.mf->nb_dof())));
          for (std::set<size_type>::iterator itt = columns.begin();
             itt != columns.end(); ++itt)
            glob_columns.insert(s + *itt);
          s += it->second.mf->nb_dof();
        } else {
          dal::bit_vector kept;
          for (std::set<size_type>::iterator itt = columns.begin();
               itt != columns.end(); ++itt)
            kept.add(*itt);
          if (it->second.filter & VDESCRFILTER_REGION)
            kept &= it->second.mf->dof_on_region(it->second.m_region);
          // kept &= mult_kept_dofs[k];
          it->second.partial_mf->adapt(kept);
          it->second.set_size(it->second.partial_mf->nb_dof());
          it->second.v_num = act_counter();
        }
      }

//         #if GETFEM_PARA_LEVEL > 1
//         if (!rk) cout << "Range basis for  multipliers for " << itbd->first << " time : " << MPI_Wtime()-tt_ref << endl;

//         #endif

      if (mults.size() > 1) {
        gmm::range_basis(MGLOB, glob_columns, 1E-12, gmm::col_major(), true);


//         #if GETFEM_PARA_LEVEL > 1
//         if (!rk) cout << "Producing partial mf for  multipliers for " << itbd->first << " time : " << MPI_Wtime()-tt_ref << endl;

//         #endif

        s = 0;
        for (size_type k = 0; k < mults.size(); ++k) {
          VAR_SET::iterator it = variables.find(mults[k]);
          dal::bit_vector kept;
          size_type nbdof = it->second.mf->nb_dof();
          for (std::set<size_type>::iterator itt = glob_columns.begin();
               itt != glob_columns.end(); ++itt)
            if (*itt >= s && *itt < s + nbdof) kept.add(*itt-s);
          if (it->second.filter & VDESCRFILTER_REGION)
            kept &= it->second.mf->dof_on_region(it->second.m_region);
          // kept &= mult_kept_dofs[k];
          it->second.partial_mf->adapt(kept);
          it->second.set_size(it->second.partial_mf->nb_dof());
          it->second.v_num = act_counter();
          s += it->second.mf->nb_dof();
        }
      }
//       #if GETFEM_PARA_LEVEL > 1
//       if (!rk) cout << "End compute size of  multipliers for " << itbd->first << " time : " << MPI_Wtime()-tt_ref << endl;

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
      for (VAR_SET::const_iterator it = variables.begin();
           it != variables.end(); ++it) {
        if (it->second.is_variable) ost << "Variable ";
        else ost << "Data     ";
        ost << std::setw(20) << std::left << it->first;
        if (it->second.n_iter == 1)
          ost << " 1 copy   ";
        else
          ost << std::setw(2) << std::right << it->second.n_iter
              << " copies ";
        if (it->second.is_fem_dofs) ost << "fem dependant ";
        else ost << "constant size ";
        size_type si = it->second.size();
        ost << std::setw(8) << std::right << si;
        if (is_complex()) ost << " complex";
        ost << " double" << ((si > 1) ? "s." : ".");
        if (it->second.is_variable &&
            is_disabled_variable(it->first)) ost << "\t (disabled)";
        else                                 ost << "\t           ";
        if (it->second.pim_data != 0) ost << "\t (is im_data)";
        if (it->second.is_affine_dependent) ost << "\t (is affine dependent)";
        ost << endl;
      }
    }
  }

  void model::listresiduals(std::ostream &ost) const {
    if (variables.size() == 0)
      ost << "Model with no variable nor data" << endl;
    else {
      bool firstvar(true);
      for (VAR_SET::const_iterator it = variables.begin();
           it != variables.end(); ++it) {
        if (it->second.is_variable) {
          const gmm::sub_interval &II = interval_of_variable(it->first);
          scalar_type res = gmm::vect_norm2(gmm::sub_vector(rrhs, II));
          if (!firstvar) cout << ", ";
          ost << "res_" << it->first << "= " << std::setw(11) << res;
          firstvar = false;
        }
      }
      ost << endl;
    }
  }

  void model::add_fixed_size_variable(const std::string &name, size_type size,
                                      size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(true, is_complex(), false, niter);
    act_size_to_be_done = true;
    variables[name].set_size(size);
  }

  void model::resize_fixed_size_variable(const std::string &name,
                                         size_type size) {
    GMM_ASSERT1(!(variables[name].is_fem_dofs),
                "Cannot explicitly resize a fem variable or data");
    GMM_ASSERT1(variables[name].pim_data == 0,
                "Cannot explicitly resize an im data");
    variables[name].set_size(size);
  }

  void resize_fixed_size_variable(const std::string &name, size_type size);


  void model::add_fixed_size_data(const std::string &name, size_type size,
                                  size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(false, is_complex(), false, niter);
    variables[name].set_size(size);
  }

  void model::add_im_data(const std::string &name, const im_data &im_data,
                          size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(false, is_complex(), false, niter);
    variables[name].pim_data = &im_data;
    variables[name].set_size(im_data.nb_filtered_index() * im_data.nb_tensor_elem());
    add_dependency(im_data);
  }

  void model::add_fem_variable(const std::string &name, const mesh_fem &mf,
                               size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_NO, &mf);
    variables[name].set_size(mf.nb_dof());
    add_dependency(mf);
    act_size_to_be_done = true;
    leading_dim = std::max(leading_dim, mf.linked_mesh().dim());
  }

  void model::add_filtered_fem_variable(const std::string &name,
                                        const mesh_fem &mf,
                                        size_type region, size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_REGION, &mf, region);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_affine_dependent_variable(const std::string &name,
                                            const std::string &org_name,
                                            scalar_type alpha) {
    check_name_validity(name);
    VAR_SET::const_iterator it = variables.find(org_name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << org_name);
    GMM_ASSERT1(it->second.is_variable && !(it->second.is_affine_dependent),
                "The original variable should be a variable");
    variables[name] = variables[org_name];
    variables[name].is_affine_dependent = true;
    variables[name].org_name = org_name;
    variables[name].alpha = alpha;
    variables[name].set_size(variables[org_name].size());
  }

  void model::add_fem_data(const std::string &name, const mesh_fem &mf,
                               dim_type qdim, size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(false, is_complex(), true, niter,
                                      VDESCRFILTER_NO, &mf, size_type(-1),
                                      qdim);
    variables[name].set_size(mf.nb_dof()*qdim);
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             const std::string &primal_name,
                             size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_CTERM, &mf, 0,
                                      1, primal_name);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             size_type region, const std::string &primal_name,
                             size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_REGION_CTERM, &mf, region,
                                      1, primal_name);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             const std::string &primal_name,const mesh_im &mim,
                             size_type region, size_type niter) {
    check_name_validity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_INFSUP, &mf, region,
                                      1, primal_name, &mim);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::disable_variable(const std::string &name) {
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
    it->second.is_disabled = true;
    for (VAR_SET::iterator itv = variables.begin();
         itv != variables.end(); ++itv) {
      if (((itv->second.filter & VDESCRFILTER_INFSUP) ||
           (itv->second.filter & VDESCRFILTER_CTERM))
          && (name.compare(itv->second.filter_var) == 0)) {
        itv->second.is_disabled = true;
      }
    }
    if (!act_size_to_be_done) resize_global_system();
  }

  void model::enable_variable(const std::string &name) {
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << name);
    it->second.is_disabled = false;
    for (VAR_SET::iterator itv = variables.begin();
         itv != variables.end(); ++itv) {
      if (((itv->second.filter & VDESCRFILTER_INFSUP) ||
           (itv->second.filter & VDESCRFILTER_CTERM))
          && (name.compare(itv->second.filter_var) == 0)) {
        itv->second.is_disabled = false;
      }
    }
    if (!act_size_to_be_done) resize_global_system();
  }

  void model::add_macro(const std::string &name, const std::string &expr)
  { check_name_validity(name); macros[name] = expr; }

  bool model::macro_exists(const std::string &name) const
  { return (macros.find(name) != macros.end()); }

  const std::string &model::get_macro(const std::string &name) const {
    std::map<std::string, std::string>::const_iterator it = macros.find(name);
    GMM_ASSERT1(it != macros.end(), "Undefined macro");
    return it->second;
  }


  void model::delete_brick(size_type ib) {
     GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
     valid_bricks.del(ib);
     active_bricks.del(ib);

     for  (size_type i = 0; i < bricks[ib].mims.size(); ++i) {
       const mesh_im *mim = bricks[ib].mims[i];
       bool found = false;
       for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
         for  (size_type j = 0; j < bricks[ibb].mims.size(); ++j)
           if (bricks[ibb].mims[j] == mim) found = true;
       }
       for(VAR_SET::iterator it2 = variables.begin();
           it2 != variables.end(); ++it2) {
         if (it2->second.is_fem_dofs &&
             (it2->second.filter & VDESCRFILTER_INFSUP) &&
             mim == it2->second.mim) found = true;
        }
       if (!found) sup_dependency(*mim);
     }

     is_linear_ = is_symmetric_ = is_coercive_ = true;
     for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
       is_linear_ = is_linear_ && bricks[ibb].pbr->is_linear();
       is_symmetric_ = is_symmetric_ &&  bricks[ibb].pbr->is_symmetric();
       is_coercive_ = is_coercive_ &&  bricks[ibb].pbr->is_coercive();
     }

     Neumann_SET::iterator it = Neumann_term_list.begin(), it2 = it;
     for (; it != Neumann_term_list.end(); it = it2) {
       it2++;
       if (it->first.second == ib) Neumann_term_list.erase(it);
     }

     bricks[ib] = brick_description();
  }

  void model::delete_variable(const std::string &varname) {
    for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
      for (size_type i = 0; i < bricks[ibb].vlist.size(); ++i)
        GMM_ASSERT1(varname.compare(bricks[ibb].vlist[i]),
                    "Cannot delete a variable which is still use by a brick");
      for (size_type i = 0; i < bricks[ibb].dlist.size(); ++i)
        GMM_ASSERT1(varname.compare(bricks[ibb].dlist[i]),
                    "Cannot delete a data which is still use by a brick");
    }

    VAR_SET::const_iterator it = variables.find(varname);
    GMM_ASSERT1(it != variables.end(), "Undefined variable " << varname);

    if (it->second.is_fem_dofs) {
      const mesh_fem *mf = it->second.mf;
      bool found = false;
      for(VAR_SET::iterator it2 = variables.begin();
          it2 != variables.end(); ++it2) {
        if (it != it2 && it2->second.is_fem_dofs && mf == it2->second.mf)
          found = true;
      }
      if (!found) sup_dependency(*mf);

      if (it->second.filter & VDESCRFILTER_INFSUP) {
        const mesh_im *mim = it->second.mim;
        found = false;
        for (dal::bv_visitor ibb(valid_bricks); !ibb.finished(); ++ibb) {
          for  (size_type j = 0; j < bricks[ibb].mims.size(); ++j)
            if (bricks[ibb].mims[j] == mim) found = true;
        }
        for(VAR_SET::iterator it2 = variables.begin();
            it2 != variables.end(); ++it2) {
          if (it != it2 && it2->second.is_fem_dofs &&
              (it2->second.filter & VDESCRFILTER_INFSUP) &&
              mim == it2->second.mim) found = true;
        }
        if (!found) sup_dependency(*mim);
      }
    }

    if (it->second.pim_data != 0) sup_dependency(*(it->second.pim_data));

    Neumann_SET::iterator itn = Neumann_term_list.begin(), itn2 = itn;
    for (; itn != Neumann_term_list.end(); itn = itn2) {
      itn2++;
      if (!(varname.compare(itn->first.first))) Neumann_term_list.erase(itn);
    }
    Neumann_terms_auxilliary_variables.erase(varname);

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


    for (size_type i = 0; i < bricks[ib].mims.size(); ++i)
      add_dependency(*(bricks[ib].mims[i]));

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

    for (size_type i=0; i < varnames.size(); ++i)
      GMM_ASSERT1(variables.find(varnames[i]) != variables.end(),
                  "Undefined model variable " << varnames[i]);
    cout << "dl == " << datanames << endl;
    for (size_type i=0; i < datanames.size(); ++i)
      GMM_ASSERT1(variables.find(datanames[i]) != variables.end(),
                  "Undefined model data or variable " << datanames[i]);

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
    for (size_type i=0; i < vl.size(); ++i)
      GMM_ASSERT1(variables.find(vl[i]) != variables.end(),
                  "Undefined model variable " << vl[i]);
  }

  void model::change_data_of_brick(size_type ib, const varnamelist &dl) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].dlist = dl;
    for (size_type i=0; i < dl.size(); ++i)
      GMM_ASSERT1(variables.find(dl[i]) != variables.end(),
                  "Undefined model variable " << dl[i]);
  }

  void model::change_mims_of_brick(size_type ib, const mimlist &ml) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].mims = ml;
    for (size_type i = 0; i < ml.size(); ++i) add_dependency(*(ml[i]));
  }

  void model::change_update_flag_of_brick(size_type ib, bool flag) {
    GMM_ASSERT1(valid_bricks[ib], "Inexistent brick");
    touch_brick(ib);
    bricks[ib].is_update_brick = flag;
  }

  void model::call_init_affine_dependent_variables(int version) {
    for (VAR_SET::iterator it = variables.begin();
         it != variables.end(); ++it)
      if (it->second.is_variable && it->second.ptsc) {
        if (version == 2)
          it->second.ptsc->init_affine_dependent_variables_precomputation(*this);
        else
          it->second.ptsc->init_affine_dependent_variables(*this);
      }
  }

  void model::shift_variables_for_time_integration(void) {
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

  void model::copy_init_time_derivative(void) {

    for (VAR_SET::iterator it = variables.begin();
         it != variables.end(); ++it)
      if (it->second.is_variable && it->second.ptsc) {

        std::string name_v, name_previous_v;
        it->second.ptsc->time_derivative_to_be_intialized(name_v,
                                                          name_previous_v);

        if (name_v.size()) {
          if (is_complex()) {
            model_complex_plain_vector v0 = this->complex_variable(name_v);
            gmm::copy(v0, this->set_complex_variable(name_previous_v));
          } else {
            const model_real_plain_vector &v0 = this->real_variable(name_v);
            gmm::copy(v0, this->set_real_variable(name_previous_v));
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

      virtual void time_derivative_to_be_intialized
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
    ptime_scheme ptsc = new first_order_theta_method_scheme(md, varname,theta);
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

    virtual void time_derivative_to_be_intialized
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
    ptime_scheme ptsc = new second_order_theta_method_scheme(md,varname,theta);
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
        scalar_type b2 = dt*(1-gamma/(scalar_type(2)*beta));

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

      virtual void time_derivative_to_be_intialized
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
    ptime_scheme ptsc = new Newmark_scheme(md, varname, beta, gamma);
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
        if (!(active_bricks[i])) ost << " (desactivated)";
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

  /**takes a list (more often it's a std::vector) of matrices
  or vectors, creates an empty copy on each thread. When the
  thread computations are done (in the destructor), accumulates
  the assembled copies into the original. Note: the matrices or
  vectors in the list are gmm::cleared, deleting the content
  in the constructor*/
  template <class CONTAINER_LIST> class list_distro
  {
    CONTAINER_LIST& original_list;
    omp_distribute<CONTAINER_LIST> distributed_list;
    typedef typename CONTAINER_LIST::value_type value_type;

    void build_distro(gmm::abstract_matrix)
    {
      //intentionally skipping thread 0, as list_distro will
      //use original_list for it
      for(size_type thread = 1; thread < num_threads(); thread++)
      {
        auto it_original = original_list.begin();
        auto it_distributed = distributed_list(thread).begin();
        for(;it_original != original_list.end(); ++it_original, ++it_distributed)
        {
          gmm::resize(*it_distributed, gmm::mat_nrows(*it_original),gmm::mat_ncols(*it_original));
        }
      }
    }

    void build_distro(gmm::abstract_vector)
    {
      //.. skipping thread 0 ..
      for(size_type thread = 1; thread < num_threads(); thread++)
      {
        auto it_original = original_list.begin();
        auto it_distributed = distributed_list(thread).begin();
        for(;it_original != original_list.end(); ++it_original, ++it_distributed)
        {
          gmm::resize(*it_distributed, gmm::vect_size(*it_original));
        }
      }
    }

    bool not_multithreaded() const { return num_threads() == 1; }

  public:

    list_distro(CONTAINER_LIST& l) : original_list(l)
    {
      if (not_multithreaded()) return;

      for(size_type thread=1; thread<num_threads(); thread++)
          distributed_list(thread).resize(original_list.size());

      build_distro(typename gmm::linalg_traits<value_type>::linalg_type());
    }

    operator CONTAINER_LIST&()
    {
      if (not_multithreaded() || this_thread() == 0) return original_list;
      else return distributed_list(this_thread());
    }

    ~list_distro()
    {
      if (not_multithreaded()) return;

      GMM_ASSERT1(!me_is_multithreaded_now(),
                  "List accumulation should not run in parallel");

      for(size_type thread = 1; thread < num_threads(); thread++)
      {
        auto it_original=original_list.begin();
        auto it_distributed = distributed_list(thread).begin();
        for(; it_original != original_list.end(); ++it_original, ++it_distributed)
                gmm::add(*it_distributed, *it_original);
        }
      }
  };


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
      { //brackets are needed because list_distro has constructor/destructor
        //semantics (as in RAII)
        list_distro<complex_matlist> cmatlist(brick.cmatlist);
        list_distro<complex_veclist> cveclist(brick.cveclist[rhs_ind]);
        list_distro<complex_veclist> cveclist_sym(brick.cveclist_sym[rhs_ind]);

        /*running the assembly in parallel*/
        gmm::standard_locale locale;
        open_mp_is_running_properly check;
        #pragma omp parallel default(shared)
        {
          brick.pbr->asm_complex_tangent_terms(*this, ib, brick.vlist,
                                               brick.dlist, brick.mims,
                                               cmatlist,
                                               cveclist,
                                               cveclist_sym,
                                               brick.region, version);
        }

      }
      brick.pbr->complex_post_assembly_in_serial(*this, ib, brick.vlist,
                                                 brick.dlist, brick.mims,
                                                 brick.cmatlist,
                                                 brick.cveclist[rhs_ind],
                                                 brick.cveclist_sym[rhs_ind],
                                                 brick.region, version);

      if (brick.is_update_brick) //contributions of pure update bricks must be deleted
      {
        for (size_type i = 0; i < brick.cmatlist.size(); ++i)
        {
          gmm::clear(brick.cmatlist[i]);
        }

        for (size_type i = 0; i < brick.cveclist.size(); ++i)
          for (size_type j = 0; j < brick.cveclist[i].size(); ++j)
          {
            gmm::clear(brick.cveclist[i][j]);
          }

        for (size_type i = 0; i < brick.cveclist_sym.size(); ++i)
          for (size_type j = 0; j < brick.cveclist_sym[i].size(); ++j)
          {
            gmm::clear(brick.cveclist_sym[i][j]);
          }
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
        list_distro<real_matlist> rmatlist(brick.rmatlist);
        list_distro<real_veclist> rveclist(brick.rveclist[rhs_ind]);
        list_distro<real_veclist> rveclist_sym(brick.rveclist_sym[rhs_ind]);

        /*running the assembly in parallel*/
        gmm::standard_locale locale;
        open_mp_is_running_properly check;
        #pragma omp parallel default(shared)
        {
          brick.pbr->asm_real_tangent_terms(*this, ib, brick.vlist,
                                            brick.dlist, brick.mims,
                                            rmatlist,
                                            rveclist,
                                            rveclist_sym,
                                            brick.region,
                                            version);
        }
      }
      brick.pbr->real_post_assembly_in_serial(*this, ib, brick.vlist,
                                              brick.dlist, brick.mims,
                                              brick.rmatlist,
                                              brick.rveclist[rhs_ind],
                                              brick.rveclist_sym[rhs_ind],
                                              brick.region, version);

      if (brick.is_update_brick) //contributions of pure update bricks must be deleted
      {
        for (size_type i = 0; i < brick.rmatlist.size(); ++i)
        {
          gmm::clear(brick.rmatlist[i]);
        }

        for (size_type i = 0; i < brick.rveclist.size(); ++i)
          for (size_type j = 0; j < brick.rveclist[i].size(); ++j)
          {
            gmm::clear(brick.rveclist[i][j]);
          }

        for (size_type i = 0; i < brick.rveclist_sym.size(); ++i)
          for (size_type j = 0; j < brick.rveclist_sym[i].size(); ++j)
          {
            gmm::clear(brick.rveclist_sym[i][j]);
          }
      }
    }
  }


  void model::set_dispatch_coeff(void) {
    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];
      if (brick.pdispatch)
        brick.pdispatch->set_dispatch_coeff(*this, ib);

    }
  }

  void model::first_iter(void) {
    context_check(); if (act_size_to_be_done) actualize_sizes();
    for (VAR_SET::iterator it = variables.begin(); it != variables.end(); ++it)
      it->second.clear_temporaries();

    set_dispatch_coeff();

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];
      bool cplx = is_complex() && brick.pbr->is_complex();
      if (brick.pdispatch) {
        if (cplx)
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

  void model::next_iter(void) {
    context_check(); if (act_size_to_be_done) actualize_sizes();
    set_dispatch_coeff();

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];
      bool cplx = is_complex() && brick.pbr->is_complex();
      if (brick.pdispatch) {
        if (cplx)
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

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
         ++it) {
      for (size_type i = 1; i < it->second.n_iter; ++i) {
        if (is_complex())
          gmm::copy(it->second.complex_value[i-1],
                    it->second.complex_value[i]);
        else
          gmm::copy(it->second.real_value[i-1],
                    it->second.real_value[i]);
      }
      if (it->second.n_iter > 1) it->second.v_num_data = act_counter();
    }
  }

  bool model::is_var_newer_than_brick(const std::string &varname,
                                      size_type ib) const {
    const brick_description &brick = bricks[ib];
    var_description &vd = variables[varname];
    return (vd.v_num > brick.v_num || vd.v_num_data > brick.v_num);
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





  void model::auxilliary_variables_of_Neumann_terms
  (const std::string &varname, std::vector<std::string> &aux_vars) const {
    std::map<std::string, std::vector<std::string> >::const_iterator
      it = Neumann_terms_auxilliary_variables.find(varname);
    if (it !=  Neumann_terms_auxilliary_variables.end())
      aux_vars = it->second;
    else
      aux_vars.resize(0);
  }

  /* Pb on this function: depend only on the variable and not on the term
     and brick. Not well managed at brick deletion.
  */
  void model::add_auxilliary_variables_of_Neumann_terms
  (const std::string &varname,
   const std::vector<std::string> &aux_vars) const {

    for (size_type i = 0; i < aux_vars.size(); ++i) {
      bool found = false;
      for (size_type j = 0;
           j < Neumann_terms_auxilliary_variables[varname].size(); ++j)
        if (!(Neumann_terms_auxilliary_variables[varname][j].compare
              (aux_vars[i])))
          found = true;
      if (!found)
        Neumann_terms_auxilliary_variables[varname].push_back(aux_vars[i]);
    }
  }

  void model::add_auxilliary_variables_of_Neumann_terms
  (const std::string &varname, const std::string &aux_var) const {
    std::vector<std::string> aux_vars(1,  aux_var);
    add_auxilliary_variables_of_Neumann_terms(varname, aux_vars);
  }

  size_type
  model::check_Neumann_terms_consistency(const std::string &varname) const {

    dal::bit_vector bnum;
    Neumann_SET::const_iterator it = Neumann_term_list.begin();
    for (; it != Neumann_term_list.end(); ++it) bnum.add(it->first.second);

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      if (bricks[ib].pbr->has_Neumann_term() && !(bnum.is_in(ib))) {
        for (size_type j = 0; j < bricks[ib].vlist.size(); ++j)
          if (!(bricks[ib].vlist[j].compare(varname))) return ib;
      }
    }
    return size_type(-1);

  }

  bool model::check_Neumann_terms_linearity(const std::string &varname) const {

    Neumann_SET::const_iterator it
      = Neumann_term_list.lower_bound(Neumann_pair(varname, 0));

    while (it != Neumann_term_list.end()
           && !(it->first.first.compare(varname))) {
      if (!(bricks[it->first.second].pbr->is_linear())) return false;
    }
    return true;
  }


  void model::compute_Neumann_terms(int version, const std::string &varname,
                                    const mesh_fem &mfvar,
                                    const model_real_plain_vector &var,
                                    fem_interpolation_context &ctx,
                                    base_small_vector &n,
                                    bgeot::base_tensor &t) const {

    // The output tensor has to have the right size. No verification.
    Neumann_SET::const_iterator it
      = Neumann_term_list.lower_bound(Neumann_pair(varname, 0));

    gmm::clear(t.as_vector());
    while (it != Neumann_term_list.end()
           && !(it->first.first.compare(varname))) {
      if (active_bricks.is_in(it->first.second))
        it->second->compute_Neumann_term(version, mfvar, var, ctx, n, t);
      ++it;
    }
  }

  void model::compute_auxilliary_Neumann_terms
  (int version, const std::string &varname,
   const mesh_fem &mfvar, const model_real_plain_vector &var,
   const std::string &aux_varname,
   fem_interpolation_context &ctx, base_small_vector &n,
   bgeot::base_tensor &t) const {

    // The output tensor has to have the right size. No verification.
    Neumann_SET::const_iterator it
      = Neumann_term_list.lower_bound(Neumann_pair(varname, 0));

    gmm::clear(t.as_vector());
    while (it != Neumann_term_list.end()
           && !(it->first.first.compare(varname))) {
      if (active_bricks.is_in(it->first.second)) {
        size_type ind = size_type(-1);
        for (size_type i = 0; i < it->second->auxilliary_variables.size(); ++i)
          if (!(aux_varname.compare(it->second->auxilliary_variables[i])))
            ind = i;
        if (ind != size_type(-1))
          it->second->compute_Neumann_term(version,mfvar,var,ctx,n,t,ind+1);
      }
      ++it;
    }
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
      if (vd.v_num_iter[ind] <= vd.v_num_data) {
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

  const model_real_sparse_matrix &model::linear_real_matrix_term
  (size_type ib, size_type iterm) {
    GMM_ASSERT1(bricks[ib].tlist[iterm].is_matrix_term,
                "Not a matrix term !");
    GMM_ASSERT1(bricks[ib].pbr->is_linear(), "Nonlinear term !");
    return bricks[ib].rmatlist[iterm];
  }

  const model_complex_sparse_matrix &model::linear_complex_matrix_term
  (size_type ib, size_type iterm) {
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

    // check variable list to test if a mesh_fem as changed.
    if (!tobecomputed ) {
      for (size_type i = 0; i < brick.vlist.size() && !tobecomputed; ++i) {
        var_description &vd = variables[brick.vlist[i]];
        if (vd.v_num > brick.v_num) tobecomputed = true;
      }
    }

    // check data list to test if a vector value of a data has changed.
    for (size_type i = 0; i < brick.dlist.size() && !tobecomputed; ++i) {
      var_description &vd = variables[brick.dlist[i]];
      if (vd.v_num > brick.v_num || vd.v_num_data > brick.v_num) {
        tobecomputed = true;
        version = build_version(version | BUILD_ON_DATA_CHANGE);
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

  void model::update_affine_dependent_variables(void) {
    for (VAR_SET::iterator it = variables.begin(); it != variables.end(); ++it)
      if (it->second.is_affine_dependent) {
        VAR_SET::iterator it2 = variables.find(it->second.org_name);
        if (it->second.size() != it2->second.size())
          it->second.set_size(it2->second.size());
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
        it->second.v_num_data = std::max(it->second.v_num_data,
                                         it2->second.v_num_data);
      }
  }

  void model::assembly(build_version version) {

#if GETFEM_PARA_LEVEL > 1
    double t_ref = MPI_Wtime();
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
        VAR_SET::iterator it1, it2;
        if (!isg) {
          it1 = variables.find(term.var1);
          GMM_ASSERT1(it1->second.is_variable, "Assembly of data not allowed");
          I1 = it1->second.I;
        }
        if (term.is_matrix_term && !isg) {
          it2 = variables.find(term.var2);
          I2 = it2->second.I;
          if (!(it2->second.is_variable)) {
            std::string vorgname = sup_previous_and_dot_to_varname(term.var2);
            VAR_SET::iterator it3 = variables.find(vorgname);
            GMM_ASSERT1(it3->second.is_variable,
                        "Assembly of data not allowed");
            I2 = it3->second.I;
            isprevious = true;
          }
          alpha *= it1->second.alpha * it2->second.alpha;
          alpha1 *= it1->second.alpha;
          alpha2 *= it2->second.alpha;
        }

        if (cplx) {
          if (term.is_matrix_term && (version & BUILD_MATRIX) && !isprevious
              && (isg || (!(it1->second.is_disabled) && !(it2->second.is_disabled)))) {
            gmm::add(gmm::scaled(brick.cmatlist[j], alpha),
                     gmm::sub_matrix(cTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first()) {
              gmm::add(gmm::scaled(gmm::transposed(brick.cmatlist[j]), alpha),
                       gmm::sub_matrix(cTM, I2, I1));
            }
          }
          if (version & BUILD_RHS) {
            if (isg || !(it1->second.is_disabled)) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.cveclist[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I1));
              }
              else {
                gmm::add(gmm::scaled(brick.cveclist[0][j],
                                     complex_type(alpha1)),
                         gmm::sub_vector(crhs, I1));
              }
            }
            if (term.is_matrix_term && brick.pbr->is_linear() && is_linear()) {
              if (it2->second.is_affine_dependent
                  && !(it1->second.is_disabled))
                gmm::mult_add(brick.cmatlist[j],
                              gmm::scaled(it2->second.affine_complex_value,
                                          complex_type(-alpha1)),
                              gmm::sub_vector(crhs, I1));
              if (term.is_symmetric && I1.first() != I2.first()
                  && it1->second.is_affine_dependent
                  && !(it2->second.is_disabled)) {
                gmm::mult_add(gmm::conjugated(brick.cmatlist[j]),
                              gmm::scaled(it1->second.affine_complex_value,
                                          complex_type(-alpha2)),
                              gmm::sub_vector(crhs, I2));
              }
            }
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
              if (!(it1->second.is_disabled))
                gmm::mult_add(brick.cmatlist[j],
                              gmm::scaled(it2->second.complex_value[0],
                                          complex_type(-alpha1)),
                              gmm::sub_vector(crhs, I1));
            }
            if (term.is_symmetric && I1.first() != I2.first() &&
                !(it2->second.is_disabled)) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.cveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I2));
              } else {
                gmm::add(gmm::scaled(brick.cveclist_sym[0][j],
                                     complex_type(alpha2)),
                         gmm::sub_vector(crhs, I2));
              }
               if (brick.pbr->is_linear()
                   && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
                 gmm::mult_add(gmm::conjugated(brick.cmatlist[j]),
                            gmm::scaled(it1->second.complex_value[0],
                                        complex_type(-alpha2)),
                            gmm::sub_vector(crhs, I2));
               }
            }
          }
        } else if (is_complex()) {
          if (term.is_matrix_term && (version & BUILD_MATRIX) && !isprevious
              && (isg || (!(it1->second.is_disabled) && !(it2->second.is_disabled)))) {
            gmm::add(gmm::scaled(brick.rmatlist[j], alpha),
                     gmm::sub_matrix(cTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first()) {
              gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), alpha),
                       gmm::sub_matrix(cTM, I2, I1));
            }
          }
          if (version & BUILD_RHS) {
            if (isg || !(it1->second.is_disabled)) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I1));
              }
              else {
                gmm::add(gmm::scaled(brick.rveclist[0][j], alpha1),
                         gmm::sub_vector(crhs, I1));
              }
            }
            if (term.is_matrix_term && brick.pbr->is_linear() && is_linear()) {
              if (it2->second.is_affine_dependent
                  && !(it1->second.is_disabled))
                gmm::mult_add(brick.rmatlist[j],
                              gmm::scaled(it2->second.affine_complex_value,
                                          complex_type(-alpha1)),
                              gmm::sub_vector(crhs, I1));
              if (term.is_symmetric && I1.first() != I2.first()
                  && it1->second.is_affine_dependent
                  && !(it2->second.is_disabled)) {
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                              gmm::scaled(it1->second.affine_complex_value,
                                          complex_type(-alpha2)),
                              gmm::sub_vector(crhs, I2));
              }
            }
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
              if (!(it1->second.is_disabled))
                gmm::mult_add(brick.rmatlist[j],
                              gmm::scaled(it2->second.complex_value[0],
                                          complex_type(-alpha1)),
                              gmm::sub_vector(crhs, I1));
            }
            if (term.is_symmetric && I1.first() != I2.first() &&
                !(it2->second.is_disabled)) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I2));
              }
              else {
                gmm::add(gmm::scaled(brick.rveclist_sym[0][j], alpha2),
                         gmm::sub_vector(crhs, I2));
              }
              if (brick.pbr->is_linear()
                  && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                             gmm::scaled(it1->second.complex_value[0],
                                          complex_type(-alpha2)),
                              gmm::sub_vector(crhs, I2));
              }
            }
          }
        } else {
          if (term.is_matrix_term && (version & BUILD_MATRIX) && !isprevious
              && (isg || (!(it1->second.is_disabled)
                          && !(it2->second.is_disabled)))) {
            gmm::add(gmm::scaled(brick.rmatlist[j], alpha),
                     gmm::sub_matrix(rTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first()) {
              gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), alpha),
                       gmm::sub_matrix(rTM, I2, I1));
            }
          }
          if (version & BUILD_RHS) {
            if (isg || !(it1->second.is_disabled)) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(rrhs, I1));
              }
              else {
                gmm::add(gmm::scaled(brick.rveclist[0][j], alpha1),
                         gmm::sub_vector(rrhs, I1));
              }
            }
            if (term.is_matrix_term && brick.pbr->is_linear() && is_linear()) {
              if (it2->second.is_affine_dependent
                  && !(it1->second.is_disabled))
                gmm::mult_add(brick.rmatlist[j],
                              gmm::scaled(it2->second.affine_real_value,
                                          -alpha1),
                              gmm::sub_vector(rrhs, I1));
              if (term.is_symmetric && I1.first() != I2.first()
                  && it1->second.is_affine_dependent
                  && !(it2->second.is_disabled)) {
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                              gmm::scaled(it1->second.affine_real_value,
                                          -alpha2),
                              gmm::sub_vector(rrhs, I2));
              }
            }
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
              if (!(it1->second.is_disabled))
                gmm::mult_add(brick.rmatlist[j],
                              gmm::scaled(it2->second.real_value[0],
                                          -alpha1),
                              gmm::sub_vector(rrhs, I1));
            }
            if (term.is_symmetric && I1.first() != I2.first() &&
                !(it2->second.is_disabled)) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(rrhs, I2));
              }
              else {
                gmm::add(gmm::scaled(brick.rveclist_sym[0][j], alpha2),
                         gmm::sub_vector(rrhs, I2));
              }
              if (brick.pbr->is_linear()
                  && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                              gmm::scaled(it1->second.real_value[0], -alpha2),
                              gmm::sub_vector(rrhs, I2));
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

    if (version & BUILD_RHS) {
      if (is_complex()) MPI_SUM_VECTOR(crhs); else MPI_SUM_VECTOR(rrhs);
    }

    // Generic expressions
    if (generic_expressions.size()) {
      GMM_TRACE2("Global generic assembly");
      ga_workspace workspace(*this);

      for (std::list<gen_expr>::iterator it = generic_expressions.begin();
             it != generic_expressions.end(); ++it) {
        workspace.add_expression(it->expr, it->mim, it->region);
      }

      if (version & BUILD_RHS) {
        if (is_complex()) {
          GMM_ASSERT1(false, "to be done");
        } else {
          model_real_plain_vector residual(gmm::vect_size(rrhs));
          workspace.set_assembled_vector(residual);
          workspace.assembly(1);
          gmm::add(gmm::scaled(residual, scalar_type(-1)), rrhs);
        }
      }

      if (version & BUILD_MATRIX) {
        if (is_complex()) {
          GMM_ASSERT1(false, "to be done");
        } else {
          workspace.set_assembled_matrix(rTM);
          workspace.assembly(2);
        }
      }
    }

    // Post simplification for dof constraints
    if ((version & BUILD_RHS) || (version & BUILD_MATRIX)) {
      if (is_complex()) {
        std::vector<size_type> dof_indices;
        std::vector<complex_type> dof_pr_values;
        std::vector<complex_type> dof_go_values;
        std::map<std::string, complex_dof_constraints_var>::const_iterator it;

        for (it = complex_dof_constraints.begin();
             it != complex_dof_constraints.end(); ++it) {
          const gmm::sub_interval &I = interval_of_variable(it->first);
          const model_complex_plain_vector &V = complex_variable(it->first);
          complex_dof_constraints_var::const_iterator itv;
          for (itv = it->second.begin(); itv != it->second.end(); ++itv) {
            dof_indices.push_back(itv->first + I.first());
            dof_go_values.push_back(itv->second);
            dof_pr_values.push_back(V[itv->first]);
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
                  model_complex_plain_vector vv(gmm::vect_size(rrhs));
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
      } else {
        std::vector<size_type> dof_indices;
        std::vector<scalar_type> dof_pr_values;
        std::vector<scalar_type> dof_go_values;
        std::map<std::string, real_dof_constraints_var>::const_iterator it;

        for (it = real_dof_constraints.begin();
             it != real_dof_constraints.end(); ++it) {
          const gmm::sub_interval &I = interval_of_variable(it->first);
          const model_real_plain_vector &V = real_variable(it->first);
          real_dof_constraints_var::const_iterator itv;
          for (itv = it->second.begin(); itv != it->second.end(); ++itv) {
            dof_indices.push_back(itv->first + I.first());
            dof_go_values.push_back(itv->second);
            dof_pr_values.push_back(V[itv->first]);
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
      approx_external_load_ = MPI_SUM_SCALAR(approx_external_load_);
    }


    #if GETFEM_PARA_LEVEL > 1
    // int rk; MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    if (MPI_IS_MASTER()) cout << "Assembly time " << MPI_Wtime()-t_ref << endl;
    #endif

  }


  const mesh_fem &model::mesh_fem_of_variable(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    return it->second.associated_mf();
  }

  const mesh_fem *model::pmesh_fem_of_variable(const std::string &name) const {
    VAR_SET::const_iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    return it->second.passociated_mf();
  }

  const model_real_plain_vector &
  model::real_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(!complex_version, "This model is a complex one");
    context_check();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.is_fem_dofs) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size(it->second.mf->nb_dof()*it->second.qdim);
    }
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.real_value[niter];
  }

  const model_complex_plain_vector &
  model::complex_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    context_check();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.is_fem_dofs) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size(it->second.mf->nb_dof()*it->second.qdim);
    }
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter  > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.complex_value[niter];
  }

  model_real_plain_vector &
  model::set_real_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(!complex_version, "This model is a complex one");
    context_check();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.is_fem_dofs) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size(it->second.mf->nb_dof()*it->second.qdim);
    }
    it->second.v_num_data = act_counter();
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.real_value[niter];
  }

  model_complex_plain_vector &
  model::set_complex_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    context_check();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.is_fem_dofs) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size(it->second.mf->nb_dof()*it->second.qdim);
    }
    it->second.v_num_data = act_counter();
    if (niter == size_type(-1)) niter = it->second.default_iter;
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
    if (act_size_to_be_done && it->second.is_fem_dofs) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size(it->second.mf->nb_dof()*it->second.qdim);
    }
    it->second.v_num_data = act_counter();
    return it->second.affine_real_value;
  }

  model_complex_plain_vector &
  model::set_complex_constant_part(const std::string &name) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    context_check();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (act_size_to_be_done && it->second.is_fem_dofs) {
      if (it->second.filter != VDESCRFILTER_NO)
        actualize_sizes();
      else
        it->second.set_size(it->second.mf->nb_dof()*it->second.qdim);
    }
    it->second.v_num_data = act_counter();
    return it->second.affine_complex_value;
  }

  void model::check_brick_stiffness_rhs(size_type ind_brick) const {

    const brick_description &brick = bricks[ind_brick];
    update_brick(ind_brick, model::BUILD_ALL);

    brick.pbr->check_stiffness_matrix_and_rhs(*this, ind_brick, brick.tlist,
                      brick.vlist, brick.dlist, brick.mims, brick.rmatlist,
                      brick.rveclist[0], brick.rveclist_sym[0], brick.region);
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

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &,
                                        const model::varnamelist &,
                                        const model::mimlist &mims,
                                        model::real_matlist &,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {
      GMM_ASSERT1(vecl.size() ==  vl_test1.size()
                  + ((directdataname.size() == 0) ? 0 : 1), "Wrong number "
                  "of terms for Generic source term assembly brick ");
      GMM_ASSERT1(mims.size() == 1, "Generic source term assembly brick "
                  "needs one and only one mesh_im");
      GMM_TRACE2("Generic source term assembly");

      gmm::clear(vecl[0]);

      if (expr.size()) {
        size_type nbgdof = md.nb_dof();
        ga_workspace workspace(md, true);
        GMM_TRACE2(name << ": generic source term assembly");
        workspace.add_expression(expr, *(mims[0]), region);  
        model::varnamelist vlmd; md.variable_list(vlmd);
        for (size_type i = 0; i < vlmd.size(); ++i)
          if (md.is_disabled_variable(vlmd[i]))
            nbgdof = std::max(nbgdof, 
                              workspace.interval_of_variable(vlmd[i]).last());
        GMM_TRACE2(name << ": generic matrix assembly");
        model_real_plain_vector V(nbgdof);
        workspace.set_assembled_vector(V);
        workspace.assembly(1);
        md.add_external_load(ib, gmm::vect_norm1(V));
        for (size_type i = 0; i < vl_test1.size(); ++i) {
          gmm::copy(gmm::sub_vector
                    (V, workspace.interval_of_variable(vl_test1[i])), vecl[i]);
        }
      }

      if (directvarname.size()) {
        gmm::copy(md.real_variable(directdataname), vecl.back());
        md.add_external_load(ib, gmm::vect_norm1(vecl.back()));
      }
    }

    gen_source_term_assembly_brick(const std::string &expr_,
                                   std::string brickname,
                                   const model::varnamelist &vl_test1_,
                                   const std::string &directvarname_,
                                   const std::string &directdataname_)
      : vl_test1(vl_test1_) {
      if (brickname.size() == 0)
        brickname = "Generic source term assembly brick";
      expr = expr_;
      set_flags(brickname, true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, false /* is complex */);
      directvarname = directvarname_; directdataname = directdataname_;
    }

  };

  size_type add_source_term_generic_assembly_brick
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   std::string brickname, std::string directvarname,
   const std::string &directdataname, bool return_if_nonlin) {

    ga_workspace workspace(md);
    size_type order = workspace.add_expression(expr, mim, region);
    GMM_ASSERT1(order <= 1, "Wrong order for a source term");
    model::varnamelist vl, vl_test1, vl_test2, dl;
    bool is_lin = workspace.used_variables(vl, vl_test1, vl_test2, dl, 1);
    if (!is_lin && return_if_nonlin) return size_type(-1);
    GMM_ASSERT1(is_lin, "Nonlinear term");

    cout << "dl = " << dl << endl;
    cout << "directvarname = " << directvarname << endl;

    if (directdataname.size()) {
      vl.push_back(directvarname);
      dl.push_back(directdataname);
    } else directvarname = "";

    cout << "dl = " << dl << endl;
    
    pbrick pbr = new gen_source_term_assembly_brick
      (expr, brickname, vl_test1, directvarname, directdataname);
    model::termlist tl;

    for (size_type i = 0; i < vl_test1.size(); ++i)
      tl.push_back(model::term_description(vl_test1[i]));
    if (directdataname.size())
      tl.push_back(model::term_description(directvarname));

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
    cout << "ok " << endl;
  }

  // ----------------------------------------------------------------------
  //
  // Linear generic assembly brick
  //
  // ----------------------------------------------------------------------

  struct gen_linear_assembly_brick : public virtual_brick {

    std::string expr;
    model::varnamelist vl_test1, vl_test2;

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
        size_type nbgdof = md.nb_dof();
        ga_workspace workspace(md, true);
        workspace.add_expression(expr, *(mims[0]), region);
        model::varnamelist vlmd; md.variable_list(vlmd);
        for (size_type i = 0; i < vlmd.size(); ++i)
          if (md.is_disabled_variable(vlmd[i]))
            nbgdof = std::max(nbgdof, 
                              workspace.interval_of_variable(vlmd[i]).last());
        GMM_TRACE2(name << ": generic matrix assembly");
        model_real_sparse_matrix R(nbgdof, nbgdof); 
        workspace.set_assembled_matrix(R);
        workspace.assembly(2);
        for (size_type i = 0; i < vl_test1.size(); ++i) {
          scalar_type alpha = scalar_type(1)
            / ( workspace.factor_of_variable(vl_test1[i]) *
                workspace.factor_of_variable(vl_test2[i]));
          gmm::copy(gmm::scaled(gmm::sub_matrix
                    (R, workspace.interval_of_variable(vl_test1[i]),
                     workspace.interval_of_variable(vl_test2[i])), alpha),
                    matl[i]);
        }
      }
    }

    gen_linear_assembly_brick(const std::string &expr_, bool is_sym,
                              bool is_coer, std::string brickname,
                              const model::varnamelist &vl_test1_,
                              const model::varnamelist &vl_test2_)
      : vl_test1(vl_test1_), vl_test2(vl_test2_) {
      if (brickname.size() == 0) brickname = "Generic linear assembly brick";
      expr = expr_;
      set_flags(brickname, true /* is linear*/,
                is_sym /* is symmetric */, is_coer /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };

  size_type add_linear_generic_assembly_brick
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   bool is_sym, bool is_coercive, std::string brickname,
   bool return_if_nonlin) {

    ga_workspace workspace(md, true);
    size_type order = workspace.add_expression(expr, mim, region);
    model::varnamelist vl, vl_test1, vl_test2, dl;
    bool is_lin = workspace.used_variables(vl, vl_test1, vl_test2, dl, 2);

    if (!is_lin && return_if_nonlin) return size_type(-1);
    GMM_ASSERT1(is_lin, "Nonlinear term");
    if (order == 0) { is_coercive = is_sym = true; }

    std::string const_expr= workspace.extract_constant_term(mim.linked_mesh());
    if (const_expr.size()) {
      add_source_term_generic_assembly_brick
        (md, mim, const_expr, region, brickname+" (source term)");
    }

    // GMM_ASSERT1(order <= 1,
    //             "This brick does not support a second order term");
    
    if (vl_test1.size()) {
      pbrick pbr = new gen_linear_assembly_brick(expr, is_sym, is_coercive,
                                                 brickname,
                                                 vl_test1, vl_test2);
      model::termlist tl;
      for (size_type i = 0; i < vl_test1.size(); ++i)
        tl.push_back(model::term_description(vl_test1[i], vl_test2[i], false));

      return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
    }
    return size_type(-1);
  }


  // ----------------------------------------------------------------------
  //
  // Nonlinear generic assembly brick
  //
  // ----------------------------------------------------------------------

  struct gen_nonlinear_assembly_brick : public virtual_brick {

    std::string expr;

    virtual void asm_real_tangent_terms(const model &md, size_type ,
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
      md.add_generic_expression(expr, *(mims[0]), region);
    }

    gen_nonlinear_assembly_brick(const std::string &expr_, bool is_sym,
                                 bool is_coer, std::string brickname = "") {
      if (brickname.size() == 0) brickname = "Generic linear assembly brick";
      expr = expr_;
      set_flags(brickname, false /* is linear*/,
                is_sym /* is symmetric */, is_coer /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };

  size_type add_nonlinear_generic_assembly_brick
  (model &md, const mesh_im &mim, const std::string &expr, size_type region,
   bool is_sym, bool is_coercive, std::string brickname) {

    ga_workspace workspace(md);
    size_type order = workspace.add_expression(expr, mim, region);
    GMM_ASSERT1(order < 2, "Order two test functions (Test2) are not allowed"
                " in assembly string for nonlinear terms");
    model::varnamelist vl, vl_test1, vl_test2, ddl, dl;
    workspace.used_variables(vl, vl_test1, vl_test2, ddl, order);

    for (size_type i = 0; i < ddl.size(); ++i)
      if (md.is_true_data(ddl[i])) dl.push_back(ddl[i]);
      else vl.push_back(ddl[i]);

    if (order == 0) { is_coercive = is_sym = true; }
    pbrick pbr = new gen_nonlinear_assembly_brick(expr, is_sym, is_coercive,
                                                  brickname);
    model::termlist tl; // No term
    // tl.push_back(model::term_description(true, is_sym));
    // TODO to be changed.

    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  // ----------------------------------------------------------------------
  //
  // Generic elliptic brick
  //
  // ----------------------------------------------------------------------

  // Deprecated
  struct generic_elliptic_Neumann_elem_term : public Neumann_elem_term {

    const mesh_fem *mf_a;
    const model_real_plain_vector *A;

    mutable fem_interpolation_context ctx_a;
    mutable base_vector coeff, val;
    mutable base_matrix grad, G;

    void compute_Neumann_term
    (int version, const mesh_fem &mfvar, const model_real_plain_vector &var,
     fem_interpolation_context& ctx, base_small_vector &n,
     base_tensor &output, size_type /*auxilliary_ind*/ = 0) const {

      if (version == 3) return;  // No contribution because the term is linear

      const mesh &m = mfvar.linked_mesh();
      size_type N = m.dim(), Q = mfvar.get_qdim(), s = 1, cv=ctx.convex_num();

      if (A) {
        s = gmm::vect_size(*A);
        if (mf_a) s = s * mf_a->get_qdim() / mf_a->nb_dof();
      }
      gmm::resize(val, s);

      if (mf_a) {
        GMM_ASSERT1(!(mf_a->is_reduced()),
                    "Sorry, to be adapted for reduced mesh fems");

        if (!(ctx_a.have_pf()) || ctx_a.convex_num() != cv
            || (ctx_a.have_pfp() != ctx.have_pfp())
            || (ctx_a.have_pfp()
                && (&(ctx.pfp()->get_point_tab())
                    != &(ctx_a.pfp()->get_point_tab())))) {

          bgeot::vectors_to_base_matrix
            (G, mf_a->linked_mesh().points_of_convex(cv));

          pfem_precomp pfp = fem_precomp(mf_a->fem_of_element(cv),
                                         &(ctx.pfp()->get_point_tab()), 0);

          if (ctx.have_pfp())
            ctx_a = fem_interpolation_context
              (mf_a->linked_mesh().trans_of_convex(cv), pfp, ctx.ii(),
               G, cv, ctx.face_num());
          else
            ctx_a = fem_interpolation_context
              (mf_a->linked_mesh().trans_of_convex(cv),
               mf_a->fem_of_element(cv), ctx.xref(), G, cv, ctx.face_num());

        } else {
          if (ctx.have_pfp())  ctx_a.set_ii(ctx.ii());
          else ctx_a.set_xref(ctx.xref());
        }

        slice_vector_on_basic_dof_of_element(*mf_a, *A, cv, coeff);
        // coeff.resize(mf_a->nb_basic_dof_of_element(cv));
        // gmm::copy(gmm::sub_vector(*A, gmm::sub_index
        //                      (mfvar.ind_basic_dof_of_element(cv))), coeff);
        ctx_a.pf()->interpolation(ctx_a, coeff, val, dim_type(s));
      } else if (A) {
        gmm::copy(*A, val);
      } else {
        val[0] = scalar_type(1);
      }

      switch (version) {
      case 1:
        gmm::resize(grad, Q, N);
        slice_vector_on_basic_dof_of_element(mfvar, var, cv, coeff);
        // coeff.resize(mfvar.nb_basic_dof_of_element(cv));
        // gmm::copy(gmm::sub_vector(var, gmm::sub_index
        //                      (mfvar.ind_basic_dof_of_element(cv))), coeff);
        ctx.pf()->interpolation_grad(ctx, coeff, grad, dim_type(Q));

        if (s == 1)
          gmm::mult_add(grad, gmm::scaled(n, val[0]), output.as_vector());
        else if (s == N*N) {
          base_vector::const_iterator it = val.begin();
          for (size_type j = 0; j < N; ++j)
            for (size_type i = 0; i < N; ++i, ++it)
              for (size_type k = 0; k < Q; ++k)
                output[k] += (*it)*grad(k,j)*n[i];
        }
        else if (s == N*N*Q*Q) {
          base_vector::const_iterator it = val.begin();
          for (size_type l = 0; l < N; ++l)
            for (size_type k = 0; k < Q; ++k)
              for (size_type j = 0; j < N; ++j)
                for (size_type i = 0; i < Q; ++i, ++it)
                  output[i] += (*it) * grad(k, l) * n[j];
        }
        break;
      case 2:
        {
          base_tensor t;
          dim_type tdim = ctx.pf()->target_dim(), qmult = dim_type(Q) / tdim;
          size_type ndof = ctx.pf()->nb_dof(cv);
          // The return tensor is t(i,j,k) with 0<=i<ndof, 0<=j<target_dim,
          // 0<=k<dim. In order to iterate on the tensor values, i should
          // increment the faster, then j, then k.
          // If target_dim == qdim, grad(phi_i)(j,k) = t(i,j,k)
          // If target_dim == 1, grad(phi_i * e_l)(l,k) = t(i,1,k)
          // General case, psi_{i*qmult+l} = phi_i * e_l  and
          //    grad(psi_{i*qmult+l})(j+tdim*l,k) = t(i,j,k)
          ctx.pf()->real_grad_base_value(ctx, t);

          if (s == 1) {
//            for (size_type l = 0; l < qmult; ++l) {
//              for (size_type p = 0; p < Q; ++p) {
//                base_tensor::const_iterator it = t.begin();
//                for (size_type k = 0; k < Q; ++k)
//                  for (size_type j = 0; j < tdim; ++j)
//                    for (size_type i = 0; i < ndof; ++i, ++it) {
//                      size_type jj = j + tdim*l;
//                      if (p == jj) output(i*qmult+l, p) += val[0]*(*it)*n[k];
//                    }
//                GMM_ASSERT1(it ==  t.end(), "Internal error");
//              }
//            }
            if (Q == 1) {
              base_tensor::const_iterator it = t.begin();
              for (size_type k = 0; k < N; ++k)
                for (size_type i = 0; i < ndof; ++i, ++it)
                  output[i] += val[0]*(*it)*n[k];
              GMM_ASSERT1(it ==  t.end(), "Internal error");
            } else {
              for (size_type l = 0; l < qmult; ++l) {
                base_tensor::const_iterator it = t.begin();
                for (size_type k = 0; k < N; ++k)
                  for (size_type j = 0; j < tdim; ++j)
                    for (size_type i = 0; i < ndof; ++i, ++it) {
                      size_type jj = j + tdim*l;
                      output(i*qmult+l, jj) += val[0]*(*it)*n[k];
                    }
                GMM_ASSERT1(it ==  t.end(), "Internal error");
              }
            }
          } else if (s == N*N) {
            if (Q == 1) {
              base_tensor::const_iterator it = t.begin();
              for (size_type k = 0; k < N; ++k)
                for (size_type i = 0; i < ndof; ++i, ++it) {
                  for (size_type q = 0; q < N; ++q)
                    output[i] += val[q+k*N]*(*it)*n[q];
                }
              GMM_ASSERT1(it ==  t.end(), "Internal error");
            } else {
              for (size_type l = 0; l < qmult; ++l) {
                base_tensor::const_iterator it = t.begin();
                for (size_type k = 0; k < N; ++k)
                  for (size_type j = 0; j < tdim; ++j)
                    for (size_type i = 0; i < ndof; ++i, ++it) {
                      size_type jj = j + tdim*l;
                      for (size_type q = 0; q < N; ++q)
                        output(i*qmult+l, jj) += val[q+k*N]*(*it)*n[q];
                    }
                GMM_ASSERT1(it ==  t.end(), "Internal error");
              }
            }
          } else if (s == N*N*Q*Q) {
            for (size_type l = 0; l < qmult; ++l) {
              for (size_type p = 0; p < Q; ++p) {
                base_tensor::const_iterator it = t.begin();
                for (size_type k = 0; k < N; ++k)
                  for (size_type j = 0; j < tdim; ++j)
                    for (size_type i = 0; i < ndof; ++i, ++it) {
                      size_type jj = j + tdim*l;
                      for (size_type q = 0; q < N; ++q)
                        output(i*qmult+l, p)
                          += val[p+q*Q+jj*N*Q+k*N*Q*Q]*(*it)*n[q];
                    }
                GMM_ASSERT1(it ==  t.end(), "Internal error");
              }
            }
          }
        }
        break;
      }
    }

    generic_elliptic_Neumann_elem_term
    (const mesh_fem *mf_a_, const model_real_plain_vector *A_)
      : mf_a(mf_a_), A(A_) {}

  };



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

    virtual void real_post_assembly_in_serial(const model &md, size_type ib,
                                              const model::varnamelist &vl,
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
      pNeumann_elem_term pNt = new generic_elliptic_Neumann_elem_term(mf_a, A);
      md.add_Neumann_term(pNt, vl[0], ib);
    }

    virtual void complex_post_assembly_in_serial(const model &md, size_type ib,
                                              const model::varnamelist &vl,
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
      pNeumann_elem_term pNt = new generic_elliptic_Neumann_elem_term(mf_a, A);
      md.add_Neumann_term(pNt, vl[0], ib);
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

    generic_elliptic_brick(void) {
      set_flags("Generic elliptic", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */);
    }

  };

  size_type add_Laplacian_brick(model &md, const mesh_im &mim,
                                const std::string &varname,
                                size_type region) {
    if (md.is_complex()) {
      pbrick pbr = new generic_elliptic_brick;
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
      return add_linear_generic_assembly_brick(md, mim, expr, region, true,
                                               true, "Laplacian", false);
    }
  }

  size_type add_generic_elliptic_brick(model &md, const mesh_im &mim,
                                       const std::string &varname,
                                       const std::string &dataname,
                                       size_type region) {
    if (md.is_complex()) {
      pbrick pbr = new generic_elliptic_brick;
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          model::varnamelist(1, dataname), tl,
                          model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
      size_type qdim = mf_u.get_qdim();
      std::string expr;
      if (qdim == 1)
        expr = "(("+dataname+")*Grad_"+varname+").Grad_"+test_varname;
      else
        expr = "(("+dataname+")*Grad_"+varname+"):Grad_"+test_varname;
      size_type ib = add_linear_generic_assembly_brick
        (md, mim, expr, region, true, true, "Generic elliptic", true);
      if (ib == size_type(-1))
        ib = add_nonlinear_generic_assembly_brick
          (md, mim, expr, region, false, false,
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

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version) const {
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

      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type ib,
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

      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }


    source_term_brick(void) {
      set_flags("Source term", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
    }


  };

  size_type add_source_term_brick(model &md, const mesh_im &mim,
                                  const std::string &varname,
                                  const std::string &dataexpr,
                                  size_type region,
                                  const std::string &directdataname) {
    if (md.is_complex()) {
      pbrick pbr = new source_term_brick;
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
        ib = add_nonlinear_generic_assembly_brick
          (md, mim, "-("+expr+")", region, false, false,
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

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
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

      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type ib,
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
      md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
    }

    normal_source_term_brick(void) {
      set_flags("Normal source term", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
    }


  };

  size_type add_normal_source_term_brick(model &md, const mesh_im &mim,
                                         const std::string &varname,
                                         const std::string &dataexpr,
                                         size_type region) {
    if (md.is_complex()) {
      pbrick pbr = new normal_source_term_brick;
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
        if (H_version) {
          if (mf_H)
            asm_real_or_complex_1_param
              (*B, mim, mf_mult, *mf_H, *H, rg, (mf_u.get_qdim() == 1) ?
               "F=data(#2);"
               "M(#1,#3)+=comp(Base(#1).Base(#3).Base(#2))(:,:,i).F(i)"
               : "F=data(qdim(#1),qdim(#1),#2);"
               "M(#1,#3)+=comp(vBase(#1).vBase(#3).Base(#2))(:,i,:,j,k).F(i,j,k);", &mf_u);
          else {
            asm_real_or_complex_1_param
              (*B, mim, mf_mult, mf_u, *H, rg, (mf_u.get_qdim() == 1) ?
               "F=data(1);"
               "M(#1,#2)+=comp(Base(#1).Base(#2).F(1))"
               : "F=data(qdim(#1),qdim(#1));"
               "M(#1,#2)+=comp(vBase(#1).vBase(#2))(:,i,:,j).F(i,j);");
          }
        }
        else if (normal_component) {
          generic_assembly assem;
          if (mf_mult.get_qdim() == 1)
            assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1).Normal())(:,:,i,i);");
          else
            assem.set("M(#2,#1)+=comp(vBase(#2).mBase(#1).Normal())(:,i,:,i,j,j);");
          assem.push_mi(mim);
          assem.push_mf(mf_u);
          assem.push_mf(mf_mult);
          assem.push_mat(*B);
          assem.assembly(rg);
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

        md.add_external_load(ib, gmm::vect_norm1(vecl[0]));

        if (penalized && (&mf_mult != &mf_u))  {
          gmm::mult(gmm::transposed(rB), rV, vecl[0]);
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
          rV = model_real_plain_vector();
        } else if (penalized)
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }

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
          if (mf_H)
            asm_real_or_complex_1_param
              (*B, mim, mf_mult, *mf_H, *H, rg, (mf_u.get_qdim() == 1) ?
               "F=data(#2);"
               "M(#1,#3)+=sym(comp(Base(#1).Base(#3).Base(#2))(:,:,i).F(i))"
               : "F=data(qdim(#1),qdim(#1),#2);"
               "M(#1,#3)+=sym(comp(vBase(#1).vBase(#3).Base(#2))(:,i,:,j,k).F(i,j,k));", &mf_u);
          else
             asm_real_or_complex_1_param
              (*B, mim, mf_mult, mf_u, *H, rg, (mf_u.get_qdim() == 1) ?
               "F=data(1);"
               "M(#1,#2)+=sym(comp(Base(#1).Base(#2)).F(1))"
               : "F=data(qdim(#1),qdim(#1));"
               "M(#1,#2)+=sym(comp(vBase(#1).vBase(#2))(:,i,:,j).F(i,j));");
        }
        else if (normal_component) {
          generic_assembly assem;
          if (mf_mult.get_qdim() == 1)
            assem.set("M(#2,#1)+=comp(Base(#2).vBase(#1).Normal())(:,:,i,i);");
          else
            assem.set("M(#2,#1)+=comp(vBase(#2).mBase(#1).Normal())(:,i,:,i,j,j);");
          assem.push_mi(mim);
          assem.push_mf(mf_u);
          assem.push_mf(mf_mult);
          assem.push_mat(gmm::real_part(*B));
          assem.assembly(rg);
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

        md.add_external_load(ib, gmm::vect_norm1(vecl[0]));

        if (penalized && (&mf_mult != &mf_u))  {
          gmm::mult(gmm::transposed(cB), cV, vecl[0]);
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
          cV = model_complex_plain_vector();
        } else if (penalized)
          gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }
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
                false /* compute each time */, false /* has a Neumann term */);
    }
  };

  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname) {
    pbrick pbr = new Dirichlet_condition_brick(false, false, false);
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
    pbrick pbr = new Dirichlet_condition_brick(true, false, false, mf_mult);
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
    pbrick pbr = new Dirichlet_condition_brick(false, false, true);
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
    pbrick pbr = new Dirichlet_condition_brick(true, false, true, mf_mult);
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
    pbrick pbr = new Dirichlet_condition_brick(false, true, false);
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
    pbrick pbr = new Dirichlet_condition_brick(true, true, false, mf_mult);
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

    virtual void asm_real_tangent_terms(const model &md, size_type /*ib*/,
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
          scalar_type val(0);
          if (A) val = (mf_data ? (*A)[i] :  (*A)[i%s]);
          md.add_real_dof_constraint(vl[0], i, val);
        }
      }
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type /*ib*/,
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

    simplification_Dirichlet_condition_brick(void) {
      set_flags("Dirichlet with simplification brick",
                true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                true /* compute each time */, false /* has a Neumann term */);
    }
  };

  size_type add_Dirichlet_condition_with_simplification
  (model &md, const std::string &varname,
   size_type region, const std::string &dataname) {
    pbrick pbr = new simplification_Dirichlet_condition_brick();
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

  // Deprecated, should be re-written
  struct dirichlet_nitsche_nonlinear_term : public nonlinear_elem_term {
    // Option:
    // 1 : matrix term H^TH/gamma
    // 2 : matrix term -(D_u G(u,lambda)[w])^TH^TH
    // 3 : matrix term theta(g-Hu)^TH(D^2_uu G(u,lambda)[w,v])
    // 4 : rhs term (H^Tg)/gamma
    // 5 : rhs term H^T((g-Hu)/gamma + HG(u,lambda))
    // 6 : rhs term -theta(g)^TH(D_uG(u, lambda)[v])
    // 7 : rhs term theta(Hu-g)^TH(DG(u, lambda)[v])
    // 8 : matrix term theta(g-Hu)^TH(D^2_{u,lambda}G(u, lambda)[w,v])
    // 9 : matrix term -(D_lambda G(u,lambda)[w])^TH^TH

    dim_type N, qdim;
    size_type option;
    const model *md;
    const std::string *varname, *auxvarname;
    bool H_version, normal_component;
    scalar_type theta, gamma0;


    base_small_vector auxg, auxn, u, g, n;
    base_tensor tp;
    scalar_type gamma;
    base_vector coeff;
    base_matrix H, HTH, auxH;
    const mesh_fem *mf_u, *mf_lambda;
    const mesh_fem *mf_data;
    const mesh_fem *mf_H;

    base_vector U;
    const base_vector &HH, &G;

    mutable bgeot::multi_index sizes_;


    void adjust_tensor_size(void) {
      switch(option) {
      case 1:
        if (qdim > 1) { sizes_.resize(2); sizes_[0] = sizes_[1] = qdim; }
        else { sizes_.resize(1); sizes_[0] = 1; }
        break;
      case 2: case 9:
        if (qdim > 1)
          { sizes_.resize(2); sizes_[0] = 0; sizes_[1] = qdim; }
        else { sizes_.resize(1); sizes_[0] = 0; }
        break;
      case 3: case 8:
        sizes_.resize(3);
        sizes_[0] = sizes_[1] = 0; sizes_[2] = 1;
        break;
      case 4: case 5:
        sizes_.resize(1); sizes_[0] = qdim;
        break;
      case 6: case 7:
        sizes_.resize(2); sizes_[0] = 0; sizes_[1] = 1;
        break;
      }

      gmm::resize(u, qdim);
      gmm::resize(auxg, 1);
      gmm::resize(auxn, qdim);
      gmm::resize(g, normal_component ? 1 : qdim);
      gmm::resize(H, qdim, qdim); gmm::resize(HTH, qdim, qdim);
      gmm::resize(auxH, qdim, 1);
    }

    const bgeot::multi_index &sizes(size_type cv) const {
      if (cv != size_type(-1))
        switch(option) {
        case 2:
          sizes_[0] = short_type(mf_u->nb_basic_dof_of_element(cv));
          break;
        case 9:
          sizes_[0] = short_type(mf_lambda->nb_basic_dof_of_element(cv));
          break;
        case 3:
          sizes_[0] = sizes_[1] = short_type(mf_u->nb_basic_dof_of_element(cv));
          break;
        case 8:
          sizes_[0] = short_type(mf_u->nb_basic_dof_of_element(cv));
          sizes_[1] = short_type(mf_lambda->nb_basic_dof_of_element(cv));
          break;
        case 6: case 7:
          sizes_[0] = short_type(mf_u->nb_basic_dof_of_element(cv));
          break;
        }
      return sizes_;
    }

    dirichlet_nitsche_nonlinear_term
    (size_type option_, const model *md_, const std::string *varname_,
     const mesh_fem *mfu_, const model_real_plain_vector *U_,
     scalar_type theta_, scalar_type gamma0_, bool H_version_,
     bool normal_component_, const mesh_fem *mf_data_ = 0,
     const model_real_plain_vector *G_ = 0, const mesh_fem *mf_H_ = 0,
     const model_real_plain_vector *H_ = 0,
     const std::string *auxvarname_ = 0, const mesh_fem *mf_lambda_ = 0
     )
      : option(option_), md(md_), varname(varname_), auxvarname(auxvarname_),
        H_version(H_version_), normal_component(normal_component_),
        theta(theta_), gamma0(gamma0_), mf_u(mfu_), mf_lambda(mf_lambda_),
        mf_data(mf_data_), mf_H(mf_H_), HH(*H_), G(*G_) {

      N = mf_u->linked_mesh().dim();
      qdim = mf_u->get_qdim();
      adjust_tensor_size();

      if (U_) {
        gmm::resize(U, mf_u->nb_basic_dof());
        mf_u->extend_vector(*U_, U);
      }

      if (mf_data) GMM_ASSERT1(!(mf_data->is_reduced()),
                               "Reduced fem not allowed for data");
      if (mf_H) GMM_ASSERT1(!(mf_H->is_reduced()),
                            "Reduced fem not allowed for data");
    }

    void compute(fem_interpolation_context &ctx, bgeot::base_tensor &t) {

      dim_type i;
      // size_type cv = ctx.convex_num();

      switch (option) {
      case 1:
        for (i = 0; i < qdim*qdim; ++i) t[i] = HTH[i]/gamma;
        break;
      case 2:
        if (qdim == 1) {
          md->compute_Neumann_terms(2, *varname, *mf_u, U, ctx, n, t);
          t *= -scalar_type(1);
        } else {
          tp.adjust_sizes(sizes_);
          md->compute_Neumann_terms(2, *varname, *mf_u, U, ctx, n, tp);
          t.mat_reduction(tp, HTH, 1);
          t *= -scalar_type(1);
        }
        break;
      case 9:
        if (qdim == 1) {
          md->compute_auxilliary_Neumann_terms(2, *varname, *mf_u, U,
                                               *auxvarname, ctx, n, t);
          t *= -scalar_type(1);
        } else {
          tp.adjust_sizes(sizes_);
          md->compute_auxilliary_Neumann_terms(2, *varname, *mf_u, U,
                                               *auxvarname,ctx, n, tp);
          t.mat_reduction(tp, HTH, 1);
          t *= -scalar_type(1);
        }
        break;
      case 3:
        sizes_[2] = qdim;
        tp.adjust_sizes(sizes_);
        sizes_[2] = 1;
        md->compute_Neumann_terms(3, *varname, *mf_u, U, ctx, n, tp);
        gmm::mult(H, gmm::scaled(u, -theta), gmm::scaled(g, theta), auxn);
        gmm::mult(gmm::transposed(H), gmm::col_vector(auxn), auxH);
        t.mat_reduction(tp, auxH, 2);
        break;
      case 8:
        sizes_[2] = qdim;
        tp.adjust_sizes(sizes_);
        sizes_[2] = 1;
        md->compute_auxilliary_Neumann_terms(3, *varname,  *mf_u, U,
                                             *auxvarname, ctx, n, tp);
        gmm::mult(H, gmm::scaled(u, -theta), gmm::scaled(g, theta), auxn);
        gmm::mult(gmm::transposed(H), gmm::col_vector(auxn), auxH);
        t.mat_reduction(tp, auxH, 2);
        break;
      case 4:
        gmm::mult(gmm::transposed(H), g, t.as_vector());
        t /= gamma;
        break;
      case 5:
        gmm::mult(H, gmm::scaled(u, -scalar_type(1)), g, auxn);
        gmm::scale(auxn, scalar_type(1)/gamma);
        tp.adjust_sizes(sizes_);
        md->compute_Neumann_terms(1, *varname, *mf_u, U, ctx, n, tp);
        gmm::mult_add(H, tp.as_vector(), auxn);
        gmm::mult(gmm::transposed(H), auxn, t.as_vector());
        break;
      case 6:
        sizes_[1] = qdim;
        tp.adjust_sizes(sizes_);
        sizes_[1] = 1;
        md->compute_Neumann_terms(2, *varname, *mf_u, U, ctx, n, tp);
        gmm::copy(gmm::scaled(g, -theta), auxn);
        gmm::mult(gmm::transposed(H), gmm::col_vector(auxn), auxH);
        t.mat_reduction(tp, auxH, 1);
        break;
      case 7:
        sizes_[1] = qdim;
        tp.adjust_sizes(sizes_);
        sizes_[1] = 1;
        md->compute_Neumann_terms(2, *varname, *mf_u, U, ctx, n, tp);
        gmm::mult(H, gmm::scaled(u, theta), gmm::scaled(g, -theta), auxn);
        gmm::mult(gmm::transposed(H), gmm::col_vector(auxn), auxH);
        t.mat_reduction(tp, auxH, 1);
        break;
      }
    }

    void prepare(fem_interpolation_context& ctx, size_type nb) {
      size_type cv = ctx.convex_num();

      switch (nb) { // last is computed first
      case 1 : // mandatory. calculate [u], [n], [gamma], [HTH]
        n = bgeot::compute_normal(ctx, ctx.face_num());
        n /= gmm::vect_norm2(n);
        if (mf_u && gmm::vect_size(U)) {
          slice_vector_on_basic_dof_of_element(*mf_u, U, cv, coeff);
          // coeff.resize(mf_u->nb_basic_dof_of_element(cv));
          // gmm::copy(gmm::sub_vector(U, gmm::sub_index
          //           (mf_u->ind_basic_dof_of_element(cv))), coeff);
          ctx.pf()->interpolation(ctx, coeff, u, qdim);
        }
        if (normal_component) {
          GMM_ASSERT1(qdim == N, "dimensions mismatch");
          for (size_type i = 0; i < qdim; ++i)
            for (size_type j = 0; j < qdim; ++j)
              HTH(i,j) = H(i,j) = n[i]*n[j];
        }
        else if (!H_version) {
          gmm::copy(gmm::identity_matrix(), HTH);
          gmm::copy(gmm::identity_matrix(), H);
        }
        else {
          GMM_ASSERT1(&HH && gmm::vect_size(HH), "Need H in this case !");
          if (!mf_H) gmm::copy(HH, H.as_vector());
          gmm::clear(HTH);
          for (size_type i = 0; i < qdim; ++i)
            for (size_type j = 0; j < qdim; ++j)
              for (size_type k = 0; k < qdim; ++k)
                HTH(i,j) += H(k,i) * H(k,j);
        }
        if (!mf_data) {
          if (&G && gmm::vect_size(G))
            if (normal_component) gmm::copy(G, auxg); else gmm::copy(G, g);
          else
            if (normal_component) gmm::clear(auxg); else gmm::clear(g);
        }
        if (normal_component) gmm::copy(gmm::scaled(n, auxg[0]), g);
        // computation of h for gamma = gamma0*h
        scalar_type emax, emin; gmm::condition_number(ctx.K(),emax,emin);
        gamma = gamma0 * emax / sqrt(scalar_type(N));
        break;

      case 2 : // calculate [g]
        if (&G && gmm::vect_size(G)) {
          size_type ndof = mf_data->nb_basic_dof_of_element(cv);
          size_type qmult = qdim / mf_data->get_qdim();
          coeff.resize(ndof * qmult);
          mesh_fem::ind_dof_ct ct = mf_data->ind_basic_dof_of_element(cv);
          for (size_type i = 0; i < ndof; ++i)
            for (size_type j = 0; j < qmult; ++j)
              coeff[i*qmult+j] = G[ct[i]*qmult+j];
          if (normal_component)
            ctx.pf()->interpolation(ctx, coeff, auxg, 1);
          else
            ctx.pf()->interpolation(ctx, coeff, g, qdim);
        }
        break;

      case 3 :// calculate [H]
        if (&HH && gmm::vect_size(HH)) {
          size_type ndof = mf_H->nb_basic_dof_of_element(cv);
          size_type qmult = qdim*qdim / mf_H->get_qdim();
          coeff.resize(ndof * qmult);
          mesh_fem::ind_dof_ct ct = mf_H->ind_basic_dof_of_element(cv);
          for (size_type i = 0; i < ndof; ++i)
            for (size_type j = 0; j < qmult; ++j)
              coeff[i*qmult+j] = HH[ct[i]*qmult+j];
          ctx.pf()->interpolation(ctx, coeff, H.as_vector(),
                                  dim_type(qdim*qdim));
        }
        break;

      default : GMM_ASSERT1(false, "Invalid option");
      }
    }
  };

  void asm_Dirichlet_Nitsche_first_tangent_term
  (model_real_sparse_matrix &M, const mesh_im &mim, const model &md,
   const std::string &varname, const mesh_fem &mfu,
   const model_real_plain_vector *U,
   scalar_type theta, scalar_type gamma0, bool H_version,
   bool normal_component, const mesh_fem *mf_H,
   const model_real_plain_vector *H, const mesh_region &rg) {

    dirichlet_nitsche_nonlinear_term nterm(2, &md, &varname, &mfu, U, theta,
                                           gamma0, H_version, normal_component,
                                           0, 0, mf_H, H);

    getfem::generic_assembly assem;

    std::string Nlinfems = mf_H ? "#1,#1,#2" : "#1";

    if (mfu.get_qdim() > 1)
      assem.set("M(#1,#1)+=comp(vBase(#1).NonLin$1(#1,"+Nlinfems+"))(:,i,:,i);");
    else
      assem.set("M(#1,#1)+=comp(Base(#1).NonLin$1(#1,#1))(:,:);");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    if (mf_H) assem.push_mf(*mf_H);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(M);
    assem.assembly(rg);
  }

  void asm_Dirichlet_Nitsche_second_tangent_term
  (model_real_sparse_matrix &M, const mesh_im &mim, const mesh_fem &mfu,
   scalar_type theta, scalar_type gamma0, bool H_version,
   bool normal_component, const mesh_fem *mf_H,
   const model_real_plain_vector *H, const mesh_region &rg) {


    dirichlet_nitsche_nonlinear_term nterm(1, 0, 0, &mfu, 0, theta, gamma0,
                                           H_version, normal_component,
                                           0, 0, mf_H, H);

    getfem::generic_assembly assem;

    std::string Nlinfems = mf_H ? "#1,#1,#2" : "#1";

    if (mfu.get_qdim() > 1)
      assem.set("M(#1,#1)+=sym(comp(NonLin$1(#1,"+Nlinfems+").vBase(#1).vBase(#1))(i,j,:,i,:,j));");
    else
      assem.set("M(#1,#1)+=sym(comp(NonLin$1(#1,#1).Base(#1).Base(#1))(i,:,:));");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    if (mf_H) assem.push_mf(*mf_H);
    assem.push_nonlinear_term(&nterm);

    assem.push_mat(M);
    assem.assembly(rg);
  }


  void asm_Dirichlet_Nitsche_third_tangent_term
  (model_real_sparse_matrix &M, const mesh_im &mim, const model &md,
   const std::string &varname, const mesh_fem &mfu,
   const model_real_plain_vector *U, scalar_type theta, scalar_type gamma0,
   bool H_version, bool normal_component,
   const mesh_fem *mf_H, const model_real_plain_vector *H,
   const mesh_fem *mf_data, const model_real_plain_vector *G,
   const mesh_region &rg) {

    dirichlet_nitsche_nonlinear_term nterm(3, &md, &varname, &mfu, U, theta,
                                           gamma0, H_version, normal_component,
                                           mf_data, G, mf_H, H);

    getfem::generic_assembly assem;

    std::string Nlinfems = "#1";
    if (mf_H && mf_data) Nlinfems = "#1,#2,#3";
    else if (mf_H) Nlinfems = "#1,#1,#2";
    else if (mf_data) Nlinfems = "#1,#2";

    assem.set("M(#1,#1)+=comp(NonLin$1(#1,"+Nlinfems+"))(:,:,i);");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    if (mf_data) assem.push_mf(*mf_data);
    if (mf_H) assem.push_mf(*mf_H);
    assem.push_nonlinear_term(&nterm);

    assem.push_mat(M);
    assem.assembly(rg);
  }


  void asm_Dirichlet_Nitsche_fourth_tangent_term
  (model_real_sparse_matrix &M, const mesh_im &mim, const model &md,
   const std::string &varname, const mesh_fem &mfu,
   const model_real_plain_vector *U,
   const std::string &auxvarname, const mesh_fem &mf_lambda,
   scalar_type theta,
   scalar_type gamma0, bool H_version, bool normal_component,
   const mesh_fem *mf_H, const model_real_plain_vector *H,
   const mesh_fem *mf_data, const model_real_plain_vector *G,
   const mesh_region &rg) {

    dirichlet_nitsche_nonlinear_term nterm(8, &md, &varname, &mfu, U, theta,
                                           gamma0, H_version, normal_component,
                                           mf_data, G, mf_H, H, &auxvarname,
                                           &mf_lambda);
    getfem::generic_assembly assem;

    std::string Nlinfems = "#1", lambdafem = "#2";
    if (mf_H && mf_data) { Nlinfems = "#1,#2,#3"; lambdafem = "#4"; }
    else if (mf_H) { Nlinfems = "#1,#1,#2"; lambdafem = "#3"; }
    else if (mf_data) { Nlinfems = "#1,#2"; lambdafem = "#3"; }

    assem.set("M(#1,"+lambdafem+")+=comp(NonLin$1(#1,"+Nlinfems+"))(:,:,i);");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    if (mf_data) assem.push_mf(*mf_data);
    if (mf_H) assem.push_mf(*mf_H);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm);

    assem.push_mat(M);
    assem.assembly(rg);
  }

  void asm_Dirichlet_Nitsche_fifth_tangent_term
  (model_real_sparse_matrix &M, const mesh_im &mim, const model &md,
   const std::string &varname, const mesh_fem &mfu,
   const model_real_plain_vector *U,
   const std::string &auxvarname, const mesh_fem &mf_lambda,
   scalar_type theta, scalar_type gamma0, bool H_version,
   bool normal_component, const mesh_fem *mf_H,
   const model_real_plain_vector *H, const mesh_region &rg) {

    dirichlet_nitsche_nonlinear_term nterm(9, &md, &varname, &mfu, U, theta,
                                           gamma0, H_version, normal_component,
                                           0, 0, mf_H, H, &auxvarname,
                                           &mf_lambda);
    getfem::generic_assembly assem;

    std::string Nlinfems = mf_H ? "#1,#1,#2" : "#1";
    std::string lambdafem = mf_H ? "#3" : "#2";

    if (mfu.get_qdim() > 1)
      assem.set("M(#1,"+lambdafem+")+=comp(vBase(#1).NonLin$1(#1,"+Nlinfems+"))(:,i,:,i);");
    else
      assem.set("M(#1,"+lambdafem+")+=comp(Base(#1).NonLin$1(#1,#1))(:,:);");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    if (mf_H) assem.push_mf(*mf_H);
    assem.push_mf(mf_lambda);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(M);
    assem.assembly(rg);
  }


  void asm_Dirichlet_Nitsche_first_rhs_term
  (model_real_plain_vector &V, const mesh_im &mim, const model &md,
   const std::string &varname, const mesh_fem &mfu,
   const model_real_plain_vector *U, scalar_type theta, scalar_type gamma0,
   bool H_version, bool normal_component,
   const mesh_fem *mf_H, const model_real_plain_vector *H,
   const mesh_fem *mf_data, const model_real_plain_vector *G, bool is_linear,
   const mesh_region &rg) {

    dirichlet_nitsche_nonlinear_term nterm(is_linear ? 4:5, &md, &varname,
                                           &mfu, U, theta, gamma0, H_version,
                                           normal_component, mf_data,
                                           G, mf_H, H);

    getfem::generic_assembly assem;
    std::string Nlinfems = "#1";
    if (mf_H && mf_data) Nlinfems = "#1,#2,#3";
    else if (mf_H) Nlinfems = "#1,#1,#2";
    else if (mf_data) Nlinfems = "#1,#2";

    if (mfu.get_qdim() > 1)
      assem.set("V(#1)+=comp(NonLin$1(#1,"+Nlinfems+").vBase(#1))(i,:,i);");
    else
      assem.set("V(#1)+=comp(NonLin$1(#1,"+Nlinfems+").Base(#1))(i,:);");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    if (mf_data) assem.push_mf(*mf_data);
    if (mf_H) assem.push_mf(*mf_H);
    assem.push_nonlinear_term(&nterm);

    assem.push_vec(V);
    assem.assembly(rg);
  }

  void asm_Dirichlet_Nitsche_second_rhs_term
  (model_real_plain_vector &V, const mesh_im &mim, const model &md,
   const std::string &varname, const mesh_fem &mfu,
   const model_real_plain_vector *U, scalar_type theta, scalar_type gamma0,
   bool H_version, bool normal_component,
   const mesh_fem *mf_H, const model_real_plain_vector *H,
   const mesh_fem *mf_data, const model_real_plain_vector *G, bool is_linear,
   const mesh_region &rg) {

    dirichlet_nitsche_nonlinear_term nterm(is_linear ? 6:7, &md, &varname,
                                           &mfu, U, theta, gamma0, H_version,
                                           normal_component, mf_data,
                                           G, mf_H, H);

    getfem::generic_assembly assem;

    std::string Nlinfems = "#1";
    if (mf_H && mf_data) Nlinfems = "#1,#2,#3";
    else if (mf_H) Nlinfems = "#1,#1,#2";
    else if (mf_data) Nlinfems = "#1,#2";

    assem.set("V(#1)+=comp(NonLin$1(#1,"+Nlinfems+"))(:,i);");
    assem.push_mi(mim);
    assem.push_mf(mfu);
    if (mf_data) assem.push_mf(*mf_data);
    if (mf_H) assem.push_mf(*mf_H);
    assem.push_nonlinear_term(&nterm);

    assem.push_vec(V);
    assem.assembly(rg);
  }


  struct Nitsche_Dirichlet_condition_brick : public virtual_brick {

    bool H_version; // The version hu = r for vector fields.
    bool normal_component; // Dirichlet on normal component for vector field.
    bool linear_version;
    scalar_type theta;

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(vecl.size() == vl.size() && matl.size() == vl.size(),
                  "Wrong number of terms for Dirichlet condition brick");
      GMM_ASSERT1(mims.size() == 1,
                  "Dirichlet condition brick need one and only one mesh_im");
      GMM_ASSERT1(vl.size() >= 1 && dl.size() >= 1 && dl.size() <= 3,
                  "Wrong number of variables for Dirichlet condition brick");


      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const model_real_plain_vector *U = &(md.real_variable(vl[0]));
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector *G = 0, *H = 0;
      const mesh_fem *mf_data = 0, *mf_H = 0;
      bool recompute_matrix = (!is_linear() && (version & model::BUILD_MATRIX))
        || (is_linear() && (!((version & model::BUILD_ON_DATA_CHANGE) != 0)
                            || md.is_var_newer_than_brick(dl[0], ib)
                            || md.is_var_newer_than_brick(dl[1], ib)));

      GMM_ASSERT1(gmm::vect_size(md.real_variable(dl[0])) == 1,
                  "Parameter gamma0 for Nitsche's method should be a scalar");
      scalar_type gamma0 = md.real_variable(dl[0])[0];

      size_type s = 0, ind = 1;
      if (dl.size() > 1 + (H_version ? 1:0)) {
        ++ind;
        G = &(md.real_variable(dl[1]));
        mf_data = md.pmesh_fem_of_variable(dl[1]);
        s = gmm::vect_size(*G);
        if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();
        size_type ss = s * ((normal_component) ? mf_u.linked_mesh().dim() : 1);
        GMM_ASSERT1(mf_u.get_qdim() == ss, dl[1] << ": bad format of "
                    "Dirichlet data. Detected dimension is " << ss
                    << " should be " << size_type(mf_u.get_qdim()));
      }

      if (H_version) {
        GMM_ASSERT1(H_version,
                    "Wrong number of data for Dirichlet condition brick");
        H = &(md.real_variable(dl[ind]));
        mf_H = md.pmesh_fem_of_variable(dl[ind]);
        s = gmm::vect_size(*H);
        if (mf_H) {
          s = s * mf_H->get_qdim() / mf_H->nb_dof();
          // GMM_ASSERT1(mf_H->get_qdim() == 1,  "Implemented only for mf_H "
          //             "a scalar finite element method");
        }
        GMM_ASSERT1(s = gmm::sqr(mf_u.get_qdim()),
                    dl[ind] << ": bad format of Dirichlet data. "
                    "Detected dimension is " << s << " should be "
                    << size_type(gmm::sqr(mf_u.get_qdim())));
      }

      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      // Test Neumann term consistency if some computation are needed
      if (recompute_matrix || (!linear_version && (version & model::BUILD_RHS))
          || (linear_version && G)) {
        size_type ifb = md.check_Neumann_terms_consistency(vl[0]);
        GMM_ASSERT1(ifb == size_type(-1),
                    "Impossible to build Nitsche's terms for Dirichlet "
                    " condition. At least '"
                    << md.brick_pointer(ifb)->brick_name() << "' is declared "
                    "after Nitsche's brick or do not declare a Neumann term.");
      }

      if (recompute_matrix) {

        gmm::clear(matl[0]);

        GMM_TRACE2("Assembly of Nitsche's tangent terms "
                   "for Dirichlet condition");
        asm_Dirichlet_Nitsche_first_tangent_term
          (matl[0], mim, md, vl[0], mf_u, U, theta, gamma0, H_version,
           normal_component, mf_H, H, rg);

        if (theta != scalar_type(0)) {
          model_real_sparse_matrix B(matl[0]);
          gmm::scale(B, theta);
          gmm::add(gmm::transposed(B), matl[0]);
        }

        asm_Dirichlet_Nitsche_second_tangent_term
          (matl[0], mim, mf_u, theta, gamma0, H_version, normal_component,
           mf_H, H, rg);

        if (theta != scalar_type(0) && !linear_version) {
          asm_Dirichlet_Nitsche_third_tangent_term
            (matl[0], mim, md, vl[0], mf_u, U, theta, gamma0, H_version,
             normal_component, mf_H, H, mf_data, G, rg);
        }

        for (size_type i = 1; i < vl.size(); ++i) { // Auxilliary variables
          gmm::clear(matl[i]);
          if (theta != scalar_type(0) && !linear_version)
            asm_Dirichlet_Nitsche_fourth_tangent_term
              (matl[i], mim, md, vl[0], mf_u, U, vl[i],
               md.mesh_fem_of_variable(vl[i]), theta, gamma0,
               H_version, normal_component, mf_H, H, mf_data, G, rg);
          asm_Dirichlet_Nitsche_fifth_tangent_term
            (matl[i], mim, md, vl[0], mf_u, U, vl[i],
             md.mesh_fem_of_variable(vl[i]), theta, gamma0,
             H_version, normal_component, mf_H, H, rg);
        }
      }

      if ((!linear_version && (version & model::BUILD_RHS))
          || (linear_version && G)) {

        GMM_TRACE2("Assembly of Nitsche's source terms "
                   "for Dirichlet condition");
        asm_Dirichlet_Nitsche_first_rhs_term
          (vecl[0], mim, md, vl[0], mf_u, U, theta, gamma0, H_version,
           normal_component, mf_H, H, mf_data, G, linear_version, rg);

        if (theta != scalar_type(0)) {
          asm_Dirichlet_Nitsche_second_rhs_term
            (vecl[0], mim, md, vl[0], mf_u, U, theta, gamma0, H_version,
             normal_component, mf_H, H, mf_data, G, linear_version, rg);
        }
      }
    }


    Nitsche_Dirichlet_condition_brick(bool H_version_,
                                      bool normal_component_,
                                      bool is_linear_,
                                      scalar_type theta_) {
      H_version = H_version_;
      normal_component = normal_component_;
      // linear_version = false;
      linear_version = is_linear_;
      theta = theta_;
      GMM_ASSERT1(!(H_version && normal_component), "Bad Dirichlet version");
      set_flags(is_linear_ ? "Dirichlet with Nitsche's method linear brick"
                : "Dirichlet with Nitsche's method nonlinear brick",
                linear_version /* is linear*/,
                (theta==scalar_type(1)) /* is symmetric */,
                (theta==scalar_type(1)) /* is coercive */,
                true /* is real */, false /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
    }
  };


  size_type add_Dirichlet_condition_with_Nitsche_method
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &gamma0name, size_type region, scalar_type theta,
   const std::string &dataname) {

    pbrick pbr = new Nitsche_Dirichlet_condition_brick
      (false, false, md.check_Neumann_terms_linearity(varname), theta);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname,
                                         theta == scalar_type(1)));
    model::varnamelist vl(1, varname);

    std::vector<std::string> aux_vars;
    md.auxilliary_variables_of_Neumann_terms(varname, aux_vars);
    for (size_type i = 0; i < aux_vars.size(); ++i) {
      vl.push_back(aux_vars[i]);
      tl.push_back(model::term_description(varname, aux_vars[i], false));
    }

    model::varnamelist dl;
    dl.push_back(gamma0name);
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }


  size_type add_normal_Dirichlet_condition_with_Nitsche_method
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &gamma0name, size_type region, scalar_type theta,
   const std::string &dataname) {
    pbrick pbr = new Nitsche_Dirichlet_condition_brick
      (false, true, md.check_Neumann_terms_linearity(varname), theta);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname,
                                         theta == scalar_type(1)));
    model::varnamelist vl(1, varname);

    std::vector<std::string> aux_vars;
    md.auxilliary_variables_of_Neumann_terms(varname, aux_vars);
    for (size_type i = 0; i < aux_vars.size(); ++i) {
      vl.push_back(aux_vars[i]);
      tl.push_back(model::term_description(varname, aux_vars[i], false));
    }

    model::varnamelist dl;
    dl.push_back(gamma0name);
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

  size_type add_generalized_Dirichlet_condition_with_Nitsche_method
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &gamma0name, size_type region, scalar_type theta,
   const std::string &dataname, const std::string &Hname) {
    pbrick pbr = new Nitsche_Dirichlet_condition_brick
      (true, false, md.check_Neumann_terms_linearity(varname), theta);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname,
                                         theta == scalar_type(1)));
    model::varnamelist vl(1, varname);

    std::vector<std::string> aux_vars;
    md.auxilliary_variables_of_Neumann_terms(varname, aux_vars);
    for (size_type i = 0; i < aux_vars.size(); ++i) {
      vl.push_back(aux_vars[i]);
      tl.push_back(model::term_description(varname, aux_vars[i], false));
    }

    model::varnamelist dl;
    dl.push_back(gamma0name);
    dl.push_back(dataname);
    dl.push_back(Hname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
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

    pointwise_constraints_brick(bool penalized) {
      set_flags(penalized ? "Pointwise cosntraints with penalization brick"
                          : "Pointwise cosntraints with multipliers brick",
                true /* is linear*/,
                true /* is symmetric */, penalized /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
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
    pbrick pbr = new pointwise_constraints_brick(true);
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
    pbrick pbr = new  pointwise_constraints_brick(false);
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

    Helmholtz_brick(void) {
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
      pbrick pbr = new Helmholtz_brick;
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

       size_type ib = add_linear_generic_assembly_brick
         (md, mim, expr, region, true, true, "Helmholtz", true);
       if (ib == size_type(-1))
         ib = add_nonlinear_generic_assembly_brick
           (md, mim, expr, region, false, false, "Helmholtz (nonlinear)");
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

    Fourier_Robin_brick(void) {
      set_flags("Fourier Robin condition", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
    }

  };

  size_type add_Fourier_Robin_brick(model &md, const mesh_im &mim,
                                    const std::string &varname,
                                    const std::string &dataexpr,
                                    size_type region) {
    if (md.is_complex()) {
      pbrick pbr = new Fourier_Robin_brick;
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      return md.add_brick(pbr, model::varnamelist(1, varname),
                          model::varnamelist(1, dataexpr), tl,
                          model::mimlist(1, &mim), region);
    } else {
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      std::string expr = "(("+dataexpr+")*"+varname+")."+test_varname;
      size_type ib = add_linear_generic_assembly_brick
        (md, mim, expr, region, true, true, "Fourier-Robin", true);
      if (ib == size_type(-1))
        ib = add_nonlinear_generic_assembly_brick
          (md, mim, expr, region, false, false, "Fourier-Robin (nonlinear)");
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

        if (penalized) {
          COEFF = &(md.real_variable(dl[0]));
          GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                      "Data for coefficient should be a scalar");

          gmm::mult(gmm::transposed(rB),
                    gmm::scaled(rL, gmm::abs((*COEFF)[0])), vecl[0]);
          gmm::mult(gmm::transposed(rB),
                    gmm::scaled(rB, gmm::abs((*COEFF)[0])), matl[0]);
        } else {
          gmm::copy(rL, vecl[0]);
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

        if (penalized) {
          COEFF = &(md.complex_variable(dl[0]));
          GMM_ASSERT1(gmm::vect_size(*COEFF) == 1,
                      "Data for coefficient should be a scalar");

          gmm::mult(gmm::transposed(cB),
                    gmm::scaled(cL, gmm::abs((*COEFF)[0])), vecl[0]);
          gmm::mult(gmm::transposed(cB),
                    gmm::scaled(cB, gmm::abs((*COEFF)[0])), matl[0]);
        } else {
          gmm::copy(cL, vecl[0]);
          gmm::copy(cB, matl[0]);
        }
      }
    }

    constraint_brick(bool penalized) {
      set_flags(penalized ? "Constraint with penalization brick"
                          : "Constraint with multipliers brick",
                true /* is linear*/,
                true /* is symmetric */, penalized /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
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
    return p->cL;
  }

  size_type add_constraint_with_penalization
  (model &md, const std::string &varname, scalar_type penalisation_coeff) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = new constraint_brick(true);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(), size_type(-1));
  }

  size_type add_constraint_with_multipliers
  (model &md, const std::string &varname, const std::string &multname) {
    pbrick pbr = new constraint_brick(false);
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

    explicit_matrix_brick(bool symmetric_, bool coercive_) {
      set_flags("Explicit matrix brick",
                true /* is linear*/,
                symmetric_ /* is symmetric */, coercive_ /* is coercive */,
                true /* is real */, true /* is complex */,
                true /* is to be computed each time */,
                false /* has a Neumann term */);
    }
  };

  size_type add_explicit_matrix
  (model &md, const std::string &varname1, const std::string &varname2,
   bool issymmetric, bool iscoercive) {
    pbrick pbr = new explicit_matrix_brick(issymmetric, iscoercive);
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

    explicit_rhs_brick(void) {
      set_flags("Explicit rhs brick",
                true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                true /* is to be computed each time */,
                false /* has a Neumann term */);
    }

  };

  size_type add_explicit_rhs
  (model &md, const std::string &varname) {
    pbrick pbr = new explicit_rhs_brick();
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

  struct iso_lin_elasticity_Neumann_elem_term : public Neumann_elem_term {

    const mesh_fem *mf_lambda;
    const model_real_plain_vector *lambda;
    const mesh_fem *mf_mu;
    const model_real_plain_vector *mu;

    mutable fem_interpolation_context ctx_mu;
    mutable base_vector coeff, val;
    mutable base_matrix grad, E, G;

    void compute_Neumann_term
    (int version, const mesh_fem &mfvar, const model_real_plain_vector &var,
     fem_interpolation_context& ctx, base_small_vector &n,
     base_tensor &output, size_type /*auxilliary_ind*/ = 0) const {

      if (version == 3) return;  // No contribution because the term is linear

      dim_type qdim = mfvar.linked_mesh().dim();
      gmm::resize(grad, qdim, qdim);
      gmm::resize(E, qdim, qdim);
      gmm::resize(val, 1);
      size_type cv = ctx.convex_num();
      scalar_type val_lambda = scalar_type(0), val_mu = scalar_type(0);

      if (mf_mu) {
        GMM_ASSERT1(!(mf_mu->is_reduced()),
                    "Sorry, to be adapted for reduced mesh fems");

        if (!(ctx_mu.have_pf()) || ctx_mu.convex_num() != cv
            || (ctx_mu.have_pfp() != ctx.have_pfp())
            || (ctx_mu.have_pfp()
                && (&(ctx.pfp()->get_point_tab())
                    != &(ctx_mu.pfp()->get_point_tab())))) {

          bgeot::vectors_to_base_matrix
            (G, mf_mu->linked_mesh().points_of_convex(cv));

          pfem_precomp pfp = fem_precomp(mf_mu->fem_of_element(cv),
                                         &(ctx.pfp()->get_point_tab()), 0);

          if (ctx.have_pfp())
            ctx_mu = fem_interpolation_context
              (mf_mu->linked_mesh().trans_of_convex(cv), pfp, ctx.ii(),
               G, cv, ctx.face_num());
          else
            ctx_mu = fem_interpolation_context
              (mf_mu->linked_mesh().trans_of_convex(cv),
               mf_mu->fem_of_element(cv), ctx.xref(), G, cv, ctx.face_num());

        } else {
          if (ctx.have_pfp())  ctx_mu.set_ii(ctx.ii());
          else ctx_mu.set_xref(ctx.xref());
        }
        slice_vector_on_basic_dof_of_element(*mf_mu, *mu, cv, coeff);
        // coeff.resize(mf_mu->nb_basic_dof_of_element(cv));
        // gmm::copy(gmm::sub_vector(*mu, gmm::sub_index
        //                      (mf_mu->ind_basic_dof_of_element(cv))), coeff);
        ctx_mu.pf()->interpolation(ctx_mu, coeff, val, 1);
        val_mu = val[0];
        slice_vector_on_basic_dof_of_element(*mf_mu, *lambda, cv, coeff);
        // gmm::copy(gmm::sub_vector(*lambda, gmm::sub_index
        //                      (mf_mu->ind_basic_dof_of_element(cv))), coeff);
        ctx_mu.pf()->interpolation(ctx_mu, coeff, val, 1);
        val_mu = val[0];
      } else {
        val_lambda = (*lambda)[0]; val_mu = (*mu)[0];
      }

      switch (version) {
      case 1:
        slice_vector_on_basic_dof_of_element(mfvar, var, cv, coeff);
        // coeff.resize(mfvar.nb_basic_dof_of_element(cv));
        // gmm::copy(gmm::sub_vector(var, gmm::sub_index
        //                      (mfvar.ind_basic_dof_of_element(cv))), coeff);
        ctx.pf()->interpolation_grad(ctx, coeff, grad, qdim);
        gmm::copy(gmm::identity_matrix(), E);
        gmm::scale(E, val_lambda * gmm::mat_trace(grad));
        gmm::add(gmm::scaled(grad, val_mu), E);
        gmm::add(gmm::scaled(gmm::transposed(grad), val_mu), E);
        gmm::mult_add(E, n, output.as_vector());
        break;
      case 2:
        {
          base_tensor t;
          dim_type tdim = ctx.pf()->target_dim(), qmult = qdim / tdim;
          size_type ndof = ctx.pf()->nb_dof(cv);
          // The return tensor is t(i,j,k) with 0<=i<ndof, 0<=j<target_dim,
          // 0<=k<dim. In order to iterate on the tensor values, i should
          // increment the faster, then j, then k.
          // If target_dim == qdim, grad(phi_i)(j,k) = t(i,j,k)
          // If target_dim == 1, grad(phi_i * e_l)(l,k) = t(i,1,k)
          // General case, psi_{i*qmult+l} = phi_i * e_l  and
          //    grad(psi_{i*qmult+l})(j+tdim*l,k) = t(i,j,k)
          ctx.pf()->real_grad_base_value(ctx, t);

          for (size_type l = 0; l < qmult; ++l) {
            for (size_type p = 0; p < qdim; ++p) {
              base_tensor::const_iterator it = t.begin();
              for (size_type k = 0; k < qdim; ++k)
                for (size_type j = 0; j < tdim; ++j)
                  for (size_type i = 0; i < ndof; ++i, ++it) {
                    size_type jj = j + tdim*l;
                    if (k == jj) output(i*qmult+l, p) += val_lambda*(*it)*n[p];
                    if (p == jj) output(i*qmult+l, p) += val_mu*(*it)*n[k];
                    if (k == p) output(i*qmult+l, p) += val_mu*(*it)*n[jj];
                  }
              GMM_ASSERT1(it ==  t.end(), "Internal error");
            }
          }
        }
        break;
      }

    }

    iso_lin_elasticity_Neumann_elem_term
    (const mesh_fem *mf_lambda_,
     const model_real_plain_vector *lambda_,
     const mesh_fem *mf_mu_, const model_real_plain_vector *mu_) :
      mf_lambda(mf_lambda_), lambda(lambda_), mf_mu(mf_mu_), mu(mu_) {
      GMM_ASSERT1(mf_lambda == mf_mu,
                  "The two coefficients should be described on the same "
                  "finite element method.");
    }

  };

  struct iso_lin_elasticity_brick : public virtual_brick {

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
                  "isotropic linearized elasticity brick has one and only "
                  "one term");
      GMM_ASSERT1(mims.size() == 1,
                  "isotropic linearized elasticity brick need one and only "
                  "one mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() >= 2 && dl.size() <= 3,
                  "Wrong number of variables for isotropic linearized "
                  "elasticity brick");

      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
        || md.is_var_newer_than_brick(dl[0], ib)
        || md.is_var_newer_than_brick(dl[1], ib);

      if (recompute_matrix) {

        const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
        const mesh &m = mf_u.linked_mesh();
        size_type N = m.dim(), Q = mf_u.get_qdim();
        GMM_ASSERT1(Q == N, "isotropic linearized elasticity brick is only "
                    "for vector field of the same dimension as the mesh");
        const mesh_im &mim = *mims[0];
        mesh_region rg(region);
        m.intersect_with_mpi_region(rg);
        const mesh_fem *mf_lambda = md.pmesh_fem_of_variable(dl[0]);
        const model_real_plain_vector *lambda = &(md.real_variable(dl[0]));
        const mesh_fem *mf_mu = md.pmesh_fem_of_variable(dl[1]);
        const model_real_plain_vector *mu = &(md.real_variable(dl[1]));

        size_type sl = gmm::vect_size(*lambda);
        if (mf_lambda) sl = sl * mf_lambda->get_qdim() / mf_lambda->nb_dof();
        size_type sm = gmm::vect_size(*mu);
        if (mf_mu) sm = sm * mf_mu->get_qdim() / mf_mu->nb_dof();

        GMM_ASSERT1(sl == 1 && sm == 1, "Bad format of isotropic linearized "
                    "elasticity brick coefficients");
        GMM_ASSERT1(mf_lambda == mf_mu,
                    "The two coefficients should be described on the same "
                    "finite element method.");

        GMM_TRACE2("Stiffness matrix assembly for isotropic linearized "
                   "elasticity");
        gmm::clear(matl[0]);
        if (mf_lambda)
          asm_stiffness_matrix_for_linear_elasticity
            (matl[0], mim, mf_u, *mf_lambda, *lambda, *mu, rg);
        else
          asm_stiffness_matrix_for_homogeneous_linear_elasticity
            (matl[0], mim, mf_u, *lambda, *mu, rg);

      }


      if  (dl.size() == 3) { // Pre-constraints given by an "initial"
        // displacement u0. Means that the computed displacement will be u - u0
        // The displacement u0 should be discribed on the same fem as the
        // variable.
        gmm::mult(matl[0],
                  gmm::scaled(md.real_variable(dl[2]), scalar_type(-1)),
                  vecl[0]);
      }
    }

    virtual void real_post_assembly_in_serial(const model &md, size_type ib,
                                              const model::varnamelist &vl,
                                              const model::varnamelist &dl,
                                              const model::mimlist &/*mims*/,
                                              model::real_matlist &/*matl*/,
                                              model::real_veclist &,
                                              model::real_veclist &,
                                              size_type /*region*/,
                                              build_version version) const
    {
      bool recompute_matrix = !((version & model::BUILD_ON_DATA_CHANGE) != 0)
        || md.is_var_newer_than_brick(dl[0], ib)
        || md.is_var_newer_than_brick(dl[1], ib);

      if (recompute_matrix)
      {
          const mesh_fem *mf_lambda = md.pmesh_fem_of_variable(dl[0]);
          const model_real_plain_vector *lambda = &(md.real_variable(dl[0]));
          const mesh_fem *mf_mu = md.pmesh_fem_of_variable(dl[1]);
          const model_real_plain_vector *mu = &(md.real_variable(dl[1]));

          pNeumann_elem_term pNt = new iso_lin_elasticity_Neumann_elem_term
            (mf_lambda, lambda, mf_mu, mu);
          md.add_Neumann_term(pNt, vl[0], ib);
      }
    }


    iso_lin_elasticity_brick(void) {
      set_flags("isotropic linearized elasticity", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };


  struct iso_lin_elasticity_new_brick : public virtual_brick {

    std::string expr, dataname3;

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
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
        size_type nbgdof = md.nb_dof();
        ga_workspace workspace(md, true);
        workspace.add_expression(expr, *(mims[0]), region);
        GMM_TRACE2(name << ": generic matrix assembly");
        model::varnamelist vlmd; md.variable_list(vlmd);
        for (size_type i = 0; i < vlmd.size(); ++i)
          if (md.is_disabled_variable(vlmd[i]))
            nbgdof = std::max(nbgdof, 
                              workspace.interval_of_variable(vlmd[i]).last());
        model_real_sparse_matrix R(nbgdof, nbgdof); 
        workspace.set_assembled_matrix(R);
        workspace.assembly(2);
        scalar_type alpha = scalar_type(1)
          / (workspace.factor_of_variable(vl[0]));
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
        md.add_external_load(ib, gmm::vect_norm1(vecl[0]));
      }

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
#if 0 // Old brick
      pbrick pbr = new iso_lin_elasticity_brick;
      model::termlist tl;
      tl.push_back(model::term_description(varname, varname, true));
      model::varnamelist dl(1, dataexpr1);
      dl.push_back(dataexpr2);
      if (dataname3.size()) dl.push_back(dataname3);
      return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
                          model::mimlist(1, &mim), region);
#else // New brick with high-level generic assembly
      std::string test_varname
        = "Test_" + sup_previous_and_dot_to_varname(varname);
      
      std::string expr1 = "("+dataexpr1+")*(Div_"+varname+"-Div_"+dataname3
        +")*Div_"+test_varname+"+("+dataexpr2+")*(Grad_"+varname+"+Grad_"
        +varname+"'-Grad_"+dataname3+"-Grad_"+dataname3+"'):Grad_"
        +test_varname;
      std::string expr2 = "("+dataexpr1+")*Div_"+varname+"*Div_"+test_varname
        + "+("+dataexpr2+")*(Grad_"+varname+"+Grad_"+varname+"'):Grad_"
        + test_varname;
      
      ga_workspace workspace(md, true);
      workspace.add_expression(expr2, mim, region);
      model::varnamelist vl, vl_test1, vl_test2, dl;
      bool is_lin = workspace.used_variables(vl, vl_test1, vl_test2, dl, 2);
      
      if (is_lin) {
        pbrick pbr = new iso_lin_elasticity_new_brick(expr2, dataname3);
        model::termlist tl;
        tl.push_back(model::term_description(varname,
                     sup_previous_and_dot_to_varname(varname), true));
        if (dataname3.size()) dl.push_back(dataname3);
        return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
      } else {
        return add_nonlinear_generic_assembly_brick
          (md, mim, dataname3.size() ? expr1 : expr2, region, false, false,
           "Linearized isotropic elasticity (with nonlinear dependance)");
      }
#endif
  }

  void compute_isotropic_linearized_Von_Mises_or_Tresca
  (model &md, const std::string &varname, const std::string &dataname_lambda,
   const std::string &dataname_mu, const mesh_fem &mf_vm,
   model_real_plain_vector &VM, bool tresca) {

    const mesh_fem &mf_u = md.mesh_fem_of_variable(varname);
    const mesh_fem *mf_lambda = md.pmesh_fem_of_variable(dataname_lambda);
    const model_real_plain_vector *lambda=&(md.real_variable(dataname_lambda));
    const mesh_fem *mf_mu = md.pmesh_fem_of_variable(dataname_mu);
    const model_real_plain_vector *mu = &(md.real_variable(dataname_mu));

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
  }

  // ----------------------------------------------------------------------
  //
  // linearized incompressibility brick  (div u = 0)
  //
  // ----------------------------------------------------------------------

  struct lin_incomp_Neumann_elem_term : public Neumann_elem_term {

    const gmm::uint64_type &var_vnum;
    const mesh_fem *mf_p;
    const model_real_plain_vector *org_P;
    mutable model_real_plain_vector P;
    mutable gmm::uint64_type vnum;

    mutable fem_interpolation_context ctx_p;
    mutable base_vector coeff, val;
    mutable base_matrix G;

    void compute_Neumann_term
    (int version, const mesh_fem &mfvar,
     const model_real_plain_vector &/* var */,
     fem_interpolation_context& ctx, base_small_vector &n,
     base_tensor &output, size_type auxilliary_ind = 0) const {

      if (version == 3) return;  // No contribution because the term is linear
      if (version == 2 && auxilliary_ind == 0) return;

      dim_type qdim = mfvar.linked_mesh().dim();
      size_type cv = ctx.convex_num();

      if (vnum != var_vnum || !(ctx_p.have_pf()) || ctx_p.convex_num() != cv
          || (ctx_p.have_pfp() != ctx.have_pfp())
          || (ctx_p.have_pfp()
              && (&(ctx.pfp()->get_point_tab())
                  != &(ctx_p.pfp()->get_point_tab())))) {

        if (vnum != var_vnum) {
          gmm::resize(P, mf_p->nb_basic_dof());
          mf_p->extend_vector(*org_P, P);
          vnum = var_vnum;
        }

        bgeot::vectors_to_base_matrix
          (G, mf_p->linked_mesh().points_of_convex(cv));

        if (ctx.have_pfp()) {
          pfem_precomp pfp = fem_precomp(mf_p->fem_of_element(cv),
                                         &(ctx.pfp()->get_point_tab()), 0);
          ctx_p = fem_interpolation_context
            (mf_p->linked_mesh().trans_of_convex(cv), pfp, ctx.ii(),
             G, cv, ctx.face_num());
        } else
          ctx_p = fem_interpolation_context
            (mf_p->linked_mesh().trans_of_convex(cv),
             mf_p->fem_of_element(cv), ctx.xref(), G, cv, ctx.face_num());
      } else {
        if (ctx_p.have_pfp()) ctx_p.set_ii(ctx.ii());
        else ctx_p.set_xref(ctx.xref());

      }

      switch (version) {
      case 1:
        slice_vector_on_basic_dof_of_element(*mf_p, P, cv, coeff);
        // coeff.resize(mf_p->nb_basic_dof_of_element(cv));
        // gmm::copy(gmm::sub_vector(P, gmm::sub_index
        //                      (mf_p->ind_basic_dof_of_element(cv))), coeff);
        ctx_p.pf()->interpolation(ctx_p, coeff, val, 1);

        for (size_type k = 0; k < qdim; ++k) output[k] -= val[0] * n[k];
        break;
      case 2:
        {
          base_tensor t;
          size_type ndof = ctx_p.pf()->nb_dof(cv);
          ctx_p.pf()->real_base_value(ctx_p, t);

          for (size_type i = 0; i < ndof; ++i)
            for (size_type k = 0; k < qdim; ++k)
              output(i, k) -= t[i]*n[k];
        }
        break;
      }

    }

    lin_incomp_Neumann_elem_term
    (const gmm::uint64_type &var_vnum_, const mesh_fem *mf_p_,
     const model_real_plain_vector *P_,
     const std::string &auxvarname)
      : var_vnum(var_vnum_), mf_p(mf_p_), org_P(P_)  {
      auxilliary_variables.push_back(auxvarname);
      gmm::resize(P, mf_p->nb_basic_dof());
      mf_p->extend_vector(*P_, P);
      vnum = var_vnum;
      gmm::resize(val, 1);
    }

  };



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


    virtual void real_post_assembly_in_serial(const model &md, size_type ib,
                                              const model::varnamelist &vl,
                                              const model::varnamelist &/*dl*/,
                                              const model::mimlist &/*mims*/,
                                              model::real_matlist &/*matl*/,
                                              model::real_veclist &,
                                              model::real_veclist &,
                                              size_type /*region*/,
                                              build_version) const
    {
        const mesh_fem &mf_p = md.mesh_fem_of_variable(vl[1]);
        pNeumann_elem_term pNt = new lin_incomp_Neumann_elem_term
          (md.version_number_of_data_variable( vl[1]), &mf_p,
                            &(md.real_variable(vl[1])), vl[1]);
        md.add_Neumann_term(pNt, vl[0], ib);
        md.add_auxilliary_variables_of_Neumann_terms(vl[0], vl[1]);
    }


    linear_incompressibility_brick(void) {
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
    pbrick pbr = new linear_incompressibility_brick();
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
    size_type ib = add_linear_generic_assembly_brick
      (md, mim, expr, region, true, true, "Linear incompressibility", true);
    if (ib == size_type(-1))
      ib = add_nonlinear_generic_assembly_brick
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

    mass_brick(void) {
      set_flags("Mass brick", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
    }

  };

  size_type add_mass_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataexpr_rho,  size_type region) {
    if (md.is_complex()) {
      pbrick pbr = new mass_brick;
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
      size_type ib = add_linear_generic_assembly_brick
        (md, mim, expr, region, true, true, "Mass matrix", true);
      if (ib == size_type(-1))
        ib = add_nonlinear_generic_assembly_brick
          (md, mim, expr, region, false, false, "Mass matrix (nonlinear)");
      return ib;
    }
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

    basic_d_on_dt_brick(void) {
      set_flags("Basic d/dt brick", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
    }

  };

  size_type add_basic_d_on_dt_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname_dt, const std::string &dataname_rho,
   size_type region) {
    pbrick pbr = new basic_d_on_dt_brick;
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

    basic_d2_on_dt2_brick(void) {
      set_flags("Basic d2/dt2 brick", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, true /* is complex */,
                false /* compute each time */, false /* has a Neumann term */);
    }

  };

  size_type add_basic_d2_on_dt2_brick
  (model &md, const mesh_im &mim, const std::string &varnameU,
   const std::string &datanameV,
   const std::string &dataname_dt,
   const std::string &dataname_alpha,
   const std::string &dataname_rho,
   size_type region) {
    pbrick pbr = new basic_d2_on_dt2_brick;
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
    pdispatcher pdispatch = new theta_method_dispatcher(THETA);
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
      md.assembly(model::BUILD_COMPLETE_RHS);

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
      md.assembly(model::BUILD_COMPLETE_RHS);

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
        clear(vectl[1]); clear(vectl_sym[1]);
      }

      if (pbr->is_linear()) { // If the problem is linear, add the term
        // coming from the previous iteration as a second rhs.
        // This rhs is only used for this.
        if (first_iter) md.update_brick(ib, model::BUILD_RHS);
        clear(vectl[1]); clear(vectl_sym[1]);
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
        clear(vectl[1]); clear(vectl_sym[1]);
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
        clear(vectl[1]); clear(vectl_sym[1]);
        md.linear_brick_add_to_rhs(ib, 1, 1);
      }

      md.reset_default_iter_of_variables(dl);
      if (!(pbr->is_linear()))
        md.reset_default_iter_of_variables(vl);
    }

    midpoint_dispatcher(void) : virtual_dispatcher(2)
    { id_num = act_counter(); }

  };

  void add_midpoint_dispatcher(model &md, dal::bit_vector ibricks) {
    pdispatcher pdispatch = new midpoint_dispatcher();
    for (dal::bv_visitor i(ibricks); !i.finished(); ++i)
      md.add_time_dispatcher(i, pdispatch);
  }


}  /* end of namespace getfem.                                             */

