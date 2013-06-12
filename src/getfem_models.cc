/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2009-2012 Yves Renard
 
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
#if GETFEM_HAVE_MUPARSER_MUPARSER_H
#include <muParser/muParser.h>
#elif GETFEM_HAVE_MUPARSER_H
#include <muParser.h>
#endif



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

  bool model::check_name_valitity(const std::string &name, bool assert) const {
    VAR_SET::const_iterator it = variables.find(name);
    if (it != variables.end()) {
      GMM_ASSERT1(!assert, "Variable " << name << " already exists");
      return false;
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
    bool valid = check_name_valitity(res_name, false);
    VAR_SET::const_iterator it = variables.find(res_name);
    GMM_ASSERT1(valid || it != variables.end(),
                "Illegal variable name : " << name);
    for (size_type i = 2; it != variables.end(); ++i) {
      std::stringstream m;
      m << name << '_' << i;
      res_name = m.str();
      it = variables.find(res_name);
    }
    return res_name;
  }

  void model::actualize_sizes(void) const {
    act_size_to_be_done = false;
    std::map<std::string, std::vector<std::string> > multipliers;
    std::map<std::string, bool > tobedone;

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
         ++it) {
      if (it->second.is_fem_dofs
	  && (it->second.filter == VDESCRFILTER_CTERM
	      || it->second.filter == VDESCRFILTER_INFSUP)) {
        VAR_SET::iterator it2 = variables.find(it->second.filter_var);
        GMM_ASSERT1(it2 != variables.end(), "The primal variable of the "
                    "multiplier does not exist");
        GMM_ASSERT1(it2->second.is_fem_dofs, "The primal variable of the "
                    "multiplier is not a fem variable");
        multipliers[it->second.filter_var].push_back(it->first);
        if (it->second.v_num < it->second.mf->version_number() ||
            it->second.v_num < it2->second.mf->version_number())
          tobedone[it->second.filter_var] = true;
      }
    }

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
          ++it) {
      if (it->second.is_fem_dofs) {
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
    }

    for (std::map<std::string, bool >::iterator itbd = tobedone.begin();
         itbd != tobedone.end(); ++itbd) {
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
      for (size_type k = 0; k < mults.size(); ++k) {
        VAR_SET::iterator it = variables.find(mults[k]);

        // This step forces the recomputation of corresponding bricks.
        // A test to check if a modification is really necessary could
        // be done first ... (difficult to coordinate with other multipliers)
        dal::bit_vector alldof; alldof.add(0, it->second.mf->nb_dof());
        it->second.partial_mf->adapt(alldof);
        it->second.set_size(it->second.partial_mf->nb_dof());

        // Obtaining the coupling matrix between the multipier and
        // the primal variable. A search is done on all the terms of the
        // model. Only the the corresponding linear terms are added.
        // If no linear term is available, a mass matrix is used.

        gmm::col_matrix< gmm::rsvector<scalar_type> >
          MM(it2->second.mf->nb_dof(), it->second.mf->nb_dof());
        bool termadded = false;

	if (it->second.filter == VDESCRFILTER_CTERM) {

	  for (size_type ib = 0; ib < bricks.size(); ++ib) {
	    const brick_description &brick = bricks[ib];
	    bool bupd = false;
	    bool cplx = is_complex() && brick.pbr->is_complex();
	    
	    for (size_type j = 0; j < brick.tlist.size(); ++j) {
	      
	      const term_description &term = brick.tlist[j];
	      
	      if (term.is_matrix_term && !mults[k].compare(term.var1) &&
		  !it2->first.compare(term.var2)) {
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
		
	      } else if (term.is_matrix_term && !mults[k].compare(term.var2) &&
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

	  if (!termadded)
	    GMM_WARNING1("No term found to filter multiplier " << it->first
			 << ". Variable is cancelled");
#if GETFEM_PARA_LEVEL > 1
	  if (termadded) {
            // we assume that all bricks take mpi_region into account but it
            // would be better if the brick itself could report if it supports
            // distributed assembly
            // This is only a reference implementation, it needs to be optimized
            // maybe by using gmm::mpi_distributed_matrix
            std::vector<scalar_type> tmpvec1(gmm::mat_nrows(MM)), tmpvec2(gmm::mat_nrows(MM));
            for (size_type k = 0; k < gmm::mat_ncols(MM); ++k) {
                gmm::copy(gmm::mat_col(MM,k),tmpvec1);
                MPI_SUM_VECTOR(tmpvec1,tmpvec2);
                gmm::copy(tmpvec2,gmm::mat_col(MM,k));
            }
          }
#endif
	} else if (it->second.filter == VDESCRFILTER_INFSUP) {
	  asm_mass_matrix(MM, *(it->second.mim), *(it2->second.mf),
			  *(it->second.mf), it->second.m_region);
	}
        //
        // filtering
        //
        std::set<size_type> columns;
        gmm::range_basis(MM, columns);
        if (mults.size() > 1) {
          gmm::copy(MM, gmm::sub_matrix
                    (MGLOB,gmm::sub_interval(0, it2->second.mf->nb_dof()),
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
          it->second.partial_mf->adapt(kept);
          it->second.set_size(it->second.partial_mf->nb_dof());
          it->second.v_num = act_counter();
        }
      }

      if (mults.size() > 1) {
        gmm::range_basis(MGLOB, glob_columns, 1E-12, gmm::col_major(), true);

        s = 0;
        for (size_type k = 0; k < mults.size(); ++k) {
          VAR_SET::iterator it = variables.find(mults[k]);
          dal::bit_vector kept;
          size_type nbdof = it->second.mf->nb_dof();
          for (std::set<size_type>::iterator itt = glob_columns.begin();
               itt != glob_columns.end(); ++itt)
            if (*itt >= s && *itt < s + nbdof) kept.add(*itt-s);
          it->second.partial_mf->adapt(kept);
          it->second.set_size(it->second.partial_mf->nb_dof());
          it->second.v_num = act_counter();
          s += it->second.mf->nb_dof();
        }
      }
    }

    size_type tot_size = 0;

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
         ++it)
      if (it->second.is_variable) {
        it->second.I = gmm::sub_interval(tot_size, it->second.size());
        tot_size += it->second.size();
      }

    if (complex_version) {
      gmm::resize(cTM, tot_size, tot_size);
      gmm::resize(crhs, tot_size);
    }
    else {
      gmm::resize(rTM, tot_size, tot_size);
      gmm::resize(rrhs, tot_size);
    }
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
	if (it->second.is_disabled) ost << "\t (disabled)";
	ost << endl;
      }
    }
  }

  void model::add_fixed_size_variable(const std::string &name, size_type size,
                                      size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(true, is_complex(), false, niter);
    act_size_to_be_done = true;
    variables[name].set_size(size);
  }

  void model::resize_fixed_size_variable(const std::string &name,
                                         size_type size) {
    GMM_ASSERT1(!(variables[name].is_fem_dofs), "Cannot explicitely resize "
                " a fem variable or data");
    variables[name].set_size(size);
  }

  void resize_fixed_size_variable(const std::string &name, size_type size);


  void model::add_fixed_size_data(const std::string &name, size_type size,
                                      size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(false, is_complex(), false, niter);
    variables[name].set_size(size);
  }

  void model::add_fem_variable(const std::string &name, const mesh_fem &mf,
                               size_type niter) {
    check_name_valitity(name);
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
    check_name_valitity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_REGION, &mf, region);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_fem_data(const std::string &name, const mesh_fem &mf,
                               dim_type qdim, size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(false, is_complex(), true, niter,
                                      VDESCRFILTER_NO, &mf, 0, qdim);
    variables[name].set_size(mf.nb_dof()*qdim);
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             const std::string &primal_name,
                             size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_CTERM, &mf, 0,
                                      1, primal_name);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  void model::add_multiplier(const std::string &name, const mesh_fem &mf,
                             const std::string &primal_name,const mesh_im &mim,
			     size_type region, size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
                                      VDESCRFILTER_INFSUP, &mf, region,
                                      1, primal_name, &mim);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    add_dependency(mf);
  }

  size_type model::add_brick(pbrick pbr, const varnamelist &varnames,
                             const varnamelist &datanames,
                             const termlist &terms,
                             const mimlist &mims, size_type region) {
    bricks.push_back(brick_description(pbr, varnames, datanames, terms,
                                       mims, region));
    size_type ib = bricks.size() - 1;
    active_bricks.add(ib);
    for  (size_type i = 0; i < bricks.back().mims.size(); ++i)
      add_dependency(*(bricks.back().mims[i]));

    GMM_ASSERT1(pbr->is_real() || is_complex(),
                "Impossible to add a complex brick to a real model");
    if (is_complex() && pbr->is_complex()) {
      bricks.back().cmatlist.resize(terms.size());
      bricks.back().cveclist[0].resize(terms.size());
      bricks.back().cveclist_sym[0].resize(terms.size());
    } else {
      bricks.back().rmatlist.resize(terms.size());
      bricks.back().rveclist[0].resize(terms.size());
      bricks.back().rveclist_sym[0].resize(terms.size());
    }
    is_linear_ = is_linear_ && pbr->is_linear();
    is_symmetric_ = is_symmetric_ && pbr->is_symmetric();
    is_coercive_ = is_coercive_ && pbr->is_coercive();

    for (size_type i=0; i < varnames.size(); ++i)
      GMM_ASSERT1(variables.find(varnames[i]) != variables.end(),
                  "Undefined model variable " << varnames[i]);
    for (size_type i=0; i < datanames.size(); ++i)
      GMM_ASSERT1(variables.find(datanames[i]) != variables.end(),
                  "Undefined model data or variable " << datanames[i]);

    return ib;
  }

  void model::add_mim_to_brick(size_type ib, const mesh_im &mim) {
    GMM_ASSERT1(ib < bricks.size(), "Inexistent brick");
    touch_brick(ib);
    bricks[ib].mims.push_back(&mim);
    add_dependency(mim);
  }

  void model::change_terms_of_brick(size_type ib, const termlist &terms) {
    GMM_ASSERT1(ib < bricks.size(), "Inexistent brick");
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
    GMM_ASSERT1(ib < bricks.size(), "Inexistent brick");
    touch_brick(ib);
    bricks[ib].vlist = vl;
    for (size_type i=0; i < vl.size(); ++i)
      GMM_ASSERT1(variables.find(vl[i]) != variables.end(),
                  "Undefined model variable " << vl[i]);
  }


  void model::add_time_dispatcher(size_type ibrick, pdispatcher pdispatch) {

    GMM_ASSERT1(ibrick < bricks.size(), "Inexistent brick");

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
    GMM_ASSERT1(ind_brick < bricks.size(), "Inexistent brick");
    GMM_ASSERT1(ind_var < bricks[ind_brick].vlist.size(),
               "Inexistent brick variable");
    return bricks[ind_brick].vlist[ind_var];
  }

  const std::string &model::dataname_of_brick(size_type ind_brick,
                                              size_type ind_data) {
    GMM_ASSERT1(ind_brick < bricks.size(), "Inexistent brick");
    GMM_ASSERT1(ind_data < bricks[ind_brick].dlist.size(),
                "Inexistent brick data");
    return bricks[ind_brick].dlist[ind_data];
  }

  void model::listbricks(std::ostream &ost, size_type base_id) const {
    if (bricks.size() == 0)
      ost << "Model with no bricks" << endl;
    else {
      ost << "List of model bricks:" << endl;
      for (size_type i = 0; i < bricks.size(); ++i) {
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


  void model::brick_call(size_type ib, build_version version,
                         size_type rhs_ind) const {
    const brick_description &brick = bricks[ib];
    bool cplx = is_complex() && brick.pbr->is_complex();

    brick_init(ib, version, rhs_ind);

    if (cplx)
      brick.pbr->asm_complex_tangent_terms(*this, ib, brick.vlist, brick.dlist,
                                           brick.mims,
                                           brick.cmatlist,
                                           brick.cveclist[rhs_ind],
                                           brick.cveclist_sym[rhs_ind],
                                           brick.region, version);
    else
      brick.pbr->asm_real_tangent_terms(*this, ib, brick.vlist, brick.dlist,
                                        brick.mims,
                                        brick.rmatlist,
                                        brick.rveclist[rhs_ind],
                                        brick.rveclist_sym[rhs_ind],
                                        brick.region, version);
  }

  void model::set_dispatch_coeff(void) {
    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];
      if (brick.pdispatch)
        brick.pdispatch->set_dispatch_coeff(*this, ib);

    }
  }

  void model::first_iter(void) {

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

  void model::auxilliary_variables_of_Neumann_terms
  (const std::string &varname, std::vector<std::string> &aux_vars) const {
    std::map<std::string, std::vector<std::string> >::const_iterator
      it = Neumann_terms_auxilliary_variables.find(varname);
    if (it !=  Neumann_terms_auxilliary_variables.end())
      aux_vars = it->second;
    else
      aux_vars.resize(0);
  }

  void model::add_auxilliary_variables_of_Neumann_terms
  (const std::string &varname, const std::vector<std::string> &aux_vars) {

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
  (const std::string &varname, const std::string &aux_var) {
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
    for (size_type i = 0; i < brick.vlist.size() && !tobecomputed; ++i) {
      var_description &vd = variables[brick.vlist[i]];
      if (vd.v_num > brick.v_num) tobecomputed = true;
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
                 gmm::scaled(V,  std::complex<scalar_type>(-1)),
                 brick.cveclist[ind_data][j]);
            } else
              gmm::mult_add
                (brick.cmatlist[j],
                 gmm::scaled(variables[term.var2].complex_value[n_iter_2],
                             std::complex<scalar_type>(-1)),
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
                             std::complex<scalar_type>(-1)),
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


  bool model::build_reduced_index(std::vector<size_type> &ind) {
    ind.resize(0);
    bool reduced = false;
    for (VAR_SET::iterator it = variables.begin(); it != variables.end(); ++it)
      if (it->second.is_variable) {
	if  (it->second.is_disabled)
	  reduced = true;
	else {
	  for (size_type i=it->second.I.first(); i < it->second.I.last(); ++i)
	    ind.push_back(i);
	}
      }
    return reduced;
  }


  void model::assembly(build_version version) {

    context_check(); if (act_size_to_be_done) actualize_sizes();
    if (is_complex()) {
      if (version & BUILD_MATRIX) gmm::clear(cTM);
      if (version & BUILD_RHS) gmm::clear(crhs);
      if (version & BUILD_PSEUDO_POTENTIAL) pseudo_potential_ = scalar_type(0);
    }
    else {
      if (version & BUILD_MATRIX) gmm::clear(rTM);
      if (version & BUILD_RHS) gmm::clear(rrhs);
      if (version & BUILD_PSEUDO_POTENTIAL) pseudo_potential_ = scalar_type(0);
    }
    clear_dof_constraints();

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {

      brick_description &brick = bricks[ib];
      
      // Disables the brick if all its variables are disabled.
      bool auto_disabled_brick = true;
      for (size_type j = 0; j < brick.vlist.size(); ++j) {
	if (!(variables[brick.vlist[j]].is_disabled))
	  auto_disabled_brick = false;
      }
      if (auto_disabled_brick) continue;

      update_brick(ib, version);

      bool cplx = is_complex() && brick.pbr->is_complex();

      scalar_type coeff0 = scalar_type(1);
      if (brick.pdispatch) coeff0 = brick.matrix_coeff;

      if (version & BUILD_PSEUDO_POTENTIAL) {

        scalar_type pseudop = scalar_type(0);
        if (cplx)
          pseudop = brick.pbr->asm_complex_pseudo_potential
            (*this, ib, brick.vlist, brick.dlist, brick.mims, brick.cmatlist,
             brick.cveclist[0], brick.cveclist_sym[0], brick.region);
        else
          pseudop = brick.pbr->asm_real_pseudo_potential
            (*this, ib, brick.vlist, brick.dlist, brick.mims, brick.rmatlist,
             brick.rveclist[0], brick.rveclist_sym[0], brick.region);

        pseudo_potential_ += pseudop * coeff0;

        GMM_ASSERT1(!(brick.pdispatch), "Pseudo potential not "
                    "supported by brick dispatcher, sorry");

      }

      // Assembly of terms

      for (size_type j = 0; j < brick.tlist.size(); ++j) {
        term_description &term = brick.tlist[j];
        bool isg = term.is_global;
        size_type nbgdof = nb_dof();
        gmm::sub_interval I1(0,nbgdof), I2(0,nbgdof);
        if (!isg) I1 = variables[term.var1].I;
        if (term.is_matrix_term && !isg) I2 = variables[term.var2].I;

        if (cplx) {
          if (term.is_matrix_term && (version & BUILD_MATRIX)) {
            gmm::add(gmm::scaled(brick.cmatlist[j], coeff0),
                     gmm::sub_matrix(cTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first()) {
              gmm::add(gmm::scaled(gmm::transposed(brick.cmatlist[j]), coeff0),
                       gmm::sub_matrix(cTM, I2, I1));
            }
          }
          if (version & BUILD_RHS) {
            if (brick.pdispatch) {
              for (size_type k = 0; k < brick.nbrhs; ++k)
                gmm::add(gmm::scaled(brick.cveclist[k][j],
                                     brick.coeffs[k]),
                         gmm::sub_vector(crhs, I1));
            }
            else
              gmm::add(brick.cveclist[0][j], gmm::sub_vector(crhs, I1));
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
              if (isg) {
                model_complex_plain_vector V(nbgdof);
                from_variables(V);
                gmm::mult_add(brick.cmatlist[j],
                            gmm::scaled(V, std::complex<scalar_type>(-coeff0)),
                            crhs);
              }
              else
                gmm::mult_add(brick.cmatlist[j],
                            gmm::scaled(variables[term.var2].complex_value[0],
                                        std::complex<scalar_type>(-coeff0)),
                            gmm::sub_vector(crhs, I1));
            }
            if (term.is_symmetric && I1.first() != I2.first()) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.cveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I2));
              }
              else
                gmm::add(brick.cveclist_sym[0][j], gmm::sub_vector(crhs, I2));
               if (brick.pbr->is_linear()
                   && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
                 gmm::mult_add(gmm::conjugated(brick.cmatlist[j]),
                            gmm::scaled(variables[term.var1].complex_value[0],
                                        std::complex<scalar_type>(-coeff0)),
                            gmm::sub_vector(crhs, I2));
               }
            }
          }
        } else if (is_complex()) {
          if (term.is_matrix_term && (version & BUILD_MATRIX)) {
            gmm::add(gmm::scaled(brick.rmatlist[j], coeff0),
                     gmm::sub_matrix(cTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first()) {
              gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), coeff0),
                       gmm::sub_matrix(cTM, I2, I1));
            }
          }
          if (version & BUILD_RHS) {
            if (brick.pdispatch) {
              for (size_type k = 0; k < brick.nbrhs; ++k)
                gmm::add(gmm::scaled(brick.rveclist[k][j],
                                     brick.coeffs[k]),
                         gmm::sub_vector(crhs, I1));
            }
            else
              gmm::add(brick.rveclist[0][j], gmm::sub_vector(crhs, I1));
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
              if (isg) {
                model_complex_plain_vector V(nbgdof);
                from_variables(V);
                gmm::mult_add(brick.rmatlist[j],
                            gmm::scaled(V, std::complex<scalar_type>(-coeff0)),
                            crhs);
              }
              else
                gmm::mult_add(brick.rmatlist[j],
                            gmm::scaled(variables[term.var2].complex_value[0],
                                        std::complex<scalar_type>(-coeff0)),
                            gmm::sub_vector(crhs, I1));
            }
            if (term.is_symmetric && I1.first() != I2.first()) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(crhs, I2));
              }
              else
                gmm::add(brick.rveclist_sym[0][j], gmm::sub_vector(crhs, I2));
              if (brick.pbr->is_linear()
                  && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                             gmm::scaled(variables[term.var1].complex_value[0],
                                          std::complex<scalar_type>(-coeff0)),
                              gmm::sub_vector(crhs, I2));
              }
            }
          }
        } else {
          if (term.is_matrix_term && (version & BUILD_MATRIX)) {
            gmm::add(gmm::scaled(brick.rmatlist[j], coeff0),
                     gmm::sub_matrix(rTM, I1, I2));
            if (term.is_symmetric && I1.first() != I2.first()) {
              gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), coeff0),
                       gmm::sub_matrix(rTM, I2, I1));
            }
          }
          if (version & BUILD_RHS) {
            if (brick.pdispatch) {
              for (size_type k = 0; k < brick.nbrhs; ++k)
                gmm::add(gmm::scaled(brick.rveclist[k][j],
                                     brick.coeffs[k]),
                         gmm::sub_vector(rrhs, I1));
            }
            else
              gmm::add(brick.rveclist[0][j], gmm::sub_vector(rrhs, I1));
            if (term.is_matrix_term && brick.pbr->is_linear()
                && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
              if (isg) {
                model_real_plain_vector V(nbgdof);
                from_variables(V);
                gmm::mult_add(brick.rmatlist[j],
                              gmm::scaled(V, -coeff0), rrhs);
              }
              gmm::mult_add(brick.rmatlist[j],
                            gmm::scaled(variables[term.var2].real_value[0],
                                        -coeff0),
                            gmm::sub_vector(rrhs, I1));
            }
            if (term.is_symmetric && I1.first() != I2.first()) {
              if (brick.pdispatch) {
                for (size_type k = 0; k < brick.nbrhs; ++k)
                  gmm::add(gmm::scaled(brick.rveclist_sym[k][j],
                                       brick.coeffs[k]),
                           gmm::sub_vector(rrhs, I2));
              }
              else
                gmm::add(brick.rveclist_sym[0][j], gmm::sub_vector(rrhs, I2));
              if (brick.pbr->is_linear()
                  && (!is_linear() || (version & BUILD_WITH_COMPLETE_RHS))) {
                gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
                              gmm::scaled(variables[term.var1].real_value[0],
                                          -coeff0),
                              gmm::sub_vector(rrhs, I2));
              }
            }
          }
        }
      }

      if (brick.pbr->is_linear())
        brick.terms_to_be_computed = false;
// Commented to allow to get this information. Should be optional ?
//       else
//         if (cplx) {
//           brick.cmatlist = complex_matlist(brick.tlist.size());
//           brick.cveclist[0] = complex_veclist(brick.tlist.size());
//         } else {
//           brick.rmatlist = real_matlist(brick.tlist.size());
//           brick.rveclist[0] = real_veclist(brick.tlist.size());
//         }

    }

    if (version & BUILD_RHS) {
      if (is_complex()) MPI_SUM_VECTOR(crhs); else MPI_SUM_VECTOR(rrhs);
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
        // In parallel, a unification of the indices and values could be done
        // in order to allow the bricks to compute dof constraints in a
        // distributed way. Not done for the moment.
        
        if (dof_indices.size()) {
          gmm::sub_index SI(dof_indices);
          gmm::sub_interval II(0, nb_dof());
          
          if (version & BUILD_RHS) {
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
    context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::const_iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.real_value[niter];
  }

  const model_complex_plain_vector &
  model::complex_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::const_iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter  > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.complex_value[niter];
  }

  model_real_plain_vector &
  model::set_real_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(!complex_version, "This model is a complex one");
    context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
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
    context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    it->second.v_num_data = act_counter();
    if (niter == size_type(-1)) niter = it->second.default_iter;
    GMM_ASSERT1(it->second.n_iter + it->second.n_temp_iter > niter,
                "Invalid iteration number "
                << niter << " for " << name);
    return it->second.complex_value[niter];
  }
    void model::check_brick_stiffness_rhs(size_type ind_brick) const
	{

	  
	  const brick_description &brick = bricks[ind_brick];
	  update_brick(ind_brick, model::BUILD_ALL);

      brick.pbr->check_stiffness_matrix_and_rhs(*this, ind_brick, brick.tlist,
		  brick.vlist, brick.dlist, brick.mims, brick.rmatlist,
            brick.rveclist[0], brick.rveclist_sym[0], brick.region);
  }


  // ----------------------------------------------------------------------
  //
  //
  // Standard bricks
  //
  //
  // ----------------------------------------------------------------------
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
        const scalar_type TINY) const {
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
            asm_real_tangent_terms(md, s, vl, dl, mims, matl, rvc1, rvc2,
                rg, model::BUILD_MATRIX);
            for(size_type iterm=0;iterm<matl.size();iterm++){

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
                asm_real_tangent_terms(md, s, vl, dl, mims, matl, rvc1, rvc2,
                    rg, model::BUILD_RHS);
                if (gmm::mat_euclidean_norm(matl[iterm])<1e-12){
                    std::cout<<"    The assembled matrix is nearly zero, skipping."<<std::endl;
                    continue;
                }
                model_real_plain_vector RHS0(rvc1[rhs_index[tlist[iterm].var1]]);

                //finite difference stiffness		
                model_real_sparse_matrix fdSM(matl[iterm].nrows(),matl[iterm].ncols());
                model_real_plain_vector&U = md.set_real_variable(tlist[iterm].var2);
                model_real_plain_vector& RHS1 =rvc1[rhs_index[tlist[iterm].var1]];
                for (size_type j=0; j < matl[iterm].ncols(); j++){
                    U[j]+=TINY;
                    gmm::fill(RHS1, 0.0);
                    asm_real_tangent_terms(md, s, vl, dl, mims, matl, rvc1, rvc2,
                        rg, model::BUILD_RHS);
                    for (size_type i=0;i<matl[iterm].nrows();i++)
                        fdSM(i,j) = (RHS0[i]-RHS1[i])/TINY;
                    U[j]-=TINY;
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
  // Generic elliptic brick
  //
  // ----------------------------------------------------------------------


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

	coeff.resize(mf_a->nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector(var, gmm::sub_index
			      (mfvar.ind_basic_dof_of_element(cv))), coeff);
	ctx_a.pf()->interpolation(ctx_a, coeff, val, dim_type(s));
      } else if (A) {
	gmm::copy(*A, val);
      } else {
	val[0] = scalar_type(1);
      }

      switch (version) {
      case 1:
	gmm::resize(grad, Q, N);
	coeff.resize(mfvar.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector(var, gmm::sub_index
			      (mfvar.ind_basic_dof_of_element(cv))), coeff);
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
// 	    for (size_type l = 0; l < qmult; ++l) {
// 	      for (size_type p = 0; p < Q; ++p) {
// 		base_tensor::const_iterator it = t.begin();
// 		for (size_type k = 0; k < Q; ++k)
// 		  for (size_type j = 0; j < tdim; ++j)
// 		    for (size_type i = 0; i < ndof; ++i, ++it) {
// 		      size_type jj = j + tdim*l;
// 		      if (p == jj) output(i*qmult+l, p) += val[0]*(*it)*n[k];
// 		    }
// 		GMM_ASSERT1(it ==  t.end(), "Internal error");
// 	      }
// 	    }
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




  struct generic_elliptic_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
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

      pNeumann_elem_term pNt = new generic_elliptic_Neumann_elem_term(mf_a, A);
      md.add_Neumann_term(pNt, vl[0], ib);
    }

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
                                                  const model::varnamelist &vl,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::real_matlist &matl,
                                                  model::real_veclist &,
                                                  model::real_veclist &,
                                                  size_type) const {
      const model_real_plain_vector &U = md.real_variable(vl[0]);
      return gmm::vect_sp(matl[0], U, U) / scalar_type(2);
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
    pbrick pbr = new generic_elliptic_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    return md.add_brick(pbr, model::varnamelist(1, varname),
                        model::varnamelist(), tl, model::mimlist(1, &mim),
                        region);
  }

  size_type add_generic_elliptic_brick(model &md, const mesh_im &mim,
                                       const std::string &varname,
                                       const std::string &dataname,
                                       size_type region) {
    pbrick pbr = new generic_elliptic_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    return md.add_brick(pbr, model::varnamelist(1, varname),
                        model::varnamelist(1, dataname), tl,
                        model::mimlist(1, &mim), region);
  }

  // ----------------------------------------------------------------------
  //
  // Source term brick
  //
  // ----------------------------------------------------------------------

  struct source_term_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
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

    }

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
                                                  const model::varnamelist &vl,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::real_matlist &,
                                                  model::real_veclist &vecl,
                                                  model::real_veclist &,
                                                  size_type) const {
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      return -gmm::vect_sp(vecl[0], u);
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type,
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

    }

    virtual scalar_type asm_complex_pseudo_potential(const model &md,size_type,
                                                 const model::varnamelist &vl,
                                                 const model::varnamelist &,
                                                 const model::mimlist &,
                                                 model::complex_matlist &,
                                                 model::complex_veclist &vecl,
                                                 model::complex_veclist &,
                                                 size_type) const {
      const model_complex_plain_vector &u = md.complex_variable(vl[0]);
      return -gmm::real(gmm::vect_hp(vecl[0], u)); /* ? */
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
                                  const std::string &dataname,
                                  size_type region,
                                  const std::string &directdataname) {
    pbrick pbr = new source_term_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname));
    model::varnamelist vdata(1, dataname);
    if (directdataname.size()) vdata.push_back(directdataname);
    return md.add_brick(pbr, model::varnamelist(1, varname),
                        vdata, tl, model::mimlist(1, &mim), region);
  }

  // ----------------------------------------------------------------------
  //
  // Normal source term brick
  //
  // ----------------------------------------------------------------------

  struct normal_source_term_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md, size_type,
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

    }

    virtual void asm_complex_tangent_terms(const model &md, size_type,
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

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
                                                  const model::varnamelist &vl,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::real_matlist &,
                                                  model::real_veclist &vecl,
                                                  model::real_veclist &,
                                                  size_type) const {
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      return -gmm::vect_sp(vecl[0], u);
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
                                         const std::string &dataname,
                                         size_type region) {
    pbrick pbr = new normal_source_term_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname));
    model::varnamelist vdata(1, dataname);
    return md.add_brick(pbr, model::varnamelist(1, varname),
                        vdata, tl, model::mimlist(1, &mim), region);
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
    mutable model_real_sparse_matrix rB;
    mutable model_real_plain_vector rV;
    mutable model_complex_sparse_matrix cB;
    mutable model_complex_plain_vector cV;

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
	size_type ss = s * ((normal_component) ? mf_u.linked_mesh().dim() : 1);
        GMM_ASSERT1(mf_u.get_qdim() == ss, dl[ind] << ": bad format of "
		    "Dirichlet data. Detected dimension is " << ss
		    << " should be " << size_type(mf_u.get_qdim()));
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
               "M(#1,#3)+=comp(Base(#1).Base(#3).Base(#2)(:,:,i).F(i))"
               : "F=data(qdim(#1),qdim(#1),#2);"
               "M(#1,#3)+=comp(vBase(#1).vBase(#3).Base(#2)(:,i,:,j,k).F(i,j,k));", &mf_u);
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
          GMM_ASSERT1(mf_data == &mf_u, "Sorry, for this brick, the data has "
                     "to be define on the same f.e.m. than the unknown");
        } else {
          s = gmm::vect_size(*A);
          GMM_ASSERT1(mf_u.get_qdim() == s, ": bad format of "
		    "Dirichlet data. Detected dimension is " << s
		    << " should be " << size_type(mf_u.get_qdim())); 
        }
      }

      mesh_region rg(region);
      // mf_u.linked_mesh().intersect_with_mpi_region(rg);

      if (mf_u.get_qdim() > 1) {
        for (mr_visitor i(rg, mf_u.linked_mesh()); !i.finished(); ++i) {
          if (mf_u.fem_of_element(i.cv()))
            GMM_ASSERT1((mf_u.fem_of_element(i.cv()))->target_dim() == 1,
                        "Intrinsically vectorial fems are not allowed");
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

    virtual void asm_complex_tangent_terms(const model &md, size_type /*ib*/,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type region,
                                           build_version /*version*/) const {
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
          GMM_ASSERT1(mf_data == &mf_u, "Sorry, for this brick, the data has "
                     "to be define on the same f.e.m. than the unknown");
        } else {
          s = gmm::vect_size(*A);
          GMM_ASSERT1(mf_u.get_qdim() == s, ": bad format of "
		    "Dirichlet data. Detected dimension is " << s
		    << " should be " << size_type(mf_u.get_qdim())); 
        }
      }

      mesh_region rg(region);
      // mf_u.linked_mesh().intersect_with_mpi_region(rg); // Not distributed
      // for the moment. To distribute, model::assembly should gather the 
      // dof constraints.

      if (mf_u.get_qdim() > 1) {
        for (mr_visitor i(rg); !i.finished(); ++i)
          GMM_ASSERT1((mf_u.fem_of_element(i.cv()))->target_dim() ==1,
                      "Intrinsically vectorial fems are not allowed");
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
	  coeff.resize(mf_u->nb_basic_dof_of_element(cv));
	  gmm::copy(gmm::sub_vector(U, gmm::sub_index
			       (mf_u->ind_basic_dof_of_element(cv))), coeff);
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
	  //	      "a scalar finite element method");
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

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type,
                                        build_version version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Pointwize constraints brick only one term");
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

    virtual void asm_complex_tangent_terms(const model &md, size_type ib,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type,
                                           build_version version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Pointwize constraints brick only one term");
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
	gmm::row_matrix<model_complex_sparse_vector> BB(nb_co*Q,mf_u.nb_dof());
	gmm::clear(cB); gmm::resize(cB, nb_co, mf_u.nb_dof());
	dal::bit_vector dof_untouched;
	getfem::mesh_trans_inv mti(mf_u.linked_mesh());
	base_node pt(N);
	for (size_type i = 0; i < nb_co; ++i) {
	  gmm::copy(gmm::real_part(gmm::sub_vector(PT,
				   gmm::sub_interval(i*N, N))), pt);
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
                                const std::string &dataname,
                                size_type region) {
    pbrick pbr = new Helmholtz_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    return md.add_brick(pbr, model::varnamelist(1, varname),
                        model::varnamelist(1, dataname), tl,
                        model::mimlist(1, &mim), region);
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
                                    const std::string &dataname,
                                    size_type region) {
    pbrick pbr = new Fourier_Robin_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    return md.add_brick(pbr, model::varnamelist(1, varname),
                        model::varnamelist(1, dataname), tl,
                        model::mimlist(1, &mim), region);
  }

  // ----------------------------------------------------------------------
  //
  // Basic nonlinear brick
  //
  // ----------------------------------------------------------------------
      
#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H

  class basic_nonlinear_term : public nonlinear_elem_term {

  public:
    dim_type N;
    const mesh_fem &mf_u;
    std::vector<scalar_type> U;
    scalar_type param;
    base_small_vector V;
    base_vector coeff;
    std::string f, dfdu, varname, paramname;
    mu::Parser parser;
    bgeot::multi_index sizes_;
    int version;
    
    template <class VECT> basic_nonlinear_term
    (const mesh_fem &mf_u_, const VECT &U_, scalar_type param_,
     const std::string &f_, const std::string &dfdu_,
     const std::string &varname_, const std::string &paramname_,
     int version_)
      : N(mf_u_.linked_mesh().dim()), mf_u(mf_u_), U(mf_u.nb_basic_dof()),
	param(param_), f(f_), dfdu(dfdu_), varname(varname_),
	paramname(paramname_), version(version_) {
      sizes_.resize(1); sizes_[0] = 1;
      V.resize(1);
      mf_u.extend_vector(U_, U);

      switch (version) {
      case 0 : parser.SetExpr(dfdu); break; // tangent matrix
      case 1 : parser.SetExpr(f); break; // rhs
      }
      parser.DefineVar(varname, &V[0]);
      if (paramname.size()) parser.DefineVar(paramname, &param);
    }

    const bgeot::multi_index &sizes(size_type) const { return sizes_; }

    virtual void compute(fem_interpolation_context &ctx,
			 bgeot::base_tensor &t) {
      size_type cv = ctx.convex_num();
      t.adjust_sizes(sizes_);
      coeff.resize(mf_u.nb_basic_dof_of_element(cv));
      gmm::copy(gmm::sub_vector
                (U, gmm::sub_index
                 (mf_u.ind_basic_dof_of_element(cv))), coeff);
      ctx.pf()->interpolation(ctx, coeff, V, 1);

      try {
	t[0] = scalar_type(parser.Eval());
      } catch (mu::Parser::exception_type &e) {
	std::cerr << "Message  : " << e.GetMsg()   << std::endl;
	std::cerr << "Formula  : " << e.GetExpr()  << std::endl;
	std::cerr << "Token    : " << e.GetToken() << std::endl;
	std::cerr << "Position : " << e.GetPos()   << std::endl;
	std::cerr << "Errc     : " << e.GetCode()  << std::endl;
	GMM_ASSERT1(false, "error in the expressions");
      }
    }

  };

  template<typename MAT, typename VECT>
  void asm_basic_nonlinear_tangent_matrix
  (const MAT &K, const mesh_im &mim, const mesh_fem &mf_u, const VECT &U,
   const std::string &f, const std::string &dfdu, const std::string &varname,
   const mesh_region &rg = mesh_region::all_convexes(),
   scalar_type param = 0., const std::string &paramname = std::string()) {

    basic_nonlinear_term nterm(mf_u, U, param, f, dfdu, 
			       varname, paramname, 0);

    generic_assembly assem;
    assem.set("M(#1,#1)+=sym(comp(NonLin(#1).Base(#1).Base(#1))(i,:,:))");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_nonlinear_term(&nterm);
    assem.push_mat(const_cast<MAT &>(K));
    assem.assembly(rg);
  }

  template<typename VECT>
  void asm_basic_nonlinear_rhs
  (const VECT &V, const mesh_im &mim, const mesh_fem &mf_u, const VECT &U,
   const std::string &f, const std::string &dfdu, const std::string &varname,
   const mesh_region &rg = mesh_region::all_convexes(),
   scalar_type param = 0., const std::string &paramname = std::string()) {

    basic_nonlinear_term nterm(mf_u, U, param, f, dfdu,
			       varname, paramname, 1);

    generic_assembly assem;
    assem.set("V(#1)+=comp(NonLin(#1).Base(#1))(i,:)");
    assem.push_mi(mim);
    assem.push_mf(mf_u);
    assem.push_nonlinear_term(&nterm);
    assem.push_vec(const_cast<VECT &>(V));
    assem.assembly(rg);
  }

#endif

  struct basic_nonlinear_brick : public virtual_brick {

    std::string f, dfdu;
    
    virtual void asm_real_tangent_terms(const model &md, size_type /* ib */,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type region,
                                        build_version version) const {
      GMM_ASSERT1(mims.size() == 1,
		  "basic nonlinear brick needs a single mesh_im");
      GMM_ASSERT1(vl.size() == 1,
		  "basic nonlinear brick needs a single variable");
      GMM_ASSERT1(dl.size() <= 1,
		  "wrong number of data for basic nonlinear brick");
      GMM_ASSERT1(matl.size() == 1,  "wrong number of terms for basic "
		  "nonlinear brick");

      const model_real_plain_vector &u = md.real_variable(vl[0]);
      const mesh_fem &mf_u = *(md.pmesh_fem_of_variable(vl[0]));
      size_type Q = mf_u.get_qdim();
      GMM_ASSERT1(Q == 1, "basic nonlinear brick is only for scalar field, "
		  "sorry");

      const model_real_plain_vector *vparam = 0;
      if (dl.size()) {
	vparam = &md.real_variable(dl[0]) ;
	size_type sl = gmm::vect_size(*vparam);
	GMM_ASSERT1(sl == 1, "the parameter in basic nonlinear brick "
		    "has to be scalar");
      }

      const mesh_im &mim = *mims[0];
      mesh_region rg(region);
      mf_u.linked_mesh().intersect_with_mpi_region(rg);

#if GETFEM_HAVE_MUPARSER_MUPARSER_H || GETFEM_HAVE_MUPARSER_H

      if (version & model::BUILD_MATRIX) {
	gmm::clear(matl[0]);
	GMM_TRACE2("basic nonlinear stiffness matrix assembly");
	if (dl.size())
	  asm_basic_nonlinear_tangent_matrix(matl[0], mim, mf_u, u, f, dfdu,
					     vl[0], rg, (*vparam)[0], dl[0]);
	else
	  asm_basic_nonlinear_tangent_matrix(matl[0], mim, mf_u, u,
					     f, dfdu, vl[0], rg);
      }

      if (version & model::BUILD_RHS) {
	if (dl.size())
	  asm_basic_nonlinear_rhs(vecl[0], mim, mf_u, u,
				  f, dfdu, vl[0], rg, (*vparam)[0], dl[0]);
	else
	  asm_basic_nonlinear_rhs(vecl[0], mim, mf_u, u, f, dfdu, vl[0], rg);
	gmm::scale(vecl[0], scalar_type(-1));
      }
      
#else

    GMM_ASSERT1(false, "Muparser is not installed, "
                    "you cannot use basic nonlinear brick");

#endif
    }

    basic_nonlinear_brick(const std::string &f_, const std::string &dfdu_)
      : f(f_), dfdu(dfdu_)
    { set_flags("basic nonlinear brick", false /* is linear*/,
		true /* is symmetric */, false /* is coercive */,
		true /* is real */, false /* is complex */,
		false /* compute each time */, false /* has a Neumann term */);
    }
    
  };

  size_type add_basic_nonlinear_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &f, const std::string &dfdu, size_type region,
   const std::string &dataname) {
    pbrick pbr = new basic_nonlinear_brick(f, dfdu);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl;
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
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

    virtual void asm_real_tangent_terms(const model &md, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type, build_version) const {
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

        gmm::mult(gmm::transposed(rB), gmm::scaled(rL, gmm::abs((*COEFF)[0])),
                  vecl[0]);
        gmm::mult(gmm::transposed(rB), gmm::scaled(rB, gmm::abs((*COEFF)[0])),
                  matl[0]);
      } else {
        gmm::copy(rL, vecl[0]);
        gmm::copy(rB, matl[0]);
      }
    }

    virtual void asm_complex_tangent_terms(const model &md, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type,
                                           build_version) const {
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

        gmm::mult(gmm::transposed(cB), gmm::scaled(cL, gmm::abs((*COEFF)[0])),
                  vecl[0]);
        gmm::mult(gmm::transposed(cB), gmm::scaled(cB, gmm::abs((*COEFF)[0])),
                  matl[0]);
      } else {
        gmm::copy(cL, vecl[0]);
        gmm::copy(cB, matl[0]);
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

    virtual void asm_real_tangent_terms(const model &, size_type,
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

    virtual void asm_complex_tangent_terms(const model &, size_type,
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

    virtual void asm_real_tangent_terms(const model &, size_type,
                                        const model::varnamelist &vl,
                                        const model::varnamelist &dl,
                                        const model::mimlist &mims,
                                        model::real_matlist &matl,
                                        model::real_veclist &vecl,
                                        model::real_veclist &,
                                        size_type, build_version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Explicit rhs has one and only one term");
      GMM_ASSERT1(mims.size() == 0, "Explicit rhs need no mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 0,
                  "Wrong number of variables for explicit rhs brick");
      gmm::copy(rL, vecl[0]);
    }

    virtual void asm_complex_tangent_terms(const model &, size_type,
                                           const model::varnamelist &vl,
                                           const model::varnamelist &dl,
                                           const model::mimlist &mims,
                                           model::complex_matlist &matl,
                                           model::complex_veclist &vecl,
                                           model::complex_veclist &,
                                           size_type,
                                           build_version) const {
      GMM_ASSERT1(vecl.size() == 1 && matl.size() == 1,
                  "Explicit rhs has one and only one term");
      GMM_ASSERT1(mims.size() == 0, "Explicit rhs need no mesh_im");
      GMM_ASSERT1(vl.size() == 1 && dl.size() == 0,
                  "Wrong number of variables for explicit rhs brick");
      gmm::copy(cL, vecl[0]);

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

	coeff.resize(mf_mu->nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector(*mu, gmm::sub_index
			      (mf_mu->ind_basic_dof_of_element(cv))), coeff);
	ctx_mu.pf()->interpolation(ctx_mu, coeff, val, 1);
	val_mu = val[0];
	gmm::copy(gmm::sub_vector(*lambda, gmm::sub_index
			      (mf_mu->ind_basic_dof_of_element(cv))), coeff);
	ctx_mu.pf()->interpolation(ctx_mu, coeff, val, 1);
	val_mu = val[0];
      } else {
	val_lambda = (*lambda)[0]; val_mu = (*mu)[0];
      }

      switch (version) {
      case 1:
	coeff.resize(mfvar.nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector(var, gmm::sub_index
			      (mfvar.ind_basic_dof_of_element(cv))), coeff);
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


	pNeumann_elem_term pNt = new iso_lin_elasticity_Neumann_elem_term
	  (mf_lambda, lambda, mf_mu, mu);
	md.add_Neumann_term(pNt, vl[0], ib);
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

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
                                                  const model::varnamelist &vl,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::real_matlist &matl,
                                                  model::real_veclist &,
                                                  model::real_veclist &,
                                                  size_type) const {
      const model_real_plain_vector &u = md.real_variable(vl[0]);
      return gmm::vect_sp(matl[0], u, u) / scalar_type(2);
    }


    iso_lin_elasticity_brick(void) {
      set_flags("isotropic linearized elasticity", true /* is linear*/,
                true /* is symmetric */, true /* is coercive */,
                true /* is real */, false /* is complex */);
    }

  };

  size_type add_isotropic_linearized_elasticity_brick
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &dataname1, const std::string &dataname2,
   size_type region, const std::string &dataname3) {
    pbrick pbr = new iso_lin_elasticity_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl(1, dataname1);
    dl.push_back(dataname2);
    if (dataname3.size()) dl.push_back(dataname3);
    return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
                        model::mimlist(1, &mim), region);
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

    GMM_ASSERT1(sl == 1 && sm == 1, "Bad format for Lamé coefficients");
    GMM_ASSERT1(mf_lambda == mf_mu,
                "The two Lamé coefficients should be described on the same "
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
        coeff.resize(mf_p->nb_basic_dof_of_element(cv));
	gmm::copy(gmm::sub_vector(P, gmm::sub_index
                                 (mf_p->ind_basic_dof_of_element(cv))), coeff);
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

    virtual void asm_real_tangent_terms(const model &md, size_type ib,
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

      pNeumann_elem_term pNt = new lin_incomp_Neumann_elem_term
	(md.version_number_of_data_variable( vl[1]), &mf_p,
	 &(md.real_variable(vl[1])), vl[1]);
      md.add_Neumann_term(pNt, vl[0], ib);

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
   const std::string &dataname) {
    pbrick pbr = new linear_incompressibility_brick();
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    if (dataname.size()) {
      dl.push_back(dataname);
      tl.push_back(model::term_description(multname, multname, true));
    }
    md.add_auxilliary_variables_of_Neumann_terms(varname, multname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
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

    virtual scalar_type asm_real_pseudo_potential(const model &md, size_type,
                                                  const model::varnamelist &vl,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::real_matlist &matl,
                                                  model::real_veclist &,
                                                  model::real_veclist &,
                                                  size_type) const {
      const model_real_plain_vector &U = md.real_variable(vl[0]);
      return gmm::vect_sp(matl[0], U, U) / scalar_type(2);
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

    virtual scalar_type asm_complex_pseudo_potential(const model &md,size_type,
                                                  const model::varnamelist &vl,
                                                  const model::varnamelist &,
                                                  const model::mimlist &,
                                                  model::complex_matlist &matl,
                                                  model::complex_veclist &,
                                                  model::complex_veclist &,
                                                  size_type) const {
      const model_complex_plain_vector &U = md.complex_variable(vl[0]);
      return gmm::real(gmm::vect_hp(matl[0], U, U) / scalar_type(2));
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
   const std::string &dataname_rho,  size_type region) {
    pbrick pbr = new mass_brick;
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist dl;
    if (dataname_rho.size())
      dl.push_back(dataname_rho);
    return md.add_brick(pbr, model::varnamelist(1, varname), dl, tl,
                        model::mimlist(1, &mim), region);
  }


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

      if (dl.size() > 4)
        recompute_matrix = recompute_matrix ||
          md.is_var_newer_than_brick(dl[4], ib);

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
      if (dl.size() > 4)
        recompute_matrix = recompute_matrix ||
          md.is_var_newer_than_brick(dl[4], ib);

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

  class theta_method_dispatcher : public virtual_dispatcher {

  public :

    typedef model::build_version build_version;

    void set_dispatch_coeff(const model &md, size_type ib) const {
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


    template <typename MATLIST, typename VECTLIST>
    inline void next_iter(const model &md, size_type ib,
                          const model::varnamelist &/* vl */,
                          const model::varnamelist &/* dl */,
                          MATLIST &/* matl */,
                          VECTLIST &vectl, VECTLIST &vectl_sym,
                          bool first_iter) const {
      if (first_iter) md.update_brick(ib, model::BUILD_RHS);

      // shift the rhs
      transfert(vectl[0], vectl[1]);
      transfert(vectl_sym[0], vectl_sym[1]);

      // add the component represented by the linear matrix terms to the
      // supplementary rhs
      md.linear_brick_add_to_rhs(ib, 1, 0);
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
     std::vector<model::real_veclist> &/* vectl */,
     std::vector<model::real_veclist> &/* vectl_sym */,
     build_version version) const
    { md.brick_call(ib, version, 0); }

    virtual void asm_complex_tangent_terms
    (const model &md, size_type ib, model::complex_matlist &/* matl */,
     std::vector<model::complex_veclist> &/* vectl */,
     std::vector<model::complex_veclist> &/* vectl_sym */,
     build_version version) const
    { md.brick_call(ib, version, 0); }

    theta_method_dispatcher(const std::string &THETA)
      : virtual_dispatcher(2) {
      param_names.push_back(THETA);
    }

  };

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
      //      typedef std::complex<scalar_type> complex_type;

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

