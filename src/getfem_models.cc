// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include <iomanip>
#include "gmm/gmm_range_basis.h"
#include "gmm/gmm_solver_cg.h"
#include "getfem/getfem_models.h"
#include "getfem/getfem_assembling.h"
#include "getfem/getfem_derivatives.h"


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
      if (it->second.is_fem_dofs && it->second.filter == VDESCRFILTER_INFSUP) {
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

	// Obtening the coupling matrix between the multipier and
	// the primal variable. A search is done on all the terms of the
	// model. The corresponding terms are added. If no term is available
	// a warning is printed and the variable is cancelled.

	gmm::col_matrix< gmm::rsvector<scalar_type> >
	  MM(it2->second.mf->nb_dof(), it->second.mf->nb_dof());
	bool termadded = false;

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
	  GMM_WARNING1("No term present to filter the multiplier " << mults[k]
		       << ". The multiplier is cancelled.");

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
	range_basis(MGLOB, glob_columns, 1E-12, gmm::col_major(), true);

	s = 0;
	for (size_type k = 0; k < mults.size(); ++k) {
	  VAR_SET::iterator it = variables.find(mults[k]);
	  dal::bit_vector kept;
	  size_type nbdof = it->second.mf->nb_dof();
	  for (std::set<size_type>::iterator itt = glob_columns.begin();
	       itt != glob_columns.end(); ++itt)
	    if (*itt > s && *itt < s + nbdof) kept.add(*itt-s);
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
	size_type d = sizeof(scalar_type);
	if (is_complex()) d *= 2;
	ost << std::setw(8) << std::right << it->second.size() * d
	    << " bytes.";
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
				      VDESCRFILTER_INFSUP, &mf, 0,
				      1, primal_name);
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


  void model::add_time_dispatcher(size_type ibrick, pdispatcher pdispatch) {
        
    GMM_ASSERT1(ibrick < bricks.size(), "Unexistent brick");

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
    GMM_ASSERT1(ind_brick < bricks.size(), "Unexistent brick");
    GMM_ASSERT1(ind_var < bricks[ind_brick].vlist.size(),
	       "Inexistent brick variable");
    return bricks[ind_brick].vlist[ind_var];
  }
  
  const std::string &model::dataname_of_brick(size_type ind_brick,
					      size_type ind_data) {
    GMM_ASSERT1(ind_brick < bricks.size(), "Unexistent brick");
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
      size_type nbd1 = variables[term.var1].size();
      size_type nbd2 = term.is_matrix_term ?
	variables[term.var2].size() : 0;
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
		  "Unexistent iteration " << ind);
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

	size_type n_iter_1 = n_iter, n_iter_2 = n_iter;
	if (n_iter == size_type(-1)) {
	  n_iter_1 = variables[term.var1].default_iter;
	  if (term.is_matrix_term)
	    n_iter_2 = variables[term.var2].default_iter;
	}
	
	if (term.is_matrix_term) {
	  if (cplx) 
	    gmm::mult_add
	      (brick.cmatlist[j], 
	       gmm::scaled(variables[term.var2].complex_value[n_iter_2],
			   std::complex<scalar_type>(-1)),
	       brick.cveclist[ind_data][j]);
	  else
	    gmm::mult_add
	      (brick.rmatlist[j],
	       gmm::scaled(variables[term.var2].real_value[n_iter_2],
			   scalar_type(-1)), brick.rveclist[ind_data][j]);
	  
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


  void model::assembly(build_version version) {

    context_check(); if (act_size_to_be_done) actualize_sizes();
    if (is_complex()) { gmm::clear(cTM); gmm::clear(crhs); }
    else { gmm::clear(rTM); gmm::clear(rrhs); }

    for (dal::bv_visitor ib(active_bricks); !ib.finished(); ++ib) {
      brick_description &brick = bricks[ib];

      update_brick(ib, version);

      bool cplx = is_complex() && brick.pbr->is_complex();
      
      // Assembly of terms
      for (size_type j = 0; j < brick.tlist.size(); ++j) {
	term_description &term = brick.tlist[j];
	gmm::sub_interval I1 = variables[term.var1].I;
	gmm::sub_interval I2(0,0);
	if (term.is_matrix_term) I2 = variables[term.var2].I;

	scalar_type coeff0 = scalar_type(1);
	if (brick.pdispatch) coeff0 = brick.matrix_coeff;

	if (cplx) {
	  if (term.is_matrix_term && (version | BUILD_MATRIX)) {
	    gmm::add(gmm::scaled(brick.cmatlist[j], coeff0),
		     gmm::sub_matrix(cTM, I1, I2));
	    if (term.is_symmetric && I1.first() != I2.first()) {
	      gmm::add(gmm::scaled(gmm::transposed(brick.cmatlist[j]), coeff0),
		       gmm::sub_matrix(cTM, I2, I1));
	    }
	  }
	  if (version | BUILD_RHS) {
	    if (brick.pdispatch) {
	      for (size_type k = 0; k < brick.nbrhs; ++k)
		gmm::add(gmm::scaled(brick.cveclist[k][j],
				     brick.coeffs[k]),
			 gmm::sub_vector(crhs, I1));
	    }
	    else
	      gmm::add(brick.cveclist[0][j], gmm::sub_vector(crhs, I1));
	    if (term.is_matrix_term && brick.pbr->is_linear()
		&& (!is_linear() || version & BUILD_WITH_COMPLETE_RHS)) {
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
		   && (!is_linear() || version & BUILD_WITH_COMPLETE_RHS)) {
		 gmm::mult_add(gmm::conjugated(brick.cmatlist[j]),
			    gmm::scaled(variables[term.var1].complex_value[0],
					std::complex<scalar_type>(-coeff0)),
			       gmm::sub_vector(crhs, I2));
	       }
	    }
	  }
	} else if (is_complex()) {
	  if (term.is_matrix_term && (version | BUILD_MATRIX)) {
	    gmm::add(gmm::scaled(brick.rmatlist[j], coeff0),
		     gmm::sub_matrix(cTM, I1, I2));
	    if (term.is_symmetric && I1.first() != I2.first()) {
	      gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), coeff0),
		       gmm::sub_matrix(cTM, I2, I1));
	    }
	  }
	  if (version | BUILD_RHS) {
	    if (brick.pdispatch) {
	      for (size_type k = 0; k < brick.nbrhs; ++k)
		gmm::add(gmm::scaled(brick.rveclist[k][j],
				     brick.coeffs[k]),
			 gmm::sub_vector(crhs, I1));
	    }
	    else
	      gmm::add(brick.rveclist[0][j], gmm::sub_vector(crhs, I1));
	    if (term.is_matrix_term && brick.pbr->is_linear()
		&& (!is_linear() || version & BUILD_WITH_COMPLETE_RHS)) {
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
		  && (!is_linear() || version & BUILD_WITH_COMPLETE_RHS)) {
		gmm::mult_add(gmm::transposed(brick.rmatlist[j]),
			     gmm::scaled(variables[term.var1].complex_value[0],
					  std::complex<scalar_type>(-coeff0)),
			      gmm::sub_vector(crhs, I2));
	      }
	    }
	  }
	} else {
	  if (term.is_matrix_term && (version | BUILD_MATRIX)) {
	    gmm::add(gmm::scaled(brick.rmatlist[j], coeff0),
		     gmm::sub_matrix(rTM, I1, I2));
	    if (term.is_symmetric && I1.first() != I2.first()) {
	      gmm::add(gmm::scaled(gmm::transposed(brick.rmatlist[j]), coeff0),
		       gmm::sub_matrix(rTM, I2, I1));
	    }
	  }
	  if (version | BUILD_RHS) {
	    if (brick.pdispatch) {
	      for (size_type k = 0; k < brick.nbrhs; ++k)
		gmm::add(gmm::scaled(brick.rveclist[k][j],
				     brick.coeffs[k]),
			 gmm::sub_vector(rrhs, I1));
	    }
	    else
	      gmm::add(brick.rveclist[0][j], gmm::sub_vector(rrhs, I1));
	    if (term.is_matrix_term && brick.pbr->is_linear()
		&& (!is_linear() || version & BUILD_WITH_COMPLETE_RHS)) {
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
		  && (!is_linear() || version & BUILD_WITH_COMPLETE_RHS)) {
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
      else 
	if (cplx) {
	  brick.cmatlist = complex_matlist(brick.tlist.size());
	  brick.cveclist[0] = complex_veclist(brick.tlist.size());
	} else {
	  brick.rmatlist = real_matlist(brick.tlist.size());
	  brick.rveclist[0] = real_veclist(brick.tlist.size());	    
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
		"Unvalid iteration number "
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
		"Unvalid iteration number "
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
		"Unvalid iteration number "
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
		"Unvalid iteration number "
		<< niter << " for " << name);
    return it->second.complex_value[niter];    
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
  // Generic elliptic brick
  //
  // ----------------------------------------------------------------------

  struct generic_elliptic_brick : public virtual_brick {

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
	GMM_ASSERT1(false,
		    "Bad format generic elliptic brick coefficient");
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

    source_term_brick(void) {
      set_flags("Source term", true /* is linear*/,
		true /* is symmetric */, true /* is coercive */,
		true /* is real */, true /* is complex */);
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

    normal_source_term_brick(void) {
      set_flags("Normal source term", true /* is linear*/,
		true /* is symmetric */, true /* is coercive */,
		true /* is real */, true /* is complex */);
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
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 2,
		  "Wrong number of variables for Dirichlet condition brick");

      bool penalized = (vl.size() == 1);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_mult = md.mesh_fem_of_variable(vl[vl.size()-1]);
      const mesh_im &mim = *mims[0];
      const model_real_plain_vector *A = 0, *COEFF = 0;
      const mesh_fem *mf_data = 0;
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
	GMM_ASSERT1(mf_u.get_qdim() == s,
		    dl[ind] << ": bad format of Dirichlet data. "
		    "Detected dimension is " << s << " should be "
		    << size_type(mf_u.get_qdim()));
      }
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (dl.size()) {
	GMM_TRACE2("Source term assembly for Dirichlet condition");	  
	if (mf_data)
	  asm_source_term(vecl[0], mim, mf_mult, *mf_data, *A, rg);
	else
	  asm_homogeneous_source_term(vecl[0], mim, mf_mult, *A, rg);
	if (penalized) gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }

      if (recompute_matrix) {
	GMM_TRACE2("Mass term assembly for Dirichlet condition");
	gmm::clear(matl[0]);
	asm_mass_matrix(matl[0], mim, mf_mult, mf_u, region);
	if (penalized) gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
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
      GMM_ASSERT1(vl.size() >= 1 && vl.size() <= 2 && dl.size() <= 2,
		  "Wrong number of variables for Dirichlet condition brick");

      bool penalized = (vl.size() == 1);
      const mesh_fem &mf_u = md.mesh_fem_of_variable(vl[0]);
      const mesh_fem &mf_mult = md.mesh_fem_of_variable(vl[vl.size()-1]);
      const mesh_im &mim = *mims[0];
      const model_complex_plain_vector *A = 0, *COEFF = 0;
      const mesh_fem *mf_data = 0;
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
	GMM_ASSERT1(mf_u.get_qdim() == s, "Bad format of Dirichlet data");
      }
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (dl.size() > ind) {
	GMM_TRACE2("Source term assembly for Dirichlet condition");
	if (mf_data)
	  asm_source_term(vecl[0], mim, mf_mult, *mf_data, *A, rg);
	else
	  asm_homogeneous_source_term(vecl[0], mim, mf_mult, *A, rg);
	if (penalized) gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }

      if (recompute_matrix) {
	GMM_TRACE2("Mass term assembly for Dirichlet condition");
	gmm::clear(matl[0]);
	asm_mass_matrix(matl[0], mim, mf_mult, mf_u, region);
	if (penalized) gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
      }
    }

    Dirichlet_condition_brick(bool penalized) {
      set_flags(penalized ? "Dirichlet with penalization brick"
		          : "Dirichlet with multipliers brick",
		true /* is linear*/,
		true /* is symmetric */, penalized /* is coercive */,
		true /* is real */, true /* is complex */);
    }
  };

  size_type add_Dirichlet_condition_with_multipliers
  (model &md, const mesh_im &mim, const std::string &varname,
   const std::string &multname, size_type region,
   const std::string &dataname) {
    pbrick pbr = new Dirichlet_condition_brick(false);
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
   const std::string &dataname) {
    std::string coeffname = md.new_name("penalization_on_" + varname);
    md.add_fixed_size_data(coeffname, 1);
    if (md.is_complex())
      md.set_complex_variable(coeffname)[0] = penalisation_coeff;
    else
      md.set_real_variable(coeffname)[0] = penalisation_coeff;
    pbrick pbr = new Dirichlet_condition_brick(true);
    model::termlist tl;
    tl.push_back(model::term_description(varname, varname, true));
    model::varnamelist vl(1, varname);
    model::varnamelist dl(1, coeffname);
    if (dataname.size()) dl.push_back(dataname);
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
		true /* is real */, true /* is complex */);
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
      }

      if (penalized) {
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
      }

      if (penalized) {
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
		true /* is real */, true /* is complex */);
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
		true /* is real */, true /* is complex */);
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
		true /* is real */, true /* is complex */);
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

  struct linear_incompressibility_brick : public virtual_brick {
    
    virtual void asm_real_tangent_terms(const model &md, size_type,
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
      const model_real_plain_vector *A = 0, *COEFF = 0;
      const mesh_fem *mf_data = 0;

      if (penalized) {
	COEFF = &(md.real_variable(dl[0]));
	mf_data = md.pmesh_fem_of_variable(dl[0]);
	size_type s = gmm::vect_size(*A);
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
	  asm_mass_matrix_param(matl[1], mim, mf_p, *mf_data, *A, rg);
	  gmm::scale(matl[1], scalar_type(-1));
	}
	else {
	  asm_mass_matrix(matl[1], mim, mf_p, rg);
	  gmm::scale(matl[1], -(*A)[0]);
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
		true /* is real */, true /* is complex */);
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
		true /* is real */, true /* is complex */);
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
		true /* is real */, true /* is complex */);
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
      typedef std::complex<scalar_type> complex_type;
      
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
    md.unable_brick(id2dt2b);
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

      pbrick pbr = md.get_brick(ib);

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
      pbrick pbr = md.get_brick(ib);
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
      pbrick pbr = md.get_brick(ib);
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

