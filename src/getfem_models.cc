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
#include "getfem/getfem_models.h"
#include "getfem/getfem_assembling.h"


namespace getfem {

  void model::var_description::set_size(size_type s) {
    if (is_complex)
      complex_value.resize(n_iter);
    else
      real_value.resize(n_iter);
    for (size_type i = 0; i < n_iter; ++i)
      if (is_complex)
	complex_value[i].resize(s);
      else
	real_value[i].resize(s);
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
    GMM_ASSERT1(!assert || valid, "Illegal variable name : " << name);
    return valid;
  }

  void model::actualize_sizes(void) {

    std::map<std::string, std::vector<std::string> > multipliers;
    std::map<std::string, bool > tobedone;

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
	case VDESCRFILTER_INFSUP:
	  VAR_SET::iterator it2 = variables.find(it->second.filter_var);
	  GMM_ASSERT1(it2 != variables.end(), "The primal variable of the "
		      "multiplier does not exist");
	  multipliers[it->second.filter_var].push_back(it->first);
	  if (it->second.v_num < it->second.mf->version_number() ||
	      it->second.v_num < it2->second.mf->version_number())
	    tobedone[it->second.filter_var] = true;
	  break;
	}
      }
    }

    for (std::map<std::string, bool >::iterator itbd = tobedone.begin();
	 itbd != tobedone.end(); ++itbd) {
      std::vector<std::string> &mults = multipliers[itbd->first];
      VAR_SET::iterator it2 = variables.find(itbd->first);
      GMM_ASSERT1(it2->second.is_fem_dofs, "The primal variable of the "
		  "multiplier is not a fem variable");

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
	
	gmm::col_matrix< gmm::rsvector<scalar_type> >
	  MM(it2->second.mf->nb_dof(), it->second.mf->nb_dof());
	asm_mass_matrix(MM, *(it->second.mim), *(it2->second.mf),
			*(it->second.mf), it->second.m_region);
	std::set<size_type> columns;
	gmm::range_basis(MM, columns);
	if (mults.size() > 1) {
	  gmm::copy(MM, gmm::sub_matrix
		    (MGLOB,gmm::sub_interval(0, it2->second.mf->nb_dof()),
		     gmm::sub_interval(s, it->second.mf->nb_dof())));
	  s += it->second.mf->nb_dof();
	  for (std::set<size_type>::iterator itt = columns.begin();
	     itt != columns.end(); ++itt)
	    glob_columns.insert(s + *itt);
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
      act_size_to_be_done = false;
    }

    size_type tot_size = 0;

    for (VAR_SET::iterator it = variables.begin(); it != variables.end();
	 ++it) {
      if (it->second.is_variable) {
	it->second.I = gmm::sub_interval(tot_size, it->second.size());
	tot_size += it->second.size();
      }
    }
      
    if (complex_version) {
      gmm::resize(cTM, tot_size, tot_size);
      gmm::resize(cstate, tot_size);
	gmm::resize(crhs, tot_size);
    }
    else {
      gmm::resize(rTM, tot_size, tot_size);
      gmm::resize(rstate, tot_size);
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

  std::string model::new_name(const std::string &name) {
    std::string res_name = name;
    check_name_valitity(res_name);
    VAR_SET::const_iterator it = variables.find(res_name);
    for (size_type i = 2; it != variables.end(); ++i) {
      std::stringstream m;
      m << name << '_' << i;
      res_name = m.str();
      it = variables.find(res_name);
    }
    return res_name;
  }


  void model::add_fixed_size_variable(const std::string &name, size_type size,
				      size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(true, is_complex(), false, niter);
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
    this->add_dependency(mf);
  }
  
  void model::add_fem_data(const std::string &name, const mesh_fem &mf,
			       dim_type qdim, size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(false, is_complex(), true, niter,
				      VDESCRFILTER_NO, &mf, 0, 0, qdim);
    variables[name].set_size(mf.nb_dof()*qdim);
    this->add_dependency(mf); 
  }

  void model::add_mult_on_region(const std::string &name, const mesh_fem &mf,
				 const mesh_im &mim,
				 const std::string &primal_name,
				 size_type region,
				 size_type niter) {
    check_name_valitity(name);
    variables[name] = var_description(true, is_complex(), true, niter,
				      VDESCRFILTER_INFSUP, &mf, &mim, region,
				      1, primal_name);
    variables[name].set_size(mf.nb_dof());
    act_size_to_be_done = true;
    this->add_dependency(mf); this->add_dependency(mim);
  }

  size_type model::add_brick(pbrick pbr, const varnamelist &varnames,
			     const varnamelist &datanames,
			     const termlist &terms,
			     const mimlist &mims, size_type region) {
    bricks.push_back(brick_description(pbr, varnames, datanames, terms,
				       mims, region));
    for  (size_type i = 0; i < bricks.back().mims.size(); ++i)
      add_dependency(*(bricks.back().mims[i]));

    GMM_ASSERT1(pbr->is_real() || is_complex(),
		"Impossible to add a complex brick to a real model");
    if (is_complex() && pbr->is_complex()) {
      bricks.back().cmatlist.resize(terms.size());
      bricks.back().cveclist.resize(terms.size());
    } else {
      bricks.back().rmatlist.resize(terms.size());
      bricks.back().rveclist.resize(terms.size());
    }
    is_linear = is_linear && pbr->is_linear();
    is_symmetric = is_symmetric && pbr->is_symmetric();
    is_coercive = is_coercive && pbr->is_coercive();

    for (size_type i=0; i < varnames.size(); ++i)
      GMM_ASSERT1(variables.find(varnames[i]) != variables.end(),
		  "Undefined model variable " << varnames[i]);
    for (size_type i=0; i < datanames.size(); ++i)
      GMM_ASSERT1(variables.find(datanames[i]) != variables.end(),
		  "Undefined model data or variable " << datanames[i]);
    
    return size_type(bricks.size() - 1);
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
    cout << "ind_brick = " << ind_brick << endl;
    cout << "ind_data = " << ind_data << " bricks[ind_brick].dlist.size() = " << bricks[ind_brick].dlist.size() << endl;
    GMM_ASSERT1(ind_brick < bricks.size(), "Inexistent brick");
    GMM_ASSERT1(ind_data < bricks[ind_brick].dlist.size(),
		"Inexistent brick data");
    return bricks[ind_brick].dlist[ind_data];
  }

  void model::listbricks(std::ostream &ost) const {
    if (bricks.size() == 0)
      ost << "Model with no bricks" << endl;
    else {
      ost << "List of model bricks:" << endl;
      for (size_type i = 0; i < bricks.size(); ++i) {
	ost << "Brick " << std::setw(3) << std::right << i
	    << " " << std::setw(20) << std::right
	    << bricks[i].pbr->brick_name() << endl;
	ost << "  concerned variables: " << bricks[i].vlist[0];
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

  void model::assembly(void) {
    context_check(); if (act_size_to_be_done) actualize_sizes();
    if (is_complex()) { gmm::clear(cTM); gmm::clear(crhs); }
    else { gmm::clear(rTM); gmm::clear(rrhs); }

    for (size_type ib = 0; ib < bricks.size(); ++ib) {
      brick_description &brick = bricks[ib];

      bool cplx = is_complex() && brick.pbr->is_complex();
      
      bool tobecomputed = brick.terms_to_be_computed
	|| !(brick.pbr->is_linear());
      
      // check variable list to test if a mesh_fem as changed. 
      for (size_type i = 0; i < brick.vlist.size() && !tobecomputed; ++i) {
	var_description &vd = variables[brick.vlist[i]];
	if (vd.v_num > brick.v_num)
	  tobecomputed = true;
      }
      
      // check data list to test if a vector value of a data has changed. 
      for (size_type i = 0; i < brick.dlist.size() && !tobecomputed; ++i) {
	var_description &vd = variables[brick.dlist[i]];
	if (vd.v_num > brick.v_num || vd.v_num_data > brick.v_num)
	  tobecomputed = true;
      }
      
      if (tobecomputed) {
	// Initialization of vector and matrices.
	for (size_type j = 0; j < brick.tlist.size(); ++j) {
	  term_description &term = brick.tlist[j];
	  size_type nbd1 = variables[term.var1].size();
	  size_type nbd2 = term.is_matrix_term ?
	    variables[term.var2].size() : 0;
	  if (term.is_matrix_term) {
	    if (cplx)
	      brick.cmatlist[j] = model_complex_sparse_matrix(nbd1, nbd2);
	    else
	      brick.rmatlist[j] = model_real_sparse_matrix(nbd1, nbd2);
	  }
	  if (!(term.is_matrix_term) || !(brick.pbr->is_linear())) {
	    if (cplx) {
	      gmm::clear(brick.cveclist[j]);
	      gmm::resize(brick.cveclist[j], nbd1);
	    } else {
	      gmm::clear(brick.rveclist[j]);
	      gmm::resize(brick.rveclist[j], nbd1);
	    }
	  }
	  brick.v_num = act_counter();
	}
	 
	// Brick call for all terms.
	if (cplx)
	  brick.pbr->asm_complex_tangent_terms(*this, brick.vlist, brick.dlist,
					       brick.mims,
					       brick.cmatlist, brick.cveclist,
					       brick.region);
	else
	  brick.pbr->asm_real_tangent_terms(*this, brick.vlist, brick.dlist,
					    brick.mims,
					    brick.rmatlist, brick.rveclist,
					    brick.region);

	if (brick.pbr->is_linear())
	  brick.terms_to_be_computed = false;
	else 
	  if (cplx) {
	    brick.cmatlist = complex_matlist(brick.tlist.size());
	    brick.cveclist = complex_veclist(brick.tlist.size());
	  } else {
	    brick.rmatlist = real_matlist(brick.tlist.size());
	    brick.rveclist = real_veclist(brick.tlist.size());	    
	  }
      }

      // Assembly of terms
      for (size_type j = 0; j < brick.tlist.size(); ++j) {
	term_description &term = brick.tlist[j];
	gmm::sub_interval I1 = variables[term.var1].I;
	gmm::sub_interval I2(0,0);
	if (term.is_matrix_term) I2 = variables[term.var2].I;

	if (cplx) {
	  if (term.is_matrix_term) {
	    gmm::add(brick.cmatlist[j], gmm::sub_matrix(cTM, I1, I2));
	    if (brick.pbr->is_linear() && !is_linear) {
	      gmm::mult(brick.cmatlist[j],
			gmm::scaled(variables[term.var1].complex_value[0],
				    std::complex<scalar_type>(-1)),
			gmm::sub_vector(crhs, I1));
	    }
	    if (term.is_symmetric && I1.first() != I2.first()) {
	      gmm::add(gmm::transposed(brick.cmatlist[j]),
		       gmm::sub_matrix(cTM, I2, I1));
	      if (brick.pbr->is_linear() && !is_linear) {
		gmm::mult(gmm::conjugated(brick.cmatlist[j]),
			  gmm::scaled(variables[term.var2].complex_value[0],
				      std::complex<scalar_type>(-1)),
			  gmm::sub_vector(crhs, I2));
	      }
	    }
	  }
	  // if (!(term.is_matrix_term) || !(brick.pbr->is_linear()))
	  gmm::add(brick.cveclist[j], gmm::sub_vector(crhs, I1));
	} else if (is_complex()) {
	  if (term.is_matrix_term) {
	    gmm::add(brick.rmatlist[j], gmm::sub_matrix(cTM, I1, I2));
	    if (brick.pbr->is_linear() && !is_linear) {
	      gmm::mult(brick.rmatlist[j],
			gmm::scaled(variables[term.var1].real_value[0],
				    scalar_type(-1)),
			gmm::sub_vector(crhs, I1));
	    }
	    if (term.is_symmetric && I1.first() != I2.first()) {
	      gmm::add(gmm::transposed(brick.rmatlist[j]),
		       gmm::sub_matrix(cTM, I2, I1));
	      if (brick.pbr->is_linear() && !is_linear) {
		gmm::mult(gmm::transposed(brick.rmatlist[j]),
			  gmm::scaled(variables[term.var2].real_value[0],
				      scalar_type(-1)),
			  gmm::sub_vector(crhs, I2));
	      }
	    }
	  }
	  // if (!(term.is_matrix_term) || !(brick.pbr->is_linear()))
	  gmm::add(brick.rveclist[j], gmm::sub_vector(crhs, I1));
	} else {
	  if (term.is_matrix_term) {
	    gmm::add(brick.rmatlist[j], gmm::sub_matrix(rTM, I1, I2));
	    if (brick.pbr->is_linear() && !is_linear) {
	      gmm::mult(brick.rmatlist[j],
			gmm::scaled(variables[term.var1].real_value[0],
				    scalar_type(-1)),
			gmm::sub_vector(rrhs, I1));
	    }
	    if (term.is_symmetric && I1.first() != I2.first()) {
	      gmm::add(gmm::transposed(brick.rmatlist[j]),
		       gmm::sub_matrix(rTM, I2, I1));
	      if (brick.pbr->is_linear() && !is_linear) {
		gmm::mult(gmm::transposed(brick.rmatlist[j]),
			  gmm::scaled(variables[term.var2].real_value[0],
				      scalar_type(-1)),
			  gmm::sub_vector(rrhs, I2));
	      }
	    }
	  }
	  // if (!(term.is_matrix_term) || !(brick.pbr->is_linear()))
	  gmm::add(brick.rveclist[j], gmm::sub_vector(rrhs, I1));
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
    // context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::const_iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    GMM_ASSERT1(it->second.n_iter > niter, "Unvalid iteration number");
    return it->second.real_value[niter];
  }
  
  const model_complex_plain_vector &
  model::complex_variable(const std::string &name, size_type niter) const {
    GMM_ASSERT1(complex_version, "This model is a real one");
    // context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::const_iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    GMM_ASSERT1(it->second.n_iter > niter, "Unvalid iteration number");
    return it->second.complex_value[niter];    
  }

  model_real_plain_vector &
  model::set_real_variable(const std::string &name, size_type niter) {
    GMM_ASSERT1(!complex_version, "This model is a complex one");
    context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    GMM_ASSERT1(it->second.n_iter > niter, "Unvalid iteration number");
    it->second.v_num_data = act_counter();
    return it->second.real_value[niter];
  }
  
  model_complex_plain_vector &
  model::set_complex_variable(const std::string &name, size_type niter) {
    GMM_ASSERT1(complex_version, "This model is a real one");
    context_check(); if (act_size_to_be_done) actualize_sizes();
    VAR_SET::iterator it = variables.find(name);
    GMM_ASSERT1(it!=variables.end(), "Undefined variable " << name);
    GMM_ASSERT1(it->second.n_iter > niter, "Unvalid iteration number");
    it->second.v_num_data = act_counter();    
    return it->second.complex_value[niter];    
  }


  // ----------------------------------------------------------------------
  //
  // Generic scalar elliptic brick
  //
  // ----------------------------------------------------------------------

  struct generic_elliptic_brick : public virtual_brick {

    virtual void asm_real_tangent_terms(const model &md,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &,
					size_type region) const {
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

    virtual void asm_complex_tangent_terms(const model &md,
					   const model::varnamelist &vl,
					   const model::varnamelist &dl,
					   const model::mimlist &mims,
					   model::complex_matlist &matl,
					   model::complex_veclist &,
					   size_type region) const {
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

    virtual void asm_real_tangent_terms(const model &md,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &,
					model::real_veclist &vecl,
					size_type region) const {
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

      size_type s = gmm::vect_size(A);
      if (mf_data) s = s * mf_data->get_qdim() / mf_data->nb_dof();

      GMM_ASSERT1(mf_u.get_qdim() == s, "Bad format of source term data");

      if (mf_data)
	asm_source_term(vecl[0], mim, mf_u, *mf_data, A, rg);
      else
	asm_homogeneous_source_term(vecl[0], mim, mf_u, A, rg);

      if (dl.size() > 1) gmm::add(md.real_variable(dl[1]), vecl[0]);

    }

    virtual void asm_complex_tangent_terms(const model &md,
					   const model::varnamelist &vl,
					   const model::varnamelist &dl,
					   const model::mimlist &mims,
					   model::complex_matlist &,
					   model::complex_veclist &vecl,
					   size_type region) const {
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
  // Dirichlet condition brick
  //
  // ----------------------------------------------------------------------
  // Two variables : with multipliers
  // One variable : penalization

  struct Dirichlet_condition_brick : public virtual_brick {
    
    virtual void asm_real_tangent_terms(const model &md,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &vecl,
					size_type region) const {
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
	GMM_ASSERT1(mf_u.get_qdim() == s, "Bad format of Dirichlet data");
      }
      mesh_region rg(region);
      mim.linked_mesh().intersect_with_mpi_region(rg);

      if (dl.size()) {
	if (mf_data)
	  asm_source_term(vecl[0], mim, mf_mult, *mf_data, *A, rg);
	else
	  asm_homogeneous_source_term(vecl[0], mim, mf_mult, *A, rg);
	if (penalized) gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }

      asm_mass_matrix(matl[0], mim, mf_mult, mf_u, region);
      if (penalized) gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
    }

    virtual void asm_complex_tangent_terms(const model &md,
					   const model::varnamelist &vl,
					   const model::varnamelist &dl,
					   const model::mimlist &mims,
					   model::complex_matlist &matl,
					   model::complex_veclist &vecl,
					   size_type region) const {
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
	if (mf_data)
	  asm_source_term(vecl[0], mim, mf_mult, *mf_data, *A, rg);
	else
	  asm_homogeneous_source_term(vecl[0], mim, mf_mult, *A, rg);
	if (penalized) gmm::scale(vecl[0], gmm::abs((*COEFF)[0]));
      }

      asm_mass_matrix(matl[0], mim, mf_mult, mf_u, region);
      if (penalized) gmm::scale(matl[0], gmm::abs((*COEFF)[0]));
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
    md.add_mult_on_region(multname, mf_mult, mim, varname, region);
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
   size_type region, scalar_type penalisation_coeff,
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

  void change_penalization_coeff_Dirichlet(model &md, size_type ind_brick,
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
  // Constraint brick
  //
  // ----------------------------------------------------------------------

  struct constraint_brick : public virtual_brick {

    model_real_sparse_matrix rB;
    model_complex_sparse_matrix cB;
    model_real_plain_vector rL;
    model_complex_plain_vector cL;
    
    virtual void asm_real_tangent_terms(const model &md,
					const model::varnamelist &vl,
					const model::varnamelist &dl,
					const model::mimlist &mims,
					model::real_matlist &matl,
					model::real_veclist &vecl,
					size_type) const {
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

    virtual void asm_complex_tangent_terms(const model &md,
					   const model::varnamelist &vl,
					   const model::varnamelist &dl,
					   const model::mimlist &mims,
					   model::complex_matlist &matl,
					   model::complex_veclist &vecl,
					   size_type) const {
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
		          : "Constraint with multiplier brick",
		true /* is linear*/,
		true /* is symmetric */, penalized /* is coercive */,
		true /* is real */, true /* is complex */);
    }


  };

  // à interfacer, dont l'acces aux contraintes stockées dans la brique ...
  // et a ranger dans une serie de briques "accès direct" ou on peut modifier
  // les sys linéaire sans intégration de m.e.f.




}  /* end of namespace getfem.                                             */

