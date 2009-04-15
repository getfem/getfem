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

#include "gmm/gmm_range_basis.h"
#include "getfem/getfem_models.h"
#include "getfem/getfem_assembling.h"


namespace getfem {


  void model::actualize_sizes(void) {
    // mets à jour la taille des variables, les filtres ...
    // à executer avant de lancer le calcul de la matrice tangente.
    
    // doit mettre aussi à jour les paramêtres en adaptant la taille
    // des vecteurs -> warning si la taille change.

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

  // Faire la doc avec les différents types de variables gérées
  //   (fem, mult ...)

  // gérer les dépendances entre l'objet et les mesh_fem, mesh_im
  //    : si un disparaît, l'objet est invalidé ...

  void model::add_fixed_size_variable(const std::string &name, size_type size,
				      size_type niter, bool is_complex) {
    variables[name] = var_description(true, is_complex, false, niter);
    variables[name].set_size(size);
  }
  
  void model::add_fixed_size_constant(const std::string &name, size_type size,
				      size_type niter, bool is_complex) {
    variables[name] = var_description(false, is_complex, false, niter);
    variables[name].set_size(size);
  }

  void model::add_fem_variable(const std::string &name, const mesh_fem &mf,
			       size_type niter, bool is_complex) {
    variables[name] = var_description(true, is_complex, true, niter,
				      VDESCRFILTER_NO, &mf);
    variables[name].set_size(mf.nb_dof());
  }
  
  void model::add_fem_constant(const std::string &name, const mesh_fem &mf,
			       dim_type qdim, size_type niter,
			       bool is_complex) {
    variables[name] = var_description(false, is_complex, true, niter,
				      VDESCRFILTER_NO, &mf, 0, 0, qdim);
    variables[name].set_size(mf.nb_dof()*qdim);
  }

  void model::add_mult_on_region(const std::string &name, const mesh_fem &mf,
				 const mesh_im &mim,
				 const std::string &primal_name, int region,
				 size_type niter, bool is_complex) {
    variables[name] = var_description(true, is_complex, true, niter,
				      VDESCRFILTER_INFSUP, &mf, &mim, region,
				      1, primal_name);
    variables[name].set_size(mf.nb_dof());
  }
  

}  /* end of namespace getfem.                                             */

