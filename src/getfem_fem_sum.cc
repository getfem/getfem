// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : GEneric Tool for Finite Element Methods (getfem)
// File    : getfem_fem_sum.cc : make the direct sum of a set of finite
//           element method.
// Date    : October 29, 2004.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2004-2005 Yves Renard
//
// This file is a part of GETFEM++
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
//========================================================================

#include <getfem_fem_sum.h>

namespace getfem {
  namespace obsolete {

  void fem_sum::update_from_context(void) const {
    dim_ = dim_type(-1);
    size_type max_dof = 0;
    const getfem_mesh *pmesh = &(mfs[0]->linked_mesh());
    for (dal::bv_visitor cv(pmesh->convex_index()); !cv.finished(); ++cv) {
      if (dim_ == dim_type(-1))
	dim_ = pmesh->structure_of_convex(cv)->dim();
      if (dim_ != pmesh->structure_of_convex(cv)->dim())
	DAL_THROW(dal::failure_error, "Convexes of different dimension"
		  ": to be done");
      max_dof = std::max(max_dof, nb_dof(cv));
    }
    std::fill(spf.begin(), spf.end(), pfem(0));
    if (max_dof == 0) return;
    
    base_node P(dim());
    std::fill(P.begin(), P.end(), scalar_type(1)/scalar_type(20));
    node_tab_.resize(max_dof);
    std::fill(node_tab_.begin(), node_tab_.end(), P);
    pspt_valid = false;
    dof_types_.resize(max_dof);
    std::fill(dof_types_.begin(), dof_types_.end(),
	      global_dof(dim()));
    
    vit.resize(target_dim()*dim());
  }
  
  void fem_sum::init(const std::vector<const mesh_fem *> &mfs_) {
    mfs = mfs_;

    const getfem_mesh *pmesh = &(mfs[0]->linked_mesh());
    for (size_type i = 0; i < mfs.size(); ++i) {
      if (&(mfs[i]->linked_mesh()) != pmesh)
	DAL_THROW(failure_error, "Meshes should be the same");
      if (mfs[i]->get_qdim() != 1) 
	DAL_THROW(dal::to_be_done_error, "fem_sum do not handle qdim != 1");
      this->add_dependency(*(mfs[i]));
    }
    is_pol = is_lag = false; es_degree = 5;
    is_equiv = real_element_defined = true;
    ntarget_dim = 1; // An extension for vectorial elements should be easy
    // The detection should be done and the multiplication of components
    // for scalar elements interpolated.
    mi2 = bgeot::multi_index(2);
    mi3 = bgeot::multi_index(3);
    update_from_context();
    spf.resize(mfs.size());
    spspt.resize(mfs.size());
    spfp.resize(mfs.size());
    taux.resize(mfs.size());
  }

  size_type fem_sum::nb_dof(size_type cv) const {
    context_check();
    size_type res = 0;
    for (size_type i = 0; i < mfs.size(); ++i)
      if (mfs[i]->convex_index()[cv])
	res += mfs[i]->nb_dof_of_element(cv);
    return res;
  }

  size_type fem_sum::index_of_global_dof(size_type cv,
						 size_type j) const {
    size_type res = 0, i;
    for (i = 0; i < mfs.size(); ++i) {
      size_type nb
	= (mfs[i]->convex_index()[cv]) ? mfs[i]->nb_dof_of_element(cv) : 0;
      if (j < nb) break;
      j -= nb;
      res += mfs[i]->nb_dof();
    }
    return res+mfs[i]->ind_dof_of_element(cv)[j];
  }
  
  bgeot::pconvex_ref fem_sum::ref_convex(size_type cv) const {
    for (size_type i = 0; i < mfs.size(); ++i)
      if (mfs[i]->convex_index()[cv])
	return mfs[i]->fem_of_element(cv)->ref_convex(cv);
    DAL_THROW(internal_error, "Internal error");
  }
  
  const bgeot::convex<base_node> &fem_sum::node_convex(size_type cv) const
  { return *(bgeot::generic_dummy_convex_ref(dim(), nb_dof(cv), mfs[0]->linked_mesh().structure_of_convex(cv)->nb_faces())); }
  
  bgeot::pstored_point_tab fem_sum::node_tab(size_type) const { 
    if (!pspt_valid)
      { pspt = bgeot::store_point_tab(node_tab_); pspt_valid = true; }
    return pspt;
  }
  
  void fem_sum::base_value(const base_node &, base_tensor &) const
  { DAL_THROW(internal_error, "No base values, real only element."); }
  
  void fem_sum::grad_base_value(const base_node &, base_tensor &) const
  { DAL_THROW(internal_error, "No grad values, real only element."); }
  
  void fem_sum::hess_base_value(const base_node &, base_tensor &) const
  { DAL_THROW(internal_error, "No hess values, real only element."); }
  
  pfem_precomp fem_sum::get_pfp(size_type i,pfem pf,
				bgeot::pstored_point_tab pspti) const {
    if (pf != spf[i] || pspti != spspt[i])
	{ spfp[i]=fem_precomp(pf, pspti); spf[i] = pf; spspt[i] = pspti; }
    return spfp[i];
  }
  
  void fem_sum::real_base_value(const fem_interpolation_context& d, 
				base_tensor &t) const {
    fem_interpolation_context& c = const_cast<fem_interpolation_context&>(d);
    size_type cv = c.convex_num(), nbdof = nb_dof(cv);
    mi2[1] = target_dim(); mi2[0] = nbdof;
    t.adjust_sizes(mi2);
    bool have_pfp = c.have_pfp();
    pfem_precomp pfp = c.pfp();
    pfem pf = c.pf();
    
    base_tensor::iterator it1 = t.begin();
    for (size_type k = 0; k < target_dim(); ++k)
      { vit[k] = it1; it1 += nbdof; }
    
    for (size_type i = 0; i < mfs.size(); ++i)
      if (mfs[i]->convex_index()[cv]) {
	pfem pfloc = mfs[i]->fem_of_element(cv);
	size_type nbdofloc = pfloc->nb_dof(cv);
	if (have_pfp) c.set_pfp(get_pfp(i, pfloc, &(pfp->get_point_tab())));
	else c.set_pf(pfloc);
	pfloc->real_base_value(c, taux[i]);
	  base_tensor::iterator it2 = taux[i].begin();
	  for (size_type k = 0; k < target_dim(); ++k)
	    for (size_type l = 0; l < nbdofloc; ++l)
	      *(vit[k]++) = *it2++;
      }
    if (have_pfp) c.set_pfp(pfp); else c.set_pf(pf);
  }
  
  void fem_sum::real_grad_base_value(const fem_interpolation_context& d, 
				     base_tensor &t) const {
    fem_interpolation_context& c = const_cast<fem_interpolation_context&>(d);
    size_type cv = c.convex_num(), nbdof = nb_dof(cv);
    mi3[2] = dim(); mi3[1] = target_dim(); mi3[0] = nbdof;
    t.adjust_sizes(mi3);
    bool have_pfp = c.have_pfp();
    pfem_precomp pfp = c.pfp();
    pfem pf = c.pf();
    base_tensor::iterator it1 = t.begin();
    for (size_type k = 0; k < size_type(target_dim() * dim()); ++k)
	{ vit[k] = it1; it1 += nbdof; }
    
    for (size_type i = 0; i < mfs.size(); ++i)
      if (mfs[i]->convex_index()[cv]) {
	pfem pfloc = mfs[i]->fem_of_element(cv);
	size_type nbdofloc = pfloc->nb_dof(cv);
	if (have_pfp) c.set_pfp(get_pfp(i, pfloc,&(pfp->get_point_tab())));
	else c.set_pf(pfloc);
	  pfloc->real_grad_base_value(c, taux[i]);
	  base_tensor::iterator it2 = taux[i].begin();
	  for (size_type k = 0; k < size_type(target_dim() * dim()); ++k)
	    for (size_type l = 0; l < nbdofloc; ++l)
	      *(vit[k]++) = *it2++;
      }
    if (have_pfp) c.set_pfp(pfp); else c.set_pf(pf);
  }
    
  void fem_sum::real_hess_base_value(const fem_interpolation_context&, 
				     base_tensor &) const
  { DAL_THROW(dal::to_be_done_error, "Sorry, to be done."); }
  
  fem_sum::fem_sum(const std::vector<const mesh_fem *> &mfs_)
  { init(mfs_); }

  DAL_SIMPLE_KEY(special_femsum_key, pfem);

  pfem new_fem_sum(const std::vector<const mesh_fem *> &mfs_) {
    pfem pf = new fem_sum(mfs_);
    dal::add_stored_object(new special_femsum_key(pf), pf);
    return pf;
  }
  
  pfem new_fem_sum(const mesh_fem &mef1, const mesh_fem &mef2) {
    std::vector<const mesh_fem *> mfs_(2);
    mfs_[0] = &mef1; mfs_[1] = &mef2;
    return new_fem_sum(mfs_);
  }
  
  pfem new_fem_sum(const mesh_fem &mef1, const mesh_fem &mef2,
		   const mesh_fem &mef3) {
    std::vector<const mesh_fem *> mfs_(3);
    mfs_[0] = &mef1; mfs_[1] = &mef2; mfs_[2] = &mef3;
    return new_fem_sum(mfs_);
  }
  
  }
}  /* end of namespace getfem.                                            */

