/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_fem_sum.h : make the direct sum of a set of finite     */
/*           element method.                                               */
/*                                                                         */
/* Date : October 29, 2004.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


#ifndef GETFEM_FEM_SUM_H__
#define GETFEM_FEM_SUM_H__

#include <getfem_fem.h>
#include <getfem_mesh_fem.h>

namespace getfem {

  class fem_sum : public virtual_fem, public context_dependencies {
    
  protected :

    std::vector<const mesh_fem *> mfs;

    mutable std::vector<base_node> node_tab_;
    mutable std::vector<base_tensor> taux;
    mutable bgeot::multi_index mi2, mi3;
    mutable std::vector<base_tensor::iterator> vit;
    mutable std::vector<pfem> spf;
    mutable std::vector<bgeot::pstored_point_tab> spspt;
    mutable std::vector<pfem_precomp> spfp;


    void build_fem(void) const {
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
		already_numerate_dof(dim()));

      vit.resize(target_dim()*dim());
      
    }

    void init(const std::vector<const mesh_fem *> &mfs_) {
      mfs = mfs_;
      const getfem_mesh *pmesh = &(mfs[0]->linked_mesh());
      for (size_type i = 0; i < mfs.size(); ++i) {
	if (&(mfs[i]->linked_mesh()) != pmesh)
	  DAL_THROW(failure_error, "Meshes should be the same");
	this->add_dependency(*(mfs[i]));
      }
      is_pol = is_lag = false; es_degree = 5;
      is_equiv = real_element_defined = true;
      ntarget_dim = 1; // An extension for vectorial elements should be easy
      // The detection should be done and the multiplication of components
      // for scalar elements interpolated.
      mi2 = bgeot::multi_index(2);
      mi3 = bgeot::multi_index(3);
      build_fem();
      spf.resize(mfs.size());
      spspt.resize(mfs.size());
      spfp.resize(mfs.size());
      taux.resize(mfs.size());
    }

  public :

    virtual size_type nb_dof(size_type cv) const {
      if (context_changed()) build_fem(); 
      size_type res = 0;
      for (size_type i = 0; i < mfs.size(); ++i)
	if (mfs[i]->convex_index()[cv])
	  res += mfs[i]->nb_dof_of_element(cv);
      return res;
    }
    virtual size_type index_of_already_numerate_dof(size_type cv,
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
    
    virtual bgeot::pconvex_ref ref_convex(size_type cv) const {
      for (size_type i = 0; i < mfs.size(); ++i)
	if (mfs[i]->convex_index()[cv])
	  return mfs[i]->fem_of_element(cv)->ref_convex(cv);
    }

    virtual const bgeot::convex<base_node> &node_convex(size_type cv) const
    { return *(bgeot::generic_dummy_convex_ref(dim(), nb_dof(cv))); }
    
    virtual bgeot::pstored_point_tab node_tab(size_type) const { 
      if (!pspt_valid)
	{ pspt = bgeot::store_point_tab(node_tab_); pspt_valid = true; }
      return pspt;
    }
    
    void base_value(const base_node &, base_tensor &) const
    { DAL_THROW(internal_error, "No base values, real only element."); }
    
    void grad_base_value(const base_node &, base_tensor &) const
    { DAL_THROW(internal_error, "No grad values, real only element."); }
    
    void hess_base_value(const base_node &, base_tensor &) const
    { DAL_THROW(internal_error, "No hess values, real only element."); }
    
    pfem_precomp get_pfp(size_type i,pfem pf,
			 bgeot::pstored_point_tab pspti) const {
      if (pf != spf[i] || pspti != spspt[i])
	{ spfp[i]=fem_precomp(pf, pspti); spf[i] = pf; spspt[i] = pspti; }
      return spfp[i];
    }

    void real_base_value(const fem_interpolation_context& d, 
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
    
    void real_grad_base_value(const fem_interpolation_context& d, 
			      base_tensor &t) const {
      fem_interpolation_context& c = const_cast<fem_interpolation_context&>(d);
      size_type cv = c.convex_num(), nbdof = nb_dof(cv);
      mi3[2] = dim(); mi3[1] = target_dim(); mi2[0] = nbdof;
      t.adjust_sizes(mi3);
      bool have_pfp = c.have_pfp();
      pfem_precomp pfp = c.pfp();
      pfem pf = c.pf();
      base_tensor::iterator it1 = t.begin();
      for (size_type k = 0; k < target_dim() * dim(); ++k)
	{ vit[k] = it1; it1 += nbdof; }

      for (size_type i = 0; i < mfs.size(); ++i)
	if (mfs[i]->convex_index()[cv]) {
	  pfem pfloc = mfs[i]->fem_of_element(cv);
	  size_type nbdofloc = pfloc->nb_dof(cv);
	  if (have_pfp) c.set_pfp(get_pfp(i, pfloc,&(pfp->get_point_tab())));
	  else c.set_pf(pfloc);
	  pfloc->real_grad_base_value(c, taux[i]);
	  base_tensor::iterator it2 = taux[i].begin();
	  for (size_type k = 0; k < target_dim() * dim(); ++k)
	    for (size_type l = 0; l < nbdofloc; ++l)
	      *(vit[k]++) = *it2++;
	}
      if (have_pfp) c.set_pfp(pfp); else c.set_pf(pf);
    }
    
    void real_hess_base_value(const fem_interpolation_context&, 
			      base_tensor &) const
    { DAL_THROW(internal_error, "Sorry, to be done."); }
    
    fem_sum(const std::vector<const mesh_fem *> &mfs_)
    { init(mfs_); }

    fem_sum(const mesh_fem &mef1, const mesh_fem &mef2) {
      std::vector<const mesh_fem *> mfs_(2);
      mfs[1] = &mef1; mfs[2] = &mef2;
      init (mfs_);
    }

    fem_sum(const mesh_fem &mef1, const mesh_fem &mef2, const mesh_fem &mef3) {
      std::vector<const mesh_fem *> mfs_(3);
      mfs[1] = &mef1; mfs[2] = &mef2; mfs[3] = &mef3;
      init (mfs_);
    }
    
  };


}  /* end of namespace getfem.                                            */

#endif
