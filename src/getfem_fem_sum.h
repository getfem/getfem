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


    void update_from_context(void) const;
    void init(const std::vector<const mesh_fem *> &mfs_);
    pfem_precomp get_pfp(size_type i,pfem pf,
			 bgeot::pstored_point_tab pspti) const;
    
  public :

    virtual size_type nb_dof(size_type cv) const;
    virtual size_type index_of_global_dof(size_type cv,
					  size_type j) const;
    
    virtual bgeot::pconvex_ref ref_convex(size_type cv) const;
    virtual const bgeot::convex<base_node> &node_convex(size_type cv) const;
    virtual bgeot::pstored_point_tab node_tab(size_type) const;
    
    void base_value(const base_node &, base_tensor &) const;
    void grad_base_value(const base_node &, base_tensor &) const;
    void hess_base_value(const base_node &, base_tensor &) const;
    
    void real_base_value(const fem_interpolation_context& d, 
			 base_tensor &t) const;
    
    void real_grad_base_value(const fem_interpolation_context& d, 
			      base_tensor &t) const;
    
    void real_hess_base_value(const fem_interpolation_context&, 
			      base_tensor &) const;
    
    fem_sum(const std::vector<const mesh_fem *> &mfs_);
    fem_sum(const mesh_fem &mef1, const mesh_fem &mef2);
    fem_sum(const mesh_fem &mef1, const mesh_fem &mef2, const mesh_fem &mef3);
    
  };


}  /* end of namespace getfem.                                            */

#endif
