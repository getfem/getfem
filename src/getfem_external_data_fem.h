/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_external_data_fem.h : definition a "fem" allowing to   */
/*           define any function.                                          */
/*                                                                         */
/* Date : June 28, 2003.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2003  Yves Renard.                                   */
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

#include <getfem_fem.h>
#include <getfem_mesh_fem.h>
#include <bgeot_geotrans_inv.h>

namespace getfem
{

  // Element to be derived. Overload at least real_base_value.
  // always 1 ddl

  class external_data_fem : public virtual_fem {
    
  public :
    
    void interpolation(const base_node &, const base_matrix &,
		       bgeot::pgeometric_trans ,
		       const base_vector &, base_node &) const
    { DAL_THROW(failure_error, "Uninstantied function"); }
    
    void interpolation(pfem_precomp , size_type ,
		       const base_matrix &,
		       bgeot::pgeometric_trans , 
		       const base_vector &, 
		       base_node &, dim_type=1) const
    { DAL_THROW(failure_error, "Uninstantied function"); }

    void interpolation_grad(const base_node &, const base_matrix &,
			    bgeot::pgeometric_trans ,
			    const base_vector &, base_matrix &) const
    { DAL_THROW(failure_error, "Uninstantied function"); }

    void base_value(const base_node &, base_tensor &) const
    { DAL_THROW(failure_error, "Uninstantied function"); }
    void grad_base_value(const base_node &, base_tensor &) const
    { DAL_THROW(failure_error, "Uninstantied function"); }
    void hess_base_value(const base_node &, base_tensor &) const
    { DAL_THROW(failure_error, "Uninstantied function"); }

    void real_base_value(pgeotrans_precomp , pfem_precomp ,
			 size_type , const base_matrix &,
			 base_tensor &, size_type ) const
    { DAL_THROW(failure_error, "Uninstantied function"); } 
    
    void real_grad_base_value(pgeotrans_precomp , pfem_precomp ,
			      size_type , const base_matrix &,
			      const base_matrix &, base_tensor &,
			      size_type ) const
    { DAL_THROW(failure_error, "Uninstantied function"); }
    
    void real_hess_base_value(pgeotrans_precomp , pfem_precomp ,
			      size_type , const base_matrix &,
			      const base_matrix &, const base_matrix &,
			      base_tensor &,
			      size_type ) const
    { DAL_THROW(failure_error, "Uninstantied function"); }

    external_data_fem(const bgeot::pconvex_ref _cvr, size_type dim = 1) {
      cvr = _cvr;
      is_equiv = real_element_defined = true;
      is_polycomp = is_pol = is_lag = false;
      es_degree = 5;
      ntarget_dim = dim;
      init_cvs_node();
      base_node pt(cvr->structure()->dim()); pt.fill(0.001);
      add_node(lagrange_dof(cvr->structure()->dim()), pt);
    }
  };


}  /* end of namespace getfem.                                            */
