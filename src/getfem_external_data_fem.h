/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_fem_gauss_points.h : definition a "fem" allowing to    */
/*           define a function only on gauss points.                       */
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

  class external_data_fem : public virtual_fem {
    
  public :
    
    size_type nb_dof(void) const { return 1; }
    
    void interpolation(const base_node &x, const base_matrix &G,
		       bgeot::pgeometric_trans pgt,
		       const base_vector &coeff, base_node &val) const
    { DAL_THROW(failure_error, "You cannot interpolate this element"); }
    
    void interpolation(pfem_precomp pfp, size_type ii,
		       const base_matrix &G,
		       bgeot::pgeometric_trans pgt, 
		       const base_vector &coeff, 
		       base_node &val, dim_type Qdim=1) const
    { DAL_THROW(failure_error, "You cannot interpolate this element"); }

    void interpolation_grad(const base_node &x, const base_matrix &G,
			    bgeot::pgeometric_trans pgt,
			    const base_vector &coeff, base_matrix &val) const
    { DAL_THROW(failure_error, "You cannot interpolate this element"); }

    void base_value(const base_node &x, base_tensor &t) const
    { DAL_THROW(failure_error, "You cannot interpolate this element"); }
    void grad_base_value(const base_node &x, base_tensor &t) const
    { DAL_THROW(failure_error, "You cannot interpolate this element"); }
    void hess_base_value(const base_node &x, base_tensor &t) const
    { DAL_THROW(failure_error, "You cannot interpolate this element"); }

    void real_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			 size_type ii, const base_matrix &G,
			 base_tensor &t) const;
    
    void real_grad_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			      size_type ii, const base_matrix &G,
			      const base_matrix &B, base_tensor &t) const;
    
    void real_hess_base_value(pgeotrans_precomp pgp, pfem_precomp pfp,
			      size_type ii, const base_matrix &G,
			      const base_matrix &B3, const base_matrix &B32,
			      base_tensor &t) const;

    gauss_points_fem(const bgeot::pconvex_ref _cvr, size_type dim = 1) {
      cvr = _cvr;
      is_equiv = real_element_defined = true;
      is_polycomp = is_pol = is_lag = false;
      es_degree = 5;
      ntarget_dim = dim;
    }

  };


}  /* end of namespace getfem.                                            */
