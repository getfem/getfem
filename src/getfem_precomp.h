/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Elements Methods (getfem)            */
/* File    :  getfem_precomp.h : pre-computations.                         */
/*     									   */
/*                                                                         */
/* Date : June 17, 2002.                                                   */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GETFEM++                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; version 2 of the License.                 */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */
/*                                                                         */
/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software Foundation, */
/* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.         */
/*                                                                         */
/* *********************************************************************** */

#ifndef __GETFEM_PRECOMP_H
#define __GETFEM_PRECOMP_H

#include <bgeot_geometric_trans.h>
#include <getfem_fem.h>

namespace getfem
{

  /* ********************************************************************* */
  /*       Precomputation on geometric transformations.                    */
  /* ********************************************************************* */

  struct _pre_geot_light;

  class _geotrans_precomp
  {
    protected :

      bgeot::pgeometric_trans pgt;
      bgeot::pstored_point_tab pspt;
      std::vector<base_matrix> pc;
      std::vector<base_matrix> hpc;

    public :

      inline const base_matrix &grad(size_type i) const { return pc[i]; }
      inline const base_matrix &hessian(size_type i) const { return hpc[i]; }

      _geotrans_precomp(const _pre_geot_light &ls);
  };
  
  typedef const _geotrans_precomp * pgeotrans_precomp;

  pgeotrans_precomp geotrans_precomp(bgeot::pgeometric_trans pg,
				     bgeot::pstored_point_tab pspt);

  /* ********************************************************************* */
  /*       Precomputation on fem.                                          */
  /* ********************************************************************* */

  class virtual_fem;
  typedef const virtual_fem * pfem;

  struct _pre_fem_light;
  
  class _fem_precomp
  {
    protected :
      
      std::vector<base_tensor> c;
      std::vector<base_tensor> pc;
      std::vector<base_tensor> hpc;

    public :

      inline const base_tensor &val(size_type i) const { return c[i]; }
      inline const base_tensor &grad(size_type i) const { return pc[i]; }
      inline const base_tensor &hess(size_type i) const { return hpc[i]; }

      _fem_precomp(const _pre_fem_light &ls);
  };
  
  typedef const _fem_precomp * pfem_precomp;

  pfem_precomp fem_precomp(pfem pf, bgeot::pstored_point_tab pspt);

}  /* end of namespace getfem.                                            */

#endif
