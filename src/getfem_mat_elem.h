/* *********************************************************************** */
/*                                                                         */
/* Library :  GEneric Tool for Finite Element Methods (getfem)             */
/* File    :  getfem_mat_elem.h : computation of elementary matrices.      */
/*     									   */
/*                                                                         */
/* Date : December 21, 2000.                                               */
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


#ifndef __GETFEM_MAT_ELEM_H
#define __GETFEM_MAT_ELEM_H

#include <getfem_mat_elem_type.h>
#include <bgeot_geometric_trans.h>
#include <getfem_fem.h>

namespace getfem
{

  class mat_elem_computation
  {
    protected : 

      bgeot::pgeometric_trans pgt;
      pmat_elem_type pme;
      base_matrix pa;

    public :

      virtual void compute(base_tensor &t, const base_matrix &a) = 0;
      virtual void compute_on_face(base_tensor &t, const base_matrix &a,
				   short_type f) = 0;
      template <class CONT> void set_pa(const CONT &a)
      {
	size_type P = (*(a.begin())).size(), NP = a.end() - a.begin();
	pa.resize(P, NP);
	typename CONT::const_iterator it = a.begin(), ite = a.end();
	base_matrix::iterator itm = pa.begin();
	for (; it != ite; ++it, itm += P)
	  std::copy((*it).begin(), (*it).end(), itm);
      }

      template <class CONT> void gen_compute(base_tensor &t, const CONT &a)
      { set_pa(a); compute(t, pa); }
      template <class CONT> void gen_compute_on_face(base_tensor &t,
						 const CONT &a, short_type f)
      { set_pa(a); compute_on_face(t, pa, f); }
  };

  typedef mat_elem_computation *pmat_elem_computation;

  struct pintegration_method
  {
    union
    {
      bgeot::ppoly_integration ppi;
      bgeot::papprox_integration pai;
    } method;
    bool is_ppi;

    const bgeot::stored_point_tab &integration_points(void) const { 
      if (is_ppi)
	return *(bgeot::org_stored_point_tab(method.pai->structure()->dim()));
      else 
	return method.pai->integration_points();
    }

    pintegration_method(bgeot::ppoly_integration p)
    { method.ppi = p; is_ppi = true; }

    pintegration_method(bgeot::papprox_integration p)
    { method.pai = p; is_ppi = false; }

    bool operator >(const pintegration_method& p) const
    { return method.ppi > p.method.ppi; } 
    bool operator <(const pintegration_method& p) const
    { return method.ppi < p.method.ppi; } 
    bool operator !=(const pintegration_method& p) const
    { return method.ppi != p.method.ppi; } 
    bool operator ==(const pintegration_method& p) const
    { return method.ppi == p.method.ppi; } 

    pintegration_method(void) {}

  };

  pmat_elem_computation mat_elem(pmat_elem_type pm, pintegration_method pi,
				 bgeot::pgeometric_trans pg);

  base_vector compute_normal(const base_matrix &G, size_type ir, bgeot::pgeometric_trans pgt, const base_node &pt);

}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_MAT_ELEM_H                                              */
