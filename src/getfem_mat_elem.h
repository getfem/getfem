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
/* Copyright (C) 2000-2002  Yves Renard.                                   */
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


#ifndef __GETFEM_MAT_ELEM_H
#define __GETFEM_MAT_ELEM_H

#include <getfem_mat_elem_type.h>
#include <bgeot_geometric_trans.h>
#include <getfem_fem.h>

namespace getfem
{

  template <class CONT> void transfert_to_G(base_matrix &G, const CONT &a) {
    size_type P = (*(a.begin())).size(), NP = a.end() - a.begin();
    G.resize(P, NP);
    typename CONT::const_iterator it = a.begin(), ite = a.end();
    base_matrix::iterator itm = G.begin();
    for (; it != ite; ++it, itm += P)
      std::copy((*it).begin(), (*it).end(), itm);
  }

  class mat_elem_computation
  {
    protected : 

      bgeot::pgeometric_trans pgt;
      pmat_elem_type pme;
      base_matrix pa;

    public :

      virtual void compute(base_tensor &t, const base_matrix &a,
			   size_type elt) = 0;
      virtual void compute_on_face(base_tensor &t, const base_matrix &a,
				   short_type f, size_type elt) = 0;
      template <class CONT>
	void gen_compute(base_tensor &t, const CONT &a,  size_type elt)
      { transfert_to_G(pa, a); compute(t, pa,elt); }
      template <class CONT>
	void gen_compute_on_face(base_tensor &t,
				 const CONT &a, short_type f, size_type elt)
      { transfert_to_G(pa, a); compute_on_face(t, pa, f, elt); }

      virtual ~mat_elem_computation() {}
      virtual size_type memsize() const = 0;
  };

  typedef mat_elem_computation *pmat_elem_computation;

  pmat_elem_computation mat_elem(pmat_elem_type pm, 
				 pintegration_method pi,
				 bgeot::pgeometric_trans pg);


}  /* end of namespace getfem.                                             */


#endif /* __GETFEM_MAT_ELEM_H                                              */
