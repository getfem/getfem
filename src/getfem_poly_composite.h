/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library : GEneric Tool for Finite Element Methods (getfem)              */
/* File    : getfem_poly_composite.h : polynomials by parts                */
/*                                                                         */
/* Date : August 26, 2002.                                                 */
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


#ifndef __GETFEM_POLY_COMPOSIT_H
#define __GETFEM_POLY_COMPOSIT_H

#include <bgeot_poly.h>
#include <bgeot_geotrans_inv.h>
#include <getfem_mesh.h>

namespace getfem
{

  struct mesh_precomposite {

    typedef dal::dynamic_tree_sorted<base_node,
      bgeot::imbricated_box_less> PTAB;

    const getfem_mesh *mesh;
    PTAB vertexes;
    std::vector<base_matrix> gtrans;
    std::vector<base_node> orgs;
    std::vector<bool> elt;
    
    const getfem_mesh &linked_mesh(void) const { return *mesh; }
    size_type nb_convex(void) const { return gtrans.size(); }
    dim_type dim(void) const { return mesh->dim(); }
    bgeot::pgeometric_trans trans_of_convex(size_type ic) const
    { return mesh->trans_of_convex(ic); }
    
    mesh_precomposite(const getfem_mesh &m);
  };

  class polynomial_composite {

  protected :
    const mesh_precomposite *mp;
    std::vector<bgeot::base_poly> polytab;

  public :
    
    template <class ITER> scalar_type eval(const ITER &it) const;
    void derivative(short_type k);

    polynomial_composite(void) {}
    polynomial_composite(const mesh_precomposite &m)
      : mp(&m), polytab(m.nb_convex()) {}

  };

  template <class ITER>
  scalar_type polynomial_composite::eval(const ITER &it) const {
    dim_type N = mp->dim();
    base_node pt(N), p0;
    std::copy(it, it + N, pt.begin());
    std::fill(mp->elt.begin(), mp->elt.end(), true);
    size_type i = size_type(-1);
    
    mesh_precomposite::PTAB::const_sorted_iterator
      it1 = mp->vertexes.sorted_ge(pt), it2 = it1;    
    size_type i1 = it1.index(), i2 = it2.index();
    --it2;

    while (i1 != size_type(-1) || i2 != size_type(-1)) {
      if (i1 != size_type(-1)) {
	size_type ii1 = md->linked_mesh().convex_to_point(i1)[0];
	if (mp->elt[ii1]) {
	  mp->elt[ii1] = false;
	  p0 = pt; p0 -= md->orgs[ii1];
	  p0 *= md->gtrans[ii1];
	  if (md->trans_of_convex(ii1).convex_ref()->is_in(p0) < 1E-10)
	    { i = ii1; break; }
	}
	++it1; i1 = it1.index();
      }
      if (i2 != size_type(-1)) {
	size_type ii2 = md->linked_mesh().convex_to_point(i2)[0];
	if (mp->elt[ii2]) {
	  mp->elt[ii2] = false;
	  p0 = pt; p0 -= md->orgs[ii2];
	  p0 *= md->gtrans[ii2];
	  if (md->trans_of_convex(ii2).convex_ref()->is_in(p0) < 1E-10)
	    { i = ii2; break; }
	}
	--it2; i2 = it2.index();
      }
    }
    if (i == size_type(-1))
      DAL_THROW(internal_error, "Element not found in composite polynomial.");
    polytab[i].eval(it);
  }


  
}  /* end of namespace getfem.                                            */


#endif
