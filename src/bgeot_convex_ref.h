/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Basic GEOmetric Tool  (bgeot)                                */
/* File    :  bgeot_convex_ref.h :  convexes of reference                  */
/*     									   */
/* Date : Septembre 28, 2001.                                              */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2001-2002  Yves Renard.                                   */
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



#ifndef BGEOT_CONVEX_REF_H__
#define BGEOT_CONVEX_REF_H__

#include <bgeot_convex.h>

namespace bgeot
{
  

  /** 
      Point tab storage.
  */
  typedef std::vector<base_node> stored_point_tab;
  typedef const stored_point_tab * pstored_point_tab;

  class comp_stored_point_tab
    : public std::binary_function<stored_point_tab, stored_point_tab, int>
  {
    public :
    int operator()(const stored_point_tab &x,
		   const stored_point_tab &y) const;
  };

  extern dal::dynamic_tree_sorted<stored_point_tab, comp_stored_point_tab>
    *stored_point_tab_tab_;
  extern bool isinit_stored_point_tab_tab;

  template<class CONT> pstored_point_tab store_point_tab(const CONT &TAB)
  { 
    if (!isinit_stored_point_tab_tab) {
      stored_point_tab_tab_ =
	new dal::dynamic_tree_sorted<stored_point_tab,
	                             comp_stored_point_tab>();
      isinit_stored_point_tab_tab = true;
    }
    typename CONT::const_iterator it = TAB.begin(), ite = TAB.end();
    size_type nb;
    for (nb = 0; it != ite; ++it, ++nb);
    stored_point_tab spt; spt.resize(nb);
    it = TAB.begin(); ite = TAB.end();
    for (nb = 0; it != ite; ++it, ++nb) spt[nb] = *it;
    return &((*stored_point_tab_tab_)[stored_point_tab_tab_->add_norepeat(spt)]);
  }

  pstored_point_tab org_stored_point_tab(size_type n);

  class mesh_structure;

  /* structures de reference.                                             */

  class convex_of_reference : public convex<base_node> {
  protected :     
    std::vector<base_small_vector> normals_;
    mutable mesh_structure *psimplexified_convex;
    convex_of_reference *basic_convex_ref_;
  public :
    convex_of_reference() : convex<base_node>(), psimplexified_convex(0), basic_convex_ref_(0) {}
    virtual scalar_type is_in(const base_node &) const = 0;
    virtual scalar_type is_in_face(short_type, const base_node &) const =0;
    const std::vector<base_small_vector> &normals(void) const { return normals_; }
    virtual ~convex_of_reference() {}

    /* returns a mesh structure composed of simplexes whose union
       is the reference convex
    */
    const mesh_structure* simplexified_convex() const;
    const convex_of_reference* basic_convex_ref() const { return basic_convex_ref_; }
    void attach_basic_convex_ref(convex_of_reference* cvr) { basic_convex_ref_ = cvr; }
  };

  typedef const convex_of_reference * pconvex_ref;

  pconvex_ref simplex_of_reference(dim_type nc, short_type k = 1);
  pconvex_ref parallelepiped_of_reference(dim_type nc);
  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b);
  pconvex_ref equilateral_simplex_of_reference(dim_type nc);

  /* fonctions en sursis ... */
  // pconvex_ref multiply_convex_of_reference(pconvex_ref a, dim_type n);

  // pconvex_ref nonconforming_triangle_ref(void);

}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_CONVEX_REF_H__                                            */
