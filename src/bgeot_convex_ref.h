// -*- c++ -*- (enables emacs c++ mode)
//========================================================================
//
// Library : Basic GEOmetric Tool  (bgeot)
// File    : bgeot_convex_ref.h : convexes of reference
//           
// Date    : Septembre 28, 2001.
// Author  : Yves Renard <Yves.Renard@insa-toulouse.fr>
//
//========================================================================
//
// Copyright (C) 2001-2005 Yves Renard
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


#ifndef BGEOT_CONVEX_REF_H__
#define BGEOT_CONVEX_REF_H__

#include <bgeot_convex.h>

namespace bgeot {

  /** Point tab storage. */
  struct stored_point_tab : public dal::static_stored_object,
			    public std::vector<base_node> {
    const base_node &operator[](size_type i) const
    { return std::vector<base_node>::operator [](i); }
    template <class IT> stored_point_tab(IT it, IT ite)
      : std::vector<base_node>(it, ite) { }
  };

  typedef boost::intrusive_ptr<const stored_point_tab> pstored_point_tab;

  pstored_point_tab store_point_tab(const stored_point_tab& spt);

  /* store a new (read-only) array of points in stored_point_tab_tab */
  template<class CONT> pstored_point_tab store_point_tab(const CONT &TAB)
  { return store_point_tab(stored_point_tab(TAB.begin(), TAB.end())); }

  class mesh_structure;

  /**
     convex_of_reference: base class for reference convexes (the order
     1 triangle (0,0)-(1,0)-(0,1), the order 2 segment
     (0)-(.5)-(1.), etc...)  stores :

       - a list of points (vertices of the reference convex, plus
       other points for reference convexes of degree > 1)

       - a normal for each face of the convex

       - a mesh structure defining the smallest simplex partition of
       the convex

       - a pointer to the "basic convex_ref": for a convex_ref of
       degree k, this is a pointer to the correspounding convex_ref of
       degree 1.
   */

  class convex_of_reference;
  typedef boost::intrusive_ptr<const convex_of_reference> pconvex_ref;
  
  
  class convex_of_reference : public dal::static_stored_object, 
			      public convex<base_node> {
  protected :     
    std::vector<base_small_vector> normals_;
    pstored_point_tab ppoints;
    mutable mesh_structure *psimplexified_convex;
    mutable const convex_of_reference *basic_convex_ref_;
  public :
    convex_of_reference() : convex<base_node>(), psimplexified_convex(0),
			    basic_convex_ref_(0) {}
    virtual scalar_type is_in(const base_node &) const = 0;
    virtual scalar_type is_in_face(short_type, const base_node &) const =0;
    const std::vector<base_small_vector> &normals(void) const
    { return normals_; }

    const stored_point_tab &points(void) const { return *ppoints; }
    const stored_point_tab &points(void) { return *ppoints; }
    
    /* returns a mesh structure composed of simplexes whose union
       is the reference convex. All simplexes have the same (direct)
       orientation.
    */
    const mesh_structure* simplexified_convex() const;
    pconvex_ref basic_convex_ref() const { return basic_convex_ref_; }
    void attach_basic_convex_ref(pconvex_ref cvr) const
    { basic_convex_ref_ = cvr.get(); }
  };

  
  /* these are the public functions for obtaining a convex of reference */

  /** returns a simplex of reference of dimension nc and degree k */
  pconvex_ref simplex_of_reference(dim_type nc, short_type k = 1);
  /** parallelepiped of reference of dimension nc (and degree 1) */
  pconvex_ref parallelepiped_of_reference(dim_type nc);
  /** tensorial product of two convex ref.
      in order to ensure unicity, it is required the a->dim() >= b->dim() */
  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b);
  /** equilateral simplex (degree 1). used only for mesh quality estimations */
  pconvex_ref equilateral_simplex_of_reference(dim_type nc);

  /** generic convex with n global nodes      */
  pconvex_ref generic_dummy_convex_ref(dim_type nc, size_type n, size_type nf);

}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_CONVEX_REF_H__                                            */
