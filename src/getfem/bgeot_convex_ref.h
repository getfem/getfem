/* -*- c++ -*- (enables emacs c++ mode) */
/*===========================================================================
 
 Copyright (C) 2001-2012 Yves Renard
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 2.1 of the License,  or
 (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 
 As a special exception, you  may use  this file  as it is a part of a free
 software  library  without  restriction.  Specifically,  if   other  files
 instantiate  templates  or  use macros or inline functions from this file,
 or  you compile this  file  and  link  it  with other files  to produce an
 executable, this file  does  not  by itself cause the resulting executable
 to be covered  by the GNU Lesser General Public License.  This   exception
 does not  however  invalidate  any  other  reasons why the executable file
 might be covered by the GNU Lesser General Public License.
 
===========================================================================*/


/**@file bgeot_convex_ref.h
   @author  Yves Renard <Yves.Renard@insa-lyon.fr>
   @date Septembre 28, 2001.
   @brief Reference convexes
*/
#ifndef BGEOT_CONVEX_REF_H__
#define BGEOT_CONVEX_REF_H__

#include "bgeot_convex.h"

namespace bgeot {

  /** @defgroup refconvexes Reference Convexes*/
  /*@{*/

  /** Point tab storage. */
  struct stored_point_tab : virtual public dal::static_stored_object,
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

  class convex_of_reference;
  typedef boost::intrusive_ptr<const convex_of_reference> pconvex_ref;
    
  /** Base class for reference convexes.

     Examples of reference convexes are (the order
     1 triangle (0,0)-(1,0)-(0,1), the order 2 segment
     (0)-(.5)-(1.), etc...). This class stores :

       - a list of points (vertices of the reference convex, plus
       other points for reference convexes of degree > 1).

       - a normal for each face of the convex.

       - a mesh structure defining the smallest simplex partition of
       the convex.

       - a pointer to the "basic convex_ref": for a convex_ref of
       degree k, this is a pointer to the correspounding convex_ref of
       degree 1.
   */
  class convex_of_reference : virtual public dal::static_stored_object, 
			      public convex<base_node> {
  protected :     
    std::vector<base_small_vector> normals_;
    pstored_point_tab ppoints;
    mutable mesh_structure *psimplexified_convex;
    mutable const convex_of_reference *basic_convex_ref_;
    convex_of_reference() : convex<base_node>(), psimplexified_convex(0),
			    basic_convex_ref_(0) {}

  public :
    /// return a negative or null number if the base_node is in the convex.
    virtual scalar_type is_in(const base_node &) const = 0;
    /** return a null (or almost zero) if pt is in the face of the convex.
     *  Does not control if the point is in the convex, but if a point
     *  supposed to be in a convex is in this face.
     */
    virtual scalar_type is_in_face(short_type, const base_node &) const =0;
    /// return the normal vector for each face.
    const std::vector<base_small_vector> &normals(void) const
    { return normals_; }
    /// return the vertices of the reference convex.
    const stored_point_tab &points(void) const { return *ppoints; }
    const stored_point_tab &points(void) { return *ppoints; }
    
    /** return a mesh structure composed of simplexes whose union
	is the reference convex. All simplexes have the same (direct)
	orientation.
    */
    const mesh_structure* simplexified_convex() const;
    /// return the associated order 1 reference convex.
    pconvex_ref basic_convex_ref() const { return basic_convex_ref_; }
    /// private function..
    void attach_basic_convex_ref(pconvex_ref cvr) const
    { basic_convex_ref_ = cvr.get(); }
  };

  
  //@name public functions for obtaining a convex of reference
  //@{

  /** returns a simplex of reference of dimension nc and degree k */
  pconvex_ref simplex_of_reference(dim_type nc, short_type k = 1);
  /** parallelepiped of reference of dimension nc (and degree 1) */
  pconvex_ref parallelepiped_of_reference(dim_type nc);
  /** prism of reference of dimension nc (and degree 1) */
  pconvex_ref prism_of_reference(dim_type nc);
  /** tensorial product of two convex ref.
      in order to ensure unicity, it is required the a->dim() >= b->dim() */
  pconvex_ref convex_ref_product(pconvex_ref a, pconvex_ref b);
  /** equilateral simplex (degree 1). used only for mesh quality estimations */
  pconvex_ref equilateral_simplex_of_reference(dim_type nc);

  /** generic convex with n global nodes      */
  pconvex_ref generic_dummy_convex_ref(dim_type nc, size_type n, size_type nf);
  //@}

  /*@}*/
}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_CONVEX_REF_H__                                            */
