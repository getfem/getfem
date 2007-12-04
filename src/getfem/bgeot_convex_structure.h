// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 1999-2008 Yves Renard
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
// As a special exception, you may use this file as part of a free software
// library without restriction.  Specifically, if other files instantiate
// templates or use macros or inline functions from this file, or you compile
// this file and link it with other files to produce an executable, this
// file does not by itself cause the resulting executable to be covered by
// the GNU General Public License.  This exception does not however
// invalidate any other reasons why the executable file might be covered by
// the GNU General Public License.
//
//===========================================================================

#ifndef BGEOT_CONVEX_STRUCTURE_H__
#define BGEOT_CONVEX_STRUCTURE_H__

/** @file bgeot_convex_structure.h
    @author  Yves Renard <Yves.Renard@insa-lyon.fr>
    @date December 20, 1999.
    @brief Definition of convex structures
 */

#include "gmm/gmm_ref.h"
#include "bgeot_config.h"
#include "bgeot_tensor.h"
#include "bgeot_poly.h"
#include "dal_static_stored_objects.h"

namespace bgeot {
  /** @defgroup convex_structure Convex Structure */
  /*@{*/

  // The number of faces for a convex is limited in certain applications
#  define MAX_FACES_PER_CV 31

  class convex_structure;

  /// Pointer on a convex structure description. 
  typedef boost::intrusive_ptr<const convex_structure> pconvex_structure;

  typedef std::vector<const convex_structure *> convex_structure_faces_ct;
  typedef std::vector<short_type>               convex_ind_ct;
  typedef gmm::tab_ref_index_ref< convex_ind_ct::const_iterator,
	         convex_ind_ct::const_iterator> ref_convex_ind_ct;

  /**  Structure of a convex.  
   *
   * This class is not to be manipulate by itself. Use
   * pconvex_structure and the functions written to produce the
   * convex structures from classicals convexes (simplexes, polygonals
   * ...). The reason is that there is no need for having more than
   * one convex structure for the same type of convex.
   */
  class convex_structure : virtual public dal::static_stored_object {
  protected :
    
    dim_type Nc;
    short_type nbpt, nbf;
    convex_structure_faces_ct  faces_struct;
    std::vector<convex_ind_ct> faces;
    convex_ind_ct              dir_points_;
    const convex_structure *basic_pcvs;
    
    pconvex_structure prod_a, prod_b; /* only filled for convex structures */
				      /* product.                          */
    public :

      /// Number of faces.
      inline short_type nb_faces(void)  const { return nbf;  }
      /// Dimension of the convex.
      inline dim_type  dim(void)        const { return Nc;   }
      /// Number of vertices.
      inline short_type nb_points(void) const { return nbpt; }
      /// Original structure (if concerned).
      pconvex_structure basic_structure(void) const 
      { return basic_pcvs; }
    /** Number of vertices of a face.
     *	@param i the face number.
     */
    inline short_type nb_points_of_face(short_type i) const
    { return faces[i].size(); }
    /** Give an array of the indexes of the vertices of a face. 
     *	The indexes are "local" to the convex.
     *	@param i the face number.
     */
    inline const convex_ind_ct &ind_points_of_face(short_type i) const
    { return faces[i]; }
    /** Return "direct" points indexes. These are the subset of points than
     *	can be used to build a direct vector basis. (rarely used)
     */
    inline const convex_ind_ct &ind_dir_points() const
    { return dir_points_; }
    /** Give a pointer array on the structures of the faces.
     *   faces_structure()[i] is a pointer on the structure of the face i.
     */
    inline const convex_structure_faces_ct &faces_structure(void) const
    { return faces_struct; }
    /** Return "direct" points indexes for a given face.
     *	@param i the face number.
     */
    inline ref_convex_ind_ct ind_dir_points_of_face(short_type i) const {
      return ref_convex_ind_ct(faces[i].begin(),
			       faces_struct[i]->ind_dir_points().begin(),
			       faces_struct[i]->ind_dir_points().end());
    }
    
    void init_for_adaptative(pconvex_structure cvs);
    void add_point_adaptative(short_type i, short_type f);
    /** Return true if the convex structure is indeed a direct product
     *	of two convex structures. 
     *	@param pprod1 the first sub-structure (optional)
     *	@param pprod2 the second sub-structure (optional)
     */
    bool is_product(pconvex_structure *pprod1=0,
		    pconvex_structure *pprod2=0) const {
      if (pprod1) *pprod1 = prod_a;
      if (pprod2) *pprod2 = prod_b;
      return prod_a ? true : false;
    }
  protected:
    convex_structure() { prod_a = prod_b = 0; }
    friend boost::intrusive_ptr<convex_structure> new_convex_structure();
  };

  inline boost::intrusive_ptr<convex_structure> new_convex_structure()
  { return boost::intrusive_ptr<convex_structure>(new convex_structure); }

  /** @name functions on convex structures
   */
  //@{

  /** Print the details of the convex structure cvs to the output stream o.
   *   For debuging purpose.
   */
  std::ostream &operator << (std::ostream &o,
				   const convex_structure &cv);

  /// Give a pointer on the structures of a simplex of dimension d.
  pconvex_structure simplex_structure(dim_type d);
  /// Give a pointer on the structures of a parallelepiped of dimension d.
  pconvex_structure parallelepiped_structure(dim_type d);
  /// Give a pointer on the structures of a polygon with n vertex.
  pconvex_structure polygon_structure(short_type);
  /** Give a pointer on the structures of a convex which is the direct
   *   product of the convexes represented by *pcvs1 and *pcvs2.
   */
  pconvex_structure convex_product_structure(pconvex_structure,
					     pconvex_structure);
  /** Give a pointer on the structures of a prism of dimension d.
   *   i.e. the direct product of a simplex of dimension d-1 and a segment.
   */
  inline pconvex_structure prism_structure(dim_type nc) { 
    return convex_product_structure(simplex_structure(nc-1),
				    simplex_structure(1));
  }

  /** Simplex structure with the Lagrange grid of degree k.
      @param n the simplex dimension.
      @param k the simplex degree.
  */
  pconvex_structure simplex_structure(dim_type n, short_type k);

  /// Generic convex with n global nodes
  pconvex_structure generic_dummy_structure(dim_type nc, size_type n,
					    size_type nf);

  //@}

  /*@}*/
}  /* end of namespace bgeot.                                             */


#endif /* BGEOT_CONVEX_STRUCTURE_H__                                      */
